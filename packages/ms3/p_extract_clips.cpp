/*
 * Copyright 2016-2017 Flatiron Institute, Simons Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "p_extract_clips.h"
#include "diskreadmda32.h"
#include "diskreadmda.h"
#include "extract_clips.h"
#include "pca.h"

#include <diskwritemda.h>

namespace P_extract_clips {
Mda32 extract_channels_from_chunk(const Mda32& chunk, const QList<int>& channels);
}

bool p_extract_clips(QStringList timeseries_list, QString event_times, const QList<int>& channels, QString clips_out, const QVariantMap& params)
{
    DiskReadMda32 X(2, timeseries_list);
    DiskReadMda ET(event_times);

    bigint M = X.N1();
    //bigint N = X.N2();
    bigint T = params["clip_size"].toInt();
    bigint L = ET.totalSize();

    if ((ET.N1()>1)&&(ET.N2()>1)) {
        qWarning() << "Error: the input (event times) must be a vector.";
        return false;
    }

    bigint M2 = M;
    if (!channels.isEmpty()) {
        M2 = channels.count();
    }

    if (!T) {
        qWarning() << "Unexpected: Clip size is zero.";
        return false;
    }

    bigint Tmid = (bigint)((T + 1) / 2) - 1;
    printf("Extracting clips (%ld,%ld,%ld) (%ld)...\n", M, T, L, M2);
    DiskWriteMda clips;
    clips.open(MDAIO_TYPE_FLOAT32, clips_out, M2, T, L);
    for (bigint i = 0; i < L; i++) {
        bigint t1 = ET.value(i) - Tmid;
        //bigint t2 = t1 + T - 1;
        Mda32 tmp;
        if (!X.readChunk(tmp, 0, t1, M, T)) {
            qWarning() << "Problem reading chunk in extract_clips";
            return false;
        }
        if (!channels.isEmpty()) {
            tmp = P_extract_clips::extract_channels_from_chunk(tmp, channels);
        }
        if (!clips.writeChunk(tmp, 0, 0, i)) {
            qWarning() << "Problem writing chunk" << i;
            return false;
        }
    }

    return true;
}

bool p_mv_extract_clips(QStringList timeseries_list, QString firings, const QList<int>& channels, QString clips_out, const QVariantMap& params)
{
    DiskReadMda32 X(2, timeseries_list);
    DiskReadMda FF(firings);

    bigint M = X.N1();
    //bigint N = X.N2();
    bigint T = params["clip_size"].toInt();
    bigint L = FF.N2();

    bigint M2 = M;
    if (!channels.isEmpty()) {
        M2 = channels.count();
    }

    if (!T) {
        qWarning() << "Unexpected: Clip size is zero.";
        return false;
    }

    bigint Tmid = (bigint)((T + 1) / 2) - 1;
    printf("Extracting clips (%ld,%ld,%ld) (%ld)...\n", M, T, L, M2);
    DiskWriteMda clips;
    clips.open(MDAIO_TYPE_FLOAT32, clips_out, M2, T, L);
    for (bigint i = 0; i < L; i++) {
        bigint t1 = FF.value(1,i) - Tmid;
        //bigint t2 = t1 + T - 1;
        Mda32 tmp;
        if (!X.readChunk(tmp, 0, t1, M, T)) {
            qWarning() << "Problem reading chunk in extract_clips";
            return false;
        }
        if (!channels.isEmpty()) {
            tmp = P_extract_clips::extract_channels_from_chunk(tmp, channels);
        }
        if (!clips.writeChunk(tmp, 0, 0, i)) {
            qWarning() << "Problem writing chunk" << i;
            return false;
        }
    }

    return true;
}

namespace P_extract_clips {
Mda32 extract_channels_from_chunk(const Mda32& X, const QList<int>& channels)
{
    //int M=X.N1();
    int T = X.N2();
    int M2 = channels.count();
    Mda32 ret(M2, T);
    for (int t = 0; t < T; t++) {
        for (int m2 = 0; m2 < M2; m2++) {
            ret.set(X.value(channels[m2] - 1, t), m2, t);
        }
    }
    return ret;
}
}




bool p_mv_extract_clips_features(QString timeseries_path, QString firings_path, QString features_out_path, int clip_size, int num_features, int subtract_mean)
{
    DiskReadMda32 X(timeseries_path);
    DiskReadMda F(firings_path);
    QVector<double> times;
    for (int i = 0; i < F.N2(); i++) {
        times << F.value(1, i);
    }
    Mda32 clips = extract_clips(X, times, clip_size);
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    Mda clips_reshaped(M * T, L);
    int NNN = M * T * L;
    for (int iii = 0; iii < NNN; iii++) {
        clips_reshaped.set(clips.get(iii), iii);
    }
    Mda CC, FF, sigma;
    pca(CC, FF, sigma, clips_reshaped, num_features, subtract_mean);
    return FF.write32(features_out_path);
}
