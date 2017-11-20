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
#include "dimension_reduce_clips.h"
#include "pca.h"

#include <QTime>

void dimension_reduce_clips(Mda32& ret, const Mda32& clips, bigint num_features_per_channel, bigint max_samples)
{
    bigint M = clips.N1();
    bigint T = clips.N2();
    bigint L = clips.N3();
    const float* clips_ptr = clips.constDataPtr();

    qDebug().noquote() << QString("Dimension reduce clips %1x%2x%3").arg(M).arg(T).arg(L);

    ret.allocate(M, num_features_per_channel, L);
    float* retptr = ret.dataPtr();
    for (bigint m = 0; m < M; m++) {
        Mda32 reshaped(T, L);
        float* reshaped_ptr = reshaped.dataPtr();
        bigint aa = 0;
        bigint bb = m;
        for (bigint i = 0; i < L; i++) {
            for (bigint t = 0; t < T; t++) {
                //reshaped.set(clips.get(bb),aa);
                reshaped_ptr[aa] = clips_ptr[bb];
                aa++;
                bb += M;
                //reshaped.setValue(clips.value(m, t, i), t, i);
            }
        }
        Mda32 CC, FF, sigma;
        QTime timerA;
        timerA.start();
        pca_subsampled(CC, FF, sigma, reshaped, num_features_per_channel, false, max_samples);
        float* FF_ptr = FF.dataPtr();
        aa = 0;
        bb = m;
        for (bigint i = 0; i < L; i++) {
            for (bigint a = 0; a < num_features_per_channel; a++) {
                //ret.setValue(FF.value(a, i), m, a, i);
                //ret.set(FF.get(aa),bb);
                retptr[bb] = FF_ptr[aa];
                aa++;
                bb += M;
            }
        }
    }
}

void dimension_reduce_clips(QString clips_path, QString reduced_clips_path, bigint num_features_per_channel, bigint max_samples)
{
    Mda32 clips(clips_path);
    Mda32 reduced_clips;
    dimension_reduce_clips(reduced_clips, clips, num_features_per_channel, max_samples);
    reduced_clips.write32(reduced_clips_path);
}
