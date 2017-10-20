/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 2/22/2017
*******************************************************/

#include "p_extract_clips.h"
#include "diskreadmda32.h"
#include "diskreadmda.h"

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
