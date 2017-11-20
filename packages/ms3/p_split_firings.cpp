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
#include "p_split_firings.h"

#include <diskreadmda32.h>
#include <mda.h>

namespace P_split_firings {

struct Session {
    bigint start_timepoint = 0;
    bigint end_timepoint = 0;
    QVector<bigint> event_indices;
};
}

bool p_split_firings(QStringList timeseries_paths, QString firings_path, QStringList firings_out_paths)
{
    if (timeseries_paths.count() == 0)
        return true;

    QVector<bigint> start_timepoints;

    QList<P_split_firings::Session> sessions;
    bigint offset = 0;
    for (bigint i = 0; i < timeseries_paths.count(); i++) {
        DiskReadMda32 X(timeseries_paths[i]);
        P_split_firings::Session S;
        S.start_timepoint = offset;
        S.end_timepoint = S.start_timepoint + X.N2() - 1;
        offset += X.N2();
        sessions << S;
    }

    Mda firings(firings_path);
    bigint current_timeseries_ind = 0;
    for (bigint i = 0; i < firings.N2(); i++) {
        double time0 = firings.value(1, i);
        while ((time0 < sessions[current_timeseries_ind].start_timepoint) && (current_timeseries_ind - 1 >= 0)) {
            current_timeseries_ind--;
        }
        while ((time0 > sessions[current_timeseries_ind].end_timepoint) && (current_timeseries_ind + 1 < sessions.count())) {
            current_timeseries_ind++;
        }
        sessions[current_timeseries_ind].event_indices << i;
    }

    for (bigint i = 0; i < sessions.count(); i++) {
        Mda firings0(firings.N1(), sessions[i].event_indices.count());
        for (bigint j = 0; j < sessions[i].event_indices.count(); j++) {
            for (bigint r = 0; r < firings.N1(); r++) {
                firings0.setValue(firings.value(r, sessions[i].event_indices[j]), r, j);
            }
            double time0 = firings0.value(1, j);
            firings0.setValue(time0 - sessions[i].start_timepoint, 1, j);
        }
        if (!firings0.write64(firings_out_paths.value(i)))
            return false;
    }

    return true;
}

namespace P_split_firings {
}

bool p_extract_firings(QString firings, const QSet<int>& clusters, QString firings_out)
{
    Mda FF1(firings);
    QVector<bigint> inds_to_use;
    for (bigint i = 0; i < FF1.N2(); i++) {
        int label = FF1.value(2, i);
        if (clusters.contains(label))
            inds_to_use << i;
    }
    Mda FF2(FF1.N1(), inds_to_use.count());
    for (bigint i = 0; i < FF2.N2(); i++) {
        for (int j = 0; j < FF2.N1(); j++) {
            FF2.setValue(FF1.value(j, inds_to_use[i]), j, i);
        }
    }
    return FF2.write64(firings_out);
}
