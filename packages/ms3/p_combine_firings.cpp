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
#include "p_combine_firings.h"

#include <diskreadmda.h>
#include <mda.h>
#include "get_sort_indices.h"

bool p_combine_firings(QStringList firings_list, QString firings_out, bool increment_labels)
{
    QVector<int> all_central_channels;
    QVector<double> all_times;
    QVector<int> all_labels;
    QVector<double> all_extra_values;

    if (firings_list.count() == 0) {
        qWarning() << "Firings list is empty.";
        return false;
    }

    DiskReadMda firings0(firings_list.value(0));
    int R = firings0.N1();

    int max_label = 0;
    int label_offset = 0;
    for (bigint i = 0; i < firings_list.count(); i++) {
        Mda F(firings_list[i]);
        for (bigint j = 0; j < F.N2(); j++) {
            int label0 = F.value(2, j);
            if (label0 > 0) {
                all_central_channels << F.value(0, j);
                all_times << F.value(1, j);
                all_labels << label_offset + label0;
                for (int r = 3; r < R; r++)
                    all_extra_values << F.value(r, j);
                if (label_offset + label0 > max_label)
                    max_label = label_offset + label0;
            }
        }
        if (increment_labels)
            label_offset = max_label;
    }

    bigint L = all_times.count();
    QList<bigint> sort_inds = get_sort_indices_bigint(all_times);

    Mda ret(R, L);
    for (bigint j = 0; j < L; j++) {
        ret.setValue(all_central_channels[sort_inds[j]], 0, j);
        ret.setValue(all_times[sort_inds[j]], 1, j);
        ret.setValue(all_labels[sort_inds[j]], 2, j);
        for (int r = 3; r < R; r++) {
            ret.setValue(all_extra_values[sort_inds[j] * (R - 3) + (r - 3)], r, j);
        }
    }
    return ret.write64(firings_out);
}
