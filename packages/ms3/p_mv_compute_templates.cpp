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
#include "p_mv_compute_templates.h"

#include <diskreadmda.h>
#include "compute_templates_0.h"

/// TODO 0.9.1 #define USE_TASK_PROGRESS, and don't use this outside of the mountainview gui -- unnecessary dependence AND risk

bool mv_compute_templates(const QString& timeseries_path, const QString& firings_path, const QString& templates_out_path, const QString& stdevs_out_path, int clip_size)
{
    DiskReadMda X(timeseries_path);
    if (X.N2() <= 1) {
        return false;
    }
    DiskReadMda firings(firings_path);
    QVector<double> times;
    QVector<int> labels;
    for (bigint i = 0; i < firings.N2(); i++) {
        times << firings.value(1, i);
        labels << (int)firings.value(2, i);
    }
    Mda templates, stdevs;
    compute_templates_stdevs(templates, stdevs, X, times, labels, clip_size);
    templates.write32(templates_out_path);
    stdevs.write32(stdevs_out_path);
    return true;
}

QList<bool> evenly_distributed_to_use(int num, int num_to_use)
{
    QList<bool> ret;
    for (int i = 0; i < num; i++) {
        ret << false;
    }
    double stride = num * 1.0 / num_to_use; //will be greater than 1
    double j = 0;
    for (int i = 0; i < num_to_use; i++) {
        ret[(int)j] = true;
        j += stride;
    }
    return ret;
}

bool mv_subfirings(QString firings_path, QString firings_out_path, QVector<int> labels, int max_per_label)
{
    QMap<int, int> counts;
    DiskReadMda F(firings_path);
    QSet<int> set;
    foreach (int label, labels) {
        set.insert(label);
    }
    QList<int> inds;
    QList<bool> to_use;
    for (int j = 0; j < F.N2(); j++) {
        int label = (int)F.value(2, j);
        if (set.contains(label)) {
            inds << j;
            to_use << true;
            counts[label]++;
        }
    }

    if (max_per_label) {
        int K = MLCompute::max<int>(labels);
        for (int k = 0; k <= K; k++) {
            if (counts[k] > max_per_label) {
                QList<bool> to_use_k = evenly_distributed_to_use(counts[k], max_per_label);
                int jj = 0;
                for (int i = 0; i < inds.count(); i++) {
                    int label = (int)F.value(2, inds[i]);
                    if (label == k) {
                        to_use[i] = to_use_k[jj];
                        jj++;
                    }
                }
            }
        }
    }

    QList<int> inds2;
    for (int i = 0; i < inds.count(); i++) {
        if (to_use[i])
            inds2 << inds[i];
    }

    Mda out(F.N1(), inds2.count());
    for (int i = 0; i < inds2.count(); i++) {
        for (int j = 0; j < F.N1(); j++) {
            out.setValue(F.value(j, inds2[i]), j, i);
        }
    }
    return out.write64(firings_out_path);
}
