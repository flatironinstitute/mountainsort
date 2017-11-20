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
#include "globaltemplatecomputer.h"

#include <mda.h>
#include "omp.h"

class GlobalTemplateComputerPrivate {
public:
    GlobalTemplateComputer* q;
    int m_clip_size = 0; //should be set

    int m_num_threads = 1;
    QVector<double> m_times;
    QVector<int> m_labels;
    int m_K = 0;

    Mda m_waveform_sums;
    QVector<bigint> m_waveform_counts;
};

GlobalTemplateComputer::GlobalTemplateComputer()
{
    d = new GlobalTemplateComputerPrivate;
    d->q = this;
}

GlobalTemplateComputer::~GlobalTemplateComputer()
{
    delete d;
}

void GlobalTemplateComputer::setNumThreads(int num_threads)
{
    d->m_num_threads = num_threads;
}

void GlobalTemplateComputer::setClipSize(int clip_size)
{
    d->m_clip_size = clip_size;
}

void GlobalTemplateComputer::setTimesLabels(const QVector<double>& times, const QVector<int>& labels)
{
    d->m_times = times;
    d->m_labels = labels;
    d->m_K = MLCompute::max(labels);
}

void GlobalTemplateComputer::processTimeChunk(bigint t, const Mda32& X, bigint padding_left, bigint padding_right)
{
    int M = X.N1();
    int T = d->m_clip_size;
    int K = d->m_K;
    int Tmid = (int)((T + 1) / 2) - 1;
    if (d->m_waveform_counts.count() != d->m_K) {
        d->m_waveform_sums.allocate(M, T, d->m_K);
        d->m_waveform_counts = QVector<bigint>(d->m_K, 0);
    }

#pragma omp parallel num_threads(d->m_num_threads)
    {
        QVector<double> local_times;
        QVector<int> local_labels;
        Mda local_sums(M, T, K);
        QVector<bigint> local_counts(K, 0);
#pragma omp critical
        {
            for (bigint a = omp_get_thread_num(); a < d->m_times.count(); a += omp_get_num_threads()) {
                local_times << d->m_times[a];
                local_labels << d->m_labels[a];
            }
        }
#pragma omp for
        for (bigint i = 0; i < local_times.count(); i++) {
            double t0 = local_times[i];
            double t1 = t0 - t + padding_left;
            if ((padding_left <= t1) && (t1 < X.N2() - padding_left - padding_right)) {
                int k0 = local_labels[i];
                if (k0 > 0) {
                    Mda32 tmp;
                    X.getChunk(tmp, 0, t1 - Tmid, M, T);
                    for (int t = 0; t < T; t++) {
                        for (int m = 0; m < M; m++) {
                            local_sums.set(local_sums.get(m, t, k0 - 1) + tmp.get(m, t), m, t, k0 - 1);
                        }
                    }
                    local_counts[k0 - 1]++;
                }
            }
        }
#pragma omp critical(set_sums_and_counts)
        for (int kk = 1; kk <= K; kk++) {
            d->m_waveform_counts[kk - 1] += local_counts[kk - 1];
            Mda tmp;
            local_sums.getChunk(tmp, 0, 0, kk - 1, M, T, 1);
            for (int t = 0; t < T; t++) {
                for (int m = 0; m < M; m++) {
                    d->m_waveform_sums.set(d->m_waveform_sums.get(m, t, kk - 1) + tmp.get(m, t), m, t, kk - 1);
                }
            }
        }
    }
}

Mda32 GlobalTemplateComputer::templates() const
{
    int M = d->m_waveform_sums.N1();
    int T = d->m_waveform_sums.N2();
    int K = d->m_waveform_sums.N3();
    Mda32 templates(M, T, K);
    for (int kk = 1; kk <= K; kk++) {
        if (d->m_waveform_counts[kk - 1]) {
            for (int t = 0; t < T; t++) {
                for (int m = 0; m < M; m++) {
                    templates.set(d->m_waveform_sums.get(m, t, kk - 1) / d->m_waveform_counts[kk - 1], m, t, kk - 1);
                }
            }
        }
    }
    return templates;
}
