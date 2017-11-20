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

#include "compute_templates_0.h"
#include <math.h>
#include "get_sort_indices.h"
#include "omp.h"

Mda compute_templates_0(const DiskReadMda& X, Mda& firings, int clip_size)
{
    QVector<double> times;
    QVector<int> labels;
    int L = firings.N2();
    for (int i = 0; i < L; i++) {
        times << firings.value(1, i);
        labels << (int)firings.value(2, i);
    }
    return compute_templates_0(X, times, labels, clip_size);
}

Mda32 compute_templates_0(const DiskReadMda32& X, Mda& firings, int clip_size)
{
    QVector<double> times;
    QVector<int> labels;
    int L = firings.N2();
    for (int i = 0; i < L; i++) {
        times << firings.value(1, i);
        labels << (int)firings.value(2, i);
    }
    return compute_templates_0(X, times, labels, clip_size);
}

Mda compute_templates_0(const DiskReadMda& X, const QVector<double>& times, const QVector<int>& labels, int clip_size)
{
    int M = X.N1();
    int T = clip_size;
    int L = times.count();

    int K = MLCompute::max<int>(labels);

    int Tmid = (int)((T + 1) / 2) - 1;

    Mda templates(M, T, K);
    QList<int> counts;
    for (int k = 0; k < K; k++)
        counts << 0;
    for (int i = 0; i < L; i++) {
        int k = labels[i];
        bigint t0 = (bigint)(times[i] + 0.5);
        if (k >= 1) {
            Mda X0;
            X.readChunk(X0, 0, t0 - Tmid, M, T);
            double* Xptr = X0.dataPtr();
            double* Tptr = templates.dataPtr(0, 0, k - 1);
            for (int i = 0; i < M * T; i++) {
                Tptr[i] += Xptr[i];
            }
            counts[k - 1]++;
        }
    }
    for (int k = 0; k < K; k++) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                if (counts[k]) {
                    templates.set(templates.get(m, t, k) / counts[k], m, t, k);
                }
            }
        }
    }

    return templates;
}

void get_sums_and_counts_for_templates(Mda& sums, Mda& counts, const Mda32& X, bigint t_offset, const QVector<double>& times, const QVector<int>& labels, int clip_size, int K)
{
    int M = X.N1();
    bigint N = X.N2();
    int T = clip_size;
    int Tmid = (int)((T + 1) / 2) - 1;
    sums.allocate(M, T, K);
    counts.allocate(1, K);
    for (bigint i = 0; i < times.count(); i++) {
        bigint t = times[i] - t_offset;
        if ((t >= clip_size) && (t < N - clip_size)) {
            int k = labels[i];
            if ((k >= 1) && (k <= K)) {
                Mda32 clip;
                X.getChunk(clip, 0, t - Tmid, M, T);
                for (int t = 0; t < T; t++) {
                    for (int m = 0; m < M; m++) {
                        sums.setValue(sums.value(m, t, k - 1) + clip.value(m, t), m, t, k - 1);
                    }
                }
                counts.set(counts.get(k - 1) + 1, k - 1);
            }
        }
    }
}

Mda32 compute_templates_in_parallel(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int clip_size)
{
    int M = X.N1();
    bigint N = X.N2();
    int T = clip_size;
    int K = MLCompute::max<int>(labels);

    Mda sums(M, T, K);
    Mda counts(1, K);

    bigint chunk_size = 1e5;
#pragma omp parallel for
    for (bigint t = 0; t < N; t += chunk_size) {
        Mda32 chunk;
#pragma omp critical(compute_templates_in_parallel1)
        {
            X.readChunk(chunk, 0, t - clip_size, M, chunk_size + 2 * clip_size);
        }
        Mda sums0;
        Mda counts0;
        get_sums_and_counts_for_templates(sums0, counts0, chunk, t - clip_size, times, labels, clip_size, K);
#pragma omp critical(compute_templates_in_parallel2)
        {
            for (bigint i = 0; i < M * T * K; i++) {
                sums.set(sums.get(i) + sums0.get(i), i);
            }
            for (int i = 0; i < K; i++) {
                counts.set(counts.get(i) + counts0.get(i), i);
            }
        }
    }
    Mda32 ret(M, T, K);
    for (int k = 0; k < K; k++) {
        if (counts.get(k)) {
            for (int t = 0; t < T; t++) {
                for (int m = 0; m < M; m++) {
                    ret.setValue(sums.value(m, t, k) / counts.get(k), m, t, k);
                }
            }
        }
    }
    return ret;
}

Mda32 compute_templates_0(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int clip_size)
{
    int M = X.N1();
    int T = clip_size;
    int L = times.count();

    int K = MLCompute::max<int>(labels);

    int Tmid = (int)((T + 1) / 2) - 1;

    Mda32 templates(M, T, K);
    QList<int> counts;
    for (int k = 0; k < K; k++)
        counts << 0;
    for (int i = 0; i < L; i++) {
        int k = labels[i];
        bigint t0 = (bigint)(times[i] + 0.5);
        if (k >= 1) {
            Mda32 X0;
            X.readChunk(X0, 0, t0 - Tmid, M, T);
            dtype32* Xptr = X0.dataPtr();
            dtype32* Tptr = templates.dataPtr(0, 0, k - 1);
            for (int i = 0; i < M * T; i++) {
                Tptr[i] += Xptr[i];
            }
            counts[k - 1]++;
        }
    }
    for (int k = 0; k < K; k++) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                if (counts[k]) {
                    templates.set(templates.get(m, t, k) / counts[k], m, t, k);
                }
            }
        }
    }

    return templates;
}

void compute_templates_stdevs(Mda& templates, Mda& stdevs, DiskReadMda& X, const QVector<double>& times, const QVector<int>& labels, int clip_size)
{
    int M = X.N1();
    int T = clip_size;
    bigint L = times.count();

    int K = MLCompute::max<int>(labels);

    int Tmid = (int)((T + 1) / 2) - 1;

    Mda sums(M, T, K);
    Mda sumsqrs(M, T, K);
    QList<bigint> counts;
    for (int k = 0; k < K; k++)
        counts << 0;
    for (bigint i = 0; i < L; i++) {
        int k = labels[i];
        bigint t0 = (bigint)(times[i] + 0.5);
        if (k >= 1) {
            Mda X0;
            X.readChunk(X0, 0, t0 - Tmid, M, T);
            double* Xptr = X0.dataPtr();
            double* sum_ptr = sums.dataPtr(0, 0, k - 1);
            double* sumsqr_ptr = sumsqrs.dataPtr(0, 0, k - 1);
            for (int i = 0; i < M * T; i++) {
                sum_ptr[i] += Xptr[i];
                sumsqr_ptr[i] += Xptr[i] * Xptr[i];
            }
            counts[k - 1]++;
        }
    }

    templates.allocate(M, T, K);
    stdevs.allocate(M, T, K);
    for (int k = 0; k < K; k++) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                if (counts[k] >= 2) {
                    double sum0 = sums.get(m, t, k);
                    double sumsqr0 = sumsqrs.get(m, t, k);
                    templates.set(sum0 / counts[k], m, t, k);
                    stdevs.set(sqrt(sumsqr0 / counts[k] - (sum0 * sum0) / (counts[k] * counts[k])), m, t, k);
                }
            }
        }
    }
}
