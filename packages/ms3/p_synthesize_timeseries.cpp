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
#include "p_synthesize_timeseries.h"

#include "diskwritemda.h"
#include "get_sort_indices.h"

#include <QTime>
#include <diskreadmda.h>
#include <random>

namespace Synthesize_timeseries {
void generate_randn(size_t N, float* X);
}

using namespace Synthesize_timeseries;

bool p_synthesize_timeseries(QString firings_in, QString waveforms_in, QString timeseries_out, P_synthesize_timeseries_opts opts)
{
    DiskReadMda firings(firings_in);
    Mda32 waveforms(waveforms_in);

    if (waveforms.N2() % opts.waveform_upsample_factor != 0) {
        qWarning() << "Waveforms dimensions not consistent with waveform_upsample_factor" << waveforms.N2() << opts.waveform_upsample_factor;
        return false;
    }

    bigint L = firings.N2();
    int M = waveforms.N1();
    int T = waveforms.N2() / opts.waveform_upsample_factor;
    int K = waveforms.N3();
    int Tmid = (int)((T + 1) / 2) - 1;

    QVector<double> times0(L);
    QVector<int> labels0(L);
    for (bigint i = 0; i < L; i++) {
        times0[i] = firings.value(1, i);
        labels0[i] = firings.value(2, i);
    }
    QList<bigint> inds = get_sort_indices_bigint(times0);
    QVector<double> times(L);
    QVector<int> labels(L);
    for (bigint i = 0; i < L; i++) {
        times[i] = times0[inds[i]];
        labels[i] = labels0[inds[i]];
    }

    if (!opts.duration) {
        opts.duration = times.value(times.count() - 1) + T * 2;
    }
    bigint NN = opts.duration;
    qDebug().noquote() << QString("Generating timeseries with %1 timepoints").arg(NN);

    DiskWriteMda Y(MDAIO_TYPE_FLOAT32, timeseries_out, M, NN);
    bigint chunk_size = 1e7;
    bigint ii = 0;
    QTime timer;
    timer.start();
    for (bigint t = 0; t < NN; t += chunk_size) {
        if (timer.elapsed() > 3000) {
            qDebug().noquote() << QString("Handling time %1 of %2 (%3%)").arg(t).arg(NN).arg(t * 100.0 / NN);
            timer.restart();
        }
        bigint N0 = chunk_size;
        if (t + N0 >= NN)
            N0 = NN - t;
        Mda32 chunk(M, N0);
        if (opts.noise_level > 0) {
            generate_randn(M * N0, chunk.dataPtr());
            for (bigint i = 0; i < M * N0; i++) {
                chunk.set(chunk.get(i) * opts.noise_level, i);
            }
        }
        while ((ii > 0) && (times[ii - 1] >= t - T))
            ii--;
        while (1) {
            if (ii >= times.count())
                break;
            double t0 = times[ii];
            if (t0 >= t + N0 + T)
                break;
            int k = labels[ii];
            if ((1 <= k) && (k <= K)) {
                Mda32 waveform0;
                waveforms.getChunk(waveform0, 0, 0, k - 1, M, T * opts.waveform_upsample_factor, 1);
                int offset = (int)((t0 - (bigint)t0) * opts.waveform_upsample_factor);
                for (int a = 0; a < T; a++) {
                    if ((a + t0 - t - Tmid >= 0) && (a + t0 - t - Tmid < chunk.N2())) {
                        for (int m = 0; m < M; m++) {
                            double val = waveform0.value(m, a * opts.waveform_upsample_factor + offset);
                            chunk.setValue(chunk.value(m, a + t0 - t - Tmid) + val, m, a + t0 - t - Tmid);
                        }
                    }
                }
            }
            ii++;
        }
        Y.writeChunk(chunk, 0, t);
    }

    return true;
}

namespace Synthesize_timeseries {
void generate_randn(size_t N, float* X)
{
    std::random_device rd; //to generate a non-deterministic integer (I believe)
    std::mt19937 e2(rd()); //the random number generator, seeded by rd() (I believe)
    std::normal_distribution<> dist(0, 1); //to generate random numbers from normal distribution
    for (size_t n = 0; n < N; n++) {
        X[n] = dist(e2);
    }
}
}
