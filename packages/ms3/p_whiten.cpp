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
#include "p_whiten.h"

#include <QTime>
#include <diskreadmda32.h>
#include <diskwritemda.h>
#include <mda.h>
#include "pca.h"
#include "omp.h"
#include <cmath>
using std::floor;

namespace P_whiten {

void scale_for_quantization(Mda32& X, double quantization_unit);

double quantize(float X, double unit)
{
    return (floor(X / unit + 0.5)) * unit;
}

void quantize(bigint N, float* X, double unit)
{
    for (bigint i = 0; i < N; i++) {
        X[i] = quantize(X[i], unit);
    }
}
Mda32 extract_channels_from_chunk(const Mda32& X, const QList<int>& channels);
}

bool p_whiten(QString timeseries, QString timeseries_out, Whiten_opts &opts)
{
    (void)opts;

    DiskReadMda32 X(timeseries);
    bigint M = X.N1();
    bigint N = X.N2();

    bigint processing_chunk_size = 1e5; //changed from 1e7 to 1e5 by jfm (do ensure we are using all threads for short datasets with large numbers of channels)

    int overhead = 3; //maybe this should be 2
    double expected_peak_ram_bytes = overhead * processing_chunk_size * M * sizeof(float);
    opts.expected_peak_ram_mb = expected_peak_ram_bytes / (1024 * 1024);

    if (opts.requirements_only)
        return true;

    qDebug().noquote() << "Expected peak RAM (MB):" << opts.expected_peak_ram_mb;

    Mda XXt(M, M);
    double* XXtptr = XXt.dataPtr();
    bigint chunk_size = processing_chunk_size;
    if (N < processing_chunk_size) {
        chunk_size = N;
    }

    {
        QTime timer;
        timer.start();
        bigint num_timepoints_handled = 0;
#pragma omp parallel for
        for (bigint timepoint = 0; timepoint < N; timepoint += chunk_size) {
            Mda32 chunk;
#pragma omp critical(lock1)
            {
                if (!X.readChunk(chunk, 0, timepoint, M, qMin(chunk_size, N - timepoint))) {
                    qWarning() << "Problem reading chunk in whiten (1)";
                }
            }
            float* chunkptr = chunk.dataPtr();
            Mda XXt0(M, M);
            double* XXt0ptr = XXt0.dataPtr();
            for (bigint i = 0; i < chunk.N2(); i++) {
                bigint aa = M * i;
                bigint bb = 0;
                for (bigint m1 = 0; m1 < M; m1++) {
                    for (bigint m2 = 0; m2 < M; m2++) {
                        XXt0ptr[bb] += chunkptr[aa + m1] * chunkptr[aa + m2];
                        bb++;
                    }
                }
            }
#pragma omp critical(lock2)
            {
                bigint bb = 0;
                for (bigint m1 = 0; m1 < M; m1++) {
                    for (bigint m2 = 0; m2 < M; m2++) {
                        XXtptr[bb] += XXt0ptr[bb];
                        bb++;
                    }
                }
                num_timepoints_handled += qMin(chunk_size, N - timepoint);
                if ((timer.elapsed() > 5000) || (num_timepoints_handled == N)) {
                    printf("%ld/%ld (%d%%)\n", num_timepoints_handled, N, (int)(num_timepoints_handled * 1.0 / N * 100));
                    timer.restart();
                }
            }
        }
    }
    if (N > 1) {
        for (bigint ii = 0; ii < M * M; ii++) {
            XXtptr[ii] /= (N - 1);
        }
    }

    /*
    qDebug().noquote() << "Debug XXt";
    {
        int aaa=0;
        for (bigint m1 = 0; m1 < M; m1++) {
            for (bigint m2 = 0; m2 < M; m2++) {
                printf("%g ",XXtptr[aaa]);
                aaa++;
            }
            printf("\n");
        }
    }
    */

    //Mda AA = get_whitening_matrix(COV);
    Mda WW;
    whitening_matrix_from_XXt(WW, XXt); // the result is symmetric (assumed below)
    double* WWptr = WW.dataPtr();

    /*
    qDebug().noquote() << "Debug WW";
    {
        int aaa=0;
        for (bigint m1 = 0; m1 < M; m1++) {
            for (bigint m2 = 0; m2 < M; m2++) {
                printf("%g ",WWptr[aaa]);
                aaa++;
            }
            printf("\n");
        }
    }
    */

    DiskWriteMda Y;
    int dtype = MDAIO_TYPE_FLOAT32;
    if (opts.quantization_unit > 0)
        dtype = MDAIO_TYPE_INT16;
    Y.open(dtype, timeseries_out, M, N);
    {
        QTime timer;
        timer.start();
        bigint num_timepoints_handled = 0;
#pragma omp parallel for
        for (bigint timepoint = 0; timepoint < N; timepoint += chunk_size) {
            Mda32 chunk_in;
#pragma omp critical(lock1)
            {
                if (!X.readChunk(chunk_in, 0, timepoint, M, qMin(chunk_size, N - timepoint))) {
                    qWarning() << "Problem reading chunk in whiten (2)";
                }
            }
            float* chunk_in_ptr = chunk_in.dataPtr();
            Mda32 chunk_out(M, chunk_in.N2());
            float* chunk_out_ptr = chunk_out.dataPtr();
            for (bigint i = 0; i < chunk_in.N2(); i++) { // explicitly do mat-mat mult ... TODO replace w/ BLAS3
                bigint aa = M * i;
                bigint bb = 0;
                for (bigint m1 = 0; m1 < M; m1++) {
                    for (bigint m2 = 0; m2 < M; m2++) {
                        chunk_out_ptr[aa + m1] += chunk_in_ptr[aa + m2] * WWptr[bb]; // actually this does dgemm w/ WW^T
                        bb++; // but since symmetric, doesn't matter.
                    }
                }
            }
#pragma omp critical(lock2)
            {
                // The following is needed to make the output deterministic, due to a very tricky floating-point problem that I honestly could not track down
                // It has something to do with multiplying by very small values of WWptr[bb]. But I truly could not pinpoint the exact problem.
                P_whiten::quantize(chunk_out.totalSize(), chunk_out.dataPtr(), 0.0001);
                if (opts.quantization_unit > 0) {
                    P_whiten::scale_for_quantization(chunk_out, opts.quantization_unit);
                }
                if (!Y.writeChunk(chunk_out, 0, timepoint)) {
                    qWarning() << "Problem writing chunk in whiten";
                }
                num_timepoints_handled += qMin(chunk_size, N - timepoint);
                if ((timer.elapsed() > 5000) || (num_timepoints_handled == N)) {
                    printf("%ld/%ld (%d%%)\n", num_timepoints_handled, N, (int)(num_timepoints_handled * 1.0 / N * 100));
                    timer.restart();
                }
            }
        }
    }
    Y.close();

    return true;
}

bool p_compute_whitening_matrix(QStringList timeseries_list, const QList<int>& channels, QString whitening_matrix_out, Whiten_opts opts)
{
    (void)opts;

    DiskReadMda32 X0(2, timeseries_list);
    bigint M = X0.N1();
    bigint N = X0.N2();
    bigint M2 = M;
    if (!channels.isEmpty()) {
        M2 = channels.count();
    }
    qDebug().noquote() << "Computing whitening matrix: M/N" << M << N;

    bigint processing_chunk_size = 1e7;

    Mda XXt(M2, M2);
    double* XXtptr = XXt.dataPtr();
    bigint chunk_size = processing_chunk_size;
    if (N < processing_chunk_size) {
        chunk_size = N;
    }

    bigint timepoint = 0;
    while (timepoint < N) {
        QList<Mda32> chunks;
        while ((timepoint < N) && (chunks.count() < omp_get_max_threads())) {
            Mda32 chunk0;
            if (!X0.readChunk(chunk0, 0, timepoint, M, qMin(chunk_size, N - timepoint))) {
                qWarning() << "Problem reading chunk in compute whiten matrix";
                return false;
            }
            if (!channels.isEmpty()) {
                chunk0 = P_whiten::extract_channels_from_chunk(chunk0, channels);
            }
            chunks << chunk0;
            timepoint += chunk_size;
        }
        bigint num_chunks = chunks.count();
#pragma omp parallel
        {
#pragma omp for
            for (bigint i = 0; i < num_chunks; i++) {
                Mda32 chunk0;
#pragma omp critical
                {
                    chunk0 = chunks[i];
                }
                Mda XXt0(M2, M2);
                double* XXt0ptr = XXt0.dataPtr();
                float* chunkptr = chunk0.dataPtr();
                for (bigint i = 0; i < chunk0.N2(); i++) {
                    bigint aa = M2 * i;
                    bigint bb = 0;
                    for (bigint m1 = 0; m1 < M2; m1++) {
                        for (bigint m2 = 0; m2 < M2; m2++) {
                            XXt0ptr[bb] += chunkptr[aa + m1] * chunkptr[aa + m2];
                            bb++;
                        }
                    }
                }
#pragma omp critical(lock2)
                {
                    bigint bb = 0;
                    for (bigint m1 = 0; m1 < M2; m1++) {
                        for (bigint m2 = 0; m2 < M2; m2++) {
                            XXtptr[bb] += XXt0ptr[bb];
                            bb++;
                        }
                    }
                }
            }
        }
    }

    if (N > 1) {
        for (bigint ii = 0; ii < M2 * M2; ii++) {
            XXtptr[ii] /= (N - 1);
        }
    }

    //Mda AA = get_whitening_matrix(COV);
    Mda WW;
    whitening_matrix_from_XXt(WW, XXt); // the result is symmetric (assumed below)

    return WW.write64(whitening_matrix_out);
}

bool p_whiten_clips(QString clips_path, QString whitening_matrix, QString clips_out_path, Whiten_opts opts)
{
    (void)opts;
    Mda WW(whitening_matrix);
    DiskReadMda32 clips(clips_path);

    double* WWptr = WW.dataPtr();

    bigint M = clips.N1();
    bigint T = clips.N2();
    bigint L = clips.N3();

    qDebug().noquote() << "Whitening..." << M << T << L;

    int dtype = MDAIO_TYPE_FLOAT32;
    if (opts.quantization_unit > 0)
        dtype = MDAIO_TYPE_INT16;
    DiskWriteMda clips_out(dtype, clips_out_path, M, T, L);

    for (bigint i = 0; i < L; i++) {
        Mda32 chunk;
        if (!clips.readChunk(chunk, 0, 0, i, M, T, 1)) {
            qWarning() << "Problem reading chunk" << i;
            return false;
        }
        float* chunk_ptr = chunk.dataPtr();

        Mda32 chunk_out(M, T);
        float* chunk_out_ptr = chunk_out.dataPtr();
        {

            for (bigint t = 0; t < T; t++) {
                bigint aa = 0;
                bigint bb = M * t;
                for (bigint m1 = 0; m1 < M; m1++) {
                    for (bigint m2 = 0; m2 < M; m2++) {
                        chunk_out_ptr[bb + m1] += chunk_ptr[bb + m2] * WWptr[aa]; // actually this does dgemm w/ WW^T
                        aa++; // but since symmetric, doesn't matter.
                    }
                }
            }
        }
        if (opts.quantization_unit > 0) {
            P_whiten::scale_for_quantization(chunk_out, opts.quantization_unit);
        }
        if (!clips_out.writeChunk(chunk_out, 0, 0, i)) {
            qWarning() << "Problem writing chunk";
            return false;
        }
    }

    /*
    {
        // The following is needed to make the output deterministic, due to a very tricky floating-point problem that I honestly could not track down
        // It has something to do with multiplying by very small values of WWptr[bb]. But I truly could not pinpoint the exact problem.
        P_whiten::quantize(clips_out.totalSize(), clips_out.dataPtr(), 0.0001);
    }
    */

    return true;
}

bool p_apply_whitening_matrix(QString timeseries, QString whitening_matrix, QString timeseries_out, Whiten_opts opts)
{
    (void)opts;

    DiskReadMda32 X(timeseries);
    bigint M = X.N1();
    bigint N = X.N2();

    bigint processing_chunk_size = 1e7;
    bigint chunk_size = processing_chunk_size;
    if (N < processing_chunk_size) {
        chunk_size = N;
    }

    Mda WW(whitening_matrix);

    double* WWptr = WW.dataPtr();

    DiskWriteMda Y;
    int dtype = MDAIO_TYPE_FLOAT32;
    if (opts.quantization_unit > 0)
        dtype = MDAIO_TYPE_INT16;
    Y.open(dtype, timeseries_out, M, N);
    {
        QTime timer;
        timer.start();
        bigint num_timepoints_handled = 0;
#pragma omp parallel for
        for (bigint timepoint = 0; timepoint < N; timepoint += chunk_size) {
            Mda32 chunk_in;
#pragma omp critical(lock1)
            {
                if (!X.readChunk(chunk_in, 0, timepoint, M, qMin(chunk_size, N - timepoint))) {
                    qWarning() << "Problem reading chunk in whiten (3)";
                }
            }
            float* chunk_in_ptr = chunk_in.dataPtr();
            Mda32 chunk_out(M, chunk_in.N2());
            float* chunk_out_ptr = chunk_out.dataPtr();
            for (bigint i = 0; i < chunk_in.N2(); i++) { // explicitly do mat-mat mult ... TODO replace w/ BLAS3
                bigint aa = M * i;
                bigint bb = 0;
                for (bigint m1 = 0; m1 < M; m1++) {
                    for (bigint m2 = 0; m2 < M; m2++) {
                        chunk_out_ptr[aa + m1] += chunk_in_ptr[aa + m2] * WWptr[bb]; // actually this does dgemm w/ WW^T
                        bb++; // but since symmetric, doesn't matter.
                    }
                }
            }
#pragma omp critical(lock2)
            {
                // The following is needed to make the output deterministic, due to a very tricky floating-point problem that I honestly could not track down
                // It has something to do with multiplying by very small values of WWptr[bb]. But I truly could not pinpoint the exact problem.
                P_whiten::quantize(chunk_out.totalSize(), chunk_out.dataPtr(), 0.0001);
                if (opts.quantization_unit > 0) {
                    P_whiten::scale_for_quantization(chunk_out, opts.quantization_unit);
                }
                if (!Y.writeChunk(chunk_out, 0, timepoint)) {
                    qWarning() << "Problem writing chunk in apply whitening matrix";
                }
                num_timepoints_handled += qMin(chunk_size, N - timepoint);
                if ((timer.elapsed() > 5000) || (num_timepoints_handled == N)) {
                    printf("%ld/%ld (%d%%)\n", num_timepoints_handled, N, (int)(num_timepoints_handled * 1.0 / N * 100));
                    timer.restart();
                }
            }
        }
    }
    Y.close();

    return true;
}

namespace P_whiten {
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

void scale_for_quantization(Mda32& X, double quantization_unit)
{
    bigint N = X.totalSize();
    for (bigint i = 0; i < N; i++) {
        X.set((int)((X.get(i) / quantization_unit) + 0.5), i);
    }
}
}
