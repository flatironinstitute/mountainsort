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
#include "p_create_multiscale_timeseries.h"
#include <QStringList>
#include <QFile>
#include <QTime>
#include "diskreadmda.h"
#include "diskwritemda.h"
#include "mlutil.h"

bigint smallest_power_of_3_larger_than(bigint N);
bool downsample_min(const DiskReadMda& X, QString out_fname, bigint N);
bool downsample_max(const DiskReadMda& X, QString out_fname, bigint N);
bool write_concatenation(QStringList input_fnames, QString output_fname);

bool p_create_multiscale_timeseries(QString path_in, QString path_out, QString tempdir)
{
    DiskReadMda X(path_in);
    X.reshape(X.N1(), X.N2() * X.N3()); //to handle the case of clips (3D array)

    bigint N = smallest_power_of_3_larger_than(X.N2());

    QStringList tmp_file_names;
    QString prev_min_fname;
    QString prev_max_fname;
    for (bigint ds_factor = 3; ds_factor <= N; ds_factor *= 3) {
        printf("ds_factor = %ld / %ld\n", ds_factor, N);
        QString rnd=MLUtil::makeRandomId(10);
        QString min_fname = tempdir+"/"+"min."+rnd;
        QString max_fname = tempdir+"/"+"max."+rnd;
        if (ds_factor == 3) {
            if (!downsample_min(X, min_fname, N)) {
                printf("Problem in downsample_min\n");
                printf("Removing temporary files\n");
                foreach (QString fname, tmp_file_names) {
                    QFile::remove(fname);
                }
                return false;
            }
            if (!downsample_max(X, max_fname, N)) {
                printf("Problem in downsample_max\n");
                printf("Removing temporary files\n");
                foreach (QString fname, tmp_file_names) {
                    QFile::remove(fname);
                }
                return false;
            }
        }
        else {
            if (!downsample_min(DiskReadMda(prev_min_fname), min_fname, (N * 3) / ds_factor)) {
                printf("Problem in downsample_min *\n");
                printf("Removing temporary files\n");
                foreach (QString fname, tmp_file_names) {
                    QFile::remove(fname);
                }
                return false;
            }
            if (!downsample_max(DiskReadMda(prev_max_fname), max_fname, (N * 3) / ds_factor)) {
                printf("Problem in downsample_max *\n");
                printf("Removing temporary files\n");
                foreach (QString fname, tmp_file_names) {
                    QFile::remove(fname);
                }
                return false;
            }
        }
        tmp_file_names << min_fname;
        tmp_file_names << max_fname;
        prev_min_fname = min_fname;
        prev_max_fname = max_fname;
        printf("...................\n");
    }

    printf("Writing concatenation...\n");
    if (!write_concatenation(tmp_file_names, path_out)) {
        printf("Problem in write_concatenation\n");
        printf("Removing temporary files\n");
        foreach (QString fname, tmp_file_names) {
            QFile::remove(fname);
        }
        return false;
    }
    printf("Removing temporary files\n");
    foreach (QString fname, tmp_file_names) {
        QFile::remove(fname);
    }

    return true;
}

bigint smallest_power_of_3_larger_than(bigint N)
{
    bigint ret = 1;
    while (ret < N) {
        ret *= 3;
    }
    return ret;
}

bool downsample_min(const DiskReadMda& X, QString out_fname, bigint N)
{
    DiskWriteMda Y;
    if (!Y.open(MDAIO_TYPE_FLOAT32, out_fname, X.N1(), N / 3)) {
        qWarning() << "Unable to open output file in downsample_min: " + out_fname;
        return false;
    }

    /// TODO choose chunk_size sensibly
    bigint chunk_size = 3 * 5000;
    QTime timer;
    timer.start();
    for (bigint ii = 0; ii < N; ii += chunk_size) {
        if (timer.elapsed() > 5000) {
            printf("downsample_min %ld/%ld (%ld%%)\n", ii, N, (bigint)(ii * 1.0 / N * 100));
            timer.restart();
        }
        bigint size0 = qMin(chunk_size, N - ii);
        Mda A;
        X.readChunk(A, 0, ii, X.N1(), size0);
        Mda B(X.N1(), size0 / 3);
        for (bigint jj = 0; jj < size0 / 3; jj++) {
            for (bigint m = 0; m < X.N1(); m++) {
                double val1 = A.value(m, jj * 3);
                double val2 = A.value(m, jj * 3 + 1);
                double val3 = A.value(m, jj * 3 + 2);
                double val = qMin(qMin(val1, val2), val3);
                B.setValue(val, m, jj);
            }
        }
        Y.writeChunk(B, 0, ii / 3);
    }

    return true;
}

bool downsample_max(const DiskReadMda& X, QString out_fname, bigint N)
{
    DiskWriteMda Y;
    if (!Y.open(MDAIO_TYPE_FLOAT32, out_fname, X.N1(), N / 3)) {
        qWarning() << "Unable to open output file in downsample_min: " + out_fname;
        return false;
    }

    /// TODO choose chunk_size sensibly
    bigint chunk_size = 3 * 5000;
    QTime timer;
    timer.start();
    for (bigint ii = 0; ii < N; ii += chunk_size) {
        if (timer.elapsed() > 5000) {
            printf("downsample_max %ld/%ld (%ld%%)\n", ii, N, (bigint)(ii * 1.0 / N * 100));
            timer.restart();
        }
        bigint size0 = qMin(chunk_size, N - ii);
        Mda A;
        X.readChunk(A, 0, ii, X.N1(), size0);
        Mda B(X.N1(), size0 / 3);
        for (bigint jj = 0; jj < size0 / 3; jj++) {
            for (bigint m = 0; m < X.N1(); m++) {
                double val1 = A.value(m, jj * 3);
                double val2 = A.value(m, jj * 3 + 1);
                double val3 = A.value(m, jj * 3 + 2);
                double val = qMax(qMax(val1, val2), val3);
                B.setValue(val, m, jj);
            }
        }
        Y.writeChunk(B, 0, ii / 3);
    }

    return true;
}

bool write_concatenation(QStringList input_fnames, QString output_fname)
{
    bigint M = 1, N = 0;
    foreach (QString fname, input_fnames) {
        DiskReadMda X(fname);
        M = X.N1();
        N += X.N2();
    }
    DiskWriteMda Y;
    if (!Y.open(MDAIO_TYPE_FLOAT32, output_fname, M, N)) {
        qWarning() << "Unable to open output file: " + output_fname;
        return false;
    }
    bigint offset = 0;
    foreach (QString fname, input_fnames) {
        DiskReadMda X(fname);

        /// TODO choose chunk_size sensibly
        bigint chunk_size = 10000;
        for (bigint ii = 0; ii < X.N2(); ii += chunk_size) {
            bigint size0 = qMin((bigint)chunk_size, (bigint)(X.N2() - ii));
            Mda tmp;
            X.readChunk(tmp, 0, ii, M, size0);
            Y.writeChunk(tmp, 0, offset);
            offset += size0;
        }
    }
    Y.close();

    return true;
}
