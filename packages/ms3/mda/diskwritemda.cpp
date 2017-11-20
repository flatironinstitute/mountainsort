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
#include "diskwritemda.h"
#include "mdaio.h"

#include <QFile>
#include <QString>
#include <mda32.h>
#include "mda.h"
#include <QDebug>

class DiskWriteMdaPrivate {
public:
    DiskWriteMda* q;
    QString m_path;
    MDAIO_HEADER m_header;
    FILE* m_file;
    bool m_requires_rename = false;

    int determine_ndims(bigint N1, bigint N2, bigint N3, bigint N4, bigint N5, bigint N6);
};

DiskWriteMda::DiskWriteMda()
{
    d = new DiskWriteMdaPrivate;
    d->q = this;
    d->m_file = 0;
}

DiskWriteMda::DiskWriteMda(int data_type, const QString& path, bigint N1, bigint N2, bigint N3, bigint N4, bigint N5, bigint N6)
{
    d = new DiskWriteMdaPrivate;
    d->q = this;
    d->m_file = 0;
    this->open(data_type, path, N1, N2, N3, N4, N5, N6);
}

DiskWriteMda::~DiskWriteMda()
{
    close();
    delete d;
}

bool DiskWriteMda::open(int data_type, const QString& path, bigint N1, bigint N2, bigint N3, bigint N4, bigint N5, bigint N6)
{
    if (d->m_file) {
        qWarning() << "Error in DiskWriteMda::open -- cannot open the file twice";
        return false; //can't open twice!
    }

    if (QFile::exists(path)) {
        if (!QFile::remove(path)) {
            qWarning() << "Unable to remove file in diskwritemda::open" << path;
            return false;
        }
    }

    d->m_path = path;

    d->m_header.data_type = data_type;
    for (int i = 0; i < MDAIO_MAX_DIMS; i++)
        d->m_header.dims[i] = 1;
    d->m_header.dims[0] = N1;
    d->m_header.dims[1] = N2;
    d->m_header.dims[2] = N3;
    d->m_header.dims[3] = N4;
    d->m_header.dims[4] = N5;
    d->m_header.dims[5] = N6;
    d->m_header.num_dims = d->determine_ndims(N1, N2, N3, N4, N5, N6);

    d->m_file = fopen((path + ".tmp").toLatin1().data(), "wb");
    d->m_requires_rename = true;

    if (!d->m_file) {
        qWarning() << "Error in DiskWriteMda::open -- problem in fopen: " + path + ".tmp";
        return false;
    }

    bigint NN = N1 * N2 * N3 * N4 * N5 * N6;
    //bigint buf_size = 1e6;

    //write the header
    mda_write_header(&d->m_header, d->m_file);

    /*
    //fill it all with zeros!
    float* zeros = (float*)malloc(sizeof(float) * buf_size);
    for (bigint i = 0; i < buf_size; i++)
        zeros[i] = 0;
    bigint i = 0;
    while (i < NN) {
        bigint num_to_write = NN - i;
        if (num_to_write > buf_size)
            num_to_write = buf_size;
        mda_write_float32(zeros, &d->m_header, num_to_write, d->m_file);
        i += buf_size;
    }
    free(zeros);
    */

    fseeko(d->m_file, d->m_header.header_size + d->m_header.num_bytes_per_entry * NN - 1, SEEK_SET);
    unsigned char zero = 0;
    fwrite(&zero, 1, 1, d->m_file);

    return true;
}

bool DiskWriteMda::open(const QString& path)
{
    if (d->m_file)
        return false; //can't open twice!

    d->m_path = path;

    d->m_requires_rename = false;
    d->m_file = fopen(path.toLatin1().data(), "r+"); //open file for update, both read and write
    mda_read_header(&d->m_header, d->m_file);

    if (!d->m_file)
        return false;

    return true;
}

void DiskWriteMda::close()
{
    if (d->m_file) {
        fclose(d->m_file);
        if (d->m_requires_rename) {
            if (!QFile::rename(d->m_path + ".tmp", d->m_path)) {
                qWarning() << "Unable to rename file in diskwritemda::open" << d->m_path + ".tmp" << d->m_path;
            }
        }
        d->m_file = 0;
    }
}

bigint DiskWriteMda::N1()
{
    if (!d->m_file)
        return 0;
    return d->m_header.dims[0];
}

bigint DiskWriteMda::N2()
{
    if (!d->m_file)
        return 0;
    return d->m_header.dims[1];
}

bigint DiskWriteMda::N3()
{
    if (!d->m_file)
        return 0;
    return d->m_header.dims[2];
}

bigint DiskWriteMda::N4()
{
    if (!d->m_file)
        return 0;
    return d->m_header.dims[3];
}

bigint DiskWriteMda::N5()
{
    if (!d->m_file)
        return 0;
    return d->m_header.dims[4];
}

bigint DiskWriteMda::N6()
{
    if (!d->m_file)
        return 0;
    return d->m_header.dims[5];
}

bigint DiskWriteMda::totalSize()
{
    return N1() * N2() * N3() * N4() * N5() * N6();
}

bool DiskWriteMda::writeChunk(Mda& X, bigint i)
{
    if (!d->m_file)
        return false;
    fseeko(d->m_file, d->m_header.header_size + d->m_header.num_bytes_per_entry * i, SEEK_SET);
    bigint size = X.totalSize();
    if (i + size > this->totalSize())
        size = this->totalSize() - i;
    if (size > 0) {
        if (!mda_write_float64(X.dataPtr(), &d->m_header, size, d->m_file))
            return false;
    }
    return true;
}

bool DiskWriteMda::writeChunk(Mda& X, bigint i1, bigint i2)
{
    if ((X.N1() == N1()) && (i1 == 0)) {
        return writeChunk(X, i1 + this->N1() * i2);
    }
    else {
        qWarning() << "dims:" << d->m_header.dims[0] << d->m_header.dims[1] << d->m_file;
        qWarning() << "This case not yet supported in 2d writeChunk" << X.N1() << X.N2() << N1() << N2() << i1 << i2;
        return false;
    }
}

bool DiskWriteMda::writeChunk(Mda& X, bigint i1, bigint i2, bigint i3)
{
    if ((i3 == 0) && (X.N3() == 1) && (X.N3() == 1))
        return writeChunk(X, i1, i2);
    if ((X.N1() == N1()) && (X.N2() == N2()) && (i1 == 0) && (i2 == 0)) {
        return writeChunk(X, i1 + this->N1() * i2 + this->N1() * this->N2() * i3);
    }
    else {
        qWarning() << "This case not yet supported in 3d writeSubArray" << X.N1() << X.N2() << X.N3() << N1() << N2() << N3() << i1 << i2 << i3;
        return false;
    }
}

bool DiskWriteMda::writeChunk(Mda32& X, bigint i)
{
    if (!d->m_file)
        return false;
    fseeko(d->m_file, d->m_header.header_size + d->m_header.num_bytes_per_entry * i, SEEK_SET);
    bigint size = X.totalSize();
    if (i + size > this->totalSize())
        size = this->totalSize() - i;
    if (size > 0) {
        return mda_write_float32(X.dataPtr(), &d->m_header, size, d->m_file);
    }
    else {
        qWarning() << "size is zero in writeChunk";
        return false;
    }
}

bool DiskWriteMda::writeChunk(Mda32& X, bigint i1, bigint i2)
{
    if ((X.N1() == N1()) && (i1 == 0)) {
        return writeChunk(X, i1 + this->N1() * i2);
    }
    else {
        qWarning() << "dims:" << d->m_header.dims[0] << d->m_header.dims[1] << d->m_file;
        qWarning() << "This case not yet supported in 2d writeSubArray" << X.N1() << X.N2() << N1() << N2() << i1 << i2;
        return false;
    }
}

bool DiskWriteMda::writeChunk(Mda32& X, bigint i1, bigint i2, bigint i3)
{
    if ((i3 == 0) && (X.N3() == 1) && (X.N3() == 1))
        return writeChunk(X, i1, i2);
    if ((X.N1() == N1()) && (X.N2() == N2()) && (i1 == 0) && (i2 == 0)) {
        return writeChunk(X, i1 + this->N1() * i2 + this->N1() * this->N2() * i3);
    }
    else {
        qWarning() << "This case not yet supported in 3d writeSubArray" << X.N1() << X.N2() << X.N3() << N1() << N2() << N3() << i1 << i2 << i3;
        return false;
    }
}

int DiskWriteMdaPrivate::determine_ndims(bigint N1, bigint N2, bigint N3, bigint N4, bigint N5, bigint N6)
{
#ifdef QT_CORE_LIB
    Q_UNUSED(N1)
    Q_UNUSED(N2)
#endif
    if (N6 > 1)
        return 6;
    if (N5 > 1)
        return 5;
    if (N4 > 1)
        return 4;
    if (N3 > 1)
        return 3;
    return 2;
}
