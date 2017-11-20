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
#include "mda32.h"
#include "mda_p.h"
#include "mdaio.h"
//#include <cachemanager.h>
#include <stdio.h>
//#include "taskprogress.h"
//#include "icounter.h"
//#include <objectregistry.h>
#include <QFile>

#define MDA_MAX_DIMS 6

class MdaDataFloat : public MdaData<float> {
};

Mda32::Mda32(bigint N1, bigint N2, bigint N3, bigint N4, bigint N5, bigint N6)
{
    d = new MdaDataFloat;
    this->allocate(N1, N2, N3, N4, N5, N6);
}

Mda32::Mda32(const QString mda_filename)
{
    d = new MdaDataFloat;
    this->read(mda_filename);
}

Mda32::Mda32(const Mda32& other)
{
    d = other.d;
}

void Mda32::operator=(const Mda32& other)
{
    d = other.d;
}

Mda32::~Mda32()
{
}

bool Mda32::allocate(bigint N1, bigint N2, bigint N3, bigint N4, bigint N5, bigint N6)
{

    return d->allocate((float)0, N1, N2, N3, N4, N5, N6);
}

bool Mda32::read(const QString& path)
{
    return read(path.toLatin1().data());
}

bool Mda32::write8(const QString& path) const
{
    return write8(path.toLatin1().data());
}

bool Mda32::write32(const QString& path) const
{
    return write32(path.toLatin1().data());
}

bool Mda32::write64(const QString& path) const
{
    return write64(path.toLatin1().data());
}

bool Mda32::writeCsv(const QString& path) const
{
    return d->write_to_text_file(path, ',');
}

bool Mda32::read(const char* path)
{
    if ((QString(path).endsWith(".txt")) || (QString(path).endsWith(".csv"))) {
        return d->read_from_text_file(path);
    }
    FILE* input_file = fopen(path, "rb");
    if (!input_file) {
        printf("Warning: Unable to open mda file for reading: %s\n", path);
        return false;
    }
    MDAIO_HEADER H;
    if (!mda_read_header(&H, input_file)) {
        qWarning() << "Problem reading mda file: " + QString(path);
        fclose(input_file);
        return false;
    }
    this->allocate(H.dims[0], H.dims[1], H.dims[2], H.dims[3], H.dims[4], H.dims[5]);
    mda_read_float32(d->data(), &H, d->totalSize(), input_file);
    //d->incrementBytesReadCounter(d->totalSize() * H.num_bytes_per_entry);
    fclose(input_file);
    return true;
}

bool Mda32::readCsv(const QString& path)
{
    return d->read_from_text_file(path);
}

bool Mda32::write8(const char* path) const
{
    if (QString(path).endsWith(".txt")) {
        return d->write_to_text_file(path, ' ');
    }
    if (QString(path).endsWith(".csv")) {
        return d->write_to_text_file(path, ',');
    }
    FILE* output_file = fopen(path, "wb");
    if (!output_file) {
        printf("Warning: Unable to open mda file for writing: %s\n", path);
        return false;
    }
    MDAIO_HEADER H;
    H.data_type = MDAIO_TYPE_BYTE;
    H.num_bytes_per_entry = 1;
    for (int i = 0; i < MDAIO_MAX_DIMS; i++)
        H.dims[i] = 1;
    for (int i = 0; i < MDA_MAX_DIMS; i++)
        H.dims[i] = d->dims(i);
    H.num_dims = d->determine_num_dims(N1(), N2(), N3(), N4(), N5(), N6());
    mda_write_header(&H, output_file);
    mda_write_float32(d->constData(), &H, d->totalSize(), output_file);
    //d->incrementBytesWrittenCounter(d->totalSize() * H.num_bytes_per_entry);
    fclose(output_file);
    return true;
}

bool Mda32::write32(const char* path) const
{
    if (QString(path).endsWith(".txt")) {
        return d->write_to_text_file(path, ' ');
    }
    if (QString(path).endsWith(".csv")) {
        return d->write_to_text_file(path, ',');
    }
    FILE* output_file = fopen(path, "wb");
    if (!output_file) {
        printf("Warning: Unable to open mda file for writing: %s\n", path);
        return false;
    }
    MDAIO_HEADER H;
    H.data_type = MDAIO_TYPE_FLOAT32;
    H.num_bytes_per_entry = 4;
    for (int i = 0; i < MDAIO_MAX_DIMS; i++)
        H.dims[i] = 1;
    for (int i = 0; i < MDA_MAX_DIMS; i++)
        H.dims[i] = d->dims(i);
    H.num_dims = d->determine_num_dims(N1(), N2(), N3(), N4(), N5(), N6());
    mda_write_header(&H, output_file);
    mda_write_float32(d->constData(), &H, d->totalSize(), output_file);
    //d->incrementBytesWrittenCounter(d->totalSize() * H.num_bytes_per_entry);
    fclose(output_file);
    return true;
}

bool Mda32::write64(const char* path) const
{
    if (QString(path).endsWith(".txt")) {
        return d->write_to_text_file(path, ' ');
    }
    if (QString(path).endsWith(".csv")) {
        return d->write_to_text_file(path, ',');
    }
    FILE* output_file = fopen(path, "wb");
    if (!output_file) {
        printf("Warning: Unable to open mda file for writing: %s\n", path);
        return false;
    }
    MDAIO_HEADER H;
    H.data_type = MDAIO_TYPE_FLOAT64;
    H.num_bytes_per_entry = 4;
    for (int i = 0; i < MDAIO_MAX_DIMS; i++)
        H.dims[i] = 1;
    for (int i = 0; i < MDA_MAX_DIMS; i++)
        H.dims[i] = d->dims(i);
    H.num_dims = d->determine_num_dims(N1(), N2(), N3(), N4(), N5(), N6());
    mda_write_header(&H, output_file);
    mda_write_float32(d->constData(), &H, d->totalSize(), output_file);
    //d->incrementBytesWrittenCounter(d->totalSize() * H.num_bytes_per_entry);
    fclose(output_file);
    return true;
}

int Mda32::ndims() const
{
    return d->determine_num_dims(N1(), N2(), N3(), N4(), N5(), N6());
}

bigint Mda32::N1() const
{
    return d->dims(0);
}

bigint Mda32::N2() const
{
    return d->dims(1);
}

bigint Mda32::N3() const
{
    return d->dims(2);
}

bigint Mda32::N4() const
{
    return d->dims(3);
}

bigint Mda32::N5() const
{
    return d->dims(4);
}

bigint Mda32::N6() const
{
    return d->dims(5);
}

bigint Mda32::totalSize() const
{
    return d->totalSize();
}

bigint Mda32::size(int dimension_index) const
{
    if (dimension_index < 0)
        return 0;
    if (dimension_index >= MDA_MAX_DIMS)
        return 1;
    return d->dims(dimension_index);
}

dtype32 Mda32::get(bigint i) const
{
    return d->at(i);
}

dtype32 Mda32::get(bigint i1, bigint i2) const
{
    return d->at(i1 + d->dims(0) * i2);
}

dtype32 Mda32::get(bigint i1, bigint i2, bigint i3) const
{
    return d->at(i1 + d->dims(0) * i2 + d->dims(0) * d->dims(1) * i3);
}

dtype32 Mda32::get(bigint i1, bigint i2, bigint i3, bigint i4, bigint i5, bigint i6) const
{
    const bigint d01 = d->dims(0) * d->dims(1);
    const bigint d02 = d01 * d->dims(2);
    const bigint d03 = d02 * d->dims(3);
    const bigint d04 = d03 * d->dims(4);

    return d->at(i1 + d->dims(0) * i2 + d01 * i3 + d02 * i4 + d03 * i5 + d04 * i6);
}

dtype32 Mda32::value(bigint i) const
{
    if (!d->safe_index(i))
        return 0;
    return get(i);
}

dtype32 Mda32::value(bigint i1, bigint i2) const
{
    if (!d->safe_index(i1, i2))
        return 0;
    return get(i1, i2);
}

dtype32 Mda32::value(bigint i1, bigint i2, bigint i3) const
{
    if (!d->safe_index(i1, i2, i3))
        return 0;
    return get(i1, i2, i3);
}

dtype32 Mda32::value(bigint i1, bigint i2, bigint i3, bigint i4, bigint i5, bigint i6) const
{
    if (!d->safe_index(i1, i2, i3, i4, i5, i6))
        return 0;
    return get(i1, i2, i3, i4, i5, i6);
}

void Mda32::setValue(dtype32 val, bigint i)
{
    if (!d->safe_index(i))
        return;
    set(val, i);
}

void Mda32::setValue(dtype32 val, bigint i1, bigint i2)
{
    if (!d->safe_index(i1, i2))
        return;
    set(val, i1, i2);
}

void Mda32::setValue(dtype32 val, bigint i1, bigint i2, bigint i3)
{
    if (!d->safe_index(i1, i2, i3))
        return;
    set(val, i1, i2, i3);
}

void Mda32::setValue(dtype32 val, bigint i1, bigint i2, bigint i3, bigint i4, bigint i5, bigint i6)
{
    if (!d->safe_index(i1, i2, i3, i4, i5, i6))
        return;
    set(val, i1, i2, i3, i4, i5, i6);
}

dtype32* Mda32::dataPtr()
{
    return d->data();
}

const dtype32* Mda32::constDataPtr() const
{
    return d->constData();
}

dtype32* Mda32::dataPtr(bigint i)
{
    return d->data() + i;
}

dtype32* Mda32::dataPtr(bigint i1, bigint i2)
{
    return d->data() + (i1 + N1() * i2);
}

dtype32* Mda32::dataPtr(bigint i1, bigint i2, bigint i3)
{
    return d->data() + (i1 + N1() * i2 + N1() * N2() * i3);
}

dtype32* Mda32::dataPtr(bigint i1, bigint i2, bigint i3, bigint i4, bigint i5, bigint i6)
{
    const bigint N12 = N1() * N2();
    const bigint N13 = N12 * N3();
    const bigint N14 = N13 * N4();
    const bigint N15 = N14 * N5();

    return d->data() + (i1 + N1() * i2 + N12 * i3 + N13 * i4 + N14 * i5 + N15 * i6);
}

void Mda32::getChunk(Mda32& ret, bigint i, bigint size) const
{
    // A lot of bugs fixed on 5/31/16
    bigint a_begin = i;
    bigint x_begin = 0;
    bigint a_end = i + size - 1;
    //    bigint x_end = size - 1;  // unused?

    if (a_begin < 0) {
        x_begin += 0 - a_begin;
        a_begin += 0 - a_begin;
    }
    if (a_end >= (bigint)d->totalSize()) {
        //        x_end += (bigint)d->totalSize() - 1 - a_end; // unused?
        a_end += (bigint)d->totalSize() - 1 - a_end;
    }

    ret.allocate(1, size);

    const float* ptr1 = this->constDataPtr();
    float* ptr2 = ret.dataPtr();

    std::copy(ptr1 + a_begin, ptr1 + a_end + 1, ptr2 + x_begin);
}

void Mda32::getChunk(Mda32& ret, bigint i1, bigint i2, bigint size1, bigint size2) const
{
    // A lot of bugs fixed on 5/31/16
    bigint a1_begin = i1;
    bigint x1_begin = 0;
    bigint a1_end = i1 + size1 - 1;
    bigint x1_end = size1 - 1;
    if (a1_begin < 0) {
        x1_begin += 0 - a1_begin;
        a1_begin += 0 - a1_begin;
    }
    if (a1_end >= N1()) {
        x1_end += N1() - 1 - a1_end;
        a1_end += N1() - 1 - a1_end;
    }

    bigint a2_begin = i2;
    bigint x2_begin = 0;
    bigint a2_end = i2 + size2 - 1;
    bigint x2_end = size2 - 1;
    if (a2_begin < 0) {
        x2_begin += 0 - a2_begin;
        a2_begin += 0 - a2_begin;
    }
    if (a2_end >= N2()) {
        x2_end += N2() - 1 - a2_end;
        a2_end += N2() - 1 - a2_end;
    }

    ret.allocate(size1, size2);

    const float* ptr1 = this->constDataPtr();
    float* ptr2 = ret.dataPtr();

    for (bigint ind2 = 0; ind2 <= a2_end - a2_begin; ind2++) {
        bigint ii_out = (ind2 + x2_begin) * size1 + x1_begin; //bug fixed on 5/31/16 by jfm
        bigint ii_in = (ind2 + a2_begin) * N1() + a1_begin; //bug fixed on 5/31/16 by jfm
        for (bigint ind1 = 0; ind1 <= a1_end - a1_begin; ind1++) {
            ptr2[ii_out] = ptr1[ii_in];
            ii_in++;
            ii_out++;
        }
    }
}

void Mda32::getChunk(Mda32& ret, bigint i1, bigint i2, bigint i3, bigint size1, bigint size2, bigint size3) const
{
    // A lot of bugs fixed on 5/31/16
    bigint a1_begin = i1;
    bigint x1_begin = 0;
    bigint a1_end = i1 + size1 - 1;
    bigint x1_end = size1 - 1;
    if (a1_begin < 0) {
        x1_begin += 0 - a1_begin;
        a1_begin += 0 - a1_begin;
    }
    if (a1_end >= N1()) {
        x1_end += N1() - 1 - a1_end;
        a1_end += N1() - 1 - a1_end;
    }

    bigint a2_begin = i2;
    bigint x2_begin = 0;
    bigint a2_end = i2 + size2 - 1;
    bigint x2_end = size2 - 1;
    if (a2_begin < 0) {
        x2_begin += 0 - a2_begin;
        a2_begin += 0 - a2_begin;
    }
    if (a2_end >= N2()) {
        x2_end += N2() - 1 - a2_end;
        a2_end += N2() - 1 - a2_end;
    }

    bigint a3_begin = i3;
    bigint x3_begin = 0;
    bigint a3_end = i3 + size3 - 1;
    bigint x3_end = size3 - 1;
    if (a3_begin < 0) {
        x3_begin += 0 - a3_begin;
        a3_begin += 0 - a3_begin;
    }
    if (a3_end >= N3()) {
        x3_end += N3() - 1 - a3_end;
        a3_end += N3() - 1 - a3_end;
    }

    ret.allocate(size1, size2, size3);

    const float* ptr1 = this->constDataPtr();
    float* ptr2 = ret.dataPtr();

    for (bigint ind3 = 0; ind3 <= a3_end - a3_begin; ind3++) {
        for (bigint ind2 = 0; ind2 <= a2_end - a2_begin; ind2++) {
            bigint ii_out = x1_begin + (ind2 + x2_begin) * size1 + (ind3 + x3_begin) * size1 * size2; //bug fixed on 5/31/16 by jfm
            bigint ii_in = a1_begin + (ind2 + a2_begin) * N1() + (ind3 + a3_begin) * N1() * N2(); //bug fixed on 5/31/16 by jfm
            for (bigint ind1 = 0; ind1 <= a1_end - a1_begin; ind1++) {
                ptr2[ii_out] = ptr1[ii_in];
                ii_in++;
                ii_out++;
            }
        }
    }
}

void Mda32::setChunk(Mda32& X, bigint i)
{
    bigint size = X.totalSize();

    bigint a_begin = i;
    bigint x_begin = 0;
    bigint a_end = i + size - 1;
    bigint x_end = size - 1;

    if (a_begin < 0) {
        a_begin += 0 - a_begin;
        x_begin += 0 - a_begin;
    }
    if (a_end >= (bigint)d->totalSize()) {
        a_end += (bigint)d->totalSize() - 1 - a_end;
        x_end += (bigint)d->totalSize() - 1 - a_end;
    }

    float* ptr1 = this->dataPtr();
    float* ptr2 = X.dataPtr();

    bigint ii = 0;
    for (bigint a = a_begin; a <= a_end; a++) {
        ptr1[a_begin + ii] = ptr2[x_begin + ii];
        ii++;
    }
}

void Mda32::setChunk(Mda32& X, bigint i1, bigint i2)
{
    bigint size1 = X.N1();
    bigint size2 = X.N2();

    bigint a1_begin = i1;
    bigint x1_begin = 0;
    bigint a1_end = i1 + size1 - 1;
    bigint x1_end = size1 - 1;
    if (a1_begin < 0) {
        a1_begin += 0 - a1_begin;
        x1_begin += 0 - a1_begin;
    }
    if (a1_end >= N1()) {
        a1_end += N1() - 1 - a1_end;
        x1_end += N1() - 1 - a1_end;
    }

    bigint a2_begin = i2;
    bigint x2_begin = 0;
    bigint a2_end = i2 + size2 - 1;
    bigint x2_end = size2 - 1;
    if (a2_begin < 0) {
        a2_begin += 0 - a2_begin;
        x2_begin += 0 - a2_begin;
    }
    if (a2_end >= N2()) {
        a2_end += N2() - 1 - a2_end;
        x2_end += N2() - 1 - a2_end;
    }

    dtype32* ptr1 = this->dataPtr();
    dtype32* ptr2 = X.dataPtr();

    for (bigint ind2 = 0; ind2 <= a2_end - a2_begin; ind2++) {
        bigint ii_out = (ind2 + x2_begin) * size1;
        bigint ii_in = (ind2 + a2_begin) * N1();
        for (bigint ind1 = 0; ind1 <= a1_end - a1_begin; ind1++) {
            ptr1[ii_in] = ptr2[ii_out];
            ii_in++;
            ii_out++;
        }
    }
}

void Mda32::setChunk(Mda32& X, bigint i1, bigint i2, bigint i3)
{
    bigint size1 = X.N1();
    bigint size2 = X.N2();
    bigint size3 = X.N3();

    bigint a1_begin = i1;
    bigint x1_begin = 0;
    bigint a1_end = i1 + size1 - 1;
    bigint x1_end = size1 - 1;
    if (a1_begin < 0) {
        a1_begin += 0 - a1_begin;
        x1_begin += 0 - a1_begin;
    }
    if (a1_end >= N1()) {
        a1_end += N1() - 1 - a1_end;
        x1_end += N1() - 1 - a1_end;
    }

    bigint a2_begin = i2;
    bigint x2_begin = 0;
    bigint a2_end = i2 + size2 - 1;
    bigint x2_end = size2 - 1;
    if (a2_begin < 0) {
        a2_begin += 0 - a2_begin;
        x2_begin += 0 - a2_begin;
    }
    if (a2_end >= N2()) {
        a2_end += N2() - 1 - a2_end;
        x2_end += N2() - 1 - a2_end;
    }

    bigint a3_begin = i3;
    bigint x3_begin = 0;
    bigint a3_end = i3 + size3 - 1;
    bigint x3_end = size3 - 1;
    if (a3_begin < 0) {
        a2_begin += 0 - a3_begin;
        x3_begin += 0 - a3_begin;
    }
    if (a3_end >= N3()) {
        a3_end += N3() - 1 - a3_end;
        x3_end += N3() - 1 - a3_end;
    }

    dtype32* ptr1 = this->dataPtr();
    dtype32* ptr2 = X.dataPtr();

    for (bigint ind3 = 0; ind3 <= a3_end - a3_begin; ind3++) {
        for (bigint ind2 = 0; ind2 <= a2_end - a2_begin; ind2++) {
            bigint ii_out = (ind2 + x2_begin) * size1 + (ind3 + x3_begin) * size1 * size2;
            bigint ii_in = (ind2 + a2_begin) * N1() + (ind3 + a3_begin) * N1() * N2();
            for (bigint ind1 = 0; ind1 <= a1_end - a1_begin; ind1++) {
                ptr1[ii_in] = ptr2[ii_out];
                ii_in++;
                ii_out++;
            }
        }
    }
}

dtype32 Mda32::minimum() const
{
    bigint NN = this->totalSize();
    const dtype32* ptr = this->constDataPtr();
    if ((!NN) || (!ptr)) {
        return 0;
    }
    return *std::min_element(ptr, ptr + NN);
}

dtype32 Mda32::maximum() const
{
    bigint NN = this->totalSize();
    const dtype32* ptr = this->constDataPtr();
    if ((!NN) || (!ptr)) {
        return 0;
    }
    return *std::max_element(ptr, ptr + NN);
}

bool Mda32::reshape(bigint N1b, bigint N2b, bigint N3b, bigint N4b, bigint N5b, bigint N6b)
{
    if (N1b * N2b * N3b * N4b * N5b * N6b != this->totalSize()) {
        qWarning() << "Unable to reshape Mda32, wrong total size";
        qWarning() << N1b << N2b << N3b << N4b << N5b << N6b;
        qWarning() << N1() << N2() << N3() << N4() << N5() << N6();
        return false;
    }
    d->setDims(N1b, N2b, N3b, N4b, N5b, N6b);
    return true;
}

void Mda32::set(dtype32 val, bigint i)
{
    d->set(val, i);
}

void Mda32::set(dtype32 val, bigint i1, bigint i2)
{
    d->set(val, i1 + d->dims(0) * i2);
}

void Mda32::set(dtype32 val, bigint i1, bigint i2, bigint i3)
{
    d->set(val, i1 + d->dims(0) * i2 + d->dims(0) * d->dims(1) * i3);
}

void Mda32::set(dtype32 val, bigint i1, bigint i2, bigint i3, bigint i4, bigint i5, bigint i6)
{
    const bigint d01 = d->dims(0) * d->dims(1);
    const bigint d02 = d01 * d->dims(2);
    const bigint d03 = d02 * d->dims(3);
    const bigint d04 = d03 * d->dims(4);

    d->set(val, i1 + d->dims(0) * i2 + d01 * i3 + d02 * i4 + d03 * i5 + d04 * i6);
}
