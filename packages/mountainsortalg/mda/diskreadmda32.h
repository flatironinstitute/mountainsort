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

#ifndef DISKREADMDA32_H
#define DISKREADMDA32_H

#include "mda32.h"
#include "mdaio.h"

class DiskReadMda32Private;
/**
 * \class DiskReadMda32
 * @brief Read-only access to a .mda file, especially useful for huge arrays that cannot be practically loaded into memory.
 *
 * See also Mda32
 */
class DiskReadMda32 {
public:
    friend class DiskReadMda32Private;
    DiskReadMda32(const QString& path = ""); ///Constructor pointing to the .mda file specified by path (file name).
    DiskReadMda32(const DiskReadMda32& other); ///Copy constructor
    DiskReadMda32(const Mda32& X); ///Constructor based on an in-memory array. This enables passing an Mda32 into a function that expects a DiskReadMda32.
    //DiskReadMda32(const QJsonObject& prv_object);
    DiskReadMda32(int concat_dimension, const QList<DiskReadMda32>& arrays); //concatenation of arrays. For now concat_dimension must be 2
    DiskReadMda32(int concat_dimension, const QStringList& array_paths); //concatenation of arrays. For now concat_dimension must be 2
    virtual ~DiskReadMda32();
    void operator=(const DiskReadMda32& other);

    ///Set the path (file name) of the .mda file to read.
    void setPath(const QString& file_path);
    //void setPrvObject(const QJsonObject& prv_object);
    void setConcatPaths(int concat_dimension, const QStringList& paths);
    void setConcatDirectory(int concat_dimension, const QString& dir_path);

    QString makePath() const; //not capturing the reshaping
    //QJsonObject toPrvObject() const;

    ///The dimensions of the array
    bigint N1() const;
    bigint N2() const;
    bigint N3() const;
    bigint N4() const;
    bigint N5() const;
    bigint N6() const;
    bigint N(int dim) const; // dim is 1-based indexing
    bigint totalSize() const; //product of N1..N6

    MDAIO_HEADER mdaioHeader() const;

    bool reshape(bigint N1b, bigint N2b, bigint N3b = 1, bigint N4b = 1, bigint N5b = 1, bigint N6b = 1);
    DiskReadMda32 reshaped(bigint N1b, bigint N2b, bigint N3b = 1, bigint N4b = 1, bigint N5b = 1, bigint N6b = 1);

    ///Retrieve a chunk of the vectorized data of size 1xN starting at position i
    bool readChunk(Mda32& X, bigint i, bigint size) const;
    ///Retrieve a chunk of the vectorized data of size N1xN2 starting at position (i1,i2)
    bool readChunk(Mda32& X, bigint i1, bigint i2, bigint size1, bigint size2) const;
    ///Retrieve a chunk of the vectorized data of size N1xN2xN3 starting at position (i1,i2,i3)
    bool readChunk(Mda32& X, bigint i1, bigint i2, bigint i3, bigint size1, bigint size2, bigint size3) const;

    ///A slow method to retrieve the value at location i of the vectorized array for example value(3+4*N1())==value(3,4). Consider using readChunk() instead
    dtype32 value(bigint i) const;
    ///A slow method to retrieve the value at location (i1,i2) of the array. Consider using readChunk() instead
    dtype32 value(bigint i1, bigint i2) const;
    ///A slow method to retrieve the value at location (i1,i2,i3) of the array. Consider using readChunk() instead
    dtype32 value(bigint i1, bigint i2, bigint i3) const;

private:
    DiskReadMda32Private* d;
};

///Unit test
void diskreadmda_unit_test();

#endif // DISKREADMDA_H
