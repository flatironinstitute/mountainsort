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

#ifndef DISKWRITEMDA_H
#define DISKWRITEMDA_H

#include "mda.h"
#include "mda32.h"

#include <QString>
#include <mdaio.h>

class DiskWriteMdaPrivate;
class DiskWriteMda {
public:
    friend class DiskWriteMdaPrivate;
    DiskWriteMda();
    DiskWriteMda(int data_type, const QString& path, bigint N1, bigint N2, bigint N3 = 1, bigint N4 = 1, bigint N5 = 1, bigint N6 = 1);
    virtual ~DiskWriteMda();
    bool open(int data_type, const QString& path, bigint N1, bigint N2, bigint N3 = 1, bigint N4 = 1, bigint N5 = 1, bigint N6 = 1);
    bool open(const QString& path);
    void close();

    bigint N1();
    bigint N2();
    bigint N3();
    bigint N4();
    bigint N5();
    bigint N6();
    bigint totalSize();

    bool writeChunk(Mda& X, bigint i);
    bool writeChunk(Mda& X, bigint i1, bigint i2);
    bool writeChunk(Mda& X, bigint i1, bigint i2, bigint i3);

    bool writeChunk(Mda32& X, bigint i);
    bool writeChunk(Mda32& X, bigint i1, bigint i2);
    bool writeChunk(Mda32& X, bigint i1, bigint i2, bigint i3);

private:
    DiskWriteMdaPrivate* d;
};

#endif // DISKWRITEMDA_H
