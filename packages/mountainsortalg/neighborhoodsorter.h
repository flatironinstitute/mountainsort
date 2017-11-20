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
#ifndef NEIGHBORHOODSORTER_H
#define NEIGHBORHOODSORTER_H

#include "mda32.h"
#include "p_mountainsort3.h"

class NeighborhoodSorterPrivate;
class NeighborhoodSorter {
public:
    friend class NeighborhoodSorterPrivate;
    NeighborhoodSorter();
    virtual ~NeighborhoodSorter();

    void setOptions(P_mountainsort3_opts opts);
    void setMaxRAM(bigint max_ram_bytes);
    void addTimeChunk(bigint t, const Mda32& X, const QList<int>& channels, bigint padding_left, bigint padding_right);
    void sort(int num_threads);
    QVector<double> times() const;
    QVector<int> labels() const;
    Mda32 templates() const;

private:
    NeighborhoodSorterPrivate* d;
};

void setDiskBackedMdaTemporaryDirectory(QString path);

class DiskBackedMda32 {
public:
    DiskBackedMda32();
    DiskBackedMda32(const Mda32& X);
    virtual ~DiskBackedMda32();
    void store(const Mda32& X);
    void retrieve(Mda32& X) const;
    void remove();

private:
    QString m_tmp_path;
};

#endif // NEIGHBORHOODSORTER_H
