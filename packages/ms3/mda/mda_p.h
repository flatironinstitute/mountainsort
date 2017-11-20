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
#ifndef MDA_P_H
#define MDA_P_H

#include <QSharedData>
//#include "icounter.h"
//#include <objectregistry.h>
#include <cstring>
#include "textfile.h"

#define MDA_MAX_DIMS 6

template <typename T>
class MdaData : public QSharedData {
public:
    typedef T value_type;
    typedef T* pointer;
    typedef T& reference;

    MdaData()
        : QSharedData()
        , m_data(0)
        , m_dims(1, 1)
        , total_size(0)
    {
        /*
        ICounterManager* manager = ObjectRegistry::getObject<ICounterManager>();
        if (manager) {
            allocatedCounter = static_cast<IIntCounter*>(manager->counter("allocated_bytes"));
            freedCounter = static_cast<IIntCounter*>(manager->counter("freed_bytes"));
            bytesReadCounter = static_cast<IIntCounter*>(manager->counter("bytes_read"));
            bytesWrittenCounter = static_cast<IIntCounter*>(manager->counter("bytes_written"));
        }
        */
    }
    MdaData(const MdaData& other)
        : QSharedData(other)
        , m_data(0)
        , m_dims(other.m_dims)
        , total_size(other.total_size)
        /*
        , allocatedCounter(other.allocatedCounter)
        , freedCounter(other.freedCounter)
        , bytesReadCounter(other.bytesReadCounter)
        , bytesWrittenCounter(other.bytesWrittenCounter)
            */
    {
        allocate(total_size);
        std::copy(other.m_data, other.m_data + other.totalSize(), m_data);
    }
    ~MdaData()
    {
        deallocate();
    }
    bool allocate(T value, bigint N1, bigint N2, bigint N3 = 1, bigint N4 = 1, bigint N5 = 1, bigint N6 = 1)
    {
        deallocate();
        setDims(N1, N2, N3, N4, N5, N6);
        if (N1 > 0 && N2 > 0 && N3 > 0 && N4 > 0 && N5 > 0 && N6 > 0)
            setTotalSize(N1 * N2 * N3 * N4 * N5 * N6);
        else
            setTotalSize(0);

        if (totalSize() > 0) {
            allocate(totalSize());
            if (!constData()) {
                qCritical() << QString("Unable to allocate Mda of size %1x%2x%3x%4x%5x%6 (total=%7)").arg(N1).arg(N2).arg(N3).arg(N4).arg(N5).arg(N6).arg(totalSize());
                exit(-1);
            }
            if (value == 0.0) {
                std::memset(data(), 0, totalSize() * sizeof(value_type));
            }
            else
                std::fill(data(), data() + totalSize(), value);
        }
        return true;
    }

    inline bigint dim(bigint idx) const { return m_dims.at(idx); }
    inline bigint N1() const { return dim(0); }
    inline bigint N2() const { return dim(1); }

    void allocate(bigint size)
    {
        //m_data = (value_type*)::allocate(size * sizeof(value_type));
        m_data = (value_type*)malloc(size * sizeof(value_type));
        if (!m_data)
            return;
        //incrementBytesAllocatedCounter(totalSize() * sizeof(value_type));
    }
    void deallocate()
    {
        if (!m_data)
            return;
        free(m_data);
        //incrementBytesFreedCounter(totalSize() * sizeof(value_type));
        m_data = 0;
    }
    inline bigint totalSize() const { return total_size; }
    inline void setTotalSize(bigint ts) { total_size = ts; }
    inline T* data() { return m_data; }
    inline const T* constData() const { return m_data; }
    inline T at(bigint idx) const { return *(constData() + idx); }
    inline T at(bigint i1, bigint i2) const { return at(i1 + dim(0) * i2); }
    inline void set(T val, bigint idx) { m_data[idx] = val; }
    inline void set(T val, bigint i1, bigint i2) { set(val, i1 + dim(0) * i2); }

    inline bigint dims(bigint idx) const
    {
        if (idx < 0 || idx >= (bigint)m_dims.size())
            return 0;
        return *(m_dims.data() + idx);
    }
    void setDims(bigint n1, bigint n2, bigint n3, bigint n4, bigint n5, bigint n6)
    {
        m_dims.resize(MDA_MAX_DIMS);
        m_dims[0] = n1;
        m_dims[1] = n2;
        m_dims[2] = n3;
        m_dims[3] = n4;
        m_dims[4] = n5;
        m_dims[5] = n6;
    }

    int determine_num_dims(bigint N1, bigint N2, bigint N3, bigint N4, bigint N5, bigint N6) const
    {
        Q_UNUSED(N1);
        Q_UNUSED(N2);
        //changed this on 2/21/17 by jfm
        //if (!(N6 > 0 && N5 > 0 && N4 > 0 && N3 > 0 && N2 > 0 && N1 > 0))
        //    return 0;
        if (N6 != 1)
            return 6;
        if (N5 != 1)
            return 5;
        if (N4 != 1)
            return 4;
        if (N3 != 1)
            return 3;
        return 2;
    }
    bool safe_index(bigint i) const
    {
        return (i < totalSize());
    }
    bool safe_index(bigint i1, bigint i2) const
    {
        return (((bigint)i1 < dims(0)) && ((bigint)i2 < dims(1)));
    }
    bool safe_index(bigint i1, bigint i2, bigint i3) const
    {
        return (((bigint)i1 < dims(0)) && ((bigint)i2 < dims(1)) && ((bigint)i3 < dims(2)));
    }
    bool safe_index(bigint i1, bigint i2, bigint i3, bigint i4, bigint i5, bigint i6) const
    {
        return (
            (0 <= i1) && (i1 < dims(0))
            && (0 <= i2) && (i2 < dims(1))
            && (0 <= i3) && (i3 < dims(2))
            && (0 <= i4) && (i4 < dims(3))
            && (0 <= i5) && (i5 < dims(4))
            && (0 <= i6) && (i6 < dims(5)));
    }

    bool read_from_text_file(const QString& path)
    {
        QString txt = TextFile::read(path);
        if (txt.isEmpty()) {
            return false;
        }
        QStringList lines = txt.split("\n", QString::SkipEmptyParts);
        QStringList lines2;
        for (int i = 0; i < lines.count(); i++) {
            QString line = lines[i].trimmed();
            if (!line.isEmpty()) {
                if (i == 0) {
                    //check whether this is a header line, if so, don't include it
                    line = line.split(",", QString::SkipEmptyParts).join(" ");
                    QList<QString> vals = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
                    bool ok;
                    vals.value(0).toDouble(&ok);
                    if (ok) {
                        lines2 << line;
                    }
                }
                else {
                    lines2 << line;
                }
            }
        }
        for (int i = 0; i < lines2.count(); i++) {
            QString line = lines2[i].trimmed();
            line = line.split(",", QString::SkipEmptyParts).join(" ");
            QList<QString> vals = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
            if (i == 0) {
                allocate(0, vals.count(), lines2.count());
            }
            for (int j = 0; j < vals.count(); j++) {
                set(vals[j].toDouble(), j, i);
            }
        }
        return true;
    }
    bool write_to_text_file(const QString& path, char sep) const
    {
        /*
        char sep = ' ';
        if (path.endsWith(".csv"))
            sep = ',';
            */
        int max_num_entries = 1e6;
        if (N1() * N2() == max_num_entries) {
            qWarning() << "mda is too large to write text file";
            return false;
        }
        QList<QString> lines;
        for (int i = 0; i < N2(); i++) {
            QStringList vals;
            for (int j = 0; j < N1(); j++) {
                vals << QString("%1").arg(at(j, i));
            }
            QString line = vals.join(sep);
            lines << line;
        }
        return TextFile::write(path, lines.join("\n"));
    }

private:
    pointer m_data;
    std::vector<bigint> m_dims;
    bigint total_size;
    /*
    mutable IIntCounter* allocatedCounter = nullptr;
    mutable IIntCounter* freedCounter = nullptr;
    mutable IIntCounter* bytesReadCounter = nullptr;
    mutable IIntCounter* bytesWrittenCounter = nullptr;
    */
};

#endif // MDA_P_H
