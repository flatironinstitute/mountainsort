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
#ifndef MDAREADER_P_H
#define MDAREADER_P_H

#include <QByteArray>
#include <QtEndian>

class QIODevice;
class Mda;
class Mda32;

class MdaIOHandler {
public:
    MdaIOHandler();
    virtual ~MdaIOHandler() {}
    virtual bool canRead() const = 0;
    virtual bool canWrite() const;
    virtual bool read(Mda*) = 0;
    virtual bool read(Mda32*);
    virtual bool write(const Mda&);
    virtual bool write(const Mda32&);
    QIODevice* device() const;
    void setDevice(QIODevice* d);
    QByteArray format() const;
    void setFormat(const QByteArray& ba);

private:
    QIODevice* dev;
    QByteArray fmt;
};

class MdaIOHandlerFactory {
public:
    virtual ~MdaIOHandlerFactory() {}
    virtual MdaIOHandler* create(QIODevice* device, const QByteArray& format = QByteArray()) const = 0;
};

class MdaIOHandlerMDA : public MdaIOHandler {
public:
    enum DataType {
        Complex = -1,
        Byte = -2,
        Float32 = -3,
        Int16 = -4,
        Int32 = -5,
        UInt16 = -6,
        Float64 = -7,
        UInt32 = -8
    };

    const int maxDims = 50;
    MdaIOHandlerMDA(QIODevice* d, const QByteArray& fmt);
    bool canRead() const;
    bool canWrite() const;
    bool read(Mda* mda);
    bool read(Mda32* mda);
    bool write(const Mda& mda);
    bool write(const Mda32& mda);

    struct Header {
        int32_t data_type;
        int32_t num_bytes_per_entry;
        int32_t num_dims;
        QVector<uint64_t> dims;
    };
    template <typename T>
    bool readLE(T* ptr)
    {
        if (device()->read((char*)ptr, sizeof(T)) != sizeof(T))
            return false;
        *ptr = qFromLittleEndian(*ptr);
        return true;
    }
    template <typename T>
    bool writeLE(T val)
    {
        val = qToLittleEndian(val);
        return (device()->write((char*)&val, sizeof(T)) == sizeof(T));
    }

    template <typename src, typename dst>
    bool readData(dst* ptr, size_t cnt = 1)
    {
        src val;
        for (size_t i = 0; i < cnt; ++i) {
            if (!readLE(&val))
                return false;
            *(ptr++) = val;
        }
        return true;
    }
    template <typename dst, typename src>
    bool writeData(src* ptr, size_t cnt = 1)
    {
        dst val;
        for (size_t i = 0; i < cnt; ++i) {
            val = *(ptr++);
            if (!writeLE(val))
                return false;
        }
        return true;
    }
};

class MdaIOHandlerMDAFactory : public MdaIOHandlerFactory {
public:
    MdaIOHandlerMDAFactory() {}
    MdaIOHandler* create(QIODevice* device, const QByteArray& format) const;
};

class MdaIOHandlerCSV : public MdaIOHandler {
public:
    MdaIOHandlerCSV(QIODevice* device, const QByteArray& format);
    bool canRead() const;
    bool canWrite() const;
    bool read(Mda* mda);
    bool read(Mda32* mda);
    bool write(const Mda& mda);
    bool write(const Mda32& mda);
};

class MdaIOHandlerCSVFactory : public MdaIOHandlerFactory {
public:
    MdaIOHandlerCSVFactory() {}
    MdaIOHandler* create(QIODevice* device, const QByteArray& format) const;
};

#endif // MDAREADER_P_H
