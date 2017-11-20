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
#include "mda/mdareader.h"
#include "mda/mdareader_p.h"
#include <QSaveFile>
#include <QFile>
#include <QFileInfo>
#include <QSharedPointer>

Q_GLOBAL_STATIC(QList<QSharedPointer<MdaIOHandlerFactory> >, _mdaReaderFactories)

class MdaReaderPrivate {
public:
    MdaReaderPrivate(MdaReader* qq, QIODevice* dev = 0)
        : q(qq)
        , device(dev)
        , ownsDevice(false)
    {
        if (_mdaReaderFactories->isEmpty()) {
            _mdaReaderFactories->append(QSharedPointer<MdaIOHandlerFactory>(new MdaIOHandlerMDAFactory));
            _mdaReaderFactories->append(QSharedPointer<MdaIOHandlerFactory>(new MdaIOHandlerCSVFactory));
        }
    }

    MdaReader* q;
    QIODevice* device;
    QByteArray format;
    bool ownsDevice;
};

MdaReader::MdaReader()
    : d(new MdaReaderPrivate(this))
{
}

MdaReader::MdaReader(QIODevice* dev, const QByteArray& format)
    : d(new MdaReaderPrivate(this, dev))
{
    d->format = format;
}

MdaReader::MdaReader(const QString& fileName, const QByteArray& format)
    : d(new MdaReaderPrivate(this))
{
    d->format = format;
    d->device = new QFile(fileName);
    d->ownsDevice = true;
}

MdaReader::~MdaReader()
{
    if (d->ownsDevice)
        d->device->deleteLater();
    delete d;
}

bool MdaReader::canRead() const
{
    for (int i = 0; i < _mdaReaderFactories->size(); ++i) {
        MdaIOHandler* handler = _mdaReaderFactories->at(i)->create(device(), format());
        if (handler->canRead()) {
            delete handler;
            return true;
        }
        delete handler;
    }
    return false;
}

QIODevice* MdaReader::device() const
{
    return d->device;
}

QString MdaReader::fileName() const
{
    if (QFile* f = qobject_cast<QFile*>(d->device)) {
        return f->fileName();
    }
    return QString();
}

QByteArray MdaReader::format() const
{
    return d->format;
}

Mda MdaReader::read()
{
    Mda mda;
    if (read(&mda)) {
        return std::move(mda);
    }
    return Mda();
}

Mda32 MdaReader::read32()
{
    Mda32 mda;
    if (read(&mda)) {
        return std::move(mda);
    }
    return Mda32();
}

bool MdaReader::read(Mda* mda)
{
    if (!device()->isOpen() && !device()->open(QIODevice::ReadOnly))
        return false;
    for (int i = 0; i < _mdaReaderFactories->size(); ++i) {
        MdaIOHandler* handler = _mdaReaderFactories->at(i)->create(device(), format());
        if (handler->canRead()) {
            bool result = handler->read(mda);
            delete handler;
            return result;
        }
        delete handler;
    }
    return false;
}

bool MdaReader::read(Mda32* mda)
{
    if (!device()->isOpen() && !device()->open(QIODevice::ReadOnly))
        return false;
    MdaIOHandlerMDA h(device(), format());
    if (!h.canRead())
        return false;
    return h.read(mda);
}

void MdaReader::setDevice(QIODevice* dev)
{
    if (d->device && d->ownsDevice) {
        d->device->deleteLater();
        d->device = 0;
        d->ownsDevice = false;
    }
    d->device = dev;
}

void MdaReader::setFileName(const QString& fileName)
{
    if (d->device && d->ownsDevice) {
        d->device->deleteLater();
        d->device = 0;
    }
    d->device = new QFile(fileName);
    d->ownsDevice = true;
}

class MdaWriterPrivate {
public:
    MdaWriterPrivate(MdaWriter* qq, QIODevice* dev = 0)
        : q(qq)
        , device(dev)
        , ownsDevice(false)
    {
    }

    MdaWriter* q;
    QIODevice* device;
    QByteArray format;
    bool ownsDevice;
};

MdaWriter::MdaWriter()
    : d(new MdaWriterPrivate(this))
{
}

MdaWriter::MdaWriter(QIODevice* dev, const QByteArray& format)
    : d(new MdaWriterPrivate(this, dev))
{
    d->format = format;
}

MdaWriter::MdaWriter(const QString& fileName, const QByteArray& format)
    : d(new MdaWriterPrivate(this))
{
    d->device = new QSaveFile(fileName);
    d->ownsDevice = true;
    d->format = format;
}

MdaWriter::~MdaWriter()
{
    if (d->ownsDevice) {
        if (d->device)
            d->device->deleteLater();
    }
    delete d;
}

bool MdaWriter::canWrite() const
{
    MdaIOHandlerMDA h(device(), format());
    return h.canWrite();
}

QIODevice* MdaWriter::device() const
{
    return d->device;
}

QString MdaWriter::fileName() const
{
    if (QFileDevice* f = qobject_cast<QFileDevice*>(d->device))
        return f->fileName();

    return QString();
}

QByteArray MdaWriter::format() const
{
    return d->format;
}

void MdaWriter::setDevice(QIODevice* dev)
{
    if (d->device && d->ownsDevice) {
        d->device->deleteLater();
        d->device = 0;
        d->ownsDevice = false;
    }
    d->device = dev;
}

void MdaWriter::setFileName(const QString& fileName)
{
    if (d->device && d->ownsDevice) {
        d->device->deleteLater();
        d->device = 0;
    }
    d->device = new QSaveFile(fileName);
    d->ownsDevice = true;
}

bool MdaWriter::write(const Mda& mda)
{
    bool closeAfterWrite = !device()->isOpen();
    if (!device()->isOpen() && !device()->open(QIODevice::WriteOnly | QIODevice::Truncate))
        return false;
    for (int i = 0; i < _mdaReaderFactories->size(); ++i) {
        MdaIOHandler* handler = _mdaReaderFactories->at(i)->create(device(), format());
        if (handler->canWrite()) {
            bool result = handler->write(mda);
            delete handler;
            if (QSaveFile* f = qobject_cast<QSaveFile*>(device()))
                f->commit();
            else if (closeAfterWrite)
                device()->close();
            return result;
        }
        delete handler;
    }
    return false;
}

bool MdaWriter::write(const Mda32& mda)
{
    bool closeAfterWrite = !device()->isOpen();
    if (!device()->isOpen() && !device()->open(QIODevice::WriteOnly | QIODevice::Truncate))
        return false;
    for (int i = 0; i < _mdaReaderFactories->size(); ++i) {
        MdaIOHandler* handler = _mdaReaderFactories->at(i)->create(device(), format());
        if (handler->canWrite()) {
            bool result = handler->write(mda);
            delete handler;
            if (QSaveFile* f = qobject_cast<QSaveFile*>(device()))
                f->commit();
            else if (closeAfterWrite)
                device()->close();
            return result;
        }
        delete handler;
    }
    return false;
}

MdaIOHandler::MdaIOHandler()
    : dev(0)
{
}

bool MdaIOHandler::canWrite() const { return false; }

bool MdaIOHandler::read(Mda32*) { return false; }

bool MdaIOHandler::write(const Mda&) { return false; }

bool MdaIOHandler::write(const Mda32&) { return false; }

QIODevice* MdaIOHandler::device() const { return dev; }

void MdaIOHandler::setDevice(QIODevice* d) { dev = d; }

QByteArray MdaIOHandler::format() const { return fmt; }

void MdaIOHandler::setFormat(const QByteArray& ba) { fmt = ba; }

MdaIOHandlerMDA::MdaIOHandlerMDA(QIODevice* d, const QByteArray& fmt)
{
    setFormat(fmt);
    setDevice(d);
}

bool MdaIOHandlerMDA::canRead() const
{
    if (!device())
        return false;
    if (format().toLower() == "mda")
        return true;
    if (format().toLower().startsWith("mda."))
        return true;
    if (format().isEmpty()) {
        // try to auto detect by file name
        if (QFileDevice* f = qobject_cast<QFileDevice*>(device())) {
            QFileInfo finfo(f->fileName());
            if (finfo.suffix().toLower() == "mda")
                return true;
        }
    }
    return false;
}

bool MdaIOHandlerMDA::canWrite() const
{
    if (!device())
        return false;
    if (format().toLower() == "mda")
        return true;
    if (format().toLower().startsWith("mda."))
        return true;
    return false;
}

bool MdaIOHandlerMDA::read(Mda* mda)
{
    bool use64bitdims = false;
    Header header;
    if (!readLE(&header.data_type))
        return false;
    if (!readLE(&header.num_bytes_per_entry))
        return false;
    if (!readLE(&header.num_dims))
        return false;
    if (header.num_dims == 0 || header.num_dims > maxDims)
        return false;
    if (header.num_dims < 0) {
        use64bitdims = true;
        header.num_dims = -header.num_dims;
    }
    header.dims.resize(header.num_dims);
    if (use64bitdims) {
        for (size_t i = 0; i < (size_t)header.num_dims; ++i) {
            if (!readLE(header.dims.data() + i))
                return false;
        }
    }
    else {
        int32_t data;
        for (size_t i = 0; i < (size_t)header.num_dims; ++i) {
            if (!readLE(&data))
                return false;
            header.dims[i] = data;
        }
    }
    while (header.dims.size() < 6)
        header.dims.append(1);
    if (!mda->allocate(header.dims[0], header.dims[1], header.dims[2], header.dims[3], header.dims[4], header.dims[5]))
        return false;
    switch (header.data_type) {
    case Byte:
        return readData<unsigned char>(mda->dataPtr(), mda->totalSize());
    case Float32:
        return readData<float>(mda->dataPtr(), mda->totalSize());
    case Int16:
        return readData<int16_t>(mda->dataPtr(), mda->totalSize());
    case Int32:
        return readData<int32_t>(mda->dataPtr(), mda->totalSize());
    case UInt16:
        return readData<uint16_t>(mda->dataPtr(), mda->totalSize());
    case Float64:
        return readData<double>(mda->dataPtr(), mda->totalSize());
    case UInt32:
        return readData<uint32_t>(mda->dataPtr(), mda->totalSize());
    default:
        return false;
    }
    return true;
}

bool MdaIOHandlerMDA::read(Mda32* mda)
{
    bool use64bitdims = false;
    Header header;
    if (!readLE(&header.data_type))
        return false;
    if (!readLE(&header.num_bytes_per_entry))
        return false;
    if (!readLE(&header.num_dims))
        return false;
    if (header.num_dims == 0 || header.num_dims > maxDims)
        return false;
    if (header.num_dims < 0) {
        use64bitdims = true;
        header.num_dims = -header.num_dims;
    }
    header.dims.resize(header.num_dims);
    if (use64bitdims) {
        for (size_t i = 0; i < (size_t)header.num_dims; ++i) {
            if (!readLE(header.dims.data() + i))
                return false;
        }
    }
    else {
        int32_t data;
        for (size_t i = 0; i < (size_t)header.num_dims; ++i) {
            if (!readLE(&data))
                return false;
            header.dims[i] = data;
        }
    }
    while (header.dims.size() < 6)
        header.dims.append(1);
    if (!mda->allocate(header.dims[0], header.dims[1], header.dims[2], header.dims[3], header.dims[4], header.dims[5]))
        return false;
    switch (header.data_type) {
    case Byte:
        return readData<unsigned char>(mda->dataPtr(), mda->totalSize());
    case Float32:
        return readData<float>(mda->dataPtr(), mda->totalSize());
    case Int16:
        return readData<int16_t>(mda->dataPtr(), mda->totalSize());
    case Int32:
        return readData<int32_t>(mda->dataPtr(), mda->totalSize());
    case UInt16:
        return readData<uint16_t>(mda->dataPtr(), mda->totalSize());
    case Float64:
        return readData<double>(mda->dataPtr(), mda->totalSize());
    case UInt32:
        return readData<uint32_t>(mda->dataPtr(), mda->totalSize());
    default:
        return false;
    }
    return true;
}

bool MdaIOHandlerMDA::write(const Mda& mda)
{
    Header header;
    if (format().toLower().startsWith("mda.")) {
        QByteArray subFormat = format().toLower().mid(4);
        if (subFormat == "byte") {
            header.data_type = Byte;
            header.num_bytes_per_entry = 1;
        }
        else if (subFormat == "float" || subFormat == "float32") {
            header.data_type = Float32;
            header.num_bytes_per_entry = 4;
        }
        else if (subFormat == "int16") {
            header.data_type = Int16;
            header.num_bytes_per_entry = 2;
        }
        else if (subFormat == "int" || subFormat == "int32") {
            header.data_type = Int32;
            header.num_bytes_per_entry = 4;
        }
        else if (subFormat == "uint16") {
            header.data_type = UInt16;
            header.num_bytes_per_entry = 2;
        }
        else if (subFormat == "double" || subFormat == "float64") {
            header.data_type = Float64;
            header.num_bytes_per_entry = 8;
        }
        else if (subFormat == "uint" || subFormat == "uint32") {
            header.data_type = UInt32;
            header.num_bytes_per_entry = 4;
        }
        else {
            return false;
        }
    }
    else {
        header.data_type = Float64;
        header.num_bytes_per_entry = 8;
    }
    header.num_dims = mda.ndims();
    header.dims.resize(header.num_dims);
    bool use64bitdims = false;
    for (int i = 0; i < header.num_dims; ++i) {
        header.dims[i] = mda.size(i);
        if (mda.size(i) > 2e9) {
            use64bitdims = true;
        }
    }
    if (use64bitdims)
        header.num_dims = -header.num_dims;
    // commit the header
    writeLE(header.data_type);
    writeLE(header.num_bytes_per_entry);
    writeLE(header.num_dims);
    if (use64bitdims) {
        for (uint64_t dim : header.dims)
            writeLE(dim);
    }
    else {
        for (int32_t dim : header.dims)
            writeLE(dim);
    }
    // commit data
    switch (header.data_type) {
    case Byte:
        return writeData<unsigned char>(mda.constDataPtr(), mda.totalSize());
    case Float32:
        return writeData<float>(mda.constDataPtr(), mda.totalSize());
    case Int16:
        return writeData<int16_t>(mda.constDataPtr(), mda.totalSize());
    case Int32:
        return writeData<int32_t>(mda.constDataPtr(), mda.totalSize());
    case UInt16:
        return writeData<uint16_t>(mda.constDataPtr(), mda.totalSize());
    case Float64:
        return writeData<double>(mda.constDataPtr(), mda.totalSize());
    case UInt32:
        return writeData<uint32_t>(mda.constDataPtr(), mda.totalSize());
    default:
        return false;
    }
}

bool MdaIOHandlerMDA::write(const Mda32& mda)
{
    Header header;
    if (format().toLower().startsWith("mda.")) {
        QByteArray subFormat = format().toLower().mid(4);
        if (subFormat == "byte") {
            header.data_type = Byte;
            header.num_bytes_per_entry = 1;
        }
        else if (subFormat == "float" || subFormat == "float32") {
            header.data_type = Float32;
            header.num_bytes_per_entry = 4;
        }
        else if (subFormat == "int16") {
            header.data_type = Int16;
            header.num_bytes_per_entry = 2;
        }
        else if (subFormat == "int" || subFormat == "int32") {
            header.data_type = Int32;
            header.num_bytes_per_entry = 4;
        }
        else if (subFormat == "uint16") {
            header.data_type = UInt16;
            header.num_bytes_per_entry = 2;
        }
        else if (subFormat == "double" || subFormat == "float64") {
            header.data_type = Float64;
            header.num_bytes_per_entry = 8;
        }
        else if (subFormat == "uint" || subFormat == "uint32") {
            header.data_type = UInt32;
            header.num_bytes_per_entry = 4;
        }
        else {
            return false;
        }
    }
    else {
        header.data_type = Float64;
        header.num_bytes_per_entry = 8;
    }
    header.num_dims = mda.ndims();
    header.dims.resize(header.num_dims);
    for (int i = 0; i < header.num_dims; ++i) {
        header.dims[i] = mda.size(i);
    }
    // commit the header
    writeLE(header.data_type);
    writeLE(header.num_bytes_per_entry);
    writeLE(header.num_dims);
    for (int32_t dim : header.dims)
        writeLE(dim);

    // commit data
    switch (header.data_type) {
    case Byte:
        return writeData<unsigned char>(mda.constDataPtr(), mda.totalSize());
    case Float32:
        return writeData<float>(mda.constDataPtr(), mda.totalSize());
    case Int16:
        return writeData<int16_t>(mda.constDataPtr(), mda.totalSize());
    case Int32:
        return writeData<int32_t>(mda.constDataPtr(), mda.totalSize());
    case UInt16:
        return writeData<uint16_t>(mda.constDataPtr(), mda.totalSize());
    case Float64:
        return writeData<double>(mda.constDataPtr(), mda.totalSize());
    case UInt32:
        return writeData<uint32_t>(mda.constDataPtr(), mda.totalSize());
    default:
        return false;
    }
}

MdaIOHandlerCSV::MdaIOHandlerCSV(QIODevice* device, const QByteArray& format)
{
    setDevice(device);
    setFormat(format);
}

bool MdaIOHandlerCSV::canRead() const
{
    if (!device())
        return false;
    if (format().toLower() == "csv")
        return true;
    if (format().isEmpty()) {
        // try to auto detect by file name
        if (QFileDevice* f = qobject_cast<QFileDevice*>(device())) {
            QFileInfo finfo(f->fileName());
            if (finfo.suffix().toLower() == "csv")
                return true;
        }
    }
    return false;
}

bool MdaIOHandlerCSV::canWrite() const
{
    if (!device())
        return false;
    if (format().toLower() == "csv")
        return true;
    return false;
}

bool MdaIOHandlerCSV::read(Mda* mda)
{
    QVector<QVector<double> > data;

    QTextStream stream(device());
    while (!device()->atEnd()) {
        QString line = stream.readLine().trimmed();
        if (line.isEmpty())
            continue;
        /// TODO: handle header
        if (line.startsWith('#'))
            continue;
        QStringList tokens = line.split(',', QString::SkipEmptyParts);
        QVector<double> row;
        foreach (const QString& token, tokens) {
            row.append(token.toDouble());
        }
        data.append(row);
    }
    if (data.isEmpty())
        return false;
    const int columnCount = data.first().count();
    const int rowCount = data.count();
    if (!mda->allocate(columnCount, rowCount))
        return false;
    for (int r = 0; r < rowCount; ++r) {
        const QVector<double>& row = data.at(r);
        for (int c = 0; c < columnCount; ++c) {
            double value = c < row.count() ? row.at(c) : 0.0;
            mda->setValue(value, c, r);
        }
    }
    return true;
}

bool MdaIOHandlerCSV::read(Mda32* mda)
{
    QVector<QVector<float> > data;

    QTextStream stream(device());
    while (!device()->atEnd()) {
        QString line = stream.readLine().trimmed();
        if (line.isEmpty())
            continue;
        /// TODO: handle header
        if (line.startsWith('#'))
            continue;
        QStringList tokens = line.split(',', QString::SkipEmptyParts);
        QVector<float> row;
        foreach (const QString& token, tokens) {
            row.append(token.toFloat());
        }
        data.append(row);
    }
    if (data.isEmpty())
        return false;
    const int columnCount = data.first().count();
    const int rowCount = data.count();
    if (!mda->allocate(columnCount, rowCount))
        return false;
    for (int r = 0; r < rowCount; ++r) {
        const QVector<float>& row = data.at(r);
        for (int c = 0; c < columnCount; ++c) {
            float value = c < row.count() ? row.at(c) : 0.0;
            mda->setValue(value, c, r);
        }
    }
    return true;
}

bool MdaIOHandlerCSV::write(const Mda& mda)
{
    QTextStream stream(device());
    for (size_t r = 0; r < (size_t)mda.N2(); ++r) {
        for (size_t c = 0; c < (size_t)mda.N1(); ++c) {
            double value = mda.value(c, r);
            if (c != 0) {
                stream << ", ";
            }
            stream << QString::number(value);
        }
        stream << endl;
    }
    return true;
}

bool MdaIOHandlerCSV::write(const Mda32& mda)
{
    QTextStream stream(device());
    for (size_t r = 0; r < (size_t)mda.N2(); ++r) {
        for (size_t c = 0; c < (size_t)mda.N1(); ++c) {
            double value = mda.value(c, r);
            if (c != 0) {
                stream << ", ";
            }
            stream << QString::number(value);
        }
        stream << endl;
    }
    return true;
}

MdaIOHandler* MdaIOHandlerCSVFactory::create(QIODevice* device, const QByteArray& format) const
{
    return new MdaIOHandlerCSV(device, format);
}

MdaIOHandler* MdaIOHandlerMDAFactory::create(QIODevice* device, const QByteArray& format) const
{
    return new MdaIOHandlerMDA(device, format);
}
