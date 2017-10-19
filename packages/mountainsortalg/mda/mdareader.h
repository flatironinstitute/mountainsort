#ifndef MDAREADER_H
#define MDAREADER_H

#include "mda/mda.h"
#include "mda/mda32.h"

class QIODevice;
class MdaReaderPrivate;
class MdaReader {
public:
    MdaReader();
    MdaReader(QIODevice*, const QByteArray& format = QByteArray());
    MdaReader(const QString& fileName, const QByteArray& format = QByteArray());
    ~MdaReader();
    bool canRead() const;
    QIODevice* device() const;
    QString fileName() const;
    QByteArray format() const;
    Mda read();
    Mda32 read32();
    bool read(Mda*);
    bool read(Mda32* mda);
    void setDevice(QIODevice*);
    void setFileName(const QString&);

private:
    MdaReaderPrivate* d;
};

class MdaWriterPrivate;
class MdaWriter {
public:
    MdaWriter();
    MdaWriter(QIODevice*, const QByteArray& format);
    MdaWriter(const QString& fileName, const QByteArray& format);
    ~MdaWriter();
    bool canWrite() const;
    QIODevice* device() const;
    QString fileName() const;
    QByteArray format() const;
    void setDevice(QIODevice*);
    void setFileName(const QString& fileName);
    bool write(const Mda&);
    bool write(const Mda32&);

private:
    MdaWriterPrivate* d;
};

#endif // MDAREADER_H
