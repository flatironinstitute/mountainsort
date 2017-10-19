/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

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
