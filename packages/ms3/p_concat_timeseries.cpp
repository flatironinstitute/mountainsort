#include "p_concat_timeseries.h"

#include <QTime>
#include <diskreadmda32.h>
#include <diskwritemda.h>

bool p_concat_timeseries(QStringList timeseries_list, QString timeseries_out)
{

    DiskReadMda32 X;
    X.setConcatPaths(2, timeseries_list);

    DiskWriteMda Y;
    Y.open(X.mdaioHeader().data_type, timeseries_out, X.N1(), X.N2());
    bigint chunk_size = 1e6;
    QTime timer;
    timer.start();
    for (bigint j = 0; j < X.N2(); j += chunk_size) {
        if (timer.elapsed() > 3000) {
            qDebug().noquote() << QString("Processed %1 of %2 timepoints (%3%)").arg(j).arg(X.N2()).arg(j * 100.0 / X.N2());
            timer.start();
        }
        bigint size0 = chunk_size;
        if (j + size0 >= X.N2())
            size0 = X.N2() - j;
        Mda32 chunk0;
        if (!X.readChunk(chunk0, 0, j, X.N1(), size0)) {
            qWarning() << "Error reading chunk.";
        }
        if (!Y.writeChunk(chunk0, 0, j)) {
            qWarning() << "Error writing chunk.";
            return false;
        }
    }

    return true;
}
