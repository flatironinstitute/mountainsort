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
