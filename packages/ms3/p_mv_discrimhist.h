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
#ifndef MV_DISCRIMHIST_H
#define MV_DISCRIMHIST_H

#include <QList>
#include <QString>
#include <QVector>
#include <diskreadmda.h>
#include <diskreadmda32.h>

struct mv_discrimhist_opts {
    QList<int> clusters;
    /// TODO clip_size is hard-coded here
    int clip_size = 80;
    QString method = "centroid"; //centroid or svm
    int num_features = 0;
};

bool mv_discrimhist(QString timeseries_path, QString firings_path, QString output_path, mv_discrimhist_opts opts);
bool get_discrimhist_data(QVector<double>& ret1, QVector<double>& ret2, const DiskReadMda32& timeseries, const DiskReadMda& firings, int k1, int k2, int clip_size, QString method);

#endif // MV_DISCRIMHIST_H
