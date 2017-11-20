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
#ifndef P_BANJOVIEW_CROSS_CORRELOGRAMS_H
#define P_BANJOVIEW_CROSS_CORRELOGRAMS_H

#include <QList>
#include <QString>

enum P_banjoview_cross_correlograms_mode {
    Autocorrelograms,
    Matrix_of_cross_correlograms
};

struct P_banjoview_cross_correlograms_opts {
    P_banjoview_cross_correlograms_mode mode;
    QList<int> clusters;
    double max_dt_msec = 100;
    double bin_size_msec = 1;
    double samplerate = 30000;
};

bool p_banjoview_cross_correlograms(QString firings, QString correlograms_out, P_banjoview_cross_correlograms_opts opts);

#endif // P_BANJOVIEW_CROSS_CORRELOGRAMS_H
