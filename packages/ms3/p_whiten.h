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
#ifndef P_WHITEN_H
#define P_WHITEN_H

#include <QString>

struct Whiten_opts {
    double quantization_unit = 0;

    bool requirements_only = false;
    double expected_peak_ram_mb = -1;
};

bool p_whiten(QString timeseries, QString timeseries_out, Whiten_opts &opts);
bool p_compute_whitening_matrix(QStringList timeseries_list, const QList<int>& channels, QString whitening_matrix_out, Whiten_opts opts);
bool p_apply_whitening_matrix(QString timeseries, QString whitening_matrix, QString timeseries_out, Whiten_opts opts);
bool p_whiten_clips(QString clips, QString whitening_matrix, QString clips_out, Whiten_opts opts);

#endif // P_WHITEN_H
