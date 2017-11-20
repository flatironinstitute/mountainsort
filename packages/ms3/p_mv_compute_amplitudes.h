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
#ifndef P_MV_COMPUTE_AMPLITUDES_H
#define P_MV_COMPUTE_AMPLITUDES_H

#include <QString>

struct p_mv_compute_amplitudes_opts {
    int clip_size = 50;
};

bool p_mv_compute_amplitudes(QString timeseries_path, QString firings_path, QString firings_out_path, p_mv_compute_amplitudes_opts opts);

#endif // P_MV_COMPUTE_AMPLITUDES_H
