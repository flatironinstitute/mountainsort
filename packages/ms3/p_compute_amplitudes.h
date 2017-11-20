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
#ifndef P_COMPUTE_AMPLITUDES_H
#define P_COMPUTE_AMPLITUDES_H

#include <QString>
#include "mlutil.h"

struct P_compute_amplitudes_opts {
    int central_channel = 0;
};

bool p_compute_amplitudes(QString timeseries, QString event_times, QString amplitudes_out, P_compute_amplitudes_opts opts);

#endif // P_COMPUTE_AMPLITUDES_H
