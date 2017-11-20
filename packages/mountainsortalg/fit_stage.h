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
#ifndef FIT_STAGE_H
#define FIT_STAGE_H

#include <QVector>
#include "mda32.h"

struct Fit_stage_opts {
    double time_channel_mask_thresh = 0.1;
};

QVector<bigint> fit_stage(Mda32& X, const QVector<double>& times, const QVector<int>& labels, Mda32& templates, Fit_stage_opts opts);

#endif // FIT_STAGE_H
