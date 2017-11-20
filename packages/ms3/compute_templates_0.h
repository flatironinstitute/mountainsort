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

#ifndef COMPUTE_TEMPLATES_0_H
#define COMPUTE_TEMPLATES_0_H

#include "diskreadmda.h"
#include "diskreadmda32.h"

Mda compute_templates_0(const DiskReadMda& X, Mda& firings, int clip_size);
Mda compute_templates_0(const DiskReadMda& X, const QVector<double>& times, const QVector<int>& labels, int clip_size);
void compute_templates_stdevs(Mda& ret_templates, Mda& ret_stdevs, DiskReadMda& X, const QVector<double>& times, const QVector<int>& labels, int clip_size);

Mda32 compute_templates_0(const DiskReadMda32& X, Mda& firings, int clip_size);
Mda32 compute_templates_0(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int clip_size);

Mda32 compute_templates_in_parallel(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int clip_size);

#endif // COMPUTE_TEMPLATES_0_H
