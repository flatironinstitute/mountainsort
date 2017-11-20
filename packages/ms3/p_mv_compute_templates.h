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
#ifndef MV_COMPUTE_TEMPLATES_H
#define MV_COMPUTE_TEMPLATES_H

#include <QString>

bool mv_compute_templates(const QString& timeseries_path, const QString& firings_path, const QString& templates_out_path, const QString& stdevs_out_path, int clip_size);
bool mv_subfirings(QString firings_path, QString firings_out_path, QVector<int> labels, int max_per_label);

#endif // MV_COMPUTE_TEMPLATES_H
