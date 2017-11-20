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
#ifndef CONSOLIDATE_CLUSTERS_H
#define CONSOLIDATE_CLUSTERS_H

#include <QVector>
#include "mda32.h"

struct Consolidate_clusters_opts {
    double consolidation_factor = 0.9;
};

void consolidate_clusters(QVector<bigint>& event_inds, QVector<double>& timepoints, QVector<int>& labels, const Mda32& templates, Consolidate_clusters_opts opts);
QMap<int, int> consolidate_clusters(const Mda32& templates, Consolidate_clusters_opts opts);

#endif // CONSOLIDATE_CLUSTERS_H
