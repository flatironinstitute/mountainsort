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

#ifndef EXTRACT_CLIPS_H
#define EXTRACT_CLIPS_H

#include "mda.h"
#include "diskreadmda.h"
#include "diskreadmda32.h"

bool extract_clips(const QString& timeseries_path, const QString& firings_path, const QString& clips_path, int clip_size, const QList<int>& channels, double t1, double t2);
//bool extract_clips_features(const QString& timeseries_path, const QString& firings_path, const QString& features_path, int clip_size, int num_features);

Mda extract_clips(const DiskReadMda& X, const QVector<double>& times, int clip_size);
Mda32 extract_clips(const DiskReadMda32& X, const QVector<double>& times, int clip_size);
Mda extract_clips(const DiskReadMda& X, const QVector<double>& times, const QVector<int>& channels, int clip_size);
Mda32 extract_clips(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& channels, int clip_size);

#endif // EXTRACT_CLIPS_H
