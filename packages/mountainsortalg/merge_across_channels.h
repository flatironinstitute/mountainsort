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
#ifndef MERGE_ACROSS_CHANNELS_H
#define MERGE_ACROSS_CHANNELS_H

#include <mda32.h>

struct Merge_across_channels_opts {
    int clip_size = 100;
    double min_peak_ratio_to_consider = 0.3; //reduced on 5/26/17 0.7->0.3
    double event_fraction_threshold = 0.3; //reduced on 5/26/17 0.5->0.3
};

void merge_across_channels(QVector<double>& times, QVector<int>& labels, QVector<int>& central_channels, Mda32& templates, Merge_across_channels_opts opts);

#endif // MERGE_ACROSS_CHANNELS_H
