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
#ifndef P_COMBINE_FIRING_SEGMENTS_H
#define P_COMBINE_FIRING_SEGMENTS_H

#include <QStringList>

struct P_combine_firing_segments_opts {
    int clip_size = 60;
    double match_score_threshold = 0.6;
    int offset_search_radius = 10;
    int num_comparison_events = 1000;
};

bool p_combine_firing_segments(QString timeseries, QStringList firings_list, QString firings_out, P_combine_firing_segments_opts opts);

#endif // P_COMBINE_FIRING_SEGMENTS_H
