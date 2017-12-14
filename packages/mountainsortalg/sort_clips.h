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
#ifndef SORT_CLIPS_H
#define SORT_CLIPS_H

#include <mda32.h>

struct Sort_clips_opts {
    int num_features = 10;
    bigint max_samples = 10000;
    double isocut_threshold = 1;
    int K_init = 200;
    QString debug_description;
};

QVector<int> sort_clips(const Mda32& clips, const Sort_clips_opts& opts);

#endif // SORT_CLIPS_H
