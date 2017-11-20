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
#ifndef DIMENSION_REDUCE_CLIPS_H
#define DIMENSION_REDUCE_CLIPS_H

#include "mda32.h"

void dimension_reduce_clips(Mda32& reduced_clips, const Mda32& clips, bigint num_features_per_channel, bigint max_samples);
void dimension_reduce_clips(QString clips_path, QString reduced_clips_path, bigint num_features_per_channel, bigint max_samples);

#endif // DIMENSION_REDUCE_CLIPS_H
