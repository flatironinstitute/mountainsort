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
#ifndef P_EXTRACT_CLIPS_H
#define P_EXTRACT_CLIPS_H

#include <QVariantMap>

bool p_extract_clips(QStringList timeseries_list, QString event_times, const QList<int>& channels, QString clips_out, const QVariantMap& params);

bool p_mv_extract_clips(QStringList timeseries_list, QString firings, const QList<int>& channels, QString clips_out, const QVariantMap& params);

bool p_mv_extract_clips_features(QString timeseries, QString firings, QString features_out, int clip_size, int num_features, int subtract_mean);


#endif // P_EXTRACT_CLIPS_H
