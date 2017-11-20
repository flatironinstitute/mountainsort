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
#ifndef P_CONCAT_FIRINGS_H
#define P_CONCAT_FIRINGS_H

#include <QString>
#include <QStringList>

bool p_concat_firings(QStringList timeseries_list, QStringList firings_list, QString timeseries_out, QString firings_out);
bool p_concat_event_times(QStringList event_times_list, QString event_times_out);

#endif // P_CONCAT_FIRINGS_H
