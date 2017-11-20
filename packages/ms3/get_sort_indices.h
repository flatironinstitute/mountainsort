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

#ifndef GET_SORT_INDICES_H
#define GET_SORT_INDICES_H

#include <QList>
#include "mlvector.h"
typedef int64_t bigint;

QList<int> get_sort_indices(const QList<int>& X);
QVector<int> get_sort_indices(const QVector<int>& X);
QList<int> get_sort_indices(const QVector<double>& X);
QList<bigint> get_sort_indices_bigint(const QVector<double>& X);

MLVector<bigint> get_sort_indices(const MLVector<double>& X);

#endif // GET_SORT_INDICES_H
