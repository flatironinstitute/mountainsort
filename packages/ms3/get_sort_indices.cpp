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
#include "get_sort_indices.h"

#include <QVector>

QList<int> get_sort_indices(const QList<int>& X)
{
    QList<int> result;
    result.reserve(X.size());
    for (int i = 0; i < X.size(); ++i)
        result << i;
    std::stable_sort(result.begin(), result.end(),
        [&X](int i1, int i2) { return X[i1] < X[i2]; });
    return result;
}

QVector<int> get_sort_indices(const QVector<int>& X)
{
    QVector<int> result(X.size());
    for (int i = 0; i < X.size(); ++i)
        result << i;
    std::stable_sort(result.begin(), result.end(),
        [&X](int i1, int i2) { return X[i1] < X[i2]; });
    return result;
}

QList<int> get_sort_indices(const QVector<double>& X)
{
    QList<int> result;
    result.reserve(X.size());
    for (int i = 0; i < X.size(); ++i)
        result << i;
    std::stable_sort(result.begin(), result.end(),
        [&X](int i1, int i2) { return X[i1] < X[i2]; });
    return result;
}

QList<bigint> get_sort_indices_bigint(const QVector<double>& X)
{
    QList<bigint> result;
    result.reserve(X.size());
    for (bigint i = 0; i < X.size(); ++i)
        result << i;
    std::stable_sort(result.begin(), result.end(),
        [&X](bigint i1, bigint i2) { return X[i1] < X[i2]; });
    return result;
}

MLVector<bigint> get_sort_indices(const MLVector<double>& X)
{
    MLVector<bigint> result;
    result.reserve(X.size());
    for (bigint i = 0; i < X.count(); ++i)
        result << i;
    std::stable_sort(result.begin(), result.end(),
        [&X](bigint i1, bigint i2) { return X[i1] < X[i2]; });
    return result;
}
