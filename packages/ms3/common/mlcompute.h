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
#ifndef MLCOMPUTE_H
#define MLCOMPUTE_H

#include <QVector>
#include <stdlib.h>
#include "mlvector.h"
typedef int64_t bigint;

namespace MLCompute {
double min(const QVector<double>& X);
double max(const QVector<double>& X);
double sum(const QVector<double>& X);
double mean(const QVector<double>& X);
double median(const QVector<double>& X);
double stdev(const QVector<double>& X);
double norm(const QVector<double>& X);
double dotProduct(const QVector<double>& X1, const QVector<double>& X2);
double correlation(const QVector<double>& X1, const QVector<double>& X2);

double norm(const QVector<float>& X);
double dotProduct(const QVector<float>& X1, const QVector<float>& X2);

template <typename T>
T max(const QVector<T>& X);
double min(bigint N, const double* X);
double max(bigint N, const double* X);
double sum(bigint N, const double* X);
double mean(bigint N, const double* X);
double dotProduct(bigint N, const double* X1, const double* X2);
double norm(bigint N, const double* X);
double min(bigint N, const float* X);
double max(bigint N, const float* X);
double sum(bigint N, const float* X);
double mean(bigint N, const float* X);
double stdev(bigint N, const float* X);
double dotProduct(bigint N, const float* X1, const float* X2);
double correlation(bigint N, const float* X1, const float* X2);
double norm(bigint N, const float* X);

double min(const MLVector<double>& X);
double max(const MLVector<double>& X);
int min(const MLVector<int>& X);
int max(const MLVector<int>& X);
}

template <typename T>
T MLCompute::max(const QVector<T>& X)
{
    return *std::max_element(X.constBegin(), X.constEnd());
}

#endif // MLCOMPUTE_H

