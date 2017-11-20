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
#include "mlcompute.h"
#include <cmath>
using std::sqrt;

double MLCompute::min(const QVector<double>& X)
{
    if (X.count() == 0)
        return 0;
    return *std::min_element(X.constBegin(), X.constEnd());
}

double MLCompute::max(const QVector<double>& X)
{
    if (X.count() == 0)
        return 0;
    return *std::max_element(X.constBegin(), X.constEnd());
}

double MLCompute::min(const MLVector<double>& X)
{
    if (X.count() == 0)
        return 0;
    return *std::min_element(X.begin(), X.end());
}

double MLCompute::max(const MLVector<double>& X)
{
    if (X.count() == 0)
        return 0;
    return *std::max_element(X.begin(), X.end());
}

int MLCompute::min(const MLVector<int>& X)
{
    if (X.count() == 0)
        return 0;
    return *std::min_element(X.begin(), X.end());
}

int MLCompute::max(const MLVector<int>& X)
{
    if (X.count() == 0)
        return 0;
    return *std::max_element(X.begin(), X.end());
}

double MLCompute::sum(const QVector<double>& X)
{
    return std::accumulate(X.constBegin(), X.constEnd(), 0.0);
}

double MLCompute::mean(const QVector<double>& X)
{
    if (X.isEmpty())
        return 0;
    double s = sum(X);
    return s / X.count();
}

double MLCompute::stdev(const QVector<double>& X)
{
    double sumsqr = std::inner_product(X.constBegin(), X.constEnd(), X.constBegin(), 0.0);
    double sum = std::accumulate(X.constBegin(), X.constEnd(), 0.0);
    int ct = X.count();
    if (ct >= 2) {
        return sqrt((sumsqr - sum * sum / ct) / (ct - 1));
    }
    else
        return 0;
}

double MLCompute::dotProduct(const QVector<double>& X1, const QVector<double>& X2)
{
    if (X1.count() != X2.count())
        return 0;
    return std::inner_product(X1.constBegin(), X1.constEnd(), X2.constBegin(), 0.0);
}

double MLCompute::norm(const QVector<double>& X)
{
    return sqrt(dotProduct(X, X));
}

double MLCompute::dotProduct(const QVector<float>& X1, const QVector<float>& X2)
{
    if (X1.count() != X2.count())
        return 0;
    return std::inner_product(X1.constBegin(), X1.constEnd(), X2.constBegin(), 0.0);
}

double MLCompute::norm(const QVector<float>& X)
{
    return sqrt(dotProduct(X, X));
}

double MLCompute::correlation(const QVector<double>& X1, const QVector<double>& X2)
{
    if (X1.count() != X2.count())
        return 0;
    int N = X1.count();
    double mean1 = mean(X1);
    double stdev1 = stdev(X1);
    double mean2 = mean(X2);
    double stdev2 = stdev(X2);
    if ((stdev1 == 0) || (stdev2 == 0))
        return 0;
    QVector<double> Y1(N);
    QVector<double> Y2(N);
    for (bigint i = 0; i < N; i++) {
        Y1[i] = (X1[i] - mean1) / stdev1;
        Y2[i] = (X2[i] - mean2) / stdev2;
    }
    return dotProduct(Y1, Y2);
}

double MLCompute::norm(bigint N, const double* X)
{
    return sqrt(dotProduct(N, X, X));
}

double MLCompute::dotProduct(bigint N, const double* X1, const double* X2)
{
    return std::inner_product(X1, X1 + N, X2, 0.0);
}

double MLCompute::dotProduct(bigint N, const float* X1, const float* X2)
{

    return std::inner_product(X1, X1 + N, X2, 0.0);
}

double MLCompute::sum(bigint N, const double* X)
{
    return std::accumulate(X, X + N, 0.0);
}

double MLCompute::mean(bigint N, const double* X)
{
    if (!N)
        return 0;
    return sum(N, X) / N;
}

double MLCompute::max(bigint N, const double* X)
{
    return N ? *std::max_element(X, X + N) : 0;
}

double MLCompute::min(bigint N, const double* X)
{
    return N ? *std::min_element(X, X + N) : 0;
}

double MLCompute::min(bigint N, const float* X)
{
    return N ? *std::min_element(X, X + N) : 0;
}

double MLCompute::max(bigint N, const float* X)
{
    return N ? *std::max_element(X, X + N) : 0;
}

double MLCompute::sum(bigint N, const float* X)
{
    return std::accumulate(X, X + N, 0.0);
}

double MLCompute::mean(bigint N, const float* X)
{
    if (!N)
        return 0;
    return sum(N, X) / N;
}

double MLCompute::norm(bigint N, const float* X)
{
    return sqrt(dotProduct(N, X, X));
}

double MLCompute::median(const QVector<double>& X)
{
    if (X.isEmpty())
        return 0;
    QVector<double> Y = X;
    qSort(Y);
    if (Y.count() % 2 == 0) {
        return (Y[Y.count() / 2 - 1] + Y[Y.count() / 2]) / 2;
    }
    else {
        return Y[Y.count() / 2];
    }
}

double MLCompute::correlation(bigint N, const float* X1, const float* X2)
{
    if (N <= 1)
        return 0;
    double mean1 = mean(N, X1);
    double stdev1 = stdev(N, X1);
    double mean2 = mean(N, X2);
    double stdev2 = stdev(N, X2);
    if ((stdev1 == 0) || (stdev2 == 0))
        return 0;
    QVector<double> Y1(N);
    QVector<double> Y2(N);
    for (bigint i = 0; i < N; i++) {
        Y1[i] = (X1[i] - mean1) / stdev1;
        Y2[i] = (X2[i] - mean2) / stdev2;
    }
    return dotProduct(Y1, Y2) / (N - 1);
}

double MLCompute::stdev(bigint N, const float* X)
{
    double sumsqr = dotProduct(N, X, X);
    double sum0 = sum(N, X);
    if (N >= 2) {
        return sqrt((sumsqr - sum0 * sum0 / N) / (N - 1));
    }
    else
        return 0;
}

