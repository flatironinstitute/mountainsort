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

