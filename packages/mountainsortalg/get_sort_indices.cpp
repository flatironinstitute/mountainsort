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
