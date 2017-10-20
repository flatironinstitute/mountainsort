#ifndef MLUTIL_H
#define MLUTIL_H

#include <QJsonValue>
#include <stdlib.h>
typedef int64_t bigint;

namespace MLUtil {
QJsonValue toJsonValue(const QByteArray& X);
QJsonValue toJsonValue(const QList<int>& X);
QJsonValue toJsonValue(const QVector<int>& X);
QJsonValue toJsonValue(const QVector<double>& X);
void fromJsonValue(QByteArray& X, const QJsonValue& val);
void fromJsonValue(QList<int>& X, const QJsonValue& val);
void fromJsonValue(QVector<int>& X, const QJsonValue& val);
void fromJsonValue(QVector<double>& X, const QJsonValue& val);
QByteArray readByteArray(const QString& path);
bool writeByteArray(const QString& path, const QByteArray& X);
QStringList toStringList(const QVariant& val); //val is either a string or a QVariantList
QString makeRandomId(int numchars);
QList<int> stringListToIntList(const QStringList& list);
QList<bigint> stringListToBigIntList(const QStringList& list);
QStringList intListToStringList(const QList<int>& list);
};

#endif // MLUTIL_H

