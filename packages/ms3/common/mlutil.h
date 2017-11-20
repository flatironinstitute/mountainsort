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

