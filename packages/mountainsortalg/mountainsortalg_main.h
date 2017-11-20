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
#ifndef MOUNTAINSORTALG_MAIN_H
#define MOUNTAINSORTALG_MAIN_H

#include <QJsonObject>
#include <QVariant>

QJsonObject get_spec();

struct ProcessorSpecFile {
    QString name;
    QString description;
    bool optional = false;

    QJsonObject get_spec();
};

struct ProcessorSpecParam {
    QString name;
    QString description;
    bool optional = false;
    QVariant default_value;

    QJsonObject get_spec();
};

struct ProcessorSpec {
    ProcessorSpec(QString name, QString version_in)
    {
        processor_name = name;
        version = version_in;
    }

    QString processor_name;
    QString version;
    QString description;
    QList<ProcessorSpecFile> inputs;
    QList<ProcessorSpecFile> outputs;
    QList<ProcessorSpecParam> parameters;

    void addInputs(QString name1, QString name2 = "", QString name3 = "", QString name4 = "", QString name5 = "");
    void addOptionalInputs(QString name1, QString name2 = "", QString name3 = "", QString name4 = "", QString name5 = "");
    void addOutputs(QString name1, QString name2 = "", QString name3 = "", QString name4 = "", QString name5 = "");
    void addOptionalOutputs(QString name1, QString name2 = "", QString name3 = "", QString name4 = "", QString name5 = "");
    void addRequiredParameters(QString name1, QString name2 = "", QString name3 = "", QString name4 = "", QString name5 = "");
    void addOptionalParameters(QString name1, QString name2 = "", QString name3 = "", QString name4 = "", QString name5 = "");
    void addRequiredParameter(QString name, QString description = "");
    void addOptionalParameter(QString name, QString description = "", QVariant default_value = QVariant());

    void addInput(QString name, QString description = "", bool optional = false);
    void addOutput(QString name, QString description = "", bool optional = false);

    QJsonObject get_spec();
};

#endif // MOUNTAINSORTALG_MAIN_H

