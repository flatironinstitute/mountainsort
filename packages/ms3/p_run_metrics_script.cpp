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
#include "p_run_metrics_script.h"
#include <QDebug>
#include <QJSEngine>
#include <QJsonObject>
#include <QJsonDocument>
#include "mlutil.h"
#include "textfile.h"

bool p_run_metrics_script(QString metrics, QString script, QString metrics_out)
{
    QJsonObject obj = QJsonDocument::fromJson(TextFile::read(metrics).toUtf8()).object();
    QJSEngine engine;

    JSConsole console;
    QJSValue consoleObj = engine.newQObject(&console);
    engine.globalObject().setProperty("console", consoleObj);

    QString program = TextFile::read(script) + "\n\n";
    program += QString("var _obj=JSON.parse('%1');\n").arg((QString)QJsonDocument(obj).toJson(QJsonDocument::Compact));
    program += QString("main(_obj);\n");
    program += QString("_obj_json=JSON.stringify(_obj);\n");
    QJSValue val = engine.evaluate(program);
    if (val.isError()) {
        qWarning() << QString("Problem executing script on line %1").arg(val.property("lineNumber").toString());
        qWarning() << val.toString();
        return false;
    }
    else {
        QString obj_json = engine.globalObject().property("_obj_json").toString();
        obj = QJsonDocument::fromJson(obj_json.toUtf8()).object();
        TextFile::write(metrics_out, QJsonDocument(obj).toJson());
        return true;
    }
}

JSConsole::JSConsole(QObject* parent)
    : QObject(parent)
{
}

void JSConsole::log(QString msg)
{
    qDebug().noquote() << msg;
}
