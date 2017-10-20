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
