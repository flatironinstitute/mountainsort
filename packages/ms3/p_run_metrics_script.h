#ifndef P_RUN_METRIC_SCRIPT_H
#define P_RUN_METRIC_SCRIPT_H

#include <QString>
#include <QObject>

bool p_run_metrics_script(QString metrics, QString script, QString metrics_out);

class JSConsole : public QObject {
    Q_OBJECT
public:
    explicit JSConsole(QObject* parent = 0);

signals:

public slots:
    void log(QString msg);
};

#endif // P_RUN_METRIC_SCRIPT_H
