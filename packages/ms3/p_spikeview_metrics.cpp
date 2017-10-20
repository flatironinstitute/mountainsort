#include "p_spikeview_metrics.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <mda.h>
#include <QJsonObject>
#include "textfile.h"

struct ClusterInfo {
    bigint num_events = 0;
    double tmin = 0;
    double tmax = 0;
};

bool p_spikeview_metrics1(QString firings, QString metrics_out, P_spikeview_metrics1_opts opts)
{
    Mda FF(firings);
    bigint L = FF.N2();
    QVector<double> times(L);
    QVector<int> labels(L);

    for (bigint i = 0; i < L; i++) {
        times[i] = FF.value(1, i);
        labels[i] = FF.value(2, i);
    }

    QMap<int, ClusterInfo> infos;
    for (bigint i = 0; i < L; i++) {
        int k = labels[i];
        infos[k].num_events++;
        if (infos[k].num_events > 1) {
            infos[k].tmin = qMin(infos[k].tmin, times[i]);
            infos[k].tmax = qMax(infos[k].tmax, times[i]);
        }
        else {
            infos[k].tmin = times[i];
            infos[k].tmax = times[i];
        }
    }

    QList<int> cluster_numbers = infos.keys();

    double samplerate = opts.samplerate;
    if (!samplerate)
        samplerate = 1;

    QJsonArray clusters;
    foreach (int k, cluster_numbers) {
        QJsonObject tmp;
        tmp["label"] = k;
        QJsonObject metrics;
        metrics["sv.num_events"] = (long long)infos[k].num_events;
        metrics["sv.tmin_sec"] = infos[k].tmin / samplerate;
        metrics["sv.tmax_sec"] = infos[k].tmax / samplerate;
        metrics["sv.dur_sec"] = (infos[k].tmax - infos[k].tmin) / samplerate;
        tmp["metrics"] = metrics;
        clusters.push_back(tmp);
    }

    {
        QJsonObject obj;
        obj["clusters"] = clusters;
        QString json = QJsonDocument(obj).toJson(QJsonDocument::Indented);
        if (!TextFile::write(metrics_out, json))
            return false;
    }

    return true;
}
