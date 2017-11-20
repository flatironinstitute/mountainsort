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
#include "p_cluster_metrics.h"
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <diskreadmda32.h>
#include <mda.h>
#include "mlutil.h"
#include "textfile.h"

namespace P_cluster_metrics {
QJsonObject get_cluster_metrics(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int k, Cluster_metrics_opts opts);
}

bool p_cluster_metrics(QString timeseries_path, QString firings_path, QString metrics_out_path, Cluster_metrics_opts opts)
{
    QJsonArray clusters;

    DiskReadMda32 X(timeseries_path);
    Mda firings(firings_path);

    QVector<double> times(firings.N2());
    QVector<int> labels(firings.N2());
    for (bigint i = 0; i < firings.N2(); i++) {
        times[i] = firings.value(1, i);
        labels[i] = firings.value(2, i);
    }

    QSet<int> used_labels_set;
    for (bigint i = 0; i < labels.count(); i++) {
        used_labels_set.insert(labels[i]);
    }
    QList<int> used_labels = used_labels_set.toList();
    qSort(used_labels);

    foreach (int k, used_labels) {
        QJsonObject tmp = P_cluster_metrics::get_cluster_metrics(X, times, labels, k, opts);
        clusters.push_back(tmp);
    }

    QJsonObject obj;
    obj["clusters"] = clusters;
    QString json = QJsonDocument(obj).toJson(QJsonDocument::Indented);

    return TextFile::write(metrics_out_path, json);
}

namespace P_cluster_metrics {
QJsonObject get_cluster_metrics(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int k, Cluster_metrics_opts opts)
{
    (void)X;
    QVector<double> times_k;
    for (bigint i = 0; i < labels.count(); i++) {
        if (labels[i] == k)
            times_k << times[i];
    }

    qSort(times_k);

    QJsonObject metrics;
    metrics["num_events"] = times_k.count();
    if (!opts.samplerate)
        opts.samplerate = 1;
    double t1_sec = times_k.value(0) / opts.samplerate;
    double t2_sec = times_k.value(times_k.count() - 1) / opts.samplerate;
    double dur_sec = t2_sec - t1_sec;
    double firing_rate = 0;
    if (dur_sec > 0) {
        firing_rate = times_k.count() / dur_sec;
    }
    metrics["t1_sec"] = t1_sec;
    metrics["t2_sec"] = t2_sec;
    metrics["dur_sec"] = dur_sec;
    metrics["firing_rate"] = firing_rate;

    QJsonObject obj;
    obj["label"] = (double)k;
    obj["metrics"] = metrics;

    return obj;
}
}

bool p_combine_cluster_metrics(QStringList metrics_list, QString metrics_out)
{
    QMap<int, QJsonObject> all_clusters_metrics;
    foreach (QString fname, metrics_list) {
        QString json = TextFile::read(fname);
        QJsonObject obj = QJsonDocument::fromJson(json.toUtf8()).object();
        QJsonArray clusters = obj["clusters"].toArray();
        for (int i = 0; i < clusters.count(); i++) {
            QJsonObject tmp = clusters[i].toObject();
            int label = tmp["label"].toInt();
            if (!all_clusters_metrics.contains(label)) {
                all_clusters_metrics[label] = QJsonObject();
            }
            QStringList names = tmp["metrics"].toObject().keys();
            foreach (QString name, names) {
                QJsonValue val = tmp["metrics"].toObject()[name];
                all_clusters_metrics[label][name] = val;
            }
        }
    }
    QJsonArray ret_clusters;
    QList<int> labels = all_clusters_metrics.keys();
    qSort(labels);
    foreach (int k, labels) {
        QJsonObject C;
        C["label"] = k;
        C["metrics"] = all_clusters_metrics[k];
        ret_clusters << C;
    }
    QJsonObject ret;
    ret["clusters"] = ret_clusters;
    return TextFile::write(metrics_out, QJsonDocument(ret).toJson(QJsonDocument::Indented));
}
