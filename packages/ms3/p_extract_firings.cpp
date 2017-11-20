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
#include "p_extract_firings.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <diskreadmda.h>
#include "mlutil.h"
#include "textfile.h"

bool p_extract_firings(QString firings, QString metrics, QString firings_out, P_extract_firings_opts opts)
{
    QSet<int> labels_to_exclude;
    QJsonParseError err;
    QJsonObject obj;

    if (!metrics.isEmpty()) {
        obj = QJsonDocument::fromJson(TextFile::read(metrics).toUtf8(), &err).object();
        if (err.error != QJsonParseError::NoError) {
            qWarning() << "Error parsing json file: " + metrics + " : " + err.errorString();
            return false;
        }
    }

    QJsonArray clusters = obj["clusters"].toArray();
    for (int i = 0; i < clusters.count(); i++) {
        QJsonObject cluster = clusters[i].toObject();
        int k = cluster["label"].toInt();
        QJsonArray tags = cluster["tags"].toArray();
        QSet<QString> tags_set;
        for (int j = 0; j < tags.count(); j++) {
            tags_set.insert(tags[j].toString());
        }
        foreach (QString str, opts.exclusion_tags) {
            if (tags_set.contains(str)) {
                labels_to_exclude.insert(k);
            }
        }
    }
    qDebug().noquote() << "Excluding clusters: " << labels_to_exclude.toList();

    QVector<bigint> inds_to_use;
    QSet<int> labels_to_include;

    DiskReadMda FF(firings);
    bigint L = FF.N2();
    for (bigint i = 0; i < L; i++) {
        double t0 = FF.value(1, i);
        if ((opts.t1 <= t0) && ((opts.t2 < 0) || (t0 <= opts.t2))) {
            int k = FF.value(2, i);
            if (!labels_to_exclude.contains(k)) {
                if (labels_to_include.contains(k)) {
                    inds_to_use << i;
                }
                else {
                    if (opts.clusters.isEmpty()) {
                        labels_to_include.insert(k);
                    }
                    else {
                        if (opts.clusters.contains(k))
                            labels_to_include.insert(k);
                        else
                            labels_to_exclude.insert(k);
                    }
                    if (labels_to_include.contains(k))
                        inds_to_use << i;
                }
            }
        }
    }

    bigint L2 = inds_to_use.count();

    qDebug().noquote() << QString("Using %1 of %2 events (%3%)").arg(L2).arg(L).arg(L2 * 100.0 / L);

    Mda FF2(FF.N1(), L2);
    for (bigint j = 0; j < L2; j++) {
        for (int r = 0; r < FF.N1(); r++) {
            FF2.set(FF.value(r, inds_to_use[j]), r, j);
        }
    }

    return FF2.write64(firings_out);
}

/*
Example .json structure
{
    "clusters": [
        {
            "label": 1,
            "metrics": {
                "dur_sec": 3599.7703666666666,
                "firing_rate": 4.681965315361643,
                "isolation": 0.9917743830787309,
                "noise_overlap": 0.3683901292596945,
                "num_events": 16854,
                "overlap_cluster": 20,
                "peak_amp": 3.2513680458068848,
                "peak_noise": 1.0339038372039795,
                "peak_snr": 3.1447489880680464,
                "t1_sec": 0.0434,
                "t2_sec": 3599.8137666666667
            },
            "tags": [
                "rejected",
                "noise_overlap"
            ]
        }
}
*/
