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
#include "p_banjoview_cross_correlograms.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <diskreadmda.h>
#include <mlvector.h>
#include "get_sort_indices.h"
#include "textfile.h"
#include <cmath>
using std::ceil;

bool p_banjoview_cross_correlograms(QString firings_path, QString correlograms_out, P_banjoview_cross_correlograms_opts opts)
{
    DiskReadMda firings(firings_path);

    //get the sorted times and labels
    MLVector<double> times(firings.N2());
    MLVector<int> labels(firings.N2());
    {
        MLVector<double> times_unsorted(firings.N2());
        MLVector<int> labels_unsorted(firings.N2());
        for (bigint i = 0; i < firings.N2(); i++) {
            times_unsorted[i] = firings.value(1, i);
            labels_unsorted[i] = firings.value(2, i);
        }
        MLVector<bigint> sort_inds = get_sort_indices(times_unsorted);
        for (bigint i = 0; i < firings.N2(); i++) {
            times[i] = times_unsorted[sort_inds[i]];
            labels[i] = labels_unsorted[sort_inds[i]];
        }
    }

    int K = MLCompute::max(labels);
    if (opts.clusters.isEmpty()) {
        QSet<int> labels_present;
        for (bigint aa = 0; aa < labels.count(); aa++) {
            labels_present.insert(labels[aa]);
        }
        for (int k = 1; k <= K; k++) {
            if (labels_present.contains(k))
                opts.clusters << k;
        }
    }

    QVector<int> k1s, k2s;
    if (opts.mode == Autocorrelograms) {
        for (int ii = 0; ii < opts.clusters.count(); ii++) {
            k1s << opts.clusters[ii];
            k2s << opts.clusters[ii];
        }
    }
    else if (opts.mode == Matrix_of_cross_correlograms) {
        for (int i1 = 0; i1 < opts.clusters.count(); i1++) {
            for (int i2 = 0; i2 < opts.clusters.count(); i2++) {
                k1s << opts.clusters[i1];
                k2s << opts.clusters[i2];
            }
        }
    }

    //set up the histogram bins
    int R = ceil((opts.max_dt_msec) / (opts.bin_size_msec));
    bigint num_bins = 2 * R;
    double tmin = -R * opts.bin_size_msec / 1000 * opts.samplerate;
    double tmax = R * opts.bin_size_msec / 1000 * opts.samplerate;
    double bin_size = opts.bin_size_msec / 1000 * opts.samplerate;

    //initialize the histograms with zeros
    QMap<bigint, QVector<double> > histogram_counts;
    QVector<double> zeros(num_bins, 0);
    for (int ii = 0; ii < k1s.count(); ii++) {
        int k1 = k1s[ii];
        int k2 = k2s[ii];
        int num = (k1 - 1) + K * (k2 - 1);
        histogram_counts[num] = zeros;
    }

    qDebug().noquote() << "k1s:" << k1s;
    qDebug().noquote() << "k2s:" << k1s;

    bigint jjj_last = 0;
    for (bigint ii = 0; ii < times.count(); ii++) {
        double t1 = times[ii];
        int k1 = labels[ii];
        bigint jjj = jjj_last;
        while ((jjj + 1 < times.count()) && (times[jjj] < t1 + tmin))
            jjj++;
        jjj_last = jjj;
        while ((jjj < times.count()) && (t1 + tmin <= times[jjj]) && (times[jjj] < t1 + tmax)) {
            if (jjj != ii) {
                double t2 = times[jjj];
                int k2 = labels[jjj];
                bigint num = (k1 - 1) + K * (k2 - 1);
                if (histogram_counts.contains(num)) {
                    double dt = t2 - t1;
                    int bin_index = (dt - tmin) / bin_size;
                    histogram_counts[num][bin_index]++;
                }
            }
            jjj++;
        }
    }

    QJsonArray correlograms;
    for (int ii = 0; ii < k1s.count(); ii++) {
        QJsonObject correlogram;
        correlogram["k1"] = k1s[ii];
        correlogram["k2"] = k2s[ii];
        correlogram["dt_min_msec"] = tmin * 1000 / opts.samplerate;
        correlogram["bin_size_msec"] = bin_size * 1000 / opts.samplerate;
        int num = (k1s[ii] - 1) + K * (k2s[ii] - 1);
        QJsonArray counts0;
        for (int a = 0; a < num_bins; a++) {
            counts0 << histogram_counts[num][a];
        }
        correlogram["counts"] = counts0;
        correlograms.push_back(correlogram);
    }

    QJsonObject output;
    output["correlograms"] = correlograms;
    QString json_output = QJsonDocument(output).toJson(QJsonDocument::Compact);

    return TextFile::write(correlograms_out, json_output);
}
