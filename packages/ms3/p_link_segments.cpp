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

#include "p_link_segments.h"

#include <mda.h>
#include "mlutil.h"
#include "get_sort_indices.h"

namespace P_link_segments {
QMap<bigint, bigint> get_label_map(const QVector<double>& times, const QVector<bigint>& labels, const QVector<double>& times_prev, const QVector<bigint>& labels_prev);
}

bool p_link_segments(QString firings_path, QString firings_prev_path, QString Kmax_prev_path, QString firings_out_path, QString firings_subset_out_path, QString Kmax_out_path, double t1, double t2, double t1_prev, double t2_prev)
{
    Mda firings(firings_path);
    Mda firings_prev(firings_prev_path);
    Mda Kmax_prev0(Kmax_prev_path);
    bigint Kmax_prev = Kmax_prev0.value(0);

    double t1_intersect = qMax(t1, t1_prev);
    double t2_intersect = qMin(t2, t2_prev);

    QVector<double> times_intersect;
    QVector<bigint> labels_intersect;
    for (bigint i = 0; i < firings.N2(); i++) {
        double t0 = firings.value(1, i);
        if ((t1_intersect <= t0) && (t0 <= t2_intersect)) {
            times_intersect << t0;
            labels_intersect << firings.value(2, i);
        }
    }
    QVector<double> times_intersect_prev;
    QVector<bigint> labels_intersect_prev;
    for (bigint i = 0; i < firings_prev.N2(); i++) {
        double t0 = firings_prev.value(1, i);
        if ((t1_intersect <= t0) && (t0 <= t2_intersect)) {
            times_intersect_prev << t0;
            labels_intersect_prev << firings_prev.value(2, i);
        }
    }

    QMap<bigint, bigint> label_map = P_link_segments::get_label_map(times_intersect, labels_intersect, times_intersect_prev, labels_intersect_prev);
    // Careful! Some labels may not have been included because the intersection may not contain them!
    bigint new_Kmax = Kmax_prev;
    for (bigint i = 0; i < firings_prev.N2(); i++) {
        //important! in case Kmax_prev was not provided
        if (firings_prev.value(2, i) > new_Kmax)
            new_Kmax = firings_prev.value(2, i);
    }

    for (bigint i = 0; i < firings.N2(); i++) {
        bigint label0 = firings.value(2, i);
        if ((!label_map.contains(label0)) || (label_map[label0] == 0)) {
            label_map[label0] = new_Kmax + 1;
            new_Kmax++;
        }
    }

    QVector<bigint> inds_subset;
    QVector<bigint> labels_out;
    for (bigint i = 0; i < firings.N2(); i++) {
        double t0 = firings.value(1, i);
        if (t0 > t2_prev) {
            inds_subset << i;
        }
        bigint label0 = firings.value(2, i);
        labels_out << label_map.value(label0);
    }

    Mda Kmax_out(1, 1);
    Kmax_out.setValue(new_Kmax, 0);
    if (!Kmax_out.write64(Kmax_out_path))
        return false;

    Mda firings_out(firings.N1(), firings.N2());
    for (bigint j = 0; j < firings.N2(); j++) {
        for (bigint r = 0; r < firings.N1(); r++) {
            firings_out.setValue(firings.value(r, j), r, j);
        }
        firings_out.setValue(labels_out[j], 2, j);
    }
    if (!firings_out.write64(firings_out_path))
        return false;

    Mda firings_subset_out(firings.N1(), inds_subset.count());
    for (bigint j = 0; j < inds_subset.count(); j++) {
        for (bigint r = 0; r < firings.N1(); r++) {
            firings_subset_out.setValue(firings.value(r, inds_subset[j]), r, j);
        }
        firings_subset_out.setValue(labels_out[inds_subset[j]], 2, j);
    }
    if (!firings_subset_out.write64(firings_subset_out_path))
        return false;

    return true;
}

namespace P_link_segments {
typedef QString LabelPair;
bigint k1(LabelPair pair) { return pair.split(",").value(0).toLong(); }
bigint k2(LabelPair pair) { return pair.split(",").value(1).toLong(); }

QVector<bigint> condense_labels(const QVector<bigint>& labels, QMap<bigint, bigint>& k_map_out)
{
    QSet<bigint> label_set;
    for (bigint i = 0; i < labels.count(); i++)
        label_set.insert(labels[i]);
    QList<bigint> label_list = label_set.toList();
    qSort(label_list);
    QMap<bigint, bigint> inv_label_map;
    for (bigint i = 0; i < label_list.count(); i++) {
        inv_label_map[label_list[i]] = i;
        k_map_out[i] = label_list[i];
    }
    QVector<bigint> ret(labels.count());
    for (bigint i = 0; i < labels.count(); i++) {
        ret[i] = inv_label_map[labels[i]];
    }
    return ret;
}

void sort_times_labels(QVector<double>& times, QVector<bigint>& labels)
{
    QList<bigint> sort_inds = get_sort_indices_bigint(times);
    QVector<double> times2(times.count());
    QVector<bigint> labels2(times.count());
    for (bigint i = 0; i < times.count(); i++) {
        times2[i] = times[sort_inds[i]];
        labels2[i] = labels[sort_inds[i]];
    }
    times = times2;
    labels = labels2;
}

bool should_link(const Mda& CM, const QVector<bigint>& counts1, const QVector<bigint>& counts2, bigint i1, bigint i2, double link_threshold)
{
    bigint best_i1 = 0;
    for (bigint j1 = 0; j1 < CM.N1(); j1++) {
        if (CM.value(j1, i2) > CM.value(best_i1, i2))
            best_i1 = j1;
    }
    bigint best_i2 = 0;
    for (bigint j2 = 0; j2 < CM.N2(); j2++) {
        if (CM.value(i1, j2) > CM.value(i1, best_i2))
            best_i2 = j2;
    }
    if (best_i1 != i1)
        return false;
    if (best_i2 != i2)
        return false;
    if (CM.value(i1, i2) < link_threshold * counts1[i1])
        return false;
    if (CM.value(i1, i2) < link_threshold * counts2[i2])
        return false;
    return true;
}

QVector<bigint> get_counts(const QVector<bigint>& labels, bigint K)
{
    QVector<bigint> ret(K);
    ret.fill(0);
    for (bigint i = 0; i < labels.count(); i++) {
        ret[labels[i]]++;
    }
    return ret;
}

Mda get_loose_confusion_matrix(const QVector<double>& times1, const QVector<bigint>& labels1, const QVector<double>& times2, const QVector<bigint>& labels2, bigint K1, bigint K2, bigint time_offset_tolerance)
{
    Mda CM(K1, K2);
    QVector<double> times1_sorted = times1, times2_sorted = times2;
    QVector<bigint> labels1_sorted = labels1, labels2_sorted = labels2;
    sort_times_labels(times1_sorted, labels1_sorted);
    sort_times_labels(times2_sorted, labels2_sorted);
    bigint i1 = 0, i2 = 0;
    while ((i1 < times1_sorted.count()) && (i2 < times2_sorted.count())) {
        double t1 = times1_sorted[i1];
        double t2 = times2_sorted[i2];
        if (t2 < t1 - time_offset_tolerance)
            i2++;
        else if (t2 > t1 + time_offset_tolerance) {
            i1++;
            if (i1 < times1_sorted.count()) {
                while ((i2 > 0) && (times2_sorted[i2] >= times1_sorted[i1]))
                    i2--;
            }
        }
        else {
            CM.setValue(CM.value(labels1[i1], labels2[i2]) + 1, labels1[i1], labels2[i2]);
            i2++;
        }
    }
    return CM;
}

QMap<bigint, bigint> get_label_map(const QVector<double>& times, const QVector<bigint>& labels, const QVector<double>& times_prev, const QVector<bigint>& labels_prev)
{
    bigint time_offset_tolerance = 10;
    double link_threshold = 0.7;

    QMap<bigint, bigint> ret;

    QMap<bigint, bigint> k1_map, k2_map;
    QVector<bigint> labels_con = condense_labels(labels, k1_map);
    QVector<bigint> labels_con_prev = condense_labels(labels_prev, k2_map);
    bigint K1 = MLCompute::max(labels_con) + 1;
    bigint K2 = MLCompute::max(labels_con_prev) + 1;

    Mda CM = get_loose_confusion_matrix(times, labels_con, times_prev, labels_con_prev, K1, K2, time_offset_tolerance);
    QVector<bigint> counts1 = get_counts(labels_con, K1);
    QVector<bigint> counts2 = get_counts(labels_con_prev, K2);

    for (bigint i1 = 0; i1 < K1; i1++) {
        for (bigint i2 = 0; i2 < K2; i2++) {
            if (should_link(CM, counts1, counts2, i1, i2, link_threshold)) {
                ret[k1_map[i1]] = k2_map[i2];
            }
        }
    }

    return ret;
}
}
