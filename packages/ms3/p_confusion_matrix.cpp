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

#include "p_confusion_matrix.h"
#include "mlutil.h"
#include "hungarian.h"
#include "get_sort_indices.h"

#include <diskreadmda.h>
#include <cmath>
using std::fabs;

namespace P_confusion_matrix {
struct MFEvent {
    int chan = -1;
    double time = -1;
    int label = -1;
};
struct MFMergeEvent {
    int chan = -1;
    double time = -1;
    int label1 = -1;
    int label2 = -1;
};

void sort_events_by_time(QList<MFEvent>& events);
void sort_events_by_time(QList<MFMergeEvent>& events);
int compute_max_label(const QList<MFEvent>& events);
}

bool p_confusion_matrix(QString firings1, QString firings2, QString confusion_matrix_out, QString matched_firings_out, QString label_map_out, QString firings2_relabeled_out, QString firings2_relabel_map_out, P_confusion_matrix_opts opts)
{
    using namespace P_confusion_matrix;

    if (opts.relabel_firings2) {
        /*
        if (firings2_relabeled_out.isEmpty()) {
            qWarning() << "firings2_relabeled_out is empty even though relabel_firings2 is true";
            return false;
        }
        */
    }
    else {
        if (!firings2_relabeled_out.isEmpty()) {
            qWarning() << "firings2_relabeled_out is not empty even though relabel_firings2 is false";
            return false;
        }
    }

    // Collect the list of events from firings1
    printf("Collecting events...\n");
    QList<MFEvent> events1;
    {
        DiskReadMda F1(firings1);
        for (bigint i = 0; i < F1.N2(); i++) {
            MFEvent evt;
            evt.chan = F1.value(0, i);
            evt.time = F1.value(1, i);
            evt.label = F1.value(2, i);
            events1 << evt;
        }
    }
    // Collect the list of events from firings2
    printf("Collecting events...\n");
    QList<MFEvent> events2;
    {
        DiskReadMda F2(firings2);
        for (bigint i = 0; i < F2.N2(); i++) {
            MFEvent evt;
            evt.chan = F2.value(0, i);
            evt.time = F2.value(1, i);
            evt.label = F2.value(2, i);
            events2 << evt;
        }
    }

    // Sort the events by time
    printf("Sorting events...\n");
    sort_events_by_time(events1);
    sort_events_by_time(events2);

    // Get K1 and K2
    int K1 = compute_max_label(events1);
    int K2 = compute_max_label(events2);

    // Count up every pair that satisfies opts.max_matching_offset -- but don't count redundantly
    printf("Counting all pairs...\n");
    bigint total_counts_12[K1 + 1][K2 + 1];
    {
        for (int k2 = 0; k2 < K2 + 1; k2++) {
            for (int k1 = 0; k1 < K1 + 1; k1++) {
                total_counts_12[k1][k2] = 0;
            }
        }

        bigint i2 = 0;
        for (bigint i1 = 0; i1 < events1.count(); i1++) {
            if (events1[i1].label > 0) {
                int present[K2 + 1];
                for (int k2 = 0; k2 < K2 + 1; k2++) {
                    present[k2] = 0;
                }
                double t1 = events1[i1].time;
                while ((i2 < events2.count()) && (events2[i2].time < t1 - opts.max_matching_offset))
                    i2++;
                bigint old_i2 = i2;
                while ((i2 < events2.count()) && (events2[i2].time <= t1 + opts.max_matching_offset)) {
                    if (events2[i2].label > 0) {
                        present[events2[i2].label] = 1;
                    }
                    i2++;
                }
                i2 = old_i2;
                for (int kk = 1; kk <= K2; kk++) {
                    if (present[kk])
                        total_counts_12[events1[i1].label][kk]++;
                }
            }
        }
    }

    // Count up every pair that satisfies opts.max_matching_offset -- but don't count redundantly
    printf("Counting all pairs...\n");
    bigint total_counts_21[K2 + 1][K1 + 1];
    {
        for (int k1 = 0; k1 < K1 + 1; k1++) {
            for (int k2 = 0; k2 < K2 + 1; k2++) {
                total_counts_21[k2][k1] = 0;
            }
        }

        bigint i1 = 0;
        for (bigint i2 = 0; i2 < events2.count(); i2++) {
            if (events2[i2].label > 0) {
                int present[K1 + 1];
                for (int k1 = 0; k1 < K1 + 1; k1++) {
                    present[k1] = 0;
                }
                double t2 = events2[i2].time;
                while ((i1 < events1.count()) && (events1[i1].time < t2 - opts.max_matching_offset))
                    i1++;
                bigint old_i1 = i1;
                while ((i1 < events1.count()) && (events1[i1].time <= t2 + opts.max_matching_offset)) {
                    if (events1[i1].label > 0) {
                        present[events1[i1].label] = 1;
                    }
                    i1++;
                }
                i1 = old_i1;
                for (int kk = 1; kk <= K1; kk++) {
                    if (present[kk])
                        total_counts_21[events2[i2].label][kk]++;
                }
            }
        }
    }

    bigint event_counts1[K1 + 1];
    for (bigint i = 0; i <= K1; i++)
        event_counts1[i] = 0;
    for (bigint i = 0; i < events1.count(); i++) {
        event_counts1[events1[i].label]++;
    }
    bigint event_counts2[K2 + 1];
    for (bigint i = 0; i <= K2; i++)
        event_counts2[i] = 0;
    for (bigint i = 0; i < events2.count(); i++) {
        event_counts2[events2[i].label]++;
    }

    std::vector<bigint> assignments1(events1.count(), -1);
    std::vector<bigint> assignments2(events2.count(), -1);
    int max_passes = 10;
    for (int pass = 1; pass <= max_passes; pass++) {
        printf("pass %d...\n", pass);
        std::vector<bigint> assignments1_thispass(events1.count(), -1);
        std::vector<bigint> assignments2_thispass(events2.count(), -1);
        { //first fill in the assignments1_thispass
            bigint i2 = 0;
            for (bigint i1 = 0; i1 < events1.count(); i1++) {
                if ((events1[i1].label > 0) && (assignments1[i1] < 0)) { // only consider it if it has a label and hasn't been assigned
                    double t1 = events1[i1].time; // this is the timepoint for event 1

                    //increase i2 until it reaches the lefthand constraint (we assume we are coming from the left)
                    while ((i2 < events2.count()) && (events2[i2].time < t1 - opts.max_matching_offset))
                        i2++;
                    bigint old_i2 = i2; //save this index for later so we can return to this spot for the next event

                    double best_match_score = 0;
                    double abs_offset_of_best_match_score = opts.max_matching_offset + 1;
                    bigint best_i2 = -1;
                    //move through the events in firings2 until we pass the righthand constraint
                    while ((i2 < events2.count()) && (events2[i2].time <= t1 + opts.max_matching_offset)) {
                        if ((events2[i2].label > 0) && (assignments2[i2] < 0)) { //only consider it if it has a label and unassigned
                            double time0 = events2[i2].time;
                            bigint numer12 = total_counts_12[events1[i1].label][events2[i2].label];
                            bigint numer21 = total_counts_21[events2[i2].label][events1[i1].label];
                            bigint denom12 = event_counts1[events1[i1].label];
                            bigint denom21 = event_counts2[events2[i2].label];
                            if (denom12 == 0)
                                denom12 = 1; //don't divide by zero
                            if (denom21 == 0)
                                denom21 = 1;
                            double match_score = qMin(numer12 * 1.0 / denom12, numer21 * 1.0 / denom21);

                            if (match_score >= best_match_score) {
                                double abs_offset = fabs(time0 - t1);
                                //in the case of a tie, use the one that is closer in offset.
                                if ((match_score > best_match_score) || ((match_score == best_match_score) && (abs_offset < abs_offset_of_best_match_score))) {
                                    best_match_score = match_score;
                                    best_i2 = i2;
                                    abs_offset_of_best_match_score = abs_offset;
                                }
                            }
                        }
                        i2++; // go to the next one
                    }
                    if (best_i2 >= 0) {
                        assignments1_thispass[i1] = best_i2;
                    }
                    i2 = old_i2; // go back so we are ready to handle the next event
                }
            }
        }
        { //next fill in the assignments2_thispass
            bigint i1 = 0;
            for (bigint i2 = 0; i2 < events2.count(); i2++) {
                if ((events2[i2].label > 0) && (assignments2[i2] < 0)) { // only consider it if it has a label and hasn't been assigned
                    double t2 = events2[i2].time; // this is the timepoint for event 2

                    //increase i1 until it reaches the lefthand constraint (we assume we are coming from the left)
                    while ((i1 < events1.count()) && (events1[i1].time < t2 - opts.max_matching_offset))
                        i1++;
                    bigint old_i1 = i1; //save this index for later so we can return to this spot for the next event

                    double best_match_score = 0;
                    double abs_offset_of_best_match_score = opts.max_matching_offset + 1;
                    bigint best_i1 = -1;
                    //move through the events in firings2 until we pass the righthand constraint
                    while ((i1 < events1.count()) && (events1[i1].time <= t2 + opts.max_matching_offset)) {
                        if ((events1[i1].label > 0) && (assignments1[i1] < 0)) { //only consider it if it has a label and unassigned
                            double time0 = events1[i1].time;
                            bigint numer12 = total_counts_12[events1[i1].label][events2[i2].label];
                            bigint numer21 = total_counts_21[events2[i2].label][events1[i1].label];
                            bigint denom12 = event_counts1[events1[i1].label];
                            bigint denom21 = event_counts2[events2[i2].label];
                            if (denom12 == 0)
                                denom12 = 1; //don't divide by zero
                            if (denom21 == 0)
                                denom21 = 1;
                            double match_score = qMin(numer12 * 1.0 / denom12, numer21 * 1.0 / denom21);

                            if (match_score >= best_match_score) {
                                double abs_offset = fabs(time0 - t2);
                                //in the case of a tie, use the one that is closer in offset.
                                if ((match_score > best_match_score) || ((match_score == best_match_score) && (abs_offset < abs_offset_of_best_match_score))) {
                                    best_match_score = match_score;
                                    best_i1 = i1;
                                    abs_offset_of_best_match_score = abs_offset;
                                }
                            }
                        }
                        i1++; // go to the next one
                    }
                    if (best_i1 >= 0) {
                        assignments2_thispass[i2] = best_i1;
                    }
                    i1 = old_i1; // go back so we are ready to handle the next event
                }
            }
        }
        //use only those where assignments1_thispass agrees with assignments2_thispass
        bool something_changed = false;
        for (bigint i1 = 0; i1 < events1.count(); i1++) {
            if (assignments1_thispass[i1] >= 0) {
                if (assignments2_thispass[assignments1_thispass[i1]] == i1) {
                    assignments1[i1] = assignments1_thispass[i1];
                    assignments2[assignments1_thispass[i1]] = i1;
                    something_changed = true;
                }
            }
        }
        if (!something_changed)
            break;
    }

    printf("Creating list of merged events...\n");
    // Create the list of matched events
    QList<MFMergeEvent> events3;
    for (bigint i1 = 0; i1 < events1.count(); i1++) {
        if (assignments1[i1] >= 0) {
            bigint i2 = assignments1[i1];
            MFMergeEvent E;
            E.chan = events1[i1].chan;
            E.time = events1[i1].time;
            E.label1 = events1[i1].label;
            E.label2 = events2[i2].label;
            events3 << E;
            events1[i1].time = -1;
            events1[i1].label = -1; //set to -1 , see below
            events2[i2].time = -1;
            events2[i2].label = -1;
        }
    }

    printf("Handling missing events...\n");
    //Now take care of all the events that did not get handled earlier ... they are unclassified
    for (bigint i1 = 0; i1 < events1.count(); i1++) {
        if (events1[i1].label > 0) { //this is how we check to see if it was not handled
            MFMergeEvent E;
            E.chan = events1[i1].chan;
            E.time = events1[i1].time;
            E.label1 = events1[i1].label;
            E.label2 = 0;
            events3 << E;
        }
    }
    for (bigint i2 = 0; i2 < events2.count(); i2++) {
        if (events2[i2].label > 0) {
            MFMergeEvent E;
            E.chan = events2[i2].chan;
            E.time = events2[i2].time;
            E.label1 = 0;
            E.label2 = events2[i2].label;
            events3 << E;
        }
    }

    // Sort the matched events by time
    printf("Sorting events...\n");
    sort_events_by_time(events3);
    bigint L = events3.count();

    // Create the confusion matrix
    printf("Writing confusion matrix...\n");
    Mda confusion_matrix(K1 + 1, K2 + 1);
    {
        for (bigint i = 0; i < L; i++) {
            int a = events3[i].label1 - 1;
            int b = events3[i].label2 - 1;
            if (a < 0)
                a = K1;
            if (b < 0)
                b = K2;
            confusion_matrix.setValue(confusion_matrix.value(a, b) + 1, a, b);
        }
    }
    if (!confusion_matrix_out.isEmpty()) {
        if (!confusion_matrix.write64(confusion_matrix_out))
            return false;
    }

    // Get the optimal label map
    Mda label_map(K1, 1);
    if ((!label_map_out.isEmpty()) || (opts.relabel_firings2)) {
        printf("Writing label_map_out...\n");
        if (K1 > 0) {
            int assignment[K1];
            Mda matrix(K1, K2);
            for (int i = 0; i < K1; i++) {
                for (int j = 0; j < K2; j++) {
                    matrix.setValue((-1) * confusion_matrix.value(i, j), i, j);
                }
            }
            double cost;
            hungarian(assignment, &cost, matrix.dataPtr(), K1, K2);
            for (int i = 0; i < K1; i++) {
                label_map.setValue(assignment[i] + 1, i);
            }
        }
        if (!label_map_out.isEmpty()) {
            if (!label_map.write64(label_map_out))
                return false;
        }
    }

    if (opts.relabel_firings2) {
        printf("Relabeling firings2...\n");
        Mda label_map_inv(K2, 1);
        for (int j = 1; j <= K1; j++) {
            if (label_map.value(j - 1) > 0) {
                label_map_inv.setValue(j, label_map.value(j - 1) - 1);
            }
        }
        int kk = K1 + 1;
        for (int j = 0; j < K2; j++) {
            if (!label_map_inv.value(j)) {
                label_map_inv.setValue(kk, j);
                kk++;
            }
        }
        for (int j = 0; j < K1; j++) {
            label_map.setValue(j + 1, j);
        }
        if (!label_map_out.isEmpty()) {
            if (!label_map.write64(label_map_out))
                return false;
        }
        Mda confusion_matrix2(K1 + 1, K2 + 1);
        for (int k1 = 1; k1 <= K1 + 1; k1++) {
            for (int k2 = 1; k2 <= K2 + 1; k2++) {
                int k2b = k2;
                if (k2 <= K2)
                    k2b = label_map_inv.value(k2 - 1);
                confusion_matrix2.setValue(confusion_matrix.value(k1 - 1, k2 - 1), k1 - 1, k2b - 1);
            }
        }
        if (!confusion_matrix_out.isEmpty()) {
            printf("Rewriting confusion matrix...\n");
            if (!confusion_matrix2.write64(confusion_matrix_out)) {
                qWarning() << "Error rewriting confusion matrix output";
                return false;
            }
        }
        if (!firings2_relabeled_out.isEmpty()) {
            printf("Writing firings_relabeled_out...\n");
            Mda firings2_relabeled(firings2);
            for (bigint i = 0; i < firings2_relabeled.N2(); i++) {
                int k = firings2_relabeled.value(2, i);
                if ((k >= 1) && (k <= K2)) {
                    firings2_relabeled.setValue(label_map_inv.value(k - 1), 2, i);
                }
            }
            if (!firings2_relabeled.write64(firings2_relabeled_out))
                return false;
        }
        if (!firings2_relabel_map_out.isEmpty()) {
            printf("Writing firings2_relabel_map...\n");
            if (!label_map_inv.write64(firings2_relabel_map_out))
                return false;
        }
        for (bigint i = 0; i < L; i++) {
            MFMergeEvent* E = &events3[i];
            int k = E->label2;
            if ((k >= 1) && (k <= K2)) {
                E->label2 = label_map_inv.value(k - 1);
            }
        }
    }

    if (!matched_firings_out.isEmpty()) {
        printf("Writing matched_firings...\n");
        // Create the matched firings file
        Mda matched_firings(4, L);
        {
            for (bigint i = 0; i < L; i++) {
                MFMergeEvent* E = &events3[i];
                matched_firings.setValue(E->chan, 0, i);
                matched_firings.setValue(E->time, 1, i);
                matched_firings.setValue(E->label1, 2, i);
                matched_firings.setValue(E->label2, 3, i);
            }
        }
        if (!matched_firings.write64(matched_firings_out))
            return false;
    }

    return true;
}

namespace P_confusion_matrix {

void sort_events_by_time(QList<MFEvent>& events)
{
    QVector<double> times(events.count());
    for (bigint i = 0; i < events.count(); i++) {
        times[i] = events[i].time;
    }
    QList<bigint> inds = get_sort_indices_bigint(times);
    QList<MFEvent> ret;
    for (bigint i = 0; i < inds.count(); i++) {
        ret << events[inds[i]];
    }
    events = ret;
}

void sort_events_by_time(QList<MFMergeEvent>& events)
{
    QVector<double> times(events.count());
    for (bigint i = 0; i < events.count(); i++) {
        times[i] = events[i].time;
    }
    QList<bigint> inds = get_sort_indices_bigint(times);
    QList<MFMergeEvent> ret;
    for (bigint i = 0; i < inds.count(); i++) {
        ret << events[inds[i]];
    }
    events = ret;
}

int compute_max_label(const QList<MFEvent>& events)
{
    int ret = 0;
    for (bigint i = 0; i < events.count(); i++) {
        if (events[i].label > ret)
            ret = events[i].label;
    }
    return ret;
}
}
