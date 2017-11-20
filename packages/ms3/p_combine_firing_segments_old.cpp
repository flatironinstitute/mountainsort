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
#include "fitstagecomputer.h"
#include "p_combine_firing_segments.h"

#include <diskreadmda.h>
#include <mda32.h>
#include "compute_templates_0.h"
#include "get_sort_indices.h"
#include "omp.h"

namespace P_combine_firing_segments {

struct SegmentInfo {
    int K;
    QVector<double> times;
    QVector<int> labels;
    QMap<int, int> label_map_with_previous;
    QMap<int, int> label_map_global;
    QVector<bool> events_to_delete;
    Mda32 templates;
    QString firings;
    double tmin, tmax;
};

void compute_segment_info(SegmentInfo& info, QString timeseries_path, QString firings_path, P_combine_firing_segments_opts opts);
int find_matching_cluster(QVector<bigint>& event_inds_to_delete_out, SegmentInfo* Sprev, SegmentInfo* S, int k2, P_combine_firing_segments_opts opts);

struct TimeChunkInfo {
    bigint t_padding; //t1 corresponds to index t_padding (because we have padding/overlap)
    bigint t1; //starting timepoint (corresponding to index t_padding)
    bigint size; //number of timepoints (excluding padding on left and right)
};

QList<TimeChunkInfo> get_time_chunk_infos(bigint M, bigint N, int num_simultaneous, bigint min_num_chunks);
QVector<double> get_subarray(const QVector<double>& X, const QVector<bigint>& inds);
QVector<int> get_subarray(const QVector<int>& X, const QVector<bigint>& inds);
QVector<bigint> get_subarray(const QVector<bigint>& X, const QVector<bigint>& inds);
}

using namespace P_combine_firing_segments;

bool p_combine_firing_segments(QString timeseries_path, QStringList firings_list, QString firings_out, P_combine_firing_segments_opts opts)
{
    qDebug().noquote() << "Computing segment info";
    QList<SegmentInfo> segments;
#pragma omp for
    for (int ii = 0; ii < firings_list.count(); ii++) {
        SegmentInfo info0;
        compute_segment_info(info0, timeseries_path, firings_list[ii], opts);
#pragma omp critical(lock1_combine_firing_segments)
        {
            segments << info0;
        }
    }

    //define label map with previous and events to delete
    qDebug().noquote() << "Defining label maps with previous";
#pragma omp for
    for (int ii = 1; ii < segments.count(); ii++) {
        SegmentInfo* Sprev, *S;
#pragma omp critical
        {
            Sprev = &segments[ii - 1];
            S = &segments[ii];
        }
        for (int k2 = 1; k2 <= S->K; k2++) {
            QVector<bigint> event_inds_to_delete;
            int k1 = find_matching_cluster(event_inds_to_delete, Sprev, S, k2, opts);
            S->label_map_with_previous[k2] = k1;
            for (bigint jj = 0; jj < event_inds_to_delete.count(); jj++) {
                S->events_to_delete[event_inds_to_delete[jj]] = true;
            }
        }
    }

    //define label map global
    qDebug().noquote() << "Defining global label maps";
    int kk = 0;
    if (segments.count() > 0) {
        SegmentInfo* S0 = &segments[0];
        for (int k = 1; k <= S0->K; k++) {
            S0->label_map_global[k] = k;
        }
        kk = S0->K + 1;
    }
    for (int ii = 1; ii < segments.count(); ii++) {
        SegmentInfo* Sprev = &segments[ii - 1];
        SegmentInfo* S = &segments[ii];
        for (int k2 = 1; k2 <= S->K; k2++) {
            int kmatch = S->label_map_with_previous.value(k2, 0);
            if (kmatch) {
                kmatch = Sprev->label_map_global.value(kmatch, 0);
            }
            if (!kmatch) {
                kmatch = kk;
                kk++;
            }
            S->label_map_global[k2] = kmatch;
        }
    }

    //define new labels
    qDebug().noquote() << "Defining new labels";
    QVector<double> new_times;
    QVector<int> new_labels;
    QVector<int> segment_numbers;
    QVector<bigint> event_indices;
    for (int ii = 0; ii < segments.count(); ii++) {
        SegmentInfo* S = &segments[ii];
        for (bigint jj = 0; jj < S->labels.count(); jj++) {
            if (!S->events_to_delete[jj]) {
                int k0 = S->labels[jj];
                int k1 = S->label_map_global.value(k0, 0);
                if (k1) {
                    new_times << S->times[jj];
                    new_labels << k1;
                    segment_numbers << ii;
                    event_indices << jj;
                }
            }
        }
    }

    DiskReadMda32 X(timeseries_path);
    Mda32 new_templates = compute_templates_in_parallel(X, new_times, new_labels, opts.clip_size);

    qDebug().noquote() << "Fit stage";
    bigint t_start = 0;
    int M = X.N1();
    FitStageComputer FSC;
    FSC.setTemplates(new_templates);
    FSC.setTimesLabels(new_times, new_labels);
    int num_threads = omp_get_max_threads();
    QList<TimeChunkInfo> time_chunk_infos = get_time_chunk_infos(M, X.N2(), num_threads, 1);
    for (bigint i = 0; i < time_chunk_infos.count(); i += num_threads) {
        QList<TimeChunkInfo> infos = time_chunk_infos.mid(i, num_threads);
        QList<Mda32> time_chunks;
        double bytes0 = 0;
        for (int j = 0; j < infos.count(); j++) {
            Mda32 time_chunk;
            X.readChunk(time_chunk, 0, t_start + infos[j].t1 - infos[j].t_padding, M, infos[j].size + 2 * infos[j].t_padding);
            time_chunks << time_chunk;
            bytes0 += time_chunk.totalSize() * sizeof(float);
        }
#pragma omp parallel for num_threads(num_threads)
        for (int j = 0; j < infos.count(); j++) {
            FSC.processTimeChunk(infos[j].t1, time_chunks[j], infos[j].t_padding, infos[j].t_padding);
        }
    }
    FSC.finalize();
    QVector<bigint> event_inds_to_use = FSC.eventIndicesToUse();
    qDebug().noquote() << QString("Using %1 of %2 events (%3%) after fit stage").arg(event_inds_to_use.count()).arg(new_times.count()).arg(event_inds_to_use.count() * 1.0 / new_times.count() * 100);
    new_times = get_subarray(new_times, event_inds_to_use);
    new_labels = get_subarray(new_labels, event_inds_to_use);
    segment_numbers = get_subarray(segment_numbers, event_inds_to_use);
    event_indices = get_subarray(event_indices, event_inds_to_use);

    qDebug().noquote() << "Creating new firings file";
    bigint R = 0;
    if (segments.count() > 0) {
        SegmentInfo* S0 = &segments[0];
        R = DiskReadMda(S0->firings).N1();
    }
    QList<DiskReadMda> all_firings;
    for (int ii = 0; ii < segments.count(); ii++) {
        all_firings << DiskReadMda(segments[ii].firings);
    }
    bigint L = new_labels.count();
    Mda FF(R, L);
    for (bigint jj = 0; jj < L; jj++) {
        int ss = segment_numbers[jj];
        bigint aa = event_indices[jj];
        for (int r = 0; r < R; r++) {
            FF.setValue(all_firings[ss].value(r, aa), r, jj);
        }
        FF.setValue(new_labels[jj], 2, jj);
    }

    qDebug().noquote() << "Writing firings file";
    return FF.write64(firings_out);
}

namespace P_combine_firing_segments {

void compute_segment_info(SegmentInfo& info, QString timeseries_path, QString firings_path, P_combine_firing_segments_opts opts)
{
    DiskReadMda32 X(timeseries_path);
    DiskReadMda firings(firings_path);
    for (bigint i = 0; i < firings.N2(); i++) {
        info.times << firings.value(1, i);
        info.labels << firings.value(2, i);
    }
    info.events_to_delete.resize(info.labels.count());
    info.events_to_delete.fill(false);
    info.K = MLCompute::max(info.labels);
    info.templates = compute_templates_0(X, info.times, info.labels, opts.clip_size);
    info.firings = firings_path;
    info.tmin = MLCompute::min(info.times);
    info.tmax = MLCompute::max(info.times);
}

double compute_noncentered_correlation(int N, const float* X1, const float* X2)
{
    double S12 = 0, S11 = 0, S22 = 0;
    for (int i = 0; i < N; i++) {
        S12 += X1[i] * X2[i];
        S11 += X1[i] * X1[i];
        S22 += X2[i] * X2[i];
    }
    if ((!S11) || (!S22))
        return 0;
    return S12 / (sqrt(S11) * sqrt(S22));
}

Mda32 time_shift_template(const Mda32& template0, int dt)
{
    int M = template0.N1();
    int T = template0.N2();
    Mda32 ret(M, T);
    for (int t = 0; t < T; t++) {
        for (int m = 0; m < M; m++) {
            ret.setValue(template0.value(m, t + dt), m, t);
        }
    }
    return ret;
}

void compute_sliding_correlation_between_templates(double& best_corr_out, int& best_offset2_out, const Mda32& template1, const Mda32& template2, P_combine_firing_segments_opts opts)
{
    int M = template1.N1();
    int T = template1.N2();
    best_corr_out = 0;
    best_offset2_out = 0;
    for (int dt = -opts.offset_search_radius; dt <= opts.offset_search_radius; dt++) {
        Mda32 template2b = time_shift_template(template2, dt);
        double val = compute_noncentered_correlation(M * T, template1.constDataPtr(), template2b.constDataPtr());
        if (val > best_corr_out) {
            best_corr_out = val;
            best_offset2_out = dt;
        }
    }
}

QVector<bigint> get_match_indices(const QVector<double>& times1, const QVector<double>& times2, int best_offset2, P_combine_firing_segments_opts opts)
{
    QVector<bigint> matching_indices;
    QList<bigint> sort_inds1 = get_sort_indices_bigint(times1);
    QList<bigint> sort_inds2 = get_sort_indices_bigint(times2);
    int a1 = best_offset2 - opts.match_tolerance;
    int a2 = best_offset2 + opts.match_tolerance;
    bigint i1 = 0;
    for (bigint i2 = 0; i2 < sort_inds2.count(); i2++) {
        double t2 = times2[sort_inds2[i2]];
        while ((i1 - 1 >= 0) && (times1[sort_inds1[i1 - 1]] >= t2 + a1))
            i1--;
        while ((i1 + 1 < sort_inds1.count()) && (times1[sort_inds1[i1 + 1]] <= t2 + a2))
            i1++;
        if ((times1[sort_inds1[i1]] >= t2 + a1) && (times1[sort_inds1[i1]] <= t2 + a2)) {
            matching_indices << sort_inds2[i2];
        }
    }
    return matching_indices;
}

int find_matching_cluster(QVector<bigint>& event_inds_to_delete_out, SegmentInfo* Sprev, SegmentInfo* S, int k2, P_combine_firing_segments_opts opts)
{
    QVector<double> times2;
    QVector<bigint> inds2;
    for (bigint i = 0; i < S->times.count(); i++) {
        if (S->times[i] <= Sprev->tmax) {
            if (S->labels[i] == k2) {
                times2 << S->times[i];
                inds2 << i;
            }
        }
    }
    if (times2.count() == 0)
        return 0;
    int M = S->templates.N1();
    int T = S->templates.N2();
    Mda32 template2;
    S->templates.getChunk(template2, 0, 0, k2 - 1, M, T, 1);
    double best_match_score = 0;
    int best_match_k1 = 0;
    QMap<int, QVector<bigint> > matching_indices;
    for (int k1 = 1; k1 <= Sprev->K; k1++) {
        Mda32 template1;
        Sprev->templates.getChunk(template1, 0, 0, k1 - 1, M, T, 1);
        double corr = 0;
        int best_offset2 = 0;
        compute_sliding_correlation_between_templates(corr, best_offset2, template1, template2, opts);
        if (corr > opts.template_correlation_threshold) {
            QVector<double> times1;
            for (bigint i = 0; i < Sprev->times.count(); i++) {
                if (Sprev->times[i] >= S->tmin) {
                    if (Sprev->labels[i] == k1)
                        times1 << Sprev->times[i];
                }
            }
            QVector<bigint> match_indices = get_match_indices(times1, times2, best_offset2, opts);
            double numer = match_indices.count();
            double denom = times1.count() + times2.count() - match_indices.count();
            double match_score = numer / denom;
            if (match_score >= opts.match_score_threshold) {
                if (match_score > best_match_score) {
                    best_match_score = match_score;
                    best_match_k1 = k1;
                }
                for (bigint bb = 0; bb < match_indices.count(); bb++) {
                    matching_indices[k1] << inds2[match_indices[bb]];
                }
            }
        }
    }
    if (best_match_score >= opts.match_score_threshold) {
        event_inds_to_delete_out = matching_indices.value(best_match_k1);
        return best_match_k1;
    }
    else
        return 0;
}

QList<TimeChunkInfo> get_time_chunk_infos(bigint M, bigint N, int num_simultaneous, bigint min_num_chunks)
{
    bigint RAM_available_for_chunks_bytes = 1e9;
    bigint chunk_size = RAM_available_for_chunks_bytes / (M * num_simultaneous * 4);
    if (chunk_size > N)
        chunk_size = N;
    if (N < chunk_size * min_num_chunks) {
        chunk_size = N / min_num_chunks;
    }
    if (chunk_size < 20000)
        chunk_size = 20000;
    bigint chunk_overlap_size = 1000;

    // Prepare the information on the time chunks
    QList<TimeChunkInfo> time_chunk_infos;
    for (bigint t = 0; t < N; t += chunk_size) {
        TimeChunkInfo info;
        info.t1 = t;
        info.t_padding = chunk_overlap_size;
        info.size = chunk_size;
        if (t + info.size > N)
            info.size = N - t;
        time_chunk_infos << info;
    }

    return time_chunk_infos;
}

QVector<double> get_subarray(const QVector<double>& X, const QVector<bigint>& inds)
{
    QVector<double> ret(inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        ret[i] = X[inds[i]];
    }
    return ret;
}

QVector<int> get_subarray(const QVector<int>& X, const QVector<bigint>& inds)
{
    QVector<int> ret(inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        ret[i] = X[inds[i]];
    }
    return ret;
}

QVector<bigint> get_subarray(const QVector<bigint>& X, const QVector<bigint>& inds)
{
    QVector<bigint> ret(inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        ret[i] = X[inds[i]];
    }
    return ret;
}
}
