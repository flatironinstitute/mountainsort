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
#include "p_combine_firing_segments.h"

#include <QTime>
#include <diskreadmda.h>
#include <mda32.h>
#include "compute_templates_0.h"
#include "get_sort_indices.h"
#include "omp.h"
#include "pca.h"
#include "kdtree.h"
#include <cmath>
using std::sqrt;

namespace P_combine_firing_segments {

struct SegmentInfo {
    int K;
    QVector<double> times;
    QVector<int> labels;
    Mda32 templates;
    QString firings;

    QMap<int, int> label_map_with_previous;
    QMap<int, int> label_map_global;
    QMap<int, double> time_offset_map_with_previous;
    QMap<int, double> time_offset_map_global;
};

void compute_segment_info(SegmentInfo& info, QString timeseries_path, QString firings_path, P_combine_firing_segments_opts opts);
void get_ending_times_labels(QVector<double>& times_out, QVector<int>& labels_out, const QVector<double>& times, const QVector<int>& labels, P_combine_firing_segments_opts opts, bool use_beginning = false);
void get_beginning_times_labels(QVector<double>& times_out, QVector<int>& labels_out, const QVector<double>& times, const QVector<int>& labels, P_combine_firing_segments_opts opts);
Mda32 extract_clips(const DiskReadMda32& X, const QVector<double>& times, int clip_size);
Mda32 compute_clips_features(const Mda32& clips, int num_features);
QVector<bigint> find_nearest_neighbors(Mda32& FF1, Mda32& FF2, int num_features);
Mda32 compute_templates_from_clips_and_labels(const Mda32& clips, const QVector<int>& labels, int K);
Mda32 align_clips(const Mda32& clips1, const QVector<int>& labels1, const Mda32& templates1, const Mda32& template_k2, int offset_search_radius);
void compute_sliding_correlation_between_templates(double& best_corr_out, int& best_offset2_out, const Mda32& template1, const Mda32& template2, int offset_search_radius);

void write_mda(const QVector<int>& X, QString fname); //for debug
void write_mda(const QVector<bigint>& X, QString fname); //for debug
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

    //define label map with previous
    qDebug().noquote() << "Defining label maps with previous";
#pragma omp parallel for
    for (int ii = 1; ii < segments.count(); ii++) {
        SegmentInfo* Sprev, *S;
#pragma omp critical
        {
            Sprev = &segments[ii - 1];
            S = &segments[ii];
        }
        DiskReadMda32 X(timeseries_path);
        QVector<double> times1, times2;
        QVector<int> labels1, labels2;
        get_ending_times_labels(times1, labels1, Sprev->times, Sprev->labels, opts);
        get_beginning_times_labels(times2, labels2, S->times, S->labels, opts);
        Mda32 clips1 = extract_clips(X, times1, opts.clip_size);
        Mda32 clips2 = extract_clips(X, times2, opts.clip_size);
        int M = clips1.N1();
        int T = clips1.N2();
        Mda32 templates1 = compute_templates_from_clips_and_labels(clips1, labels1, Sprev->K);
        Mda32 templates2 = compute_templates_from_clips_and_labels(clips2, labels2, S->K);
        Mda neighbor_counts(Sprev->K, S->K);
        Mda counts(S->K, 1);
        for (int k2 = 1; k2 <= S->K; k2++) {
            QVector<bigint> inds_k2;
            for (bigint a = 0; a < labels2.count(); a++) {
                if (labels2[a] == k2)
                    inds_k2 << a;
            }
            Mda32 template_k2;
            templates2.getChunk(template_k2, 0, 0, k2 - 1, M, T, 1);
            Mda32 clips1_aligned = align_clips(clips1, labels1, templates1, template_k2, opts.offset_search_radius);
            Mda32 clips_k2(M, T, inds_k2.count());
            for (bigint a = 0; a < inds_k2.count(); a++) {
                Mda32 tmp;
                clips2.getChunk(tmp, 0, 0, inds_k2[a], M, T, 1);
                clips_k2.setChunk(tmp, 0, 0, a);
            }
            Mda32 clips1_aligned_reshaped(M * T, clips1_aligned.N3());
            memcpy(clips1_aligned_reshaped.dataPtr(), clips1_aligned.dataPtr(), sizeof(float) * clips1_aligned.totalSize());
            Mda32 clips_k2_reshaped(M * T, clips_k2.N3());
            memcpy(clips_k2_reshaped.dataPtr(), clips_k2.dataPtr(), sizeof(float) * clips_k2.totalSize());
            QTime timer;
            timer.start();
            int num_features = 10;
            QVector<bigint> neighbor_inds = find_nearest_neighbors(clips1_aligned_reshaped, clips_k2_reshaped, num_features);
            for (bigint a = 0; a < neighbor_inds.count(); a++) {
                int k1 = labels1[neighbor_inds[a]];
                if (k1 > 0) {
                    neighbor_counts.set(neighbor_counts.get(k1 - 1, k2 - 1) + 1, k1 - 1, k2 - 1);
                }
            }
            counts.set(inds_k2.count(), k2 - 1);
        }
        Mda match_scores(Sprev->K, S->K);
        for (int k2 = 1; k2 <= S->K; k2++) {
            for (int k1 = 1; k1 <= Sprev->K; k1++) {
                double numer = neighbor_counts.get(k1 - 1, k2 - 1);
                double denom = counts.get(k2 - 1);
                if (denom)
                    match_scores.set(numer / denom, k1 - 1, k2 - 1);
            }
        }

        /*
        Mda32 clips12(M,T,clips1.N3()+clips2.N3());
        clips12.setChunk(clips1,0,0,0);
        clips12.setChunk(clips2,0,0,clips1.N3());
        int num_features=10;
        Mda32 FF12=compute_clips_features(clips12,num_features);
        Mda32 FF1,FF2;
        FF12.getChunk(FF1,0,0,num_features,clips1.N3());
        FF12.getChunk(FF2,0,clips1.N3(),num_features,clips2.N3());
        QVector<bigint> neighbors12=find_nearest_neighbors(FF1,FF2);
        QVector<bigint> neighbors21=find_nearest_neighbors(FF2,FF1);
        QVector<bigint> counts1(Sprev->K);
        QVector<bigint> counts2(S->K);
        counts1.fill(0);
        counts2.fill(0);
        Mda match_counts(Sprev->K,S->K);
        for (bigint i=0; i<times1.count(); i++) {
            int k1=labels1[i];
            int k2=labels2[neighbors21[i]];
            if ((k1>0)&&(k2>0)) {
                match_counts.set(match_counts.get(k1-1,k2-1)+1,k1-1,k2-1);
                counts1[k1-1]++;
            }
        }
        for (bigint i=0; i<times2.count(); i++) {
            int k1=labels1[neighbors12[i]];
            int k2=labels2[i];
            if ((k1>0)&&(k2>0)) {
                match_counts.set(match_counts.get(k1-1,k2-1)+1,k1-1,k2-1);
                counts2[k2-1]++;
            }
        }
        */
        /*
        clips1.write32(QString("/home/magland/tmp/clips1_%1.mda").arg(ii));
        clips2.write32(QString("/home/magland/tmp/clips2_%1.mda").arg(ii));
        write_mda(labels1,QString("/home/magland/tmp/labels1_%1.mda").arg(ii));
        write_mda(labels2,QString("/home/magland/tmp/labels2_%1.mda").arg(ii));
        match_counts.write64(QString("/home/magland/tmp/match_counts_%1.mda").arg(ii));
        FF1.write64(QString("/home/magland/tmp/FF1_%1.mda").arg(ii));
        FF2.write64(QString("/home/magland/tmp/FF2_%1.mda").arg(ii));
        write_mda(neighbors12,QString("/home/magland/tmp/neighbors_12_%1.mda").arg(ii));
        write_mda(neighbors21,QString("/home/magland/tmp/neighbors_21_%1.mda").arg(ii));
        */
        /*
        Mda match_scores(Sprev->K,S->K);
        for (int k2=1; k2<=S->K; k2++) {
            for (int k1=1; k1<=Sprev->K; k1++) {
                double numer=2*match_counts.get(k1-1,k2-1);
                double denom=counts1[k1-1]+counts2[k2-1];
                if (denom) match_scores.set(numer/denom,k1-1,k2-1);
            }
        }
        */

        int num_matches = 0;
        while (true) {
            int best_k1 = 0, best_k2 = 0;
            double best_score = 0;
            for (int k2 = 1; k2 <= S->K; k2++) {
                for (int k1 = 1; k1 <= Sprev->K; k1++) {
                    double score0 = match_scores.value(k1 - 1, k2 - 1);
                    if (score0 > best_score) {
                        best_score = score0;
                        best_k1 = k1;
                        best_k2 = k2;
                    }
                }
            }
            if (best_score > opts.match_score_threshold) {
                for (int k1 = 1; k1 <= Sprev->K; k1++) {
                    match_scores.setValue(0, k1 - 1, best_k2 - 1);
                }
                for (int k2 = 1; k2 <= S->K; k2++) {
                    match_scores.setValue(0, best_k1 - 1, k2 - 1);
                }
                qDebug().noquote() << QString("Matching %1 to %2 in segment %3 (score=%4)").arg(best_k1).arg(best_k2).arg(ii).arg(best_score);
                num_matches++;
                S->label_map_with_previous[best_k2] = best_k1;
                Mda32 template1, template2;
                templates1.getChunk(template1, 0, 0, best_k1 - 1, M, T, 1);
                templates2.getChunk(template2, 0, 0, best_k2 - 1, M, T, 1);
                double corr;
                int offset2;
                compute_sliding_correlation_between_templates(corr, offset2, template1, template2, opts.offset_search_radius);
                S->time_offset_map_with_previous[best_k2] = offset2;
            }
            else
                break;
        }
        qDebug().noquote() << QString("Matched %1 of %2 clusters in segment %3").arg(num_matches).arg(S->K).arg(ii);
    }

    //define label map global
    qDebug().noquote() << "Defining global label maps";
    int kk = 0;
    if (segments.count() > 0) {
        SegmentInfo* S0 = &segments[0];
        for (int k = 1; k <= S0->K; k++) {
            S0->label_map_global[k] = k;
            S0->time_offset_map_global[k] = 0;
        }
        kk = S0->K + 1;
    }
    for (int ii = 1; ii < segments.count(); ii++) {
        SegmentInfo* Sprev = &segments[ii - 1];
        SegmentInfo* S = &segments[ii];
        for (int k2 = 1; k2 <= S->K; k2++) {
            int kmatch = S->label_map_with_previous.value(k2, 0);
            double offset0 = S->time_offset_map_with_previous.value(k2, 0);
            if (kmatch) {
                kmatch = Sprev->label_map_global.value(kmatch, 0);
                offset0 = offset0 + Sprev->time_offset_map_global.value(kmatch, 0);
            }
            if (!kmatch) {
                kmatch = kk;
                kk++;
            }
            S->label_map_global[k2] = kmatch;
            S->time_offset_map_global[k2] = offset0;
        }
    }

    //define new labels
    qDebug().noquote() << "Defining new labels";
    bigint L = 0;
    for (int ii = 0; ii < segments.count(); ii++) {
        L += segments[ii].times.count();
    }
    qDebug().noquote() << QString("Total number of events: %1").arg(L);
    Mda new_times(L, 1);
    Mda new_labels(L, 1);
    Mda segment_numbers(L, 1);
    Mda event_indices(L, 1);
    bigint jjjj = 0;
    for (int ii = 0; ii < segments.count(); ii++) {
        SegmentInfo* S = &segments[ii];
        for (bigint jj = 0; jj < S->labels.count(); jj++) {
            int k0 = S->labels[jj];
            int k1 = S->label_map_global.value(k0, 0);
            double offset0 = S->time_offset_map_global.value(k0, 0);
            new_times.set(S->times[jj] + offset0, jjjj);
            new_labels.set(k1, jjjj);
            segment_numbers.set(ii, jjjj);
            event_indices.set(jj, jjjj);
            jjjj++;
        }
    }

    //DiskReadMda32 X(timeseries_path);
    //Mda32 new_templates = compute_templates_in_parallel(X, new_times, new_labels, opts.clip_size);

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
    Mda FF(R, L);
    for (bigint jj = 0; jj < L; jj++) {
        int ss = segment_numbers.get(jj);
        bigint aa = event_indices.get(jj);
        for (int r = 0; r < R; r++) {
            FF.setValue(all_firings[ss].value(r, aa), r, jj);
        }
        FF.setValue(new_times.get(jj), 1, jj);
        FF.setValue(new_labels.get(jj), 2, jj);
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
    info.K = MLCompute::max(info.labels);
    info.templates = compute_templates_0(X, info.times, info.labels, opts.clip_size);
    info.firings = firings_path;
    //info.tmin = MLCompute::min(info.times);
    //info.tmax = MLCompute::max(info.times);
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

void compute_sliding_correlation_between_templates(double& best_corr_out, int& best_offset2_out, const Mda32& template1, const Mda32& template2, int offset_search_radius)
{
    int M = template1.N1();
    int T = template1.N2();
    best_corr_out = 0;
    best_offset2_out = 0;
    for (int dt = -offset_search_radius; dt <= offset_search_radius; dt++) {
        Mda32 template2b = time_shift_template(template2, dt);
        double val = compute_noncentered_correlation(M * T, template1.constDataPtr(), template2b.constDataPtr());
        if (val > best_corr_out) {
            best_corr_out = val;
            best_offset2_out = dt;
        }
    }
}

Mda32 extract_clips(const DiskReadMda32& X, const QVector<double>& times, int clip_size)
{
    bigint M = X.N1();
    bigint N = X.N2();
    bigint T = clip_size;
    bigint L = times.count();
    bigint Tmid = (bigint)((T + 1) / 2) - 1;
    Mda32 clips(M, T, L);
    for (bigint i = 0; i < L; i++) {
        bigint t1 = (bigint)times[i] - Tmid;
        bigint t2 = t1 + T - 1;
        if ((t1 >= 0) && (t2 < N)) {
            Mda32 tmp;
            X.readChunk(tmp, 0, t1, M, T);
            for (bigint t = 0; t < T; t++) {
                for (bigint m = 0; m < M; m++) {
                    clips.set(tmp.get(m, t), m, t, i);
                }
            }
        }
    }
    return clips;
}

/*
double compute_match_score(const DiskReadMda32 &X,const QVector<double> &times1_in,const QVector<double> &times2_in, P_combine_firing_segments_opts opts) {
    Mda32 clips1=extract_clips(X,times1,opts.clip_size);
    Mda32 clips2=extract_clips(X,times2,opts.clip_size);
    int M=clips1.N1();
    int T=clips1.N2();
    Mda32 clips(M,T,clips1.N3()+clips2.N3());
    for (bigint j=0; j<clips1.N3(); j++) {
        for (int t=0; t<T; t++) {
            for (int m=0; m<M; m++) {
                clips.set(clips1.get(m,t,j),m,t,j);
            }
        }
    }
    for (bigint j=0; j<clips2.N3(); j++) {
        for (int t=0; t<T; t++) {
            for (int m=0; m<M; m++) {
                clips.set(clips2.get(m,t,j),m,t,clips1.N3()+j);
            }
        }
    }

    int num_features=10;
    bigint L0=clips.N3();
    Mda32 FF;
    {
        // do this inside a code block so memory gets released
        Mda32 clips_reshaped(M * T, L0);
        bigint iii = 0;
        for (bigint ii = 0; ii < L0; ii++) {
            for (bigint t = 0; t < T; t++) {
                for (bigint m = 0; m < M; m++) {
                    clips_reshaped.set(clips.value(m, t, ii), iii);
                    iii++;
                }
            }
        }

        Mda32 CC, sigma;
        bigint max_samples=10000;
        pca_subsampled(CC, FF, sigma, clips_reshaped, num_features, false, max_samples); //should we subtract the mean?
    }

    isosplit5_opts i5_opts;
    i5_opts.isocut_threshold = 1;
    i5_opts.K_init = 200;

    QVector<int> labels0(L0);
    {
        QTime timer;
        timer.start();
        if (!isosplit5(labels0.data(), num_features, L0, FF.dataPtr(), i5_opts)) {
            qWarning() << "Isosplit5 returned with an error. Aborting";
            //clips.write32("/home/magland/tmp/debug_clips.mda");
            //FF.write32("/home/magland/tmp/debug_FF.mda");
            abort();
        }
        qDebug().noquote() << QString("Time elapsed for isosplit (%1x%2) - K=%3: %4 sec").arg(FF.N1()).arg(FF.N2()).arg(MLCompute::max(labels0)).arg(timer.elapsed() * 1.0 / 1000);
    }
    int KK=MLCompute::max(labels0);
    QVector<bigint> counts1(KK),counts2(KK);
    counts1.fill(0);
    counts2.fill(0);
    for (bigint j=0; j<clips1.N3(); j++) {
        int k=labels0[j];
        if ((1<=k)&&(k<=KK)) {
            counts1[k-1]++;
        }
    }
    for (bigint j=0; j<clips2.N3(); j++) {
        int k=labels0[clips1.N3()+j];
        if ((1<=k)&&(k<=KK)) {
            counts2[k-1]++;
        }
    }
    bigint numer=0;
    for (bigint k=1; k<=KK; k++) {
        int count1=counts1[k-1];
        int count2=counts2[k-1];
        if (count1>count2) {
            numer+=count1;
        }
        else {
            numer+=count2;
        }
    }
    bigint denom=L0;
    double frac=0;
    if (denom) frac=numer*1.0/denom;
    // if perfect split, frac=1
    // if perfect join, frac=0.5
    return 1 - (frac-0.5)/0.5;
}
*/

/*
int find_matching_cluster(QString timeseries_path,SegmentInfo* Sprev, SegmentInfo* S, int k2, P_combine_firing_segments_opts opts)
{
    DiskReadMda32 X(timeseries_path);

    QVector<double> times2;
    QVector<bigint> inds2;
    for (bigint i = 0; i < S->times.count(); i++) {
        //if (S->times[i] <= Sprev->tmax) {
            if (S->labels[i] == k2) {
                times2 << S->times[i];
                inds2 << i;
            }
        //}
    }
    if (times2.count() == 0)
        return 0;
    int M = S->templates.N1();
    int T = S->templates.N2();
    Mda32 template2;
    S->templates.getChunk(template2, 0, 0, k2 - 1, M, T, 1);
    double best_match_score = 0;
    int best_match_k1 = 0;
    for (int k1 = 1; k1 <= Sprev->K; k1++) {
        Mda32 template1;
        Sprev->templates.getChunk(template1, 0, 0, k1 - 1, M, T, 1);
        double corr = 0;
        int best_offset2 = 0;
        compute_sliding_correlation_between_templates(corr, best_offset2, template1, template2, opts);
        if (corr > opts.template_correlation_threshold) {
            QVector<double> times1;
            for (bigint i = 0; i < Sprev->times.count(); i++) {
                //if (Sprev->times[i] >= S->tmin) {
                    if (Sprev->labels[i] == k1)
                        times1 << Sprev->times[i];
                //}
            }
            QVector<double> times2_adjusted(times2.count());
            for (bigint j=0; j<times2.count(); j++) {
                times2_adjusted[j]=times2[j]+best_offset2;
            }
            double match_score=compute_match_score(X,times1,times2_adjusted,opts);
            if (match_score >= opts.match_score_threshold) {
                if (match_score > best_match_score) {
                    best_match_score = match_score;
                    best_match_k1 = k1;
                }
            }
        }
    }
    if (best_match_score >= opts.match_score_threshold) {
        return best_match_k1;
    }
    else
        return 0;
}
*/

Mda32 compute_mean_clip(const Mda32& clips)
{
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    Mda32 ret;
    ret.allocate(M, T);
    int aaa = 0;
    for (int i = 0; i < L; i++) {
        int bbb = 0;
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                ret.set(ret.get(bbb) + clips.get(aaa), bbb);
                aaa++;
                bbb++;
            }
        }
    }
    if (L) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                ret.set(ret.get(m, t) / L, m, t);
            }
        }
    }
    return ret;
}

double compute_template_similarity_score(const DiskReadMda32& X, SegmentInfo* Sprev, SegmentInfo* S, int k1, int k2, P_combine_firing_segments_opts opts)
{
    QVector<double> times1, times2;
    for (bigint i = 0; i < Sprev->times.count(); i++) {
        if (Sprev->labels[i] == k1)
            times1 << Sprev->times[i];
    }
    for (bigint i = 0; i < S->times.count(); i++) {
        if (S->labels[i] == k2)
            times2 << S->times[i];
    }
    qSort(times1);
    qSort(times2);
    bigint num_events = opts.num_comparison_events;
    if (times1.count() < num_events)
        num_events = times1.count();
    if (times2.count() < num_events)
        num_events = times2.count();
    times1 = times1.mid(times1.count() - num_events);
    times2 = times2.mid(0, num_events);
    Mda32 clips1 = extract_clips(X, times1, opts.clip_size);
    Mda32 clips2 = extract_clips(X, times2, opts.clip_size);
    Mda32 template1 = compute_mean_clip(clips1);
    Mda32 template2 = compute_mean_clip(clips2);
    double corr;
    int offset2;
    compute_sliding_correlation_between_templates(corr, offset2, template1, template2, opts.offset_search_radius);
    return corr;
}

void get_ending_times_labels(QVector<double>& times_out, QVector<int>& labels_out, const QVector<double>& times, const QVector<int>& labels, P_combine_firing_segments_opts opts, bool use_beginning)
{
    QVector<double> times_ret;
    QVector<int> labels_ret;
    int K = MLCompute::max(labels);
    for (int k = 1; k <= K; k++) {
        QVector<double> times_k;
        for (bigint i = 0; i < times.count(); i++) {
            if (labels[i] == k) {
                times_k << times[i];
            }
        }
        qSort(times_k);
        if (times_k.count() > opts.num_comparison_events) {
            if (use_beginning)
                times_k = times_k.mid(0, opts.num_comparison_events);
            else
                times_k = times_k.mid(times_k.count() - opts.num_comparison_events);
        }
        for (bigint aa = 0; aa < times_k.count(); aa++) {
            times_ret << times_k[aa];
            labels_ret << k;
        }
    }
    QList<bigint> sort_inds = get_sort_indices_bigint(times_ret);
    times_out.clear();
    labels_out.clear();
    for (bigint aa = 0; aa < sort_inds.count(); aa++) {
        times_out << times_ret[sort_inds[aa]];
        labels_out << labels_ret[sort_inds[aa]];
    }
}
void get_beginning_times_labels(QVector<double>& times_out, QVector<int>& labels_out, const QVector<double>& times, const QVector<int>& labels, P_combine_firing_segments_opts opts)
{
    get_ending_times_labels(times_out, labels_out, times, labels, opts, true);
}

Mda32 compute_clips_features(const Mda32& clips, int num_features)
{
    int M = clips.N1();
    int T = clips.N2();
    bigint L0 = clips.N3();
    Mda32 FF;
    {
        // do this inside a code block so memory gets released
        Mda32 clips_reshaped(M * T, L0);
        bigint iii = 0;
        for (bigint ii = 0; ii < L0; ii++) {
            for (bigint t = 0; t < T; t++) {
                for (bigint m = 0; m < M; m++) {
                    clips_reshaped.set(clips.value(m, t, ii), iii);
                    iii++;
                }
            }
        }

        Mda32 CC, sigma;
        bigint max_samples = 10000;
        pca_subsampled(CC, FF, sigma, clips_reshaped, num_features, false, max_samples); //should we subtract the mean?
    }
    return FF;
}

QVector<bigint> find_nearest_neighbors(Mda32& FF1_in, Mda32& FF2_in, int num_features)
{
    Mda32 FF1, FF2;
    if (num_features) {
        Mda32 FFF(FF1_in.N1(), FF1_in.N2() + FF2_in.N2());
        FFF.setChunk(FF1_in, 0, 0);
        FFF.setChunk(FF2_in, 0, FF1_in.N2());
        Mda32 CC, sigma, features;
        bigint max_samples = 10000;
        pca_subsampled(CC, features, sigma, FFF, num_features, false, max_samples);
        features.getChunk(FF1, 0, 0, FF1_in.N1(), FF1_in.N2());
        features.getChunk(FF2, 0, FF1_in.N2(), FF1_in.N1(), FF2_in.N2());
    }
    else {
        FF1 = FF1_in;
        FF2 = FF2_in;
    }
    KdTree tree;
    tree.create(FF1);
    QVector<bigint> ret;
    for (bigint i = 0; i < FF2.N2(); i++) {
        QVector<float> V;
        for (int a = 0; a < FF1.N1(); a++)
            V << FF2.get(a, i);
        QList<int> inds = tree.findApproxKNearestNeighbors(FF1, V, 1, 100);
        ret << inds[0];
    }
    return ret;
}

void write_mda(const QVector<int>& X, QString fname)
{
    Mda A(X.count(), 1);
    for (bigint i = 0; i < X.count(); i++) {
        A.set(X[i], i);
    }
    A.write64(fname);
}
void write_mda(const QVector<bigint>& X, QString fname)
{
    Mda A(X.count(), 1);
    for (bigint i = 0; i < X.count(); i++) {
        A.set(X[i], i);
    }
    A.write64(fname);
}

Mda32 compute_templates_from_clips_and_labels(const Mda32& clips, const QVector<int>& labels, int K)
{
    int M = clips.N1();
    int T = clips.N2();
    int L = labels.count();

    Mda32 templates(M, T, K);
    QVector<int> counts(K);
    counts.fill(0);
    for (int i = 0; i < L; i++) {
        int k = labels[i];
        if ((k >= 1) && (k <= K)) {
            Mda32 X0;
            clips.getChunk(X0, 0, 0, i, M, T, 1);
            dtype32* Xptr = X0.dataPtr();
            dtype32* Tptr = templates.dataPtr(0, 0, k - 1);
            for (int i = 0; i < M * T; i++) {
                Tptr[i] += Xptr[i];
            }
            counts[k - 1]++;
        }
    }
    for (int k = 1; k <= K; k++) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                if (counts[k - 1]) {
                    templates.set(templates.get(m, t, k - 1) / counts[k - 1], m, t, k - 1);
                }
            }
        }
    }
    return templates;
}

QVector<int> find_optimal_alignment_shifts(const Mda32& templates, const Mda32& template_ref, int offset_search_radius)
{
    int M = templates.N1();
    int T = templates.N2();
    QVector<int> ret(templates.N3());
    ret.fill(0);
    for (int i = 0; i < templates.N3(); i++) {
        Mda32 template1;
        templates.getChunk(template1, 0, 0, i, M, T, 1);
        double corr;
        int offset2;
        compute_sliding_correlation_between_templates(corr, offset2, template1, template_ref, offset_search_radius);
        ret[i] = -offset2;
    }
    return ret;
}

Mda32 shift_clip(const Mda32& clip, int shift)
{
    if (!clip.N2())
        return clip;
    Mda32 ret(clip.N1(), clip.N2());
    for (int t = 0; t < clip.N2(); t++) {
        int t2 = t - shift;
        while (t2 < 0)
            t2 += clip.N2();
        while (t2 >= clip.N2())
            t2 -= clip.N2();
        for (int m = 0; m < clip.N1(); m++) {
            ret.set(clip.get(m, t2), m, t);
        }
    }
    return ret;
}

Mda32 align_clips(const Mda32& clips1, const QVector<int>& labels1, const Mda32& templates1, const Mda32& template_k2, int offset_search_radius)
{
    QVector<int> offsets1 = find_optimal_alignment_shifts(templates1, template_k2, offset_search_radius);
    int M = clips1.N1();
    int T = clips1.N2();
    Mda32 clips1_aligned(M, T, clips1.N3());
    for (bigint i = 0; i < clips1.N3(); i++) {
        Mda32 clip;
        clips1.getChunk(clip, 0, 0, i, M, T, 1);
        int k = labels1[i];
        if ((1 <= k) && (k <= offsets1.count())) {
            int offset = offsets1[k - 1];
            clip = shift_clip(clip, offset);
        }
        clips1_aligned.setChunk(clip, 0, 0, i);
    }
    return clips1_aligned;
}
}
