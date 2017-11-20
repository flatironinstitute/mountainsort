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
#include "fit_stage.h"

typedef QList<bigint> BigintList;

QList<bigint> get_time_channel_mask(const Mda32& template0, double thresh);
bool is_dirty(float* dirty_ptr, const QList<bigint>& tchmask);
double compute_score(bigint M, bigint T, float* X_ptr, float* template0, const QList<bigint>& tchmask);
QVector<bigint> find_events_to_use(const QVector<double>& times, const QVector<double>& scores, const Fit_stage_opts& opts, int clip_size);
void subtract_scaled_template(bigint N, double* X, double* template0, double scale_min, double scale_max);
void subtract_scaled_template(bigint M, bigint T, float* X_ptr, float* dirty_ptr, float* template0, const QList<bigint>& tchmask, double scale_min, double scale_max);

QVector<bigint> fit_stage(Mda32& X, const QVector<double>& times, const QVector<int>& labels, Mda32& templates, Fit_stage_opts opts)
{
    int M = X.N1();
    //bigint N=X.N2();
    int T = templates.N2();
    int K = templates.N3();
    bigint L = times.count();

    QVector<bool> ret(L, false);

    QList<BigintList> time_channel_mask;
    QList<int> mask_sizes_for_display;
    for (bigint i = 0; i < K; i++) {
        Mda32 template0, template0_stdev;
        templates.getChunk(template0, 0, 0, i, M, T, 1);
        time_channel_mask << get_time_channel_mask(template0, opts.time_channel_mask_thresh); //use only the channels with highest maxval
        mask_sizes_for_display << time_channel_mask[i].count();
    }
    //qDebug().noquote() << QString("Mask sizes:") << mask_sizes_for_display;

    bigint Tmid = (bigint)((T + 1) / 2) - 1; //the center timepoint in a clip (zero-indexed)

    //compute the L2-norms of the templates ahead of time
    QVector<double> template_norms;
    template_norms << 0;
    for (bigint k = 1; k <= K; k++) {
        template_norms << MLCompute::norm(M * T, templates.dataPtr(0, 0, k - 1));
    }

    QVector<double> scores(L);
    Mda32 dirty(X.N1(), X.N2()); //timepoints/channels for which we need to recompute the score if an event overlaps
    for (bigint ii = 0; ii < dirty.totalSize(); ii++) {
        dirty.setValue(1, ii);
    }

    QVector<bigint> event_inds_to_consider;
    for (bigint i = 0; i < L; i++)
        event_inds_to_consider << i;
    QVector<bigint> event_inds_to_use;

    bigint num_score_computes = 0;
    bigint num_to_use = 0;
    //keep passing through the data until nothing changes anymore
    bool something_changed = true;
    bigint num_passes = 0;
    //while ((something_changed)&&(num_passes<2)) {
    while (something_changed) {
        num_passes++;

        QList<bigint> inds_to_try; //indices of the events to try on this pass
        QVector<double> scores_to_try;
        QVector<double> times_to_try;
        QVector<bigint> labels_to_try;

        //QVector<double> template_norms_to_try;
        QVector<bigint> new_event_inds_to_consider;
        for (bigint kk = 0; kk < event_inds_to_consider.count(); kk++) { //loop through the events
            bigint i = event_inds_to_consider[kk];

            double t0 = times[i];
            bigint k0 = labels[i];
            if (k0 > 0) { //make sure we have a positive label (don't know why we wouldn't)
                bigint tt = (bigint)(t0 - Tmid + 0.5); //start time of clip
                double score0 = 0;
                if ((tt >= 0) && (tt + T <= X.N2())) { //make sure we are in range
                    BigintList tchmask = time_channel_mask[k0 - 1];
                    if (!is_dirty(dirty.dataPtr(0, tt), tchmask)) {
                        // we don't need to recompute the score
                        score0 = scores[i];
                    }
                    else {
                        //we do need to recompute it.

                        //The score will be how much something like the L2-norm is decreased
                        score0 = compute_score(M, T, X.dataPtr(0, tt), templates.dataPtr(0, 0, k0 - 1), tchmask);
                        num_score_computes++;
                        /*
                        if (score0 < template_norms[k0] * template_norms[k0] * 0.1)
                            score0 = 0; //the norm of the improvement needs to be at least 0.5 times the norm of the template
                            */
                        scores[i] = score0;
                    }
                }
                //the score needs to be at least as large as neglogprior in order to accept the spike
                //double neglogprior = 30;
                double neglogprior = 0;
                if (score0 > neglogprior) {
                    //we are not committing to using this event yet... see below for next step
                    scores_to_try << score0;
                    times_to_try << t0;
                    labels_to_try << k0;
                    inds_to_try << i;
                }
                else {
                    //not gonna use it
                }
            }
        }
        //Look at those events to try and see if we should use them
        QVector<bigint> to_use = find_events_to_use(times_to_try, scores_to_try, opts, templates.N2());

        //at this point, nothing is dirty
        for (bigint i = 0; i < dirty.totalSize(); i++) {
            dirty.set(0, i);
        }

        //amplitude scaling for template subtraction
        double scale_min = 1;
        double scale_max = 1;

        //for all those we are going to "use", we want to subtract out the corresponding templates from the timeseries data
        something_changed = false;
        bigint num_added = 0;
        for (bigint aa = 0; aa < to_use.count(); aa++) {
            if (to_use[aa] == 1) {
                BigintList tchmask = time_channel_mask[labels_to_try[aa] - 1];
                something_changed = true;
                num_added++;
                bigint tt = (bigint)(times_to_try[aa] - Tmid + 0.5);
                subtract_scaled_template(M, T, X.dataPtr(0, tt), dirty.dataPtr(0, tt), templates.dataPtr(0, 0, labels_to_try[aa] - 1), tchmask, scale_min, scale_max);
                event_inds_to_use << inds_to_try[aa];
                num_to_use++;
            }
            else {
                //consider it for next time
                new_event_inds_to_consider << inds_to_try[aa];
            }
        }
        event_inds_to_consider = new_event_inds_to_consider;
    }

    //qDebug().noquote() << QString("Processed %1 events in time chunk (using %2%), %3 score computes per event, %4 passes").arg(L).arg(num_to_use * 100.0 / L).arg(num_score_computes * 1.0 / L).arg(num_passes);

    return event_inds_to_use;
}

QList<bigint> get_time_channel_mask(const Mda32& template0, double thresh)
{
    bigint M = template0.N1();
    bigint T = template0.N2();
    QList<bigint> ret;
    int ii = 0;
    for (bigint t = 0; t < T; t++) {
        for (bigint m = 0; m < M; m++) {
            double val = qAbs(template0.value(m, t));
            if (val > thresh)
                ret << ii;
            ii++;
        }
    }
    return ret;
}

bool is_dirty(float* dirty_ptr, const QList<bigint>& tchmask)
{
    for (bigint i = 0; i < tchmask.count(); i++) {
        if (dirty_ptr[i])
            return true;
    }
    return false;
}

double compute_score(bigint M, bigint T, float* X_ptr, float* template0, const QList<bigint>& tchmask)
{
    (void)M;
    (void)T;
    double before_sumsqr = 0;
    double after_sumsqr = 0;
    for (bigint j = 0; j < tchmask.count(); j++) {
        bigint i = tchmask[j];
        double val = X_ptr[i];
        before_sumsqr += val * val;
        val -= template0[i];
        after_sumsqr += val * val;
    }
    return before_sumsqr - after_sumsqr;
}

QVector<bigint> find_events_to_use(const QVector<double>& times, const QVector<double>& scores, const Fit_stage_opts& opts, int clip_size)
{
    (void)opts;
    QVector<bigint> to_use;
    bigint L = times.count();
    for (bigint i = 0; i < L; i++)
        to_use << 0; //start out not using any
    for (bigint i = 0; i < L; i++) {
        if (scores[i] > 0) { //score has to at least be positive
            to_use[i] = 1; //for now we say we are using it
            {
                // but let's check nearby events that may have a larger score

                bigint j = i;
                while ((j >= 0) && (times[j] >= times[i] - clip_size)) {
                    if ((i != j) && (scores[j] >= scores[i])) {
                        to_use[i] = 0; //actually not using it because there is something bigger to the left
                    }
                    j--;
                }
            }
            {
                bigint j = i;
                while ((j < times.count()) && (times[j] <= times[i] + clip_size)) {
                    if (scores[j] > scores[i]) {
                        to_use[i] = 0; //actually not using it because there is something bigger to the right
                    }
                    j++;
                }
            }
        }
    }
    return to_use;
}

void subtract_scaled_template(bigint N, double* X, double* template0, double scale_min, double scale_max)
{
    double S12 = 0, S22 = 0;
    for (bigint i = 0; i < N; i++) {
        S22 += template0[i] * template0[i];
        S12 += X[i] * template0[i];
    }
    double alpha = 1;
    if (S22)
        alpha = S12 / S22;
    alpha = qMin(scale_max, qMax(scale_min, alpha));
    for (bigint i = 0; i < N; i++) {
        X[i] -= alpha * template0[i];
    }
}

void subtract_scaled_template(bigint M, bigint T, float* X_ptr, float* dirty_ptr, float* template0, const QList<bigint>& tchmask, double scale_min, double scale_max)
{
    (void)M;
    (void)T;
    double S12 = 0, S22 = 0;
    for (bigint j = 0; j < tchmask.count(); j++) {
        bigint i = tchmask[j];
        S22 += template0[i] * template0[i];
        S12 += X_ptr[i] * template0[i];
    }
    double alpha = 1;
    if (S22)
        alpha = S12 / S22;
    alpha = qMin(scale_max, qMax(scale_min, alpha));
    for (bigint j = 0; j < tchmask.count(); j++) {
        bigint i = tchmask[j];
        X_ptr[i] -= alpha * template0[i];
        dirty_ptr[i] = 1;
    }
}
