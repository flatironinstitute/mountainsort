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
#include "p_mv_compute_amplitudes.h"

#include "compute_templates_0.h"
#include "extract_clips.h"
#include "mlutil.h"

#include <QFile>

QVector<int> find_label_inds_00(const QVector<int>& labels, int k);
double compute_max_amp_of_template(int* max_amp_index, const Mda32& template0);

bool p_mv_compute_amplitudes(QString timeseries_path, QString firings_path, QString firings_out_path, p_mv_compute_amplitudes_opts opts)
{
    DiskReadMda32 X(timeseries_path);
    DiskReadMda firings(firings_path);
    QVector<double> times;
    QVector<int> labels;
    for (int i = 0; i < firings.N2(); i++) {
        times << firings.value(1, i);
        labels << firings.value(2, i);
    }
    Mda32 templates = compute_templates_0(X, times, labels, opts.clip_size);
    int K = MLCompute::max(labels);
    int A = qMax((long long)firings.N1(), 4LL);
    Mda firings_out(A, firings.N2());
    for (int i = 0; i < firings.N2(); i++) {
        for (int a = 0; a < A; a++) {
            firings_out.setValue(firings.value(a, i), a, i);
        }
        firings_out.setValue(0, 3, i);
    }
    for (int k = 1; k <= K; k++) {
        QVector<int> inds_k = find_label_inds_00(labels, k);
        QVector<double> times_k(inds_k.count());
        for (int i = 0; i < inds_k.count(); i++) {
            times_k[i] = times[inds_k[i]];
        }
        Mda32 template0;
        templates.getChunk(template0, 0, 0, k - 1, templates.N1(), templates.N2(), 1);
        int max_amp_index;
        double max_amp = compute_max_amp_of_template(&max_amp_index, template0);
        //double template_norm=MLCompute::norm(template0.totalSize(),template0.constDataPtr());
        if (max_amp) {
            Mda32 clips_k = extract_clips(X, times_k, opts.clip_size);
            for (int i = 0; i < inds_k.count(); i++) {
                Mda32 clip0;
                clips_k.getChunk(clip0, 0, 0, i, clips_k.N1(), clips_k.N2(), 1);
                /*
                //this was the tricky way to do it, but the simpler way is much better
                double tmp = MLCompute::dotProduct(clip0.N1() * clip0.N2(), clip0.dataPtr(), template0.dataPtr());
                double amp0=tmp/(template_norm*template_norm)*max_amp;
                */
                /// TODO: This involves a lot of extra unnecessary computations
                double amp0 = clip0.value(max_amp_index);
                firings_out.setValue(amp0, 3, inds_k[i]);
            }
        }
    }

    return firings_out.write64(firings_out_path);
}

QVector<int> find_label_inds_00(const QVector<int>& labels, int k)
{
    QVector<int> ret;
    for (int i = 0; i < labels.count(); i++) {
        if (labels[i] == k)
            ret << i;
    }
    return ret;
}

double compute_max_amp_of_template(int* max_amp_index, const Mda32& template0)
{
    int N = template0.totalSize();
    const float* ptr = template0.constDataPtr();
    double ret = 0;
    *max_amp_index = 0;
    for (int i = 0; i < N; i++) {
        if (qAbs(ptr[i]) > qAbs(ret)) {
            ret = ptr[i];
            *max_amp_index = i;
        }
    }
    return ret;
}
