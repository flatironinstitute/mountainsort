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
#include "p_mv_discrimhist.h"

#include "diskreadmda.h"
#include "extract_clips.h"
#include "mlutil.h"

struct discrimhist_data {
    int k1, k2;
    QVector<double> data1;
    QVector<double> data2;
};

/// TODO parallelize mv_distrimhist
bool mv_discrimhist(QString timeseries_path, QString firings_path, QString output_path, mv_discrimhist_opts opts)
{
    DiskReadMda32 timeseries(timeseries_path);
    DiskReadMda firings(firings_path);

    QList<discrimhist_data> datas;
    for (int i1 = 0; i1 < opts.clusters.count(); i1++) {
        for (int i2 = i1 + 1; i2 < opts.clusters.count(); i2++) {
            int k1 = opts.clusters[i1];
            int k2 = opts.clusters[i2];
            discrimhist_data DD;
            DD.k1 = k1;
            DD.k2 = k2;
            if (!get_discrimhist_data(DD.data1, DD.data2, timeseries, firings, k1, k2, opts.clip_size, opts.method)) {
                return false;
            }
            datas << DD;
        }
    }

    int total_count = 0;
    for (int i = 0; i < datas.count(); i++) {
        total_count += datas[i].data1.count();
        total_count += datas[i].data2.count();
    }

    Mda output(4, total_count);
    //first two rows are k1/k2, third row is k1 or k2, fourth row is the value
    int jj = 0;
    for (int i = 0; i < datas.count(); i++) {
        int k1 = datas[i].k1;
        int k2 = datas[i].k2;
        for (int k = 0; k < datas[i].data1.count(); k++) {
            output.setValue(k1, 0, jj);
            output.setValue(k2, 1, jj);
            output.setValue(k1, 2, jj);
            output.setValue(datas[i].data1[k], 3, jj);
            jj++;
        }
        for (int k = 0; k < datas[i].data2.count(); k++) {
            output.setValue(k1, 0, jj);
            output.setValue(k2, 1, jj);
            output.setValue(k2, 2, jj);
            output.setValue(datas[i].data2[k], 3, jj);
            jj++;
        }
    }

    output.write32(output_path);

    return true;
}

double compute_dot_product(int N, float* v1, float* v2)
{
    double ip = 0;
    for (int m = 0; m < N; m++)
        ip += v1[m] * v2[m];
    return ip;
}

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


bool get_discrimhist_data(QVector<double>& ret1, QVector<double>& ret2, const DiskReadMda32& timeseries, const DiskReadMda& firings, int k1, int k2, int clip_size, QString method)
{
    //Assemble the times where neuron k1 fires (times1) and where neuron k2 fires (times2)
    QVector<double> times1, times2;
    for (int i = 0; i < firings.N2(); i++) {
        int label = (int)firings.value(2, i);
        if (label == k1) {
            times1 << firings.value(1, i);
        }
        if (label == k2) {
            times2 << firings.value(1, i);
        }
    }

    //extract the clips
    Mda32 clips1 = extract_clips(timeseries, times1, clip_size);
    Mda32 clips2 = extract_clips(timeseries, times2, clip_size);
    float* ptr_clips1 = clips1.dataPtr();
    float* ptr_clips2 = clips2.dataPtr();

    //the direction to project onto and the cutoff to separate the two clusters
    Mda32 discrim_direction;
    double cutoff = 0;

    //compute the centroids
    Mda32 centroid1 = compute_mean_clip(clips1);
    Mda32 centroid2 = compute_mean_clip(clips2);
    int M = centroid1.N1();
    int T = centroid1.N2();
    int MT = M * T;
    //int L1 = clips1.N3();
    //int L2 = clips2.N3();

    //allocate the discrim direction
    discrim_direction.allocate(centroid1.N1(), centroid1.N2());
    float* ptr_dd = discrim_direction.dataPtr();

    if (method == "centroid") {
        //the direction will be the vector connecting the two centroids
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                discrim_direction.setValue(centroid2.value(m, t) - centroid1.value(m, t), m, t);
            }
        }
    }
    /*
    else if (method == "svm") {
        //the direction and cutoff are determined using support vector machine
        QVector<double> direction0;
        Mda32 X0(MT, L1 + L2); //concatenation (and reshaped) of clips1 and clips2
        QVector<int> labels0;
        for (int j = 0; j < L1; j++) {
            Mda32 tmp;
            clips1.getChunk(tmp, 0, 0, j, M, T, 1);
            for (int k = 0; k < MT; k++) {
                X0.setValue(tmp.value(k), k, j);
            }
            labels0 << 1;
        }
        for (int j = 0; j < L2; j++) {
            Mda32 tmp;
            clips2.getChunk(tmp, 0, 0, j, M, T, 1);
            for (int k = 0; k < MT; k++) {
                X0.setValue(tmp.value(k), k, L2 + j);
            }
            labels0 << 2;
        }
        get_svm_discrim_direction(cutoff, direction0, X0, labels0);
        for (int i = 0; i < MT; i++) {
            discrim_direction.setValue(direction0[i], i);
        }
    }
    */
    else {
        qWarning() << "Unsupported discrimhist method: " + method;
        return false;
    }

    //normalize the discrim direction
    double norm0 = MLCompute::norm(MT, ptr_dd);
    if (!norm0)
        norm0 = 1;

    ret1.clear();
    for (int i = 0; i < clips1.N3(); i++) {
        ret1 << (compute_dot_product(MT, ptr_dd, &ptr_clips1[MT * i]) - cutoff) / (norm0 * norm0);
    }
    ret2.clear();
    for (int i = 0; i < clips2.N3(); i++) {
        ret2 << (compute_dot_product(MT, ptr_dd, &ptr_clips2[MT * i]) - cutoff) / (norm0 * norm0);
    }

    return true;
}
