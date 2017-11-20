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
#include "p_compute_templates.h"

#include <QTime>
#include <diskreadmda.h>
#include <diskreadmda32.h>
#include "mlutil.h"

bool p_compute_templates(QStringList timeseries_list, QString firings_path, QString templates_out, int clip_size, const QList<int>& clusters_in)
{
    QList<int> clusters = clusters_in;
    DiskReadMda32 X(2, timeseries_list);
    DiskReadMda firings(firings_path);

    bigint M = X.N1();
    //bigint N = X.N2();
    bigint T = clip_size;
    QVector<double> times;
    QVector<bigint> labels;

    if (!T) {
        qWarning() << "Unexpected: Clip size is zero.";
        return false;
    }

    if (clusters.isEmpty()) {
        int Kmax = 0;
        for (bigint i = 0; i < firings.N2(); i++) {
            int label0 = firings.value(2, i);
            if (label0 > Kmax)
                Kmax = label0;
        }
        for (int kk = 1; kk <= Kmax; kk++)
            clusters << kk;
    }
    if (clusters.isEmpty()) {
        //there must really be no clusters. so we should return no templates
        Mda32 templates00;
        return templates00.write32(templates_out);
    }

    bigint Kmax = MLCompute::max(clusters.toVector());
    QVector<int> label_map(Kmax + 1);
    label_map.fill(-1);
    for (int c = 0; c < clusters.count(); c++) {
        label_map[clusters[c]] = c;
    }

    QSet<int> clusters_set = clusters.toSet();

    for (bigint i = 0; i < firings.N2(); i++) {
        int label0 = firings.value(2, i);
        if (clusters_set.contains(label0)) {
            times << firings.value(1, i);
            labels << label0;
        }
    }

    bigint K0 = clusters.count();

    Mda templates(M, T, K0);
    QVector<double> counts(K0);
    counts.fill(0);

    double* templates_ptr = templates.dataPtr();

    bigint Tmid = (bigint)((T + 1) / 2) - 1;
    printf("computing templates (M=%ld,T=%ld,K=%ld,L=%d)...\n", M, T, K0, times.count());
    QTime timer;
    timer.start();
    for (bigint i = 0; i < times.count(); i++) {
        bigint t1 = times[i] - Tmid;
        //bigint t2 = t1 + T - 1;
        int k = label_map[labels[i]];
        if (timer.elapsed() > 3000) {
            qDebug().noquote() << QString("Compute templates: processing event %1 of %2").arg(i).arg(times.count());
            timer.restart();
        }
        {
            Mda32 tmp;
            tmp.allocate(M, T);
            if (!X.readChunk(tmp, 0, t1, M, T)) {
                qWarning() << "Problem reading chunk of timeseries list:" << t1;
                return false;
            }
            float* tmp_ptr = tmp.dataPtr();
            bigint offset = M * T * k;
            for (bigint aa = 0; aa < M * T; aa++) {
                templates_ptr[offset + aa] += tmp_ptr[aa];
            }
            counts[k]++;
        }
    }
    bigint bb = 0;
    for (bigint k = 0; k < K0; k++) {
        for (bigint aa = 0; aa < M * T; aa++) {
            if (counts[k])
                templates_ptr[bb] /= counts[k];
            bb++;
        }
    }

    return templates.write32(templates_out);
}
