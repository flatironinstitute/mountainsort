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

#include "omp.h"
#include "fit_stage.h"

class FitStageComputerPrivate {
public:
    FitStageComputer* q;

    int m_num_threads = 1;
    QVector<double> m_times;
    QVector<int> m_labels;
    Mda32 m_templates;

    QVector<bigint> m_event_inds_to_use;
};

FitStageComputer::FitStageComputer()
{
    d = new FitStageComputerPrivate;
    d->q = this;
}

FitStageComputer::~FitStageComputer()
{
    delete d;
}

void FitStageComputer::setTemplates(const Mda32& templates)
{
    d->m_templates = templates;
}

void FitStageComputer::setTimesLabels(const QVector<double>& times, const QVector<int>& labels)
{
    d->m_times = times;
    d->m_labels = labels;
}

void FitStageComputer::processTimeChunk(bigint t, const Mda32& X, bigint padding_left, bigint padding_right)
{
    //int M = X.N1();
    //int T = d->m_templates.N2();
    //int Tmid = (int)((T + 1) / 2) - 1;

    {
        QVector<bigint> local_event_inds;
        QVector<double> local_times;
        QVector<int> local_labels;
        Mda32 local_templates;
        Mda32 local_X;
#pragma omp critical(fit_stage_computer)
        {
            for (bigint a = 0; a < d->m_times.count(); a++) {
                double t0 = d->m_times[a];
                double t1 = t0 - t + padding_left;
                if ((0 <= t1) && (t1 < X.N2())) {
                    local_event_inds << a;
                    local_times << t1;
                    local_labels << d->m_labels[a];
                }
            }
            local_templates = d->m_templates;
            local_X = X;
        }
        Fit_stage_opts oo;
        QVector<bigint> inds_to_use_00 = fit_stage(local_X, local_times, local_labels, local_templates, oo);
//debug
#pragma omp critical(fit_stage_computer)
        {
            for (bigint j = 0; j < inds_to_use_00.count(); j++) {
                double t1 = local_times[inds_to_use_00[j]];
                if ((padding_left <= t1) && (t1 < X.N2() - padding_left - padding_right)) {
                    d->m_event_inds_to_use << local_event_inds[inds_to_use_00[j]];
                }
            }
        }
    }
}

void FitStageComputer::finalize()
{
    qSort(d->m_event_inds_to_use);
}

QVector<bigint> FitStageComputer::eventIndicesToUse() const
{
    return d->m_event_inds_to_use;
}
