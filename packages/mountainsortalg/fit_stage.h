#ifndef FIT_STAGE_H
#define FIT_STAGE_H

#include <QVector>
#include "mda32.h"

struct Fit_stage_opts {
    double time_channel_mask_thresh = 0.1;
};

QVector<bigint> fit_stage(Mda32& X, const QVector<double>& times, const QVector<int>& labels, Mda32& templates, Fit_stage_opts opts);

#endif // FIT_STAGE_H
