#include "p_mv_compute_templates.h"

#include <diskreadmda.h>
#include "compute_templates_0.h"

/// TODO 0.9.1 #define USE_TASK_PROGRESS, and don't use this outside of the mountainview gui -- unnecessary dependence AND risk

bool mv_compute_templates(const QString& timeseries_path, const QString& firings_path, const QString& templates_out_path, const QString& stdevs_out_path, int clip_size)
{
    DiskReadMda X(timeseries_path);
    if (X.N2() <= 1) {
        return false;
    }
    DiskReadMda firings(firings_path);
    QVector<double> times;
    QVector<int> labels;
    for (bigint i = 0; i < firings.N2(); i++) {
        times << firings.value(1, i);
        labels << (int)firings.value(2, i);
    }
    Mda templates, stdevs;
    compute_templates_stdevs(templates, stdevs, X, times, labels, clip_size);
    templates.write32(templates_out_path);
    stdevs.write32(stdevs_out_path);
    return true;
}
