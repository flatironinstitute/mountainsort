#ifndef P_MV_COMPUTE_AMPLITUDES_H
#define P_MV_COMPUTE_AMPLITUDES_H

#include <QString>

struct p_mv_compute_amplitudes_opts {
    int clip_size = 50;
};

bool p_mv_compute_amplitudes(QString timeseries_path, QString firings_path, QString firings_out_path, p_mv_compute_amplitudes_opts opts);

#endif // P_MV_COMPUTE_AMPLITUDES_H
