#ifndef P_COMPUTE_AMPLITUDES_H
#define P_COMPUTE_AMPLITUDES_H

#include <QString>
#include "mlutil.h"

struct P_compute_amplitudes_opts {
    int central_channel = 0;
};

bool p_compute_amplitudes(QString timeseries, QString event_times, QString amplitudes_out, P_compute_amplitudes_opts opts);

#endif // P_COMPUTE_AMPLITUDES_H
