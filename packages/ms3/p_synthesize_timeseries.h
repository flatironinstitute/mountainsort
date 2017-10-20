#ifndef P_SYNTHESIZE_TIMESERIES_H
#define P_SYNTHESIZE_TIMESERIES_H

#include <QString>
#include "mdaio.h"

struct P_synthesize_timeseries_opts {
    double noise_level = 0;
    bigint duration = 0;
    int waveform_upsample_factor = 13;
};

bool p_synthesize_timeseries(QString firings, QString waveforms, QString timeseries_out, P_synthesize_timeseries_opts opts);

#endif // P_SYNTHESIZE_TIMESERIES_H
