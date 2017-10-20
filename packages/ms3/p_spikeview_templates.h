#ifndef P_SPIKEVIEW_TEMPLATES_H
#define P_SPIKEVIEW_TEMPLATES_H

#include <QString>
#include "mdaio.h"

struct P_spikeview_templates_opts {
    int clip_size = 100;
    bigint max_events_per_template = 1000;
    int filt_padding = 0;
    double freq_min = 0;
    double freq_max = 0;
    double freq_wid = 1000;
    double samplerate = 0;
    bool subtract_temporal_mean = false;
};

bool p_spikeview_templates(QString timeseries, QString firings, QString templates_out, P_spikeview_templates_opts opts);

#endif // P_SPIKEVIEW_TEMPLATES_H
