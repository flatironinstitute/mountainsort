#ifndef P_SPIKEVIEW_METRICS_H
#define P_SPIKEVIEW_METRICS_H

#include <QString>

struct P_spikeview_metrics1_opts {
    double samplerate = 0;
};

bool p_spikeview_metrics1(QString firings, QString metrics_out, P_spikeview_metrics1_opts opts);

#endif // P_SPIKEVIEW_METRICS_H
