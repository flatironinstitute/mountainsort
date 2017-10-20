#ifndef P_BANDPASS_FILTER_H
#define P_BANDPASS_FILTER_H

#include <QString>

struct Bandpass_filter_opts {
    double samplerate = 0;
    double freq_min = 0;
    double freq_max = 0;
    double freq_wid = 0;
    double quantization_unit = 0;
    int subsample_factor = 1;
};

bool p_bandpass_filter(QString timeseries, QString timeseries_out, Bandpass_filter_opts opts);

#endif // P_BANDPASS_FILTER_H
