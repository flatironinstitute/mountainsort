#ifndef MERGE_ACROSS_CHANNELS_H
#define MERGE_ACROSS_CHANNELS_H

#include <mda32.h>

struct Merge_across_channels_opts {
    int clip_size = 100;
    double min_peak_ratio_to_consider = 0.3; //reduced on 5/26/17 0.7->0.3
    double event_fraction_threshold = 0.3; //reduced on 5/26/17 0.5->0.3
};

void merge_across_channels(QVector<double>& times, QVector<int>& labels, QVector<int>& central_channels, Mda32& templates, Merge_across_channels_opts opts);

#endif // MERGE_ACROSS_CHANNELS_H
