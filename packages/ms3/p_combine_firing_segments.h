#ifndef P_COMBINE_FIRING_SEGMENTS_H
#define P_COMBINE_FIRING_SEGMENTS_H

#include <QStringList>

struct P_combine_firing_segments_opts {
    int clip_size = 60;
    double match_score_threshold = 0.6;
    int offset_search_radius = 10;
    int num_comparison_events = 1000;
};

bool p_combine_firing_segments(QString timeseries, QStringList firings_list, QString firings_out, P_combine_firing_segments_opts opts);

#endif // P_COMBINE_FIRING_SEGMENTS_H
