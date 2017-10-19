#ifndef P_MOUNTAINSORT3_H
#define P_MOUNTAINSORT3_H

#include <QString>
typedef int64_t bigint;

struct P_mountainsort3_opts {
    double adjacency_radius = 0;
    int clip_size = 50;

    double detect_threshold = 5;
    int detect_interval = 10;
    int detect_sign = 0;

    int num_features = 10;
    int num_features_per_channel = 10;
    bigint max_pca_samples = 10000;

    bool consolidate_clusters = true;
    double consolidation_factor = 0.9;

    bool merge_across_channels = true;

    bool fit_stage = true;

    double t1 = -1;
    double t2 = -1;
};

bool p_mountainsort3(QString timeseries, QString geom, QString firings_out, QString temp_path, const P_mountainsort3_opts& opts);

#endif // P_MOUNTAINSORT3_H
