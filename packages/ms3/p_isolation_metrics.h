#ifndef P_ISOLATION_METRICS_H
#define P_ISOLATION_METRICS_H

#include <QStringList>

struct P_isolation_metrics_opts {
    QList<int> cluster_numbers;
    int clip_size = 50;
    int num_features = 10;
    int K_nearest = 6;
    int exhaustive_search_num = 100;
    int max_num_to_use = 500;
    int min_num_to_use = 100;
    bool do_not_compute_pair_metrics = false;
    bool compute_bursting_parents = false;
    double bursting_parent_waveform_correlation_threshold = 0.8;
    double bursting_parent_window = 50 * 30; //in timepoints
    double bursting_parent_factor = 2;
    double bursting_parent_z_threshold = 8;
};

bool p_isolation_metrics(QStringList timeseries_list, QString firings, QString metrics_out, QString pair_metrics_out, P_isolation_metrics_opts opts);

#endif // P_ISOLATION_METRICS_H
