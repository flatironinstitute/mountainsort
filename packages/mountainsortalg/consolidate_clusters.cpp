#include "consolidate_clusters.h"

bool should_use_template(const Mda32& template0, Consolidate_clusters_opts opts);

void consolidate_clusters(QVector<bigint>& event_inds, QVector<double>& timepoints, QVector<int>& labels, const Mda32& templates, Consolidate_clusters_opts opts)
{
    bigint L = labels.count();
    int M = templates.N1();
    int T = templates.N2();
    int K = MLCompute::max(labels);

    QVector<int> to_use(K + 1, 0);

    for (int k = 1; k <= K; k++) {
        Mda32 template0;
        templates.getChunk(template0, 0, 0, k - 1, M, T, 1);
        if (should_use_template(template0, opts)) {
            to_use[k] = 1;
        }
    }

    QMap<int, int> label_map;
    label_map[0] = 0;
    int knext = 1;
    for (int k = 1; k <= K; k++) {
        if (to_use[k]) {
            label_map[k] = knext;
            knext++;
        }
        else {
            label_map[k] = 0;
        }
    };

    event_inds.clear();
    QVector<int> new_labels;
    QVector<double> new_timepoints;
    for (int i = 0; i < L; i++) {
        int k0 = label_map.value(labels[i]);
        if (k0 > 0) {
            event_inds << i;
            new_timepoints << timepoints[i];
            new_labels << k0;
        }
    }
    timepoints = new_timepoints;
    labels = new_labels;
}

bool should_use_template(const Mda32& template0, Consolidate_clusters_opts opts)
{
    int central_channel = 0;

    int peak_location_tolerance = 10;

    bigint M = template0.N1();
    bigint T = template0.N2();
    bigint Tmid = (int)((T + 1) / 2) - 1;

    double abs_peak_on_central_channel = 0;
    bigint abs_peak_t_on_central_channel = Tmid;
    {
        for (bigint t = 0; t < T; t++) {
            double val = template0.value(central_channel, t);
            if (fabs(val) > abs_peak_on_central_channel) {
                abs_peak_on_central_channel = fabs(val);
                abs_peak_t_on_central_channel = t;
            }
        }
    }

    double abs_peak = 0;
    //bigint abs_peak_t = Tmid;
    for (bigint t = 0; t < T; t++) {
        for (bigint m = 0; m < M; m++) {
            double val = template0.value(m, t);
            if (fabs(val) > abs_peak) {
                abs_peak = fabs(val);
                //abs_peak_t = t;
            }
        }
    }
    {
        if (abs_peak_on_central_channel < abs_peak * opts.consolidation_factor)
            return false;
        if (fabs(abs_peak_t_on_central_channel - Tmid) > peak_location_tolerance)
            return false;
    }

    return true;
}

QMap<int, int> consolidate_clusters(const Mda32& templates, Consolidate_clusters_opts opts)
{
    int M = templates.N1();
    int T = templates.N2();
    int K = templates.N3();
    QVector<int> to_use(K + 1, 0);

    for (int k = 1; k <= K; k++) {
        Mda32 template0;
        templates.getChunk(template0, 0, 0, k - 1, M, T, 1);
        if (should_use_template(template0, opts)) {
            to_use[k] = 1;
        }
    }

    QMap<int, int> label_map;
    label_map[0] = 0;
    int knext = 1;
    for (int k = 1; k <= K; k++) {
        if (to_use[k]) {
            label_map[k] = knext;
            knext++;
        }
        else {
            label_map[k] = 0;
        }
    };

    return label_map;
}
