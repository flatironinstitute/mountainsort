#ifndef CONSOLIDATE_CLUSTERS_H
#define CONSOLIDATE_CLUSTERS_H

#include <QVector>
#include "mda32.h"

struct Consolidate_clusters_opts {
    double consolidation_factor = 0.9;
};

void consolidate_clusters(QVector<bigint>& event_inds, QVector<double>& timepoints, QVector<int>& labels, const Mda32& templates, Consolidate_clusters_opts opts);
QMap<int, int> consolidate_clusters(const Mda32& templates, Consolidate_clusters_opts opts);

#endif // CONSOLIDATE_CLUSTERS_H
