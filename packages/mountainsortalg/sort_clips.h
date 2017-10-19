#ifndef SORT_CLIPS_H
#define SORT_CLIPS_H

#include <mda32.h>

struct Sort_clips_opts {
    int num_features = 10;
    bigint max_samples = 10000;
    double isocut_threshold = 1;
    int K_init = 200;
};

QVector<int> sort_clips(const Mda32& clips, const Sort_clips_opts& opts);

#endif // SORT_CLIPS_H
