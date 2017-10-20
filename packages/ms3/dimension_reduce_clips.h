#ifndef DIMENSION_REDUCE_CLIPS_H
#define DIMENSION_REDUCE_CLIPS_H

#include "mda32.h"

void dimension_reduce_clips(Mda32& reduced_clips, const Mda32& clips, bigint num_features_per_channel, bigint max_samples);
void dimension_reduce_clips(QString clips_path, QString reduced_clips_path, bigint num_features_per_channel, bigint max_samples);

#endif // DIMENSION_REDUCE_CLIPS_H
