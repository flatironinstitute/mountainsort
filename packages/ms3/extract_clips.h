/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef EXTRACT_CLIPS_H
#define EXTRACT_CLIPS_H

#include "mda.h"
#include "diskreadmda.h"
#include "diskreadmda32.h"

bool extract_clips(const QString& timeseries_path, const QString& firings_path, const QString& clips_path, int clip_size, const QList<int>& channels, double t1, double t2);
//bool extract_clips_features(const QString& timeseries_path, const QString& firings_path, const QString& features_path, int clip_size, int num_features);

Mda extract_clips(const DiskReadMda& X, const QVector<double>& times, int clip_size);
Mda32 extract_clips(const DiskReadMda32& X, const QVector<double>& times, int clip_size);
Mda extract_clips(const DiskReadMda& X, const QVector<double>& times, const QVector<int>& channels, int clip_size);
Mda32 extract_clips(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& channels, int clip_size);

#endif // EXTRACT_CLIPS_H
