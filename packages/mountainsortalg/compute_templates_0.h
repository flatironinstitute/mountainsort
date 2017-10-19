/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 3/29/2016
*******************************************************/

#ifndef COMPUTE_TEMPLATES_0_H
#define COMPUTE_TEMPLATES_0_H

#include "diskreadmda.h"
#include "diskreadmda32.h"

Mda compute_templates_0(const DiskReadMda& X, Mda& firings, int clip_size);
Mda compute_templates_0(const DiskReadMda& X, const QVector<double>& times, const QVector<int>& labels, int clip_size);
void compute_templates_stdevs(Mda& ret_templates, Mda& ret_stdevs, DiskReadMda& X, const QVector<double>& times, const QVector<int>& labels, int clip_size);

Mda32 compute_templates_0(const DiskReadMda32& X, Mda& firings, int clip_size);
Mda32 compute_templates_0(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int clip_size);

Mda32 compute_templates_in_parallel(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int clip_size);

#endif // COMPUTE_TEMPLATES_0_H
