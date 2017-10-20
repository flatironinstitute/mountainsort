#ifndef P_COMPUTE_TEMPLATES_H
#define P_COMPUTE_TEMPLATES_H

#include <QStringList>
#include "mlutil.h"

bool p_compute_templates(QStringList timeseries_list, QString firings, QString templates_out, int clip_size, const QList<int>& clusters);

#endif // P_COMPUTE_TEMPLATES_H
