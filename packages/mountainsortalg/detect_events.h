#ifndef DETECT_EVENTS_H
#define DETECT_EVENTS_H

#include <QVector>
#include "mlcompute.h"

QVector<double> detect_events(const QVector<double>& X, double detect_threshold, double detect_interval, int sign);

#endif // DETECT_EVENTS_H
