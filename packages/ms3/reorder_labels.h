#ifndef REORDER_LABELS_H
#define REORDER_LABELS_H

#include <QMap>
#include <mda32.h>

QMap<int, int> reorder_labels(const Mda32& templates);

#endif // REORDER_LABELS_H
