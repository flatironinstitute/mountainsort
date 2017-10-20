#ifndef P_SPLIT_FIRINGS_H
#define P_SPLIT_FIRINGS_H

#include <QStringList>

bool p_split_firings(QStringList timeseries_list, QString firings, QStringList firings_out_list);
bool p_extract_firings(QString firings, const QSet<int>& clusters, QString firings_out);

#endif // P_SPLIT_FIRINGS_H
