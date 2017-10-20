#ifndef P_CONCAT_FIRINGS_H
#define P_CONCAT_FIRINGS_H

#include <QString>
#include <QStringList>

bool p_concat_firings(QStringList timeseries_list, QStringList firings_list, QString timeseries_out, QString firings_out);
bool p_concat_event_times(QStringList event_times_list, QString event_times_out);

#endif // P_CONCAT_FIRINGS_H
