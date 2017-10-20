#ifndef P_EXTRACT_FIRINGS_H
#define P_EXTRACT_FIRINGS_H

#include <QString>
#include <QStringList>

struct P_extract_firings_opts {
    QStringList exclusion_tags;
    QList<int> clusters;
    double t1, t2;
};

bool p_extract_firings(QString firings, QString metrics, QString firings_out, P_extract_firings_opts opts);

#endif // P_EXTRACT_FIRINGS_H
