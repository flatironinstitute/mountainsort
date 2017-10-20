#ifndef P_BANJOVIEW_CROSS_CORRELOGRAMS_H
#define P_BANJOVIEW_CROSS_CORRELOGRAMS_H

#include <QList>
#include <QString>

enum P_banjoview_cross_correlograms_mode {
    Autocorrelograms,
    Matrix_of_cross_correlograms
};

struct P_banjoview_cross_correlograms_opts {
    P_banjoview_cross_correlograms_mode mode;
    QList<int> clusters;
    double max_dt_msec = 100;
    double bin_size_msec = 1;
    double samplerate = 30000;
};

bool p_banjoview_cross_correlograms(QString firings, QString correlograms_out, P_banjoview_cross_correlograms_opts opts);

#endif // P_BANJOVIEW_CROSS_CORRELOGRAMS_H
