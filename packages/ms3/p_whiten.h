#ifndef P_WHITEN_H
#define P_WHITEN_H

#include <QString>

struct Whiten_opts {
    double quantization_unit = 0;
};

bool p_whiten(QString timeseries, QString timeseries_out, Whiten_opts opts);
bool p_compute_whitening_matrix(QStringList timeseries_list, const QList<int>& channels, QString whitening_matrix_out, Whiten_opts opts);
bool p_apply_whitening_matrix(QString timeseries, QString whitening_matrix, QString timeseries_out, Whiten_opts opts);
bool p_whiten_clips(QString clips, QString whitening_matrix, QString clips_out, Whiten_opts opts);

#endif // P_WHITEN_H
