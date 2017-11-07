#ifndef MV_DISCRIMHIST_H
#define MV_DISCRIMHIST_H

#include <QList>
#include <QString>
#include <QVector>
#include <diskreadmda.h>
#include <diskreadmda32.h>

struct mv_discrimhist_opts {
    QList<int> clusters;
    /// TODO clip_size is hard-coded here
    int clip_size = 80;
    QString method = "centroid"; //centroid or svm
    int num_features = 0;
};

bool mv_discrimhist(QString timeseries_path, QString firings_path, QString output_path, mv_discrimhist_opts opts);
bool get_discrimhist_data(QVector<double>& ret1, QVector<double>& ret2, const DiskReadMda32& timeseries, const DiskReadMda& firings, int k1, int k2, int clip_size, QString method);

#endif // MV_DISCRIMHIST_H
