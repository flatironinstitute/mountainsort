#ifndef P_CONFUSION_MATRIX_H
#define P_CONFUSION_MATRIX_H

#include <QString>

struct P_confusion_matrix_opts {
    int max_matching_offset = 30;
    bool relabel_firings2 = false;
};

bool p_confusion_matrix(QString firings1, QString firings2, QString confusion_matrix_out, QString matched_firings_out, QString label_map_out, QString firings2_relabeled_out, QString firings2_relabel_map_out, P_confusion_matrix_opts opts);

#endif // P_CONFUSION_MATRIX_H
