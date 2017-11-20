/*
 * Copyright 2016-2017 Flatiron Institute, Simons Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef P_CONFUSION_MATRIX_H
#define P_CONFUSION_MATRIX_H

#include <QString>

struct P_confusion_matrix_opts {
    int max_matching_offset = 30;
    bool relabel_firings2 = false;
};

bool p_confusion_matrix(QString firings1, QString firings2, QString confusion_matrix_out, QString matched_firings_out, QString label_map_out, QString firings2_relabeled_out, QString firings2_relabel_map_out, P_confusion_matrix_opts opts);

#endif // P_CONFUSION_MATRIX_H
