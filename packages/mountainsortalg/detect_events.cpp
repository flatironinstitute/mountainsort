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
#include "detect_events.h"
#include <cmath>

using std::fabs;

/*
QVector<double> detect_events(const QVector<double>& X, double detect_threshold, double detect_interval, int sign)
{
    double mean = MLCompute::mean(X);
    double stdev = MLCompute::stdev(X);
    double threshold2 = detect_threshold * stdev;

    bigint N = X.count();
    QVector<bigint> to_use(N);
    to_use.fill(0);
    bigint last_best_ind = 0;
    double last_best_val = 0;
    for (bigint n = 0; n < N; n++) {
        double val = (X[n] - mean);
        if (sign < 0)
            val = -val;
        else if (sign == 0)
            val = fabs(val);
        if (n - last_best_ind > detect_interval)
            last_best_val = 0;
        if (val >= threshold2) {
            if (last_best_val > 0) {
                if (val > last_best_val) {
                    to_use[n] = 1;
                    to_use[last_best_ind] = 0;
                    last_best_ind = n;
                    last_best_val = val;
                }
            }
            else {
                if (val > 0) {
                    to_use[n] = 1;
                    last_best_ind = n;
                    last_best_val = val;
                }
            }
        }
    }
    QVector<double> times;
    for (bigint n = 0; n < N; n++) {
        if (to_use[n]) {
            times << n;
        }
    }
    return times;
}
*/

QVector<double> detect_events(const QVector<double>& X, double detect_threshold, double detect_interval, int sign)
{
    double mean = MLCompute::mean(X);
    double stdev = MLCompute::stdev(X);
    double threshold2 = detect_threshold * stdev;
    //double threshold2=detect_threshold;

    QVector<bigint> timepoints_to_consider;
    QVector<double> vals_to_consider;
    for (bigint n = 0; n < X.count(); n++) {
        //double val=X[n];
        double val = X[n] - mean;
        if (sign < 0)
            val = -val;
        else if (sign == 0)
            val = fabs(val);
        if (val > threshold2) {
            timepoints_to_consider << n;
            vals_to_consider << val;
        }
    }
    QVector<bool> to_use(timepoints_to_consider.count(), false);
    bigint last_best_ind = -1;
    for (bigint i = 0; i < timepoints_to_consider.count(); i++) {
        if (last_best_ind >= 0) {
            if (timepoints_to_consider[last_best_ind] < timepoints_to_consider[i] - detect_interval) {
                last_best_ind = -1;
                //last best ind is not within range. so update it
                if ((i > 0) && (timepoints_to_consider[i - 1] >= timepoints_to_consider[i] - detect_interval)) {
                    last_best_ind = i - 1;
                    bigint jj = last_best_ind;
                    while ((jj - 1 >= 0) && (timepoints_to_consider[jj - 1] >= timepoints_to_consider[i] - detect_interval)) {
                        if (vals_to_consider[jj - 1] > vals_to_consider[last_best_ind]) {
                            last_best_ind = jj - 1;
                        }
                        jj--;
                    }
                }
                else {
                    last_best_ind = -1;
                }
            }
        }
        if (last_best_ind >= 0) {
            if (vals_to_consider[i] > vals_to_consider[last_best_ind]) {
                to_use[i] = true;
                to_use[last_best_ind] = false;
                last_best_ind = i;
            }
            else {
                to_use[i] = false;
            }
        }
        else {
            to_use[i] = true;
            last_best_ind = i;
        }
    }

    QVector<double> times;
    for (bigint i = 0; i < timepoints_to_consider.count(); i++) {
        if (to_use[i]) {
            times << timepoints_to_consider[i];
        }
    }
    return times;
}
