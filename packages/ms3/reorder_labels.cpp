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
#include "reorder_labels.h"

#include <QVector>
#include <mda32.h>
#include <cmath>
using std::sqrt;


struct template_comparer_struct {
    int channel;
    double template_peak;
    int index;
};
struct template_comparer {
    bool operator()(const template_comparer_struct& a, const template_comparer_struct& b) const
    {
        //put the zeros at the end
        if ((a.template_peak == 0) && (b.template_peak > 0))
            return false;
        if ((a.template_peak > 0) && (b.template_peak == 0))
            return true;

        if (a.channel < b.channel)
            return true;
        else if (a.channel == b.channel) {
            if (a.template_peak > b.template_peak)
                return true;
            else if (a.template_peak == b.template_peak)
                return (a.index < b.index);
            else
                return false;
        }
        else
            return false;
    }
};
template_comparer_struct compute_comparer(const Mda32& template0, int index)
{
    //QList<double> abs_peak_values;
    template_comparer_struct ret;
    int peak_channel = 1;
    double abs_peak_value = 0;
    for (int t = 0; t < template0.N2(); t++) {
        for (int m = 0; m < template0.N1(); m++) {
            double val = template0.get(m, t);
            if (fabs(val) > abs_peak_value) {
                abs_peak_value = fabs(val);
                peak_channel = m + 1;
            }
        }
    }
    ret.index = index;
    ret.channel = peak_channel;
    ret.template_peak = abs_peak_value;
    return ret;
}

//returns a label map
QMap<int, int> reorder_labels(const Mda32& templates)
{
    QList<template_comparer_struct> list;
    for (int i = 0; i < templates.N3(); i++) {
        Mda32 template0;
        templates.getChunk(template0, 0, 0, i, templates.N1(), templates.N2(), 1);
        list << compute_comparer(template0, i);
    }
    qSort(list.begin(), list.end(), template_comparer());
    //QList<int> sort_indices;
    QMap<int, int> label_map;
    for (int i = 0; i < list.count(); i++) {
        label_map[list[i].index + 1] = i + 1;
    }
    return label_map;
}
