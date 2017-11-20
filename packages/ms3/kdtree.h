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
#ifndef KDTREE_H
#define KDTREE_H

#include "mda32.h"

class KdTreePrivate;
class KdTree {
public:
    friend class KdTreePrivate;
    KdTree();
    virtual ~KdTree();
    void create(const Mda32& X);
    QList<int> allIndices() const;
    QList<int> findApproxKNearestNeighbors(const Mda32& X, const QVector<float>& p, int K, int exhaustive_search_num) const;

private:
    KdTreePrivate* d;
};

#endif // KDTREE_H
