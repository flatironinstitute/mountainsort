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
#include "kdtree.h"
#include "mlutil.h"
#include "get_sort_indices.h"

class KdTreePrivate {
public:
    KdTree* q;

    QList<int> m_leaf_indices;
    KdTree* m_left_tree = 0;
    KdTree* m_center_tree = 0;
    KdTree* m_right_tree = 0;
    float m_cutoff_1 = 0;
    float m_cutoff_2 = 0;
    QVector<float> m_projection_direction;
    int m_num_datapoints = 0;

    void create(const Mda32& X, const QList<int>& indices);
    static QVector<float> get_projection_direction(const Mda32& X);
    static double compute_distsqr(int M, const float* x, const float* y);
    static QList<int> find_closest_K_candidates(const Mda32& X, const QVector<float>& p, int K, const QList<int>& indices);
};

KdTree::KdTree()
{
    d = new KdTreePrivate;
    d->q = this;
}

KdTree::~KdTree()
{
    if (d->m_left_tree)
        delete d->m_left_tree;
    if (d->m_center_tree)
        delete d->m_center_tree;
    if (d->m_right_tree)
        delete d->m_right_tree;
    delete d;
}

void KdTree::create(const Mda32& X)
{
    QList<int> indices;
    for (int i = 0; i < X.N2(); i++)
        indices << i;
    d->create(X, indices);
}

QList<int> KdTree::allIndices() const
{
    if ((!d->m_left_tree) && (!d->m_center_tree) && (!d->m_right_tree)) {
        //leaf node
        return d->m_leaf_indices;
    }
    else {
        QList<int> ret;
        if (d->m_left_tree)
            ret.append(d->m_left_tree->allIndices());
        if (d->m_center_tree)
            ret.append(d->m_center_tree->allIndices());
        if (d->m_right_tree)
            ret.append(d->m_right_tree->allIndices());
        return ret;
    }
}

void KdTreePrivate::create(const Mda32& X, const QList<int>& indices)
{
    if (indices.isEmpty())
        return;

    m_num_datapoints = indices.count();

    int M = X.N1();
    //define the projection direction
    m_projection_direction = get_projection_direction(X);
    //compute the inner product of the points with the projection direction
    const float* ptr = X.constDataPtr();
    QVector<float> vals(indices.count());
    for (int i = 0; i < indices.count(); i++) {
        vals[i] = MLCompute::dotProduct(M, m_projection_direction.data(), &ptr[indices[i] * M]);
    }
    //sort the values and create the cutoff as the median value
    QVector<float> vals_sorted = vals;
    qSort(vals_sorted);
    m_cutoff_1 = vals_sorted.value(vals_sorted.count() * 3 / 7);
    m_cutoff_2 = vals_sorted.value(vals_sorted.count() * 4 / 7);
    //decide which points go to the left and right
    QList<int> indices_left;
    QList<int> indices_center;
    QList<int> indices_right;
    for (int i = 0; i < indices.count(); i++) {
        if (vals.value(i) < m_cutoff_1)
            indices_left << indices[i];
        else if (vals.value(i) < m_cutoff_2)
            indices_center << indices[i];
        else
            indices_right << indices[i];
    }
    //if one of these is the empty set, then we are done -- effectively a leaf node
    if ((indices_left.isEmpty()) || (indices_center.isEmpty()) || (indices_right.isEmpty())) {
        m_leaf_indices = indices;
        return;
    }
    //create the left and right trees
    m_left_tree = new KdTree;
    m_left_tree->d->create(X, indices_left);
    m_center_tree = new KdTree;
    m_center_tree->d->create(X, indices_center);
    m_right_tree = new KdTree;
    m_right_tree->d->create(X, indices_right);
}

QList<int> KdTreePrivate::find_closest_K_candidates(const Mda32& X, const QVector<float>& p, int K, const QList<int>& indices)
{
    int M = X.N1();
    const float* ptr = X.constDataPtr();
    int num = indices.count(); //should be same as m_num_datapoints
    QVector<double> distsqrs(num);
    for (int i = 0; i < num; i++) {
        distsqrs[i] = compute_distsqr(M, p.data(), &ptr[indices[i] * M]);
    }
    QList<int> inds = get_sort_indices(distsqrs);
    QList<int> ret;
    for (int i = 0; (i < inds.count()) && (i < K); i++) {
        ret << indices[inds[i]];
    }
    return ret;
}

QList<int> KdTree::findApproxKNearestNeighbors(const Mda32& X, const QVector<float>& p, int K, int exhaustive_search_num) const
{
    int M = X.N1();
    //const float *ptr=X.constDataPtr();
    if ((d->m_num_datapoints <= exhaustive_search_num) || (!d->m_left_tree) || (!d->m_right_tree) || (!d->m_center_tree)) {
        //we are now going to do an exaustive search
        QList<int> indices = allIndices();
        return d->find_closest_K_candidates(X, p, K, indices);
    }
    else {
        float val = MLCompute::dotProduct(M, p.data(), d->m_projection_direction.data());
        QList<int> candidates;
        if (val < d->m_cutoff_1) {
            candidates.append(d->m_left_tree->findApproxKNearestNeighbors(X, p, K, exhaustive_search_num));
            candidates.append(d->m_center_tree->findApproxKNearestNeighbors(X, p, K, exhaustive_search_num));
        }
        else if (val < d->m_cutoff_2) {
            candidates.append(d->m_left_tree->findApproxKNearestNeighbors(X, p, K, exhaustive_search_num));
            candidates.append(d->m_center_tree->findApproxKNearestNeighbors(X, p, K, exhaustive_search_num));
            candidates.append(d->m_right_tree->findApproxKNearestNeighbors(X, p, K, exhaustive_search_num));
        }
        else {
            candidates.append(d->m_right_tree->findApproxKNearestNeighbors(X, p, K, exhaustive_search_num));
            candidates.append(d->m_right_tree->findApproxKNearestNeighbors(X, p, K, exhaustive_search_num));
        }
        return d->find_closest_K_candidates(X, p, K, candidates);
    }
}

QVector<float> KdTreePrivate::get_projection_direction(const Mda32& X)
{
    int M = X.N1();
    QVector<float> ret(M);
    for (int i = 0; i < M; i++) {
        double randval = (qrand() % RAND_MAX) * 1.0 / RAND_MAX;
        ret[i] = randval * 2 - 1;
    }
    return ret;
    /*
    Mda32 CC, FF, sig;
    pca(CC, FF, sig, X, 1, false);
    QVector<float> ret;
    for (int i = 0; i < X.N1(); i++)
        ret << CC.value(i);
    return ret;
    */
}
double KdTreePrivate::compute_distsqr(int M, const float* x, const float* y)
{
    double ret = 0;
    for (int i = 0; i < M; i++) {
        ret += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return ret;
}
