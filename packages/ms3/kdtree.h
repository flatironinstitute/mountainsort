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
