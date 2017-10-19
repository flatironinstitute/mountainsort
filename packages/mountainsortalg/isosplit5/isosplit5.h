#ifndef ISOSPLIT5_H
#define ISOSPLIT5_H

//#include "mlcommon.h"
#include "isocut5.h"

struct isosplit5_opts {
    float isocut_threshold = 1.0;
    int min_cluster_size = 10;
    int K_init = 200;
    bool refine_clusters = false;
    int max_iterations_per_pass = 500;
};

bool isosplit5(int* labels_out, bigint M, bigint N, float* X, isosplit5_opts opts);

/*
 * MCWRAP [ labels_out[1,N] ] = isosplit5_mex(X[M,N])
 * SET_INPUT M = size(X,1)
 * SET_INPUT N = size(X,2)
 * SOURCES isosplit5.cpp isocut5.cpp jisotonic5.cpp
 * HEADERS isosplit5.h isocut5.h jisotonic5.h
 */
void isosplit5_mex(double* labels_out, int M, int N, double* X);

#endif // ISOSPLIT5_H
