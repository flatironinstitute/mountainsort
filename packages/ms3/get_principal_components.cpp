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
#include "get_principal_components.h"
#include <QDebug>
#include <math.h>
#include "get_sort_indices.h"
#include <stdlib.h>
#include <iostream>
#include <QTime>

/*
void make_random_vector(int M, double* v)
{
    for (int i = 0; i < M; i++)
        v[i] = (qrand() * 1.0 / RAND_MAX) * 2 - 1;
}

void XXT_vector_mult(int M, int N, double* A, double* v)
{
    double* v2 = (double*)malloc(sizeof(double) * N);
    for (int n = 0; n < N; n++)
        v2[n] = 0;
    int ii = 0;
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            v2[n] += A[ii] * v[m];
            ii++;
        }
    }
    ii = 0;
    for (int m = 0; m < M; m++)
        v[m] = 0;
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            v[m] += A[ii] * v2[n];
            ii++;
        }
    }
    free(v2);
}

void normalize_vector(int M, double* v)
{
    double sumsqr = 0;
    for (int m = 0; m < M; m++)
        sumsqr += v[m] * v[m];
    if (sumsqr > 0) {
        double norm = sqrt(sumsqr);
        for (int m = 0; m < M; m++)
            v[m] /= norm;
    }
}

void get_top_component(int M, int N, double* comp, double* data, int num_iterations = 30)
{
    make_random_vector(M, comp);
    for (int ii = 0; ii < num_iterations; ii++) {
        XXT_vector_mult(M, N, data, comp);
        normalize_vector(M, comp);
    }
}

void subtract_component_from_vector(int M, double* v, double* data)
{
    double ip = 0;
    for (int m = 0; m < M; m++)
        ip += v[m] * data[m];
    for (int m = 0; m < M; m++)
        data[m] -= ip * v[m];
}

void subtract_component_from_data(int M, int N, double* v, double* data)
{
    for (int n = 0; n < N; n++) {
        subtract_component_from_vector(M, v, &data[M * n]);
    }
}

double compute_inner_product(int M, double* v1, double* v2)
{
    double ip = 0;
    for (int m = 0; m < M; m++)
        ip += v1[m] * v2[m];
    return ip;
}

void do_pca(int M, int N, int ncomp, double* components, double* features, double* data)
{
    double* working_data = (double*)malloc(sizeof(double) * M * N);
    for (int ii = 0; ii < M * N; ii++)
        working_data[ii] = data[ii];
    double* v = (double*)malloc(sizeof(double) * M);
    for (int k = 0; k < ncomp; k++) {
        get_top_component(M, N, v, working_data);
        subtract_component_from_data(M, N, v, working_data);
        memcpy(&components[M * k], v, M * sizeof(double));
    }
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < ncomp; k++) {
            features[k + ncomp * n] = compute_inner_product(M, &components[M * k], &data[M * n]);
        }
    }
}

void get_pca_components(int M, int N, int ncomp, float* components, float* data)
{
    double* data2 = (double*)malloc(sizeof(double) * M * N);
    double* features2 = (double*)malloc(sizeof(double) * N * ncomp);
    double* components2 = (double*)malloc(sizeof(double) * M * ncomp);
    for (int ii = 0; ii < M * N; ii++)
        data2[ii] = data[ii];
    do_pca(M, N, ncomp, components2, features2, data2);
    for (int ii = 0; ii < M * ncomp; ii++)
        components[ii] = components2[ii];
    free(data2);
    free(features2);
    free(components2);
}

void get_pca_features_2(int M, int N, int ncomp, float* features, float* data)
{
    double* data2 = (double*)malloc(sizeof(double) * M * N);
    double* features2 = (double*)malloc(sizeof(double) * N * ncomp);
    double* components2 = (double*)malloc(sizeof(double) * M * ncomp);
    for (int ii = 0; ii < M * N; ii++)
        data2[ii] = data[ii];
    do_pca(M, N, ncomp, components2, features2, data2);
    for (int ii = 0; ii < N * ncomp; ii++)
        features[ii] = features2[ii];
    free(data2);
    free(features2);
    free(components2);
}
void get_pca_features_2(int M, int N, int ncomp, double* features, double* data)
{
    double* components2 = (double*)malloc(sizeof(double) * M * ncomp);
    do_pca(M, N, ncomp, components2, features, data);
    free(components2);
}
*/
