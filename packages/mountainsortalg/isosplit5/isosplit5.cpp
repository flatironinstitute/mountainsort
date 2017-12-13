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
#include "isosplit5.h"
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <math.h>
#include "isocut5.h"

typedef std::vector<std::vector<bigint> > intarray2d;
void alloc(intarray2d& X, bigint N1, bigint N2)
{
    X.resize(N1);
    for (bigint i = 0; i < N1; i++) {
        X[i].resize(N2);
    }
}

namespace ns_isosplit5 {
struct kmeans_opts {
    bigint num_iterations = 0;
};

bigint compute_max(bigint N, int* labels);
bigint compute_max(bigint N, bigint* inds);
void kmeans_multistep(int* labels, bigint M, bigint N, float* X, bigint K1, bigint K2, bigint K3, kmeans_opts opts);
void kmeans_maxsize(int* labels, bigint M, bigint N, float* X, bigint maxsize, kmeans_opts opts);
void compare_clusters(double* dip_score, std::vector<bigint>* new_labels1, std::vector<bigint>* new_labels2, bigint M, bigint N1, bigint N2, float* X1, float* X2, float* centroid1, float* centroid2);
void compute_centroids(float* centroids, bigint M, bigint N, bigint Kmax, float* X, int* labels, std::vector<bigint>& clusters_to_compute_vec);
void compute_covmats(float* covmats, bigint M, bigint N, bigint Kmax, float* X, int* labels, float* centroids, std::vector<bigint>& clusters_to_compute_vec);
void get_pairs_to_compare(std::vector<bigint>* inds1, std::vector<bigint>* inds2, bigint M, bigint K, float* active_centroids, const intarray2d& active_comparisons_made);
void compare_pairs(std::vector<bigint>* clusters_changed, bigint* total_num_label_changes, bigint M, bigint N, float* X, int* labels, const std::vector<bigint>& inds1, const std::vector<bigint>& inds2, const isosplit5_opts& opts, float* centroids, float* covmats); //the labels are updated
}

namespace smi {
bool get_inverse_via_lu_decomposition(int M, float* out, float* in);
}

void isosplit5_mex(double* labels_out, int M, int N, double* X)
{
    float* Xf = (float*)malloc(sizeof(float) * M * N);
    int* labelsi = (int*)malloc(sizeof(int) * N);
    for (bigint i = 0; i < M * N; i++)
        Xf[i] = X[i];
    isosplit5_opts opts;
    //opts.refine_clusters=true;
    isosplit5(labelsi, M, N, Xf, opts);
    for (bigint i = 0; i < N; i++)
        labels_out[i] = labelsi[i];
    free(Xf);
    free(labelsi);
}

struct parcelate2_opts {
    bool final_reassign = false; //not yet implemented
};

struct p2_parcel {
    std::vector<bigint> indices;
    std::vector<float> centroid;
    double radius;
};

void print_matrix(bigint M, bigint N, float* A)
{
    for (bigint m = 0; m < M; m++) {
        for (bigint n = 0; n < N; n++) {
            float val = A[m + M * n];
            printf("%g ", val);
        }
        printf("\n");
    }
}

std::vector<float> p2_compute_centroid(bigint M, float* X, const std::vector<bigint>& indices)
{
    std::vector<double> ret(M);
    double count = 0;
    for (bigint m = 0; m < M; m++) {
        ret[m] = 0;
    }
    for (bigint i = 0; i < (bigint)indices.size(); i++) {
        for (bigint m = 0; m < M; m++) {
            ret[m] += X[m + M * indices[i]];
        }
        count++;
    }
    if (count) {
        for (bigint m = 0; m < M; m++) {
            ret[m] /= count;
        }
    }
    std::vector<float> retf(M);
    for (bigint m = 0; m < M; m++)
        retf[m] = ret[m];
    return retf;
}

double p2_compute_max_distance(const std::vector<float>& centroid, bigint M, float* X, const std::vector<bigint>& indices)
{
    double max_dist = 0;
    for (bigint i = 0; i < (bigint)indices.size(); i++) {
        double dist = 0;
        for (bigint m = 0; m < M; m++) {
            double val = centroid[m] - X[m + M * indices[i]];
            dist += val * val;
        }
        dist = sqrt(dist);
        if (dist > max_dist)
            max_dist = dist;
    }
    return max_dist;
}

std::vector<bigint> p2_randsample(bigint N, bigint K)
{
    (void)N;
    // Not we are not actually randomizing here. There's a reason, I believe.
    std::vector<bigint> inds;
    for (bigint a = 0; a < K; a++)
        inds.push_back(a);
    return inds;
    /*
    if (K>N) K=N;std::vector<bigint> inds;
    std::vector<bigint> used(N);
    for (bigint i=0; i<N; i++)
        used[i]=0;
    for (bigint k=0; k<K; k++)
        used[k]=1;
    for (bigint k=0; k<K; k++) {
        bigint ii=rand()%N;
        bigint tmp=used[k];
        used[k]=used[ii];
        used[ii]=tmp;
    }
    std::vector<bigint> inds;
    for (bigint i=0; i<N; i++) {
        if (used[i])
            inds.push_back(i);
    }
    return inds;
    */
}

bool parcelate2(int* labels, bigint M, bigint N, float* X, bigint target_parcel_size, bigint target_num_parcels, const parcelate2_opts& p2opts)
{
    std::vector<p2_parcel> parcels;

    for (bigint i = 0; i < N; i++)
        labels[i] = 1;

    p2_parcel P;
    P.indices.resize(N);
    for (bigint i = 0; i < N; i++)
        P.indices[i] = i;
    P.centroid = p2_compute_centroid(M, X, P.indices);
    P.radius = p2_compute_max_distance(P.centroid, M, X, P.indices);
    parcels.push_back(P);

    bigint split_factor = 3; // split factor around 2.71 is in a sense ideal

    double target_radius;
    while ((bigint)parcels.size() < target_num_parcels) {
        bool candidate_found = false;
        for (bigint i = 0; i < (bigint)parcels.size(); i++) {
            std::vector<bigint>* indices = &parcels[i].indices;
            if ((bigint)indices->size() > target_parcel_size) {
                if (parcels[i].radius > 0)
                    candidate_found = true;
            }
        }
        if (!candidate_found) {
            // nothing else will ever be split
            break;
        }

        target_radius = 0;
        for (bigint i = 0; i < (bigint)parcels.size(); i++) {
            if ((bigint)parcels[i].indices.size() > target_parcel_size) {
                double tmp = parcels[i].radius * 0.95;
                if (tmp > target_radius)
                    target_radius = tmp;
            }
        }
        if (target_radius == 0) {
            printf("Unexpected target radius of zero.\n");
            break;
        }

        bigint p_index = 0;
        while (p_index < (bigint)parcels.size()) {
            std::vector<bigint> inds = parcels[p_index].indices;
            double rad = parcels[p_index].radius;
            bigint sz = parcels[p_index].indices.size();
            if ((sz > target_parcel_size) && (rad >= target_radius)) {
                std::vector<bigint> assignments(inds.size());
                std::vector<bigint> iii = p2_randsample(sz, split_factor);
                // iii[1] = iii[0]; //force failure for testing debug
                for (bigint i = 0; i < (bigint)inds.size(); i++) {
                    bigint best_pt = -1;
                    double best_dist = 0;
                    for (bigint j = 0; j < (bigint)iii.size(); j++) {
                        double dist = 0;
                        for (bigint m = 0; m < M; m++) {
                            double val = X[m + M * inds[iii[j]]] - X[m + M * inds[i]];
                            dist += val * val;
                        }
                        dist = sqrt(dist);
                        if ((best_pt < 0) || (dist < best_dist)) {
                            best_dist = dist;
                            best_pt = j;
                        }
                    }
                    assignments[i] = best_pt;
                }
                parcels[p_index].indices.clear();
                for (bigint i = 0; i < (bigint)inds.size(); i++) {
                    if (assignments[i] == 0) {
                        parcels[p_index].indices.push_back(inds[i]);
                        labels[inds[i]] = p_index + 1;
                    }
                }
                parcels[p_index].centroid = p2_compute_centroid(M, X, parcels[p_index].indices);
                parcels[p_index].radius = p2_compute_max_distance(parcels[p_index].centroid, M, X, parcels[p_index].indices);
                for (bigint jj = 1; jj < (bigint)iii.size(); jj++) {
                    p2_parcel PP;
                    for (bigint i = 0; i < (bigint)inds.size(); i++) {
                        if (assignments[i] == jj) {
                            PP.indices.push_back(inds[i]);
                            labels[inds[i]] = parcels.size() + 1;
                        }
                    }
                    PP.centroid = p2_compute_centroid(M, X, PP.indices);
                    PP.radius = p2_compute_max_distance(PP.centroid, M, X, PP.indices);
                    if (PP.indices.size() > 0)
                        parcels.push_back(PP);
                    else {
                        for (int aa=0; aa<split_factor; aa++) {
                            printf("DEBUG: iii[%d]=%ld: ",aa,iii[aa]);
                            for (int mm=0; mm<M; mm++) {
                                printf("%g ",X[mm+M*inds[iii[aa]]]);
                            }
                            printf("\n");
                        }
                        printf("Unexpected problem. New parcel has no points -- perhaps dataset contains duplicate points? -- original size = %ld.\n", sz);
                        return false;
                    }
                }
                if ((bigint)parcels[p_index].indices.size() == sz) {
                    printf("Warning: Size did not change after splitting parcel.\n");
                    p_index++;
                }
            }
            else {
                p_index++;
            }
        }
    }

    //final reassign not yet implemented
    if (p2opts.final_reassign) {
        //centroids=get_parcel_centroids(parcels);
        //labels=knnsearch(centroids',X','K',1)';
    }

    return true;
}

bool isosplit5(int* labels, bigint M, bigint N, float* X, isosplit5_opts opts)
{

    // compute the initial clusters
    bigint target_parcel_size = opts.min_cluster_size;
    bigint target_num_parcels = opts.K_init;
    // !! important not to do a final reassign because then the shapes will not be conducive to isosplit iterations -- hexagons are not good for isosplit!
    parcelate2_opts p2opts;
    p2opts.final_reassign = false;
    if (!parcelate2(labels, M, N, X, target_parcel_size, target_num_parcels, p2opts))
        return false;
    int Kmax = ns_isosplit5::compute_max(N, labels);

    float* centroids = (float*)malloc(sizeof(float) * M * Kmax);
    float* covmats = (float*)malloc(sizeof(float) * M * M * Kmax);
    std::vector<bigint> clusters_to_compute_vec;
    for (bigint k = 0; k < Kmax; k++)
        clusters_to_compute_vec.push_back(1);
    ns_isosplit5::compute_centroids(centroids, M, N, Kmax, X, labels, clusters_to_compute_vec);
    ns_isosplit5::compute_covmats(covmats, M, N, Kmax, X, labels, centroids, clusters_to_compute_vec);

    // The active labels are those that are still being used -- for now, everything is active
    int active_labels_vec[Kmax];
    for (bigint i = 0; i < Kmax; i++)
        active_labels_vec[i] = 1;
    std::vector<int> active_labels;
    for (bigint i = 0; i < Kmax; i++)
        active_labels.push_back(i + 1);

    // Repeat while something has been merged in the pass
    bool final_pass = false; // plus we do one final pass at the end
    intarray2d comparisons_made; // Keep a matrix of comparisons that have been made in this pass
    alloc(comparisons_made, Kmax, Kmax);
    for (bigint i1 = 0; i1 < Kmax; i1++)
        for (bigint i2 = 0; i2 < Kmax; i2++)
            comparisons_made[i1][i2] = 0;
    while (true) { //passes
        bool something_merged = false; //Keep track of whether something has merged in this pass. If not, do a final pass.
        std::vector<bigint> clusters_changed_vec_in_pass(Kmax); //Keep track of the clusters that have changed in this pass so that we can update the comparisons_made matrix at the end
        for (bigint i = 0; i < Kmax; i++)
            clusters_changed_vec_in_pass[i] = 0;
        bigint iteration_number = 0;
        while (true) { //iterations

            std::vector<bigint> clusters_changed_vec_in_iteration(Kmax); //Keep track of the clusters that have changed in this iteration so that we can update centroids and covmats
            for (bigint i = 0; i < Kmax; i++)
                clusters_changed_vec_in_iteration[i] = 0;

            iteration_number++;
            if (iteration_number > opts.max_iterations_per_pass) {
                printf("Warning: max iterations per pass exceeded.\n");
                break;
            }

            if (active_labels.size() > 0) {
                // Create an array of active centroids and comparisons made, for determining the pairs to compare
                float* active_centroids = (float*)malloc(sizeof(float) * M * active_labels.size());
                for (bigint i = 0; i < (bigint)active_labels.size(); i++) {
                    for (bigint m = 0; m < M; m++) {
                        active_centroids[m + M * i] = centroids[m + M * (active_labels[i] - 1)];
                    }
                }
                intarray2d active_comparisons_made;
                alloc(active_comparisons_made, active_labels.size(), active_labels.size());
                for (bigint i1 = 0; i1 < (bigint)active_labels.size(); i1++) {
                    for (bigint i2 = 0; i2 < (bigint)active_labels.size(); i2++) {
                        active_comparisons_made[i1][i2] = comparisons_made[active_labels[i1] - 1][active_labels[i2] - 1];
                    }
                }

                // Find the pairs to compare on this iteration
                // These will be closest pairs of active clusters that have not yet
                // been compared in this pass
                std::vector<bigint> inds1, inds2;
                ns_isosplit5::get_pairs_to_compare(&inds1, &inds2, M, active_labels.size(), active_centroids, active_comparisons_made);
                std::vector<bigint> inds1b, inds2b; //remap the clusters to the original labeling
                for (bigint i = 0; i < (bigint)inds1.size(); i++) {
                    inds1b.push_back(active_labels[inds1[i] - 1]);
                    inds2b.push_back(active_labels[inds2[i] - 1]);
                }

                // If we didn't find any, break from this iteration
                if (inds1b.size() == 0) {
                    break;
                }

                // Actually compare the pairs -- in principle this operation could be parallelized
                std::vector<bigint> clusters_changed;
                bigint total_num_label_changes = 0;
                ns_isosplit5::compare_pairs(&clusters_changed, &total_num_label_changes, M, N, X, labels, inds1b, inds2b, opts, centroids, covmats); //the labels are updated
                for (bigint i = 0; i < (bigint)clusters_changed.size(); i++) {
                    clusters_changed_vec_in_pass[clusters_changed[i] - 1] = 1;
                    clusters_changed_vec_in_iteration[clusters_changed[i] - 1] = 1;
                }

                // Update which comparisons have been made
                for (bigint j = 0; j < (bigint)inds1b.size(); j++) {
                    comparisons_made[inds1b[j] - 1][inds2b[j] - 1] = 1;
                    comparisons_made[inds2b[j] - 1][inds1b[j] - 1] = 1;
                }

                // Recompute the centers for those that have changed in this iteration
                ns_isosplit5::compute_centroids(centroids, M, N, Kmax, X, labels, clusters_changed_vec_in_iteration);
                ns_isosplit5::compute_covmats(covmats, M, N, Kmax, X, labels, centroids, clusters_changed_vec_in_iteration);

                // For diagnostics
                //printf ("total num label changes = %d\n",total_num_label_changes);

                // Determine whether something has merged and update the active labels
                for (bigint i = 0; i < Kmax; i++)
                    active_labels_vec[i] = 0;
                for (bigint i = 0; i < N; i++)
                    active_labels_vec[labels[i] - 1] = 1;
                std::vector<int> new_active_labels;
                for (bigint i = 0; i < Kmax; i++)
                    if (active_labels_vec[i])
                        new_active_labels.push_back(i + 1);
                if (new_active_labels.size() < active_labels.size())
                    something_merged = true;
                active_labels = new_active_labels;

                free(active_centroids);
            }
        }

        // zero out the comparisons made matrix only for those that have changed in this pass
        for (bigint i = 0; i < Kmax; i++) {
            if (clusters_changed_vec_in_pass[i]) {
                for (bigint j = 0; j < Kmax; j++) {
                    comparisons_made[i][j] = 0;
                    comparisons_made[j][i] = 0;
                }
            }
        }

        if (something_merged)
            final_pass = false;
        if (final_pass)
            break; // This was the final pass and nothing has merged
        if (!something_merged)
            final_pass = true; // If we are done, do one last pass for final redistributes
    }

    // We should remap the labels to occupy the first natural numbers
    bigint labels_map[Kmax];
    for (bigint i = 0; i < Kmax; i++)
        labels_map[i] = 0;
    for (bigint i = 0; i < (bigint)active_labels.size(); i++) {
        labels_map[active_labels[i] - 1] = i + 1;
    }
    for (bigint i = 0; i < N; i++) {
        labels[i] = labels_map[labels[i] - 1];
    }

    // If the user wants to refine the clusters, then we repeat isosplit on each
    // of the new clusters, recursively. Unless we only found only one cluster.
    bigint K = ns_isosplit5::compute_max(N, labels);

    if ((opts.refine_clusters) && (K > 1)) {
        int* labels_split = (int*)malloc(sizeof(int) * N);
        isosplit5_opts opts2 = opts;
        opts2.refine_clusters = true; // Maybe we should provide an option on whether to do recursive refinement
        bigint k_offset = 0;
        for (bigint k = 1; k <= K; k++) {
            std::vector<bigint> inds_k;
            for (bigint i = 0; i < N; i++)
                if (labels[i] == k)
                    inds_k.push_back(i);
            if (inds_k.size() > 0) {
                float* X_k = (float*)malloc(sizeof(float) * M * inds_k.size()); //Warning: this may cause memory problems -- especially for recursive case
                int* labels_k = (int*)malloc(sizeof(int) * inds_k.size());
                for (bigint i = 0; i < (bigint)inds_k.size(); i++) {
                    for (bigint m = 0; m < M; m++) {
                        X_k[m + M * i] = X[m + M * inds_k[i]];
                    }
                }
                isosplit5(labels_k, M, inds_k.size(), X_k, opts2);
                for (bigint i = 0; i < (bigint)inds_k.size(); i++) {
                    labels_split[inds_k[i]] = k_offset + labels_k[i];
                }
                k_offset += ns_isosplit5::compute_max(inds_k.size(), labels_k);
                free(labels_k);
                free(X_k);
            }
        }
        for (bigint i = 0; i < N; i++)
            labels[i] = labels_split[i];
        free(labels_split);
    }

    free(centroids);
    free(covmats);

    return true;
}

/*



*/

/*
void isosplit5_old(bigint *labels_out,bigint M, bigint N,float *X,isosplit5_opts opts) {
    for (bigint i=0; i<N; i++) {
        labels_out[i]=1;
    }

    isosplit5_data DD(M,N,X);

    DD.initialize_labels();
    DD.compute_all_centroids();

    bigint max_iterations=500;
    bigint max_iterations_without_merges=5;

    bigint iteration_number=1;
    bigint num_iterations_without_merges=0;
    while (true) {
        iteration_number++;
        if (iteration_number>max_iterations) {
            printf ("isosplit5: Exceeded maximum number of iterations. Breaking.");
            break;
        }

        printf ("Number of active labels: %d\n",DD.get_active_labels().size());

        std::vector<bigint> k1s,k2s;
        DD.get_pairs_to_compare(&k1s,&k2s);

        printf ("compare %d pairs\n",k1s.size());
        std::vector<bigint> old_active_labels=DD.get_active_labels();
        bigint num_changes=DD.compare_pairs(k1s,k2s,opts.isocut_threshold);
        std::vector<bigint> new_active_labels=DD.get_active_labels();
        printf ("  %d changes\n",num_changes);

        if (new_active_labels.size()==old_active_labels.size())
            num_iterations_without_merges++;
        else
            num_iterations_without_merges=0;

        if (num_iterations_without_merges>=max_iterations_without_merges)
            break;
    }

    for (bigint pass=1; pass<=2; pass++) {
        std::vector<bigint> active_labels=DD.get_active_labels();
        for (bigint i1=0; i1<(bigint)active_labels.size(); i1++) {
            for (bigint i2=i1+1; i2<(bigint)active_labels.size(); i2++) {
                bigint k1=active_labels[i1];
                bigint k2=active_labels[i2];
                if ((DD.active_labels_vec[k1-1])&&(DD.active_labels_vec[k2-1])) {
                    printf ("Number of active labels: %d\n",DD.get_active_labels().size());
                    printf ("compare %d/%d (pass %d)\n",k1,k2,pass);
                    std::vector<bigint> k1s,k2s;
                    k1s.push_back(k1);
                    k2s.push_back(k2);
                    bigint num_changes=DD.compare_pairs(k1s,k2s,opts.isocut_threshold);
                    printf ("  %d changes\n",num_changes);
                }
            }
        }
    }

    std::vector<bigint> active_labels=DD.get_active_labels();
    std::vector<bigint> labels_map(ns_isosplit5::compute_max(N,DD.labels)+1);
    for (bigint i=0; i<(bigint)active_labels.size(); i++) {
        labels_map[active_labels[i]]=i+1;
    }
    for (bigint i=0; i<N; i++) {
        labels_out[i]=labels_map[DD.labels[i]];
    }
}
*/

namespace ns_isosplit5 {
bigint compute_max(bigint N, int* labels)
{
    if (N == 0)
        return 0;
    bigint ret = labels[0];
    for (bigint i = 0; i < N; i++) {
        if (labels[i] > ret)
            ret = labels[i];
    }
    return ret;
}

/*
bigint compute_max(bigint N, bigint* inds)
{
    if (N == 0)
        return 0;
    bigint ret = inds[0];
    for (bigint i = 0; i < N; i++) {
        if (inds[i] > ret)
            ret = inds[i];
    }
    return ret;
}
*/

void kmeans_initialize(double* centroids, bigint M, bigint N, bigint K, float* X)
{
    std::vector<bigint> used(N);
    for (bigint i = 0; i < N; i++)
        used[i] = 0;
    for (bigint k = 0; k < K; k++)
        used[k] = 1;
    for (bigint k = 0; k < K; k++) {
        bigint ii = rand() % N;
        bigint tmp = used[k];
        used[k] = used[ii];
        used[ii] = tmp;
    }
    std::vector<bigint> inds;
    for (bigint i = 0; i < N; i++) {
        if (used[i])
            inds.push_back(i);
    }
    for (bigint k = 0; k < (bigint)inds.size(); k++) {
        for (bigint m = 0; m < M; m++) {
            centroids[m + M * k] = X[m + M * inds[k]];
        }
    }
}
double compute_dist(bigint M, float* X, double* Y)
{
    double sumsqr = 0;
    for (bigint m = 0; m < M; m++) {
        double val = X[m] - Y[m];
        sumsqr += val * val;
    }
    return sqrt(sumsqr);
}

bigint kmeans_assign2(bigint M, bigint K, float* X0, double* centroids)
{
    bigint ret = 0;
    double best_dist = 0;
    for (bigint k = 1; k <= K; k++) {
        double dist = compute_dist(M, X0, &centroids[M * (k - 1)]);
        if ((ret == 0) || (dist < best_dist)) {
            best_dist = dist;
            ret = k;
        }
    }
    return ret;
}
void kmeans_assign(int* labels, bigint M, bigint N, bigint K, float* X, double* centroids)
{
    for (bigint i = 0; i < N; i++) {
        labels[i] = kmeans_assign2(M, K, &X[M * i], centroids);
    }
}
void kmeans_centroids(double* centroids, bigint M, bigint N, bigint K, float* X, int* labels)
{
    std::vector<bigint> counts(K);
    for (bigint k = 1; k <= K; k++) {
        counts[k - 1] = 0;
        for (bigint m = 0; m < M; m++) {
            centroids[m + (k - 1) * M] = 0;
        }
    }
    for (bigint i = 0; i < N; i++) {
        bigint k = labels[i];
        for (bigint m = 0; m < M; m++) {
            centroids[m + (k - 1) * M] += X[m + i * M];
        }
        counts[k - 1]++;
    }
    for (bigint k = 1; k <= K; k++) {
        if (counts[k - 1]) {
            for (bigint m = 0; m < M; m++) {
                centroids[m + (k - 1) * M] /= counts[k - 1];
            }
        }
    }
}

void kmeans(int* labels, bigint M, bigint N, float* X, bigint K, kmeans_opts opts)
{
    if (K > N)
        K = N;
    double* centroids = (double*)malloc(sizeof(double) * M * K);
    kmeans_initialize(centroids, M, N, K, X);
    for (bigint it = 1; it <= opts.num_iterations; it++) {
        kmeans_assign(labels, M, N, K, X, centroids);
        kmeans_centroids(centroids, M, N, K, X, labels);
    }
    kmeans_assign(labels, M, N, K, X, centroids);
    free(centroids);
}

void extract_subarray(float* X_sub, bigint M, float* X, const std::vector<bigint>& inds)
{
    for (bigint i = 0; i < (bigint)inds.size(); i++) {
        for (bigint m = 0; m < M; m++) {
            X_sub[m + i * M] = X[m + inds[i] * M];
        }
    }
}

void kmeans_maxsize(int* labels, bigint M, bigint N, float* X, bigint maxsize, kmeans_opts opts)
{
    if (N <= maxsize) {
        for (bigint i = 0; i < N; i++)
            labels[i] = 1;
        return;
    }
    bigint K = ceil(N * 1.0 / maxsize);
    int* labels1 = (int*)malloc(sizeof(int) * N);
    kmeans(labels1, M, N, X, K, opts);
    bigint L1 = compute_max(N, labels1);
    bigint current_max_k = 0;
    for (bigint k = 1; k <= L1; k++) {
        std::vector<bigint> inds_k;
        for (bigint i = 0; i < N; i++) {
            if (labels1[i] == k)
                inds_k.push_back(i);
        }
        if (inds_k.size() > 0) {
            float* X2 = (float*)malloc(sizeof(float) * M * inds_k.size());
            int* labels2 = (int*)malloc(sizeof(int) * inds_k.size());
            extract_subarray(X2, M, X, inds_k);
            kmeans_maxsize(labels2, M, inds_k.size(), X2, maxsize, opts);
            for (bigint j = 0; j < (bigint)inds_k.size(); j++) {
                labels[inds_k[j]] = current_max_k + labels2[j];
            }
            current_max_k += compute_max(inds_k.size(), labels2);
            free(X2);
            free(labels2);
        }
    }
    free(labels1);
}

void kmeans_multistep(int* labels, bigint M, bigint N, float* X, bigint K1, bigint K2, bigint K3, kmeans_opts opts)
{
    if (K2 > 1) {
        int* labels1 = (int*)malloc(sizeof(int) * N);
        for (bigint i = 0; i < N; i++)
            labels1[i] = 0;
        kmeans_multistep(labels1, M, N, X, K2, K3, 0, opts);
        bigint L1 = compute_max(N, labels1);
        bigint current_max_k = 0;
        for (bigint k = 1; k <= L1; k++) {
            std::vector<bigint> inds_k;
            for (bigint i = 0; i < N; i++) {
                if (labels1[i] == k)
                    inds_k.push_back(i);
            }
            if (inds_k.size() > 0) {
                float* X2 = (float*)malloc(sizeof(float) * M * inds_k.size());
                int* labels2 = (int*)malloc(sizeof(int) * inds_k.size());
                extract_subarray(X2, M, X, inds_k);
                kmeans_multistep(labels2, M, inds_k.size(), X2, K1, 0, 0, opts);
                for (bigint j = 0; j < (bigint)inds_k.size(); j++) {
                    labels[inds_k[j]] = current_max_k + labels2[j];
                }
                current_max_k += compute_max(inds_k.size(), labels2);
                free(X2);
                free(labels2);
            }
        }
        free(labels1);
    }
    else {
        kmeans(labels, M, N, X, K1, opts);
    }
}

double dot_product(bigint N, float* X, float* Y)
{
    double ret = 0;
    for (bigint i = 0; i < N; i++) {
        ret += X[i] * Y[i];
    }
    return ret;
}

void normalize_vector(bigint N, float* V)
{
    double norm = sqrt(dot_product(N, V, V));
    if (!norm)
        return;
    for (bigint i = 0; i < N; i++)
        V[i] /= norm;
}

void compare_clusters(double* dip_score, std::vector<bigint>* new_labels1, std::vector<bigint>* new_labels2, bigint M, bigint N1, bigint N2, float* X1, float* X2, double* centroid1, double* centroid2)
{
    float* V = (float*)malloc(sizeof(float) * M);
    float* projection = (float*)malloc(sizeof(float) * (N1 + N2));
    for (bigint m = 0; m < M; m++) {
        V[m] = centroid2[m] - centroid1[m];
    }
    normalize_vector(M, V);
    for (bigint i = 0; i < N1; i++) {
        projection[i] = dot_product(M, V, &X1[M * i]);
    }
    for (bigint i = 0; i < N2; i++) {
        projection[N1 + i] = dot_product(M, V, &X2[M * i]);
    }
    isocut5_opts icopts;
    icopts.already_sorted = false;
    double cutpoint;
    isocut5(dip_score, &cutpoint, N1 + N2, projection, icopts);
    new_labels1->resize(N1);
    new_labels2->resize(N2);
    for (bigint i = 0; i < N1; i++) {
        if (projection[i] < cutpoint)
            (*new_labels1)[i] = 1;
        else
            (*new_labels1)[i] = 2;
    }
    for (bigint i = 0; i < N2; i++) {
        if (projection[N1 + i] < cutpoint)
            (*new_labels2)[i] = 1;
        else
            (*new_labels2)[i] = 2;
    }
    free(projection);
    free(V);
}

void compute_centroids(float* centroids, bigint M, bigint N, bigint Kmax, float* X, int* labels, std::vector<bigint>& cluster_to_compute_vec)
{
    std::vector<double> C(M * Kmax);
    for (bigint jj = 0; jj < M * Kmax; jj++)
        C[jj] = 0;
    std::vector<double> counts(Kmax);
    for (bigint k = 0; k < Kmax; k++)
        counts[k] = 0;
    for (bigint i = 0; i < N; i++) {
        bigint k0 = labels[i];
        bigint i0 = k0 - 1;
        if (cluster_to_compute_vec[i0]) {
            for (bigint m = 0; m < M; m++) {
                C[m + M * i0] += X[m + M * i];
            }
            counts[i0]++;
        }
    }
    for (bigint k = 0; k < Kmax; k++) {
        if (cluster_to_compute_vec[k]) {
            if (counts[k]) {
                for (bigint m = 0; m < M; m++) {
                    C[m + M * k] /= counts[k];
                }
            }
        }
    }
    for (bigint k = 0; k < Kmax; k++) {
        if (cluster_to_compute_vec[k]) {
            for (bigint m = 0; m < M; m++) {
                centroids[m + k * M] = C[m + k * M];
            }
        }
    }
}

void compute_covmats(float* covmats, bigint M, bigint N, bigint Kmax, float* X, int* labels, float* centroids, std::vector<bigint>& clusters_to_compute_vec)
{
    std::vector<double> C(M * M * Kmax);
    for (bigint jj = 0; jj < M * M * Kmax; jj++)
        C[jj] = 0;
    std::vector<double> counts(Kmax);
    for (bigint k = 0; k < Kmax; k++)
        counts[k] = 0;
    for (bigint i = 0; i < N; i++) {
        bigint i0 = labels[i] - 1;
        if (clusters_to_compute_vec[i0]) {
            for (bigint m1 = 0; m1 < M; m1++) {
                for (bigint m2 = 0; m2 < M; m2++) {
                    C[m1 + M * m2 + M * M * i0] += (X[m1 + M * i] - centroids[m1 + i0 * M]) * (X[m2 + M * i] - centroids[m2 + i0 * M]);
                }
            }
            counts[i0]++;
        }
    }
    for (bigint k = 0; k < Kmax; k++) {
        if (clusters_to_compute_vec[k]) {
            if (counts[k]) {
                for (bigint m1 = 0; m1 < M; m1++) {
                    for (bigint m2 = 0; m2 < M; m2++) {
                        C[m1 + m2 * M + M * M * k] /= counts[k];
                    }
                }
            }
        }
    }
    for (bigint k = 0; k < Kmax; k++) {
        if (clusters_to_compute_vec[k]) {
            for (bigint mm = 0; mm < M * M; mm++) {
                covmats[mm + k * M * M] = C[mm + k * M * M];
            }
        }
    }
}

void get_pairs_to_compare(std::vector<bigint>* inds1, std::vector<bigint>* inds2, bigint M, bigint K, float* active_centroids, const intarray2d& active_comparisons_made)
{
    inds1->clear();
    inds2->clear();
    double dists[K][K];
    for (bigint k1 = 0; k1 < K; k1++) {
        for (bigint k2 = 0; k2 < K; k2++) {
            if ((active_comparisons_made[k1][k2]) || (k1 == k2))
                dists[k1][k2] = -1;
            else {
                double dist = 0;
                for (bigint m = 0; m < M; m++) {
                    double val = active_centroids[m + M * k1] - active_centroids[m + M * k2];
                    dist += val * val;
                }
                dist = sqrt(dist);
                dists[k1][k2] = dist;
            }
        }
    }
    // important to only take the mutal closest pairs -- unlike how we originally did it
    //bool something_changed = true;
    //while (something_changed) {
    //something_changed = false;
    std::vector<bigint> best_inds(K);
    for (bigint k = 0; k < K; k++) {
        bigint best_ind = -1;
        double best_distance = -1;
        for (bigint k2 = 0; k2 < K; k2++) {
            if (dists[k][k2] >= 0) {
                if ((best_distance < 0) || (dists[k][k2] < best_distance)) {
                    best_distance = dists[k][k2];
                    best_ind = k2;
                }
            }
        }
        best_inds[k] = best_ind;
    }
    for (bigint j = 0; j < K; j++) {
        if (best_inds[j] > j) {
            if (best_inds[best_inds[j]] == j) { //mutual
                if (dists[j][best_inds[j]] >= 0) {
                    inds1->push_back(j + 1);
                    inds2->push_back(best_inds[j] + 1);
                    for (bigint aa = 0; aa < K; aa++) {
                        dists[j][aa] = -1;
                        dists[aa][j] = -1;
                        dists[best_inds[j]][aa] = -1;
                        dists[aa][best_inds[j]] = -1;
                    }
                    //something_changed = true;
                }
            }
        }
    }
    //}
}

std::vector<float> compute_centroid(bigint M, bigint N, float* X)
{
    std::vector<double> ret(M);
    double count = 0;
    for (bigint m = 0; m < M; m++) {
        ret[m] = 0;
    }
    for (bigint i = 0; i < N; i++) {
        for (bigint m = 0; m < M; m++) {
            ret[m] += X[m + M * i];
        }
        count++;
    }
    if (count) {
        for (bigint m = 0; m < M; m++) {
            ret[m] /= count;
        }
    }
    std::vector<float> retf(M);
    for (bigint m = 0; m < M; m++)
        retf[m] = ret[m];
    return retf;
}

bool matinv(bigint M, float* out, float* in)
{
    return smi::get_inverse_via_lu_decomposition(M, out, in);
}

void matvec(bigint M, bigint N, float* out, float* mat, float* vec)
{
    for (bigint m = 0; m < M; m++) {
        float val = 0;
        for (bigint n = 0; n < N; n++) {
            val += mat[m + M * n] * vec[n];
        }
        out[m] = val;
    }
}

double dbg_compute_mean(const std::vector<float>& X)
{
    double ret = 0;
    for (bigint i = 0; i < (bigint)X.size(); i++)
        ret += X[i];
    return ret / X.size();
}

double dbg_compute_var(const std::vector<float>& X)
{
    double mu = dbg_compute_mean(X);
    double ret = 0;
    for (bigint i = 0; i < (bigint)X.size(); i++)
        ret += (X[i] - mu) * (X[i] - mu);
    return ret / X.size();
}

bool merge_test(std::vector<bigint>* L12, bigint M, bigint N1, bigint N2, float* X1, float* X2, const isosplit5_opts& opts, float* centroid1, float* centroid2, float* covmat1, float* covmat2)
{
    L12->resize(N1 + N2);
    for (bigint i = 0; i < N1 + N2; i++)
        (*L12)[i] = 1;
    if ((N1 == 0) || (N2 == 0)) {
        printf("Error in merge test: N1 or N2 is zero.\n");
        return true;
    }

    //std::vector<float> centroid1 = compute_centroid(M, N1, X1);
    //std::vector<float> centroid2 = compute_centroid(M, N2, X2);

    std::vector<float> V(M);
    for (bigint m = 0; m < M; m++) {
        V[m] = centroid2[m] - centroid1[m];
    }

    std::vector<float> avg_covmat;
    avg_covmat.resize(M * M);
    for (bigint rr = 0; rr < M * M; rr++) {
        avg_covmat[rr] = (covmat1[rr] + covmat2[rr]) / 2;
    }
    std::vector<float> inv_avg_covmat;
    inv_avg_covmat.resize(M * M);
    if (!matinv(M, inv_avg_covmat.data(), avg_covmat.data())) {
        printf("Unable to invert matrix. This may be due to the fact that you have duplicate events. Contact Jeremy if this is not the case, or if you would prefer the program to continue in this case. Aborting.\n");
        abort();
        return false;
    }

    std::vector<float> V2(M);
    matvec(M, M, V2.data(), inv_avg_covmat.data(), V.data());
    //matvec(M,M,V.data(),inv_avg_covmat.data(),V2.data());
    for (bigint i = 0; i < M; i++)
        V[i] = V2[i];

    /*
    printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
    print_matrix(M,M,avg_covmat.data());
    printf("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\n");
    print_matrix(M,M,inv_avg_covmat.data());
    printf("\n\n");
    */

    double sumsqr = 0;
    for (bigint m = 0; m < M; m++) {
        sumsqr += V[m] * V[m];
    }
    if (sumsqr) {
        for (bigint m = 0; m < M; m++)
            V[m] /= sqrt(sumsqr);
    }

    std::vector<float> projection1(N1), projection2(N2), projection12(N1 + N2);
    for (bigint i = 0; i < N1; i++) {
        double tmp = 0;
        for (bigint m = 0; m < M; m++)
            tmp += V[m] * X1[m + i * M];
        projection1[i] = tmp;
        projection12[i] = tmp;
    }
    for (bigint i = 0; i < N2; i++) {
        double tmp = 0;
        for (bigint m = 0; m < M; m++)
            tmp += V[m] * X2[m + i * M];
        projection2[i] = tmp;
        projection12[N1 + i] = tmp;
    }

    bool do_merge;
    isocut5_opts oo;
    oo.already_sorted = false;
    double dipscore, cutpoint;
    isocut5(&dipscore, &cutpoint, N1 + N2, projection12.data(), oo);

    if (dipscore < opts.isocut_threshold) {
        do_merge = true;
    }
    else {
        do_merge = false;
    }
    for (bigint i = 0; i < N1 + N2; i++) {
        if (projection12[i] < cutpoint)
            (*L12)[i] = 1;
        else
            (*L12)[i] = 2;
    }

    return do_merge;
}

void compare_pairs(std::vector<bigint>* clusters_changed, bigint* total_num_label_changes, bigint M, bigint N, float* X, int* labels, const std::vector<bigint>& k1s, const std::vector<bigint>& k2s, const isosplit5_opts& opts, float* centroids, float* covmats)
{
    bigint Kmax = ns_isosplit5::compute_max(N, labels);
    std::vector<bigint> clusters_changed_vec(Kmax);
    for (bigint i = 0; i < Kmax; i++)
        clusters_changed_vec[i] = 0;
    int* new_labels = (int*)malloc(sizeof(bigint) * N);
    *total_num_label_changes = 0;
    for (bigint i = 0; i < N; i++)
        new_labels[i] = labels[i];
    for (bigint i1 = 0; i1 < (bigint)k1s.size(); i1++) {
        int k1 = k1s[i1];
        int k2 = k2s[i1];
        std::vector<bigint> inds1, inds2;
        for (bigint i = 0; i < N; i++) {
            if (labels[i] == k1)
                inds1.push_back(i);
            if (labels[i] == k2)
                inds2.push_back(i);
        }
        if ((inds1.size() > 0) && (inds2.size() > 0)) {
            std::vector<bigint> inds12;
            inds12.insert(inds12.end(), inds1.begin(), inds1.end());
            inds12.insert(inds12.end(), inds2.begin(), inds2.end());
            std::vector<bigint> L12_old(inds12.size());
            for (bigint i = 0; i < (bigint)inds1.size(); i++)
                L12_old[i] = 1;
            for (bigint i = 0; i < (bigint)inds2.size(); i++)
                L12_old[inds1.size() + i] = 2;
            std::vector<bigint> L12(inds12.size());

            bool do_merge;
            if (((bigint)inds1.size() < opts.min_cluster_size) || ((bigint)inds2.size() < opts.min_cluster_size)) {
                do_merge = true;
            }
            else {
                float* X1 = (float*)malloc(sizeof(float) * M * inds1.size());
                float* X2 = (float*)malloc(sizeof(float) * M * inds2.size());
                extract_subarray(X1, M, X, inds1);
                extract_subarray(X2, M, X, inds2);
                do_merge = merge_test(&L12, M, inds1.size(), inds2.size(), X1, X2, opts, &centroids[(k1 - 1) * M], &centroids[(k2 - 1) * M], &covmats[(k1 - 1) * M * M], &covmats[(k2 - 1) * M * M]);
                free(X1);
                free(X2);
            }
            if (do_merge) {
                for (bigint i = 0; i < (bigint)inds2.size(); i++) {
                    new_labels[inds2[i]] = k1;
                }
                *total_num_label_changes += inds2.size();
                clusters_changed_vec[k1 - 1] = 1;
                clusters_changed_vec[k2 - 1] = 1;
            }
            else {
                //redistribute
                bool something_was_redistributed = false;
                for (bigint i = 0; i < (bigint)inds1.size(); i++) {
                    if (L12[i] == 2) {
                        new_labels[inds1[i]] = k2;
                        (*total_num_label_changes)++;
                        something_was_redistributed = true;
                    }
                }
                for (bigint i = 0; i < (bigint)inds2.size(); i++) {
                    if (L12[inds1.size() + i] == 1) {
                        new_labels[inds2[i]] = k1;
                        (*total_num_label_changes)++;
                        something_was_redistributed = true;
                    }
                }
                if (something_was_redistributed) {
                    clusters_changed_vec[k1 - 1] = 1;
                    clusters_changed_vec[k2 - 1] = 1;
                }
            }
        }
    }
    clusters_changed->clear();
    for (int k = 0; k < Kmax; k++)
        if (clusters_changed_vec[k])
            clusters_changed->push_back(k + 1);
    for (bigint i = 0; i < N; i++)
        labels[i] = new_labels[i];
    free(new_labels);
}
}

void get_pairs_to_compare3(std::vector<bigint>* i1s, std::vector<bigint>* i2s, bigint M, bigint N, double* centroids)
{
    float distances[N][N];
    bigint used[N];
    for (bigint i = 0; i < N; i++)
        used[i] = 0;
    for (bigint i = 0; i < N; i++) {
        for (bigint j = i; j < N; j++) {
            double sumsqr = 0;
            for (bigint m = 0; m < M; m++) {
                float diff0 = centroids[m + M * i] - centroids[m + M * j];
                sumsqr += diff0 * diff0;
            }
            distances[i][j] = distances[j][i] = sqrt(sumsqr);
        }
    }

    bool something_changed = true;
    while (something_changed) {
        something_changed = false;
        bigint closest[N];
        for (bigint i = 0; i < N; i++)
            closest[i] = -1;
        for (bigint i = 0; i < N; i++) {
            double best_distance = -1;
            if (!used[i]) {
                for (bigint j = 0; j < N; j++) {
                    if ((!used[j]) && (j != i)) {
                        double dist = distances[i][j];
                        if ((best_distance < 0) || (dist < best_distance)) {
                            best_distance = dist;
                            closest[i] = j;
                        }
                    }
                }
            }
        }
        for (bigint i = 0; i < N; i++) {
            if (!used[i]) {
                if ((closest[i] >= 0) && (closest[i] > i) && (closest[closest[i]] == i)) { //mutual nearest neighbor among those that have not yet been used
                    i1s->push_back(i);
                    i2s->push_back(closest[i]);
                    used[i] = 1;
                    used[closest[i]] = 1;
                    something_changed = true;
                }
            }
        }
    }
}

void get_pairs_to_compare2(std::vector<bigint>* i1s, std::vector<bigint>* i2s, bigint M, bigint N, double* centroids)
{
    float* centroidsf = (float*)malloc(sizeof(float) * M * N);
    for (bigint i = 0; i < M * N; i++)
        centroidsf[i] = centroids[i];
    int* groups = (int*)malloc(sizeof(int) * N);
    bigint maxsize = 1000;
    ns_isosplit5::kmeans_opts oo;
    ns_isosplit5::kmeans_maxsize(groups, M, N, centroidsf, maxsize, oo);
    bigint num_groups = ns_isosplit5::compute_max(N, groups);

    for (bigint group_number = 1; group_number <= num_groups; group_number++) {
        std::vector<bigint> inds_group;
        for (bigint i = 0; i < N; i++)
            if (groups[i] == group_number)
                inds_group.push_back(i);
        bigint N0 = inds_group.size();
        if (N0 > 0) {
            float* centroids0f = (float*)malloc(sizeof(float) * M * N0);
            double* centroids0 = (double*)malloc(sizeof(double) * M * N0);
            ns_isosplit5::extract_subarray(centroids0f, M, centroidsf, inds_group);
            for (bigint i = 0; i < M * N0; i++)
                centroids0[i] = centroids0f[i];
            std::vector<bigint> i1s0, i2s0;
            get_pairs_to_compare3(&i1s0, &i2s0, M, N0, centroids0);
            for (bigint jj = 0; jj < (bigint)i1s0.size(); jj++) {
                i1s->push_back(inds_group[i1s0[jj]]);
                i2s->push_back(inds_group[i2s0[jj]]);
            }
            free(centroids0f);
            free(centroids0);
        }
    }

    free(groups);
    free(centroidsf);
}

////////////// SQUARE MATRIX INVERSION ////////////////////

namespace smi {

// calculate the cofactor of element (row,col)
void get_minor(bigint M, float* out, float* in, bigint row, bigint col)
{
    bigint rr = 0;
    for (bigint i = 0; i < M; i++) {
        if (i != row) {
            bigint cc = 0;
            for (bigint j = 0; j < M; j++) {
                if (j != col) {
                    out[rr + (M - 1) * cc] = in[i + M * j];
                    cc++;
                }
            }
            rr++;
        }
    }
}

// Calculate the determinant recursively.
double determinant(bigint M, float* A)
{
    if (M < 1) {
        printf("Error. Cannot take determinant when M<1.\n");
        return 0;
    }
    if (M == 1)
        return A[0];

    double ret = 0;

    std::vector<float> minor;
    minor.resize((M - 1) * (M - 1));
    for (bigint i = 0; i < M; i++) {
        get_minor(M, minor.data(), A, 0, i);
        double sgn = 1;
        if (i % 2 == 1)
            sgn = -1;
        ret += A[0 + M * i] * determinant(M - 1, minor.data()) * sgn;
    }

    return ret;
}

// matrix inversion
void get_inverse_via_formula(bigint M, float* out, float* in)
{
    if (M == 1) {
        if (in[0])
            out[0] = 1 / in[0];
        else
            out[0] = 0;
        return;
    }
    // get the determinant of a
    double det = determinant(M, in);
    if (det == 0) {
        for (bigint m = 0; m < M * M; m++)
            out[m] = 0;
        return;
    }
    det = 1 / det;

    std::vector<float> minor;
    minor.resize((M - 1) * (M - 1));
    for (bigint j = 0; j < M; j++) {
        for (bigint i = 0; i < M; i++) {
            // get the co-factor (matrix) of A(j,i)
            get_minor(M, minor.data(), in, j, i);
            out[i + M * j] = det * determinant(M - 1, minor.data());
            if ((i + j) % 2 == 1)
                out[i + M * j] = -out[i + M * j];
        }
    }
}

/* Copyright 2015 Chandra Shekhar (chandraiitk AT yahoo DOT co DOT in).
  Homepage: https://sites.google.com/site/chandraacads
 * * */

/* This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
 * * */

/* This program computes inverse of a square matrix, based on LUP decomposition.
 *
 * Tested with GCC-4.8.3 on 64 bit Linux (Fedora-20).
 *
 * Compilation:     "gcc -O2 Matrix_inverse_LUP.c -o Mat_inv_LUP.exe -lm -Wall"
 * Execution:     "./Mat_inv_LUP.exe"
 * * */

/* This function decomposes the matrix 'A' into L, U, and P. If successful,
 * the L and the U are stored in 'A', and information about the pivot in 'P'.
 * The diagonal elements of 'L' are all 1, and therefore they are not stored. */
bigint LUPdecompose(int M, float* A, int* P)
{
    bigint i, j, k, kd = 0, T;
    float p, t;

    /* Finding the pivot of the LUP decomposition. */
    for (i = 1; i < M; i++)
        P[i] = i; //Initializing.

    for (k = 1; k < M - 1; k++) {
        p = 0;
        for (i = k; i < M; i++) {
            t = A[i + M * k];
            if (t < 0)
                t *= -1; //Abosolute value of 't'.
            if (t > p) {
                p = t;
                kd = i;
            }
        }

        if (p == 0) {
            printf("\nLUPdecompose(): ERROR: A singular matrix is supplied.\n"
                   "\tRefusing to proceed any further.\n");
            return -1;
        }

        /* Exchanging the rows according to the pivot determined above. */
        T = P[kd];
        P[kd] = P[k];
        P[k] = T;
        for (i = 1; i < M; i++) {
            t = A[kd + M * i];
            A[kd + M * i] = A[k + M * i];
            A[k + M * i] = t;
        }

        for (i = k + 1; i < M; i++) //Performing substraction to decompose A as LU.
        {
            A[i + M * k] = A[i + M * k] / A[k + M * k];
            for (j = k + 1; j < M; j++)
                A[i + M * j] -= A[i + M * k] * A[k + M * j];
        }
    } //Now, 'A' contains the L (without the diagonal elements, which are all 1)
    //and the U.

    return 0;
}

/* This function calculates the inverse of the LUP decomposed matrix 'LU' and pivoting
 * information stored in 'P'. The inverse is returned through the matrix 'LU' itselt.
 * 'B', X', and 'Y' are used as temporary spaces. */
bigint LUPinverse(bigint M, int* P, float* LU,
    float* B, float* X, float* Y)
{
    bigint i, j, n, m;
    float t;

    //Initializing X and Y.
    for (n = 1; n < M; n++)
        X[n] = Y[n] = 0;

    /* Solving LUX = Pe, in order to calculate the inverse of 'A'. Here, 'e' is a column
 * vector of the identity matrix of size 'M-1'. Solving for all 'e'. */
    for (i = 1; i < M; i++) {
        //Storing elements of the i-th column of the identity matrix in i-th row of 'B'.
        for (j = 1; j < M; j++)
            B[i + M * j] = 0;
        B[i + M * i] = 1;

        //Solving Ly = Pb.
        for (n = 1; n < M; n++) {
            t = 0;
            for (m = 1; m <= n - 1; m++)
                t += LU[n + M * m] * Y[m];
            Y[n] = B[i + M * P[n]] - t;
        }

        //Solving Ux = y.
        for (n = M - 1; n >= 1; n--) {
            t = 0;
            for (m = n + 1; m < M; m++)
                t += LU[n + M * m] * X[m];
            X[n] = (Y[n] - t) / LU[n + M * n];
        } //Now, X contains the solution.

        for (j = 1; j < M; j++)
            B[i + M * j] = X[j]; //Copying 'X' into the same row of 'B'.
    } //Now, 'B' the transpose of the inverse of 'A'.

    /* Copying transpose of 'B' into 'LU', which would the inverse of 'A'. */
    for (i = 1; i < M; i++)
        for (j = 1; j < M; j++)
            LU[i + M * j] = B[j + M * i];

    return 0;
}

bool get_inverse_via_lu_decomposition(int M, float* out, float* in)
{
    std::vector<float> A((M + 1) * (M + 1));
    std::vector<int> P((M + 1));
    for (bigint i = 0; i < M; i++) {
        for (bigint j = 0; j < M; j++) {
            A[(1 + i) + (M + 1) * (1 + j)] = in[i + M * j];
        }
    }
    int ret = LUPdecompose(M + 1, A.data(), P.data());
    if (ret < 0) {
        //handle case reported by Alex Morley. Don't proceed because inverse will crash (I believe)
        return false;
    }
    std::vector<float> B((M + 1) * (M + 1));
    std::vector<float> X(M + 1);
    std::vector<float> Y(M + 1);
    LUPinverse(M + 1, P.data(), A.data(), B.data(), X.data(), Y.data());
    for (bigint i = 0; i < M; i++) {
        for (bigint j = 0; j < M; j++) {
            out[i + M * j] = A[(1 + i) + (M + 1) * (1 + j)];
        }
    }
    return true;
}
}
