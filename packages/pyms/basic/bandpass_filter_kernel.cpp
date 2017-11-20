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
#include "omp.h"
#include "fftw3.h"
#include "bandpass_filter_kernel.h"
#include <math.h>
#include <stdlib.h>

void multiply_by_factor(bigint N, float* X, double factor)
{
    for (bigint i = 0; i < N; i++)
        X[i] *= factor;
}

bool do_fft_1d_r2c(bigint M, bigint N, float* out, float* in)
{
    bigint MN = M * N;

    fftw_complex* in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * MN);
    fftw_complex* out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * MN);
    for (bigint ii = 0; ii < MN; ii++) {
        in2[ii][0] = in[ii];
        in2[ii][1] = 0;
    }

    /*
     * From FFTW docs:
     * howmany is the number of transforms to compute.
     * The resulting plan computes howmany transforms,
     * where the input of the k-th transform is at
     * location in+k*idist (in C pointer arithmetic),
     * and its output is at location out+k*odist.
     * Plans obtained in this way can often be faster
     * than calling FFTW multiple times for the individual
     * transforms. The basic fftw_plan_dft interface corresponds
     * to howmany=1 (in which case the dist parameters are ignored).
     *
     * Each of the howmany transforms has rank rank
     * and size n, as in the basic interface.
     * In addition, the advanced interface allows the
     * input and output arrays of each transform to be
     * row-major subarrays of larger rank-rank arrays,
     * described by inembed and onembed parameters,
     * respectively. {i,o}nembed must be arrays of length
     * rank, and n should be elementwise less than or equal
     * to {i,o}nembed. Passing NULL for an nembed parameter
     * is equivalent to passing n (i.e. same physical and
     * logical dimensions, as in the basic interface.)
     *
     * The stride parameters indicate that the j-th element
     * of the input or output arrays is located at j*istride
     * or j*ostride, respectively. (For a multi-dimensional array,
     * j is the ordinary row-major index.) When combined with
     * the k-th transform in a howmany loop, from above, this
     * means that the (j,k)-th element is at j*stride+k*dist.
     * (The basic fftw_plan_dft interface corresponds to a stride
     * of 1.)
     */
    fftw_plan p;
    bigint rank = 1;
    int n[] = { (int)N };
    bigint howmany = M;
    int* inembed = n;
    bigint istride = M;
    bigint idist = 1;
    int* onembed = n;
    bigint ostride = M;
    bigint odist = 1;
    bigint sign = FFTW_FORWARD;
    unsigned flags = FFTW_ESTIMATE;
#pragma omp critical
    p = fftw_plan_many_dft(rank, n, howmany, in2, inembed, istride, idist, out2, onembed, ostride, odist, sign, flags);
    //p=fftw_plan_dft_1d(N,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE);

    fftw_execute(p);
    for (bigint ii = 0; ii < MN; ii++) {
        out[ii * 2] = out2[ii][0];
        out[ii * 2 + 1] = out2[ii][1];
    }
    fftw_free(in2);
    fftw_free(out2);

/*
    if (num_threads>1) {
        fftw_cleanup_threads();
    }
    */

#pragma omp critical
    fftw_destroy_plan(p);

    return true;
}

bool do_ifft_1d_c2r(bigint M, bigint N, float* out, float* in, bigint start_write, bigint end_write)
{
    /*
    if (num_threads>1) {
        fftw_init_threads();
        fftw_plan_with_nthreads(num_threads);
    }
    */

    bigint MN = M * N;

    fftw_complex* in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * MN);
    fftw_complex* out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * MN);
    for (bigint ii = 0; ii < MN; ii++) {
        in2[ii][0] = in[ii * 2];
        in2[ii][1] = in[ii * 2 + 1];
    }

    fftw_plan p;
    bigint rank = 1;
    int n[] = { (int)N };
    bigint howmany = M;
    int* inembed = n;
    bigint istride = M;
    bigint idist = 1;
    int* onembed = n;
    bigint ostride = M;
    bigint odist = 1;
    bigint sign = FFTW_BACKWARD;
    unsigned flags = FFTW_ESTIMATE;
#pragma omp critical
    p = fftw_plan_many_dft(rank, n, howmany, in2, inembed, istride, idist, out2, onembed, ostride, odist, sign, flags);
    //p=fftw_plan_dft_1d(N,in2,out2,FFTW_BACKWARD,FFTW_ESTIMATE);

    fftw_execute(p);

    for (bigint ii = start_write*M; ii < (end_write+1)*M; ii++) {
        out[ii] = out2[ii][0];
    }
    fftw_free(in2);
    fftw_free(out2);

/*
    if (num_threads>1) {
        fftw_cleanup_threads();
    }
    */

#pragma omp critical
    fftw_destroy_plan(p);

    return true;
}

void multiply_complex_by_real_kernel(bigint M, bigint N, float* Y, double* kernel)
{
    bigint bb = 0;
    bigint aa = 0;
    for (bigint i = 0; i < N; i++) {
        for (bigint j = 0; j < M; j++) {
            Y[bb * 2] *= kernel[aa];
            Y[bb * 2 + 1] *= kernel[aa];
            bb++;
        }
        aa++;
    }
    //#endif
}

void define_kernel(bigint N, double* kernel, double samplefreq, double freq_min, double freq_max, double freq_wid)
{
    // Matches ahb's code /matlab/processors/ms_bandpass_filter.m
    // improved ahb, changing tanh to erf, correct -3dB pts  6/14/16
    double T = N / samplefreq; // total time
    double df = 1 / T; // frequency grid
    double relwid = 3.0; // relative bottom-end roll-off width param, kills low freqs by factor 1e-5.

    //printf("filter params: %.15g %.15g %.15g \n", freq_min, freq_max, freq_wid); // debug
    //freq_wid = 1000.0; // *** why not correctly read in? override hack

    for (bigint i = 0; i < N; i++) {
        const double fgrid = (i <= (N + 1) / 2) ? df * i : df * (i - N); // why const? (ahb)
        const double absf = fabs(fgrid);
        double val = 1.0;
        if (freq_min != 0) { // (suggested by ahb) added on 3/3/16 by jfm
            if (i == 0)
                val = 0.0; // kill DC part exactly - ahb
            else
                val *= (1 + erf(relwid * (absf - freq_min) / freq_min)) / 2;
        }
        if (freq_max != 0) { // added on 3/3/16 by jfm
            val *= (1 - erf((absf - freq_max) / freq_wid)) / 2;
        }
        kernel[i] = sqrt(val); // note sqrt of filter func to apply to spectral intensity not ampl
    }
}

void bandpass_filter_kernel(bigint M,bigint N,float *X, double samplerate, double freq_min, double freq_max, double freq_wid, bigint start_write, bigint end_write)
{
    /*
    //subtract channel means
    double means[M];
    for (int m=0; m<M; m++) means[m]=0;
    for (bigint t=0; t<N; t++) {
        for (int m=0; m<M; m++) {
            means[m]+=X[m+M*t];
        }
    }
    for (int m=0; m<M; m++) means[m]/=N;
    for (bigint t=0; t<N; t++) {
        for (int m=0; m<M; m++) {
            X[m+M*t]-=means[m];
        }
    }
    */

    bigint MN = M * N;
    
    double* kernel0 = (double*)malloc(sizeof(double) * N);
    float* Xhat = (float*)malloc(sizeof(float) * MN * 2);
    define_kernel(N, kernel0, samplerate, freq_min, freq_max, freq_wid);

    do_fft_1d_r2c(M, N, Xhat, X);
    multiply_complex_by_real_kernel(M, N, Xhat, kernel0);
    do_ifft_1d_c2r(M, N, X, Xhat, start_write, end_write);

    //multiply_by_factor(MN, X, 1.0 / N);
    multiply_by_factor(M*(end_write-start_write+1), &X[M*start_write], 1.0 / N);

    free(kernel0);
    free(Xhat);
}

/*
Mda32 subsample(const Mda32& timeseries, int subsample_factor)
{
    int M = timeseries.N1();
    bigint N = timeseries.N2();
    Mda32 Y(M, N / subsample_factor);
    for (bigint n2 = 0; n2 < N / subsample_factor; n2++) {
        for (int m = 0; m < M; m++) {
            Y.set(timeseries.get(m, n2 * subsample_factor), m, n2);
        }
    }
    return Y;
}
*/

void bandpass_filter_kernel_multithread(bigint M,bigint N,float *X, double samplerate, double freq_min, double freq_max, double freq_wid) {
    bigint chunk_size =   200000;
    bigint overlap_size = 20000;

    //subtract channel means
    double means[M];
    for (int m=0; m<M; m++) means[m]=0;
    for (bigint t=0; t<N; t++) {
        for (int m=0; m<M; m++) {
            means[m]+=X[m+M*t];
        }
    }
    for (int m=0; m<M; m++) means[m]/=N;
    for (bigint t=0; t<N; t++) {
        for (int m=0; m<M; m++) {
            X[m+M*t]-=means[m];
        }
    }

    #pragma omp parallel for
    for (bigint t=0; t<N; t+=chunk_size) {
        if (t==0) printf("Using %d threads ---\n",omp_get_num_threads());
        bigint t1=t,t2=t+chunk_size-1;
        if (t2>N-1) t2=N-1;
        bigint s1=t1-overlap_size,s2=t2+overlap_size;
        if (s1<0) s1=0;
        if (s2>N-1) s2=N-1;
        bandpass_filter_kernel(M,s2-s1+1,&X[M*s1],samplerate,freq_min,freq_max,freq_wid,t1-s1,t2-s1);
    }
}