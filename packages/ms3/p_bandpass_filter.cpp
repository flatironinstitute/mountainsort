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
#include "p_bandpass_filter.h"

#include <QTime>
#include <diskreadmda32.h>
#include <diskwritemda.h>
#include "omp.h"
#include "fftw3.h"
#include <QFile>
#include <QFileInfo>
#include <QCoreApplication>
#include <cmath>
using std::fabs;
using std::erf;
using std::sqrt;

namespace P_bandpass_filter {
void define_kernel(bigint N, double* kernel, double samplefreq, double freq_min, double freq_max, double freq_wid);
void multiply_by_factor(bigint N, float* X, double factor);
struct Kernel_runner {
    Kernel_runner()
    {
    }

    ~Kernel_runner()
    {
        fftwf_free(data_in);
        fftwf_free(data_out);
        free(kernel0);
        //delete p_fft;
        //delete p_ifft;
    }
    void init(bigint M_in, bigint N_in, double samplerate, double freq_min, double freq_max, double freq_wid)
    {
        M = M_in;
        N = N_in;
        MN = M * N;

        /*
        p_fft=new fftwf_plan; //this nonsense is necessary because we cannot instantiate fftw plans in multiple threads simultaneously
        p_ifft=new fftwf_plan;
        */

        data_in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * MN);
        data_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * MN);
        kernel0 = (double*)malloc(sizeof(double) * N);

        define_kernel(N, kernel0, samplerate, freq_min, freq_max, freq_wid);

        bigint rank = 1;
        int n[] = { (int)N };
        bigint howmany = M;
        int* inembed = n;
        bigint istride = M;
        bigint idist = 1;
        int* onembed = n;
        bigint ostride = M;
        bigint odist = 1;
        unsigned flags = FFTW_ESTIMATE;
        p_fft = fftwf_plan_many_dft(rank, n, howmany, data_in, inembed, istride, idist, data_out, onembed, ostride, odist, FFTW_FORWARD, flags);
        p_ifft = fftwf_plan_many_dft(rank, n, howmany, data_out, inembed, istride, idist, data_in, onembed, ostride, odist, FFTW_BACKWARD, flags);
    }
    void apply(Mda32& chunk)
    {
        //set input data
        for (bigint i = 0; i < MN; i++) {
            data_in[i][0] = chunk.get(i);
            data_in[i][1] = 0;
        }
        //fft
        fftwf_execute(p_fft);
        //multiply by kernel
        double factor = 1.0 / N;
        bigint aa = 0;
        for (bigint i = 0; i < N; i++) {
            for (bigint m = 0; m < M; m++) {
                data_out[aa][0] *= kernel0[i] * factor;
                data_out[aa][1] *= kernel0[i] * factor;
                aa++;
            }
        }
        fftwf_execute(p_ifft);
        //set the output data
        for (bigint i = 0; i < MN; i++) {
            chunk.set(data_in[i][0], i);
        }
    }

    bigint M;
    bigint N, MN;
    fftwf_complex* data_in;
    fftwf_complex* data_out;
    double* kernel0;
    fftwf_plan p_fft;
    fftwf_plan p_ifft;
};
Mda32 bandpass_filter_kernel(Mda32& X, double samplerate, double freq_min, double freq_max, double freq_wid);
Mda32 subsample(const Mda32& timeseries, int subsample_factor);
}

bool p_bandpass_filter(QString timeseries, QString timeseries_out, Bandpass_filter_opts &opts)
{
    if (opts.freq_max == 0) {
        return QFile::copy(timeseries, timeseries_out);
    }

    if (!opts.subsample_factor)
        opts.subsample_factor = 1;

    bool do_write = true;
    //if (opts.testcode.split(",").contains("nowrite"))
    //    do_write = false;

    DiskReadMda32 X;
    if (QFileInfo(timeseries).isDir())
        X.setConcatDirectory(2, timeseries);
    else
        X.setPath(timeseries);

    const bigint M = X.N1();
    const bigint Naa = X.N2();
    bigint N2 = Naa / opts.subsample_factor;

    bigint dtype = MDAIO_TYPE_FLOAT32;
    if (opts.quantization_unit) {
        dtype = MDAIO_TYPE_INT16;
    }
    DiskWriteMda Ybb(dtype, timeseries_out, M, N2);

    QTime timer_status;
    timer_status.start();

    bigint num_threads = omp_get_max_threads();

    int overhead = 6; // complex in, complex out, chunk, chunk2, extra
    bigint overlap_size = 60000;
    bigint min_chunk_size = overlap_size*2;
    double target_ram_bytes = 1.0 * (1024*1024*1024);
    double target_chunk_size = target_ram_bytes / (overhead * sizeof(float) * M * num_threads) - 2 * overlap_size;
    bigint chunk_size = (bigint)(target_chunk_size);
    if (chunk_size < min_chunk_size) chunk_size = min_chunk_size;

    double expected_ram_bytes = overhead * sizeof(float)*M*num_threads*(chunk_size+2*overlap_size);
    opts.expected_peak_ram_mb = expected_ram_bytes / (1024 * 1024);

    if (opts.requirements_only) {
        return true;
    }

    printf("************+++ Using chunk size / overlap size: %ld / %ld (num threads=%ld)\n", chunk_size, overlap_size, num_threads);
    qDebug().noquote() << "Expected peak RAM usage (MB):" << opts.expected_peak_ram_mb;
    qDebug().noquote() << "samplerate/freq_min/freq_max/freq_wid:" << opts.samplerate << opts.freq_min << opts.freq_max << opts.freq_wid;

    bool ret = true;
    bigint num_timepoints_handled = 0;
#pragma omp parallel
    {
        // one kernel runner for each parallel thread so they don't intersect
        P_bandpass_filter::Kernel_runner KR;
#pragma omp critical(lock1)
        {
            KR.init(M, chunk_size + 2 * overlap_size, opts.samplerate, opts.freq_min, opts.freq_max, opts.freq_wid);
        }
#pragma omp for
        for (bigint timepoint = 0; timepoint < Naa; timepoint += chunk_size) {
            Mda32 chunk;
#pragma omp critical(lock1)
            {
                if (!X.readChunk(chunk, 0, timepoint - overlap_size, M, chunk_size + 2 * overlap_size)) {
                    qWarning() << "Error reading chunk";
                    ret = false;
                }
            }
            //if (!opts.testcode.split(",").contains("nokernel")) {
            QTime kernel_timer;
            kernel_timer.start();
            KR.apply(chunk);
            //chunk = P_bandpass_filter::bandpass_filter_kernel(chunk, opts.samplerate, opts.freq_min, opts.freq_max, opts.freq_wid);
            //qDebug().noquote() << "Kernel timer elapsed: " << kernel_timer.elapsed() << " for chunk at " << timepoint << " of " << N;
            //}

            Mda32 chunk2;
            {
                chunk.getChunk(chunk2, 0, overlap_size, M, chunk_size);
            }
#pragma omp critical(lock1)
            {
                {
                    if (do_write) {
                        if (opts.subsample_factor > 1) {
                            chunk2 = P_bandpass_filter::subsample(chunk2, opts.subsample_factor);
                        }
                        if (opts.quantization_unit) {
                            P_bandpass_filter::multiply_by_factor(chunk2.totalSize(), chunk2.dataPtr(), 1.0 / opts.quantization_unit);
                        }
                        if (!Ybb.writeChunk(chunk2, 0, timepoint / opts.subsample_factor)) {
                            qWarning() << "Error writing chunk";
                            ret = false;
                        }
                    }
                }
                num_timepoints_handled += qMin((bigint)chunk_size, Naa - timepoint);
                if ((timer_status.elapsed() > 5000) || (num_timepoints_handled == Naa) || (timepoint == 0)) {
                    printf("%ld/%ld (%d%%) -- using %d threads.\n",
                        num_timepoints_handled, Naa,
                        (int)(num_timepoints_handled * 1.0 / Naa * 100),
                        omp_get_num_threads());
                    timer_status.restart();
                }
            }
        }
    }

    return ret;
}

namespace P_bandpass_filter {

void multiply_by_factor(bigint N, float* X, double factor)
{
    /*bigint start = 0;
#ifdef USE_SSE2
    __m128d factor_m128 = _mm_load_pd1(&factor);
    for (; start < (N / 2) * 2; start += 2) {
        double* chunk = X + start;
        __m128d x = _mm_load_pd(chunk);
        __m128d result = _mm_mul_pd(x, factor_m128);
        _mm_store_pd(chunk, result);
    }
#endif
*/
    for (bigint i = 0; i < N; i++)
        X[i] *= factor;
}

bool do_fft_1d_r2c(bigint M, bigint N, float* out, float* in)
{
    /*
    if (num_threads>1) {
        fftwf_init_threads();
        fftwf_plan_with_nthreads(num_threads);
    }
    */

    bigint MN = M * N;

    fftwf_complex* in2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * MN);
    fftwf_complex* out2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * MN);
    for (bigint ii = 0; ii < MN; ii++) {
        //in2[ii][0]=in[ii*2];
        //in2[ii][1]=in[ii*2+1];
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
     * transforms. The basic fftwf_plan_dft interface corresponds
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
     * (The basic fftwf_plan_dft interface corresponds to a stride
     * of 1.)
     */
    fftwf_plan p;
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
    p = fftwf_plan_many_dft(rank, n, howmany, in2, inembed, istride, idist, out2, onembed, ostride, odist, sign, flags);
    //p=fftwf_plan_dft_1d(N,in2,out2,fftwf_FORWARD,fftwf_ESTIMATE);

    fftwf_execute(p);
    for (bigint ii = 0; ii < MN; ii++) {
        out[ii * 2] = out2[ii][0];
        out[ii * 2 + 1] = out2[ii][1];
    }
    fftwf_free(in2);
    fftwf_free(out2);

/*
    if (num_threads>1) {
        fftwf_cleanup_threads();
    }
    */

#pragma omp critical
    fftwf_destroy_plan(p);

    return true;
}

bool do_ifft_1d_c2r(bigint M, bigint N, float* out, float* in)
{
    /*
    if (num_threads>1) {
        fftwf_init_threads();
        fftwf_plan_with_nthreads(num_threads);
    }
    */

    bigint MN = M * N;

    fftwf_complex* in2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * MN);
    fftwf_complex* out2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * MN);
    for (bigint ii = 0; ii < MN; ii++) {
        in2[ii][0] = in[ii * 2];
        in2[ii][1] = in[ii * 2 + 1];
    }

    fftwf_plan p;
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
    p = fftwf_plan_many_dft(rank, n, howmany, in2, inembed, istride, idist, out2, onembed, ostride, odist, sign, flags);
    //p=fftwf_plan_dft_1d(N,in2,out2,fftwf_BACKWARD,fftwf_ESTIMATE);

    fftwf_execute(p);
    for (bigint ii = 0; ii < MN; ii++) {
        out[ii] = out2[ii][0];
    }
    fftwf_free(in2);
    fftwf_free(out2);

/*
    if (num_threads>1) {
        fftwf_cleanup_threads();
    }
    */

#pragma omp critical
    fftwf_destroy_plan(p);

    return true;
}

void multiply_complex_by_real_kernel(bigint M, bigint N, float* Y, double* kernel)
{
    bigint bb = 0;
    bigint aa = 0;
    /*
#ifdef USE_SSE2
    for (bigint i = 0; i < N; i++) {
        __m128d kernel_m128 = _mm_load_pd1(kernel + aa);
        for (bigint j = 0; j < M; j++) {
            __m128d Y_m128 = _mm_load_pd(Y + bb * 2);
            __m128d result = _mm_mul_pd(Y_m128, kernel_m128);
            _mm_store_pd(Y + bb * 2, result);
            bb++;
        }
        aa++;
    }
#else
*/
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

Mda32 bandpass_filter_kernel(Mda32& X, double samplerate, double freq_min, double freq_max, double freq_wid)
{
    QTime timer;
    timer.start();
    bigint M = X.N1();
    bigint N = X.N2();
    bigint MN = M * N;
    Mda32 Y(M, N);
    float* Xptr = X.dataPtr();
    float* Yptr = Y.dataPtr();

    double* kernel0 = (double*)allocate(sizeof(double) * N);
    float* Xhat = (float*)allocate(sizeof(float) * MN * 2);
    define_kernel(N, kernel0, samplerate, freq_min, freq_max, freq_wid);

    do_fft_1d_r2c(M, N, Xhat, Xptr);
    multiply_complex_by_real_kernel(M, N, Xhat, kernel0);
    do_ifft_1d_c2r(M, N, Yptr, Xhat);

    multiply_by_factor(MN, Yptr, 1.0 / N);

    free(kernel0);
    free(Xhat);

    double sec = timer.elapsed() * 1.0 / 1000;
    double rate = X.totalSize() / sec;
    printf("bandpass filter rate: %g numbers/sec\n", rate);

    return Y;
}

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
}
