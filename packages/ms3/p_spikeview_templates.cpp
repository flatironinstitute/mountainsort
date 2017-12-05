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
#include "p_spikeview_templates.h"

#include <QTime>
#include <diskreadmda32.h>
#include <mda.h>
#include <cmath>
using std::fabs;
using std::erf;
using std::sqrt;
#include "fftw3.h"


namespace P_spikeview_templates {
struct ClusterData {
    QVector<double> times;
    Mda32 template0;
};

struct Bandpass_filter_runner {
    Bandpass_filter_runner();
    ~Bandpass_filter_runner();
    void init(bigint M_in, bigint N_in, double samplerate, double freq_min, double freq_max, double freq_wid);
    void apply(Mda32& chunk) const;
    void define_kernel(bigint N, double* kernel, double samplefreq, double freq_min, double freq_max, double freq_wid);
    bigint M;
    bigint N, MN;
    fftwf_complex* data_in = 0;
    fftwf_complex* data_out = 0;
    double* kernel0 = 0;
    fftwf_plan p_fft;
    fftwf_plan p_ifft;
};

QVector<double> subsample(const QVector<double>& X, bigint max_elments);
void compute_template(Mda32& template0, const DiskReadMda32& X, const QVector<double>& times, P_spikeview_templates_opts opts, const Bandpass_filter_runner& BP);
}

using namespace P_spikeview_templates;

bool p_spikeview_templates(QString timeseries, QString firings, QString templates_out, P_spikeview_templates_opts opts)
{
    Mda FF(firings);
    bigint L = FF.N2();
    QVector<double> times(L);
    QVector<int> labels(L);

    qDebug().noquote() << "Reading firings file...";
    for (bigint i = 0; i < L; i++) {
        times[i] = FF.value(1, i);
        labels[i] = FF.value(2, i);
    }

    QMap<int, ClusterData> cluster_data;
    for (bigint i = 0; i < L; i++) {
        int k = labels[i];
        cluster_data[k].times << times[i];
    }

    qDebug().noquote() << "Subsampling as needed" << opts.max_events_per_template;
    QList<int> keys = cluster_data.keys();
    foreach (int key, keys) {
        if ((opts.max_events_per_template > 0) && (cluster_data[key].times.count() > opts.max_events_per_template)) {
            cluster_data[key].times = subsample(cluster_data[key].times, opts.max_events_per_template);
        }
    }

    DiskReadMda32 X(timeseries);
    int M = X.N1();
    int T = opts.clip_size;
#pragma omp parallel
    {
        QTime timer;
        timer.start();
        bool first = true;
        Bandpass_filter_runner BP; //one per parallel thread
#pragma omp critical(cc1)
        {
            if ((opts.samplerate != 0) && ((opts.freq_min != 0) || (opts.freq_max != 0))) {
                BP.init(X.N1(), opts.clip_size + 2 * opts.filt_padding, opts.samplerate, opts.freq_min, opts.freq_max, opts.freq_wid);
            }
        }
#pragma omp for
        for (int ii = 0; ii < keys.count(); ii++) {
            ClusterData* CD;
#pragma omp critical(cc2)
            {
                int key = keys[ii];
                if ((timer.elapsed() > 3000) || (first)) {
                    qDebug().noquote() << QString("Computing template for cluster %1/%2").arg(key).arg(MLCompute::max(keys.toVector()));
                    timer.restart();
                    first = false;
                }
                CD = &cluster_data[key];
            }
            /// Witold -- this is a bad pitfall. I needed to use a new instance of DiskReadMda32 here because accessing the same instance across multiple threads was causing WRONG results --very bad -- please help make DiskReadMda32 thread safe
            DiskReadMda32 X0(timeseries);
            compute_template(CD->template0, X0, CD->times, opts, BP);
        }
    }

    qDebug().noquote() << "Assembling templates";
    int K = MLCompute::max(labels);
    Mda32 templates(M, T, K);
    for (int k = 1; k <= K; k++) {
        if (cluster_data.contains(k)) {
            templates.setChunk(cluster_data[k].template0, 0, 0, k - 1);
        }
    }

    templates.write32(templates_out);

    return true;
}

namespace P_spikeview_templates {
QVector<double> subsample(const QVector<double>& X, bigint max_elements)
{
    bigint increment = X.count() / max_elements;
    if (increment < 1)
        increment = 1;
    QVector<double> Y;
    for (bigint j = 0; (j < max_elements) && (j * increment < X.count()); j++) {
        Y << X[j * increment];
    }
    return Y;
}

void subtract_temporal_mean(Mda32& clip)
{
    int M = clip.N1();

    QVector<double> means(M);
    means.fill(0);
    for (int t = 0; t < clip.N2(); t++) {
        for (int m = 0; m < M; m++) {
            means[m] += clip.get(m, t);
        }
    }
    for (int m = 0; m < M; m++) {
        means[m] /= clip.N2();
    }
    for (int t = 0; t < clip.N2(); t++) {
        for (int m = 0; m < M; m++) {
            clip.set(clip.get(m, t) - means[m], m, t);
        }
    }
}

void compute_template(Mda32& template0, const DiskReadMda32& X, const QVector<double>& times, P_spikeview_templates_opts opts, const Bandpass_filter_runner& BP)
{
    int M = X.N1();
    int T = opts.clip_size;
    Mda sum(M, T + 2 * opts.filt_padding);

    int Tmid = (int)((T + 1) / 2) - 1;
    for (bigint i = 0; i < times.count(); i++) {
        bigint t1 = times[i] - Tmid - opts.filt_padding;
        Mda32 tmp;
        X.readChunk(tmp, 0, t1, M, T + opts.filt_padding * 2);
        for (int t = 0; t < T + 2 * opts.filt_padding; t++) {
            for (int m = 0; m < M; m++) {
                double val = tmp.get(m, t);
                sum.set(sum.get(m, t) + val, m, t);
            }
        }
    }

    Mda32 template1;
    template1.allocate(M, T + 2 * opts.filt_padding);

    if (times.count() > 0) {
        for (bigint j = 0; j < M * (T + 2 * opts.filt_padding); j++) {
            double sum0 = sum.get(j);
            template1.set(sum0 / times.count(), j);
        }
    }

    if (opts.subtract_temporal_mean) {
        subtract_temporal_mean(template1);
    }
    if ((opts.samplerate != 0) && ((opts.freq_min != 0) || (opts.freq_max != 0))) {
        BP.apply(template1);
    }

    template0.allocate(M, T);
    for (int t = 0; t < T; t++) {
        for (int m = 0; m < M; m++) {
            template0.set(template1.get(m, t + opts.filt_padding), m, t);
        }
    }
}

Bandpass_filter_runner::Bandpass_filter_runner()
{
}

Bandpass_filter_runner::~Bandpass_filter_runner()
{
    if (data_in)
        fftwf_free(data_in);
    if (data_out)
        fftwf_free(data_out);
    if (kernel0)
        free(kernel0);
    //delete p_fft;
    //delete p_ifft;
}
void Bandpass_filter_runner::init(bigint M_in, bigint N_in, double samplerate, double freq_min, double freq_max, double freq_wid)
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
void Bandpass_filter_runner::apply(Mda32& chunk) const
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
void Bandpass_filter_runner::define_kernel(bigint N, double* kernel, double samplefreq, double freq_min, double freq_max, double freq_wid)
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
}
