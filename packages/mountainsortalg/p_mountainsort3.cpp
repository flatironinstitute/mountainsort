#include "fitstagecomputer.h"
#include "p_mountainsort3.h"

#include <QTime>
#include <diskreadmda32.h>
#include <mda.h>
#include <mda32.h>
#include "omp.h"
#include "neighborhoodsorter.h"
#include "get_sort_indices.h"
#include "merge_across_channels.h"
#include "globaltemplatecomputer.h"
#include "fit_stage.h"
#include "compute_templates_0.h"

namespace MountainSort3 {

struct TimeChunkInfo {
    bigint t_padding; //t1 corresponds to index t_padding (because we have padding/overlap)
    bigint t1; //starting timepoint (corresponding to index t_padding)
    bigint size; //number of timepoints (excluding padding on left and right)
};

QList<int> get_channels_from_geom(const Mda& geom, bigint m, double adjacency_radius);
//void extract_channels(Mda32& ret, const Mda32& X, const QList<int>& channels);
QList<TimeChunkInfo> get_time_chunk_infos(bigint M, bigint N, int num_threads, bigint min_num_chunks);

QList<QList<int> > get_neighborhood_batches(int M, int batch_size)
{
    QList<QList<int> > batches;
    for (int m = 1; m <= M; m += batch_size) {
        QList<int> batch;
        for (int m2 = m; (m2 < m + batch_size) && (m2 <= M); m2++)
            batch << m2;
        batches << batch;
    }
    return batches;
}

QVector<double> reorder(const QVector<double>& X, const QList<bigint>& inds);
QVector<int> reorder(const QVector<int>& X, const QList<bigint>& inds);
QVector<double> get_subarray(const QVector<double>& X, const QVector<bigint>& inds);
QVector<int> get_subarray(const QVector<int>& X, const QVector<bigint>& inds);

struct ProgressReporter {
    void startProcessingStep(QString name, double total_expected_bytes_to_process = 0)
    {
        qDebug().noquote() << "";
        qDebug().noquote() << QString("Starting  [%1]...").arg(name);
        m_processing_step_name = name;
        m_processing_step_timer.start();
        m_progress_timer.start();
        m_bytes_processed = 0;
        m_bytes_read = 0;
        m_read_time_sec = 0;
        m_total_expected_bytes_to_process = total_expected_bytes_to_process;
    }
    void addBytesProcessed(double bytes)
    {
        double elapsed_sec = m_progress_timer.elapsed() * 1.0 / 1000;
        double mb = bytes / 1e6;
        m_bytes_processed += bytes;
        QString status;
        if (m_total_expected_bytes_to_process) {
            status = QString("%1%").arg(m_bytes_processed * 100.0 / m_total_expected_bytes_to_process);
        }
        qDebug().noquote() << QString("Processed %2 MB in %3 sec (%4 MB/sec): Total %5 MB [%1] (%6)...").arg(m_processing_step_name).arg(mb).arg(elapsed_sec).arg(mb / elapsed_sec).arg(m_bytes_processed / 1e6).arg(status);
        m_progress_timer.restart();
    }
    void addBytesRead(double bytes)
    {
        double elapsed_sec = m_progress_timer.elapsed() * 1.0 / 1000;
        double mb = bytes / 1e6;
        m_read_time_sec += elapsed_sec;
        m_bytes_read += bytes;
        qDebug().noquote() << QString("Read %2 MB in %3 sec (%4 MB/sec): Total %5 MB [%1]...").arg(m_processing_step_name).arg(mb).arg(elapsed_sec).arg(mb / elapsed_sec).arg(m_bytes_read / 1e6);
        m_progress_timer.restart();
    }

    void endProcessingStep()
    {
        double elapsed_sec = m_processing_step_timer.elapsed() * 1.0 / 1000;
        double mb = m_bytes_processed / 1e6;
        qDebug().noquote() << QString("Elapsed time -- %2 sec (%3 MB/sec) [%1]...\n").arg(m_processing_step_name).arg(elapsed_sec).arg(mb / elapsed_sec);
        qDebug().noquote() << "";
        m_summary_text += QString("%2 sec (includes %3 sec reading) [%1]\n").arg(m_processing_step_name).arg(elapsed_sec).arg(m_read_time_sec);
    }
    void printSummaryText()
    {
        qDebug().noquote() << "";
        qDebug().noquote() << m_summary_text;
        qDebug().noquote() << "";
    }

private:
    QString m_processing_step_name;
    QTime m_processing_step_timer;
    QTime m_progress_timer;
    double m_total_expected_bytes_to_process = 0;
    double m_bytes_processed = 0;
    double m_bytes_read = 0;
    double m_read_time_sec = 0;
    QString m_summary_text;
};
}

using namespace MountainSort3;

bool p_mountainsort3(QString timeseries, QString geom, QString firings_out, QString temp_path, const P_mountainsort3_opts& opts)
{
    (void)temp_path;
    int tot_threads = omp_get_max_threads();
    omp_set_nested(1);

    ProgressReporter PR;

    DiskReadMda32 X(timeseries);
    bigint M = X.N1();
    bigint N = X.N2();

    bigint t_start = opts.t1;
    bigint t_end = opts.t2;
    if (t_start < 0)
        t_start = 0;
    if (t_end < 0)
        t_end = N - 1;
    N = t_end - t_start + 1;

    qDebug().noquote() << QString("Starting sort: %1 channels, %2 timepoints").arg(M).arg(N);

    Mda Geom(geom);

    int num_simultaneous_neighborhoods = qMin(tot_threads, (int)M);
    int num_threads_within_neighborhoods = qMax(1, tot_threads / num_simultaneous_neighborhoods);
    QList<QList<int> > neighborhood_batches = get_neighborhood_batches(M, num_simultaneous_neighborhoods);
    QList<TimeChunkInfo> time_chunk_infos;

    // The amount of clips data (in the neighborhood sorters) that is permitted to stay in memory
    bigint clips_RAM = 20e9; //how to set this?
    bigint clips_RAM_per_neighborhood = clips_RAM / M;

    QMap<int, NeighborhoodSorter*> neighborhood_sorters;
    QMap<int, QList<int> > neighborhood_channels;
    for (bigint m = 1; m <= M; m++) {
        neighborhood_sorters[m] = new NeighborhoodSorter;
        neighborhood_sorters[m]->setOptions(opts);
        neighborhood_channels[m] = get_channels_from_geom(Geom, m, opts.adjacency_radius);
        neighborhood_sorters[m]->setMaxRAM(clips_RAM_per_neighborhood); //allow each neighborhood sorter to store a certain amount of clips in memory -- after that they store it temporarily on the disk
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Read data and add to neighborhood sorters
    PR.startProcessingStep("Read data and add to neighborhood sorters", M * N * sizeof(float));
    time_chunk_infos = get_time_chunk_infos(M, N, neighborhood_batches.count(), 1);
    for (bigint i = 0; i < time_chunk_infos.count(); i++) {
        TimeChunkInfo TCI = time_chunk_infos[i];
        Mda32 time_chunk;
        X.readChunk(time_chunk, 0, t_start + TCI.t1 - TCI.t_padding, M, TCI.size + 2 * TCI.t_padding);
        PR.addBytesRead(time_chunk.totalSize() * sizeof(float));
        for (int j = 0; j < neighborhood_batches.count(); j++) {
            QList<int> neighborhoods = neighborhood_batches[j];
            double bytes0 = 0;
#pragma omp parallel for num_threads(num_simultaneous_neighborhoods)
            for (int k = 0; k < neighborhoods.count(); k++) {
                int m;
#pragma omp critical(m1)
                {
                    m = neighborhoods[k]; //I don't know why this needs to be in a critical section
                }
                //extract_channels(neighborhood_time_chunk, time_chunk, neighborhood_channels[m]);
                neighborhood_sorters[m]->addTimeChunk(TCI.t1, time_chunk, neighborhood_channels[m], TCI.t_padding, TCI.t_padding);
#pragma omp critical(a1)
                {
                    bytes0 += time_chunk.N2() * sizeof(float);
                }
            }
            PR.addBytesProcessed(bytes0);
        }
    }
    PR.endProcessingStep();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sort in each neighborhood sorter
    PR.startProcessingStep("Sort in each neighborhood", M * N * sizeof(float));
    for (int j = 0; j < neighborhood_batches.count(); j++) {
        QList<int> neighborhoods = neighborhood_batches[j];
        double bytes0 = 0;
#pragma omp parallel for num_threads(num_simultaneous_neighborhoods)
        for (int k = 0; k < neighborhoods.count(); k++) {
            int m;
#pragma omp critical(m2)
            {
                m = neighborhoods[k]; //I don't know why this needs to be in a critical section
            }
            neighborhood_sorters[m]->sort(num_threads_within_neighborhoods);
#pragma omp critical(a1)
            {
                bytes0 += N * sizeof(float);
            }
        }
        PR.addBytesProcessed(bytes0);
    }
    PR.endProcessingStep();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Collect the events and sort them by time
    PR.startProcessingStep("Collect events and sort them by time");
    QVector<double> times;
    QVector<int> labels;
    QVector<int> central_channels;
    int K_offset = 0;
    for (bigint m = 1; m <= M; m++) {
        QVector<double> times0 = neighborhood_sorters[m]->times();
        QVector<int> labels0 = neighborhood_sorters[m]->labels();
        for (bigint a = 0; a < labels0.count(); a++)
            labels0[a] += K_offset;
        times.append(times0);
        labels.append(labels0);
        K_offset = qMax(K_offset, MLCompute::max(labels0));
        central_channels.append(QVector<int>(neighborhood_sorters[m]->labels().count(), m));
    }
    QList<bigint> sort_inds = get_sort_indices_bigint(times);
    times = reorder(times, sort_inds);
    labels = reorder(labels, sort_inds);
    central_channels = reorder(central_channels, sort_inds);
    PR.endProcessingStep();

    qDeleteAll(neighborhood_sorters);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute global templates
    PR.startProcessingStep("Compute global templates", M * N * sizeof(float));
    time_chunk_infos = get_time_chunk_infos(M, N, 1, 1);
    Mda32 templates;
    {
        GlobalTemplateComputer GTC;
        GTC.setNumThreads(tot_threads);
        GTC.setClipSize(opts.clip_size);
        GTC.setTimesLabels(times, labels);
        for (bigint i = 0; i < time_chunk_infos.count(); i++) {
            TimeChunkInfo TCI = time_chunk_infos[i];
            Mda32 time_chunk;
            X.readChunk(time_chunk, 0, t_start + TCI.t1 - TCI.t_padding, M, TCI.size + 2 * TCI.t_padding);
            PR.addBytesRead(time_chunk.totalSize() * sizeof(float));
            GTC.processTimeChunk(TCI.t1, time_chunk, TCI.t_padding, TCI.t_padding);
            PR.addBytesProcessed(time_chunk.totalSize() * sizeof(float));
        }
        templates = GTC.templates();
    }
    PR.endProcessingStep();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Merge across channels
    if (opts.merge_across_channels) {
        PR.startProcessingStep("Merge across channels");
        {
            Merge_across_channels_opts oo;
            oo.clip_size = opts.clip_size;
            merge_across_channels(times, labels, central_channels, templates, oo);
        }
        PR.endProcessingStep();
    }

    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute global templates
    PR.startProcessingStep("Compute global templates");
    time_chunk_infos = get_time_chunk_infos(M, N, 1, 1);
    {
        GlobalTemplateComputer GTC;
        GTC.setNumThreads(tot_threads);
        GTC.setClipSize(opts.clip_size);
        GTC.setTimesLabels(times, labels);
        for (bigint i = 0; i < time_chunk_infos.count(); i++) {
            TimeChunkInfo TCI = time_chunk_infos[i];
            Mda32 time_chunk;
            X.readChunk(time_chunk, 0, t_start + TCI.t1 - TCI.t_padding, M, TCI.size + 2 * TCI.t_padding);
            PR.addBytesRead(time_chunk.totalSize() * sizeof(float));
            GTC.processTimeChunk(TCI.t1,time_chunk, TCI.t_padding, TCI.t_padding);
            PR.addBytesProcessed(time_chunk.totalSize() * sizeof(float));
        }
        templates = GTC.templates();
    }
    PR.endProcessingStep();
    */

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    if (opts.fit_stage) {
        PR.startProcessingStep("Fit stage", M * N * sizeof(float));
        time_chunk_infos = get_time_chunk_infos(M, N, tot_threads, tot_threads);
        FitStageComputer FSC;
        FSC.setTemplates(templates);
        FSC.setTimesLabels(times, labels);
        for (bigint i = 0; i < time_chunk_infos.count(); i += tot_threads) {
            QList<TimeChunkInfo> infos = time_chunk_infos.mid(i, tot_threads);
            QList<Mda32> time_chunks;
            double bytes0 = 0;
            for (int j = 0; j < infos.count(); j++) {
                Mda32 time_chunk;
                X.readChunk(time_chunk, 0, t_start + infos[j].t1 - infos[j].t_padding, M, infos[j].size + 2 * infos[j].t_padding);
                time_chunks << time_chunk;
                bytes0 += time_chunk.totalSize() * sizeof(float);
            }
            PR.addBytesRead(bytes0);
#pragma omp parallel for num_threads(tot_threads)
            for (int j = 0; j < infos.count(); j++) {
                FSC.processTimeChunk(infos[j].t1, time_chunks[j], infos[j].t_padding, infos[j].t_padding);
            }
            PR.addBytesProcessed(bytes0);
        }
        FSC.finalize();
        QVector<bigint> event_inds_to_use = FSC.eventIndicesToUse();
        qDebug().noquote() << QString("Using %1 of %2 events (%3%) after fit stage").arg(event_inds_to_use.count()).arg(times.count()).arg(event_inds_to_use.count() * 1.0 / times.count() * 100);
        times = get_subarray(times, event_inds_to_use);
        labels = get_subarray(labels, event_inds_to_use);
        central_channels = get_subarray(central_channels, event_inds_to_use);
        PR.endProcessingStep();
    }

    /*
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// FIT STAGE
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        qDebug().noquote() << "******* Fit stage...";
        QTime timer;
        timer.start();
        QTime progress_timer;
        progress_timer.start();
        int num_time_chunks_read = 0;

        QVector<bigint> inds_to_use;

        QTime timerA;
        timerA.start();
        int num_time_threads=1;
        int chunk_overlap_size=1000;
        int progress_msec=2000;
        for (int it = 0; it < time_chunk_infos.count(); it += num_time_threads) {
            QList<TimeChunkInfo> infos = time_chunk_infos.mid(it, num_time_threads);
            QVector<Mda32> time_chunks(infos.count());
            double num_bytes_read = 0;
            for (int j = 0; j < infos.count(); j++) {
                TimeChunkInfo TCI = infos[j];
                X.readChunk(time_chunks[j], 0, t_start + TCI.t1 - TCI.t_padding, M, TCI.size + 2 * TCI.t_padding);
                num_time_chunks_read++;
                num_bytes_read += time_chunks[j].totalSize() * 4;
            }
            {
                double elapsed_sec = timerA.elapsed() * 1.0 / 1000;
                qDebug().noquote() << QString("Elapsed for reading: %1 (%2 MB/sec)").arg(timerA.restart()).arg(num_bytes_read / 1e6 / (elapsed_sec));
            }
#pragma omp parallel for
            for (int j = 0; j < infos.count(); j++) {
                TimeChunkInfo TCI = infos[j];
                Mda32 time_chunk = time_chunks[j];

                QVector<bigint> local_inds;
                QVector<double> local_times;
                QVector<int> local_labels;

                bigint ii = 0;
                while ((ii < times.count()) && (times[ii] < TCI.t1))
                    ii++;
                while ((ii < times.count()) && (times[ii] < TCI.t1 + TCI.size)) {
                    local_inds << ii;
                    local_times << times[ii] - TCI.t1 + chunk_overlap_size;
                    local_labels << labels[ii];
                    ii++;
                }
                Fit_stage_opts oo;
                QVector<bigint> local_inds_to_use = fit_stage(time_chunk, local_times, local_labels, templates, oo);
#pragma omp critical(fit_stage_set_inds_to_use1)
                {
                    for (bigint a = 0; a < local_inds_to_use.count(); a++) {
                        inds_to_use << local_inds[local_inds_to_use[a]];
                    }
                }
            }
            {
                double elapsed_sec = timerA.elapsed() * 1.0 / 1000;
                qDebug().noquote() << QString("Elapsed for fit stage: %1 (%2 MB/sec)").arg(timerA.restart()).arg(num_bytes_read / 1e6 / (elapsed_sec));
            }
            if (progress_timer.elapsed() > progress_msec) {
                qDebug().noquote() << QString("Fit stage: Processed %1 time chunks of %2...").arg(num_time_chunks_read).arg(time_chunk_infos.count());
                progress_timer.restart();
            }
        }

        qSort(inds_to_use);
        qDebug().noquote() << QString("Fit stage: Using %1 of %2 events").arg(inds_to_use.count()).arg(times.count());

        QVector<double> times2(inds_to_use.count());
        QVector<int> labels2(inds_to_use.count());
        QVector<int> central_channels2(inds_to_use.count());
        for (bigint a = 0; a < inds_to_use.count(); a++) {
            times2[a] = times[inds_to_use[a]];
            labels2[a] = labels[inds_to_use[a]];
            central_channels2[a] = central_channels[inds_to_use[a]];
        }
        times = times2;
        labels = labels2;
        central_channels = central_channels2;

        qDebug().noquote() << "Computing global templates...";
        //templates = compute_templates_in_parallel(X, times, labels, opts.clip_size);

        qDebug().noquote() << "Elapsed (Fit stage): " << timer.elapsed() * 1.0 / 1000;
    }
    */

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Create and write firings output
    PR.startProcessingStep("Create and write firings output");
    {
        bigint L = times.count();
        Mda firings(3, L);
        for (bigint i = 0; i < L; i++) {
            firings.set(central_channels[i], 0, i);
            firings.set(t_start + times[i], 1, i);
            firings.set(labels[i], 2, i);
        }

        firings.write64(firings_out);
    }
    PR.endProcessingStep();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    PR.printSummaryText();

    return true;
}

namespace MountainSort3 {
QList<int> get_channels_from_geom(const Mda& geom, bigint m, double adjacency_radius)
{
    bigint M = geom.N2();
    QList<int> ret;
    ret << m;
    for (int m2 = 1; m2 <= M; m2++) {
        if (m2 != m) {
            double sumsqr = 0;
            for (int i = 0; i < geom.N1(); i++) {
                double tmp = geom.value(i, m - 1) - geom.value(i, m2 - 1);
                sumsqr += tmp * tmp;
            }
            double dist = sqrt(sumsqr);
            if (dist <= adjacency_radius)
                ret << m2;
        }
    }
    return ret;
}

/*
void extract_channels(Mda32& ret, const Mda32& X, const QList<int>& channels)
{
    bigint M2 = channels.count();
    bigint N = X.N2();
    ret.allocate(M2, N);
    for (bigint t = 0; t < N; t++) {
        for (bigint ii = 0; ii < M2; ii++) {
            ret.set(X.get(channels[ii] - 1, t), ii, t);
        }
    }
}
*/

QList<TimeChunkInfo> get_time_chunk_infos(bigint M, bigint N, int num_simultaneous, bigint min_num_chunks)
{
    bigint RAM_available_for_chunks_bytes = 1e9;
    bigint chunk_size = RAM_available_for_chunks_bytes / (M * num_simultaneous * 4);
    if (chunk_size > N)
        chunk_size = N;
    if (N < chunk_size * min_num_chunks) {
        chunk_size = N / min_num_chunks;
    }
    if (chunk_size < 20000)
        chunk_size = 20000;
    bigint chunk_overlap_size = 1000;

    // Prepare the information on the time chunks
    QList<TimeChunkInfo> time_chunk_infos;
    for (bigint t = 0; t < N; t += chunk_size) {
        TimeChunkInfo info;
        info.t1 = t;
        info.t_padding = chunk_overlap_size;
        info.size = chunk_size;
        if (t + info.size > N)
            info.size = N - t;
        time_chunk_infos << info;
    }

    return time_chunk_infos;
}

QVector<double> reorder(const QVector<double>& X, const QList<bigint>& inds)
{
    QVector<double> ret(inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        ret[i] = X[inds[i]];
    }
    return ret;
}

QVector<int> reorder(const QVector<int>& X, const QList<bigint>& inds)
{
    QVector<int> ret(inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        ret[i] = X[inds[i]];
    }
    return ret;
}

QVector<double> get_subarray(const QVector<double>& X, const QVector<bigint>& inds)
{
    QVector<double> ret(inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        ret[i] = X[inds[i]];
    }
    return ret;
}

QVector<int> get_subarray(const QVector<int>& X, const QVector<bigint>& inds)
{
    QVector<int> ret(inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        ret[i] = X[inds[i]];
    }
    return ret;
}
}
