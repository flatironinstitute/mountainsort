#include "neighborhoodsorter.h"
#include "detect_events.h"
#include "pca.h"
#include "sort_clips.h"
#include "consolidate_clusters.h"

#include <QFile>
#include <QCoreApplication>
//#include <cachemanager.h>

class NeighborhoodSorterPrivate {
public:
    NeighborhoodSorter* q;
    P_mountainsort3_opts m_opts;

    bigint m_max_ram_bytes = 1e9;
    bigint m_M = 0;
    QVector<double> m_times;
    QVector<DiskBackedMda32> m_accumulated_clips;
    QVector<Mda32> m_accumulated_clips_buffer;
    double m_num_bytes = 0;

    QVector<int> m_labels;
    Mda32 m_templates;

    bigint size_of_accumulated_clips_buffer();
    void clear_accumulated_clips_buffer();

    void dimension_reduce_clips(Mda32& ret, const Mda32& clips, bigint num_features_per_channel, bigint max_samples);
    static Mda32 compute_templates_from_clips(const Mda32& clips, const QVector<int>& labels, int num_threads);
    static QVector<double> get_subarray(const QVector<double>& X, const QVector<bigint>& inds);
    static QVector<int> get_subarray(const QVector<int>& X, const QVector<bigint>& inds);
    static Mda32 get_subclips(const Mda32& clips, const QVector<bigint>& inds);
    static void get_clip(Mda32& clip, const Mda32& X, const QList<int>& channels, bigint t0, int T);
};

NeighborhoodSorter::NeighborhoodSorter()
{
    d = new NeighborhoodSorterPrivate;
    d->q = this;
}

NeighborhoodSorter::~NeighborhoodSorter()
{
    delete d;
}

void NeighborhoodSorter::setOptions(P_mountainsort3_opts opts)
{
    d->m_opts = opts;
}

void NeighborhoodSorter::setMaxRAM(bigint max_ram_bytes)
{
    d->m_max_ram_bytes = max_ram_bytes;
}

void NeighborhoodSorter::addTimeChunk(bigint t, const Mda32& X, const QList<int>& channels, bigint padding_left, bigint padding_right)
{
    d->m_M = channels.count();
    if (channels.count() == 0)
        return;
    int central_channel = channels[0];
    int T = d->m_opts.clip_size;
    //int Tmid = (int)((T + 1) / 2) - 1;
    QVector<double> X0;
    for (bigint i = 0; i < X.N2(); i++) {
        X0 << X.value(central_channel - 1, i);
    }

    QVector<double> times0 = detect_events(X0, d->m_opts.detect_threshold, d->m_opts.detect_interval, d->m_opts.detect_sign);
    QVector<double> times1; // only the times that fit within the proper time chunk (excluding padding)
    for (bigint i = 0; i < times0.count(); i++) {
        double t0 = times0[i];
        if ((0 <= t0 - padding_left) && (t0 - padding_left < X.N2() - padding_left - padding_right)) {
            times1 << t0;
        }
    }

    Mda32 clips0(X.N1(), T, times1.count());
    for (bigint i = 0; i < times1.count(); i++) {
        double t0 = times1[i];
        Mda32 clip0;
        d->get_clip(clip0, X, channels, t0, T);
        clips0.setChunk(clip0, 0, 0, i);
        d->m_times << t0 - padding_left + t;
    }
    d->m_accumulated_clips_buffer << clips0;
    if (d->size_of_accumulated_clips_buffer() > d->m_max_ram_bytes)
        d->clear_accumulated_clips_buffer();
}

void NeighborhoodSorter::sort(int num_threads)
{
    int T = d->m_opts.clip_size;

    // get the clips
    Mda32 clips;
    clips.allocate(d->m_M, T, d->m_times.count());

    if (d->m_accumulated_clips.count() > 0) {
        d->clear_accumulated_clips_buffer(); //important!

        bigint ii = 0;
        for (bigint j = 0; j < d->m_accumulated_clips.count(); j++) {
            Mda32 clips0;
            d->m_accumulated_clips[j].retrieve(clips0);
            d->m_accumulated_clips[j].remove(); //remove the temporary file
            clips.setChunk(clips0, 0, 0, ii);
            ii += clips0.N3();
        }
        if (ii != d->m_times.count()) {
            qWarning() << "Unexpected error putting together the clips in NeighborhoodSorter::sort()." << ii << d->m_times.count();
            abort();
        }
    }
    else {
        bigint ii = 0;
        for (bigint j = 0; j < d->m_accumulated_clips_buffer.count(); j++) {
            Mda32 clips0 = d->m_accumulated_clips_buffer[j];
            clips.setChunk(clips0, 0, 0, ii);
            ii += clips0.N3();
        }
        if (ii != d->m_times.count()) {
            qWarning() << "Unexpected error putting together the clips in NeighborhoodSorter::sort() (*)." << ii << d->m_times.count();
            abort();
        }
        d->m_accumulated_clips_buffer.clear();
    }

    // dimension reduce clips
    Mda32 reduced_clips;
    d->dimension_reduce_clips(reduced_clips, clips, d->m_opts.num_features_per_channel, d->m_opts.max_pca_samples);

    // Sort
    Sort_clips_opts ooo;
    ooo.max_samples = d->m_opts.max_pca_samples;
    ooo.num_features = d->m_opts.num_features;
    d->m_labels = sort_clips(reduced_clips, ooo);
    qDebug().noquote() << QString("Sorted %1 clips and found %2 clusters").arg(reduced_clips.N3()).arg(MLCompute::max(d->m_labels));

    // Compute templates
    d->m_templates = d->compute_templates_from_clips(clips, d->m_labels, num_threads);

    // Consolidate clusters
    if (d->m_opts.consolidate_clusters) {
        Consolidate_clusters_opts oo;
        QMap<int, int> label_map = consolidate_clusters(d->m_templates, oo);
        QVector<bigint> inds_to_keep;
        for (bigint i = 0; i < d->m_labels.count(); i++) {
            int k0 = d->m_labels[i];
            if ((k0 > 0) && (label_map[k0] > 0))
                inds_to_keep << i;
        }
        qDebug().noquote() << QString("Consolidate. Keeping %1 of %2 events").arg(inds_to_keep.count()).arg(d->m_labels.count());
        d->m_times = d->get_subarray(d->m_times, inds_to_keep);
        d->m_labels = d->get_subarray(d->m_labels, inds_to_keep);
        clips = d->get_subclips(clips, inds_to_keep);
        reduced_clips = d->get_subclips(reduced_clips, inds_to_keep);
        for (bigint i = 0; i < d->m_labels.count(); i++) {
            int k0 = label_map[d->m_labels[i]];
            d->m_labels[i] = k0;
        }
    }

    // Compute templates
    d->m_templates = d->compute_templates_from_clips(clips, d->m_labels, num_threads);
}

QVector<double> NeighborhoodSorter::times() const
{
    return d->m_times;
}

QVector<int> NeighborhoodSorter::labels() const
{
    return d->m_labels;
}

Mda32 NeighborhoodSorter::templates() const
{
    return d->m_templates;
}

bigint NeighborhoodSorterPrivate::size_of_accumulated_clips_buffer()
{
    bigint ret = 0;
    for (bigint i = 0; i < m_accumulated_clips_buffer.count(); i++) {
        ret += m_accumulated_clips_buffer[i].totalSize() * sizeof(float);
    }
    return ret;
}

void NeighborhoodSorterPrivate::clear_accumulated_clips_buffer()
{
    bigint L = 0;
    for (bigint i = 0; i < m_accumulated_clips_buffer.count(); i++) {
        L += m_accumulated_clips_buffer[i].N3();
    }
    if (L > 0) {
        Mda32 clips0(m_M, m_opts.clip_size, L);
        bigint jj = 0;
        for (bigint i = 0; i < m_accumulated_clips_buffer.count(); i++) {
            clips0.setChunk(m_accumulated_clips_buffer[i], 0, 0, jj);
            jj += m_accumulated_clips_buffer[i].N3();
        }
        DiskBackedMda32 tmp(clips0);
        m_accumulated_clips << tmp;
        m_accumulated_clips_buffer.clear();
    }
}

void NeighborhoodSorterPrivate::dimension_reduce_clips(Mda32& ret, const Mda32& clips, bigint num_features_per_channel, bigint max_samples)
{
    bigint M = clips.N1();
    bigint T = clips.N2();
    bigint L = clips.N3();
    const float* clips_ptr = clips.constDataPtr();

    qDebug().noquote() << QString("Dimension reduce clips %1x%2x%3").arg(M).arg(T).arg(L);

    ret.allocate(M, num_features_per_channel, L);
    float* retptr = ret.dataPtr();
    for (bigint m = 0; m < M; m++) {
        Mda32 reshaped(T, L);
        float* reshaped_ptr = reshaped.dataPtr();
        bigint aa = 0;
        bigint bb = m;
        for (bigint i = 0; i < L; i++) {
            for (bigint t = 0; t < T; t++) {
                //reshaped.set(clips.get(bb),aa);
                reshaped_ptr[aa] = clips_ptr[bb];
                aa++;
                bb += M;
                //reshaped.setValue(clips.value(m, t, i), t, i);
            }
        }
        Mda32 CC, FF, sigma;
        pca_subsampled(CC, FF, sigma, reshaped, num_features_per_channel, false, max_samples);
        float* FF_ptr = FF.dataPtr();
        aa = 0;
        bb = m;
        for (bigint i = 0; i < L; i++) {
            for (bigint a = 0; a < num_features_per_channel; a++) {
                //ret.setValue(FF.value(a, i), m, a, i);
                //ret.set(FF.get(aa),bb);
                retptr[bb] = FF_ptr[aa];
                aa++;
                bb += M;
            }
        }
    }
}

Mda32 NeighborhoodSorterPrivate::compute_templates_from_clips(const Mda32& clips, const QVector<int>& labels, int num_threads)
{
    int M = clips.N1();
    int T = clips.N2();
    bigint L = clips.N3();
    int K = MLCompute::max(labels);
    Mda sums(M, T, K);
    QVector<bigint> counts(K, 0);
#pragma omp parallel num_threads(num_threads)
    {
        Mda local_sums(M, T, K);
        QVector<bigint> local_counts(K, 0);
#pragma omp for
        for (bigint i = 0; i < L; i++) {
            int k0 = labels[i];
            if (k0 > 0) {
                Mda32 tmp;
                clips.getChunk(tmp, 0, 0, i, M, T, 1);
                for (int t = 0; t < T; t++) {
                    for (int m = 0; m < M; m++) {
                        local_sums.set(local_sums.get(m, t, k0 - 1) + tmp.get(m, t), m, t, k0 - 1);
                    }
                }
                local_counts[k0 - 1]++;
            }
        }
#pragma omp critical(set_sums_and_counts)
        for (int kk = 1; kk <= K; kk++) {
            counts[kk - 1] += local_counts[kk - 1];
            Mda tmp;
            local_sums.getChunk(tmp, 0, 0, kk - 1, M, T, 1);
            for (int t = 0; t < T; t++) {
                for (int m = 0; m < M; m++) {
                    sums.set(sums.get(m, t, kk - 1) + tmp.get(m, t), m, t, kk - 1);
                }
            }
        }
    }
    Mda32 templates(M, T, K);
    for (int kk = 1; kk <= K; kk++) {
        if (counts[kk - 1]) {
            for (int t = 0; t < T; t++) {
                for (int m = 0; m < M; m++) {
                    templates.set(sums.get(m, t, kk - 1) / counts[kk - 1], m, t, kk - 1);
                }
            }
        }
    }
    return templates;
}

QVector<double> NeighborhoodSorterPrivate::get_subarray(const QVector<double>& X, const QVector<bigint>& inds)
{
    QVector<double> ret(inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        ret[i] = X[inds[i]];
    }
    return ret;
}

QVector<int> NeighborhoodSorterPrivate::get_subarray(const QVector<int>& X, const QVector<bigint>& inds)
{
    QVector<int> ret(inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        ret[i] = X[inds[i]];
    }
    return ret;
}

Mda32 NeighborhoodSorterPrivate::get_subclips(const Mda32& clips, const QVector<bigint>& inds)
{
    Mda32 ret(clips.N1(), clips.N2(), inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        Mda32 tmp;
        clips.getChunk(tmp, 0, 0, inds[i], clips.N1(), clips.N2(), 1);
        ret.setChunk(tmp, 0, 0, i);
    }
    return ret;
}

DiskBackedMda32::DiskBackedMda32()
{
}

DiskBackedMda32::DiskBackedMda32(const Mda32& X)
{
    this->store(X);
}

DiskBackedMda32::~DiskBackedMda32()
{
    //don't remove here -- user must call remove() explicitly
}

void DiskBackedMda32::store(const Mda32& X)
{
    if (m_tmp_path.isEmpty()) {
        qWarning() << "Unable to use DiskBackedMda32 without a temporary path specified";
        abort();
        //m_tmp_path = CacheManager::globalInstance()->makeLocalFile(MLUtil::makeRandomId() + ".DiskBackedMda32.mda");
        //CacheManager::globalInstance()->setTemporaryFileExpirePid(m_tmp_path, QCoreApplication::applicationPid());
    }
    X.write32(m_tmp_path);
}

void DiskBackedMda32::retrieve(Mda32& X) const
{
    X.read(m_tmp_path);
}

void DiskBackedMda32::remove()
{
    if (!m_tmp_path.isEmpty()) {
        QFile::remove(m_tmp_path);
    }
}

void NeighborhoodSorterPrivate::get_clip(Mda32& clip, const Mda32& X, const QList<int>& channels, bigint t0, int T)
{
    int M = channels.count();
    int Tmid = (int)((T + 1) / 2) - 1;
    clip.allocate(M, T);
    for (int t = 0; t < T; t++) {
        for (int m = 0; m < M; m++) {
            clip.set(X.value(channels[m] - 1, t0 + t - Tmid), m, t);
        }
    }
}
