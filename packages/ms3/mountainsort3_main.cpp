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

#include "mdaio.h"
#include "clparams.h"
#include "mlutil.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QCoreApplication>
#include "mountainsort3_main.h"
//#include "p_multineighborhood_sort.h"
//#include "p_preprocess.h"
#include "p_run_metrics_script.h"
#include "p_spikeview_metrics.h"
#include "p_spikeview_templates.h"
#include "omp.h"
//#include "p_synthesize_timeseries.h"
#include "p_combine_firing_segments.h"
#include "p_extract_firings.h"
#include "p_concat_timeseries.h"
#include "p_banjoview_cross_correlograms.h"
#include "p_create_multiscale_timeseries.h"
#include "p_bandpass_filter.h"
#include "p_whiten.h"
#include "p_extract_clips.h"
#include "p_compute_templates.h"
#include "p_create_firings.h"
#include "p_combine_firings.h"
#include "p_apply_timestamp_offset.h"
#include "p_link_segments.h"
#include "p_cluster_metrics.h"
#include "p_isolation_metrics.h"
#include "p_concat_firings.h"
#include "p_concat_timeseries.h"
#include "p_split_firings.h"
#include "p_load_test.h"
#include "p_compute_amplitudes.h"
#include "p_confusion_matrix.h"
#include "p_reorder_labels.h"
#include "p_mask_out_artifacts.h"
#include "p_mv_compute_templates.h"
#include "p_mv_compute_amplitudes.h"
#include "p_mv_discrimhist.h"

QJsonObject get_spec()
{
    QJsonArray processors;

    {
        ProcessorSpec X("ms3.run_metrics_script", "0.1");
        X.addInputs("metrics", "script");
        X.addOutputs("metrics_out");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("spikeview.metrics1", "0.13");
        X.addInputs("firings");
        X.addOutputs("metrics_out");
        X.addRequiredParameter("samplerate");
        processors.push_back(X.get_spec());
    }
#ifndef NO_FFTW3
    {
        ProcessorSpec X("spikeview.templates", "0.16");
        X.addInputs("timeseries", "firings");
        X.addOutputs("templates_out");
        X.addOptionalParameter("clip_size", "", 100);
        X.addOptionalParameter("max_events_per_template", "", 500);
        X.addOptionalParameters("samplerate", "freq_min", "freq_max");
        X.addOptionalParameter("subtract_temporal_mean", "", "false");
        processors.push_back(X.get_spec());
    }
#endif
    {
        ProcessorSpec X("banjoview.cross_correlograms", "0.11");
        X.addInputs("firings");
        X.addOutputs("correlograms_out");
        X.addRequiredParameters("samplerate", "max_dt_msec", "bin_size_msec");
        X.addOptionalParameter("mode", "autocorrelograms or matrix_of_cross_correlograms", "autocorrelograms");
        X.addOptionalParameter("clusters", "", "");
        processors.push_back(X.get_spec());
    }
    /*
    {
        ProcessorSpec X("ms3.synthesize_timeseries", "0.12");
        X.addInputs("firings", "waveforms");
        X.addOutputs("timeseries_out");
        X.addOptionalParameter("noise_level", "", 0);
        X.addOptionalParameter("duration", "", 0);
        X.addOptionalParameter("waveform_upsample_factor", "", 13);
        processors.push_back(X.get_spec());
    }
    */
    {
        ProcessorSpec X("ms3.combine_firing_segments", "0.13");
        X.addInputs("timeseries", "firings_list");
        X.addOutputs("firings_out");
        X.addOptionalParameter("clip_size", "", 60);
        X.addOptionalParameter("match_score_threshold", "", 0.6);
        X.addOptionalParameter("offset_search_radius", "", 10);
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.extract_firings", "0.2");
        X.addInputs("firings");
        X.addOptionalInputs("metrics");
        X.addOutputs("firings_out");
        X.addOptionalParameter("exclusion_tags", "", "");
        X.addOptionalParameter("clusters", "", "");
        X.addOptionalParameter("t1", "", "");
        X.addOptionalParameter("t2", "", "");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.concat_timeseries", "0.11");
        X.addInputs("timeseries_list");
        X.addOutputs("timeseries_out");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.create_multiscale_timeseries", "0.1");
        X.addInputs("timeseries");
        X.addOutputs("timeseries_out");
        processors.push_back(X.get_spec());
    }




#ifndef NO_FFTW3
    {
        ProcessorSpec X("ms3.bandpass_filter", "0.20");
        X.addInputs("timeseries");
        X.addOutputs("timeseries_out");
        X.addRequiredParameters("samplerate", "freq_min", "freq_max");
        X.addOptionalParameter("freq_wid", "", 1000);
        X.addOptionalParameter("quantization_unit", "", 0);
        X.addOptionalParameter("subsample_factor", "", 1);
        X.can_return_requirements = true;
        processors.push_back(X.get_spec());
    }
#endif
    {
        ProcessorSpec X("ms3.whiten", "0.1");
        X.addInputs("timeseries");
        X.addOutputs("timeseries_out");
        //X.addRequiredParameters();
        X.addOptionalParameter("quantization_unit", "", 0);
        X.can_return_requirements = true;
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.compute_whitening_matrix", "0.11");
        X.addInputs("timeseries_list");
        X.addOutputs("whitening_matrix_out");
        X.addOptionalParameter("channels");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.whiten_clips", "0.1");
        X.addInputs("clips", "whitening_matrix");
        X.addOutputs("clips_out");
        //X.addRequiredParameters();
        X.addOptionalParameter("quantization_unit", "", 0);
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.apply_whitening_matrix", "0.1");
        X.addInputs("timeseries", "whitening_matrix");
        X.addOutputs("timeseries_out");
        //X.addRequiredParameters();
        X.addOptionalParameter("quantization_unit", "", 0);
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.extract_clips", "0.11");
        X.addInput("timeseries");
        X.addInput("event_times","1xL array of event timestamps");
        X.addOutputs("clips_out");
        X.addRequiredParameters("clip_size");
        X.addOptionalParameter("channels");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.mv_extract_clips", "0.1");
        X.addInput("timeseries");
        X.addInput("firings","3xL array");
        X.addOutputs("clips_out");
        X.addRequiredParameters("clip_size");
        X.addOptionalParameter("channels");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.mv_extract_clips_features", "0.1");
        X.addInputs("timeseries", "firings");
        X.addOutputs("features_out");
        X.addRequiredParameters("clip_size","num_features","subtract_mean");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.compute_templates", "0.11");
        X.addInputs("timeseries", "firings");
        X.addOutputs("templates_out");
        X.addRequiredParameters("clip_size");
        X.addOptionalParameter("clusters", "Comma-separated list of clusters to inclue", "");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.reorder_labels", "0.11");
        X.addInputs("templates", "firings");
        X.addOutputs("firings_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.create_firings", "0.1");
        X.addInputs("event_times", "labels");
        X.addOptionalInputs("amplitudes");
        X.addOutputs("firings_out");
        X.addRequiredParameters("central_channel");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.combine_firings", "0.1");
        X.addInputs("firings_list");
        X.addOutputs("firings_out");
        //X.addRequiredParameters();
        X.addOptionalParameter("increment_labels", "", "true");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.apply_timestamp_offset", "0.1");
        X.addInputs("firings");
        X.addOutputs("firings_out");
        X.addRequiredParameters("timestamp_offset");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.link_segments", "0.1");
        X.addInputs("firings", "firings_prev", "Kmax_prev");
        X.addOutputs("firings_out", "Kmax_out");
        X.addOutputs("firings_subset_out", "Kmax_out");
        X.addRequiredParameters("t1", "t2", "t1_prev", "t2_prev");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.cluster_metrics", "0.11");
        X.addInputs("timeseries", "firings");
        X.addOutputs("cluster_metrics_out");
        X.addRequiredParameters("samplerate");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.isolation_metrics", "0.15j");
        X.addInputs("timeseries", "firings");
        X.addOutputs("metrics_out");
        X.addOptionalOutputs("pair_metrics_out");
        X.addOptionalParameter("compute_bursting_parents", "", "false");
        processors.push_back(X.get_spec());
    }
    /*
    {
        ProcessorSpec X("ms3.extract_firings", "0.11");
        X.addInputs("firings");
        X.addOutputs("firings_out");
        X.addInputs("clusters");
        processors.push_back(X.get_spec());
    }
    */
    {
        ProcessorSpec X("ms3.combine_cluster_metrics", "0.1");
        X.addInputs("metrics_list");
        X.addOutputs("metrics_out");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.split_firings", "0.15");
        X.addInputs("timeseries_list", "firings");
        X.addOutputs("firings_out_list");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.mv_discrimhist", "0.1");
        X.addInputs("timeseries", "firings");
        X.addOutputs("output");
        X.addRequiredParameters("clusters");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.concat_firings", "0.13");
        X.addInputs("firings_list");
        X.addOptionalInputs("timeseries_list");
        X.addOutputs("firings_out");
        X.addOptionalOutputs("timeseries_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.concat_event_times", "0.1");
        X.addInputs("event_times_list");
        X.addOutputs("event_times_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.load_test", "0.1");
        X.addOutputs("stats_out");
        X.addOptionalParameters("num_read_bytes", "num_write_bytes", "num_cpu_ops");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.compute_amplitudes", "0.1");
        X.addInputs("timeseries", "event_times");
        X.addOutputs("amplitudes_out");
        X.addRequiredParameters("central_channel");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.confusion_matrix", "0.17");
        X.addInputs("firings1", "firings2");
        X.addOutputs("confusion_matrix_out");
        X.addOptionalOutputs("matched_firings_out", "label_map_out", "firings2_relabeled_out", "firings2_relabel_map_out");
        X.addOptionalParameter("max_matching_offset", "", 30);
        X.addOptionalParameter("relabel_firings2", "", "false");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.mask_out_artifacts", "0.1");
        X.addInputs("timeseries");
        X.addOutputs("timeseries_out");
        X.addOptionalParameter("threshold","",6);
        X.addOptionalParameter("interval_size","",2000);
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.mv_compute_templates", "0.1");
        X.addInputs("timeseries","firings");
        X.addOutputs("templates_out","stdevs_out");
        X.addOptionalParameter("clip_size","",200);
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.mv_subfirings", "0.1");
        X.addInputs("firings");
        X.addOutputs("firings_out");
        X.addRequiredParameter("labels");
        X.addOptionalParameter("max_per_label","",0);
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("ms3.mv_compute_amplitudes", "0.1");
        X.addInputs("timeseries", "firings");
        X.addOutputs("firings_out");
        processors.push_back(X.get_spec());
    }


    QJsonObject ret;
    ret["processors"] = processors;
    return ret;
}

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);

    CLParams CLP(argc, argv);

    QString arg1 = CLP.unnamed_parameters.value(0);
    QString arg2 = CLP.unnamed_parameters.value(1);

    QString pname;
    bool requirements_only = false;
    if (arg1 == "spec") {
        QJsonObject spec = get_spec();
        if (arg2.isEmpty()) {
            QString json = QJsonDocument(spec).toJson(QJsonDocument::Indented);
            printf("%s\n", json.toUtf8().data());
            return 0;
        }
        else {
            pname = arg2;
            QJsonArray P = spec["processors"].toArray();
            bool found=false;
            QJsonObject spec0;
            for (int i=0; i < P.count(); i++) {
                QJsonObject spec1 = P[i].toObject();
                if (spec1["name"] == pname) {
                    spec0 = spec1;
                    found=true;
                }
            }
            if (!found) {
                printf("{}\n");
                return -1;
            }
            QString json = QJsonDocument(spec0).toJson(QJsonDocument::Indented);
            printf("%s\n", json.toUtf8().data());
            return 0;
        }
    }
    else if (arg1 == "requirements") {
        pname = arg2;
        requirements_only = true;
    }
    else {
        pname = arg1;
    }

    bool ret = false;

    if (CLP.named_parameters.contains("_request_num_threads")) {
        int num_threads = CLP.named_parameters.value("_request_num_threads", 0).toInt();
        if (num_threads) {
            qDebug().noquote() << "Setting num threads:" << num_threads;
            omp_set_num_threads(num_threads);
        }
    }

    if (pname == "ms3.run_metrics_script") {
        QString metrics = CLP.named_parameters["metrics"].toString();
        QString script = CLP.named_parameters["script"].toString();
        QString metrics_out = CLP.named_parameters["metrics_out"].toString();
        ret = p_run_metrics_script(metrics, script, metrics_out);
    }
    else if (pname == "spikeview.metrics1") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString metrics_out = CLP.named_parameters["metrics_out"].toString();
        P_spikeview_metrics1_opts opts;
        opts.samplerate = CLP.named_parameters["samplerate"].toDouble();
        ret = p_spikeview_metrics1(firings, metrics_out, opts);
    }
#ifndef NO_FFTW3
    else if (pname == "spikeview.templates") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString templates_out = CLP.named_parameters["templates_out"].toString();
        P_spikeview_templates_opts opts;
        opts.clip_size = CLP.named_parameters["clip_size"].toInt();
        opts.max_events_per_template = CLP.named_parameters["max_events_per_template"].toDouble();
        opts.samplerate = CLP.named_parameters.value("samplerate", 0).toDouble();
        opts.freq_min = CLP.named_parameters.value("freq_min", 0).toDouble();
        opts.freq_max = CLP.named_parameters.value("freq_max", 0).toDouble();
        if ((opts.samplerate != 0) && ((opts.freq_min != 0) || (opts.freq_max != 0)))
            opts.filt_padding = 50;
        opts.subtract_temporal_mean = (CLP.named_parameters.value("subtract_temporal_mean") == "true");
        ret = p_spikeview_templates(timeseries, firings, templates_out, opts);
    }
#endif
    else if (pname == "banjoview.cross_correlograms") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString correlograms_out = CLP.named_parameters["correlograms_out"].toString();
        P_banjoview_cross_correlograms_opts opts;
        opts.samplerate = CLP.named_parameters["samplerate"].toDouble();
        opts.max_dt_msec = CLP.named_parameters["max_dt_msec"].toDouble();
        opts.bin_size_msec = CLP.named_parameters["bin_size_msec"].toDouble();
        QString mode = CLP.named_parameters.value("mode", "autocorrelograms").toString();
        if (mode == "autocorrelograms")
            opts.mode = Autocorrelograms;
        else if (mode == "matrix_of_cross_correlograms")
            opts.mode = Matrix_of_cross_correlograms;
        else {
            qWarning() << "Unexpected mode: " + mode;
            return -1;
        }
        ret = p_banjoview_cross_correlograms(firings, correlograms_out, opts);
    }
    /*
    else if (pname == "ms3.synthesize_timeseries") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString waveforms = CLP.named_parameters["waveforms"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        P_synthesize_timeseries_opts opts;
        opts.noise_level = CLP.named_parameters["noise_level"].toDouble();
        opts.duration = CLP.named_parameters["duration"].toDouble();
        opts.waveform_upsample_factor = CLP.named_parameters.value("waveform_upsample_factor", 13).toDouble();
        ret = p_synthesize_timeseries(firings, waveforms, timeseries_out, opts);
    }
    */
    else if (pname == "ms3.combine_firing_segments") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QStringList firings_list = MLUtil::toStringList(CLP.named_parameters["firings_list"]);
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        P_combine_firing_segments_opts opts;
        if (CLP.named_parameters.contains("clip_size"))
            opts.clip_size = CLP.named_parameters["clip_size"].toInt();
        if (CLP.named_parameters.contains("match_score_threshold"))
            opts.match_score_threshold = CLP.named_parameters["match_score_threshold"].toDouble();
        if (CLP.named_parameters.contains("offset_search_radius"))
            opts.offset_search_radius = CLP.named_parameters["offset_search_radius"].toDouble();
        ret = p_combine_firing_segments(timeseries, firings_list, firings_out, opts);
    }
    else if (pname == "ms3.extract_firings") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString metrics = CLP.named_parameters["metrics"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        P_extract_firings_opts opts;
        opts.exclusion_tags = MLUtil::toStringList(CLP.named_parameters["exclusion_tags"]);
        {
            QStringList clusters_str = CLP.named_parameters["clusters"].toString().split(",");
            opts.clusters = MLUtil::stringListToIntList(clusters_str);
        }
        {
            QString t1_str = CLP.named_parameters["t1"].toString();
            if (t1_str.isEmpty())
                t1_str = "-1";
            QString t2_str = CLP.named_parameters["t2"].toString();
            if (t2_str.isEmpty())
                t2_str = "-1";
            opts.t1 = t1_str.toDouble();
            opts.t2 = t2_str.toDouble();
        }
        ret = p_extract_firings(firings, metrics, firings_out, opts);
    }
    else if (pname == "ms3.concat_timeseries") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        ret = p_concat_timeseries(timeseries_list, timeseries_out);
    }
    else if (pname == "ms3.create_multiscale_timeseries") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        QString tempdir = CLP.named_parameters["_tempdir"].toString();
        if (tempdir.isEmpty()) {
            qWarning() << "Required _tempdir parameter not provided";
            return -1;
        }
        ret = p_create_multiscale_timeseries(timeseries,timeseries_out,tempdir);
    }
#ifndef NO_FFTW3
    else if (pname == "ms3.bandpass_filter") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        Bandpass_filter_opts opts;
        opts.samplerate = CLP.named_parameters["samplerate"].toDouble();
        opts.freq_min = CLP.named_parameters["freq_min"].toDouble();
        opts.freq_max = CLP.named_parameters["freq_max"].toDouble();
        opts.freq_wid = CLP.named_parameters.value("freq_wid", 1000).toDouble();
        opts.quantization_unit = CLP.named_parameters.value("quantization_unit").toDouble();
        opts.subsample_factor = CLP.named_parameters.value("subsample_factor", 1).toInt();
        if (requirements_only) opts.requirements_only = true;
        ret = p_bandpass_filter(timeseries, timeseries_out, opts);
        if (requirements_only) {
            if (ret) {
                QJsonObject requirements0;
                requirements0["peak_ram_mb"] = opts.expected_peak_ram_mb;
                QString json = QJsonDocument(requirements0).toJson(QJsonDocument::Indented);
                printf("%s\n", json.toUtf8().data());
            }
            else {
                printf("Error in processor while trying to retrieve requirements.\n");
            }
        }
    }
#endif
    else if (pname == "ms3.whiten") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        Whiten_opts opts;
        opts.quantization_unit = CLP.named_parameters["quantization_unit"].toDouble();
        if (requirements_only)
            opts.requirements_only = true;
        ret = p_whiten(timeseries, timeseries_out, opts);
        if (requirements_only) {
            if (ret) {
                QJsonObject requirements0;
                requirements0["peak_ram_mb"] = opts.expected_peak_ram_mb;
                QString json = QJsonDocument(requirements0).toJson(QJsonDocument::Indented);
                printf("%s\n", json.toUtf8().data());
            }
            else {
                printf("Error in processor while trying to retrieve requirements.\n");
            }
        }
    }
    else if (pname == "ms3.compute_whitening_matrix") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QString whitening_matrix_out = CLP.named_parameters["whitening_matrix_out"].toString();
        QStringList channels_str = CLP.named_parameters["channels"].toString().split(",", QString::SkipEmptyParts);
        QList<int> channels = MLUtil::stringListToIntList(channels_str);
        Whiten_opts opts;
        ret = p_compute_whitening_matrix(timeseries_list, channels, whitening_matrix_out, opts);
    }
    else if (pname == "ms3.whiten_clips") {
        QString clips = CLP.named_parameters["clips"].toString();
        QString whitening_matrix = CLP.named_parameters["whitening_matrix"].toString();
        QString clips_out = CLP.named_parameters["clips_out"].toString();
        Whiten_opts opts;
        opts.quantization_unit = CLP.named_parameters["quantization_unit"].toDouble();
        ret = p_whiten_clips(clips, whitening_matrix, clips_out, opts);
    }
    else if (pname == "ms3.apply_whitening_matrix") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString whitening_matrix = CLP.named_parameters["whitening_matrix"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        Whiten_opts opts;
        opts.quantization_unit = CLP.named_parameters["quantization_unit"].toDouble();
        ret = p_apply_whitening_matrix(timeseries, whitening_matrix, timeseries_out, opts);
    }
    else if (pname == "ms3.extract_clips") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries"]);
        QString event_times = CLP.named_parameters["event_times"].toString();
        QString clips_out = CLP.named_parameters["clips_out"].toString();
        QStringList channels_str = CLP.named_parameters["channels"].toString().split(",", QString::SkipEmptyParts);
        QList<int> channels = MLUtil::stringListToIntList(channels_str);
        ret = p_extract_clips(timeseries_list, event_times, channels, clips_out, CLP.named_parameters);
    }
    else if (pname == "ms3.mv_extract_clips") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries"]);
        QString firings = CLP.named_parameters["firings"].toString();
        QString clips_out = CLP.named_parameters["clips_out"].toString();
        QStringList channels_str = CLP.named_parameters["channels"].toString().split(",", QString::SkipEmptyParts);
        QList<int> channels = MLUtil::stringListToIntList(channels_str);
        ret = p_mv_extract_clips(timeseries_list, firings, channels, clips_out, CLP.named_parameters);
    }
    else if (pname == "ms3.mv_extract_clips_features") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString features_out = CLP.named_parameters["features_out"].toString();
        int clip_size = CLP.named_parameters["clip_size"].toInt();
        int num_features = CLP.named_parameters["num_features"].toInt();
        int subtract_mean = CLP.named_parameters["subtract_mean"].toInt();
        ret = p_mv_extract_clips_features(timeseries, firings, features_out, clip_size, num_features, subtract_mean);
    }
    else if (pname == "ms3.compute_templates") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries"]);
        QString firings = CLP.named_parameters["firings"].toString();
        QString templates_out = CLP.named_parameters["templates_out"].toString();
        int clip_size = CLP.named_parameters["clip_size"].toInt();
        QList<int> clusters = MLUtil::stringListToIntList(CLP.named_parameters["clusters"].toString().split(",", QString::SkipEmptyParts));
        ret = p_compute_templates(timeseries_list, firings, templates_out, clip_size, clusters);
    }
    else if (pname == "ms3.reorder_labels") {
        QString templates = CLP.named_parameters["templates"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        ret = p_reorder_labels(templates, firings, firings_out);
    }
    else if (pname == "ms3.create_firings") {
        QString event_times = CLP.named_parameters["event_times"].toString();
        QString labels = CLP.named_parameters["labels"].toString();
        QString amplitudes = CLP.named_parameters["amplitudes"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        int central_channel = CLP.named_parameters["central_channel"].toInt();
        ret = p_create_firings(event_times, labels, amplitudes, firings_out, central_channel);
    }
    else if (pname == "ms3.combine_firings") {
        QStringList firings_list = MLUtil::toStringList(CLP.named_parameters["firings_list"]);
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        QString tmp = CLP.named_parameters.value("increment_labels", "").toString();
        bool increment_labels = (tmp == "true");
        ret = p_combine_firings(firings_list, firings_out, increment_labels);
    }
    else if (pname == "ms3.apply_timestamp_offset") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        double timestamp_offset = CLP.named_parameters["timestamp_offset"].toDouble();
        ret = p_apply_timestamp_offset(firings, firings_out, timestamp_offset);
    }
    else if (pname == "ms3.link_segments") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_prev = CLP.named_parameters["firings_prev"].toString();
        QString Kmax_prev = CLP.named_parameters["Kmax_prev"].toString();

        QString firings_out = CLP.named_parameters["firings_out"].toString();
        QString firings_subset_out = CLP.named_parameters["firings_subset_out"].toString();
        QString Kmax_out = CLP.named_parameters["Kmax_out"].toString();

        double t1 = CLP.named_parameters["t1"].toDouble();
        double t2 = CLP.named_parameters["t2"].toDouble();
        double t1_prev = CLP.named_parameters["t1_prev"].toDouble();
        double t2_prev = CLP.named_parameters["t2_prev"].toDouble();
        ret = p_link_segments(firings, firings_prev, Kmax_prev, firings_out, firings_subset_out, Kmax_out, t1, t2, t1_prev, t2_prev);
    }
    else if (pname == "ms3.cluster_metrics") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString cluster_metrics_out = CLP.named_parameters["cluster_metrics_out"].toString();
        Cluster_metrics_opts opts;
        opts.samplerate = CLP.named_parameters["samplerate"].toDouble();
        ret = p_cluster_metrics(timeseries, firings, cluster_metrics_out, opts);
    }
    else if (pname == "ms3.isolation_metrics") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries"]);
        QString firings = CLP.named_parameters["firings"].toString();
        QString metrics_out = CLP.named_parameters["metrics_out"].toString();
        QString pair_metrics_out = CLP.named_parameters["pair_metrics_out"].toString();
        P_isolation_metrics_opts opts;
        opts.compute_bursting_parents = (CLP.named_parameters["compute_bursting_parents"].toString() == "true");
        ret = p_isolation_metrics(timeseries_list, firings, metrics_out, pair_metrics_out, opts);
    }
    else if (pname == "ms3.combine_cluster_metrics") {
        QStringList metrics_list = MLUtil::toStringList(CLP.named_parameters["metrics_list"]);
        QString metrics_out = CLP.named_parameters["metrics_out"].toString();
        ret = p_combine_cluster_metrics(metrics_list, metrics_out);
    }
    else if (pname == "ms3.split_firings") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QString firings = CLP.named_parameters["firings"].toString();
        QStringList firings_out_list = MLUtil::toStringList(CLP.named_parameters["firings_out_list"]);
        ret = p_split_firings(timeseries_list, firings, firings_out_list);
    }
    else if (pname == "ms3.mv_discrimhist") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString output = CLP.named_parameters["output"].toString();
        mv_discrimhist_opts opts;
        {
            QStringList clusters_str = CLP.named_parameters["clusters"].toString().split(",");
            opts.clusters = MLUtil::stringListToIntList(clusters_str);
        }
        ret = mv_discrimhist(timeseries, firings, output, opts);
    }
    /*
    else if (pname == "ms3.extract_firings") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        QStringList clusters_str = MLUtil::toStringList(CLP.named_parameters["clusters"]);
        QSet<int> clusters = MLUtil::stringListToIntList(clusters_str).toSet();
        ret = p_extract_firings(firings, clusters, firings_out);
    }
    */
    else if (pname == "ms3.concat_firings") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QStringList firings_list = MLUtil::toStringList(CLP.named_parameters["firings_list"]);
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        ret = p_concat_firings(timeseries_list, firings_list, timeseries_out, firings_out);
    }
    else if (pname == "ms3.concat_event_times") {
        QStringList event_times_list = MLUtil::toStringList(CLP.named_parameters["event_times_list"]);
        QString event_times_out = CLP.named_parameters["event_times_out"].toString();
        ret = p_concat_event_times(event_times_list, event_times_out);
    }
    else if (pname == "ms3.load_test") {
        QString stats_out = CLP.named_parameters["stats_out"].toString();
        P_load_test_opts opts;
        opts.num_cpu_ops = CLP.named_parameters["num_cpu_ops"].toDouble();
        opts.num_read_bytes = CLP.named_parameters["num_read_bytes"].toDouble();
        opts.num_write_bytes = CLP.named_parameters["num_write_bytes"].toDouble();
        ret = p_load_test(stats_out, opts);
    }
    else if (pname == "ms3.compute_amplitudes") {
        P_compute_amplitudes_opts opts;
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString event_times = CLP.named_parameters["event_times"].toString();
        opts.central_channel = CLP.named_parameters["central_channel"].toInt();
        QString amplitudes_out = CLP.named_parameters["amplitudes_out"].toString();
        ret = p_compute_amplitudes(timeseries, event_times, amplitudes_out, opts);
    }
    else if (pname == "ms3.confusion_matrix") {
        P_confusion_matrix_opts opts;
        QString firings1 = CLP.named_parameters["firings1"].toString();
        QString firings2 = CLP.named_parameters["firings2"].toString();
        QString confusion_matrix_out = CLP.named_parameters["confusion_matrix_out"].toString();
        QString matched_firings_out = CLP.named_parameters.value("matched_firings_out").toString();
        QString label_map_out = CLP.named_parameters.value("label_map_out").toString();
        QString firings2_relabeled_out = CLP.named_parameters.value("firings2_relabeled_out").toString();
        QString firings2_relabel_map_out = CLP.named_parameters.value("firings2_relabel_map_out").toString();
        if (CLP.named_parameters.contains("max_matching_offset")) {
            opts.max_matching_offset = CLP.named_parameters.value("max_matching_offset").toInt();
        }
        opts.relabel_firings2 = (CLP.named_parameters.value("relabel_firings2", "false").toString() == "true");
        ret = p_confusion_matrix(firings1, firings2, confusion_matrix_out, matched_firings_out, label_map_out, firings2_relabeled_out, firings2_relabel_map_out, opts);
    }
    else if (pname == "ms3.mask_out_artifacts") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        double threshold = CLP.named_parameters["threshold"].toDouble();
        bigint interval_size = CLP.named_parameters["interval_size"].toDouble();
        ret = p_mask_out_artifacts(timeseries,timeseries_out,threshold,interval_size);
    }
    else if (pname == "ms3.mv_compute_templates") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString templates_out = CLP.named_parameters["templates_out"].toString();
        QString stdevs_out = CLP.named_parameters["stdevs_out"].toString();
        int clip_size = CLP.named_parameters["clip_size"].toInt();
        ret = mv_compute_templates(timeseries,firings,templates_out,stdevs_out,clip_size);
    }
    else if (pname == "ms3.mv_subfirings") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        QStringList labels_str = MLUtil::toStringList(CLP.named_parameters["labels"]);
        QList<int> labels = MLUtil::stringListToIntList(labels_str);
        bigint max_per_label = CLP.named_parameters["max_per_label"].toDouble();
        ret = mv_subfirings(firings,firings_out,labels.toVector(),max_per_label);
    }
    else if (pname == "ms3.mv_compute_amplitudes") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        p_mv_compute_amplitudes_opts opts;
        ret = p_mv_compute_amplitudes(timeseries, firings, firings_out, opts);
    }
    else {
        qWarning() << "Unexpected processor name: " + pname;
        return -1;
    }

    if (!ret)
        return -1;

    return 0;
}

QJsonObject ProcessorSpecFile::get_spec()
{
    QJsonObject ret;
    ret["name"] = name;
    ret["description"] = description;
    ret["optional"] = optional;
    return ret;
}

QJsonObject ProcessorSpecParam::get_spec()
{
    QJsonObject ret;
    ret["name"] = name;
    ret["description"] = description;
    ret["optional"] = optional;
    if (default_value.isValid())
        ret["default_value"] = QJsonValue::fromVariant(default_value);
    return ret;
}

void ProcessorSpec::addInputs(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addInput(name1);
    if (!name2.isEmpty())
        addInput(name2);
    if (!name3.isEmpty())
        addInput(name3);
    if (!name4.isEmpty())
        addInput(name4);
    if (!name5.isEmpty())
        addInput(name5);
}

void ProcessorSpec::addOptionalInputs(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addInput(name1, "", true);
    if (!name2.isEmpty())
        addInput(name2, "", true);
    if (!name3.isEmpty())
        addInput(name3, "", true);
    if (!name4.isEmpty())
        addInput(name4, "", true);
    if (!name5.isEmpty())
        addInput(name5, "", true);
}

void ProcessorSpec::addOutputs(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addOutput(name1);
    if (!name2.isEmpty())
        addOutput(name2);
    if (!name3.isEmpty())
        addOutput(name3);
    if (!name4.isEmpty())
        addOutput(name4);
    if (!name5.isEmpty())
        addOutput(name5);
}

void ProcessorSpec::addOptionalOutputs(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addOutput(name1, "", true);
    if (!name2.isEmpty())
        addOutput(name2, "", true);
    if (!name3.isEmpty())
        addOutput(name3, "", true);
    if (!name4.isEmpty())
        addOutput(name4, "", true);
    if (!name5.isEmpty())
        addOutput(name5, "", true);
}

void ProcessorSpec::addRequiredParameters(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addRequiredParameter(name1);
    if (!name2.isEmpty())
        addRequiredParameter(name2);
    if (!name3.isEmpty())
        addRequiredParameter(name3);
    if (!name4.isEmpty())
        addRequiredParameter(name4);
    if (!name5.isEmpty())
        addRequiredParameter(name5);
}

void ProcessorSpec::addOptionalParameters(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addOptionalParameter(name1);
    if (!name2.isEmpty())
        addOptionalParameter(name2);
    if (!name3.isEmpty())
        addOptionalParameter(name3);
    if (!name4.isEmpty())
        addOptionalParameter(name4);
    if (!name5.isEmpty())
        addOptionalParameter(name5);
}

void ProcessorSpec::addInput(QString name, QString description, bool optional)
{
    ProcessorSpecFile X;
    X.name = name;
    X.description = description;
    X.optional = optional;
    inputs.append(X);
}

void ProcessorSpec::addOutput(QString name, QString description, bool optional)
{
    ProcessorSpecFile X;
    X.name = name;
    X.description = description;
    X.optional = optional;
    outputs.append(X);
}

void ProcessorSpec::addRequiredParameter(QString name, QString description)
{
    ProcessorSpecParam X;
    X.name = name;
    X.description = description;
    X.optional = false;
    parameters.append(X);
}

void ProcessorSpec::addOptionalParameter(QString name, QString description, QVariant default_value)
{
    ProcessorSpecParam X;
    X.name = name;
    X.description = description;
    X.optional = true;
    X.default_value = default_value;
    parameters.append(X);
}

QJsonObject ProcessorSpec::get_spec()
{
    QJsonObject ret;
    ret["name"] = processor_name;
    ret["version"] = version;
    ret["description"] = description;
    QJsonArray inputs0;
    for (int i = 0; i < inputs.count(); i++) {
        inputs0.push_back(inputs[i].get_spec());
    }
    ret["inputs"] = inputs0;
    QJsonArray outputs0;
    for (int i = 0; i < outputs.count(); i++) {
        outputs0.push_back(outputs[i].get_spec());
    }
    ret["outputs"] = outputs0;
    QJsonArray parameters0;
    for (int i = 0; i < parameters.count(); i++) {
        parameters0.push_back(parameters[i].get_spec());
    }
    ret["parameters"] = parameters0;
    ret["exe_command"] = qApp->applicationFilePath() + " " + processor_name + " $(arguments)";
    if (this->can_return_requirements)
        ret["requirements_command"] = qApp->applicationFilePath() + " requirements " + processor_name + " $(arguments)";
    return ret;
}
