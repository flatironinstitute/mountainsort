/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 2/22/2017
*******************************************************/

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
#include "p_synthesize_timeseries.h"
#include "p_combine_firing_segments.h"
#include "p_extract_firings.h"
#include "p_concat_timeseries.h"
#include "p_banjoview_cross_correlograms.h"
#include "p_create_multiscale_timeseries.h"

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
    {
        ProcessorSpec X("ms3.synthesize_timeseries", "0.12");
        X.addInputs("firings", "waveforms");
        X.addOutputs("timeseries_out");
        X.addOptionalParameter("noise_level", "", 0);
        X.addOptionalParameter("duration", "", 0);
        X.addOptionalParameter("waveform_upsample_factor", "", 13);
        processors.push_back(X.get_spec());
    }
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
        ProcessorSpec X("ms3.extract_firings", "0.1");
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

    QJsonObject ret;
    ret["processors"] = processors;
    return ret;
}

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);

    CLParams CLP(argc, argv);

    QString arg1 = CLP.unnamed_parameters.value(0);

    if (arg1 == "spec") {
        QJsonObject spec = get_spec();
        QString json = QJsonDocument(spec).toJson(QJsonDocument::Indented);
        printf("%s\n", json.toUtf8().data());
        return 0;
    }

    bool ret = false;

    if (CLP.named_parameters.contains("_request_num_threads")) {
        int num_threads = CLP.named_parameters.value("_request_num_threads", 0).toInt();
        if (num_threads) {
            qDebug().noquote() << "Setting num threads:" << num_threads;
            omp_set_num_threads(num_threads);
        }
    }

    if (arg1 == "ms3.run_metrics_script") {
        QString metrics = CLP.named_parameters["metrics"].toString();
        QString script = CLP.named_parameters["script"].toString();
        QString metrics_out = CLP.named_parameters["metrics_out"].toString();
        ret = p_run_metrics_script(metrics, script, metrics_out);
    }
    else if (arg1 == "spikeview.metrics1") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString metrics_out = CLP.named_parameters["metrics_out"].toString();
        P_spikeview_metrics1_opts opts;
        opts.samplerate = CLP.named_parameters["samplerate"].toDouble();
        ret = p_spikeview_metrics1(firings, metrics_out, opts);
    }
#ifndef NO_FFTW3
    else if (arg1 == "spikeview.templates") {
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
    else if (arg1 == "banjoview.cross_correlograms") {
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
    else if (arg1 == "ms3.synthesize_timeseries") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString waveforms = CLP.named_parameters["waveforms"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        P_synthesize_timeseries_opts opts;
        opts.noise_level = CLP.named_parameters["noise_level"].toDouble();
        opts.duration = CLP.named_parameters["duration"].toDouble();
        opts.waveform_upsample_factor = CLP.named_parameters.value("waveform_upsample_factor", 13).toDouble();
        ret = p_synthesize_timeseries(firings, waveforms, timeseries_out, opts);
    }
    else if (arg1 == "ms3.combine_firing_segments") {
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
    else if (arg1 == "ms3.extract_firings") {
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
    else if (arg1 == "ms3.concat_timeseries") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        ret = p_concat_timeseries(timeseries_list, timeseries_out);
    }
    else if (arg1 == "ms3.create_multiscale_timeseries") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        QString tempdir = CLP.named_parameters["_tempdir"].toString();
        if (tempdir.isEmpty()) {
            qWarning() << "Required _tempdir parameter not provided";
            return -1;
        }
        ret = p_create_multiscale_timeseries(timeseries,timeseries_out,tempdir);
    }
    else {
        qWarning() << "Unexpected processor name: " + arg1;
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
    return ret;
}
