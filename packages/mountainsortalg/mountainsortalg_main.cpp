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

#include "mountainsortalg_main.h"
#include "omp.h"
#include "p_mountainsort3.h"

#include <QCoreApplication>
#include <QJsonArray>
#include <QJsonDocument>
#include <clparams.h>
#include <QDebug>
#include <QDir>

QJsonObject get_spec()
{
    QJsonArray processors;
    {
        ProcessorSpec X("mountainsortalg.ms3alg", "0.12");
        X.addInput("timeseries","Preprocessed timeseries");
        X.addOptionalInputs("geom");
        X.addOutput("firings_out","Firings array (times/labels)");
        X.addOptionalParameter("adjacency_radius", "Radius of local sorting neighborhood, corresponding to the geometry file. 0 means each channel is sorted independently. -1 means all channels are included in every neighborhood.", -1);
        X.addOptionalParameter("consolidate_clusters", "", "true");
        X.addOptionalParameter("consolidation_factor", "", 0.9);
        X.addOptionalParameter("clip_size", "", 50);
        X.addOptionalParameter("detect_interval", "", 10);
        X.addOptionalParameter("detect_threshold", "", 3);
        X.addOptionalParameter("detect_sign", "", 0);
        X.addOptionalParameter("merge_across_channels", "", "true");
        X.addOptionalParameter("fit_stage", "", "true");
        X.addOptionalParameter("t1", "Start timepoint to do the sorting (default 0 means start at beginning)", 0);
        X.addOptionalParameter("t2", "End timepoint for sorting (default -1 means go to end of timeseries)", -1);
        QJsonObject spec0=X.get_spec();
        processors.push_back(spec0);
        spec0["name"]="mountainsortalg.ms3"; //backward compatibility
        spec0["description"]="WARNING! Please use mountainsortalg.ms3alg rather than mountainsortalg.ms3. Right now they are equivalent, but the latter will be removed in the future.";
        QString str=spec0["exe_command"].toString();
        str=str.replace(".ms3alg",".ms3");
        spec0["exe_command"]=str;
        processors.push_back(spec0);
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

    if ((arg1 == "mountainsortalg.ms3alg")||(arg1 == "mountainsortalg.ms3")) {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString geom = CLP.named_parameters["geom"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        P_mountainsort3_opts opts;
        opts.adjacency_radius = CLP.named_parameters.value("adjacency_radius").toDouble();
        opts.consolidate_clusters = (CLP.named_parameters.value("consolidate_clusters").toString() == "true");
        opts.consolidation_factor = CLP.named_parameters.value("consolidation_factor").toDouble();
        opts.clip_size = CLP.named_parameters.value("clip_size").toDouble();
        opts.detect_interval = CLP.named_parameters.value("detect_interval").toDouble();
        opts.detect_threshold = CLP.named_parameters.value("detect_threshold").toDouble();
        opts.detect_sign = CLP.named_parameters.value("detect_sign").toInt();
        opts.merge_across_channels = (CLP.named_parameters.value("merge_across_channels").toString() == "true");
        opts.fit_stage = (CLP.named_parameters.value("fit_stage").toString() == "true");
        opts.t1 = CLP.named_parameters.value("t1","-1").toDouble();
        opts.t2 = CLP.named_parameters.value("t2","-1").toDouble();
        QString temp_path = CLP.named_parameters.value("_tempdir").toString();
        if (temp_path.isEmpty()) {
            temp_path = QDir::currentPath();
        }
        ret = p_mountainsort3(timeseries, geom, firings_out, temp_path, opts);

        if (arg1 == "mountainsortalg.ms3") {
            qWarning() << "********************************************************************************************************************************";
            qWarning() << "WARNING! Please use mountainsortalg.ms3alg rather than mountainsortalg.ms3. Right now they are equivalent, but the latter will be removed in the future.";
            qWarning() << "********************************************************************************************************************************";
            qWarning() << "";
            qWarning() << "";
        }
    }
    /*
    else if (arg1 == "mountainsort.run_metrics_script") {
        QString metrics = CLP.named_parameters["metrics"].toString();
        QString script = CLP.named_parameters["script"].toString();
        QString metrics_out = CLP.named_parameters["metrics_out"].toString();
        ret = p_run_metrics_script(metrics, script, metrics_out);
    }
    */
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
