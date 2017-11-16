Processing pipelines
====================

A processing pipeline is a collection of processing steps that are related by their inputs and outputs. In spike sorting, for example, a typical pipeline would comprise

* Timeseries data extraction
* Bandpass filter
* Whitening
* Spike sorting
* Metric derivation
* Automatic cluster annotation
* Extraction of final spike sorting output (firings)
* Computation of summary data for visualization

Each step of the pipeline can either be a single low level processor (see :doc:`processing system`) or can be another pipeline (sub-pipeline). Just like processors, each pipeline contains a specification object describing the inputs, outputs, and parameters. As with processors, inputs and outputs always represent data files (or lists of a data files), whereas parameters are strings, integers, or lists thereof.

Pipelines versus processors
---------------------------

While pipelines and processors are similar in that they both operate on files, the difference is that processors are installed and registered directly on the processing server, whereas pipelines can simply be represented by a JSON object. Therefore, a pipeline may run on any computer, or even on a web browser, whereas the processors that get executed must always run on a specific computer that has direct access to the underlying data. Of course, pipelines are made up of processors, so the individual processing steps will ultimately need to run on the computer with the data. Therefore, running a pipeline job involves specifying a server where the individual processor steps will get executed.

Formal definition of a pipeline
-------------------------------

A pipeline consists of (a) a specification object and (b) a list of pipeline steps.

First, the pipeline specification (spec) is a JSON object as in the following example:

.. code :: JSON

	{
		"name": "sort_ms3_v1",
		"description": "",
    	"inputs": [
      		{"name": "raw"},
      		{"name": "geom"},
      		{"name": "annotation_script"}
    	],
    	"outputs": [
      		{"name": "pre"},
      		{"name": "filt"},
      		{"name": "firings"},
      		{"name": "firings1"}
    	],
    	"parameters": [
      		{"name": "detect_sign"},
      		{"name": "detect_threshold"},
      		{"name": "adjacency_radius"},
      		{"name": "samplerate"}
    	]
	}

Thus, the pipeline has a name, an optional description, a collection of named inputs, outputs, and parameters. Note that this is very similar to the spec object of processors.

The second component of a pipeline is a list of processing steps. In our example we have

.. code :: JSON

  [
    {
      "processor_name": "mountainsort.bandpass_filter",
      "inputs": {
        "timeseries": "raw"
      },
      "outputs": {
        "timeseries_out": "filt"
      },
      "parameters": {
        "samplerate": "${samplerate}",
        "freq_min": "300",
        "freq_max": "6000",
        "freq_wid": "1000",
        "quantization_unit": "",
        "subsample_factor": ""
      },
      "step_type": "processor"
    },
    {
      "processor_name": "mountainsort.whiten",
      "inputs": {
        "timeseries": "filt"
      },
      "outputs": {
        "timeseries_out": "pre"
      },
      "parameters": {
        "quantization_unit": ""
      },
      "step_type": "processor"
    },
    {
      "processor_name": "mountainsort.mountainsort3",
      "inputs": {
        "timeseries": "pre",
        "geom": "geom"
      },
      "outputs": {
        "firings_out": "firings1"
      },
      "parameters": {
        "adjacency_radius": "${adjacency_radius}",
        "consolidate_clusters": "",
        "consolidation_factor": "",
        "clip_size": "",
        "detect_interval": "",
        "detect_threshold": "",
        "detect_sign": "${detect_sign}",
        "merge_across_channels": "",
        "fit_stage": "",
        "t1": "",
        "t2": ""
      },
      "step_type": "processor"
    },
    {
      "processor_name": "mountainsort.cluster_metrics",
      "inputs": {
        "timeseries": "pre",
        "firings": "firings1"
      },
      "outputs": {
        "cluster_metrics_out": "metrics1"
      },
      "parameters": {
        "samplerate": "${samplerate}"
      },
      "step_type": "processor"
    },
    {
      "processor_name": "mountainsort.isolation_metrics",
      "inputs": {
        "timeseries": "pre",
        "firings": "firings1"
      },
      "outputs": {
        "metrics_out": "metrics2",
        "pair_metrics_out": ""
      },
      "parameters": {
        "compute_bursting_parents": "true"
      },
      "step_type": "processor"
    },
    {
      "processor_name": "mountainsort.combine_cluster_metrics",
      "inputs": {
        "metrics_list": "metrics1,metrics2"
      },
      "outputs": {
        "metrics_out": "metrics3"
      },
      "parameters": {},
      "step_type": "processor"
    },
    {
      "processor_name": "mountainsort.run_metrics_script",
      "inputs": {
        "metrics": "metrics3",
        "script": "annotation_script"
      },
      "outputs": {
        "metrics_out": "metrics_annotated"
      },
      "parameters": {},
      "step_type": "processor"
    },
    {
      "processor_name": "mountainsort.extract_firings",
      "inputs": {
        "firings": "firings1",
        "metrics": "metrics_annotated"
      },
      "outputs": {
        "firings_out": "firings"
      },
      "parameters": {
        "exclusion_tags": "rejected",
        "clusters": "",
        "t1": "",
        "t2": ""
      },
      "step_type": "processor"
    }
  ]