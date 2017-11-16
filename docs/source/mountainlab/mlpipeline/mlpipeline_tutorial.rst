MLPipeline Tutorial
===================

MLPipeline is a tool for creating, running, and sharing analysis pipelines for use in MountainLab. It can either be used on your local computer (desktop) or via web browser. These two methods each have advantages. The goal of MLPipeline is to provide a consistent interface so that the same processing pipelines can be used on the local machine and on the web.

Running on the desktop
----------------------

To start the MLPipeline desktop pipeline, first make sure the software is installed properly [link needed]. Then type the following from a terminal

.. code:: bash

	mlpipeline

If the larinet server is not already running on your machine, you will get instructions on starting that service which provides access to the mountainlab processing commands. Elsewhere you will see how larinet (along with a second server called kulelepoller) can be used to enable (password-protected) remote access to your computer (or processing server) via web browser from anywhere on the internet.

Running from a web browser
--------------------------

Alternatively you can use MLPipeline from any web browser connected to the internet. Navigate to

.. code:: bash

	https://mlp.herokuapp.com

Constructing a processing pipeline
----------------------------------

A MLPipeline document comprises a collection of pipelines with one special pipeline called "main". By default you start with an empty main pipeline. If your pipeline does not appear to be empty, that means somebody has used your computer already and saved something in browser storage under the name "default.mlp". To start with a new .mlp document, simply select File->New from the menu.

.. image:: https://www.dropbox.com/s/o4hydnw3rt6ewiy/mlpipeline1.png?dl=1
	:height: 150
	:alt: MLPipeline screenshot...

A pipeline is a list of pipeline steps together with an input-output-parameter specification. To add a step to the main pipeline, click "Add step". As of today there are three types of processing steps: (1) Processor steps; (2) Pipeline steps; and (3) JSON output steps. Let's choose a processor step that does not require any input files. From the dropdown menu, choose "pyms.synthesize_timeseries".

The next dialog box allows you to specify the inputs, outputs, and parameters for this step. For this processor, the two inputs (firings and waveforms) are optional. If these are left blank, the outputsut will be a pure noise dataset. The output (timeseries_out) is mandatory, so type in a name like "raw". The parameters are all optional, but why not specify noise_level=1, samplerate=30000 and duration=10 (seconds).

The main pipeline now has a single step, which we can run by clicking the single-step play button. If you are running the desktop program, this will (fingers crossed) run the process on the local computer. Effectively it will run the mp-queue-process command with the specified parameters. The actual paths of the input and output files are managed internally.

Once completed the file "raw" should appear in the list of output files. That can either be downloaded directly (if not too large) or can be used in subsequent processing steps.

Next we'll add a second step to do a bandpass filter on the synthesized timeseries. Once again click "Add step", and select the processor "ms3.bandpass_filter". For the input "timeseries" enter "raw" (matching the output from the first step) and for "timeseries_out" enter "filt". For parameters, put in samplerate=30000 (Hz), freq_min=300 (Hz) and freq_max=6000 (Hz). Just enter the number, not the units. Once again click the single-step play button and the process should run on the local computer, and a new output named "filt" will appear in the list of output files.

.. image:: https://www.dropbox.com/s/0h541cxzn6lxyj7/mlpipeline2.png?dl=1
	:height: 400
	:alt: MLPipeline screenshot...

The order of the processing steps does not matter because the system will automatically determine the execution order based on the input/output dependencies, although it is good practice to list the steps in logical order.

Using sub-pipelines
-------------------

The second type of processing step is a "pipeline" step. This is similar to a subroutine in a programming language. 




