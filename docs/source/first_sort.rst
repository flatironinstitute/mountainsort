Sorting your own data
==========================

Sorting your own data will require a few steps that will be outlined on this page:
    1. Convert raw data into the .mda file format
    2. Create geometry configuration file (geom.csv)
    3. Specify parameters (ie. sampling rate)
    4. Select the sorting pipeline
    5. Sort the data
    6. Visualize the output

Preparing raw data
==================

The first step of spike sorting using MountainSort is to prepare a raw $M\times N$ timeseries dataset in .mda format. Here $M$ is the number of electrode channels and $N$ is the number of timepoints. Note that if your electrode array can be split into multiple independent channel subsets, then you should sort each of these subsets separately.

Converting from raw binary
--------------------------

If your data is in raw binary format then you should use the pyms.extract_timeseries processor. For example, if you had a tetrode dataset located at /path/to/raw.dat in raw int16 format, then you could run the following to generate a new file called raw.mda:

.. code :: bash

	mp-run-process pyms.extract_timeseries --timeseries=/path/to/raw.dat --timeseries_out=raw.mda --timeseries_dtype=int16 --timeseries_num_channels=4

There are also options for extracting a subset of channels (1-based indexing) or timepoints (0-based indexing). For more information on the various capabilities of this processor, run the following:

.. code :: bash

	mp-spec pyms.extract_timeseries

*Could some users please contribute some examples/tutorials of how to convert their specific data into .mda format. Ultimately we'd love python processor plugins to do specific data conversions.*

Using matlab
------------

If your data is in a format that can be loaded into matlab, then you can create the raw.mda file using the writemda16i or writemda32 matlab functions, as follows:

.. code :: matlab

	cd mountainlab/matlab
	run mountainlab_setup.m

	% prepare your MxN array called X %

	writemda16i(X,'raw.mda');

	% If you need to save as float32 type, use writemda32(X,'raw.mda') %

However, this may not be suitable for very large datasets. Therefore the command-line procedure is generally preferred.

Using python3/numpy
-------------------

If your data is in a format that can be loaded using python3/numpy, then you can create the raw.mda using the functions provided in packages/pymountainsort/mlpy

You will want to use the writemda functions

.. code :: python

	import numpy as np
        cd mountainlab/packages/mountainsort/packages/pyms

	from mlpy import writemda16i, writemda32

	# prepare your MxN array called X

	writemda16i(X,'raw.mda')

	# If you need to save as float32 type, use writemda32(X,'raw.mda')

*If this is how you convert your data, please consider sharing the code, and ultimately we will create a mountainlab processor for this task.*

Here are some more specifics on the :doc:`MDA file format <mda_file_format>`


Specifying the electrode array geometry
=======================================

The geometry of the electrode array in relation to how it is stored in the raw.mda file is done in a geom.csv file, containing the 2D or 3D locations of the electrodes. It is a comma-separated text file where each line represents an electrode channel, and the columns are the geometric coordinates. These coordinates can be in any unit so long as they correspond to the adjacency_radius sorting parameter. For complex geometries, it is encouraged to use microns. 

For example a tetrode might have the following geom.csv:

.. code ::

    0,0
    -25,25
    25,25
    0,50
    
This file is used in conjunction with the adjacency_radius sorting parameter and determines the local electrode neighborhoods. 

If adjacency_radius=-1, or the geom.csv is not present, then there is only one electrode neighborhood containing all the channels. 

If the adjacency_radius=0, then each channel is sorted independently. 

This file is also used by the viewer for display.

Specifying the recording and sorting parameters
===============================================

params.json contains sorting parameters that are specific to the dataset. At a minimum it should contain the sample rate in Hz. This file would appear as follows:

.. code ::

    {"samplerate":30000}

You can also specify whether to look for positive spike peaks (detect_sign=1), negative (detect_sign=-1), or both (detect_sign=0). Typically, most electrophysiology datasets are made up almost entirely of negative spikes. This is specified as follows:


.. code ::

    {"samplerate":30000, "detect_sign":-1}


Select the sorting pipeline
===========================

With MountainSort, there is the mountainsort3 pipeline included, and you can also build your own pipeline :doc:`processing_pipelines`.

The mountainsort3 pipeline is found in 'mountainlab/packages/mountainsort/pipelines'

Sort the data
=============

You will now call the sorting pipeline, passing it the paths to the timeseries, geometry information, and parameters files. Assuming that you are running it from the directory where all the files are (including the .mlp file), and mountainlab was installed in your home directory:

.. code:: bash

  mlp-run mountainsort3.mlp sort --raw=raw.mda --geom=geom.csv --firings_out=firings2.mda --_params=params.json --curate=true

View the output
===============

You can launch the sorting results in the MountainView GUI using:

.. code:: bash

    mountainview --raw=raw.mda --geom=geom.csv --firings=firings.mda --samplerate=30000

Other arguments can be passed to mountainview, allowing for other timeseries (filtered and preprocessed/whitened data) and metrics to be viewed.

All arguments are the paths to the relevant file

.. code:: bash

    mountainview --raw=raw.mda --filt=filt.mda --pre=pre.mda --geom=geom.csv --firings=firings.mda --metrics=metrics.json


Accessing the output
====================
The main output from the sorting is the firings.mda, containing times and labels

Format of the firings.mda
-------------------------

"firings.mda" is the output file containing the times (sample number or index, NOT in seconds) and corresponding labels.

The output of a sorting run is provided in a 2D array usually named "firings.mda". The dimensions are RxL where L is the number of events and R is at least 3.

Each column is a firing event.

The first row contains the integer channels corresponding to the primary channels for each firing event. It is important to note that the channel identification number is relative. In other words, if you only sort channels 61-64, the channel identifications will be 1-4.
This primary identification channel information is optional and can be filled with zeros. It is especially useful for algorithms that sort on (neighborhoods of) individual channels and then consolidate the spike types.

The second row contains the integer time points (1-based indexing) of the firing events

The third row contains the integer labels, or the sorted spike types.

The fourth row (optional and not currently exported by default) contains the peak amplitudes for the firing events.

Further rows may be used in the future for providing reliability metrics for individual events, or quantities that identify outliers.

The curated firings
-------------------

If a pipeline is used that contains a curation step, the firings would have been automatically curated, (ie. having putative noise clusters removed). This will typically have the same name "firings.mda", and will take the same form as described above, but will typically have clusters removed.

Note that curated firings can be empty if all clusters are removed for being of low quality.

Other outputs
-------------

Intermediate files are available after processing, such as the filtered and preprocessed timeseries, metrics, and label map (if curation was used).

Beyond these intermediate files, there are no other expected outputs from pipelines at this time. For example, the clips and templates must be extracted from the timeseries using the information from the firings.mda

Depending on what you are trying to extract from the timeseries based upon the firings, there may already be a processor for what you are looking for or you may have to write it on your own. See :doc:`processing_pipelines` for more details.
