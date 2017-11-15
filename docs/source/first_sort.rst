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

The first step of spike sorting using MountainSort is to prepare a raw $M\times N$ timeseries dataset in .mda format. Here M is the number of electrode channels and N is the number of timepoints. 

Converting from raw binary
--------------------------

If your data is in raw binary format then you should use the pyms.extract_timeseries processor. For example, if you had a tetrode dataset located at /path/to/raw.dat in raw int16 format, then you could run the following to generate a new file called raw.mda:

.. code :: bash

	mp-run-process pyms.extract_timeseries --timeseries=/path/to/raw.dat --timeseries_out=raw.mda --timeseries_dtype=int16 --timeseries_num_channels=4

There are also options for extracting a subset of channels (1-based indexing) or timepoints (0-based indexing). For more information on the various capabilities of this processor, run the following:

.. code :: bash

	mp-spec pyms.extract_timeseries

Using matlab
------------

If your data is in a format that can be loaded into matlab, then you can create the raw.mda file using the writemda16i or writemda32 matlab functions, as follows:

.. code :: matlab

	cd mountainlab/matlab
	run mountainlab_setup.m

	% prepare your MxN array called X %

	writemda16i(X,'raw.mda');

	% If you need to save as float32 type, use writemda32(X,'raw.mda') %

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

MDA file format
===============

Principles of the .mda format
-----------------------------

The .mda file format was created as a simple method for storing multi-dimensional arrays of numbers. Of course the simplest way would be to store the array as a raw binary file, but the problem with this is that fundamental information required to read the data is missing – specifically,

* the data type (e.g., float32, int16, byte, complex float, etc).
* the number of dimensions
* the size of the dimensions (e.g., number of rows and columns in a matrix)

How should this information be included? There are many strategies, but we choose to include these in a minimal binary header.

In contrast to file formats that can hold multiple data entitities, each .mda file is guaranteed to contain one and only one multi-dimensional array of byte, integer, or floating point numbers. The .mda file contains a small well-defined header containing only the minimal information required to read the array, namely the number and size of the dimensions as well as the data format of the entries. Immediately following the header, the data of the multi-dimensional array is stored in raw binary format.

File format description
-----------------------

The .mda file format has evolved slightly over time (for example the first version only supported complex numbers), so please forgive the few arbitrary choices.

The first four bytes contains a 32-bit signed integer containing a negative number representing the data format:

.. code ::

  -1 is complex float32
  -2 is byte
  -3 is float32
  -4 is int16
  -5 is int32
  -6 is uint16
  -7 is double
  -8 is uint32

The next four bytes contains a 32-bit signed integer representing the number of bytes in each entry (okay a bit redundant, I know).

The next four bytes contains a 32-bit signed integer representing the number of dimensions (num_dims should be between 1 and 50).

The next 4*num_dims bytes contains a list of signed 32-bit integers representing the size of each of the dimensions.

That's it! Next comes the raw data.

Reading and writing .mda files
------------------------------

The easiest way to read and write .mda files is by using the readmda and writemda* functions available in matlab or python, or by using the C++ classes for mda i/o.

For example, in matlab you can do the following after setting up the appropriate paths:

.. code :: matlab

  > X=readmda('myfile.mda');
  > writemda32(X,'newfile.mda');
  > writemda16i(X,'newfile_16bit_integer.mda');

The python functions are available by importing the mlpy library (see mountainlab/packages/pymountainsort)

Examples of C++ usage are found in the mountainsortalg package: mountainlab/packages/mountainsortalg

Reading the .mda file header from the command-line
--------------------------------------------------

You can get information about the datatype and dimensions of a .mda file using the "mda" commandline utility as follows:

.. code :: bash

  > mda myfile.mda

Specifying the electrode array geometry
==========================

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
==========================

params.json contains sorting parameters that are specific to the dataset. At a minimum it should contain the sample rate in Hz. This file would appear as follows:

.. code ::

    {"samplerate":30000}

You can also specify whether to look for positive spike peaks (detect_sign=1), negative (detect_sign=-1), or both (detect_sign=0). Typically, most electrophysiology datasets are made up almost entirely of negative spikes. This is specified as follows:


.. code ::

    {"samplerate":30000, "detect_sign":-1}


Select the sorting pipeline
==========================

With MountainSort, there is the mountainsort3 pipeline included, and you can also build your own pipeline :doc:`processing_pipelines`.

The mountainsort3 pipeline is found in 'mountainlab/packages/mountainsort/pipelines'

Sort the data
==========================

You will now call the sorting pipeline, passing it the paths to the timeseries, geometry information, and parameters files. Assuming that you are running it from the directory where all the files are, and mountainlab was installed in your home directory:

.. code:: bash

  mlp-run ~/mountainlab/packages/mountainsort/pipelines/mountainsort3.mlp sort --raw=raw.mda --geom=geom.csv --firings_out=firings2.mda --_params=params.json --curate=true

View the output
==========================

You can launch the sorting results in the MountainView GUI using:

.. code ::

    mountainview --raw=raw.mda --geom=geom.csv --firings=firings.mda

Other arguments can be passed to mountainview, allowing for other timeseries (filtered and preprocessed/whitened data) and metrics to be viewed.

All arguments are the paths to the relevant file

.. code ::

    mountainview --raw=raw.mda --filt=filt.mda --pre=pre.mda --geom=geom.csv --firings=firings.mda --metrics=metrics.json


Accessing the output
==========================
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
-------------------------

If a pipeline is used that contains a curation step, the firings would have been automatically curated, (ie. having putative noise clusters removed). This will typically have the same name "firings.mda", and will take the same form as described above, but will typically have clusters removed.

Note that curated firings can be empty if all clusters are removed for being of low quality.

Other outputs
-------------

Intermediate files are available after processing, such as the filtered and preprocessed timeseries, metrics, and label map (if curation was used).

Beyond these intermediate files, there are no other expected outputs from pipelines at this time. For example, the clips and templates must be extracted from the timeseries using the information from the firings.mda

Depending on what you are trying to extract from the timeseries based upon the firings, there may already be a processor for what you are looking for or you may have to write it on your own. See :doc:`processing_pipelines` for more details.
