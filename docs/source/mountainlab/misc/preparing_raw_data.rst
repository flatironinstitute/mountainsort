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

[Request for help: could somebody write more details on this]