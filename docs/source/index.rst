MountainSort Documentation
==========================

MountainSort is spike sorting software developed by Jeremy Magland, Alex Barnett, and Leslie Greengard at the Center for Computational Biology, Flatiron Institute in close collaboration with Jason Chung and Loren Frank at UCSF department of Physiology. The algorithm was featured in

`Chung, Jason E., Jeremy F. Magland, Alex H. Barnett, Vanessa M. Tolosa, Angela C. Tooker, Kye Y. Lee, Kedar G. Shah, Sarah H. Felix, Loren M. Frank, and Leslie F. Greengard. "A Fully Automated Approach to Spike Sorting." Neuron 95, no. 6 (2017): 1381-1394. <http://www.cell.com/neuron/fulltext/S0896-6273(17)30745-6>`_

Installation and Getting started
--------------------------------

There are various ways to install and/or use MountainSort. The best choice will depend on how you plan to interact with the program. You can use MountainSort...

* as a plugin package to MountainLab (ML)
* as a standalone program
* from the web interface (cloud computing)

Here we will describe installation as a plugin to ML (recommended), and the remarks below will indicate how it could be used as a standalone program. MountainLab is a general framework for scientific data analysis, sharing, and visualization. See `<https://mountainlab.readthedocs.org>`_ for instructions on installing MountainLab.

After that, you can either install MountainSort via docker as described below, or by cloning the mountainsort repository into the mountainlab/packages directory and then following the compilation instructions below. In that location (mountainlab/packages), ML will automatically detect it as a plugin.

To install as a docker package, do the following (after installing mountainlab and a recent version of docker of course):

.. code:: bash

  mldock install https://github.com/flatironinstitute/mountainsort.git#master:packages/pyms pyms
  mldock install https://github.com/flatironinstitute/mountainsort.git#master:packages/mountainsortalg mountainsortalg
  mldock install https://github.com/flatironinstitute/mountainsort.git#master:packages/ms3 ms3

Note that you (and any users of the software) should be in the "docker" group in order to use docker. (Note that this essentially means that you will have root access to the machine.)

If you want to install it without docker, you must first install the following prerequisites:

* `Qt5 <http://mountainlab.readthedocs.io/en/latest/installation/qt5_installation.html>`_
* Python3 with the needed packages

On Debian (e.g., Ubuntu 16.04) you can do the following to install python and the required packages

.. code:: bash

  sudo apt-get install python3 python3-pip

  # Note: you may want to use a virtualenv or other system to manage your python packages
  pip3 install numpy scipy matplotlib pybind11 cppimport sklearn

Then you must compile mountainsortalg and ms3 using Qt5/C++:

.. code:: bash
  
  cd mountainlab/packages
  git clone https://github.com/flatironinstitute/mountainsort

  cd mountainsort/packages/mountainsortalg
  qmake
  make -j

  cd ../../mountainsort/packages/ms3
  qmake
  make -j

Note that there is also C++ code in the python part (pyms), but that will get compiled on the fly by cppimport and pybind11.

Installing MountainView
-----------------------

You will probably want to visualize your spike sorting results. It is recommended that you install `MountainView <https://github.com/flatironinstitute/mountainview.git>`_.

Testing the installation
------------------------

If you installed MountainSort as a plugin package to MountainLab, then you should see that the processors have been properly installed by running

.. code:: bash

  mp-list-processors

At the time of writing these docs, I have the following processors:

.. code:: bash

	magland@dub:~/dev/mountainsort/docs$ mp-list-processors 
	mountainsortalg.ms3
	pyms.bandpass_filter
	pyms.compute_templates
	pyms.concatenate_firings
	pyms.extract_clips
	pyms.extract_geom
	pyms.extract_timeseries
	pyms.handle_drift_in_segment
	pyms.join_segments
	pyms.normalize_channels
	pyms.synthesize_drifting_timeseries
	pyms.synthesize_random_firings
	pyms.synthesize_random_waveforms
	pyms.synthesize_timeseries

To see the inputs/outputs for each of these registered processors, use the mp-spec command as described in the MountainLab documentation.

The following command will give me a synthetic (pure noise) dataset

.. code:: bash

	mp-run-process pyms.synthesize_timeseries --timeseries_out=sim.mda --duration=10 --samplerate=30000

If successful, then we can check the dimensions and datatype using the "mda" command:

.. code:: bash

	> mda sim.mda
	{
	    "data_type": -3,
	    "data_type_string": "float32",
	    "dims": [4,300000],
	    "header_size": 20,
	    "num_bytes_per_entry": 4,
	    "num_dims": 2
	}

All arrays are stored in the `.mda file format <http://mountainlab.readthedocs.io/en/latest/mda_file_format.html>`_. If you have installed mountainview, you can visualize this pure noise dataset by running

.. code:: bash

	> mountainview --raw=raw.mda --samplerate=30000

We can then filter the timeseries using the pyms.bandpass_filter processor (use mp-spec to determine the proper inputs/outputs).

If you are not using MountainLab, you can still run these commands with a bit more effort because you will not have the assistance of tools such as mp-spec, mp-list-processors, and mda:

.. code:: bash

	packages/pyms/basic/basic.mp pyms.synthesize_timeseries --timeseries_out=sim.mda --duration=10 --samplerate=30000

You can also plunge into the python code itself to use these tools from within your python programs. However, note that the processors operate on files rather than taking numpy arrays as arguments.

If you are more comfortable in Matlab, or if your raw data is loadable into Matlab, ML has utilities for reading and writing .mda files and for wrapping ML processors. For example, the to generate the above data one could also execute (from within matlab):

.. code:: matlab

	cd mountainlab/matlab
	mlsetup

	inputs=struct();
	outputs=struct('timeseries_out','tmp_raw.mda');
	params=struct('duration',10,'samplerate',30000);
	opts=struct;
	mp_run_process('pyms.synthesize_timeseries',inputs,outputs,params,opts);
	X=readmda('tmp_raw.mda');
	disp(size(X));
