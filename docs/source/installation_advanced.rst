Installation for developers
==========================

`Installing using the debbian packages <https://mountainlab.readthedocs.org>`_ is recommended for most users. However, if you wish to be able to compile MountainLab, MLPipelines, MountainSort, and MountainView, you will need to install several dependencies first. Again, Linux/Ubuntu and Debian are the currently supported development platforms. Other Linux flavors should also work. Mac and Windows are not currently not supported.

Prerequisites
-------------

If you are on Ubuntu 16.04 or later, you can get away with using package managers to install the prerequisites:

.. code :: bash

	# Note: Run the apt-get commands as root, or using sudo

	# Install qt5
	apt-get update
	apt-get install software-properties-common
	apt-add-repository ppa:ubuntu-sdk-team/ppa
	apt-get update
	apt-get install qtdeclarative5-dev qt5-default qtbase5-dev qtscript5-dev make g++

	# For MLPipeline, you will also need to install the webkit module for Qt
	sudo apt-get install libqt5webkit5-dev

	# Install nodejs and npm
	apt-get update
	apt-get install nodejs npm nodejs-legacy

	# Install python3, python3-pip, and packages
	apt-get update
	apt-get install python3 python3-pip
	pip3 install numpy scipy pybind11 cppimport numpydoc
	# Note: you may want to use a virtualenv or other system to manage your python packages

	# Install docker (optional) for using mldock
	apt-get update
	apt-get install docker.io

	# If you are going to install the mountainsort plugin package, install fftw and sklearn
	apt-get update
	apt-get install libfftw3-dev
	pip3 install sklearn

	# Optionally, you can install matlab or octave
	apt-get update
	apt-get install octave

Otherwise, if you are on a different operating system, use the following links for installing the prequisites:

* :doc:`Qt5 (version 5.5 or later) <qt5_installation>` 
* :doc:`NodeJS <nodejs_installation>`
* :doc:`python3 together with some packages <python installation>` (see above)
* :doc:`FFTW <fftw_installation>`
* Optional: matlab or octave

Compilation
-----------

Important: You should be a regular user when you perform this step -- do not use sudo here or your files will be owned by root.

First time:

.. code :: bash

	git clone https://github.com/flatironinstitute/mountainlab.git
	cd mountainlab
	
	./compile_components.sh

Subsequent updates:

.. code :: bash

	cd mountainlab
	git pull
	./compile_components.sh


You must add mountainlab/bin to your PATH environment variable. For example append the following to your ~/.bashrc file, and open a new terminal (or, source .bashrc):

.. code :: bash

	export PATH=[/path/to/mountainlab]/bin:$PATH

Installing MLPipeline
---------------------

Do the following (after following the prerequisite installation instructions above)

.. code :: bash

	git clone https://github.com/flatironinstitute/mlpipeline.git
	cd mlpipeline
	./compile_components.sh

You must add mlpipeline/bin to your PATH environment variable.
Also add mlpipeline/utils/mlp to your PATH environment variable.

.. code :: bash

	# Then test to see if this opens the GUI:
	mlpipeline

The first time you run this program, some configuration instructions will appear on the window.

Installing the MountainSort plugin package
------------------------------------------

MountainLab packages can be added in one of two ways. They can be added using docker via the "mldock" command, or (preferred for now), by cloning the package repository into the packages/ directory and compiling them there.

For MountainSort, simply do the following (after following the prerequisite installation instructions above)

.. code :: bash
	
	cd mountainlab/packages
	git clone https://github.com/flatironinstitute/mountainsort.git
	cd mountainsort
	./compile_components.sh

	# Then test to see if we have the mountainsort processors
	mp-list-processors

Subsequently, to update the package periodically:

.. code :: bash

	cd mountainlab/packages/mountainsort
	git pull
	./compile_components.sh

Installing MountainView (spike sorting visualization)
-----------------------------------------------------

Do the following (after following the prerequisite installation instructions above)

.. code :: bash

	git clone https://github.com/flatironinstitute/mountainview.git
	cd mountainview
	./compile_components.sh

You must add mountainview/bin to your PATH environment variable.

.. code :: bash

	# Then test to see if this opens the GUI:
	mountainview


Testing the installation
------------------------

If you installed MountainSort as a plugin package to MountainLab, then you should see that the processors have been properly installed by running

.. code:: bash

  mp-list-processors

At the time of writing these docs, I have the following processors:

.. code:: bash

	magland@dub:~/dev/mountainsort/docs$ mp-list-processors 
	banjoview.cross_correlograms
	mountainsortalg.ms3
	ms3.apply_timestamp_offset
	ms3.apply_whitening_matrix
	ms3.bandpass_filter
	ms3.cluster_metrics
	ms3.combine_cluster_metrics
	ms3.combine_firing_segments
	ms3.combine_firings
	ms3.compute_amplitudes
	ms3.compute_templates
	ms3.compute_whitening_matrix
	ms3.concat_event_times
	ms3.concat_firings
	ms3.concat_timeseries
	ms3.confusion_matrix
	ms3.create_firings
	ms3.create_multiscale_timeseries
	ms3.extract_clips
	ms3.extract_firings
	ms3.isolation_metrics
	ms3.link_segments
	ms3.load_test
	ms3.mask_out_artifacts
	ms3.mv_compute_amplitudes
	ms3.mv_compute_templates
	ms3.mv_extract_clips
	ms3.mv_subfirings
	ms3.reorder_labels
	ms3.run_metrics_script
	ms3.split_firings
	ms3.synthesize_timeseries
	ms3.whiten
	ms3.whiten_clips
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
	spikeview.metrics1
	spikeview.templates

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
