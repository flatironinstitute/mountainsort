MountainSort Documentation
==========================

MountainSort is spike sorting software developed by Jeremy Magland, Alex Barnett, and Leslie Greengard at the Center for Computational Biology, Flatiron Institute in close collaboration with Jason Chung and Loren Frank at UCSF department of Physiology. The algorithm was featured in

`Chung, Jason E.*, Jeremy F. Magland*, Alex H. Barnett, Vanessa M. Tolosa, Angela C. Tooker, Kye Y. Lee, Kedar G. Shah, Sarah H. Felix, Loren M. Frank, and Leslie F. Greengard. "A Fully Automated Approach to Spike Sorting." Neuron 95, no. 6 (2017): 1381-1394. <http://www.cell.com/neuron/fulltext/S0896-6273(17)30745-6>`_

MountainSort is a plugin package to :doc:`MountainLab <mountainlab/mountainlab>` which is under development by Jeremy Magland and Witold Wysota. MountainLab is a general framework for scientific data analysis, sharing, and visualization.

Getting started with MountainSort
---------------------------------

There are various ways to install and/or use MountainSort. The best choice will depend on how you plan to interact with the program. You can use MountainSort...

* as a plugin package to :doc:`MountainLab <mountainlab/mountainlab>` (ML)
* from the web interface (cloud computing)
* as a standalone program

Below, we will describe installation as a plugin to ML (recommended), and the remarks below will indicate how it could be used as a standalone program.

Supported operating systems
---------------------------

Ubuntu 16.04 is the currently supported development platform. Other Linux flavors should also work. Currently, Mac and Windows are not supported.

If you do not have a linux machine available, we recommend setting up an `Ubuntu virtual machine. <https://help.ubuntu.com/community/VirtualMachines>`_

Installation
------------

The following instructions are for installing MountainSort on Ubuntu 16.04 (recommended). Installation instructions requiring compilation can be :doc:`found here <installation_advanced>`. 

*Note: The following packages have not quite been updated yet. Please wait until Nov 21st*

.. code:: bash

	add-apt-repository -y ppa:magland/mountainlab
	apt-get install mountainlab
	apt-get install mlpipeline
	apt-get install mountainsort
	apt-get install mountainview

Once installed, run the following to choose a temporary directory path. This is where MountainSort will store large intermediate files during processing. Put it somewhere with space.

.. code:: bash

  mlconfig


Testing the installation
------------------------

The first thing to try is

.. code:: bash

  mp-list-processors

This will list the mountainlab processors installed on your system. For example, you should see "ms3.bandpass_filter", "ms3.whiten", and "mountainsortalg.ms3alg". These are among the core steps of the MountainSort spike sorting pipeline.

Next, to get an idea for how processors work, try

.. code:: bash

  mp-spec ms3.bandpass_filter

This will give the specification (inputs/outputs/parameters) for this particular processor.

Next, try the examples in the mountainsort_examples repository

**1. Clone the examples repo:**

.. code:: bash

  git clone https://github.com/flatironinstitute/mountainsort_examples
  cd mountainsort_examples/example1_mlp

**2. Simulate data for the test:**

.. code:: bash

  mlp-run synthesize_v1.mlp synthesize --samplerate=30000 --duration=600 --timeseries=raw.mda --geom=geom.csv --waveforms_true=waveforms_true.mda --num_channels=10 --num_units=50

This will generate test raw data 'raw.mda', geometry data 'geom.csv', and waveform data 'waveforms_true.mda' in the current directory

**3. Sort the test data**

You will now call the mountainsort3 sort pipeline, passing it the newly-created raw data 'raw.mda' and geometry data 'geom.csv'. You will also tell it what to call the output firings, 'firings.mda'. Finally, you will pass it parameters, already in the directory, 'params.json'.

.. code:: bash

  mlp-run mountainsort3.mlp sort --raw=raw.mda --geom=geom.csv --firings_out=firings.mda --_params=params.json

**4. View the test sorting**

The GUI only requires a timeseries, in this case raw data, 'raw.mda', and the firings information (times/labels), 'firings.mda'. We can also pass it the geometry information and samplerate.

.. code:: bash

  mountainview --raw=raw.mda --firings=firings.mda --geom=geom.csv --samplerate=30000

**5. Re-sort the data with automated curation (masking of low-quality clusters and bursting-related merging)**

This time, you will add the automated curation option, '--curate=true'. This will mask out low-quality clusters and do bursting-related merging.

.. code:: bash

  mlp-run mountainsort3.mlp sort --raw=raw.mda --geom=geom.csv --firings_out=firings2.mda --_params=params.json --curate=true

**6. View the curated test sorting**

.. code:: bash

  mountainview --raw=raw.mda --firings=firings2.mda --geom=geom.csv --samplerate=30000

Note that sorting low signal-to-noise ratio data with relabeling may result in there being no apparent clusters (all clusters are of low quality). For this reason, we suggest first sorting your data without curation.
 
You are now ready to sort your own data :doc:`first_sort`
