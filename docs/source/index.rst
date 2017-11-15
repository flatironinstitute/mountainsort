MountainSort Documentation
==========================

MountainSort is spike sorting software developed by Jeremy Magland, Alex Barnett, and Leslie Greengard at the Center for Computational Biology, Flatiron Institute in close collaboration with Jason Chung and Loren Frank at UCSF department of Physiology. The algorithm was featured in

`Chung, Jason E.*, Jeremy F. Magland*, Alex H. Barnett, Vanessa M. Tolosa, Angela C. Tooker, Kye Y. Lee, Kedar G. Shah, Sarah H. Felix, Loren M. Frank, and Leslie F. Greengard. "A Fully Automated Approach to Spike Sorting." Neuron 95, no. 6 (2017): 1381-1394. <http://www.cell.com/neuron/fulltext/S0896-6273(17)30745-6>`_

Getting started with MountainSort
------------

There are various ways to install and/or use MountainSort. The best choice will depend on how you plan to interact with the program. You can use MountainSort...

* as a plugin package to MountainLab (ML)
* from the web interface (cloud computing)
* as a standalone program

Below, we will describe installation as a plugin to ML (recommended), and the remarks below will indicate how it could be used as a standalone program. MountainLab is a general framework for scientific data analysis, sharing, and visualization.

Supported operating systems
---------------------------

Ubuntu 16.04 is the currently supported development platform. Other Linux flavors should also work. Currently, Mac and Windows are not supported.

If you do not have a linux machine available, we recommend setting up an `Ubuntu virtual machine. <https://help.ubuntu.com/community/VirtualMachines>`_

Installation
------------------------

The following instructions are for installing MountainSort on Ubuntu 16.04 (recommended). Installation instructions for developers (requiring compilation) can be found here :doc:`installation_advanced`. 

1. Install MountainLab
2. Install MLPipeline
3. Install MountainView
4. Install MountainSort
5. Configure the directory structure. From the terminal, run:

.. code:: bash

  mlconfig

Follow the directions to choose a temporary directory path. This is where MountainSort will store large intermediate files during processing. Put it somewhere with space.

Testing the installation
------------------------

To test the installation, run the following commands from the terminal:

1. Change directory to the example1_mlp directory (depends on where you installed mountainlab and mountainsort).

.. code:: bash

  cd mountainlab/packages/mountainsort/examples/example1_mlp

2. Simulate data for the test

.. code:: bash

  mlp-run synthesize_v1.mlp synthesize --samplerate=30000 --duration=60 --raw=raw.mda --geom=geom.csv --waveforms_true=waveforms_true.mda

This will generate test raw data 'raw.mda', geometry data 'geom.csv', and waveform data 'waveforms_true.mda' in the current directory

3. Sort the test data

You will now call the mountainsort3 sort pipeline, passing it the newly-created raw data 'raw.mda' and geometry data 'geom.csv'. You will also tell it what to call the output firings, 'firings.mda'. Finally, you will pass it parameters, already in the directory, 'params.json'.

.. code:: bash

  mlp-run ../../pipelines/mountainsort3.mlp sort --raw=raw.mda --geom=geom.csv --firings_out=firings.mda --_params=params.json

3. View the test sorting

The GUI only requires a timeseries, in this case raw data, 'raw.mda', and the firings information (times/labels), 'firings.mda'.
.. code:: bash

  mountainview --raw=raw.mda --firings=firings.mda

4. Re-sort the data with automated curation (masking of low-quality clusters and bursting-related merging)

This time, you will add the automated curation option, '--curate=true'. This will mask out low-quality clusters and do bursting-related merging.

.. code:: bash

  mlp-run ../../pipelines/mountainsort3.mlp sort --raw=raw.mda --geom=geom.csv --firings_out=firings2.mda --_params=params.json --curate=true

5. View the curated test sorting

.. code:: bash

  mountainview --raw=raw.mda --firings=firings2.mda

Note that sorting low signal-to-noise ratio data with relabeling may result in there being no apparent clusters (all clusters are of low quality). For this reason, we suggest first sorting your data without curation.
 
You are now ready to sort your own data :doc:`first_sort`
