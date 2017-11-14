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

Linux/Ubuntu and Debian are the currently supported development platforms. Other Linux flavors should also work. Mac and Windows are not currently not supported.

If you do not have a linux machine available, we recommend setting up an `Ubuntu virtual machine. <https://help.ubuntu.com/community/VirtualMachines>`_

Installation
------------------------

The following instructions are for using the debbian packages (recommended). Installation instructions for developers (requiring compilation) can be found here :doc: `installation_advanced>`. 

1. Install MountainLab debbian package
2. Install the MLPipeline package into mountainlab
3. Install MountainView package
4. Install MountainSort package

Testing the installation
------------------------

To test the installation:

1. Change directory to the example1_mlp directory (depends on where you installed mountainlab and mountainsort)

.. code:: bash

  cd mountainlab/packages/mountainsort/examples/example1_mlp

2. Simulate data for the test

.. code:: bash

  mlp-run synthesize_v1.mlp synthesize --samplerate=30000 --duration=60 --timeseries=raw.mda --geom=geom.csv --waveforms_true=waveforms_true.mda

3. Sort the test data

.. code:: bash

  mlp-run ../../pipelines/mountainsort3.mlp sort --timeseries=raw.mda --geom=geom.csv --firings_out=firings.mda params.json

3. View the test sorting

.. code:: bash

  view raw.mda geom.csv params.json firings.mda

You are now ready to sort your own data :doc: `first_sort`
