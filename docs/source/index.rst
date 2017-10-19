MountainSort Documentation
==========================

MountainSort is spike sorting software developed by Jeremy Magland, Alex Barnett, and Leslie Greengard at the Center for Computational Biology, Flatiron Institute in close collaboration with Jason Chung and Loren Frank at UCSF department of Physiology. The algorithm was featured in

`Chung, Jason E., Jeremy F. Magland, Alex H. Barnett, Vanessa M. Tolosa, Angela C. Tooker, Kye Y. Lee, Kedar G. Shah, Sarah H. Felix, Loren M. Frank, and Leslie F. Greengard. "A Fully Automated Approach to Spike Sorting." Neuron 95, no. 6 (2017): 1381-1394. <http://www.cell.com/neuron/fulltext/S0896-6273(17)30745-6>`_

Installation and Getting started
--------------------------------

There are various ways to install and/or use MountainSort. The best choice will depend on how you plan to interact with the program. You can use MountainSort...

* via MountainLab as a plugin package
* as a standalone program
* from the web interface (cloud computing)

Installing as a plugin package to MountainLab
---------------------------------------------

MountainLab (ML) is a general framework for scientific data analysis, sharing, and visualization. To install MountainSort as a ML plugin package, you must first install Mountainlab (see `<https://github.com/flatironinstitute/mountainlab.git>`_).

After that, you can either install MountainSort via docker or by cloning the mountainsort repository into the mountainlab/packages directory and then following the instructions below for installing MountainSort as a standalone program.

To install as a docker package, do the following (after installing mountainlab and a recent version of docker of course):

.. code:: bash

  mldock install https://github.com/flatironinstitute/mountainsort.git#master:packages/pyms pyms
  mldock install https://github.com/flatironinstitute/mountainsort.git#master:packages/mountainsortalg mountainsortalg


Installing as a standalone program
----------------------------------

Prerequesites: Qt5, python3 with numpy, pybind11, cppimport, scilab, and sklearn

After installing the above prerequisites, you must compile the mountainsortalg component:

.. code:: bash
  
  cd mountainsort/packages/mountainsortalg
  qmake
  make -j

Note that there is also C++ code in the python part (pyms), but that will get compiled on the fly by cppimport and pybind11.

Testing the installation
------------------------

