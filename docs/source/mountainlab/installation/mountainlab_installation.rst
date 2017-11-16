MountainLab installation
========================

Supported operating systems
---------------------------

Linux/Ubuntu and Debian are the currently supported development platforms. Other Linux flavors should also work. There are plans to support Mac, but Windows is not currently supported. 

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
	pip3 install numpy scipy pybind11 cppimport
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

Installing MLPipeline
---------------------

Do the following (after following the prerequisite installation instructions above)

.. code :: bash

	git clone https://github.com/flatironinstitute/mlpipeline.git
	cd mlpipeline
	./compile_components.sh

You must add mlpipeline/bin to your PATH environment variable.

.. code :: bash

	# Then test to see if this opens the GUI:
	mlpipeline

The first time you run this program, some configuration instructions will appear on the window.

If you get stuck
----------------

If necessary, contact Jeremy. I'm happy to help, and we can improve the docs. I'm also happy to invite you to the slack team for troubleshooting, feedback, etc.

