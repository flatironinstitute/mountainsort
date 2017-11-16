Processing System
=================

The processing system is the lowest level component of MountainLab and forms the foundation of the entire project. It allows you to:

* Install packages of processors

* Execute those processors from a Linux terminal

* Create and register custom processors in any language (C++, python, matlab/octave, etc.)

* Queue processing jobs for batch analyses

Getting started
===============

First `install MountainLab <https://github.com/magland/mountainlab/blob/master/old/doc/installation.md>`_ on your Linux machine.

The following command lists all the processors that are registered in the system

.. code:: bash

	> mp-list-processors

To support multiple languages (python, matlab, C++, etc) all processors operate on files. So they are arbitrary executables. Each registered processor has a specification, or a spec. Use the mp-spec command to see the specification for a particular processor:


.. code :: bash

	> mp-spec pyms.extract_clips
	{
	    "description": "Extract clips corresponding to spike events",
	    "exe_command": "python3 /home/magland/dev/mountainlab/packages/pymountainsort/basic/basic.py pyms.extract_clips $(arguments)",
	    "has_test": true,
	    "inputs": [
	        {
	            "description": "Path of timeseries mda file (MxN) from which to draw the event clips (snippets)",
	            "name": "timeseries",
	            "optional": false
	        },
	        {
	            "description": "Path of firings mda file (RxL) where R>=2 and L is the number of events. Second row are timestamps.",
	            "name": "firings",
	            "optional": false
	        }
	    ],
	    "name": "pyms.extract_clips",
	    "outputs": [
	        {
	            "description": "Path of clips mda file (MxTxL). T=clip_size",
	            "name": "clips_out",
	            "optional": false
	        }
	    ],
	    "parameters": [
	        {
	            "datatype": "int",
	            "default_value": 100,
	            "description": "(Optional) clip size, aka snippet size, aka number of timepoints in a single clip",
	            "name": "clip_size",
	            "optional": true
	        }
	    ],
	    "version": "0.1"
	}


This tells us what are the inputs / outputs / parameters that can be passed to the processor. Inputs and outputs are always files, and parameters are essentially always strings. Elsewhere I will discuss inputs having a .prv extension in which case the system will search and find the actual file in a different location.

To call a processor, one would for example type:

.. code :: bash

	> mp-run-process pyms.extract_clips --timeseries=raw.mda --firings=firings.mda --clips_out=output.mda --clip_size=123


There are 3 such functions

.. code :: bash

	mp-exec-process
	mp-run-process
	mp-queue-process

The first just calls the processor, plain and simple. The second calls it and remembers the result, caching information about the checksums of the output files, so that if it is run a second time it does not need to recompute. The third queues the process to run at a later time when resources become available.

You can add your own processors (in essentially any language). More on that later.
