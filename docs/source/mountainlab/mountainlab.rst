MountainLab Documentation
=========================

MountainLab is data processing, sharing and visualization software for scientists. It is built around MountainSort, a spike sorting algorithm, but is designed to more generally applicable.

MountainLab is built in layers in order to maintain flexibility and simplicity. The bottom layer allows users to run individual processing commands from a Linux terminal. The top layers allow cloud-based data processing and sharing of analysis pipelines and results through web browser interfaces.

Installation
------------

:doc:`Installation instructions for MountainLab and optionally MountainSort/MountainView<installation/mountainlab_installation>`

Software components
-------------------

MountainLab comprises the following components:

* :doc:`processing_system` (mproc)

  - Manages libraries of registered processors
  - Custom processors can be written in any lanuage (C++, python, matlab/octave, etc.)

..

* :doc:`prv_system` (prv)

  - Large data files are represented using tiny .prv text files that serve as universal handles
  - Separation of huge data files from analysis workspace
  - Facilitates sharing of processing results

..

* :doc:`mda_file_format` (mda)

  - Simple binary file format for n-dimensional arrays

..

* :doc:`MLPipeline <mlpipeline/mlpipeline_tutorial>`

  - Create and execute processing pipelines on local machine or via web browser
  - Register your computer as a processing servers

Components specific to spike sorting
------------------------------------

* `MountainSort <http://mountainsort.readthedocs.io>`_

  - Plugin package to MountainLab
  - Provides spike sorting processors

..

* MountainView

  - Desktop visualization of processing results
  - (Currently specific to spike sorting)

..

Miscellaneous topics
--------------------

.. toctree::
   :maxdepth: 2

  misc/preparing_raw_data
  misc/waveform_drift
  misc/cordion_plans

  


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

