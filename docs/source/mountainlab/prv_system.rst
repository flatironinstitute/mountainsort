PRV system
==========

Universal file pointers (.prv files)
------------------------------------

A .prv file is a very small text file that serves as a universal file pointer to another (possibly very large) file. The original file may be located on the same machine, a local server, or on a remote server somewhere in the cloud.

For example, here is an example file named raw.mda.prv:


.. code:: bash

  {
    "original_checksum": "9a42d2fb7d0e2cfd8c206c2ddc30640472ffab3d",
    "original_fcs": "head1000-da39a3ee5e6b4b0d3255bfef95601890afd80709",
    "original_path": "20160426_kanye_02_r1.nt16.mda",
    "original_size": 666616580,
    "prv_version": 0.1
  }

The fields in this JSON document represent information about the file raw.mda at the time the .prv file was created using the prv-create command-line utility as described below. 

* original_checksum: The sha-1 hash of the entire raw file.

* original_fcs: (optional) a code that can be used to quickly determine whether a given file is a candidate for a match.

* original_path: Only used for information in case the user needs a reminder on the name and location of the raw file at the time the .prv file was created.

* original_size: The size in bytes of the raw file.

* prv_version: Always 0.1 for now.

This small text file can be used as a convienent substitute for the original file on a file system or in a document. If the original file is located somewhere on the local computer, the prv-locate utility may be used to efficiently find the path to that file.

The prv command-line utility
----------------------------

The prv command-line utility can be used to create .prv files, locate the original files at a later date, or download the original files from a remote server.

Use prv-create to create a new .prv file that serves as a pointer to the original file:

.. code :: bash

  > prv-create [source_file_name] [destination_file_name]
  > prv-create raw.mda raw.mda.prv
  
Later, use prv-locate to find the original file, even if it had moved:


.. code :: bash

  > prv-locate [prv_file_name]
  > prv-locate raw.mda.prv

The full path to the file is written to stdout.

Use the mlconfig utility to set up the search paths for prv-locate on your system:

.. code :: bash

  > mlconfig

The prv utility can also be used to locate and/or download files stored on a remote computer that is running a processing node server:

.. code :: bash

  > prv-locate [prv_file_name] --server=[name of processing server]

If found, this utility writes to stdout a url path to the remote file. The file can be downloaded using:

.. code :: bash

  > prv-download raw.mda.prv /path/to/destination/raw.mda --server=[name of processing server]

Using prv files in processing
-----------------------------

It is easy to use .prv files in the :doc:`processing_system`. Simply use the path of the .prv file intead of the orignial in the appropriate fields:

.. code :: bash

	> mp-run-process pyms.extract_clips --timeseries=raw.mda.prv --firings=firings.mda.prv --clips_out=output.mda --clip_size=123

The system will automatically search the local computer for the corresponding file prior to running the processing.
