MDA file format
===============

Principles of the .mda format
-----------------------------

The .mda file format was created as a simple method for storing multi-dimensional arrays of numbers. Of course the simplest way would be to store the array as a raw binary file, but the problem with this is that fundamental information required to read the data is missing – specifically,

* the data type (e.g., float32, int16, byte, complex float, etc).
* the number of dimensions
* the size of the dimensions (e.g., number of rows and columns in a matrix)

How should this information be included? There are many strategies, but we choose to include these in a minimal binary header.

In contrast to file formats that can hold multiple data entitities, each .mda file is guaranteed to contain one and only one multi-dimensional array of byte, integer, or floating point numbers. The .mda file contains a small well-defined header containing only the minimal information required to read the array, namely the number and size of the dimensions as well as the data format of the entries. Immediately following the header, the data of the multi-dimensional array is stored in raw binary format.

File format description
-----------------------

The .mda file format has evolved slightly over time (for example the first version only supported complex numbers), so please forgive the few arbitrary choices.

The first four bytes contains a 32-bit signed integer containing a negative number representing the data format:

.. code ::

  -1 is complex float32
  -2 is byte
  -3 is float32
  -4 is int16
  -5 is int32
  -6 is uint16
  -7 is double
  -8 is uint32

The next four bytes contains a 32-bit signed integer representing the number of bytes in each entry (okay a bit redundant, I know).

The next four bytes contains a 32-bit signed integer representing the number of dimensions (num_dims should be between 1 and 50).

The next 4*num_dims bytes contains a list of signed 32-bit integers representing the size of each of the dimensions.

That's it! Next comes the raw data.

Reading and writing .mda files
------------------------------

The easiest way to read and write .mda files is by using the readmda and writemda* functions available in matlab or python, or by using the C++ classes for mda i/o.

For example, in matlab you can do the following after setting up the appropriate paths:

.. code :: matlab

  > X=readmda('myfile.mda');
  > writemda32(X,'newfile.mda');
  > writemda16i(X,'newfile_16bit_integer.mda');

The python functions are available by importing the mlpy library (see mountainlab/packages/pymountainsort)

Examples of C++ usage are found in the mountainsortalg package: mountainlab/packages/mountainsortalg

Reading the .mda file header from the command-line
--------------------------------------------------

You can get information about the datatype and dimensions of a .mda file using the "mda" commandline utility as follows:

.. code :: bash

  > mda myfile.mda

