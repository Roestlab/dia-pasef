Data Conversion
===============

The generated .tdf files from DIA pasef runs can be converted to standard 
formats (mzML) with diapysef.

Using the script ``convertTDFtoMzML.py`` can convert the .tdf file to
a single mzML file. It allows merging of frames for the same precursors, 
filtering of range of frames, splitting of files by overlapping window 
settings, and compression of data with PyMSNumpress.

Inputs
------
-- ``-a``: Analysis directory of the raw data (.d) (Required)
-- ``-o``: Output filename (Required)
-- ``-m``: Number of frames for merging
-- ``-overlap``: Number of overlapping windows
-- ``-r``: Range of frames to convert

For detailed options and descriptions, simple type:

.. code:: bash

   convertTDFtoMzML.py --help


Example
-------

.. code:: bash

   data_dir=diaPasef_run.d
   output_file='diaPasef_run.mzML'

   convertTDFtoMzML.py -a=$data_dir -o=$output_file


The converted mzML files can be processed with the assay library in ``OpenSwathWorkflow``.

