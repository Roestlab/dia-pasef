Library Generation
==================

For most data-independent acquisition (DIA) analysis, a well-represented 
spectral library is required for precursors, peptide, and protein 
identifications. Currently, we acquire our library through using deeply
fractionated sample runs that are analyzed in MaxQuant.

The general steps for generating the spectral library are, annotating ion
mobility values in the MQ output, correcting retention time and ion mobilty 
values against iRT peptides, formatting it to OpenSwath read-able format, 
generate the spectral library formats (.TraML, .pqp, .tsv) with OpenSwath.


Annotating Ion Mobility
^^^^^^^^^^^^^^^^^^^^^^^
Generating the library requires the DDA output files from MaxQuant.
* msms.txt
* allPeptides.txt
* evidence.txt

Before generating the assay library, the ion mobility values have to be 
converted from scan numbers in MaxQuant output to the standardized 1/K0 
units. 

Follow the instructions at :doc: `convenientfunctions` for annotating
and writing out files with standardized units.


Correcting RT and IM values
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For correcting the retention time and ion mobility values, an iRT peptide 
library file including ion mobility dimension is required. Having the 
MaxQuant output directory (annotated with ion mobiltiy 1/K0) and iRT file,
we can generate an OpenSwath readable format of library.

.. code:: python

   import diapysef as dp
   import pandas as pd
   
   # read in the MQOUT files
   msms = pd.read_csv("mqout_dir/msms.txt", sep = "\t")
   # read in the annotated evidence file
   evidence = pd.read_csv("mqout_dir/annotated_evidence.csv")
   irt = "irt_file.tsv"
   
   # Alignment of RT and IM values
   ptsv = dp.pasef_to_tsv(evidence, msms, irt_file = irt, \
   im_column="IonMobilityIndexK0", rt_alignment="nonlinear", im_alignment="linear")
   
   # save the library table
   ptsv.to_csv(ptsv, "mqout.tsv", sep = "\t", index = False)


The function allows a few options for the alignment of retention time and ion
mobility values. The input ``rt_alignment`` has options for ``nonlinear`` for a 
lowess alignment, ``linear`` alignment, and ``None``. For ``im_alignment``, it
contains options of ``linear`` or ``None``. For standard LC-MS/MS analysis, we 
suggest alignment using ``nonlinear`` for ``rt_alignment`` and ``linear`` for 
``im_alignment``.

Alternatively, we also have a script as a commandline tool that can be called with:

.. code:: bash

   python create_library.py

For details and options of the script, simply type:

.. code:: bash

   python create_library.py --help


Generating Assay Library
^^^^^^^^^^^^^^^^^^^^^^^^

After generating the library with corrected retention time and ion mobility values,
a standardized assay library can be generated through OpenSwath applications.

Standardize Assay Library
-------------------------

First, the library table can be converted into standard assay formats (.TraML, .pqp, .tsv)
using ``OpenSwathAssayGenerator``.

Inputs
------
``OpenSwathAssayGenerator`` requires several inputs parameters.
-- ``in`` : Input file of library transitions
-- ``out``: Output file with valid extensions
-- ``swath_windows_file``: File contains isolation window information


An example would be this:

.. code:: bash

   OpenSwathAssayGenerator -in mqout.tsv -out hela_assaylib.TraML \ 
                           -swath_windows_file window_setting.txt
   

Generate Decoy
--------------

The assay library can add decoy peptides for statistical validation of scores 
and identification. It can be done through ``OpenSwathDecoyGenerator``.

Inputs
------
``OpenSwathDecoyGenerator`` requires several input parameters.
-- ``in``: Input file of the target assay library
-- ``out``: Output file of the target-decoy library
-- ``method``: Method of generating the decoy
-- ``switchKR``: Boolean of switching the termini of the decoy peptides

An example would be this:

.. code:: bash

   OpenSwathDecoyGenerator -in hela_assaylib.TraML \
                           -out hela_target_decoy_assaylib.TraML \
                           -method pseudo-reverse \
                           -switchKR true


After generating the target decoy library, the assay library file is ready for the
input for ``OpenSwathWorkflow`` parameter, ``tr``.













