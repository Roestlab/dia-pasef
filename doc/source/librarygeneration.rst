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



Correcting RT and IM values
^^^^^^^^^^^^^^^^^^^^^^^^^^^


For correcting the retention time and ion mobility values, an iRT peptide 
library file including ion mobility dimension is required. Having the 
MaxQuant output directory (annotated with ion mobiltiy 1/K0) and 
