Quantification and Identification
=================================

Using the assay library and the mzML files, identification and quantification of peptides
can be performed with ``OpenSwathWorkflow`` and ``PyProphet``. For detailed description 
and documentation of the downstream analysis, please refer to `their documentation website 
<http://openswath.org/en/latest/docs/binaries.html>`_. The newest ``OpenMS`` version 2.4.0 
includes functionalities in handling ion mobility informations. Here are some of the 
input parameters that are additional to the regular parameters.

Inputs
------
-- ``-ion_mobility_window``: Ion mobility extraction window of precursor
-- ``-im_extraction_window_ms1``: Use ion mobility on MS1 level
-- ``-irt_im_extraction_window``: iRT extraction of the ion mobilty correction values
-- ``-use_ms1_ion_mobility``: Performs extraction on MS1 level ion mobility level
-- ``-Calibration:ms1_im_calibration``: Use ms1 for ion mobility calibration
-- ``-Calibration:im_correction_function``: Choose im correction function
-- ``-Calibration:debug_im_file``: Record the ion mobility correction data
-- ``-Scoring:Scores:use_ion_mobility_scores``: Add ion mobility for scoring


Output
------

``OpenSwathWorkflow`` can generate ``.tsv``, ``.osw`` for identification and scorign output. It 
is also capable of generating the chromatogram files with extension ``.sqmass``. The quantified 
output ``.tsv`` and ``.osw`` can be statistically validated with ``PyProphet``. 

Statistical Validation
----------------------

``PyProphet`` can take the scores generated from ``OpenSwathWorkflow`` and statistically validate 
the precursor identifications. For detailed documentation, please refer to `the website 
<http://openswath.org/en/latest/docs/binaries.html#pyprophet>`_. 

