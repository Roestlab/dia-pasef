.. diapysef documentation master file, created by
   sphinx-quickstart on Sat Mar  9 16:49:47 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Summary
========
diapysef is a convenience package for working with DIA-PASEF data. It has 
functionalities to convert Bruker raw files into a format that OpenMS can 
understand. Thus OpenSwath can be used to analyze the data and TOPPView can
 be used to visualize. 

The diapysef package contains python bindings for the processing and analysis
of mass spectrometry data coupled to TIMS-TOF based proteomics. It provides
open-source access to algorithms specifically designed for DIA (data-
independent acquisition) experiments. These bindings include algorithms that
allow conversions of raw data (Bruker tdf.d), generation of spectral library, 
basic signal-processing (compression, filtering), and simple visualization 
of the raw data.

Note: The current documentation relates to the 0.3.3 version of diapysef.
Note: Please acquire the Bruker tdf-sdk through official distributions. 

.. toctree::
   :maxdepth: 2
   :caption: Installation
   
    installation

.. toctree::
   :maxdepth: 2
   :caption: Data Access and Visualization
   
   convenientfunctions
   datavisualization

.. toctree::
   :maxdepth: 2
   :caption: Library Generation
   
   librarygeneration

.. toctree::
   :maxdepth: 2
   :caption: Data Conversion
 
   dataconversion

.. toctree::
   :maxdepth: 2
   :caption: Quantification with OpenSwath
    
   openswath


Indices and tables
==================

* :ref:`genindex`

Acknowledgements
----------------
The tools and workflows are collaboratively developed with the data acquired 
by the `Mann Group at MPI, Bremen <https://www.biochem.mpg.de/en/rd/mann>`_ , 
the `Aebersold Group at IMSB, ETH Zurich <http://www.imsb.ethz.ch/research/aebersold.html>`_,
and Bruker Daltonics. The core pipeline is also referenced from `OpenMS <http://www.openms.org>`_ 
, `PyProphet <https://github.com/PyProphet>`_, and `msproteomicstools <https://github.com/msproteomicstools>`_.
