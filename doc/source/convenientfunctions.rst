Data Access
============

After installation of diapysef and all of its dependencies, make sure that the 
Bruker sdk file is in your current working directory, and you should be able
to import the package or call the scripts.

.. code:: python

   import diapysef


should import all the functions in diapysef, which allows access to the tdf data.

Succefully importing diapysef should result in the following output:

.. code:: python

   Found Bruker sdk. Access to the raw data is possible.


Otherwise, if it results in:

.. code:: python

   Bruker sdk not found. Some functionalities that need access to raw data will 
   not be available. To activate that functionality place libtimsdata.so (Linux) 
   or timsdata.dll in the src folder.


Please put a the Bruker sdk file in the working directory or make a symlink.

Getting Help
------------
To get more detailed description of the available functions in diapysef, you can use 
the python function ``help``.

.. code:: python

   >>> import diapysef
   >>> help(diapysef)
   ...
   ...


Data Structure
==============

The functions of diapysef can either be used in a python shell, or directly executed
as a commandline tool.

Window Layout
^^^^^^^^^^^^^

diapysef is able to read directly into the sqlite-based .tdf raw files and store
information about the acquisition methods.

.. code:: python

   import diapysef as dp
   dia = pd.TimsData(pasefdata)
   win = dia.get_windows()
   win.to_csv("window_layout.csv")
   print("File Written")


Alternatively, the function can be called directly in command line.

.. code:: bash

   get_dia_windows.py pasefdata outputfile


Annotation of Ion Mobility Values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
With our pipeline using MaxQuant DDA data for our spectral library generation, we need
to convert the ion mobility units to 1/K0 values for standardized measure. The MaxQuant
output files evidence and all_peptides have ion mobility values in scan numbers, which
correspond to the ion mobility 1/K0 values. 

Converting the scan numbers to 1/K0 values require a raw .tdf file that contains the 
calibration information.

.. code:: python

   import diapysef as dp
   pas = dp.PasefData(pasefdata)
   mq = dp.PasefMQData(MQout_directory)
   
   mq.get_all_peptides()
   mq.annotate_ion_mobilitty(pas)
   all_pep = mq.all_peptides
   all_pep.to_csv("all_peptides_1K0.csv")

   mq.get_evidence
   mq.annotate_ion_mobility(pas)
   ev = mq.evidence
   ev.to_csv("evidence_1K0.csv")


Alternative, these functions can be called in command line:

.. code:: bash

   annotate_mq_ionmobility.py MQOUT_director Pasefdata output_prefix



Visualization of window placement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The above function writes a simple csv file that contains information of the pasef run
including isolation window width, isolation window starts and ends, ion mobility drift 
time starts and ends, collision energy, etc. With these information, we can plot the 
window placements in the dimensions of m/z and ion mobility.

.. code:: python

   import diapysef as dp
   import pandas as pd

   dia = dp.TimsData(pasefdata)
   win = pd.read_csv("window_layout.csv")
   dp.plot_window_layout(windows = win)


If you have a MaxQuant output of the library precursors, you can also map the precursors
with the windows.

.. code:: python

   import diapysef as dp
   import pandas as pd
   
   dia = dp.TimsData(pasefdata)
   win = pd.read_csv("window_layout.csv")
   precursors = pd.read_csv("all_peptides.csv")
   dp.plot_window_layout(windows = win, precursor_map = precursors)


Or alternatively, the function can be called in command line.

.. code:: bash

   plot_dia_windows.py window_layout.csv annotated_all_peptides_1K0.csv

.. note::

   Precursor map is the allPeptides.csv file generated from MaxQuant and needs to be 
   properly annotated with 1/K0 values.























   
