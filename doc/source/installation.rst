diapysef Installation
====================

diapysef functionalities require the tdf-sdk library file distributed by 
Bruker. Please acquire the tdf-sdk file from Bruker or download ProteoWizard
from their official websites. It is recommended to have the tdf-sdk libtimsdata.so
(linux) or timsdata.dll (windows) in the current working directory.

Binaries
********

The binary packages can be downloaded at `this shared dropbox link 
<https://www.dropbox.com/sh/elmaubry6274ay5/AACTRyA2ixLJ5-ozLN5rv_J5a?dl=0>`_.
It supports python versions 3.6 and 3.7. The package is recommended to be installed 
in a python virtual environment. After downloading the file, you can install it by 
typing

.. code-block:: bash

   pip install diapysef-0.3.3-py3-none-any.whl


Downstream data analysis with ion mobility is processed with OpenMS version 2.4.0 with
OpenSwath bindings. Detailed documentation of ``OpenSWath`` can be found ` here 
<http://openswath.org/en/latest/docs/openswath.html>`_.

Source
******

To download diapysef from :index:`source` and the developmental version, you have to
clone the github repository (not published yet) and compile for the .whl file.

.. code-block:: bash

   git clone git@github.com:Roestlab/dia-pasef.git
   pip install setup.py
   python setup.py sdist bdist_wheel
   pip install diapysef-*.whl


