#!/usr/bin/env python
from diapysef.timsdata import TimsData
import sys
import pandas as pd
"""
Extracts the SWATH windows from a Bruker raw file.
The resulting windows contain the m/z center and widths and the Ion mobility
start and end coordinates (in scan number and 1/K0 units)"""
experiment = sys.argv[1]
outfile = sys.argv[2]

dia = TimsData(experiment)
windows = dia.get_windows()

windows.to_csv(outfile)
