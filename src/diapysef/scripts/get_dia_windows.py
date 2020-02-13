#!/usr/bin/env python
from __future__ import print_function
from diapysef.timsdata import DiaPasefData
import sys
import pandas as pd
"""
Extracts the SWATH windows from a Bruker raw file.
The resulting windows contain the m/z center and widths and the Ion mobility
start and end coordinates (in scan number and 1/K0 units)
"""
if len(sys.argv) < 3:
    print("Usage: get_dia_windows.py tims_analysis_directory output_file")
    sys.exit()

experiment = sys.argv[1]
outfile = sys.argv[2]

dia = DiaPasefData(experiment)
windows = dia.get_windows()
windows = windows.loc[:, ~windows.columns.duplicated()]

windows.to_csv(outfile)
