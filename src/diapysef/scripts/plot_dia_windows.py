#!/usr/bin/env python
from __future__ import print_function
from diapysef.plotting import plot_window_layout
import sys
import pandas as pd
"""
Plots the window layout of a tims experiment with precursor density of a pasef experiment.
"""
if len(sys.argv) < 3:
    print("Usage: plot_dia_windows.py window_file precursor_map")
    sys.exit()

windows = pd.read_csv(sys.argv[1])
precursors = pd.read_csv(sys.argv[2])

plot_window_layout(windows = windows, precursor_map = precursors)
