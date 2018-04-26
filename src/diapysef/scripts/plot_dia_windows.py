#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

"""
Plots the window layout of a tims experiment with precursor density of a pasef experiment.
"""
import sys
if len(sys.argv) < 3:
    print("Usage: plot_dia_windows.py window_file precursor_map")
    sys.exit()

from diapysef.plotting import plot_window_layout
import pandas as pd

windows = pd.read_csv(sys.argv[1])
precursors = pd.read_csv(sys.argv[2])

plot_window_layout(windows = windows, precursor_map = precursors)
