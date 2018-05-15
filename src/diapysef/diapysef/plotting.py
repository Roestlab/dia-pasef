#!/usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os.path
import pickle

def plot_window_layout (windows, precursor_map = None, display_sc = False):
    """Plots the windows with an optional background of ms1 features"""
    # if type(precursor_map) == 'PasefMQData':
    #     if(not hasattr(precursor_map, all_peptides)):
    #         precursor_map.get_all_peptides()
    #         precursor_map.annotate_ion_mobility()
    #     mq = precursor_map.all_peptides
    # elif type(precursor_map) == 'DataFrame':
    w = windows

    if precursor_map is None:
        mq = pd.read_pickle(os.path.join(os.path.dirname(__file__), 'data/evidence_example.pickle'))
        # ax = pickle.load(open(os.path.join(os.path.dirname(__file__), 'data/all_peptides_density.pickle'), 'rb'))
    else:
        mq = precursor_map
    if not display_sc:
        mq = mq[mq.Charge > 1]

    f, ax = plt.subplots(figsize = (8,6))
    ax.hist2d(mq['m/z'], mq['IonMobilityIndexK0'], bins = [1000,1000], norm = LogNorm())

    for p in [
            patches.Rectangle(
                (w['IsolationMz'][i]-w['IsolationWidth'][i]/2, w['IMend'][i]),
                w['IsolationWidth'][i],
                w['IMstart'][i] - w['IMend'][i],
                alpha = 0.4) for i in range(len(w))
    ]:
        ax.add_patch(p)

    ax.set(xlim = (min(mq['m/z'].min(), (w['IsolationMz'] - w['IsolationWidth']/2).min()),
                   max(mq['m/z'].max(), (w['IsolationMz'] + w['IsolationWidth']/2).max())),
           ylim = (min(mq['IonMobilityIndexK0'].min(), w['IMend'].min()),
                   max(mq['IonMobilityIndexK0'].max(), w['IMstart'].max())))
    plt.xlabel('m/z')
    plt.ylabel('Ion Mobility (1/K0)')

    plt.show()
