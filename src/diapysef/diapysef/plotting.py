#!/usr/bin/env python
from __future__ import print_function
import pickle
import os.path
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
from tqdm import tqdm
# Logging
import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)

# Plotting
matplotlib.use('Agg')

# Data


def plot_window_layout(windows, precursor_map=None, display_sc=False):
    """Plots the windows with an optional background of ms1 features"""
    # if type(precursor_map) == 'PasefMQData':
    #     if(not hasattr(precursor_map, all_peptides)):
    #         precursor_map.get_all_peptides()
    #         precursor_map.annotate_ion_mobility()
    #     mq = precursor_map.all_peptides
    # elif type(precursor_map) == 'DataFrame':
    w = windows

    if precursor_map is None:
        mq = pd.read_pickle(os.path.join(os.path.dirname(
            __file__), 'data/evidence_example.pickle'))
        # ax = pickle.load(open(os.path.join(os.path.dirname(__file__), 'data/all_peptides_density.pickle'), 'rb'))
    else:
        mq = precursor_map
    if not display_sc:
        mq = mq[mq.Charge > 1]

    f, ax = plt.subplots(figsize=(8, 6))
    ax.hist2d(mq['m/z'], mq['IonMobilityIndexK0'],
              bins=[1000, 1000], norm=LogNorm())

    for p in [
            patches.Rectangle(
                (w['IsolationMz'][i]-w['IsolationWidth'][i]/2, w['IMend'][i]),
                w['IsolationWidth'][i],
                w['IMstart'][i] - w['IMend'][i],
                alpha=0.4) for i in range(len(w))
    ]:
        ax.add_patch(p)

    ax.set(xlim=(min(mq['m/z'].min(), (w['IsolationMz'] - w['IsolationWidth']/2).min()),
                 max(mq['m/z'].max(), (w['IsolationMz'] + w['IsolationWidth']/2).max())),
           ylim=(min(mq['IonMobilityIndexK0'].min(), w['IMend'].min()),
                 max(mq['IonMobilityIndexK0'].max(), w['IMstart'].max())))
    plt.xlabel('m/z')
    plt.ylabel('Ion Mobility (1/K0)')

    plt.show()


def get_2d_heatmap_data(data_long, x_col='rt', y_col='im', z_col='int'):
    """
    Get X, Y and Z data for plotting a heatmap

    Params:
        data_long: (data.frame) data.frame containing data for plotting, must have at least 3 columns to extract data for X, Y and Z
        x_col: (str) the column that should be used for the X-axis
        y_col: (str) the column that should by used for the Y-axis
        z_col: (str) the column that should be used for the Z-axis

    Returns:
        Returns X, Y and Z values of equal dimensions for plotting to grid on a heatmap
    """

    x = data_long[x_col].to_numpy()
    y = data_long[y_col].to_numpy()
    z = data_long[z_col].to_numpy()

    # print(f"X: {x.shape} | Y: {y.shape} | Z: {z.shape}")

    # Pivot table to grid
    pdata = pd.DataFrame(data={'x': x, 'y': y, 'z': z})
    pdata = pdata.pivot_table(index='y', columns='x', values='z')

    X = pdata.columns.to_numpy()
    Y = pdata.index.to_numpy()
    Z = pdata.to_numpy()

    return X, Y, Z


def plot_2d_rt_im_heatmap(data, current_peptide, plot_contours=False, fig=plt.figure(1)):
    """
    Plot a Heatmap of RT and IM

    Params:
        data: (data.frame) data.frame containing ms_level information, retention time, ion mobility and intensity
        current_peptide: (str) string of current peptide being plotting
        plot_contours: (bool) Should contour lines be plotted?
        fig: (matplotlib figure object) figure object from matplotlib pyplot

    Returns:
        None
    """
    plt.suptitle(f"Peptide(z): {current_peptide}")
    # Plot MS1 Data if available
    data_sub = data.loc[(data.ms_level == 1)]
    if data_sub.shape[0] != 0:
        ax1 = plt.subplot(1, 2, 1)
        # Get axis data
        X, Y, Z = get_2d_heatmap_data(data_sub)
        # Plot heatmap
        c = plt.pcolormesh(X, Y, Z, cmap='afmhot_r')
        if plot_contours:
            ax1.contour(X, Y, Z, colors='red', alpha=0.55)
        ax1.set_title('MS-Level = 1')
        fig.colorbar(c, ax=ax1)
    # Plot MS2 Data if available
    data_sub = data.loc[(data.ms_level == 2)]
    if data_sub.shape[0] != 0:
        ax2 = plt.subplot(1, 2, 2)
        # Get axis data
        X, Y, Z = get_2d_heatmap_data(data_sub)
        # Plot heatmap
        c1 = plt.pcolormesh(X, Y, Z, cmap='afmhot_r')
        if plot_contours:
            ax2.contour(X, Y, Z, colors='red', alpha=0.55)
        ax2.set_title('MS-Level = 2')
        fig.colorbar(c1, ax=ax2)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False,
                    bottom=False, left=False, right=False)
    plt.xlabel("Retention Time [sec]")
    plt.ylabel("Ion Mobility [1/K0]")


def save_report_2d_rt_im_heatmap(infile, outpdf="diapasef_rt_im_heatmap.pdf", plot_contours=False):
    """
    Save a report of 2D RT and IM heatmap plots to pdf

    params:
        infile: (str) tsv file containing data from targeted extracted export experiment. Should contain information for peptide, charge, ms_level, retention time, ion mobility and intensity
        outpdf: (str) the pdf file name to save the plots to
        plot_contours: (bool) Should contour lines be plotted?

    Returns:
        None
    """
    data = pd.read_csv(infile, sep="\t")
    # Add a group_id to group peptide and charge
    data['group_id'] = data['peptide'] + '_' + data['charge'].astype(str)

    unique_peptide_groups = np.unique(data[['group_id']])
    with PdfPages(outpdf) as pdf:
        pbar = tqdm(range(len(unique_peptide_groups)))
        pbar_desc = "INFO: Plotting"
        fig = plt.figure(1, figsize=(12.75, 8.25))
        for peptide_group in pbar:
            current_peptide = unique_peptide_groups[peptide_group]

            data_peptide_sub = data.loc[(
                data.group_id.isin([current_peptide]))]

            pbar_desc = f"INFO: Plotting..{current_peptide}"
            pbar.set_description(pbar_desc)
            plot_2d_rt_im_heatmap(
                data_peptide_sub, current_peptide, plot_contours, fig)

            # Save to PDF
            pdf.savefig()  # save on the fly
            plt.clf()  # clear figure once saved
        plt.close()
