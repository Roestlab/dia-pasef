#!/usr/bin/env python
from __future__ import print_function
import pickle
import os.path
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.offsetbox import AnchoredText
import matplotlib
import seaborn as sns
from tqdm import tqdm
# Logging
import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)
from pyopenms import *
from scipy.signal import savgol_filter
import cv2

from scipy.signal import find_peaks, peak_widths
# from two_dimension_gradient_ascent_peak_finder import get_one_dimension_peaks

# Plotting
matplotlib.use('Agg')

# Data

def get_one_dimension_peaks(arr):  
  # Get mean of 1D arrays for im and rt
  im_acr_rt = np.mean(arr, axis=1)
  rt_acc_im = np.mean(arr, axis=0)
  # Smooth 1D arrays for peak picking
  im_acr_rt = savgol_filter(im_acr_rt, window_length=9, polyorder=3 )
  rt_acc_im = savgol_filter(rt_acc_im, window_length=5, polyorder=3 )
  # Get peaks and boundaries for im array
  im_peaks, im_pk_prop = find_peaks(im_acr_rt, prominence=1)
  im_pk_widths, im_pk_width_heights, im_pk_left, im_peak_right = peak_widths(im_acr_rt, im_peaks, rel_height=0.98)
  # Get peaks and boundaries for rt array
  rt_peaks, rt_pk_prop = find_peaks(rt_acc_im, prominence=1)
  rt_pk_widths, rt_pk_width_heights, rt_pk_left, rt_peak_right = peak_widths(rt_acc_im, rt_peaks, rel_height=0.98)
  
  return (im_peaks, im_pk_left, im_peak_right, im_pk_width_heights), (rt_peaks, rt_pk_left, rt_peak_right, rt_pk_width_heights)

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

def plot_2d_qaunt_results_check(arr: np.array, arr_blur: np.array, i: list, j: list, masked: np.array, mask: np.array, label_mask: any=None, fname: str="2D_Quantification_check.png", print_plot: bool=False):
    """
    Plot results from 2D quantification to check performance

    Parameters:
    arr: (np.array) A (n x m) 2D array of intensities
    arr_blur: (np.array) A (n x m) 2D array, same size as arr, but is a smooth/blurred array.
    i: (list) list of x values for seeds. Obtained from get_max method.
    j: (list) list of y values for seeds. Obtained from get_max method.
    masked: (np.masked_array) A (n x m) 2D  masked array
    mask: (np.array) A (n x m) 2D  array of trues and falses. Essentially the same as masked, but not of class masked_array.
    label_mask: (np,array) A (n x m) 2D array of labels
    fname: (str) Filename to save plot as.

    Return:
    None
    """
    # Get unique feature labels
    if label_mask is not None:
        unique_feature_labels = np.unique(label_mask)[np.logical_not(np.isnan(np.unique(label_mask)))]
        # Generate a palette of n unique feature label colors
        palette = sns.color_palette("Dark2", len(unique_feature_labels))
    plt.close()
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, figsize=(10, 10), sharex=True, sharey=True)
    fig.set_tight_layout(True)
    ax1.imshow(arr, aspect='auto', cmap="afmhot_r")
    ax1.set_title('Original')
    ax2.imshow(arr_blur, aspect='auto', cmap="afmhot_r")
    ax2.set_title('Smoothed')
    ax3.imshow(arr, aspect='auto', cmap="afmhot_r")
    ax3.set_title('Seeds')
    # plt.imshow(image, origin='lower', aspect='auto')
    ax3.plot(i,j,'go', markersize=10, alpha=0.8)
    ax4.imshow(masked, 'jet', interpolation='none', alpha=0.7, aspect='auto')
    ax4.set_title('Smoothed\nMask')
    if label_mask is not None:
        ax5.imshow(mask, aspect='auto')
        for label, col in zip(unique_feature_labels, palette):
            ax5.scatter(np.where(label_mask==label)[1], np.where(label_mask==label)[0],color=col, s=15)
        ax5.set_title('Labeled\nMask')
        ax6.imshow(arr, aspect='auto', cmap="afmhot_r")
        for label, col in zip(unique_feature_labels, palette):
            ax6.scatter(np.where(label_mask==label)[1], np.where(label_mask==label)[0],color=col, s=15)
        ax6.set_title('Labeled\nFeature')
        # for y, x, label, col in zip(i, j, unique_feature_labels, palette):
            # ax6.text(y, x, str(int(label)), color="white", fontsize=10, bbox=dict(fill="gray", edgecolor="green", linewidth=2))
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.2, hspace=0.2)
    if print_plot:
      plt.show()
    # Save Figure
    fig.savefig(fname, dpi=fig.dpi)

def plot_2d_qaunt_results_check_watershed(arr, arr_blur, im_arr, rt_arr, im_prominence, rt_prominence, pep_coord, mask, label_mask, i, j, fname: str="2D_Quantification_check.png", print_plot: bool=False):

    
    
    plt.close("all")
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3,2, figsize=(10, 15), sharex=False, sharey=False)
    if pep_coord is not None and 'decoy' in pep_coord.keys():
      fig.suptitle(f"{pep_coord['peptide']}_{pep_coord['charge']}_(decoy={pep_coord['decoy']})")
    elif pep_coord is not None:
      fig.suptitle(f"{pep_coord['peptide']}_{pep_coord['charge']}")
    fig.supxlabel('Retention Time (s)')
    fig.supylabel('Ion Mobility', y=0.35)
    fig.set_tight_layout(True)
    
    # Get 1D data
    im_peak_data, rt_peak_data = get_one_dimension_peaks(arr_blur, im_prominence=im_prominence, rt_prominence=rt_prominence)
    
    # XIC
    palette = sns.color_palette("Dark2", len(rt_peak_data[0]))
    # plt.close("all")
    # fig, ax1 = plt.subplots()
    ax1.plot(rt_arr, np.mean(arr, axis=0))
    oned_quant_str = ""
    label_id = 1
    for left, right, apex, height, col in zip(rt_peak_data[1].astype(int), rt_peak_data[2].astype(int), rt_peak_data[0].astype(int), rt_peak_data[3], palette):
      # print(col)
      ax1.vlines(rt_arr[left], ymin=0, ymax=np.max(np.mean(arr_blur, axis=0)), color=col)
      ax1.vlines(rt_arr[right], ymin=0, ymax=np.max(np.mean(arr_blur, axis=0)), color=col)
      ax1.hlines(height, rt_arr[left], rt_arr[right], color=col)
      ax1.plot(rt_arr[apex], np.mean(arr, axis=0)[apex], 'x')
      ax1.text(rt_arr[apex], np.mean(arr, axis=0)[apex], str(int(label_id)), color="white", fontsize=6, bbox=dict(fill="gray", edgecolor="green", linewidth=1))
      # np.sum(np.sum(arr, axis=0)[left:right])
      if oned_quant_str=="":
        oned_quant_str = f"F{label_id} Int: {np.round(np.sum(np.mean(arr, axis=0)[left:right]))}"
      else:
        oned_quant_str = oned_quant_str + "\n" + f"F{label_id} Int: {np.round(np.sum(np.mean(arr, axis=0)[left:right]))}"
      label_id+=1
    at = AnchoredText(oned_quant_str, prop=dict(size=8), frameon=True, loc='upper right')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax1.add_artist(at)
    if pep_coord is not None and 'intensity' in pep_coord.keys():
      ax1.axvline(pep_coord['rt_apex'], color='r')
      at = AnchoredText(f"EG Int:\n{np.round(pep_coord['intensity'])}", prop=dict(size=8), frameon=True, loc='upper left')
      at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
      ax1.add_artist(at)
    ax1.set_xlabel('Retention time (s)')
    ax1.set_ylabel('Mean Intensity')
    ax1.set_title('XIC')
    
    # XIM
    palette = sns.color_palette("Dark2", len(im_peak_data[0]))
    # plt.close("all")
    # fig, ax2 = plt.subplots()
    ax2.plot(im_arr, np.mean(arr, axis=1))
    oned_quant_str = ""
    label_id = 1
    for left, right, apex, height, col in zip(im_peak_data[1].astype(int), im_peak_data[2].astype(int), im_peak_data[0].astype(int), im_peak_data[3], palette):
      # print(col)
      ax2.vlines(im_arr[left], ymin=0, ymax=np.max(np.mean(arr_blur, axis=1)), color=col)
      ax2.vlines(im_arr[right], ymin=0, ymax=np.max(np.mean(arr_blur, axis=1)), color=col)
      ax2.hlines(height, im_arr[left], im_arr[right], color=col)
      ax2.plot(im_arr[apex], np.mean(arr, axis=1)[apex], 'x')
      ax2.text(im_arr[apex], np.mean(arr, axis=1)[apex], str(int(label_id)), color="white", fontsize=6, bbox=dict(fill="gray", edgecolor="green", linewidth=1))
      # np.sum(np.sum(arr, axis=0)[left:right])
      if oned_quant_str=="":
        oned_quant_str = f"F{label_id} Int: {np.round(np.sum(np.mean(arr, axis=0)[left:right]))}"
      else:
        oned_quant_str = oned_quant_str + "\n" + f"F{label_id} Int: {np.round(np.sum(np.mean(arr, axis=0)[left:right]))}"
      label_id+=1
    at = AnchoredText(oned_quant_str, prop=dict(size=8), frameon=True, loc='upper right')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax2.add_artist(at)
    if pep_coord is not None:
      ax2.axvline(pep_coord['im_apex'], color='r')
    ax2.set_xlabel('Ion Mobility')
    ax2.set_ylabel('Mean Intensity')
    ax2.set_title('XIM')
    
    # 2D data
    unique_feature_labels = np.unique(label_mask[label_mask!=0])
    palette = sns.color_palette("Dark2", len(unique_feature_labels))
    
    # Original data
    ax3.imshow(arr, aspect='auto', cmap="afmhot_r", extent=[np.min(rt_arr), np.max(rt_arr), np.max(im_arr), np.min(im_arr)])
    ax3.set_title('Original')
    
    # Mask
    ax4.imshow(mask, aspect='auto', cmap="afmhot_r", extent=[np.min(rt_arr), np.max(rt_arr), np.max(im_arr), np.min(im_arr)])
    ax4.plot(rt_arr[i], im_arr[j], 'go', markersize=6, alpha=0.8)
    ax4.set_title('Mask and Seeds')
    
    # Colred Label
    ax5.imshow(label_mask, aspect='auto', cmap="jet", extent=[np.min(rt_arr), np.max(rt_arr), np.max(im_arr), np.min(im_arr)])
    ax5.set_title('Labelled Mask')
    
    # Labeled data
    ax6.imshow(arr, aspect='auto', cmap="afmhot_r", extent=[np.min(rt_arr), np.max(rt_arr), np.max(im_arr), np.min(im_arr)])
    
    for label, col in zip(unique_feature_labels, palette):
        # ax5.scatter(rt_arr[np.where(label_mask==label)[1]], im_arr[np.where(label_mask==label)[0]],color=col, s=15, alpha=0.6)
        contours, _  = cv2.findContours(np.array(label_mask==label).astype(np.uint8), cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
        contour = contours[0]
        # contour = np.array([[[rt_arr[v[0][0]], im_arr[v[0][1]]]] for v in contour])
        M = cv2.moments(contour)
        # x = int(M["m10"] / M["m00"])
        # y = int(M["m01"] / M["m00"])
        xs = [rt_arr[v[0][0]] for v in contour]
        ys = [im_arr[(v[0][1])] for v in contour]
        ax6.plot(xs, ys, color='g', linewidth=6)
        # ax5.contour(label_mask==label)
    for y, x, label, col in zip(i, j, unique_feature_labels, palette):
        ax6.text(rt_arr[y], im_arr[x], str(int(label)), color="white", fontsize=6, bbox=dict(fill="gray", edgecolor="green", linewidth=1))
    twod_quant_str = ""
    for label_id in np.unique(label_mask[label_mask!=0]):
      if twod_quant_str=="":
        twod_quant_str = f"F{label_id} Int: {np.round(np.sum(arr[label_mask==label_id]))}"
      else:
        twod_quant_str = twod_quant_str + "\n" + f"F{label_id} Int: {np.round(np.sum(arr[label_mask==label_id]))}"
    at = AnchoredText(twod_quant_str, prop=dict(size=8), frameon=True, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax6.add_artist(at)
    ax6.set_title('Labelled Features')
    if print_plot:
      plt.show()
    # Save Figure
    fig.savefig(fname, dpi=fig.dpi)

def plot_xic_xim_subplot():
    
    # from pyopenms import *

    plt.close('all')
    fig, ((ax1, ax2)) = plt.subplots(1,2, figsize=(10, 10), sharex=False, sharey=True)
    fig.set_tight_layout(True)
    
    
    uni_ims = np.unique(data_filt.im)
    # exp_im_indx = np.argmin(np.abs(np.unique(uni_ims) - 0.971191683025569))
    exp_im_indx = np.argmin(np.abs(np.unique(uni_ims) - 0.831202491539261))
    exp_im_indx = np.argmin(np.abs(np.unique(uni_ims) - 1.146577))
    exp_im_indx = np.argmin(np.abs(np.unique(uni_ims) - pep_coord['im_apex']))
    uni_ims[exp_im_indx-3:exp_im_indx+3]
    
    exp = MSExperiment()
    # Set raw data (RT and intensity)
    for im_val in uni_ims[exp_im_indx-3:exp_im_indx+3]:
        # Create new chromatogram
        chromatogram = MSChromatogram()
        data_sub = data_filt[['rt', 'im', 'int']][data_filt.im==im_val]
        rt = data_sub.rt.to_numpy()
        intensity = data_sub.int.to_numpy()
        chromatogram.set_peaks([rt, intensity])
        # Sort the peaks according to ascending retention time
        chromatogram.sortByPosition()
        # Add meta information to the chromatogram
        # chromatogram.setNativeID(filter_peptides)
        # Add chromatogram to experiment
        exp.addChromatogram(chromatogram)
    
    
    for chrom in exp.getChromatograms():
        retention_times, intensities = chrom.get_peaks()
        ax1.plot(retention_times, intensities, label = chrom.getNativeID())
    
    # ax1.axvline(1221.90222167969, color='r')
    # ax1.axvline(1246.90295410156, color='r')
    # ax1.axvline(298.951477050781, color='r')
    # ax1.axvline(322.169372558594, color='r')
    # ax1.axvline(13.23825*60, color='r')
    ax1.axvline(pep_coord['rt_apex'], color='r')
    ax1.set_title('XIC')
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('intensity (cps)')
    ax1.legend()
    
    uni_ims = np.unique(data_filt.rt)
    # exp_im_indx = np.argmin(np.abs(np.unique(uni_ims) - 1234.55))
    # exp_im_indx = np.argmin(np.abs(np.unique(uni_ims) - 308.921))
    exp_im_indx = np.argmin(np.abs(np.unique(uni_ims) - pep_coord['rt_apex']))
    
    exp = MSExperiment()
    # Set raw data (RT and intensity)
    for im_val in uni_ims[exp_im_indx-3:exp_im_indx+3]:
        # Create new chromatogram
        chromatogram = MSChromatogram()
        data_sub = data_filt[['rt', 'im', 'int']][data_filt.rt==im_val]
        rt = data_sub.im.to_numpy()
        intensity = data_sub.int.to_numpy()
        chromatogram.set_peaks([rt, intensity])
        # Sort the peaks according to ascending retention time
        chromatogram.sortByPosition()
        # Add meta information to the chromatogram
        # chromatogram.setNativeID(filter_peptides)
        # Add chromatogram to experiment
        exp.addChromatogram(chromatogram)
    
    for chrom in exp.getChromatograms():
        retention_times, intensities = chrom.get_peaks()
        ax2.plot(retention_times, intensities, label = chrom.getNativeID())
    
    # plt.axvline(1221.90222167969, color='r')
    # plt.axvline(1246.90295410156, color='r')
    # ax2.axvline(0.971191683025569, color='r')
    # ax2.axvline(0.831202491539261, color='r')
    ax2.axvline(pep_coord['im_apex'], color='r')
    ax2.set_title('XIM')
    ax2.set_xlabel('time (s)')
    ax2.set_ylabel('intensity (cps)')
    ax2.legend()
    # plt.tight_layout()
    fig.suptitle(filter_peptides)
    # fig.tight_layout()
    plt.show()
    plt.close('all')


class two_dimension_plotter:
  """
  Class to produce RT and IM 2D plots
  """
  def __init__(self, arr, im_arr=None, rt_arr=None, pep_coord=None, labelled_mask=None, precursor_id=None, fname="two_dimension_report.png", print_plot=False):
    """
    Initialize object

    params:
        arr: (numpy.ndarray) (im, rt) 2D array to plot
        im_arr: (numpy.ndarray) (im,) 1D array of ion mobility values. (Optional) Will use indexes if None
        rt_arr: (numpy.ndarray) (rt,) 1D array of retention values. (Optional) Will use indexes if None
        pep_coord: (dict) information for current precursor to plot. i.e. {'peptide': 'T(UniMod:21)ELISVSEVHPSR', 'precursor_mz': 767.3691, 'charge': 2, 'rt_apex': 1730.08, 'im_apex': 1.026132868499893, 'qvalue': 0.0, 'product_mz': [496.2627, 811.4057, 910.4741, 997.5061, 1110.5902, 1223.6743], 'product_charge': [1, 1, 1, 1, 1, 1], 'product_annotation': ['y4^1', 'y7^1', 'y8^1', 'y9^1', 'y10^1', 'y11^1'], 'product_detecting': [1, 1, 1, 1, 1, 1], 'rt_boundaries': [1718.036865234375, 1751.983642578125]}
        labelled_mask: (numpy.ndarray) (im, rt) labelled feature mask. 0's are background (not a feature), 1, 2, 3, ... n are labelled features
        precursor_id: (str) string of current precursor
        fname: (str) name to save file as
        print_plot: (bool) print plot if using an interactive session
        
    Returns:
        None
    """
    
    self.arr = arr
    self.using_im_shape_indexes = False
    self.using_rt_shape_indexes = False
    if im_arr is None:
      im_arr = np.arange(arr.shape[0])
      self.using_im_shape_indexes = True
    self.im_arr = im_arr
    if rt_arr is None:
      rt_arr = np.arange(arr.shape[1])
      self.using_rt_shape_indexes = True
    self.rt_arr = rt_arr
    self.pep_coord = pep_coord
    self.labelled_mask = labelled_mask
    self.precursor_id = precursor_id
    self.fname = fname
    self.print_plot = print_plot
    
    # Figure stuff
    self.fig = plt.figure(figsize=(10,10))
    
  def plot(self):
    """
    Plot Generation
    """
    # Generate single plot of arr
    
    gs = gridspec.GridSpec(3, 3)
    if self.pep_coord is not None and 'decoy' in self.pep_coord.keys():
      self.fig.suptitle(f"{self.pep_coord['peptide']}_{self.pep_coord['charge']}_(decoy={self.pep_coord['decoy']})")
    elif self.pep_coord is not None:
      self.fig.suptitle(f"{self.pep_coord['peptide']}_{self.pep_coord['charge']}")
    elif self.precursor_id is not None:
      self.fig.suptitle(self.precursor_id)
    self.fig.set_tight_layout(True)
    
    # Set gridspace
    main_plot = plt.subplot(gs[1:3, :2])
    xic_plot = plt.subplot(gs[0, :2],sharex=main_plot)
    xim_plot = plt.subplot(gs[1:3, 2],sharey=main_plot)
    cbar_ax = plt.subplot(gs[0, 2])
    
    # two dimension arr
    main_plot_img = main_plot.imshow(self.arr, aspect='auto', cmap="afmhot_r", extent=[np.min(self.rt_arr), np.max(self.rt_arr), np.max(self.im_arr), np.min(self.im_arr)])
    ## Add labelled Contours if available
    if self.labelled_mask is not None:
      unique_feature_labels = np.unique(self.labelled_mask[self.labelled_mask!=0])
      palette = sns.color_palette("Dark2", len(unique_feature_labels))
      twod_quant_str = ""
      for label, col in zip(unique_feature_labels, palette):
        tmp_mask = np.array(self.labelled_mask==label).astype(np.uint8)
        # blur mask to try help generate better contours.
        tmp_mask = cv2.medianBlur(tmp_mask, 5)
        contours, _  = cv2.findContours(tmp_mask, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
        contour = contours[0]
        M = cv2.moments(contour)
        xs = [self.rt_arr[v[0][0]] for v in contour]
        ys = [self.im_arr[(v[0][1])] for v in contour]
        main_plot.plot(xs, ys, color=col, linewidth=3)
        tmp_arr_feature = np.zeros(self.arr.shape)
        tmp_arr_feature[tmp_mask.astype(bool)] = self.arr[tmp_mask.astype(bool)]
        y, x = np.unravel_index(np.argmax(tmp_arr_feature, axis=None), self.arr.shape)
        main_plot.text(self.rt_arr[x], self.im_arr[y], str(int(label)), color="white", fontsize=6, bbox=dict(fill="gray", edgecolor="green", linewidth=1))
        if twod_quant_str=="":
          twod_quant_str = f"F{label} Int: {np.round(np.sum(tmp_arr_feature)):,}"
        else:
          twod_quant_str = twod_quant_str + "\n" + f"F{label} Int: {np.round(np.sum(tmp_arr_feature)):,}"
      at = AnchoredText(twod_quant_str, prop=dict(size=8), frameon=True, loc='upper left')
      at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
      main_plot.add_artist(at)
    main_plot.set_xlabel('Retention Time (s)')
    main_plot.set_ylabel('Ion Mobility')
    
    # one dimension arr 1
    xic_plot.plot(self.rt_arr, np.sum(self.arr, axis=0))
    if self.pep_coord is not None and not self.using_rt_shape_indexes:
      xic_plot.vlines(self.pep_coord['rt_boundaries'][0], ymin=0, ymax=np.max(np.sum(self.arr, axis=0)), color="r")
      xic_plot.vlines(self.pep_coord['rt_boundaries'][1], ymin=0, ymax=np.max(np.sum(self.arr, axis=0)), color="r")
      xic_plot.plot(self.pep_coord['rt_apex'], np.max(np.sum(self.arr, axis=0)), 'x')
      # Add Summed Intensity
      left_index = np.argmin(np.abs(self.rt_arr - self.pep_coord['rt_boundaries'][0]))
      right_index = np.argmin(np.abs(self.rt_arr - self.pep_coord['rt_boundaries'][1]))
      oned_quant_str = r"$\sum_{i=left}^{right} Int$: " + f"{np.round(np.sum(np.sum(self.arr, axis=0)[left_index:right_index])):,}"
      oned_quant_str = oned_quant_str + f"\n Apex Int: {np.round(np.max(np.sum(self.arr, axis=0))):,}"
      at = AnchoredText(oned_quant_str, prop=dict(size=8), frameon=True, loc='upper left')
      at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
      xic_plot.add_artist(at)
    xic_plot.tick_params(bottom=False, labelbottom=False)
    xic_plot.set_ylabel('Summed Intensity')
    
    # one dimension arr 2
    xim_plot.plot(np.sum(self.arr, axis=1), self.im_arr)
    if self.pep_coord is not None and not self.using_im_shape_indexes:
      xim_plot.plot(np.max(np.sum(self.arr, axis=1)), self.pep_coord['im_apex'], 'x')
      oned_quant_str = f"Apex Int: {np.round(np.max(np.sum(self.arr, axis=1))):,}"
      at = AnchoredText(oned_quant_str, prop=dict(size=8), frameon=True, loc='upper right')
      at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
      xim_plot.add_artist(at)
    xim_plot.tick_params(left=False, labelleft=False)
    xim_plot.set_xlabel('Summed Intensity')
    
    # Add color bar
    self.fig.colorbar(main_plot_img, cax=cbar_ax)
    
    if self.print_plot:
      plt.show()
    
    self.fig.savefig(self.fname, dpi=self.fig.dpi)

  
