#!/usr/bin/env python3
import ast
import click
import sys
import pickle as pkl
from datetime import datetime

from .util import setup_logger, argument_value_log, check_bruker_sdk
from .targeted_data_extraction import TargeteddiaPASEFExperiment, generate_coordinates
from .convert_tdf_to_mzml import convert_diapasef_tdf_to_mzml
from .plotting import save_report_2d_rt_im_heatmap

# Main Command Line Interface


@click.group(chain=True)
@click.version_option()
def cli():
    '''
    Mobi-DIK (Ion Mobility DIA Tool-Kit) is a package for analysis of DIA data coupled to ion mobility. 

    Visit http://openswath.org/en/latest/docs/mobi-dik.html for usage instructions and help
    '''
    check_bruker_sdk()


# https://stackoverflow.com/a/47730333
class PythonLiteralOption(click.Option):
    def type_cast_value(self, ctx, value):
        if not isinstance(value, str):  # required for Click>=8.0.0
            return value
        try:
            return ast.literal_eval(value)
        except Exception:
            raise click.BadParameter(value)

# Main Targeted DIA-PASEF Data Extraction


@cli.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='Raw data, diaPASEF mzML file.')
@click.option('--coords', 'target_coordinates', required=True, type=click.Path(exists=True), help="""File that contains target coordinates to extract data for, or a file that can be used to generate target coordinates from.\b\n\nCan be one of:\n\npickle (.pkl) - a pickle file that contains a python dictionary of coordinates, i.e. { "TELISVSEVHPS(UniMod:21)R" : {'precursor_mz':767.3691,'charge':2,'product_mz': [311.0639, 312.6357, 322.1867, ..., 11432.6832],'rt_apex':1736.98,'rt_boundaries':[1719.823, 1759.129],'im_apex':0.9884788}, "T(UniMod:21)ELISVSEVHPSR" : {'precursor_mz':767.3691, 'charge':2,'product_mz': [311.0639, 312.6357, 322.1867, ..., 11432.6832],'rt_apex':1730.08,'rt_boundaries':[1718.037, 1751.984],'im_apex':1.0261329}}""")
@click.option('--out', 'outfile', default=('targed_data_extraction.mzML'), show_default=True, type=str, help='Filename to save extracted data to. Must be type mzML')
@click.option('--mz_tol', default=20, show_default=True, type=int, help='The m/z tolerance toget get upper and lower bounds arround target mz. Must be in ppm.')
@click.option('--rt_window', default=50, show_default=True, type=int, help='The total window range of RT, i.e. a window of 50 would be 25 points to either side of the target RT.')
@click.option('--im_window', default=0.06, show_default=True, type=float, help='The total window range of IM, i.e. a window of 0.06 would be 0.03 points to either side of the target IM.')
@click.option('--mslevel', default='[1]', show_default=True, cls=PythonLiteralOption, help='list of mslevel(s) to extract data for. i.e. [1,2] would extract data for MS1 and MS2.')
@click.option('--readOptions', 'readoptions', default='ondisk', show_default=True, type=click.Choice(['ondisk', 'cached']), help='Context to estimate gene-level FDR control.')
@click.option('--verbose', default=0, show_default=True, type=int, help='Level of verbosity. 0 - just displays info, 1 - display some debug info, 10 displays a lot of debug info.')
@click.option('--log_file', default='mobidik_data_extraction.log', show_default=True, type=str, help='Log file to save console messages.')
@click.option('--threads', default=1, show_default=True, type=int, help='Number of threads to parallelize filtering of spectrums across threads.')
def targeted_extraction(infile, target_coordinates, outfile, mz_tol, rt_window, im_window, mslevel, readoptions, verbose, log_file, threads):
    '''
    Extract from the raw data given a set of target coordinates to extract for.
    '''

    # Initialise logger
    setup_logger(log_file, verbose)

    if verbose == 10:
        args_dict = locals()
        argument_value_log(args_dict)

    if target_coordinates.endswith('.pkl'):
        pickle_file = open(target_coordinates, 'rb')
        peptides = pkl.load(pickle_file)
    else:
        raise AssertionError(
            f"Wrong input type ({target_coordinates}) for target coordinates (--coords)! Has to be a pickle file, see --help for example.")

    exp = TargeteddiaPASEFExperiment(
        infile, peptides, mz_tol, rt_window, im_window, mslevel, readoptions, verbose, None, threads)
    click.echo(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Loading data...")
    exp.load_data()
    click.echo(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Reducing spectra using targeted coordinates...")
    exp.reduce_spectra(outfile)
    click.echo(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Finished extracting targeted spectra!")


# Export filtered mzML to tabular data for plotting
@cli.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='A filtered targeted diaPASEF mzML file.')
@click.option('--out', 'outfile', default=('diapasef_extracted_data.tsv'), show_default=True, type=str, help='Filename to save extracted data to.')
@click.option('--mslevel', default='[1]', show_default=True, cls=PythonLiteralOption, help='list of mslevel(s) to extract data for. i.e. [1,2] would extract data for MS1 and MS2.')
@click.option('--verbose', default=0, show_default=True, type=int, help='Level of verbosity. 0 - just displays info, 1 - display some debug info, 10 displays a lot of debug info.')
@click.option('--log_file', default='mobidik_export.log', show_default=True, type=str, help='Log file to save console messages.')
def export(infile, outfile, mslevel, verbose, log_file):
    '''
    Export a reduced targeted mzML file to a tsv file
    '''
    # Initialise logger
    setup_logger(log_file, verbose)

    if verbose == 10:
        args_dict = locals()
        argument_value_log(args_dict)

    exp = TargeteddiaPASEFExperiment(
        infile, None, None, None, None, mslevel, verbose, None, None)
    click.echo(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Loading data...")
    exp.load_data(is_filtered=True)
    exp.save_filtered_tsv(mslevel, outfile)
    click.echo(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Finished exporting data!")

# Generate a pickle file containing a dictionary of peptide coordinates for targeted data extraction


@cli.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='An OSW file post Pyprophet statiscal validation.')
@click.option('--out', 'outfile', default=('peptides_coordinates.pkl'), show_default=True, type=str, help='Filename to save pickle of peptide coordinates.')
@click.option('--run_id', default=None, show_default=True, type=int, help='Run id in OSW file corresponding to DIA-PASEF run. Required if your OSW file is a merged OSW file.')
@click.option('--target_peptides', default=None, show_default=True, type=str, help='list of peptides (ModifiedPeptideSequence) to generate coordinates for. i.e. ["T(UniMod:21)ELISVSEVHPSR", "TELIS(UniMod:21)VSEVHPSR"]. You can alternatively pass a text file containing comma separate peptide sequences.')
@click.option('--m_score', default=0.05, show_default=True, type=float, help='QValue to filter for peptides below QValue, generate coordinates for these peptides.')
@click.option('--use_transition_peptide_mapping/--no-use_transition_peptide_mapping', default=False, show_default=True, help='Use the TRANSITION_PEPTIDE_MAPPING when getting PRODUCT MZ, instead of joining on TRANSITION_PRECURSOR_MAPPING.')
@click.option('--use_only_detecting_transitions/--no-use_only_detecting_transitions', default=True, show_default=True, help='Only include product m/z of detecting transitions. i.e do not use identifying transitions.')
@click.option('--verbose', default=0, show_default=True, type=int, help='Level of verbosity. 0 - just displays info, 1 - display some debug info, 10 displays a lot of debug info.')
@click.option('--log_file', default='mobidik_peptide_coordinate_generation.log', show_default=True, type=str, help='Log file to save console messages.')
def prepare_coordinates(infile, outfile, run_id, target_peptides, m_score, use_transition_peptide_mapping, use_only_detecting_transitions, verbose, log_file):
    '''
    Generate peptide coordinates for targeted extraction of DIA-PASEF data
    '''
    # Initialise logger
    setup_logger(log_file, verbose)

    if verbose == 10:
        args_dict = locals()
        argument_value_log(args_dict)

    click.echo(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Generating coordinates...")
    generate_coordinates(infile, outfile, run_id, target_peptides, m_score,
                         use_transition_peptide_mapping, use_only_detecting_transitions, verbose)
    click.echo(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Finished generating coordinates!")

# Conversion program to convert a Bruker TIMS .d data file to mzML format


@cli.command()
@click.option('--in', 'analysis_dir', required=True, type=click.Path(exists=True), help='The location of the directory containing raw data (usually .d).')
@click.option('--out', 'output_fname', required=True, type=str, help='The name of the output file (mzML).')
@click.option('--merge', 'merge_scans', default=-1, show_default=True, type=int, help='Number of consecutive frames to sum up (squash). This is useful to boost S/N if exactly repeated frames are measured.')
@click.option('--keep_frames/--no-keep_frames', 'keep_frames', default=False, show_default=True, help='Whether to store frames exactly as measured or split them into individual spectra by precursor isolation window (default is to split them - this is almost always what you want).')
@click.option('--verbose', 'verbosity', default=-1, show_default=True, type=int, help='Verbosity.')
@click.option('--overlap', 'overlap_scans', default=-1, show_default=True, type=int, help='How many overlapping windows were recorded for the same m/z window - will split the output into N output files.')
@click.option('--framerange', 'frame_limit', default='[-1, -1]', show_default=True, cls=PythonLiteralOption, help='The minimum and maximum Frames to convert. Useful to only convert a part of a file.')
def convertTDFtoMzML(analysis_dir, output_fname, merge_scans, keep_frames, verbosity, overlap_scans, frame_limit):
    '''
    Conversion program to convert a Bruker TIMS .d data file to mzML format
    '''
    click.echo(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Converting {analysis_dir}...")
    convert_diapasef_tdf_to_mzml(
        analysis_dir, output_fname, merge_scans, keep_frames, verbosity, overlap_scans, frame_limit)
    click.echo(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Finished converting TDF data to mzML!")

# Generate a report for spectic type of plots


@cli.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='Data tsv file that contains data to be plotting. i.e peptide sequence, charge state, m/z, MS level, retention time, ion mobility, and intensity')
@click.option('--out', 'outpdf', required=True, type=str, help='The pdf file name to save the plots to.')
@click.option('--type', 'plot_type', default='rt_im_heatmap', show_default=True, type=click.Choice(['rt_im_heatmap']), help='Type of plot to generate.')
# Plot Type: rt_im_heatmap Parameters
@click.option('--plot_contours/--no-plot_contours', 'plot_contours', default=False, show_default=True, help='Should contour lines be plotted? Arg for type rt_im_heatmap')
@click.option('--verbose', default=0, show_default=True, type=int, help='Level of verbosity. 0 - just displays info, 1 - display some debug info, 10 displays a lot of debug info.')
@click.option('--log_file', default='mobidik_report.log', show_default=True, type=str, help='Log file to save console messages.')
def report(infile, outpdf, plot_type, plot_contours, verbose, log_file):
    '''
    Generate a report for a specfific type of plot
    '''
    # Initialise logger
    setup_logger(log_file, verbose)
    if verbose == 10:
        args_dict = locals()
        argument_value_log(args_dict)
    if plot_type == 'rt_im_heatmap':
        click.echo(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Generating a report of plots for a Retention Time and Ion Mobility Heatmaps...")
        save_report_2d_rt_im_heatmap(infile, outpdf, plot_contours)
        click.echo(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Finished generating report!")
    else:
        raise click.ClickException(f'plot type {plot_type} is not supported')


if __name__ == '__main__':
    cli(obj={})
