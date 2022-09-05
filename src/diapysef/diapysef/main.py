#!/usr/bin/env python3
import ast
import click
import sys
import pickle as pkl
from datetime import datetime

from .targeted_data_extraction import setup_logger, TargeteddiaPASEFExperiment, generate_coordinates

# Main Command Line Interface
@click.group(chain=True)
@click.version_option()
def cli( ):
    '''
    Mobi-DIK (Ion Mobility DIA Tool-Kit) is a package for analysis of DIA data coupled to ion mobility. 

    Visit http://openswath.org/en/latest/docs/mobi-dik.html for usage instructions and help
    '''


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
@click.option('--out', 'outfile', default=('diapasef_extracted_data.tsv'), show_default=True, type=str, help='Filename to save extracted data to.')
@click.option('--mz_tol', default=20, show_default=True, type=int, help='The m/z tolerance toget get upper and lower bounds arround target mz. Must be in ppm.')
@click.option('--rt_window', default=50, show_default=True, type=int, help='The total window range of RT, i.e. a window of 50 would be 25 points to either side of the target RT.')
@click.option('--im_window', default=0.06, show_default=True, type=float, help='The total window range of IM, i.e. a window of 0.06 would be 0.03 points to either side of the target IM.')
@click.option('--mslevel', default='[1]', show_default=True, cls=PythonLiteralOption, help='list of mslevel(s) to extract data for. i.e. [1,2] would extract data for MS1 and MS2.')
@click.option('--verbose', default=0, show_default=True, type=int, help='Level of verbosity. 0 - just displays info, 1 - display some debug info, 10 displays a lot of debug info.')
@click.option('--log_file', default='diapasef_data_extraction.log', show_default=True, type=str, help='Log file to save console messages.')
@click.option('--threads', default=1, show_default=True, type=int, help='Number of threads to parallelize filtering of spectrums across threads.')
def extract_target( infile, target_coordinates, outfile, mz_tol, rt_window, im_window, mslevel, verbose, log_file, threads ):
  '''
  Extract from the raw data given a set of target coordinates to extract for.
  '''
  
  if target_coordinates.endswith('.pkl'):
    pickle_file = open(target_coordinates, 'rb')
    peptides = pkl.load(pickle_file)
  else:
    raise AssertionError(f"Wrong input type ({target_coordinates}) for target coordinates (--coords)! Has to be a pickle file, see --help for example.")
  
  exp = TargeteddiaPASEFExperiment(infile, peptides, mz_tol, rt_window, im_window, mslevel, verbose, log_file, threads)
  click.echo(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO Loading data...")
  exp.load_data()
  click.echo(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Reducing spectra using targeted coordinates...")
  exp.reduce_spectra("targed_data_extraction.mzML")
  click.echo(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Finished extracting targeted spectra!")
#   exp.save_filtered_tsv(outfile)

# Export filtered mzML to tabular data for plotting
@cli.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='A filtered targeted diaPASEF mzML file.')
@click.option('--out', 'outfile', default=('diapasef_extracted_data.tsv'), show_default=True, type=str, help='Filename to save extracted data to.')
@click.option('--mslevel', default='[1]', show_default=True, cls=PythonLiteralOption, help='list of mslevel(s) to extract data for. i.e. [1,2] would extract data for MS1 and MS2.')
@click.option('--verbose', default=0, show_default=True, type=int, help='Level of verbosity. 0 - just displays info, 1 - display some debug info, 10 displays a lot of debug info.')
@click.option('--log_file', default='mobidik_export.log', show_default=True, type=str, help='Log file to save console messages.')
def export( infile, outfile, mslevel, verbose, log_file ):
    '''
    Export a reduced targeted mzML file to a tsv file
    '''
    exp = TargeteddiaPASEFExperiment(infile, None, None, None, None, mslevel, verbose, log_file, None)
    click.echo(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO Loading data...")
    exp.load_data(is_filtered=True)
    exp.save_filtered_tsv(mslevel, outfile)
    click.echo(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Finished exporting data!")

# Generate a pickle file containing a dictionary of peptide coordinates for targeted data extraction
@cli.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='An OSW file post Pyprophet statiscal validation.')
@click.option('--out', 'outfile', default=('diapasef_extracted_data.tsv'), show_default=True, type=str, help='Filename to save pickle of peptide coordinates.')
@click.option('--run_id', default=None, show_default=True, type=int, help='Run id in OSW file corresponding to DIA-PASEF run. Required if your OSW file is a merged OSW file.')
@click.option('--target_peptides', default=None, show_default=True, cls=PythonLiteralOption, help='list of peptides (ModifiedPeptideSequence) to generate coordinates for. i.e. ["T(UniMod:21)ELISVSEVHPSR", "TELIS(UniMod:21)VSEVHPSR"]')
@click.option('--m_score', default=0.05, show_default=True, type=float, help='QValue to filter for peptides below QValue, generate coordinates for these peptides.')
@click.option('--use_transition_peptide_mapping/--no-use_transition_peptide_mapping', default=False, show_default=True, help='Use the TRANSITION_PEPTIDE_MAPPING when getting PRODUCT MZ, instead of joining on TRANSITION_PRECURSOR_MAPPING.')
@click.option('--use_only_detecting_transitions/--no-use_only_detecting_transitions', default=False, show_default=True, help='Only include product m/z of detecting transitions. i.e do not use identifying transitions.')
@click.option('--verbose', default=0, show_default=True, type=int, help='Level of verbosity. 0 - just displays info, 1 - display some debug info, 10 displays a lot of debug info.')
@click.option('--log_file', default='mobidik_peptide_coordinate_generation.log', show_default=True, type=str, help='Log file to save console messages.')
def generate_coordinates( infile, outfile, run_id, target_peptides, m_score, use_transition_peptide_mapping, use_only_detecting_transitions, verbose, log_file):
    '''
    Generate peptide coordinates for targeted extraction of DIA-PASEF data
    '''
    # Initialise logger
    setup_logger(log_file, verbose)
    click.echo(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO Generating coordinates...")
    generate_coordinates(infile, outfile, run_id, target_peptides, m_score, use_transition_peptide_mapping, use_only_detecting_transitions, verbose)
    click.echo(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] INFO: Finished generating coordinates!")


if __name__ == '__main__':
    cli(obj={})
