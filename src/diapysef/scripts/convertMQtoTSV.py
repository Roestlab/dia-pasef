#!python
from __future__ import print_function
from diapysef.pasefdata import PasefData
from diapysef.pasefdata import PasefMQData
from diapysef.conversions import pasef_to_tsv
import pandas as pd
import sys

"""
Convert a MaxQuant output to a spectral library. It can take a combined output from multiple runs. If
multiple runs are present, the best run (lowest PEP) is selected per precursor for the library. Note:
This script does not perform any selection of the transitions, this has to be done in a second step,
e.g. with the OpenSwathAssayGenerator. Same goes for decoy generation (see OpenSwathDecoyGenerator).
"""

parser = argparse.ArgumentParser(description ="""
Usage:
convertMQtoTSV --maxquant_dir=path/to/mqoutput/
               --output_name=path/to/outfile.tsv
               --irt_file=path/to/irtfile.tsv
               --pasef_dir=path/to/pasef/raw/data/
               --pdf_output_name=rtcalibration_plot
Convert a MaxQuant output to a spectral library. It can take a combined output from multiple runs. If
multiple runs are present, the best run (lowest PEP) is selected per precursor for the library. Note:
This script does not perform any selection of the transitions, this has to be done in a second step,
e.g. with the OpenSwathAssayGenerator. Same goes for decoy generation (see OpenSwathDecoyGenerator).
""")
parser.add_argument("-i", "--maxquant_dir",
                    help = """The directory that contains the MaxQuant files. At least evidence.txt
                    and msms.txt must be present in that directory""",
                    dest = 'mq',
                    required = True)
parser.add_argument("-o", "--output_name",
                    help = "The name of the output file",
                    dest = "output_fname",
                    required = True)
parser.add_argument("-r", "--irt_file",
                    help = "Location of the irt file. Must be .tsv",
                    dest = "irt_file",
                    required = True)
parser.add_argument("-p", "--pasef_dir",
                    help = """Location of the raw Pasef data directory. Will be used to
                    calibrate the ion mobility""",
                    dest = "pasef_data",
                    required = True)
parser.add_argument("-v", "--pdf_output_name",
                    help = "The name of the pdf file that contains the RT calibration visualization",
                    default = "rtcalibration",
                    dest = "pdf_output_fname")
args = parser.parse_args()


pas = PasefData(pasef_data)
mq = PasefMQData(mq)
mq.get_evidence()
mq.get_msms()
mq.annotate_ion_mobility(pas)
msms = mq.msms
evidence = mq.evidence
irt = pd.read_csv(irt_file, sep="\t")
ptsv = pasef_to_tsv(evidence, msms, irt, pdfout = pdf_output_fname, im_column = "IonMobilityIndexK0")

ptsv.to_csv(output_fname, sep = "\t", index=False)
