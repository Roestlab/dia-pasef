#!/usr/bin/env python
from __future__ import print_function
from diapysef.pasefdata import PasefData
from diapysef.pasefdata import PasefMQData
import sys

"""
Annotate the output of a Maxquant analysis with ion mobility in 1/K0.
Usage: annotate_mq_ionmobility.py mqout_dir pasef_analysis_dir output_prefix
"""

if len(sys.argv) < 4:
    print("Usage: annotate_mq_ionmobility.py mqout_dir pasef_analysis_dir output_prefix")
    sys.exit()

mqout = sys.argv[1]
pasefdata = sys.argv[2]
outfile = sys.argv[3]

pas = PasefData(pasefdata)
mq = PasefMQData(mqout)

mq.get_all_peptides(pas)
mq.get_evidence(pas)

mq.all_peptides.to_csv(outfile + "_all_peptides.csv")
mq.evidence.to_csv(outfile + "_evidence.csv")
