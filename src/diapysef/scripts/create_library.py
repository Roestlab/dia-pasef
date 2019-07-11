#!/usr/bin/env python

# Python script to run retention time and ion mobility alignments

#################################################################
#		Import Modules 
################################################################
import diapysef as dp
import pandas as pd
import sys
import os
import click

from diapysef.pasefdata import PasefData
from diapysef.pasefdata import PasefMQData
from diapysef.pasefdata import MQData

###############################################################
#		Load files
##############################################################

@click.command()
#@click.option('--ion_mobility', envvar='ion_mobility', required=True, default=True, help='MaxQuant output contains ion mobility information (set to false for regular MaxQuant output without ion mobility).')
@click.option('--pasefdata', envvar='pasefdata', type=click.Path(exists=True), help='Pasef data file path.')
@click.option('--mqout' , envvar='mqout', required=True, type=click.Path(exists=True), help='Directory from MaxQuant output.')
@click.option('--irt', envvar='irt', type=click.Path(exists=True), help='Full path file of iRT assay library.')
@click.option('--outfile', envvar='outfile', default='mqout.tsv', show_default=True, type=click.Path(exists=False), help='Output assay library file name.')

# Library values alignment
@click.option('--rt_alignment', envvar='rt_alignment', default='nonlinear', show_default=True, type=click.Choice(['linear', 'nonlinear']) , help='Retention time alignment option.')
@click.option('--im_alignment', envvar='im_alignment',default='linear', show_default=True, type=click.Choice(['linear', 'None']), help='Ion Mobility alignment options.')
@click.option('--all_peptides_out', envvar='all_peptides_out', default='annotate1K0_all_peptides.csv', type=click.Path(exists=False), help='Ion mobility annotated all_peptides.csv output file')
@click.option('--evidence_out', envvar='evidence_out', default='annotated1K0_evidence.csv', type=click.Path(exists=False) , help='Ouput evidence.csv file with annotated ion mobility values.')

def main(pasefdata, mqout, irt, outfile, rt_alignment, im_alignment, all_peptides_out, evidence_out):
    mq = MQData(mqout)
    mq.get_evidence()
    ev = mq.evidence
    mq.get_all_peptides()
    all_pep = mq.all_peptides
    if 'Ion mobility index' in ev.columns:
        ion_mobility = True
        mq = PasefMQData(mqout)
        mq.get_all_peptides()
        all_pep = mq.all_peptides
        if 'Calibrated 1/K0' not in ev.columns:
            if pasefdata is None:
                print('WARNING: Ion mobility conversion required. Please specify --pasefdata for ion mobility information')
                sys.exit()
            elif os.path.exists(pasefdata):
                pas = PasefData(pasefdata)
                print('Annotating ion mobility on MQout all_peptides file ...')
                mq.annotate_ion_mobility(pas)
                all_pep = mq.all_peptides
                print('Writing out annotated all_peptides to %s' % all_peptides_out)
                all_pep.to_csv(all_peptides_out)
                mq.get_evidence()
                mq.annotate_ion_mobility(pas)
                print('Annotating ion mobility values for MQout evidence file ...')
                ev = mq.evidence
                print('Writing out annotated evidence file to %s' % evidence_out)
                ev.to_csv(evidence_out)
            else:
                print('WARNING: Ion mobility conversion required. Please specify --pasefdata')
    else:
        ion_mobility = False

    mq.get_msms()
    msms = mq.msms
    pirt = pd.read_csv(irt, sep="\t")

    #############################################################
    #		Perform alignment
    ############################################################
    if ion_mobility is True:
        im_column ='IonMobilityIndexK0'
    else:
        im_column = None
        im_alignment = None
    ptsv = dp.pasef_to_tsv(ev, msms, irt_file=pirt, ion_mobility=ion_mobility, im_column=im_column, rt_alignment=rt_alignment, im_alignment = im_alignment)

    ptsv.to_csv(outfile,sep="\t", index = False)
    print('mqout library generated, proceed to OpenSwathAssayGenerator.')



if __name__ == "__main__":
    main()

