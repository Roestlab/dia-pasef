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

###############################################################
#		Load files
##############################################################

@click.command()
@click.option('--pasefdata', envvar='pasefdata', required=True, type=click.Path(exists=True), help='Pasef data file path.')
@click.option('--mqout' , envvar='mqout', required=True, type=click.Path(exists=True), help='Directory from MaxQuant output.')
@click.option('--irt', envvar='irt',required=True, type=click.Path(exists=True), help='Full path file of iRT assay library.')
@click.option('--outfile', envvar='outfile', default='mqout.tsv', show_default=True, type=click.Path(exists=False), help='Output assay library file name.')

# Library values alignment
@click.option('--rt_alignment', envvar='rt_alignment', default='nonlinear', show_default=True, type=click.Choice(['linear', 'nonlinear']) , help='Retention time alignment option.')
@click.option('--im_alignment', envvar='im_alignment',default='linear', show_default=True, type=click.Choice(['linear', 'None']), help='Ion Mobility alignment options.')
@click.option('--all_peptides_out', envvar='all_peptides_out', default='annotate1K0_all_peptides.csv', type=click.Path(exists=False), help='Ion mobility annotated all_peptides.csv output file')
@click.option('--evidence_out', envvar='evidence_out', default='annotated1K0_evidence.csv', type=click.Path(exists=False) , help='Ouput evidence.csv file with annotated ion mobility values.')

def main(pasefdata, mqout, irt, outfile, rt_alignment, im_alignment, all_peptides_out,evidence_out):
    
    pas = PasefData(pasefdata)
    mq = PasefMQData(mqout)
    mq.get_all_peptides()
    mq.annotate_ion_mobility(pas)
    print('Annotating ion mobility on MQout all_peptides file ...')
    all_pep = mq.all_peptides

    if all_peptides_out is not None:
        all_pep.to_csv(all_peptides_out)
        print('Writing out annotated all_peptides to %s' % all_peptides_out)

    mq.get_evidence()
    mq.annotate_ion_mobility(pas)
    print('Annotating ion mobility values for MQout evidence file ...')
    ev = mq.evidence
    if evidence_out is not None:
        ev.to_csv(evidence_out)
        print('Writing out annotated evidence file to %s' % evidence_out)

    mq.get_msms()
    msms = mq.msms
    pirt = pd.read_csv(irt, sep="\t")

    #############################################################
    #		Perform alignment
    ############################################################
    ptsv = dp.pasef_to_tsv(ev, msms, irt_file=pirt, im_column="IonMobilityIndexK0", rt_alignment=rt_alignment, im_alignment=im_alignment)

    ptsv.to_csv(outfile,sep="\t")
    print('mqout library generated, proceed to OpenSwathAssayGenerator.')



if __name__ == "__main__":
    main()

