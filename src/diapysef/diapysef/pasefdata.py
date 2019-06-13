#!/usr/bin/env python
from __future__ import print_function
import os, sys
import pandas as pd
import numpy as np
import sqlite3
from diapysef.timsdata import TimsData

class MQData:
    """Read outputs of maxquant for library generation functionalities."""

    def __init__ (self, maxquant_directory):
        if(not os.path.exists(maxquant_directory)):
            raise ValueError("Directory: %s not found. Please make sure that you specified the right directory." % maxquant_directory)
        self.evidence_data = os.path.join(maxquant_directory, "evidence.txt")
        self.msms_data = os.path.join(maxquant_directory, "msms.txt")
        self.all_peptides_data = os.path.join(maxquant_directory, "allPeptides.txt")

    def get_evidence (self):
        """Reads the evidence output of maxquant as pandas dataframe."""
        try:
            self.evidence = pd.read_csv(self.evidence_data, sep = '\t')
        except:
            print("Data File %s not found. Make sure you specified the right directory." % self.evidence_data)

    def get_msms (self):
        """Reads the msms output of maxquant as pandas dataframe."""
        try:
            self.msms = pd.read_csv(self.msms_data, sep = '\t')
        except:
            print("Data File %s not found. Make sure you specified the right directory." % self.msms_data)

    def get_all_peptides (self):
        """Reads the allPeptides output of maxquant as pandas dataframe."""
        try:
            self.all_peptides = pd.read_csv(self.all_peptides_data, sep ='\t')
        except:
            print("Data File %s not found. Make sure you specified the right directory." % self.all_peptides_data)

class PasefMQData(MQData):
    """Reads outputs of maxquant with ion mobility functionalities."""

    def get_evidence (self, timsdata = None):
        """Reads the evidence output of maxquant as pandas dataframe."""
        try:
            self.evidence = pd.read_table(self.evidence_data)
        except:
            print("Data File %s not found. Make sure you specified the right directory." % self.evidence_data)
        if timsdata is not None:
            self.annotate_ion_mobility(timsdata)

    def get_all_peptides (self, timsdata = None):
        """Reads the allPeptides output of maxquant as pandas dataframe."""
        try:
            self.all_peptides = pd.read_table(self.all_peptides_data)
        except:
            print("Data File %s not found. Make sure you specified the right directory." % self.all_peptides_data)
        if timsdata is not None:
            self.annotate_ion_mobility(timsdata)

    def annotate_ion_mobility (self, pasef_data):
        """Adds ion mobility colums to the maxquant results that are currently loaded."""
        if hasattr(self, 'all_peptides'):
            self.all_peptides['IonMobilityIndexK0'] = pasef_data.scanNumToOneOverK0(1,self.all_peptides['Ion mobility index'])

        if hasattr(self, 'evidence'):
            self.evidence['IonMobilityIndexK0'] = pasef_data.scanNumToOneOverK0(1,self.evidence['Ion mobility index'])
    def convert_to_lib(self, irt_file):
        if hasattr(self, 'msms') & hasattr(self, 'evidence'):
            pasef_to_lib(self.evidence, self.msms, irt_file)
        else:
            print("Msms and evidence need to be present for library generation")
            sys.exit()

class PasefData(TimsData):

    # def __init__ (self, analysis_directory):
    #     if sys.version_info.major == 2:
    #         if not isinstance(analysis_directory, unicode):
    #             raise ValueError("analysis_directory must be a Unicode string.")
    #     if sys.version_info.major == 3:
    #         if not isinstance(analysis_directory, str):
    #             raise ValueError("analysis_directory must be a string.")

    #     self.conn = sqlite3.connect(os.path.join(analysis_directory, "analysis.tdf"))

    # def get_conversion_func(self):
    #     q = self.conn.execute("SELECT * FROM TimsCalibration")
    #     calib = q.fetchone()
    #     def convert_im(im):
    #         im = np.array(im)
    #         return(1/(calib[8]+calib[9]/(calib[4]+((calib[5]-calib[4])/calib[3])*(im-calib[6]-calib[2]))))
    #     # Mobility[1/k0] = 1/(c6+c7/(c2+((c3-c2)/c1)*(scanno-c4-c0)))
    #     return convert_im

    # def scanNumToOneOverK0 (self, mzs):
    #     convert_scan_num = self.get_conversion_func()
    #     return(convert_scan_num(mzs))

    def getMetadata(self):
        metadata = pd.read_sql(" SELECT * FROM Frames JOIN PasefFrameMsMsInfo ON Frames.Id = PasefFrameMsMsInfo.Frame JOIN Precursors ON PasefFrameMsMsInfo.Precursor = Precursors.Id" , self.conn)
        self.metadata = metadata

    def filterMetadata(self,
                       imstart=None,
                       imend=None,
                       mzstart=None,
                       mzend=None,
                       rtstart=None,
                       rtend=None,
                       charge=None,
                       metadata=None):
        if metadata is None:
            if self.metadata is None:
                self.getMetadata()
            metadata = self.metadata
            # metadata = pd.read_sql(" SELECT * FROM Frames JOIN PasefFrameMsMsInfo ON Frames.Id = PasefFrameMsMsInfo.Frame JOIN Precursors ON PasefFrameMsMsInfo.Precursor = Precursors.Id" , conn)

        if mzstart is not None:
            scans_filt = metadata[metadata['MonoisotopicMz'] >= mzstart]
        else:
            scans_filt = metadata
        if mzend is not None:
            scans_filt = scans_filt[scans_filt['MonoisotopicMz'] <= mzend]
        if imstart is not None:
            scans_filt = scans_filt[scans_filt['ScanNumBegin'] <= imstart]
        if imend is not None:
            scans_filt = scans_filt[scans_filt['ScanNumEnd'] >= imend]
        if rtstart is not None:
            scans_filt = scans_filt[scans_filt['Time'] >= rtstart]
        if rtend is not None:
            scans_filt = scans_filt[scans_filt['Time'] <= rtend]
        if charge is not None:
            scans_filt = scans_filt[scans_filt['Charge'] == charge]
        return(scans_filt)

        

