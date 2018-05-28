#!/usr/bin/env python
from __future__ import print_function

"""Conversion program to convert a Bruker TIMS file to a single mzML

"""
import argparse
import sys
import sqlite3
import time
import pyopenms
import numpy as np
import matplotlib.pyplot as plt
from ctypes import cdll

import diapysef.timsdata

try:
    if sys.platform[:5] == "win32" or sys.platform[:5] == "win64":
        libname = "timsdata.dll"
        dll = cdll.LoadLibrary(libname)
    elif sys.platform[:5] == "linux":
        libname = "libtimsdata.so"
        dll = cdll.LoadLibrary(libname)
    else:
        raise Exception("Unsupported platform.")
except OSError as e:
    print("This functionality can only be carried out if the bruker sdk is present. Please install it first. The sdk can be installed by installing proteowizard(version >=3, http://proteowizard.sourceforge.net), or by placing the a library file in your path (For windows this will be timsdata.dll and for Linux libtimsdata.so).\n")
    sys.exit()

def store_frame(frame_id, td, conn, exp, verbose=False, compressFrame=True):
    """
    Store a single frame as an individual mzML file

    Note that there are two ways to store the data:
      (i) Multiple spectra per frame (for visualization), compressFrame is False. This is the
      easiest way to visualize and process the data but involves a few hacks,
      namely storing the IM axis as the RT of each spectrum.

      (ii) One spectrum per frame, compressFrame is True. This puts all peaks
      into a single spectrum (while storing the IM data in an extra array).
      This is more efficient for storage and allows analysis that is ignorant
      of the IM dimension.

      Note that msms = 2 means that we have an MS2 scan whereas msms = 8 stands
      for pasef scan.
    """
    # Get a projected mass spectrum:
    q = conn.execute("SELECT NumScans, Time, Polarity, MsMsType FROM Frames WHERE Id={0}".format(frame_id))
    tmp = q.fetchone()
    num_scans = tmp[0]
    time = tmp[1]
    pol = tmp[2]
    msms = int(tmp[3])

    center = -1
    width = -1
    next_scan_switch = -1
    mslevel = 1
    scan_data = []
    scan_data_it = 0

    # Check whether we have a MS2 or a PASEF scan
    if msms == 2:
        q = conn.execute("SELECT TriggerMass, IsolationWidth, PrecursorCharge, CollisionEnergy FROM FrameMsMsInfo WHERE Frame={0}".format(frame_id))
        tmp = q.fetchone()
        center = float(tmp[0])
        width = float(tmp[1])
        mslevel = 2
    elif msms == 8:
        q = conn.execute("SELECT IsolationMz, IsolationWidth, ScanNumBegin, ScanNumEnd, CollisionEnergy FROM PasefFrameMsMsInfo WHERE Frame={0} ORDER BY IsolationMz ASC".format(frame_id))
        scandata = q.fetchall()
        tmp = scandata[scan_data_it]
        center = float(tmp[0])
        width = float(tmp[1])
        scan_start = int(tmp[2])
        scan_end = int(tmp[3])
        next_scan_switch = scan_start
        mslevel = 2

    if verbose:
        print("mslevel", mslevel, msms)

    # Get the mapping of the ion mobility axis
    scan_number_axis = np.arange(num_scans, dtype=np.float64)
    ook0_axis = td.scanNumToOneOverK0(frame_id, scan_number_axis)

    allmz = []
    allint = []
    allim = []

    # Traverse in reversed order to get low ion mobilities first
    for k, scan in reversed(list(enumerate(td.readScans(frame_id, 0, num_scans)))):
        index = np.array(scan[0], dtype=np.float64)
        mz = td.indexToMz(frame_id, index)
        intens = scan[1]
        drift_time = ook0_axis [ k ] 
        if compressFrame:
            allmz.append(mz)
            allint.append(intens)
            allim.append( [drift_time for dr_time in mz] )

            # We have multiple MS2 spectra in each frame, we need to separate
            # them based on the information from PasefFrameMsMsInfo which
            # indicates the switch scan and the isolation parameter for each
            # quadrupole isolation.
            if next_scan_switch != -1 and next_scan_switch == k:

                if verbose:
                    print("Switch to new scan at", k, "with mapping", scandata)

                sframe = handle_compressed_frame(allmz, allint, allim, mslevel, time, center, width)
                exp.consumeSpectrum(sframe)
                allmz = []
                allint = []
                allim = []
                if k == 0: continue

                scan_data_it += 1
                tmp = scandata[scan_data_it]
                center = float(tmp[0])
                width = float(tmp[1])
                scan_start = int(tmp[2])
                scan_end = int(tmp[3])
                next_scan_switch = scan_start
            continue

        # Store data in OpenMS Spectrum file -> each TOF push is an individual
        # spectrum and we store the ion mobility in the precursor. The frame
        # can be reconstructed by grouping all spectra with the same RT.
        s = pyopenms.MSSpectrum()
        s.setMSLevel(mslevel)
        s.set_peaks( (mz, intens) ) 
        s.setRT(time)
        p = pyopenms.Precursor()
        p.setDriftTime(drift_time)
        if mslevel == 2:
            p.setMZ(center)
            p.setIsolationWindowUpperOffset(width / 2.0)
            p.setIsolationWindowLowerOffset(width / 2.0)
        s.setPrecursors([p])
        exp.consumeSpectrum(s)

    if compressFrame and next_scan_switch == -1:
        sframe = handle_compressed_frame(allmz, allint, allim, mslevel, time, center, width)
        exp.consumeSpectrum(sframe)

def handle_compressed_frame(allmz, allint, allim, mslevel, rtime, center, width):
    mz = np.concatenate(allmz)
    intens = np.concatenate(allint)
    ims = np.concatenate(allim)

    fda = pyopenms.FloatDataArray()
    fda.setName("Ion Mobility")
    fda.resize(len(mz))
    for k, val in enumerate(ims):
        fda[k] = val

    sframe = pyopenms.MSSpectrum()
    sframe.setMSLevel(mslevel)
    sframe.setRT(rtime)
    sframe.setFloatDataArrays([fda])
    p = pyopenms.Precursor()
    if mslevel == 2:
        p.setMZ(center)
        p.setIsolationWindowUpperOffset(width / 2.0)
        p.setIsolationWindowLowerOffset(width / 2.0)
    sframe.setPrecursors([p])
    sframe.set_peaks( (mz, intens) )
    sframe.sortByPosition()
    return sframe

class MergeConsumer():
    """
        Merging consumer that merges MS2 spectra with the same precursor. The
        number of consecutive spectra to be merged is a user parameter.

        This class merges m/z and intensity coordinates as well as the first
        FloatDataArray.

        Does not implement chromatogram consuming.
    """

    def __init__(self, consumer, merge_nr):
        self._internal_consumer = consumer
        self._spectrum_storage = {}
        self._merge_nr = merge_nr

    def __del__(self):
        """
        Cleanup: write all stored spectra to disk
        """
        for k, v in self._spectrum_storage.items():
            if v:
                merge_spec = self._mergeSpectra(v)
                self._internal_consumer.consumeSpectrum(merge_spec)

    def consumeSpectrum(self, s):
        """
            consume individual spectra:
                - write MS1 spectra directly to disk
                - collect MS2 spectra and write afer merging n spectra together
        """
        if s.getMSLevel() == 1:
            self._internal_consumer.consumeSpectrum(s)
        else:

            mz = int(s.getPrecursors()[0].getMZ()*10)
            tmp = self._spectrum_storage.get(mz, [])
            tmp.append(s)

            if len(tmp) >= self._merge_nr:
                merge_spec = self._mergeSpectra(tmp)
                self._internal_consumer.consumeSpectrum(merge_spec)
                tmp = []

            self._spectrum_storage[mz] = tmp

    def _mergeSpectra(self, tmp):
        """
        Perform spectral merging
         - we pick one spectrum (the first one) as our reference and add
           all other data from the other spectra to it
         - we check that merged spectra have equal precursor m/z and same float array
        """

        merge_spec = pyopenms.MSSpectrum(tmp[0])
        fda = tmp[0].getFloatDataArrays()[0]
        fda.clear()
        allmz = []
        allint = []
        for q in tmp:
            m, i = q.get_peaks()
            allmz.append(m)
            allint.append(i)

            # Sanity checks, precursors of merged spectra and float arrays need to match
            assert q.getPrecursors()[0].getMZ() - merge_spec.getPrecursors()[0].getMZ() < 1e-5
            assert len(q.getFloatDataArrays()) == len(merge_spec.getFloatDataArrays())
            assert q.getFloatDataArrays()[0].getName() == merge_spec.getFloatDataArrays()[0].getName()

            # TODO this is not very efficient, fix in pyOpenMS!
            for d in q.getFloatDataArrays()[0]:
                fda.push_back(d)

        mz = np.concatenate(allmz)
        intens = np.concatenate(allint)

        # create merge spec
        merge_spec.set_peaks( (mz, intens) )
        merge_spec.setFloatDataArrays([fda])
        merge_spec.sortByPosition()
        return merge_spec

def main():

    parser = argparse.ArgumentParser(description ="Conversion program to convert a Bruker TIMS file to a single mzML")
    parser.add_argument("-a", "--analysis_dir",
                        help = "The location of the directory containing raw data (usually .d)",
                        dest = 'analysis_dir',
                        required = True)
    parser.add_argument("-o", "--output_name",
                        help = "The name of the output file",
                        dest = "output_fname",
                        required = True)
    parser.add_argument("-m", "--merge",
                        help = "How many frames to sum up into one",
                        type = int,
                        default = -1,
                        dest = "merge_scans")
    parser.add_argument("-r", "--framerange",
                        help = "The minimum and maximum Frames to convert.",
                        type = int,
                        nargs = 2,
                        default = [-1, -1],
                        dest = "frame_limit")
    args = parser.parse_args()
    print("Running conversion with these parameters:\n ")
    analysis_dir = args.analysis_dir
    print("Raw file directory: ", analysis_dir)
    output_fname = args.output_fname
    print("Output name: ", output_fname)
    merge_scans = args.merge_scans
    print("Scans to merge: ", merge_scans)
    frame_limit = args.frame_limit
    print("Frame limits: ", frame_limit)

    # if len(sys.argv) < 3:
    #     raise RuntimeError("need arguments: tdf_directory output.mzML")

    # analysis_dir = sys.argv[1]
    # output_fname = sys.argv[2]

    if sys.version_info.major == 2:
        analysis_dir = unicode(analysis_dir)

    # merge_scans = -1
    # if len(sys.argv) > 3:
    #     merge_scans = int(sys.argv[3])
    #     print("Will merge", merge_scans, "scans")

    td = diapysef.timsdata.TimsData(analysis_dir)
    conn = td.conn

    # Get total frame count:
    q = conn.execute("SELECT COUNT(*) FROM Frames")
    row = q.fetchone()
    N = row[0]
    print("Analysis has {0} frames.".format(N))

    # Store output
    if output_fname.lower().endswith("mzml"):
        consumer = pyopenms.PlainMSDataWritingConsumer(output_fname)

        # Compress output
        try:
            opt = consumer.getOptions()
            cfg = pyopenms.NumpressConfig()
            cfg.estimate_fixed_point = True
            cfg.numpressErrorTolerance = -1.0 # skip check, faster
            cfg.setCompression(b"linear");
            cfg.linear_fp_mass_acc = -1; # set the desired RT accuracy in seconds
            opt.setNumpressConfigurationMassTime(cfg)
            cfg = pyopenms.NumpressConfig()
            cfg = pyopenms.NumpressConfig()
            cfg.estimate_fixed_point = True
            cfg.numpressErrorTolerance = -1.0 # skip check, faster
            cfg.setCompression(b"slof");
            opt.setNumpressConfigurationIntensity(cfg)
            opt.setCompression(True) # zlib compression
            consumer.setOptions(opt)
        except Exception:
            print("Your version of pyOpenMS does not support any compression, your files may get rather large")
            pass

    if output_fname.lower().endswith("sqmass"):
        consumer = pyopenms.MSDataSqlConsumer(output_fname)

    if merge_scans != -1:
        consumer = MergeConsumer(consumer, merge_scans)

    if frame_limit[0] > -1 and frame_limit[0] <= N:
        lower_frame = frame_limit[0]
    elif frame_limit[0] == -1:
        lower_frame = 0
    else:
        raise ValueError("Lower Frame limit is not in the permitted range of frames")
    if frame_limit[1] > frame_limit[0] and frame_limit[1] <= N:
        upper_frame = frame_limit[1]
    elif frame_limit[1] == -1:
        upper_frame = N
    else:
        raise ValueError("Upper Frame limit is not in the permitted range of frames")

    for frame_id in range(lower_frame, upper_frame):
        store_frame(frame_id+1, td, conn, consumer, compressFrame=True, verbose=False)

if __name__ == "__main__":
    main()

