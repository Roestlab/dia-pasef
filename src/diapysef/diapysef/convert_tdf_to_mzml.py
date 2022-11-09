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
#import matplotlib.pyplot as plt
from ctypes import cdll
from tqdm import tqdm

from .timsdata import TimsData
from .merge_consumer import MergeConsumer
from .splitting_consumer import SplittingConsumer
from .util import setCompressionOptions


try:
    if __name__ == "__main__":
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

def store_frame(frame_id, td, conn, exp, verbose=False, compressFrame=True, keep_frames=False):
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
      for pasef scan. (New tdf 5.1 has msms = 9 for pasef scan)
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
    scanBoundariesk0 = (0.0, 0.0)

    in_scan = False
    scandata = None

    # Check whether we have a MS2 or a PASEF scan
    if msms == 2:
        q = conn.execute("SELECT TriggerMass, IsolationWidth, PrecursorCharge, CollisionEnergy FROM FrameMsMsInfo WHERE Frame={0}".format(frame_id))
        tmp = q.fetchone()
        center = float(tmp[0])
        width = float(tmp[1])
        mslevel = 2
    # new tdf 5.1 has pasef scan msms = 9 
    elif msms == 9:
        q = conn.execute("SELECT IsolationMz, IsolationWidth, ScanNumBegin, ScanNumEnd, CollisionEnergy, Frame FROM DiaFrameMsMsWindows INNER JOIN DiaFrameMsMsInfo ON DiaFrameMsMsWindows.WindowGroup = DiaFrameMsMsInfo.WindowGroup WHERE Frame={0} ORDER BY  ScanNumBegin DESC".format(frame_id))
        scandata = q.fetchall()
        tmp = scandata[scan_data_it]
        center = float(tmp[0])
        width = float(tmp[1])
        scan_start = int(tmp[2])
        scan_end = int(tmp[3])
        next_scan_switch = scan_end
        scanBoundariesk0 = td.scanNumToOneOverK0(frame_id, np.array([scan_end, scan_start])) #spectrum 1/k0 boundaries for the first swath in the frame
        # Check if we already are in the new scan (if there is no
        # gap between scans, happens for diaPASEF):
        if next_scan_switch == num_scans:
            next_scan_switch = scan_start
            in_scan = True

        mslevel = 2
    elif msms == 8:
        q = conn.execute("SELECT IsolationMz, IsolationWidth, ScanNumBegin, ScanNumEnd, CollisionEnergy FROM PasefFrameMsMsInfo WHERE Frame={0} ORDER BY ScanNumBegin DESC".format(frame_id))
        scandata = q.fetchall()
        tmp = scandata[scan_data_it]
        center = float(tmp[0])
        width = float(tmp[1])
        scan_start = int(tmp[2])
        scan_end = int(tmp[3])
        next_scan_switch = scan_end
        scanBoundariesk0 = td.scanNumToOneOverK0(frame_id, np.array([scan_end, scan_start]))
        # Check if we already are in the new scan (if there is no
        # gap between scans, happens for diaPASEF):
        if next_scan_switch == num_scans:
            next_scan_switch = scan_start
            in_scan = True

        mslevel = 2
    else:
        # MS1 
        q = conn.execute("select NumScans from Frames where Id={0}".format(frame_id))
        tmp = q.fetchone()
        scan_start = 0
        scan_end = int(tmp[0])
        scanBoundariesk0 = td.scanNumToOneOverK0(frame_id, np.array([scan_end, scan_start])) #1/k0 boundaries for ms1 spectrum 
        # print(scanBoundariesk0[1])

    if verbose:
        print("Frame", frame_id, "mslevel", mslevel, msms, "contains nr scans:", num_scans, "and nr pasef scans", len(scandata) if scandata else -1)
        print("Scandata for PASEF:", scandata)
    if keep_frames:
        next_scan_switch = -1

    # Get the mapping of the ion mobility axis
    scan_number_axis = np.arange(num_scans, dtype=np.float64)
    ook0_axis = td.scanNumToOneOverK0(frame_id, scan_number_axis)

    nr_scans_created = 0
    allmz = []
    allint = []
    allim = []

    # Traverse in reversed order to get low ion mobilities first (and high scan times first)
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
            if next_scan_switch >= 0 and next_scan_switch >= k:

                if verbose:
                    print("Switch to new scan at", k, "/", next_scan_switch, "store scan of size", len(allmz))

                if in_scan:
                    # Only store spectrum when actually inside a scan, skip the "between scan" pushes
                    sframe = handle_compressed_frame(allmz, allint, allim, mslevel, time, center, width, scanBoundariesk0)
                    sframe.setNativeID("frame=%s_scan=%s" % (frame_id, next_scan_switch) )
                    exp.consumeSpectrum(sframe)
                    nr_scans_created += 1

                allmz = []
                allint = []
                allim = []
                if k == 0: continue

                if in_scan:
                    scan_data_it += 1

                    if scan_data_it >= len(scandata):
                        if verbose: print("LEFT the last scan, nothing else to do here")
                        next_scan_switch = -2
                        continue

                    # Already prepare for next scan
                    tmp = scandata[scan_data_it]
                    center = float(tmp[0])
                    width = float(tmp[1])
                    scan_start = int(tmp[2])
                    scan_end = int(tmp[3])
                    scanBoundariesk0 = td.scanNumToOneOverK0(frame_id, np.array([scan_end, scan_start]))

                    in_scan = False
                    next_scan_switch = scan_end

                    if verbose: print("LEAVING scan now, next scan starts at:", next_scan_switch)

                    # Check if we already are in the new scan (if there is no
                    # gap between scans, happens for diaPASEF):
                    if k == next_scan_switch:
                        if verbose: print("STARTING new scan immediately at", k, ":",  center - width/2.0, center + width/2.0, "scan will end at:", next_scan_switch)
                        next_scan_switch = scan_start
                        in_scan = True
                else:
                    in_scan = True
                    next_scan_switch = scan_start
                    if verbose: print("STARTING new scan at", k, ":",  center - width/2.0, center + width/2.0, "scan will end at:", next_scan_switch)

            continue

        # Store data in OpenMS Spectrum file -> each TOF push is an individual
        # spectrum and we store the ion mobility in the precursor. The frame
        # can be reconstructed by grouping all spectra with the same RT.
        s = pyopenms.MSSpectrum()
        s.setMSLevel(mslevel)
        s.set_peaks( (mz, intens) ) 
        s.setRT(time)
        s.setNativeID("frame=%s spec %s" % (frame_id, k) )
        p = pyopenms.Precursor()
        p.setDriftTime(drift_time)
        if mslevel == 2:
            p.setMZ(center)
            p.setIsolationWindowUpperOffset(width / 2.0)
            p.setIsolationWindowLowerOffset(width / 2.0)
        s.setPrecursors([p])
        exp.consumeSpectrum(s)

    # Store data compressed for cases where the whole frame represents a single spectrum (e.g. MS1)
    if compressFrame and next_scan_switch == -1:
        sframe = handle_compressed_frame(allmz, allint, allim, mslevel, time, center, width, scanBoundariesk0)
        sframe.setNativeID("frame=%s" % frame_id)
        exp.consumeSpectrum(sframe)
        nr_scans_created += 1

    if scandata is not None and (nr_scans_created != len(scandata)):
        raise Exception("Something went quite wrong here, we expected", len(scandata), "scans, but only created", nr_scans_created)

def handle_compressed_frame(allmz, allint, allim, mslevel, rtime, center, width, scanBoundariesk0):
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

    p = pyopenms.Precursor()

    if mslevel == 2:
        p.setMZ(center)
        p.setIsolationWindowUpperOffset(width / 2.0)
        p.setIsolationWindowLowerOffset(width / 2.0)

    # set IM boundaries of spectrum, set here to be consistsent with proteowizard
    sframe.setMetaValue('ion mobility lower limit', scanBoundariesk0[0])
    sframe.setMetaValue('ion mobility upper limit', scanBoundariesk0[1])

    sframe.setPrecursors([p])
    sframe.set_peaks( (mz, intens) )
    sframe.setFloatDataArrays([fda])
    sframe.sortByPosition()

    return sframe

def get_consumer(output_fname):
    # Store output
    if output_fname.lower().endswith("mzml"):
        consumer = pyopenms.PlainMSDataWritingConsumer(output_fname)

        # Compress output
        try:
            opt = consumer.getOptions()
            setCompressionOptions(opt)
            consumer.setOptions(opt)
        except Exception as e:
            print(e)
            print("Your version of pyOpenMS does not support any compression, your files may get rather large")
            pass

    elif output_fname.lower().endswith("sqmass"):
        consumer = pyopenms.MSDataSqlConsumer(output_fname)

    else:
        raise Exception("Supported filenames: mzML and sqMass.")

    return consumer

def convert_diapasef_tdf_to_mzml(analysis_dir, output_fname, merge_scans, keep_frames, verbosity, overlap_scans, frame_limit):
    """
    Main method to perform conversion form bruker tdf data to mzML format
    """

    if sys.version_info.major == 2:
        analysis_dir = unicode(analysis_dir)

    # merge_scans = -1
    # if len(sys.argv) > 3:
    #     merge_scans = int(sys.argv[3])
    #     print("Will merge", merge_scans, "scans")

    td = TimsData(analysis_dir)
    conn = td.conn

    # Get total frame count:
    q = conn.execute("SELECT COUNT(*) FROM Frames")
    row = q.fetchone()
    N = row[0]
    print("Analysis has {0} frames.".format(N))

    consumer = get_consumer(output_fname)

    if merge_scans != -1:
        consumer = MergeConsumer(consumer, merge_scans)

    if overlap_scans > 1:
        
        # For overlapping scans, we need to create N different consumers (N files on disk) 
        # with different file names which will then contain the overlapped spectra 
        
        consumers = []
        for k in range(overlap_scans):
            outspl = output_fname.rsplit(".", 1)
            outname = outspl[0] + "_" + str(k) + "." + outspl[1]
            consumer = get_consumer(outname)
            consumers.append(consumer)

        if merge_scans != -1:
            m_consumers = []
            for c in consumers:
                consumer = MergeConsumer(c, merge_scans)
                m_consumers.append(consumer)
            consumers = m_consumers # use the merge consumers in the overlap consumer

        consumer = SplittingConsumer(consumers)

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

    for frame_id in tqdm(range(lower_frame, upper_frame)):
        store_frame(frame_id+1, td, conn, consumer, compressFrame=True, verbose=verbosity > 1, keep_frames=keep_frames)
    


def main():

    parser = argparse.ArgumentParser(description ="Conversion program to convert a Bruker raw data file from a timsTOF Pro instrument into a single mzML.")
    parser.add_argument("-a", "--analysis_dir",
                        help = "The location of the directory containing raw data (usually .d)",
                        dest = 'analysis_dir',
                        required = True)
    parser.add_argument("-o", "--output_name",
                        help = "The name of the output file (mzML)",
                        dest = "output_fname",
                        required = True)
    parser.add_argument("-m", "--merge",
                        help = "Number of consecutive frames to sum up (squash). This is useful to boost S/N if exactly repeated frames are measured.",
                        type = int,
                        default = -1,
                        dest = "merge_scans")
    parser.add_argument("--keep_frames",
                        help = "Whether to store frames exactly as measured or split them into individual spectra by precursor isolation window (default is to split them - this is almost always what you want).",
                        type = bool,
                        default = False,
                        dest = "keep_frames")
    parser.add_argument("--verbose",
                        help = "Verbosity",
                        type = int,
                        default = -1,
                        dest = "verbosity")
    parser.add_argument("--overlap",
                        help = "How many overlapping windows were recorded for the same m/z window - will split the output into N output files.",
                        type = int,
                        default = -1,
                        dest = "overlap_scans")
    parser.add_argument("-r", "--framerange",
                        help = "The minimum and maximum Frames to convert. Useful to only convert a part of a file.",
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
    overlap_scans = args.overlap_scans
    print("Overlapping scans: ", overlap_scans)
    frame_limit = args.frame_limit
    print("Frame limits: ", frame_limit)

    # if len(sys.argv) < 3:
    #     raise RuntimeError("need arguments: tdf_directory output.mzML")

    # analysis_dir = sys.argv[1]
    # output_fname = sys.argv[2]
    
    convert_diapasef_tdf_to_mzml(analysis_dir, output_fname, merge_scans, args.keep_frames, args.verbosity, overlap_scans, frame_limit)

    print("Conversion completed, press Enter to continue.")

if __name__ == "__main__":
    main()

