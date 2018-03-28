# -*- coding: utf-8 -*-
"""Conversion program to convert a Bruker TIMS file to a single mzML

"""

import sys

import timsdata, sqlite3, sys, time
import pyopenms
import numpy as np, matplotlib.pyplot as plt

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
    mslevel = 1
    if msms == 2:
        q = conn.execute("SELECT TriggerMass, IsolationWidth, PrecursorCharge, CollisionEnergy FROM FrameMsMsInfo WHERE Frame={0}".format(frame_id))
        tmp = q.fetchone()
        center = float(tmp[0])
        width = float(tmp[1])
        mslevel = 2

    if verbose:
        print "mslevel", mslevel, msms

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
        drift_time = ook0_axis [k] 
        if compressFrame:
            allmz.append(mz)
            allint.append(intens)
            allim.append([drift_time for k in mz])
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
        if msms == 2:
            p.setMZ(center)
            p.setIsolationWindowUpperOffset(width / 2.0)
            p.setIsolationWindowLowerOffset(width / 2.0)
        s.setPrecursors([p])
        exp.consumeSpectrum(s)


    if compressFrame:
        mz = np.concatenate(allmz)
        intens = np.concatenate(allint)
        ims = np.concatenate(allim)
        # print "   leeen", len(mz), len(intens)

        fda = pyopenms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.resize(len(mz))
        for k,val in enumerate(ims):
            fda[k] = val

        sframe = pyopenms.MSSpectrum()
        sframe.setMSLevel(mslevel)
        sframe.setRT(time)
        sframe.setFloatDataArrays([fda])
        p = pyopenms.Precursor()
        if msms == 2:
            p.setMZ(center)
            p.setIsolationWindowUpperOffset(width / 2.0)
            p.setIsolationWindowLowerOffset(width / 2.0)
        sframe.setPrecursors([p])
        sframe.set_peaks( (mz, intens) )
        sframe.sortByPosition()
        exp.consumeSpectrum(sframe)

def main():

    if len(sys.argv) < 3:
        raise RuntimeError("need arguments: tdf_directory output.mzML")

    analysis_dir = sys.argv[1]
    output_fname = sys.argv[2]

    if sys.version_info.major == 2:
        analysis_dir = unicode(analysis_dir)

    td = timsdata.TimsData(analysis_dir)
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
            cfg.setCompression("linear");
            cfg.linear_fp_mass_acc = -1; # set the desired RT accuracy in seconds
            opt.setNumpressConfigurationMassTime(cfg)
            cfg = pyopenms.NumpressConfig()
            cfg = pyopenms.NumpressConfig()
            cfg.estimate_fixed_point = True
            cfg.numpressErrorTolerance = -1.0 # skip check, faster
            cfg.setCompression("slof");
            opt.setNumpressConfigurationIntensity(cfg)
            opt.setCompression(True) # zlib compression
            consumer.setOptions(opt)
        except Exception:
            pass

    if output_fname.lower().endswith("sqmass"):
        consumer = pyopenms.MSDataSqlConsumer(output_fname)

    for frame_id in range(N):
        store_frame(frame_id+1, td, conn, consumer, compressFrame=True)

if __name__ == "__main__":
    main()

