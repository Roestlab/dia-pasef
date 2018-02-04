# -*- coding: utf-8 -*-
"""Conversion program (testing only) to convert each frame in TIMS to a separate mzML

It maps the ion mobility dimension to the RT axis, for easy visualization in TOPPView.

DO NOT USE THIS IN PRODUCTION
"""

import sys

if len(sys.argv) < 2:
    raise RuntimeError("need arguments: tdf_directory")

analysis_dir = sys.argv[1]

if sys.version_info.major == 2:
    analysis_dir = unicode(analysis_dir)

import timsdata, sqlite3, sys, time
import pyopenms
import numpy as np, matplotlib.pyplot as plt

def store_frame(frame_id, filename, q, verbose=False):
    """
    Store a single frame as an individual mzML file

    Note that this is the easiest way to visualize and process the data but
    involves a few hacks, namely storing the IM axis as the RT of each
    spectrum. 
    """
    # Get a projected mass spectrum:
    q = conn.execute("SELECT NumScans FROM Frames WHERE Id={0}".format(frame_id))
    num_scans = q.fetchone()[0]

    # Get the mapping of the ion mobility axis
    scan_number_axis = np.arange(num_scans, dtype=np.float64)
    ook0_axis = td.scanNumToOneOverK0(frame_id, scan_number_axis)

    # Get a new MSExperiment from OpenMS
    e = pyopenms.MSExperiment()

    for k, scan in enumerate(td.readScans(frame_id, 0, num_scans)):
        index = np.array(scan[0], dtype=np.float64)
        mz = td.indexToMz(frame_id, index)
        intens = scan[1]

        # Store data in OpenMS Spectrum file -> each TOF push is an individual
        # spectrum and we set the "RT" to be the ion mobility dimension
        s = pyopenms.MSSpectrum()
        s.set_peaks( (mz, intens) ) 
        s.setRT( ook0_axis [ k ] )
        e.addSpectrum(s)

        if verbose: print "scan ", k
        if len(mz) > 0 and verbose:

            print "  len: ", len(mz), len(intens)
            print "  at ook0", ook0_axis[k]
            print "  len: ", mz, intens
            for p in s:
                print p.getMZ(), p.getIntensity()

    # Store file at designated position
    pyopenms.MzMLFile().store(filename, e)


td = timsdata.TimsData(analysis_dir)
conn = td.conn

# Get total frame count:
q = conn.execute("SELECT COUNT(*) FROM Frames")
row = q.fetchone()
N = row[0]
print("Analysis has {0} frames.".format(N))

# For testing
## store_frame(10, "/tmp/test2.mzML")

for frame_id in range(N):
    store_frame(frame_id+1, "/tmp/test_%s.mzML" % (frame_id+1), q)

