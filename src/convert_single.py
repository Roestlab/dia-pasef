# -*- coding: utf-8 -*-
"""Conversion program to convert a Bruker TIMS file to a single mzML

"""

import sys

if len(sys.argv) < 3:
    raise RuntimeError("need arguments: tdf_directory output.mzML")

analysis_dir = sys.argv[1]
output_fname = sys.argv[2]

if sys.version_info.major == 2:
    analysis_dir = unicode(analysis_dir)

import timsdata, sqlite3, sys, time
import pyopenms
import numpy as np, matplotlib.pyplot as plt

def store_frame(frame_id, filename, q, exp, verbose=False):
    """
    Store a single frame as an individual mzML file

    Note that this is the easiest way to visualize and process the data but
    involves a few hacks, namely storing the IM axis as the RT of each
    spectrum. 
    """
    # Get a projected mass spectrum:
    q = conn.execute("SELECT NumScans, Time, Polarity  FROM Frames WHERE Id={0}".format(frame_id))
    tmp = q.fetchone()
    num_scans = tmp[0]
    time = tmp[1]
    pol = tmp[2]

    # Get the mapping of the ion mobility axis
    scan_number_axis = np.arange(num_scans, dtype=np.float64)
    ook0_axis = td.scanNumToOneOverK0(frame_id, scan_number_axis)

    # Traverse in reversed order to get low ion mobilities first
    for k, scan in reversed(list(enumerate(td.readScans(frame_id, 0, num_scans)))):
        index = np.array(scan[0], dtype=np.float64)
        mz = td.indexToMz(frame_id, index)
        intens = scan[1]

        # Store data in OpenMS Spectrum file -> each TOF push is an individual
        # spectrum and we store the ion mobility in the precursor. The frame
        # can be reconstructed by grouping all spectra with the same RT.
        s = pyopenms.MSSpectrum()
        s.set_peaks( (mz, intens) ) 
        s.setRT(time)
        p = pyopenms.Precursor()
        p.setDriftTime( ook0_axis [ k ] )
        s.setPrecursors([p])
        exp.consumeSpectrum(s)

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
if output_fname.lower().endswith("sqmass"):
    consumer = pyopenms.MSDataSqlConsumer(output_fname)

for frame_id in range(N):
    store_frame(frame_id+1, "null", q, consumer)


