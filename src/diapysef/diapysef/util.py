#!/usr/bin/env python
from __future__ import print_function
import pyopenms

def setCompressionOptions(opt):
    """
    Adds suitable compression options for an object of type
    pyopenms.PeakFileOptions
        - compresses mass / time arrays with numpress linear
        - compresses intensity with slof (log integer)
        - compresses ion mobility with slof (log integer)
    """
    cfg = pyopenms.NumpressConfig()
    cfg.estimate_fixed_point = True
    cfg.numpressErrorTolerance = -1.0 # skip check, faster
    cfg.setCompression(b"linear");
    cfg.linear_fp_mass_acc = -1; # set the desired RT accuracy in seconds
    opt.setNumpressConfigurationMassTime(cfg)
    cfg = pyopenms.NumpressConfig()
    cfg.estimate_fixed_point = True
    cfg.numpressErrorTolerance = -1.0 # skip check, faster
    cfg.setCompression(b"slof");
    opt.setNumpressConfigurationIntensity(cfg)
    opt.setCompression(True) # zlib compression

    # Now also try to compress float data arrays (this is not enabled in all
    # versions of pyOpenMS).
    try:
        cfg = pyopenms.NumpressConfig()
        cfg.estimate_fixed_point = True
        cfg.numpressErrorTolerance = -1.0 # skip check, faster
        cfg.setCompression(b"slof");
        opt.setNumpressConfigurationFloatDataArray(cfg)
    except Exception:
        pass

