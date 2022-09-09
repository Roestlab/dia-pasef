#!/usr/bin/env python
from __future__ import print_function
import pyopenms
import sys
import click
import ast

# Logging and performance modules
from functools import wraps
import contextlib
import traceback
from time import time
import logging
from datetime import datetime

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

def method_timer(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        logging.debug('method:%r args:[%r, %r] took: %2.4f sec' % (
            f.__name__, args, kw, te-ts))
        return result
    return wrap


@contextlib.contextmanager
def code_block_timer(ident, log_type=logging.info):
    tstart = time()
    yield
    elapsed = time() - tstart
    log_type("{0}: Elapsed {1} ms".format(ident, elapsed))


def setup_logger(log_file, verbose):
    '''
    Setup logger
    '''
    if verbose == 0:
        use_verbosity = logging.INFO
    else:
        use_verbosity = logging.DEBUG

    root = logging.getLogger(__name__)
    root.setLevel(use_verbosity)
    logging.basicConfig(level=use_verbosity, filename=log_file, filemode='w',
                        format='[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(use_verbosity)
    root.addHandler(handler)


def check_sqlite_table(con, table):
    '''
    Check if a table exists in an SQL database
    '''
    table_present = False
    c = con.cursor()
    c.execute(
        'SELECT count(name) FROM sqlite_master WHERE type="table" AND name="%s"' % table)
    if c.fetchone()[0] == 1:
        table_present = True
    else:
        table_present = False
    c.fetchall()

    return (table_present)

def check_im_array(im_array):
    '''
    Check Ion Mobility Array to make sure values are valid
    '''
    
    if any(im_array<0):
      raise click.ClickException(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] ERROR: There are values below 0 in the input ion mobility array! {im_array[im_array<0]}. Most likely a pyopenms memory view issue.")
    elif 'e' in str(im_array):
      raise click.ClickException(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] ERROR: The input ion mobility array seems to be in scientific notation! {im_array}. Most likely a pyopenms memory view issue.")

def type_cast_value(value):
    '''
    Convert a string value to python literal
    '''
    if not isinstance(value, str):  # required for Click>=8.0.0
        return value
    try:
        return ast.literal_eval(value)
    except Exception:
        raise click.BadParameter(value)

def argument_value_log(args_dict):
    '''
    Print argument and value
    '''
    logging.debug("---------------- Input Parameters --------------------------")
    for key, item in args_dict.items():
        logging.debug(f"Parameter: {key} = {item}")
    logging.debug("------------------------------------------------------------\n")
