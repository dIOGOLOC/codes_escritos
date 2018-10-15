"""
Module that parses global parameters from a configuration
file at first import, to make them available to the other
parts of the program.
"""

import configparser
import os
import glob
import json


def select_and_parse_config_file(basedir='.', ext='cnf', verbose=True):
    """
    Reads a configuration file and returns an instance of ConfigParser:
    First, looks for files in *basedir* with extension *ext*.
    Asks user to select a file if several files are found,
    and parses it using ConfigParser module.
    @rtype: L{ConfigParser.ConfigParser}
    """
    config_files = glob.glob(os.path.join(basedir, u'*.{}'.format(ext)))


    if not config_files:
        raise Exception("No configuration file found!")

    if len(config_files) == 1:
        # only one configuration file
        config_file = config_files[0]
    else:
        print("Select a configuration file:")
        for i, f in enumerate(config_files, start=1):
            print("{} - {}".format(i, f))
        res = int(input(''))
        config_file = config_files[res - 1]

    if verbose:
        print("Reading configuration file: {}".format(config_file))

    conf = configparser.ConfigParser(allow_no_value=True)
    conf.read(config_file)

    return conf

# ==========================
# parsing configuration file
# ==========================

config = select_and_parse_config_file(basedir='.', ext='cnf', verbose=True)

# -----
# paths
# -----

# directory of raw sac files
DIR_SAC = config.get('paths', 'DIR_SAC')

# directory of selected RF files
DIR_SEL_SAC = config.get('paths', 'DIR_SEL_SAC')

#Directory to save JSON Files
OUTPUT_JSON_FILE_DIR =  config.get('paths', 'OUTPUT_JSON_FILE_DIR')

#Stations CSV FILE path
STA_CSV_FILE  =  config.get('paths', 'STA_CSV_FILE')

#Gaussian Filter
GAUSSIAN_FILTER = config.getfloat('paths', 'GAUSSIAN_FILTER')

#RADIAL RF EXTENSION
RADIAL_EXT = config.get('paths', 'RADIAL_EXT')

#TRANSVERSAL RF EXTENSION
TRANSVERSAL_EXT = config.get('paths', 'TRANSVERSAL_EXT')

# --------
# quality
# --------

#Trace check (minimium amplitude of trace)
CODA_TRACE_CHECK_AMP_MIN = config.getfloat('quality', 'CODA_TRACE_CHECK_AMP_MIN')

#Trace check (maximum amplitude of trace)
CODA_TRACE_CHECK_AMP_MAX = config.getfloat('quality', 'CODA_TRACE_CHECK_AMP_MAX')

#Trace check (minimium amplitude of trace)
ZERO_AMP_MIN = config.getfloat('quality', 'ZERO_AMP_MIN')

#Trace check (coda amplitude starts in (seconds))
CODA_TRACE_CHECK = config.getfloat('quality', 'CODA_TRACE_CHECK')

#Trace check (coda stardart deviation/anomalous amplitude ratio)
CODA_RATIO_AMP = config.getfloat('quality', 'CODA_RATIO_AMP')

#Percent recoveries of the observed radial component
RF_PERCENT = config.getfloat('quality', 'RF_PERCENT')

#Minimum event distance 
EV_GCARC_MIN = config.getfloat('quality', 'EV_GCARC_MIN')

#Maximum event distance 
EV_GCARC_MAX = config.getfloat('quality', 'EV_GCARC_MAX')

#Minimum event magnitude 
EV_MAGNITUDE_MB = config.getfloat('quality', 'EV_MAGNITUDE_MB')

#Time to trim the seismogram before P-wave arrival
CUT_BEFORE_P = config.getfloat('quality', 'CUT_BEFORE_P')

#Time to trim the seismogram after P-wave arrival
CUT_AFTER_P = config.getfloat('quality', 'CUT_AFTER_P')

#Sampling Rate of the traces
SAMPLING_RATE = config.getfloat('quality', 'SAMPLING_RATE')

# -----
# plot
# -----

#Plot X limits
X_LIM_MIN = config.getfloat('plot', 'X_LIM_MIN')

X_LIM_MAX = config.getfloat('plot', 'X_LIM_MAX')

