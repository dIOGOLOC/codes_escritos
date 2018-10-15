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

# directory of raw data
DIR_RAW_DATA = config.get('paths', 'DIR_RAW_DATA')

# directory of raw sac files
DIR_SAC = config.get('paths', 'DIR_SAC')

#Directory to save seismograms
DIR_EVENT = config.get('paths', 'DIR_EVENT')

#Directory to save seismograms
DIR_EVENT_NO_PP = config.get('paths', 'DIR_EVENT_NO_PP')

#Stations CSV FILE path
STA_CSV_FILE  =  config.get('paths', 'STA_CSV_FILE')

#Directory to save JSON Files
OUTPUT_JSON_FILE_DIR =  config.get('paths', 'OUTPUT_JSON_FILE_DIR')


# -----
# copy
# -----
#File format
FILE_TYPE = config.get('copy', 'FILE_TYPE')

#File suffix
FILE_SUFFIX = config.get('copy', 'FILE_SUFFIX')

#Do want to remove the mean of your data? (information in https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.detrend.html#obspy.core.trace.Trace.detrend)
FILTERS_RMEAN = config.getboolean('copy', 'FILTERS_RMEAN')
RMEAN_TYPE = config.get('copy', 'RMEAN_TYPE') 

#Do want to remove the trend of your data? (information in https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.detrend.html#obspy.core.trace.Trace.detrend)
FILTERS_DETREND = config.getboolean('copy', 'FILTERS_DETREND')
DETREND_TYPE = config.get('copy', 'DETREND_TYPE') 

#Do want to use taper in your data? (information in https://docs.obspy.org/master/packages/autogen/obspy.core.trace.Trace.taper.html)
FILTERS_TAPER = config.getboolean('copy', 'FILTERS_TAPER')
TAPER_TYPE = config.get('copy', 'TAPER_TYPE') 
TAPER_MAX_PERCENTAGE = config.getfloat('copy', 'TAPER_MAX_PERCENTAGE')

#Do want to use filter in your data? (help in https://docs.obspy.org/tutorial/code_snippets/filtering_seismograms.html)
FILTER = config.getboolean('copy', 'FILTER')
LOWPASS_FREQ = config.getfloat('copy', 'LOWPASS_FREQ')
LOWPASS_CORNER = config.getfloat('copy', 'LOWPASS_CORNER')
LOWPASS_ZEROPHASE = config.getboolean('copy', 'LOWPASS_ZEROPHASE')

HIGHPASS_FREQ = config.getfloat('copy', 'HIGHPASS_FREQ')
HIGHPASS_CORNER = config.getfloat('copy', 'HIGHPASS_CORNER')
HIGHPASS_ZEROPHASE = config.getboolean('copy', 'HIGHPASS_ZEROPHASE')


#Do want to interpolate your data?
#This operation is performed in place on the actual data arrays. The raw data is not accessible anymore afterwards. 
#Be careful when downsampling data and make sure to apply an appropriate anti-aliasing lowpass filter before interpolating in case it’s necessary.
INTERPOLATE = config.getboolean('copy', 'INTERPOLATE')
SAMPLING_RATE = config.getfloat('copy', 'SAMPLING_RATE')

# -----
# trim
# -----

# how many concurrent processes at the multiprocessing?
MP_PROCESSES = config.getint('trim', 'MP_PROCESSES') 

#Taup_time model to calculate travel times
TAUPY_MODEL = config.get('trim', 'TAUPY_MODEL') 

#Date event initial  (exemple: "2008-01-01")
INITIAL_DATE_EVENT = config.get('trim', 'INITIAL_DATE_EVENT') 

#Date event final  (exemple: "2008-01-01")
FINAL_DATE_EVENT = config.get('trim', 'FINAL_DATE_EVENT') 

#Minimum event distance 
EV_GCARC_MIN = config.getfloat('trim', 'EV_GCARC_MIN')

#Maximum event distance 
EV_GCARC_MAX = config.getfloat('trim', 'EV_GCARC_MAX')

#Minimum event magnitude 
EV_MAGNITUDE_MB = config.getfloat('trim', 'EV_MAGNITUDE_MB')

#Time to trim the seismogram before P-wave arrival
CUT_BEFORE_P = config.getfloat('trim', 'CUT_BEFORE_P')

#Time to trim the seismogram after P-wave arrival
CUT_AFTER_P = config.getfloat('trim', 'CUT_AFTER_P')

#Component name N
KCMPNM_N = config.get('trim', 'KCMPNM_N')

#Component name E
KCMPNM_E = config.get('trim', 'KCMPNM_E')

#Component name Z
KCMPNM_Z = config.get('trim', 'KCMPNM_Z')

#Network name
knetwk = config.get('trim', 'knetwk')

#Component N sac file name suffix
NAME_SUFFIX_N = config.get('trim', 'NAME_SUFFIX_N')

#Component E sac file name suffix
NAME_SUFFIX_E = config.get('trim', 'NAME_SUFFIX_E')

#Component Z sac file name suffix
NAME_SUFFIX_Z = config.get('trim', 'NAME_SUFFIX_Z')


# -----
# merge
# -----

#Do want to use filter in your data? (help in https://docs.obspy.org/tutorial/code_snippets/filtering_seismograms.html)
FILTERS = config.getboolean('merge', 'FILTERS')