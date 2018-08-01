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

# directory of raw sac files.
DIR_SAC = config.get('paths', 'DIR_SAC')

#Directory to save seismograms
DIR_EVENT = config.get('paths', 'DIR_EVENT')

#NEIC CSV FILE path
NEIC_CSV_FILE  = config.get('paths', 'NEIC_CSV_FILE')

#Stations CSV FILE path
STA_CSV_FILE  =  config.get('paths', 'STA_CSV_FILE')

#Directory to save JSON Files
OUTPUT_JSON_FILE_DIR =  config.get('paths', 'OUTPUT_JSON_FILE_DIR')

# -----
# trim
# -----

#Taup_time model to calculate travel times
TAUPY_MODEL = config.get('trim', 'TAUPY_MODEL') 

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
KCMPNM_E = config.get('trim', 'KCMPNM_N')

#Component name Z
KCMPNM_Z = config.get('trim', 'KCMPNM_N')

#Network name
knetwk = config.get('trim', 'KCMPNM_N')

#Component N sac file name suffix
NAME_SUFFIX_N = config.get('trim', 'NAME_SUFFIX_N')

#Component E sac file name suffix
NAME_SUFFIX_E = config.get('trim', 'NAME_SUFFIX_N')

#Component Z sac file name suffix
NAME_SUFFIX_Z = config.get('trim', 'NAME_SUFFIX_N')
