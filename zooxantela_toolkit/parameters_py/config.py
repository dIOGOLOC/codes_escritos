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
DIR_DATA = config.get('paths', 'DIR_DATA')

#Directory to save JSON Files
OUTPUT_JSON_FILE_DIR =  config.get('paths', 'OUTPUT_JSON_FILE_DIR')

#Directory to save Figures
OUTPUT_FIGURE_DIR = config.get('paths', 'OUTPUT_FIGURE_DIR')

#Directory to save PSD
OUTPUT_PSD_DIR = config.get('paths', 'OUTPUT_PSD_DIR')

#Directory to save EVENT data
OUTPUT_EV_DIR = config.get('paths', 'OUTPUT_EV_DIR')

#Directory to save XML File
OUTPUT_XML_FILE_DIR = config.get('paths', 'OUTPUT_XML_FILE_DIR')

#Stations CSV FILE path
STA_CSV_FILE  =  config.get('paths', 'STA_CSV_FILE')

#Local events FILE path
LOCAL_CSV_FILE  = config.get('paths', 'LOCAL_CSV_FILE')

# ----------
# local_evt
# ----------

#Local event start date
LOCAL_EVENT_START_DATE  = config.get('local_evt', 'LOCAL_EVENT_START_DATE')

#Local event final date
LOCAL_EVENT_FINAL_DATE  = config.get('local_evt', 'LOCAL_EVENT_FINAL_DATE')

#Minimum event magnitude 
LOCAL_EV_MAGNITUDE_MIN = config.getfloat('local_evt', 'LOCAL_EV_MAGNITUDE_MIN')

#Shapefile to filter local event data
SHP_AREA_DELIMITER = config.get('local_evt', 'SHP_AREA_DELIMITER')

#Time to trim the seismogram before P-wave arrival
CUT_BEFORE_P_LOCAL = config.getfloat('local_evt', 'CUT_BEFORE_P_LOCAL')

#Time to trim the seismogram after P-wave arrival
CUT_AFTER_P_LOCAL = config.getfloat('local_evt', 'CUT_AFTER_P_LOCAL')

#Choose an available host and port to download data (https://docs.obspy.org/packages/obspy.clients.arclink.html):
USER = config.get('local_evt', 'USER')
HOST = config.get('local_evt', 'HOST')
PORT = config.get('local_evt', 'PORT')


# ------
# event
# ------

#Taup_time model to calculate travel times
TAUPY_MODEL = config.get('event', 'TAUPY_MODEL') 

#Minimum event distance 
EV_GCARC_MIN = config.getfloat('event', 'EV_GCARC_MIN')

#Maximum event distance 
EV_GCARC_MAX = config.getfloat('event', 'EV_GCARC_MAX')


#Minimum event magnitude 
EV_MAGNITUDE_MB = config.getfloat('event', 'EV_MAGNITUDE_MB')

#Date event initial  (exemple: 2008-01-01T00:00:00)
INITIAL_DATE_EVENT = config.get('event', 'INITIAL_DATE_EVENT')

#Date event final  (exemple: 2008-01-01T00:00:00)
FINAL_DATE_EVENT = config.get('event', 'FINAL_DATE_EVENT')

#Time to trim the seismogram before P-wave arrival
CUT_BEFORE_P = config.getfloat('event', 'CUT_BEFORE_P')

#Time to trim the seismogram after P-wave arrival
CUT_AFTER_P = config.getfloat('event', 'CUT_AFTER_P')