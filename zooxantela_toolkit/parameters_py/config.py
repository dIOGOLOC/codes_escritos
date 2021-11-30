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

# directory of raw status files
DIR_STATUS = config.get('paths', 'DIR_STATUS')

#Directory to save JSON Files
OUTPUT_JSON_FILE_DIR =  config.get('paths', 'OUTPUT_JSON_FILE_DIR')

#Directory to save Figures
OUTPUT_FIGURE_DIR = config.get('paths', 'OUTPUT_FIGURE_DIR')

#Directory to save PSD
OUTPUT_PSD_DIR = config.get('paths', 'OUTPUT_PSD_DIR')

#Directory to save EVENT data
OUTPUT_EV_DIR = config.get('paths', 'OUTPUT_EV_DIR')

#XML File
XML_FILE = config.get('paths', 'XML_FILE')

#Stations CSV FILE path
STA_CSV_FILE  =  config.get('paths', 'STA_CSV_FILE')

#Local events FILE path
LOCAL_CSV_FILE  = config.get('paths', 'LOCAL_CSV_FILE')

#Shapefile  boundary states
BOUNDARY_STATES_SHP = config.get('paths', 'BOUNDARY_STATES_SHP')

# -----
# ppsd
# -----

#Number worker processes
NUM_PROCESS = config.getint('ppsd', 'NUM_PROCESS')

#PPSD start date
INITIAL_DATE = config.get('ppsd', 'INITIAL_DATE')

#PPSD final date
FINAL_DATE = config.get('ppsd', 'FINAL_DATE')

#Percentage fo the days to process and plot the PPSD?
DAY_PERCENTAGE = config.getint('ppsd', 'DAY_PERCENTAGE')

#Restricts the data that is included in the stack by time of day and weekday.
#Monday is 1, Sunday is 7, -1 for any day of week.
#For example, using time_of_weekday=[(-1, 22, 24)]
#only individual spectra that have a starttime in between 10pm and 12am are used in the stack for all days of week
#time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
TIME_OF_WEEKDAY_DAY = config.getint('ppsd', 'TIME_OF_WEEKDAY_DAY')
TIME_OF_WEEKDAY_START_HOUR = config.getint('ppsd', 'TIME_OF_WEEKDAY_START_HOUR')
TIME_OF_WEEKDAY_FINAL_HOUR = config.getint('ppsd', 'TIME_OF_WEEKDAY_FINAL_HOUR')

#Extracting PSD values for given period in seconds.
#Selects the period bin whose center period is closest to the specified period.
PERIOD_PSD = config.getfloat('ppsd', 'PERIOD_PSD')

#Maximum and minimum amplitude of the PSD.
AMP_PSD_MIN = config.getfloat('ppsd', 'AMP_PSD_MIN')
AMP_PSD_MAX = config.getfloat('ppsd', 'AMP_PSD_MAX')

#Maximum and minimum amplitude of the PSD (hydrophones).
AMP_PSD_HYDROPHONE_MIN = config.getfloat('ppsd', 'AMP_PSD_HYDROPHONE_MIN')
AMP_PSD_HYDROPHONE_MAX = config.getfloat('ppsd', 'AMP_PSD_HYDROPHONE_MAX')

# ----------
# local_evt
# ----------

#Local event start date
LOCAL_EVENT_START_DATE  = config.get('local_evt', 'LOCAL_EVENT_START_DATE')

#Local event final date
LOCAL_EVENT_FINAL_DATE  = config.get('local_evt', 'LOCAL_EVENT_FINAL_DATE')

#Minimum event magnitude
LOCAL_EV_MAGNITUDE_MIN = config.getfloat('local_evt', 'LOCAL_EV_MAGNITUDE_MIN')

#Minimum event distance
LOCAL_EV_DISTANCE_MIN = config.getfloat('local_evt', 'LOCAL_EV_DISTANCE_MIN')

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
