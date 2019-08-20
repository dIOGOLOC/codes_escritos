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
# PATH
# -----

#-------
#INPUTs

# directory of raw files
DIR_DATA = config.get('PATH', 'DIR_DATA')

# XML File
XML_FILE = config.get('PATH', 'XML_FILE')

#-------
#OUTPUTs

#Directory to save JSON Files
OUTPUT_JSON_FILE_DIR = config.get('PATH', 'OUTPUT_JSON_FILE_DIR')

#Directory to save Figures
OUTPUT_FIGURE_DIR = config.get('PATH', 'OUTPUT_FIGURE_DIR')

#Directory to save PSD
OUTPUT_PSD_DIR = config.get('PATH', 'OUTPUT_PSD_DIR')

# --------
# PPSD
# --------

#Initial date of the Data
INITIAL_DATE = config.get('PPSD', 'INITIAL_DATE')

#Final date of the data
FINAL_DATE = config.get('PPSD', 'FINAL_DATE')

#Restricts the data that is included in the stack by time of day and weekday. 
#Monday is 1, Sunday is 7, -1 for any day of week. 
#For example, using time_of_weekday=[(-1, 22, 24)] 
#only individual spectra that have a starttime in between 10pm and 12am are used in the stack for all days of week
#time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
TIME_OF_WEEKDAY_DAY = config.getint('PPSD', 'TIME_OF_WEEKDAY_DAY')
TIME_OF_WEEKDAY_START_HOUR = config.getfloat('PPSD', 'TIME_OF_WEEKDAY_START_HOUR')
TIME_OF_WEEKDAY_FINAL_HOUR = config.getfloat('PPSD', 'TIME_OF_WEEKDAY_FINAL_HOUR')

# --------
# CLIENT
# --------

#The user name is used for identification with the ArcLink server.
USER = config.get('CLIENT', 'USER')

#Host name of the remote ArcLink server (default host is 'webdc.eu').
HOST = config.get('CLIENT', 'HOST')

#Port of the remote ArcLink server (default port is 18002).
PORT = config.getint('CLIENT', 'PORT') 

#A string containing the name of the institution of the requesting person (default is an 'Anonymous').
INSTITUTION = config.get('CLIENT', 'HOST')