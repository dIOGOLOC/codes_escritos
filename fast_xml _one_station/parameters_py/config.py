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

#Directory to save JSON Files
OUTPUT_JSON_FILE_DIR =  config.get('PATH', 'OUTPUT_JSON_FILE_DIR')

#Directory to save XML File
OUTPUT_XML_FILE_DIR = config.get('PATH', 'OUTPUT_XML_FILE_DIR')

#Stations CSV FILE path
STA_CSV_FILE  =  config.get('PATH', 'STA_CSV_FILE')


# ---
# XML
# ---

#Name of the Network
NETWORK_CODE = config.get('XML', 'NETWORK_CODE')

#Name of the deployer
SOURCE = config.get('XML', 'SOURCE')

#Lotation code
LOCATION = config.get('XML', 'LOCATION')

#Description of the Network
NETWORK_DESCRIPTION = config.get('XML', 'NETWORK_DESCRIPTION')

#Start date of the Network
START_DATE = config.get('XML', 'START_DATE')

#Sampling Rate of the seismogram
SAMPLING_RATE = config.getfloat('XML', 'SAMPLING_RATE')