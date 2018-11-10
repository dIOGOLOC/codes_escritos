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

#receiver functions dir
RF_DIR = config.get('paths', 'RF_DIR')

# extension of receiver functions
RF_EXT = config.get('paths', 'RF_EXT') 

#Directory of the model of IASP91 model of 10 km thick
MODEL_FILE_NPZ = config.get('paths', 'MODEL_FILE_NPZ')

# dir of the stations data
STA_DIR = config.get('paths', 'STA_DIR')

# dir of the Pds results
PdS_DIR = config.get('paths', 'PdS_DIR')

# dir of the piercing points results
PP_DIR = config.get('paths', 'PP_DIR')

# dir of Ray Path Figures
RAY_PATH_FIGURE = config.get('paths', 'RAY_PATH_FIGURE')

# dir of piercing points Figures
PP_FIGURE = config.get('paths', 'PP_FIGURE')

# dir of selected piercing points txt
PP_SELEC_DIR = config.get('paths', 'PP_SELEC_DIR')

# ---------------
# maps parameters
# ---------------

# paths to shapefiles (boundaries, tectonic provinces and labels)
BOUNDARY_1_SHP = config.get('maps', 'BOUNDARY_1_SHP')
BOUNDARY_1_SHP_NAME = config.get('maps', 'BOUNDARY_1_SHP_NAME')

BOUNDARY_2_SHP = config.get('maps', 'BOUNDARY_2_SHP')
BOUNDARY_2_SHP_NAME = config.get('maps', 'BOUNDARY_2_SHP_NAME')

TECTO_SHP = config.get('maps', 'TECTO_SHP')
TECTO_SHP_name = config.get('maps', 'TECTO_SHP_name')

# Lat and Lon of large maps and small maps
# center (lat/lon) of the map

PROJECT_LAT = config.getfloat('maps', 'PROJECT_LAT')
PROJECT_LON = config.getfloat('maps', 'PROJECT_LON')

#lower and upper corner (lat/lon) of the large map

LLCRNRLON_LARGE= config.getfloat('maps', 'LLCRNRLON_LARGE')
LLCRNRLAT_LARGE= config.getfloat('maps', 'LLCRNRLAT_LARGE')
URCRNRLON_LARGE= config.getfloat('maps', 'URCRNRLON_LARGE')
URCRNRLAT_LARGE= config.getfloat('maps', 'URCRNRLAT_LARGE')

#lower and upper corner (lat/lon) of the small map

LLCRNRLON_SMALL= config.getfloat('maps', 'LLCRNRLON_SMALL')
LLCRNRLAT_SMALL= config.getfloat('maps', 'LLCRNRLAT_SMALL')
URCRNRLON_SMALL= config.getfloat('maps', 'URCRNRLON_SMALL')
URCRNRLAT_SMALL= config.getfloat('maps', 'URCRNRLAT_SMALL')

#figures extention and dpi
EXT_FIG= config.get('maps', 'EXT_FIG')
DPI_FIG= config.getfloat('maps', 'DPI_FIG')

#Date event initial  (exemple: 2008-01-01T00:00:00)
INITIAL_DATE_EVENT = config.get('maps', 'INITIAL_DATE_EVENT')

#Date event final  (exemple: 2008-01-01T00:00:00)
FINAL_DATE_EVENT = config.get('maps', 'FINAL_DATE_EVENT')

# ---------------
# time
# ---------------

# how many concurrent processes at the multiprocessing?
MP_PROCESSES = config.getint('time', 'MP_PROCESSES')

MIN_DEPTH = config.getint('time', 'MIN_DEPTH')           
MAX_DEPTH = config.getint('time', 'MAX_DEPTH')
INTER_DEPTH = config.getint('time', 'INTER_DEPTH')

# ---------------
# migration
# ---------------

#Rotate grid?
ROTATE_GRID = config.getboolean('migration', 'ROTATE_GRID') 

#Filter by shapefile the grid?
FILTER_BY_SHAPEFILE = config.getboolean('migration', 'FILTER_BY_SHAPEFILE') 

#Shapefile containing lines or polygons representing grid boundary
SHAPEFILE_GRID = config.get('migration', 'SHAPEFILE_GRID') 

#What is the grid rotation angle in degrees? (clockwise rotation)
ROTATE_ANGLE = config.getfloat('migration', 'ROTATE_ANGLE') 

#GRID POINTS multiplyer
GRID_PP_MULT = config.getfloat('migration', 'GRID_PP_MULT') 

# distance between grid points and piercing points for each depth (degree)
DIST_GRID_PP = config.getfloat('migration', 'DIST_GRID_PP') 

# mean distance between grid points and piercing points (degree)
DIST_GRID_PP_MED = config.getfloat('migration', 'DIST_GRID_PP_MED') 

# number of piercing points per bin 
NUMBER_PP_PER_BIN = config.getint('migration', 'NUMBER_PP_PER_BIN')  

#False or True for depth estimation
DEPTH_ESTIMATION = config.getboolean('migration', 'DEPTH_ESTIMATION')  

#False or True for depth estimation with bootstrap
BOOTSTRAP_DEPTH_ESTIMATION = config.getboolean('migration', 'BOOTSTRAP_DEPTH_ESTIMATION') 

#Number of interations to compute bootstrapping
BOOTSTRAP_INTERATOR = config.getint('migration', 'BOOTSTRAP_INTERATOR')  

#GAMMA Number
GAMMA = config.getfloat('migration', 'GAMMA') 

#Number that multiply the difference between the cross-section points
SECTION_NUM = config.getfloat('migration', 'SECTION_NUM') 

#Number that multiply the difference between the cross-section points
DEPTH_RANGE = config.getfloat('migration', 'DEPTH_RANGE') 