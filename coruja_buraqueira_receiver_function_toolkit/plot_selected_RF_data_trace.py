#!/usr/bin/python -u
"""
This script reads sac data from a set of stations and merge the data per day.
"""
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
from multiprocessing import Pool
import pandas as pd

# =====================================
# Importing RF data script_py 
# =====================================

from visual_py.plot_raw_data import plot_station_raw_RF_TRACE

# ==================================================
#  Importing some parameters from configuration file 
# ==================================================

from parameters_py.config import (
					OUTPUT_FEATHER_FILE_DIR,GAUSSIAN_FILTER
				   )


# ============================================
#  Importing station dictionary from JSON file 
# ============================================

print('\n')
print('Looking for STATIONS data in FEATHER file in '+OUTPUT_FEATHER_FILE_DIR)
print('\n')

FEATHER_FILE = 'RF_dic_'+str(GAUSSIAN_FILTER)

FEATHER_FILE_NAME = FEATHER_FILE.replace(".", "_")

filename_RF = OUTPUT_FEATHER_FILE_DIR+FEATHER_FILE_NAME+'.feather'
dic_RF = pd.read_feather(filename_RF)

filename_STA = OUTPUT_FEATHER_FILE_DIR+'STA_dic.feather'

dic_STA = pd.read_feather(filename_STA)
kstnm_STA = dic_STA['KSTNM'].values

# ==============
#  Plotting data
# ==============

print('Plotting data')
print('\n')

for i,j in enumerate(kstnm_STA):
	print('Plotting data to: '+kstnm_STA[i])
	print('\n')

	df_kstnm = dic_RF[dic_RF['kstnm'] == j]

	plot_station_raw_RF_TRACE(df_kstnm)

print('Plotting finished!')

