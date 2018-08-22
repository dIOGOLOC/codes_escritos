#!/usr/bin/python -u
"""
This script reads sac data from a set of stations and merge the data per day.
"""
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool


# =====================================
# Importing RF data script_py 
# =====================================

from visual_py.plot_raw_data import plot_station_raw_RF

# ==================================================
#  Importing some parameters from configuration file 
# ==================================================

from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,GAUSSIAN_FILTER
				   )


# ============================================
#  Importing station dictionary from JSON file 
# ============================================

print('\n')
print('Looking for STATIONS data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

JSON_FILE = 'RF_dic_'+str(GAUSSIAN_FILTER)

JSON_FILE_NAME = JSON_FILE.replace(".", "_")

filename_RF = OUTPUT_JSON_FILE_DIR+JSON_FILE_NAME+'.json'
dic_RF = json.load(open(filename_RF))

dataR = dic_RF['dataR']
dataT = dic_RF['dataT']

npts = dic_RF['npts']
kstnm = dic_RF['kstnm']
nzyear = dic_RF['nzyear']
nzjday = dic_RF['nzjday']
nzhour = dic_RF['nzhour']
nzmin = dic_RF['nzmin']
nzmsec = dic_RF['nzmsec']
evla = dic_RF['evla']
evlo = dic_RF['evlo']
evdp = dic_RF['evdp']
mag = dic_RF['mag']
stla = dic_RF['stla']
stlo = dic_RF['stlo']
user0 = dic_RF['user0']
user5 = dic_RF['user5']
user8 = dic_RF['user8']
dist = dic_RF['dist']
az = dic_RF['az']
baz = dic_RF['baz']
gcarc = dic_RF['gcarc']
b = dic_RF['b']
e = dic_RF['e']

filename_STA = OUTPUT_JSON_FILE_DIR+'STA_dic.json'

dic_STA = json.load(open(filename_STA))
kstnm_STA = dic_STA['kstnm']

# ==============================
#  Creating stations Input lists
# ==============================
input_list = [[]]*len(kstnm_STA)
print('Creating stations input lists')
print('\n')

for i,j in enumerate(kstnm_STA):
	print('Creating input list: '+j)
	print('\n')
	input_list[i] = [dataR[k] for k,l in enumerate(kstnm)  if l == j]
print('\n')

# ==============
#  Plotting data
# ==============

print('Plotting data')
print('\n')

for i,j in enumerate(input_list):
	print('Plotting data to: '+kstnm_STA[i])
	print('\n')
	plot_station_raw_RF(j,kstnm_STA[i])

print('Plotting finished!')

