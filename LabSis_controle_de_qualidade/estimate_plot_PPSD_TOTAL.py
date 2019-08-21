#!/usr/bin/python -u
"""
Script to estimate and save PPSD data of the raw data
"""
import os
import json
import obspy
import datetime
import random



# ==============================
# Generating DATA availability
# ==============================

from PSD_py.PSD_save import calc_PSD

from visual_py.plot_PSD_DATA import plot_PPSD_TOTAL_data



from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,OUTPUT_PSD_DIR,
				   )

# ========================================
# Importing data from raw files directory 
# ========================================

data_lista = []

for root, dirs, files in os.walk(DIR_DATA):
	for name in files:
		data_lista.append(os.path.join(root, name))

data_lista = sorted(data_lista)

# ==============
#  Get PPSD Data 
# ==============

print('Get Data for calculating PPSD for each station')
print('\n')

calc_PSD_result = calc_PSD(data_lista)

# ===========================
# Finding stations PPSD data
# ===========================

print('\n')
print('Looking for PPSD STATIONS data in'+OUTPUT_PSD_DIR)
print('\n')

datafile_lst = [] 
for root, dirs, files in os.walk(OUTPUT_PSD_DIR):
	for directories in dirs:
		datafile_name = os.path.join(root, directories)
		if '.PPSD' in datafile_name:
			datafile_lst.append(datafile_name)
datafile_lstS = sorted(datafile_lst)

# ================
#  plot PPSD Data 
# ================

for i,j in enumerate(datafile_lstS):
	plot_PPSD_TOTAL_data(j)
