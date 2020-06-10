#!/usr/bin/python -u
"""
Script reads data from a set of stations and trim by event data.
"""

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool


# =====================================
# Importing trim data script_py 
# =====================================

from visual_py.event_plot import cut_data_by_event,plot_event_data,plot_event_dataset


# ==================================================
#  Importing some parameters from configuration file 
# ==================================================

from parameters_py.config import (
										DIR_DATA,TAUPY_MODEL,EV_GCARC_MIN,EV_GCARC_MAX,OUTPUT_JSON_FILE_DIR,MP_PROCESSES,OUTPUT_EV_DIR
				   )


# ============================================
#  Importing station dictionary from JSON file 
# ============================================

print('\n')
print('Looking for STATIONS data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'STA_dic.json'

sta_dic = json.load(open(filename_STA))

kstnm = sta_dic['KSTNM']
stla = sta_dic['STLA']
stlo = sta_dic['STLO']

for i in kstnm:
	print('Station = '+i)
print('\n')

# ============================================
#  Importing Event dictionary from JSON file 
# ============================================


print('\n')
print('Looking for Events data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'EVENT_dic.json'

event_dic = json.load(open(filename_STA))

ev_year = event_dic['ev_year']
ev_month = event_dic['ev_month']
ev_julday = event_dic['ev_julday']
ev_day = event_dic['ev_day']
ev_hour = event_dic['ev_hour']
ev_minute = event_dic['ev_minute']
ev_second = event_dic['ev_second']
ev_microsecond = event_dic['ev_microsecond']
ev_timeUTC = event_dic['ev_timeUTC']
evla = event_dic['evla']
evlo = event_dic['evlo']
evdp = event_dic['evdp']
mag = event_dic['mag']

print('Number of events = '+str(len(mag)))
print('\n')

# ==============================
#  Creating stations Input lists
# ==============================

print('Creating stations input lists')
print('\n')

input_list = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	print('Creating input list: '+j)
	print('\n')
	input_list[i] = [
			[kstnm[i],stla[i],stlo[i],ev_timeUTC[k],ev_julday[k],ev_year[k],ev_month[k],ev_day[k],ev_hour[k],ev_minute[k],ev_second[k],ev_microsecond[k],evla[k],evlo[k],evdp[k],mag[k]]
			 for k,l in enumerate(ev_year)
			]
print('\n')

# ================
# trim event data
# ================
'''
print('Events for each station')
print('\n')

for i,j in enumerate(input_list):
	print('Station: '+kstnm[i])
	[cut_data_by_event(kstnm=kstnm[i],stla=stla[i],stlo=stlo[i],ev_timeUTC=ev_timeUTC[k],ev_julday=ev_julday[k],ev_year=ev_year[k],ev_month=ev_month[k],
						ev_day=evla[k],ev_hour=ev_hour[k],ev_minute=ev_minute[k],ev_second=ev_second[k],ev_microsecond=ev_microsecond[k],
						ev_lat=evla[k],ev_long=evlo[k],ev_depth=evdp[k],ev_mag=mag[k]) 
						for k,l in enumerate(ev_year)]

print('Cutting finished!')


# ================
#  plot event data
# ================

datafile_lst = [] 
for root, dirs, files in os.walk(OUTPUT_EV_DIR):
	for directories in dirs:
		datafile_name = os.path.join(root, directories)
		datafile_lst.append(datafile_name)
datafile_lstS = sorted(datafile_lst)

for i,j in enumerate(datafile_lstS):
	plot_event_data(j)

'''
# ====================
#  plot event dataset
# ====================

datafile_lst = [] 
for root, dirs, files in os.walk(OUTPUT_EV_DIR):
	for directories in dirs:
		datafile_name = os.path.join(root, directories)
		datafile_lst.append(datafile_name)
datafile_lstS = sorted(datafile_lst)

for i,j in enumerate(datafile_lstS):
	plot_event_dataset(j)

print('Plotting finished!')
