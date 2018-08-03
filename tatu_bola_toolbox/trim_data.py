#!/usr/bin/python -u
"""
This script reads sac data from a set of stations and trim the data.
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

from pre_processing_py.trim_seismogram import cut_data_by_event

# ==================================================
#  Importing some parameters from configuration file 
# ==================================================

from parameters_py.config import (
					DIR_SAC,DIR_EVENT,NEIC_CSV_FILE,STA_CSV_FILE,OUTPUT_JSON_FILE_DIR,TAUPY_MODEL,
					EV_GCARC_MIN,EV_GCARC_MAX,EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,KCMPNM_N,
					KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,NAME_SUFFIX_Z,MP_PROCESSES,
					DETREND,DETREND_TYPE,TAPER,TAPER_TYPE,TAPER_MAX_PERCENTAGE,LOWPASS_FREQ,LOWPASS_CORNER,
					LOWPASS_ZEROPHASE,HIGHPASS_FREQ,HIGHPASS_CORNER,HIGHPASS_ZEROPHASE,RMEAN,RMEAN_TYPE,
					INTERPOLATE,SAMPLING_RATE
				   )



# ==================================
#  Function to call cut data script
# ==================================

def parallel_trim_data(kstnm,stla,stlo,ev_timeUTC,ev_julday,ev_year,ev_month,ev_day,ev_hour,ev_minute,ev_second,ev_microsecond,	ev_lat,ev_long,ev_depth,ev_mag):
			
	station_data_result = cut_data_by_event(kstnm=kstnm,stla=stla,stlo=stlo,ev_timeUTC=ev_timeUTC,ev_julday=ev_julday,ev_year=ev_year,ev_month=ev_month,
						ev_day=ev_day,ev_hour=ev_hour,ev_minute=ev_minute,ev_second=ev_second,ev_microsecond=ev_microsecond,
						ev_lat=ev_lat,ev_long=ev_long,ev_depth=ev_depth,ev_mag=ev_mag)

	return station_data_result


# ============================================
#  Importing station dictionary from JSON file 
# ============================================

print('\n')
print('Looking for STATIONS data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'STA_dic.json'

sta_dic = json.load(open(filename_STA))

kstnm = sta_dic['kstnm']
stla = sta_dic['stla']
stlo = sta_dic['stlo']

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

# ==============
#  Cutting data
# ==============

print('Cutting data and saving events folder for each station')
print('\n')

pool_trim = Pool(MP_PROCESSES)
error_station_data_dic = {'kstnm':[],'error_dir':[]}
for i,j in enumerate(input_list[8:9]):
	trim_data_result = pool_trim.starmap(parallel_trim_data, j)
	
	error_station_data_dic['error_dir'].append(pool_trim.starmap(parallel_trim_data, j))
	error_station_data_dic['kstnm'].append(kstnm[i])

pool_trim.close()


#Saving stations data error in JSON file
print('Saving stations data error in JSON file')
with open(OUTPUT_JSON_FILE_DIR+'STA_data_error.json', 'w') as fp:
	json.dump(error_station_data_dic, fp)

print('Cutting finished!')
