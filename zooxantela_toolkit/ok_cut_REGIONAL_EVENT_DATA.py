#!/usr/bin/python -u
'''
--------------------------------------------------------------------------------
  Function to trim/plot the regional dataset according to Regional events time
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2021


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
This code will trim and plot the local datase according to a given an event time
and a list of stations.

More information in:
https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.trim.html


Inputs:
JSON file with event description:
    ev_timeUTC: event time in UTC (str)
    ev_year: year of the event
    ev_month: month of the event
    ev_day: day of the event
    ev_julday: julian day of the event
    ev_hour: hour of the event
    ev_minute: minute of the event
    ev_second: second of the event
    ev_microsecond: microsecond of the event
    evla: latitude of the event
    evlo: longitude of the event
    evdp: depth of the event
    mag: magnitude of the event

Outputs:
Event traces (format: SAC)

Examples of Usage (in command line):
   >> python cut_REGIONAL_EVENT_DATA.py

'''

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

from visual_py.event_plot import cut_data_by_event


# ==================================================
#  Importing some parameters from configuration file
# ==================================================

from parameters_py.config import (
								 DIR_DATA,TAUPY_MODEL,EV_GCARC_MIN,EV_GCARC_MAX,OUTPUT_JSON_FILE_DIR,OUTPUT_EV_DIR
								 )


# ============================================
#  Importing station dictionary from JSON file
# ============================================

print('\n')
print('Looking for STATIONS data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'STA_dic.json'

sta_dic = json.load(open(filename_STA))

knetwk = sta_dic['KNETWK']
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

# ================
# trim event data
# ================

print('Events for each station')
print('\n')

for i,j in enumerate(kstnm):
	print('Station: '+kstnm[i])
	[cut_data_by_event(knetwk=knetwk[i],kstnm=kstnm[i],stla=stla[i],stlo=stlo[i],ev_timeUTC=ev_timeUTC[k],ev_julday=ev_julday[k],ev_year=ev_year[k],ev_month=ev_month[k],
						ev_day=evla[k],ev_hour=ev_hour[k],ev_minute=ev_minute[k],ev_second=ev_second[k],ev_microsecond=ev_microsecond[k],
						ev_lat=evla[k],ev_long=evlo[k],ev_depth=evdp[k],ev_mag=mag[k])
						for k,l in enumerate(ev_year)]

print('Cutting finished!')
