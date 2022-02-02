#!/usr/bin/python -u
'''
--------------------------------------------------------------------------------
   Downloading or trimming local the dataset according to local events time
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2021


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
This code will download from a data center (f:cut_download_data_by_event) or
trim the local dataset (f:cut_data_by_local_event) according to a given an event
time and a list of stations.

More information in:
https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.trim.html
or
https://docs.obspy.org/packages/obspy.clients.fdsn.html


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
Event traces (format:SAC) and Images (format:PNG)


Examples of Usage (in command line):
   >> python download_cut_LOCAL_EVENT_DATA.py
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

from visual_py.local_event_download_plot import cut_download_data_by_event

# ==================================================
#  Importing some parameters from configuration file
# ==================================================

from parameters_py.config import (OUTPUT_JSON_FILE_DIR,
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
stel = sta_dic['STEL']

for i in kstnm:
	print('Station = '+i)
print('\n')

# ==================================================
#  Importing Local Event dictionary from JSON file
# ==================================================


print('\n')
print('Looking for Events data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'LOCAL_EVENT_dic.json'

event_dic = json.load(open(filename_STA))

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

print('Trim local event data for each station')
print('\n')

for i,j in enumerate(kstnm):
	print('Station: '+kstnm[i])
	[cut_download_data_by_event(knetwk=knetwk[i],kstnm=kstnm[i],stla=stla[i],stlo=stlo[i],evla=evla[k],evlo=evlo[k],evdp=evdp[k],evmag=mag[k],ev_timeUTC=ev_timeUTC[k])
						for k,l in enumerate(mag)]

print('Finished!')
