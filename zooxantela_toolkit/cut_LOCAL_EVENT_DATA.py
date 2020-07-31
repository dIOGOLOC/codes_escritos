#!/usr/bin/python -u
"""
Scritp to download and cut local events data
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

from visual_py.local_event_download_plot import cut_data_by_local_event

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
for i,j in enumerate(ev_timeUTC):
	print(j)
print('\n')

# ================
# trim event data
# ================

print('Trim local event data for each station')
print('\n')

for i,j in enumerate(kstnm):
	print('Station: '+kstnm[i])
	[cut_data_by_local_event(knetwk=knetwk[i],kstnm=kstnm[i],stla=stla[i],stlo=stlo[i],evla=evla[k],evlo=evlo[k],evdp=evdp[k],evmag=mag[k],ev_timeUTC=ev_timeUTC[k]) 
						for k,l in enumerate(mag)]

print('Finished!')