'''
Script to collect information from stations and create the XML file

The IRIS DMC Nominal Response Library

IRIS DMC began to collect an “authoritative” set of manufacturers’ 
recommended nominal instrument responses in SEED RESP format and
 publish these on the web at http://ds.iris.edu/NRL. 

The goal behind the Library is to make it easier for the seismological 
community to both share and create metadata for common instrumentation,
and to improve response accuracy for users of the data. Because the wide
range of features available in modern instrumentation can make locating 
the appropriate response a challenge, the Library organizes responses 
based on hardware and/or acquisition choices, elucidating these choices 
as needed. “Nominal” responses may also have limits to their validity, 
so the Library notes these as well.

Nominal Response Library (NRL) includes:
	o 146 sensor responses (representing 10 sensor manufacturers); and
	o 4705 datalogger responses (from 7 datalogger manufacturers).

An example of STA_CSV_FILE is shown bellow:

NAME;LAT;LON;ELEV;SENSOR_KEYS;DATALOGGER_KEYS;ACCER_KEYS;HYDROPHONE_KEYS
ACJC;-5.5843;-35.7861;0;Sprengnether (now Eentec)*S6000*640 and 213 Ohms;REF TEK*RT 130 & 130-SMA*1*100;
NBAN;-9.6687;-36.2749;0;REF TEK*RT 151*A*120;REF TEK*RT 130 & 130-SMA*1*100;REF TEK*RT 131 (also 130-SMA)*131A-02 (also 130-SMA)*SF1500S
NBBC;-5.5222;-45.2844;0;REF TEK*RT 151*A*120;REF TEK*RT 130 & 130-SMA*1*100;REF TEK*RT 131 (also 130-SMA)*131A-02 (also 130-SMA)*SF1500S

SENSOR_KEYS and DATALOGGER_KEYS are located in http://ds.iris.edu/NRL, and are used to create the XML file:
response = nrl.get_response( # doctest: +SKIP
    sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])

'''
#Importing modules

import numpy as np
import os
import json

from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,STA_CSV_FILE
				   )



print('Get Station Information')
print('\n')

sta_lat_lon = np.genfromtxt(STA_CSV_FILE,skip_header=1,usecols=[1,2,3],delimiter=';')

sta_name =  np.genfromtxt(STA_CSV_FILE,dtype='str',skip_header=1,usecols=[0,4,5,6,7],delimiter=';')

KSTNM = sta_name[0]
STLA = sta_lat_lon[0]
STLO = sta_lat_lon[1]
STEL = sta_lat_lon[2]

if sta_name[1]:
	SENSOR_KEYS = sta_name[1]
else:
	SENSOR_KEYS = []

if sta_name[2]:
	DATALOGGER_KEYS = sta_name[2]
else:
	DATALOGGER_KEYS = []

if sta_name[3]:
	ACCER_KEYS = sta_name[3]
else:
	ACCER_KEYS = []

if sta_name[4]:
	HYDROPHONE_KEYS = sta_name[4]
else:
	HYDROPHONE_KEYS = [] 


sta_event = {
		'KSTNM':KSTNM,
		'STLA':STLA,
		'STLO':STLO,
		'STEL':STEL,
		'SENSOR_KEYS':SENSOR_KEYS,
		'DATALOGGER_KEYS':DATALOGGER_KEYS,
		'ACCER_KEYS':ACCER_KEYS,
		'HYDROPHONE_KEYS':HYDROPHONE_KEYS
			}

print('\n')
print('Saving Station Information in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'STA_dic.json', 'w') as fp:
	json.dump(sta_event, fp)