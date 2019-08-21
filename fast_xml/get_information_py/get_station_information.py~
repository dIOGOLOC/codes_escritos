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

NAME;LAT;LON;ELEV;SENSOR_KEYS;DATALOGGER_KEYS;ACCER_KEYS
PFBR;-4.3939;-41.4457;0;Geotech,KS2000,2000,120;Geotech,Smart24,1,20,100,100,Linear Phase;
NBMO;-3.3107;-40.0413;0;REF TEK,RT 151,A,120;REF TEK,RT 130 & 130-SMA,1,100;
NBPB;-5.5459;-39.5836;0;REF TEK,RT 151,A,120;REF TEK,RT 130 & 130-SMA,1,100;
NBPS;-4.3939;-41.4457;0;REF TEK,RT 151,A,120;REF TEK,RT 130 & 130-SMA,1,100;

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

sta_name =  np.genfromtxt(STA_CSV_FILE,dtype='str',skip_header=1,usecols=[0,4,5,6],delimiter=';')

sta_event = {
		'KSTNM':[],
		'STLA':[],
		'STLO':[],
		'STEL':[],
		'SENSOR_KEYS':[],
		'DATALOGGER_KEYS':[],
		'ACCER_KEYS':[]
	    }

for i,j in enumerate(sta_name):
	sta_event['STLA'].append(sta_lat_lon[i][0])
	sta_event['STLO'].append(sta_lat_lon[i][1])
	sta_event['STEL'].append(sta_lat_lon[i][2])
	sta_event['KSTNM'].append(j[0])
	sta_event['SENSOR_KEYS'].append(j[1])
	sta_event['DATALOGGER_KEYS'].append(j[2])
	sta_event['ACCER_KEYS'].append(j[3])

print('Number of Stations: '+str(len(sta_event['KSTNM'])))
for i,j in enumerate(sta_event['KSTNM']):
	print('Station: '+j)
	print('SENSOR_KEYS: '+sta_event['SENSOR_KEYS'][i])
	print('ACCER_KEYS: '+sta_event['ACCER_KEYS'][i])
	print('DATALOGGER_KEYS: '+sta_event['DATALOGGER_KEYS'][i])
	print('\n')


print('\n')

print('Saving Station Information in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'STA_dic.json', 'w') as fp:
	json.dump(sta_event, fp)

