'''
Script to collect information from stations

An example of STA_CSV_FILE is shown bellow:

knetwk;kstnm;lat;lon;elev
ON;ABR01;-17.9646;-38.6959;38.00
ON;ALF01;-20.6169;-40.7252;22.00
ON;ANA01;-14.5310;-40.7452;740.70
ON;CAJ01;-24.8501;-48.2433;372.00
ON;CAM01;-21.8257;-41.6574;31.00
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

sta_lat_lon = np.genfromtxt(STA_CSV_FILE,skip_header=1,usecols=[2,3,4],delimiter=';')

sta_name =  np.genfromtxt(STA_CSV_FILE,dtype='str',skip_header=1,usecols=[0,1],delimiter=';')

sta_event = {

		'KNETWK':[],
		'KSTNM':[],
		'STLA':[],
		'STLO':[],
		'STEL':[],
	    }

for i,j in enumerate(sta_name):
	sta_event['STLA'].append(sta_lat_lon[i][0])
	sta_event['STLO'].append(sta_lat_lon[i][1])
	sta_event['STEL'].append(sta_lat_lon[i][2])
	sta_event['KNETWK'].append(j[0])
	sta_event['KSTNM'].append(j[1])


print('Number of Stations: '+str(len(sta_event['KSTNM'])))
for i,j in enumerate(sta_event['KSTNM']):
	print('Station: '+sta_event['KNETWK'][i]+'.'+j)
	print('\n')


print('\n')

print('Saving Station Information in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'STA_dic.json', 'w') as fp:
	json.dump(sta_event, fp)