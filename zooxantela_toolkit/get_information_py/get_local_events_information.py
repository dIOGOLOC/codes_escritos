'''
Script to get local events informations from some text file:

An example of LOCAL_CSV_FILE is shown bellow:

 YEAR MMDD HHMMSS  LAT. LONG.  ERR ST DEPTH MAG. T CAT Io  AREA LOCALITY   COMMENTS 
 2001 0107 035015  -17.70 -44.70  10 MG   0.  3.4  1  I   -       Pirapora     (UnB)
 2001 0123 092131  -05.28 -39.42  50 CE   0.  3.3  1  I   -       Quixeramobim (IAG,UFRN)
 2001 0221 152021  -11.28 -74.51  10 PU  33.  5.5  0  I   2       Central Peru (IRIS)AC-IIMM
 2001 0226 204200  -04.41 -38.29  05 CE   0.  3.7  1  I   -       Cascavel     (IAG,UnB)
'''

import numpy as np
from obspy import UTCDateTime
import os
import json
from obspy import read_events
from obspy.clients.fdsn import Client
from datetime import datetime
import shapefile
from shapely.geometry import Point, shape


from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,LOCAL_CSV_FILE,SHP_AREA_DELIMITER,
					LOCAL_EVENT_START_DATE,LOCAL_EVENT_FINAL_DATE,LOCAL_EV_MAGNITUDE_MIN
				   )



print('Get Local Events Parameters')
print('\n')

local_event_info_txt = np.genfromtxt(LOCAL_CSV_FILE,dtype='str',skip_header=1,usecols=[0,1,2,3,4,7,8])

dic_local_event = {
		'ev_timeUTC':[],
		'evla':[],
		'evlo':[],
		'evdp':[],
		'mag':[]}

for i,j in enumerate(local_event_info_txt):
	try:
		time_str = j[0]+'-'+j[1][:2]+'-'+j[1][2:]+'T'+j[2][:2]+':'+j[2][2:4]+':'+j[2][4:]
		time = datetime.fromisoformat(j[0]+'-'+j[1][:2]+'-'+j[1][2:]+'T'+j[2][:2]+':'+j[2][2:4]+':'+j[2][4:])
		
		point = Point(float(j[4]), float(j[3])) # an x,y tuple
		shp = shapefile.Reader(SHP_AREA_DELIMITER) #open the shapefile
		if datetime.fromisoformat(LOCAL_EVENT_START_DATE) <= time <= datetime.fromisoformat(LOCAL_EVENT_FINAL_DATE) and point.within(shape(shp.shapes())) == True and float(j[6]) >= LOCAL_EV_MAGNITUDE_MIN:
			dic_local_event['ev_timeUTC'].append(time_str)
			print('Event '+str(i+1)+' - '+str(time_str))
			dic_local_event['evla'].append(float(j[3]))
			dic_local_event['evlo'].append(float(j[4]))
			dic_local_event['evdp'].append(float(j[5]))
			dic_local_event['mag'].append(float(j[6]))
	except:
		print('Error in Local event time '+' - line '+str(i+1)+' or Point out of the study area')


print('\n')
print('Number of Events: '+str(len(dic_local_event['mag'])))
print('\n')
		
print('Saving Event Parameters in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)

with open(OUTPUT_JSON_FILE_DIR+'LOCAL_EVENT_dic.json', 'w') as fp:
	json.dump(dic_local_event, fp)