'''
--------------------------------------------------------------------------------
            Function to collect information of a local events
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2021


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
Given a starttime and endtime, this code will return a JSON file with a list of
events downloaded from Preliminary Seismic Boletim by Centro de Sismologia da
USP that are located inside a determined area given by a shapefile.

More information in:
http://moho.iag.usp.br/eq/bulletin/

Dataset in:
http://moho.iag.usp.br/boletim/boletim_txt/boletim2000.txt
and
http://moho.iag.usp.br/boletim/boletim_txt/boletim2001p.txt


Inputs:
LOCAL_CSV_FILE

An example of LOCAL_CSV_FILE download in:
http://moho.iag.usp.br/boletim/boletim_txt/boletim2001p.txt is shown bellow:

 YEAR MMDD HHMMSS  LAT. LONG.  ERR ST DEPTH MAG. T CAT Io  AREA LOCALITY   COMMENTS
 2001 0107 035015  -17.70 -44.70  10 MG   0.  3.4  1  I   -       Pirapora     (UnB)
 2001 0123 092131  -05.28 -39.42  50 CE   0.  3.3  1  I   -       Quixeramobim (IAG,UFRN)
 2001 0221 152021  -11.28 -74.51  10 PU  33.  5.5  0  I   2       Central Peru (IRIS)AC-IIMM
 2001 0226 204200  -04.41 -38.29  05 CE   0.  3.7  1  I   -       Cascavel     (IAG,UnB)


Data description:
    YEAR: year of the event
    MMDD: month and day of the event
    HHMMSS: hour,minute second of the event
    LAT: latitude of the event
    LONG: longitude of the event
    DEPTH: depth of the event
    MAG: magnitude of the event
    ...


Outputs:
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

        #check if the event (point[x,y]) is inside the area (shapefile)
        point = Point(float(j[4]), float(j[3])) # an x,y tuple
        shp = shapefile.Reader(SHP_AREA_DELIMITER) #open the shapefile
        polygon = shape(shp.shapeRecords()[0].shape.__geo_interface__) # 1 polygon

        if datetime.fromisoformat(LOCAL_EVENT_START_DATE) <= time <= datetime.fromisoformat(LOCAL_EVENT_FINAL_DATE) and polygon.contains(point) == True and float(j[6]) >= LOCAL_EV_MAGNITUDE_MIN:
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
