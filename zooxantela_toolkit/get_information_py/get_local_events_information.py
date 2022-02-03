'''
--------------------------------------------------------------------------------
            Function to collect information of a local events
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 02/2022


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
Given a starttime and endtime, this code will return a JSON file with a list of
events downloaded from Preliminary Seismic Boletim by Centro de Sismologia da
USP that are located inside a determined area given by a shapefile.

More information in:
http://moho.iag.usp.br/eq/latest

Inputs:
LOCAL_CSV_FILE


An example of LOCAL_CSV_FILE download in:
http://moho.iag.usp.br/eq/latest

evid;origin;longitude;latitude;depth;magnitude;magnitudet;region;author;mode
usp2020mums;2020-06-30T23:35:40.31Z;-66.721;-23.861;203.5;4.1;mb;"Jujuy Province, Argentina";jroberto;M
usp2020muij;2020-06-30T21:21:44.662Z;-69.448;-25.673;10.0;4.7;mb;"Northern Chile";cleusa;M
usp2020mtqe;2020-06-30T12:10:54.829Z;-66.554;-23.669;208.1;4.5;mb;"Jujuy Province, Argentina";cleusa;M
usp2020mtgy;2020-06-30T07:31:46.66Z;-66.742;-23.558;231.0;3.9;mb;"Jujuy Province, Argentina";cleusa;M
usp2020mskk;2020-06-29T20:07:40.43Z;-40.434;-3.896;0.0;1.9;MLv;"Groairas/CE";jroberto;M
usp2020mrlp;2020-06-29T07:36:22.4Z;-40.199;-3.409;0.0;1.4;mR;"Santana do Acaraú/CE";jroberto;M
usp2020mrkj;2020-06-29T06:59:14.529Z;-40.278;-3.350;0.0;1.3;mR;"Santana do Acaraú/CE";jroberto;M
usp2020mrhl;2020-06-29T05:29:47.46Z;-63.906;-18.941;78.3;3.9;mb;"Central Bolivia";jroberto;M
usp2020mrgr;2020-06-29T05:06:31.226Z;-63.787;-18.674;10.0;4.6;mb;"Central Bolivia";jroberto;M
usp2020mrbx;2020-06-29T02:42:26.453Z;-71.789;-15.723;124.6;4.6;mb;"Southern Peru";cleusa;M

Data description:
    evid: event name
    origin: year-month-dayThour:minute:second of the event
    longitude: latitude of the event
    latitude: longitude of the event
    depth: depth of the event
    magnitude: magnitude of the event

    see http://moho.iag.usp.br/eq/bulletin and http://moho.iag.usp.br/eq/latest
    for more informations
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

local_event_info_txt = np.genfromtxt(LOCAL_CSV_FILE,dtype='str',skip_header=1,usecols=[1,2,3,4,5],delimiter=";")
dic_local_event = {
		'ev_timeUTC':[],
		'evla':[],
		'evlo':[],
		'evdp':[],
		'mag':[]}

for i,j in enumerate(local_event_info_txt):
    try:
        time_str =  j[0]
        time = datetime.fromisoformat(UTCDateTime(time_str).isoformat())

        #check if the event (point[x,y]) is inside the area (shapefile)
        point = Point(float(j[1]), float(j[2])) # an x,y tuple
        shp = shapefile.Reader(SHP_AREA_DELIMITER) #open the shapefile
        polygon = shape(shp.shapeRecords()[0].shape.__geo_interface__) # 1 polygon

        if datetime.fromisoformat(LOCAL_EVENT_START_DATE) <= time <= datetime.fromisoformat(LOCAL_EVENT_FINAL_DATE) and polygon.contains(point) == True and float(j[4]) >= LOCAL_EV_MAGNITUDE_MIN:
            dic_local_event['ev_timeUTC'].append(time_str)
            print('Event '+str(i+1)+' - '+str(time_str))
            dic_local_event['evla'].append(float(j[2]))
            dic_local_event['evlo'].append(float(j[1]))
            dic_local_event['evdp'].append(float(j[3]))
            dic_local_event['mag'].append(float(j[4]))

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
