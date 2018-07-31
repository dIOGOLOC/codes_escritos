# coding: utf-8

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json


from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG
				   )

print('Looking for receiver functions data in JSON file in '+STA_DIR)

filename_STA = STA_DIR+'sta_dic.json'

sta_dic = json.load(open(filename_STA))

event_depth = sta_dic['event_depth']
event_lat = sta_dic['event_lat']
event_long = sta_dic['event_long']
event_dist = sta_dic['event_dist']
event_gcarc = sta_dic['event_gcarc']
event_sta = sta_dic['event_sta']
event_ray = sta_dic['event_ray']
sta_lat = sta_dic['sta_lat']
sta_long = sta_dic['sta_long']
sta_data = sta_dic['sta_data']
sta_time = sta_dic['sta_time']

print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
print('\n')
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)


print('Ppds conversion PHASES')
phase_lst = ['PPv'+str(i)+'s' for i in range(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)]
print(phase_lst)
print('\n')

#Creating dir for JSON file
os.makedirs(PdS_DIR,exist_ok=True)
#Creating dic for JSON file
PPvs_dic = {'dist':[],'depth':[],'time':[]}

dist_event = []
depth = []
time = []
for i,j in enumerate(event_depth):
	print('Event '+str(i+1)+' of '+str(len(event_depth))+' events')
	arrivalsP = model_10_km.get_travel_times_geo(source_depth_in_km=j, source_latitude_in_deg=event_lat[i], 
                                    source_longitude_in_deg=event_long[i], receiver_latitude_in_deg=sta_lat[i], 
                                    receiver_longitude_in_deg=sta_long[i], phase_list='P')

	arrivals = model_10_km.get_travel_times_geo(source_depth_in_km=j, source_latitude_in_deg=event_lat[i], 
		                            source_longitude_in_deg=event_long[i], receiver_latitude_in_deg=sta_lat[i], 
		                            receiver_longitude_in_deg=sta_long[i], phase_list=phase_lst)


	print('source_depth_in_km = '+str(j))
	print('source_latitude_in_deg = '+str(event_lat[i]))
	print('source_longitude_in_deg = '+str(event_long[i]))
	print('receiver_latitude_in_deg = '+str(sta_lat[i]))
	print('receiver_longitude_in_deg = '+str(sta_long[i]))
	
	print('P time = '+str(arrivalsP[0].time))
	print('PPds conversion time = '+str(arrivals[0].time))
	print('\n')
	#Saving times in JSON file
	PPvs_dic['time'].append(arrivals[0].time - arrivalsP[0].time)
	PPvs_dic['depth'].append(float(arrivals[0].name[3:-1]))
	PPvs_dic['dist'].append(arrivals[0].distance)


print('Saving Ppds times in JSON file')
#Saving JSON file
with open(PdS_DIR+'PPvs_dic.json', 'w') as fp:
	json.dump(PPvs_dic, fp)

print("Ppds estimated!")
