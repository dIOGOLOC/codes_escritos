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
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)

print('P-wave time calculation')
arrivalsP = []
for i,j in enumerate(event_depth):
	print('Event '+str(i)+' of '+str(len(event_depth))+' events')
	arrivalsP.append(model_10_km.get_travel_times_geo(
                                    source_depth_in_km=j, source_latitude_in_deg=event_lat[i], 
                                    source_longitude_in_deg=event_long[i], receiver_latitude_in_deg=sta_lat[i], 
                                    receiver_longitude_in_deg=sta_long[i], phase_list='P'))



print('PPvs conversion time calculation')
phase_lst = ['PPv'+str(i)+'s' for i in range(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)]
print(phase_lst)



arrivals = []
for i,j in enumerate(event_depth):
	print('Event '+str(i)+' of '+str(len(event_depth))+' events')
	arrivals.append(model_10_km.get_travel_times_geo(
		                            source_depth_in_km=j, source_latitude_in_deg=event_lat[i], 
		                            source_longitude_in_deg=event_long[i], receiver_latitude_in_deg=sta_lat[i], 
		                            receiver_longitude_in_deg=sta_long[i], phase_list=phase_lst))




dist_event = [[]]*len(arrivals)
depth = [[]]*len(arrivals)
time = [[]]*len(arrivals)

for i,j in enumerate(arrivals):
            time[i] = [l.time - arrivalsP[i][0].time for k,l in enumerate(j)]
            depth[i] = [float(l.name[3:-1]) for k,l in enumerate(j)]
            dist_event[i] = [l.distance for k,l in enumerate(j)]


#Saving times in JSON file
print('Saving PPvs times in JSON file')

os.makedirs(PdS_DIR,exist_ok=True)

PPvs_dic = {'dist':[],'depth':[],'time':[]}
for i,j in enumerate(dist_event):
	PPvs_dic['dist'].append(j)
	PPvs_dic['depth'].append(depth[i])
	PPvs_dic['time'].append(time[i])

with open(PdS_DIR+'PPvs_dic.json', 'w') as fp:
	json.dump(PPvs_dic, fp)

print("PPvs estimated!")

