# coding: utf-8

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from scipy.signal import triang

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG
				   )

print('Looking for receiver functions files in '+RF_DIR)
print('\n')
ev_list = []
ev_listS = []
for root, dirs, files in os.walk(RF_DIR):
    for datafile in files:
        if datafile.endswith(RF_EXT):
            ev_list.append(os.path.join(root, datafile))
ev_listS = sorted(ev_list)

ev = obspy.Stream()
for i,j in enumerate(ev_listS):
	ev += obspy.read(j)


print('Get Header Parameters')
print('\n')
sta_dic = {
	'event_depth':[],
	'event_lat':[],
	'event_long':[],
	'event_mag':[],
	'event_dist':[],
	'event_gcarc':[],
	'event_sta':[],
	'event_ray':[],
	'sta_lat':[],
	'sta_long':[],
	'sta_data':[],
	'sta_time':[]
	}

for i,j in enumerate(ev):
	if j.stats.sac.gcarc > 30:
		#check if the event depth is in km 
		if j.stats.sac.evdp > 1000:
			sta_dic['event_depth'].append(float(j.stats.sac.evdp/1000))
		else:
			sta_dic['event_depth'].append(float(j.stats.sac.evdp))
		sta_dic['event_mag'].append(float(j.stats.sac.mag))
		sta_dic['event_lat'].append(float(j.stats.sac.evla))
		sta_dic['event_long'].append(float(j.stats.sac.evlo))
		sta_dic['event_dist'].append(float(j.stats.sac.dist))
		sta_dic['event_gcarc'].append(float(j.stats.sac.gcarc))
		sta_dic['event_sta'].append(j.stats.station)
		sta_dic['event_ray'].append(float(j.stats.sac.user8))
		sta_dic['sta_lat'].append(float(j.stats.sac.stla))
		sta_dic['sta_long'].append(float(j.stats.sac.stlo))
		sta_dic['sta_data'].append(j.data[100:2700].tolist())
		sta_dic['sta_time'].append((j.times()[100:2700]-10).tolist())

print('Saving RF Header data in JSON file')
print('\n')

os.makedirs(STA_DIR,exist_ok=True)

with open(STA_DIR+'sta_dic.json', 'w') as fp:
	json.dump(sta_dic, fp)
