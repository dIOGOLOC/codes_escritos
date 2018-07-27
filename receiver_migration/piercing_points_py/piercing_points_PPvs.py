
# coding: utf-8

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,											
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,STA_DIR,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,BOUNDARY_1_SHP,DEPTH_1,DEPTH_2,
					BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG
				   )
print('\n')
print('Looking for receiver functions data in JSON file in '+STA_DIR)
print('\n')

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
print('\n')
model_THICKNESS_km = TauPyModel(model=MODEL_FILE_NPZ)

PHASES = 'PPv'+"{0:.0f}".format(DEPTH_1)+'s','PPv'+"{0:.0f}".format((DEPTH_1+DEPTH_2)/2)+'s','PPv'+"{0:.0f}".format(DEPTH_2)+'s'

for fase in PHASES:
	print('Phase '+fase+' calculation')
	print('\n')
	arrivalsP410s = []
	for i,j in enumerate(event_depth):
		print('Event '+str(i)+' of '+str(len(event_depth))+' events')
		arrivalsP410s.append(model_THICKNESS_km.get_pierce_points_geo(
		                            source_depth_in_km=j, source_latitude_in_deg=event_lat[i], 
		                            source_longitude_in_deg=event_long[i], receiver_latitude_in_deg=sta_lat[i], 
		                            receiver_longitude_in_deg=sta_long[i], phase_list=[fase]))


	dist_P410s = [[]]*len(arrivalsP410s)
	depth_P410s = [[]]*len(arrivalsP410s)
	time_P410s = [[]]*len(arrivalsP410s)
	lat_P410s = [[]]*len(arrivalsP410s)
	lon_P410s = [[]]*len(arrivalsP410s)

	for i,j in enumerate(arrivalsP410s):
		time_P410s[i] = [l.pierce['time'] for k,l in enumerate(j)][0].tolist()
		depth_P410s[i] = [l.pierce['depth'] for k,l in enumerate(j)][0].tolist()
		dist_P410s[i] = [l.pierce['dist'] for k,l in enumerate(j)][0].tolist()
		lat_P410s[i] = [l.pierce['lat'] for k,l in enumerate(j)][0].tolist()
		lon_P410s[i] = [l.pierce['lon'] for k,l in enumerate(j)][0].tolist()

	#Saving Piercing Points in JSON file
	print('Saving Piercing Points in JSON file')

	os.makedirs(PP_DIR,exist_ok=True)

	PP_dic = {'dist':[],'depth':[],'time':[],'lat':[],'lon':[]}
	for i,j in enumerate(dist_P410s):
		PP_dic['dist'].append(j)
		PP_dic['depth'].append(depth_P410s[i])
		PP_dic['time'].append(time_P410s[i])
		PP_dic['lat'].append(lat_P410s[i])
		PP_dic['lon'].append(lon_P410s[i])

	with open(PP_DIR+'PP_'+fase+'_dic.json', 'w') as fp:
		json.dump(PP_dic, fp)

	print('Piercing Points to '+fase+' estimated!')
	print('\n')
