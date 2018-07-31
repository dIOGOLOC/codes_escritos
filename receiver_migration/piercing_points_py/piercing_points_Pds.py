
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
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,DEPTH_1,DEPTH_2,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG
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

print('Creating the earth layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

dist_med_camada_terra = [abs(c - ((DEPTH_1+DEPTH_2)/2)) for x,c in enumerate(camadas_terra_10_km)]

DEPTH_MED = camadas_terra_10_km[dist_med_camada_terra.index(min(dist_med_camada_terra))]

PHASES = 'P'+"{0:.0f}".format(DEPTH_1)+'s','P'+"{0:.0f}".format(DEPTH_MED)+'s','P'+"{0:.0f}".format(DEPTH_2)+'s'

#Creating Piercing Points in JSON file
os.makedirs(PP_DIR,exist_ok=True)
PP_dic = {'dist':[],'depth':[],'time':[],'lat':[],'lon':[]}

for fase in PHASES:
	print('Phase '+fase+' calculation')
	print('\n')
	for i,j in enumerate(event_depth):
		print('Event '+str(i)+' of '+str(len(event_depth))+' events')
		arrivalsP410s = model_THICKNESS_km.get_pierce_points_geo(source_depth_in_km=j, source_latitude_in_deg=event_lat[i], 
		                            source_longitude_in_deg=event_long[i], receiver_latitude_in_deg=sta_lat[i], 
		                            receiver_longitude_in_deg=sta_long[i], phase_list=[fase])


		PP_dic['time'].append(arrivalsP410s[0].pierce['time'].tolist())
		PP_dic['depth'].append(arrivalsP410s[0].pierce['depth'].tolist())
		PP_dic['dist'].append(arrivalsP410s[0].pierce['dist'].tolist())
		PP_dic['lat'].append(arrivalsP410s[0].pierce['lat'].tolist())
		PP_dic['lon'].append(arrivalsP410s[0].pierce['lon'].tolist())

	#Saving Piercing Points in JSON file
	print('Saving Piercing Points in JSON file')

	with open(PP_DIR+'PP_'+fase+'_dic.json', 'w') as fp:
		json.dump(PP_dic, fp)

	print('Piercing Points to '+fase+' estimated!')
	print('\n')
