''
# coding: utf-8

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool


from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG
				   )

# ==================================================
# Importing earth model from obspy.taup.TauPyModel 
# ==================================================

print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
print('\n')

model_THICKNESS_km = TauPyModel(model=MODEL_FILE_NPZ)


# ==================
# Creating Function
# ==================

def arrivals_calculation(number,fase,ev_depth,ev_lat,ev_long,st_lat,st_long,phase_folder):
	print('Phase '+fase+' calculation')
	print('\n')

	piercing_points = model_THICKNESS_km.get_pierce_points_geo(source_depth_in_km=ev_depth, 
														source_latitude_in_deg=ev_lat, 
														source_longitude_in_deg=ev_long, 
														receiver_latitude_in_deg=st_lat, 
														receiver_longitude_in_deg=st_long, 
														phase_list=[fase])

	print('source_depth_in_km = '+str(ev_depth))
	print('source_latitude_in_deg = '+str(ev_lat))
	print('source_longitude_in_deg = '+str(ev_long))
	print('receiver_latitude_in_deg = '+str(st_lat))
	print('receiver_longitude_in_deg = '+str(st_long))
	print('Phase = '+fase)

	PP_dic = {
				'time':piercing_points[0].pierce['time'].tolist(),
				'depth':piercing_points[0].pierce['depth'].tolist(),
				'lat':piercing_points[0].pierce['lat'].tolist(),
				'lon':piercing_points[0].pierce['lon'].tolist()
			}

	print('Saving Piercing Points in JSON file')
	print('\n')

	with open(phase_folder+'PP_dic'+str(number)+'.json', 'w') as fp:
		json.dump(PP_dic, fp)