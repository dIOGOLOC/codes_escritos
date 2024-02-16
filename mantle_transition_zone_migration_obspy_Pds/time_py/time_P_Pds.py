import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json


from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG,OUTPUT_DIR
				   )


print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
print('\n')
model_THICKNESS_km = TauPyModel(model=MODEL_FILE_NPZ)


print('Creating the Earth layered model')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)
print(camadas_terra_10_km)
print('\n')


# ====================
# Creating Pds list
# ====================

print('Creating Pds list')
print('\n')

PHASES = ['P'+str(i)+'s' for i in range(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)]
PHASES.insert(0,'P')
s = ','
PHASES_lst = s.join(PHASES)
print(PHASES_lst)
print('\n')


# ======================================
# Function to estimate Pds travel times  
# ======================================


def travel_time_calculation_Pds(input):

	number = input[0]
	ev_depth = input[1]
	ev_lat = input[2]
	ev_long = input[3]
	st_lat = input[4]
	st_long = input[5]
	JSON_FOLDER = input[6]
	
	arrivals = model_THICKNESS_km.get_travel_times_geo(source_depth_in_km=ev_depth, 
														source_latitude_in_deg=ev_lat, 
														source_longitude_in_deg=ev_long, 
														receiver_latitude_in_deg=st_lat, 
														receiver_longitude_in_deg=st_long, 
														phase_list=PHASES)

	Pds_dic = {"arrivals": []}
	for j in arrivals:
		phase_dic = {"phase": j.name,"time": j.time, "rayparam": j.ray_param,'ev_lat': ev_lat,'ev_long': ev_long,'st_lat': st_lat,'st_long': st_long}
		Pds_dic["arrivals"].append(phase_dic)

	with open(JSON_FOLDER+'Pds_dic_'+str(number)+'.json', 'w') as fp:
		json.dump(Pds_dic, fp)
