import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json


from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,MODEL_FILE_TAUP,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,					
					PP_FIGURE,EXT_FIG,DPI_FIG
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


def travel_time_calculation_Pds(number,ev_depth,ev_lat,ev_long,st_lat,st_long,JSON_FOLDER):

	arrivals = model_THICKNESS_km.get_travel_times_geo(source_depth_in_km=ev_depth, 
														source_latitude_in_deg=ev_lat, 
														source_longitude_in_deg=ev_long, 
														receiver_latitude_in_deg=st_lat, 
														receiver_longitude_in_deg=st_long, 
														phase_list=PHASES)

	print('source_depth_in_km = '+str(ev_depth))
	print('source_latitude_in_deg = '+str(ev_lat))
	print('source_longitude_in_deg = '+str(ev_long))
	print('receiver_latitude_in_deg = '+str(st_lat))
	print('receiver_longitude_in_deg = '+str(st_long))

	Pds_dic = {"arrivals": []}
	for i, j in enumerate(arrivals):
		phase_dic = {"phase": j.name,"time": j.time, "rayparam": j.ray_param}
		Pds_dic["arrivals"].append(phase_dic)
	
	#Saving Piercing Points in JSON file
	print('Saving Pds Travel Times in JSON file')
	print('\n')

	with open(JSON_FOLDER+'Pds_dic_'+str(number)+'.json', 'w') as fp:
		json.dump(Pds_dic, fp)

print('Creating Ppds list')
print('\n')

PHASES_Ppds = ['PPv'+str(i)+'s' for i in range(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)]
PHASES_Ppds.insert(0,'P')
s = ','
PHASES_Ppds_lst = s.join(PHASES_Ppds)
print(PHASES_Ppds_lst)
print('\n')


# ======================================
# Function to estimate Ppds travel times  
# ======================================

def travel_time_calculation_Ppds(number,ev_depth,ev_lat,ev_long,st_lat,st_long,JSON_FOLDER):

	arrivals = model_THICKNESS_km.get_travel_times_geo(source_depth_in_km=ev_depth, 
														source_latitude_in_deg=ev_lat, 
														source_longitude_in_deg=ev_long, 
														receiver_latitude_in_deg=st_lat, 
														receiver_longitude_in_deg=st_long, 
														phase_list=PHASES_Ppds)

	print('source_depth_in_km = '+str(ev_depth))
	print('source_latitude_in_deg = '+str(ev_lat))
	print('source_longitude_in_deg = '+str(ev_long))
	print('receiver_latitude_in_deg = '+str(st_lat))
	print('receiver_longitude_in_deg = '+str(st_long))

	Pds_dic = {"arrivals": []}
	for i, j in enumerate(arrivals):
		phase_dic = {"phase": j.name,"time": j.time, "rayparam": j.ray_param}
		Pds_dic["arrivals"].append(phase_dic)
	
	#Saving Piercing Points in JSON file
	print('Saving Pds Travel Times in JSON file')
	print('\n')

	with open(JSON_FOLDER+'Pds_dic_'+str(number)+'.json', 'w') as fp:
		json.dump(Pds_dic, fp)