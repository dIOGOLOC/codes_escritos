# coding: utf-8

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json


from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG
				   )


print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
print('\n')
model_THICKNESS_km = TauPyModel(model=MODEL_FILE_NPZ)

# ====================
# Creating Pds list
# ====================

print('Creating Pds list')
print('\n')

PHASES = ['P'+str(i)+'s' for i in range(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)]
print(PHASES)
print('\n')

# ======================================
# Function to estimate Pds travel times  
# ======================================


def travel_time_calculation_Pds(ev_depth,ev_lat,ev_long,st_lat,st_long):
		arrivalsP = model_THICKNESS_km.get_travel_times_geo(source_depth_in_km=ev_depth, source_latitude_in_deg=ev_lat, 
					            source_longitude_in_deg=ev_long, receiver_latitude_in_deg=st_lat, 
					            receiver_longitude_in_deg=st_long,phase_list='P')

		arrivals = model_THICKNESS_km.get_travel_times_geo(source_depth_in_km=ev_depth, source_latitude_in_deg=ev_lat, 
					            source_longitude_in_deg=ev_long, receiver_latitude_in_deg=st_lat, 
					            receiver_longitude_in_deg=st_long, phase_list=PHASES)


		
		time = []
		depth = []
		dist = []

		for k,l in enumerate(arrivals):
			print('source_depth_in_km = '+str(ev_depth))
			print('source_latitude_in_deg = '+str(ev_lat))
			print('source_longitude_in_deg = '+str(ev_long))
			print('receiver_latitude_in_deg = '+str(st_lat))
			print('receiver_longitude_in_deg = '+str(st_long))
			print('P time = '+str(arrivalsP[0].time))
			print('Ps conversion time = '+str(l.time))
			print('Pds - P time = '+str(l.time - arrivalsP[0].time))
			print('Depth = '+l.name[1:-1])
			print('\n')

			time.append(l.time - arrivalsP[0].time)
			depth.append(float(l.name[1:-1]))
			dist.append(l.distance)
			
		return [time,depth,dist]


print('Creating Ppds list')
print('\n')

PHASES_Ppds = ['PPv'+str(i)+'s' for i in range(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)]
print(PHASES_Ppds)
print('\n')


# ======================================
# Function to estimate Ppds travel times  
# ======================================

def travel_time_calculation_Ppds(ev_depth,ev_lat,ev_long,st_lat,st_long):
	arrivalsP = model_THICKNESS_km.get_travel_times_geo(source_depth_in_km=ev_depth, source_latitude_in_deg=ev_lat, 
			                    source_longitude_in_deg=ev_long, receiver_latitude_in_deg=st_lat, 
			                    receiver_longitude_in_deg=st_long,phase_list='P')

	arrivals = model_THICKNESS_km.get_travel_times_geo(source_depth_in_km=ev_depth, source_latitude_in_deg=ev_lat, 
			                    source_longitude_in_deg=ev_long, receiver_latitude_in_deg=st_lat, 
			                    receiver_longitude_in_deg=st_long, phase_list=PHASES_Ppds)

	time = []
	depth = []
	dist = []

	for k,l in enumerate(arrivals):
		print('source_depth_in_km = '+str(ev_depth))
		print('source_latitude_in_deg = '+str(ev_lat))
		print('source_longitude_in_deg = '+str(ev_long))
		print('receiver_latitude_in_deg = '+str(st_lat))
		print('receiver_longitude_in_deg = '+str(st_long))
		print('P time = '+str(arrivalsP[0].time))
		print('Ppds conversion time = '+str(l.time))
		print('Ppds - P time = '+str(l.time - arrivalsP[0].time))
		print('Depth = '+l.name[3:-1])
		print('\n')

		time.append(l.time - arrivalsP[0].time)
		depth.append(float(l.name[3:-1]))
		dist.append(l.distance)

	return [time,depth,dist]
