
# coding: utf-8

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool


from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,STA_DIR,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG
				   )

print('Importing earth model from obspy.taup.TauPyModel')
print('\n')
model_THICKNESS_km = TauPyModel(model=MODEL_FILE_NPZ)

# ==================
# Creating Function
# ==================

def arrivals_calculation(fase,ev_depth,ev_lat,ev_long,st_lat,st_long):
	print('Phase '+fase+' calculation')
	print('\n')
	arrivals = model_THICKNESS_km.get_pierce_points_geo(source_depth_in_km=ev_depth, source_latitude_in_deg=ev_lat, 
			                    source_longitude_in_deg=ev_long, receiver_latitude_in_deg=st_lat, 
			                    receiver_longitude_in_deg=st_long, phase_list=[fase])

	time = arrivals[0].pierce['time'].tolist()
	depth = arrivals[0].pierce['depth'].tolist()
	dist = arrivals[0].pierce['dist'].tolist()
	lat = arrivals[0].pierce['lat'].tolist()
	lon = arrivals[0].pierce['lon'].tolist()

	return [time,depth,dist,lat,lon]