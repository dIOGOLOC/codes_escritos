''
# coding: utf-8

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
from multiprocessing import Pool
import pyarrow.feather as feather
import pandas as pd
import gc

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG,Ps_OR_Sp_PHASE
				   )

# ==================================================
# Importing earth model from obspy.taup.TauPyModel 
# ==================================================

print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
print('\n')

model_THICKNESS_km = TauPyModel(model=MODEL_FILE_NPZ)


# ===================
# Creating Phase list
# ===================

print('Creating '+Ps_OR_Sp_PHASE+' list')
print('\n')

if Ps_OR_Sp_PHASE == 'Ps':
	PHASES = ['P'+str(i)+'s' for i in range(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)]
	PHASES.insert(0,'P')
	s = ','

if Ps_OR_Sp_PHASE == 'Sp':
	PHASES = ['S'+str(i)+'p' for i in range(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)]
	PHASES.insert(0,'S')
	s = ','

print('\n')

# ==================
# Creating Function
# ==================

def get_pierce_points_geo_calculation(input):

	number = input[0]
	ev_depth = input[1]
	ev_lat = input[2]
	ev_long = input[3]
	st_lat = input[4]
	st_long = input[5]
	phase_folder = input[6]

	piercing_points = model_THICKNESS_km.get_pierce_points_geo(source_depth_in_km=ev_depth, 
														source_latitude_in_deg=ev_lat, 
														source_longitude_in_deg=ev_long, 
														receiver_latitude_in_deg=st_lat, 
														receiver_longitude_in_deg=st_long, 
														phase_list=PHASES)

	PP_dic = {
				'time':piercing_points[0].pierce['time'].tolist(),
				'depth':piercing_points[0].pierce['depth'].tolist(),
				'lat':piercing_points[0].pierce['lat'].tolist(),
				'lon':piercing_points[0].pierce['lon'].tolist()
			}

	PP_df = pd.DataFrame.from_dict(PP_dic)

	file_feather_name = phase_folder+'PP_dic'+str(number)+'.feather'
	feather.write_feather(PP_df, file_feather_name)

	#deleting references
	del PP_dic
	del PP_df

	#triggering collection
	gc.collect()