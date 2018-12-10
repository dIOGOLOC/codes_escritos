
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
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,MODEL_FILE_TAUP,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,STA_DIR,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,PP_FIGURE,EXT_FIG,DPI_FIG
				   )

# ==================
# Creating Function
# ==================

def arrivals_calculation(number,fase,ev_depth,ev_lat,ev_long,st_lat,st_long,phase_folder):
	print('Phase '+fase+' calculation')
	print('\n')

	os.system('taup_pierce -mod '+MODEL_FILE_TAUP+' -h '+str(ev_depth)+' -ph '+fase+' -station '+str(st_lat)+' '+str(st_long)+' -evt '+str(ev_lat)+' '+str(ev_long)+' -o '+phase_folder+'PP_'+fase+'_dic_'+str(number))

	print('source_depth_in_km = '+str(ev_depth))
	print('source_latitude_in_deg = '+str(ev_lat))
	print('source_longitude_in_deg = '+str(ev_long))
	print('receiver_latitude_in_deg = '+str(st_lat))
	print('receiver_longitude_in_deg = '+str(st_long))
	print('Phase = '+fase)