#!/usr/bin/python -u
"""
This script reads receiver function data from a set of stations, and
calculates the ray path between all pairs station-event. The data MUST be 
located in some root folder (RF_DIR). The station-event informations 
(coordinates) MUST be in the header of the files. In the current version 
of the program, files should be organized inside their directory as you prefer, 
but they need to finish with some extension name (RF_EXT), e.g.:
*.mseed, *.sac, *.itr, *.eqr
The implemented algorithm follows the lines of Gao and Liu (2014),
"Imaging mantle discontinuities using multiply-reflected P-to-S
conversions", Earth and Planetary Science Letters 402 (2014) 99–106.
Here we utilize both the P-to-S converted phase (Pds) and the multiply reflected 
and converted phase (Ppds) at the discontinuities to simultaneously determine the 
depth of mantle discontinuities and velocity anomalies in the overlying layer. 
Note that all the parameters are defined in the configuration file.
"""
import warnings
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool

# =====================================
# Importing piercing points scritp_py 
# =====================================

from piercing_points_py.piercing_points import arrivals_calculation

# ==================================================
# Importing some parameters from configuration file 
# ==================================================

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,STA_DIR,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					PP_FIGURE,EXT_FIG,DPI_FIG,MP_PROCESSES
				   )
# =====================================
# Function to estimate piercing points  
# =====================================

def parallel_piercing_points(number,PHASE,ev_depth,ev_lat,ev_long,st_lat,st_long):
	PP_dic = {'dist':[],'depth':[],'time':[],'lat':[],'lon':[],'number':[]}
	arrivals_calculation_lst = arrivals_calculation(fase=PHASE,ev_depth=ev_depth,ev_lat=ev_lat,ev_long=ev_long,st_lat=st_lat,st_long=st_long)

	PP_dic['time'].append(arrivals_calculation_lst[0])
	PP_dic['depth'].append(arrivals_calculation_lst[1])
	PP_dic['dist'].append(arrivals_calculation_lst[2])
	PP_dic['lat'].append(arrivals_calculation_lst[3])
	PP_dic['lon'].append(arrivals_calculation_lst[4])
	print(number)
	PP_dic['number'].append(number)

	return PP_dic

# ============================================
# Importing station dictionary from JSON file 
# ============================================

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

# ==================================================
# Importing earth model from obspy.taup.TauPyModel 
# ==================================================

print('Importing earth model from obspy.taup.TauPyModel')
print('\n')
model_THICKNESS_km = TauPyModel(model=MODEL_FILE_NPZ)


# =========================
# Creating the earth layers
# =========================

print('Creating the earth layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

# =========================
# Creating the middle DEPTH
# =========================

print('Creating the middle DEPTH')
print('\n')

dist_med_camada_terra = [abs(c - ((410+660)/2)) for x,c in enumerate(camadas_terra_10_km)]

DEPTH_MED = camadas_terra_10_km[dist_med_camada_terra.index(min(dist_med_camada_terra))]

# ====================
# Creating Pds list
# ====================

print('Creating Pds list')
print('\n')

PHASES = 'P410s','P'+"{0:.0f}".format(DEPTH_MED)+'s','P660s'

# =========================
# Creating Pds Input lists
# ========================

print('Creating Pds Input lists')
print('\n')


input_list_410 = [[i,PHASES[0],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i]] for i,j in enumerate(event_depth)]

input_list_DEPTH_MED = [[i,PHASES[1],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i]] for i,j in enumerate(event_depth)]

input_list_660 = [[i,PHASES[2],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i]] for i,j in enumerate(event_depth)]

# ================================
# Calculating Pds Piercing Points
# ================================

print('Calculating Piercing Points to '+PHASES[0])
print('\n')

pool_410 = Pool(MP_PROCESSES)
PP_dic_410 = pool_410.starmap(parallel_piercing_points, input_list_410)
pool_410.close()
#Saving Piercing Points in JSON file
print('Saving Piercing Points in JSON file')
os.makedirs(PP_DIR,exist_ok=True)
with open(PP_DIR+'PP_'+PHASES[0]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_410, fp)

print('Piercing Points to '+PHASES[0]+' estimated!')
print('\n')


# ============================
print('Calculating Piercing Points to '+PHASES[1])
print('\n')

pool_DEPTH_MED = Pool(MP_PROCESSES)
PP_dic_DEPTH_MED = pool_DEPTH_MED.starmap(parallel_piercing_points, input_list_DEPTH_MED)
pool_DEPTH_MED.close()
#Saving Piercing Points in JSON file
print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES[1]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_DEPTH_MED, fp)

print('Piercing Points to '+PHASES[1]+' estimated!')
print('\n')


# ============================

print('Calculating Piercing Points to '+PHASES[2])
print('\n')

pool_660 = Pool(MP_PROCESSES)
PP_dic_660 = pool_660.starmap(parallel_piercing_points, input_list_660)
pool_660.close()
#Saving Piercing Points in JSON file
print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES[2]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_660, fp)

print('Piercing Points to '+PHASES[2]+' estimated!')
print('\n')

# ============================================================================================================================


# ====================
# Creating Ppds list
# ====================

print('Creating Ppds list')
print('\n')

PHASES_Ppds = 'PPv410s','PPv'+"{0:.0f}".format(DEPTH_MED)+'s','PPv660s'

# =========================
# Creating Pds Input lists
# ========================

print('Creating Pds Input lists')
print('\n')


input_list_410_Ppds = [[i,PHASES_Ppds[0],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i]] for i,j in enumerate(event_depth)]

input_list_DEPTH_MED_Ppds = [[i,PHASES_Ppds[1],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i]] for i,j in enumerate(event_depth)]

input_list_660_Ppds = [[i,PHASES_Ppds[2],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i]] for i,j in enumerate(event_depth)]

# ================================
# Calculating Pds Piercing Points
# ================================

print('Calculating Piercing Points to '+PHASES_Ppds[0])
print('\n')

pool_410_Ppds = Pool(MP_PROCESSES)
PP_dic_410_Ppds = pool_410_Ppds.starmap(parallel_piercing_points, input_list_410_Ppds)
pool_410_Ppds.close()
#Saving Piercing Points in JSON file
print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES_Ppds[0]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_410_Ppds, fp)

print('Piercing Points to '+PHASES_Ppds[0]+' estimated!')
print('\n')


# ============================
print('Calculating Piercing Points to '+PHASES_Ppds[1])
print('\n')

pool_DEPTH_MED_Ppds = Pool(MP_PROCESSES)
PP_dic_DEPTH_MED_Ppds = pool_DEPTH_MED_Ppds.starmap(parallel_piercing_points, input_list_DEPTH_MED_Ppds)
pool_DEPTH_MED_Ppds.close()
#Saving Piercing Points in JSON file
print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES_Ppds[1]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_DEPTH_MED_Ppds, fp)

print('Piercing Points to '+PHASES_Ppds[1]+' estimated!')
print('\n')


# ============================

print('Calculating Piercing Points to '+PHASES[2])
print('\n')

pool_660_Ppds = Pool(MP_PROCESSES)
PP_dic_660_Ppds = pool_660_Ppds.starmap(parallel_piercing_points, input_list_660_Ppds)
pool_660_Ppds.close()
#Saving Piercing Points in JSON file
print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES_Ppds[2]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_660_Ppds, fp)

print('Piercing Points to '+PHASES_Ppds[2]+' estimated!')
print('\n')

# ============================================================================================================================