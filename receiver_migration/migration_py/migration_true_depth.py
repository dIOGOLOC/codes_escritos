
# coding: utf-8

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import copy
import matplotlib
from matplotlib.cm import get_cmap
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
import shapefile
from fatiando import gridder, utils
import scipy.io
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import json
import random
from matplotlib.colors import Normalize
from numpy import ma
from matplotlib import cbook




from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,PHASES_LST,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,DEPTH_1,DEPTH_2,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP_MED,PHASES_PPvs_LST,DIST_GRID_PP,
					LINEAR_STACKING,DEPTH_ESTIMATION,DEPTH_RANGE,BOOTSTRAP_INTERATOR,BOOTSTRAP_DEPTH_ESTIMATION,
					GAMMA
				   )

print('Starting Receiver Functions migration code to estimate true depths of the discontinuities')
print('\n')


print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)

for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[0] == DEPTH_1:
		Vp_depth_1 = j[2]
		Vs_depth_1 = j[4]
		
for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[0] == DEPTH_2:
		Vp_depth_2 = j[2]
		Vs_depth_2 = j[4]
		



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

print('Creating the earth layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

print('Importing piercing points to each PHASE')
print('\n')

PHASES = PHASES_LST.split(',')

print('Importing Piercing Points '+PHASES[0])
print('\n')

filename_1 = PP_DIR+'PP_'+PHASES[0]+'_dic.json'

PP_1_dic = json.load(open(filename_1))



PP_dist_1 = PP_1_dic['dist']
PP_time_1 = PP_1_dic['time']
PP_lat_1 = PP_1_dic['lat']
PP_lon_1 = PP_1_dic['lon']
PP_depth_1 = PP_1_dic['depth']


print('Importing Piercing Points '+PHASES[1])
print('\n')

filename_2 = PP_DIR+'PP_'+PHASES[1]+'_dic.json'

PP_2_dic = json.load(open(filename_2))



PP_depth_2 = PP_2_dic['dist']
PP_time_2 = PP_2_dic['time']
PP_lat_2 = PP_2_dic['lat']
PP_lon_2 = PP_2_dic['lon']
PP_depth_2 = PP_2_dic['depth']

print('Piercing Points - '+str(DEPTH_1))
print('\n')

pp_1_lat  = [[]]*len(PP_lon_1)
pp_1_long  = [[]]*len(PP_lon_1)


for i,j in enumerate(PP_lon_1):
    for k,l in enumerate(j):
        if LLCRNRLON_SMALL<= l <= URCRNRLON_SMALL and PP_depth_1[i][k] == DEPTH_1:
                pp_1_lat[i] = PP_lat_1[i][k] 
                pp_1_long[i] = l


print('Piercing Points - '+str(DEPTH_2))
print('\n')


pp_2_lat  = [[]]*len(PP_lon_2)
pp_2_long  = [[]]*len(PP_lon_2)


for i,j in enumerate(PP_lon_2):
    for k,l in enumerate(j):
        if LLCRNRLON_SMALL <= l <= URCRNRLON_SMALL and PP_depth_2[i][k] == DEPTH_2:
                pp_2_lat[i] = PP_lat_2[i][k] 
                pp_2_long[i] = l


print('Creating GRID POINTS')
print('\n')

area = (LLCRNRLON_SMALL,URCRNRLON_SMALL, LLCRNRLAT_SMALL, URCRNRLAT_SMALL)

shape = (abs(abs(URCRNRLON_SMALL) - abs(LLCRNRLON_SMALL))*3, abs(abs(URCRNRLAT_SMALL) - abs(LLCRNRLAT_SMALL))*3)

grdx, grdy = gridder.regular(area, shape)


print('Filtering grid points')
print('\n')


dist_pp_grid_min = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_min[i] = [np.sqrt((j - pp_1_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_1_lat)]
    
dist_pp_grid_max = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_max[i] = [np.sqrt((j - pp_2_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_2_lat)]



dist_grid_piercing_points = DIST_GRID_PP
number_piercing_points_per_bin = NUMBER_PP_PER_BIN 

grid_sel_min = []
grid_sel_min_data = []
for i,j in enumerate(dist_pp_grid_min):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[number_piercing_points_per_bin]] < dist_grid_piercing_points:
        grid_sel_min.append((grdx[i],grdy[i]))
        
grid_sel_max = []
grid_sel_min_data = []

for i,j in enumerate(dist_pp_grid_max):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[number_piercing_points_per_bin]] < dist_grid_piercing_points:
        grid_sel_max.append((grdx[i],grdy[i]))


grid_sel = grid_sel_min+grid_sel_max


grid_selected = set(map(tuple,grid_sel))


grid_sel_x = []
grid_sel_y = []

for i,j in enumerate(grid_selected):
    grid_sel_x.append(j[0])
    grid_sel_y.append(j[1])


print('Importing depths and times to the Ps conversion to each event for all stations')
print('\n')

filename_Pds = PdS_DIR+'Pds_dic.json'

PdS_Dic = json.load(open(filename_Pds))

Pds_dist = PdS_Dic['dist']
Pds_time = PdS_Dic['time']
Pds_depth = PdS_Dic['depth']


print('Importing depths and times to the Ppds conversion to each event for all stations')
print('\n')

filename_Ppds = PdS_DIR+'PPvs_dic.json'

Ppds_Dic = json.load(open(filename_Ppds))

Ppds_dist = Ppds_Dic['dist']
Ppds_time = Ppds_Dic['time']
Ppds_depth = Ppds_Dic['depth']

print('Migrating Pds data...')
print('\n')

RF_amplitude_time_Pds = [[]]*len(Pds_depth)
RF_amplitude_depth_Pds = [[]]*len(Pds_depth)
for i,j in enumerate(Pds_depth):
    sta_t_Pds = j
    RF_t_Pds = camadas_terra_10_km
    RF_amplitude_time_Pds[i] = [Pds_time[i][sta_t_Pds.index(l)] if l in sta_t_Pds else -1 for k,l in enumerate(RF_t_Pds)]
    RF_amplitude_depth_Pds[i] = [sta_t_Pds[sta_t_Pds.index(l)] for k,l in enumerate(RF_t_Pds) if l in sta_t_Pds]

RF_amplitude_Pds = [[]]*len(RF_amplitude_time_Pds)

for i,j in enumerate(RF_amplitude_time_Pds):
    sta_t_Pds = [round(l,1) for k,l in enumerate(sta_time[i])]
    RF_t_Pds = [round(l,1) for k,l in enumerate(j)]
    RF_amplitude_Pds[i] = [sta_data[i][sta_t_Pds.index(l)] if l != -1 else 0 for k,l in enumerate(RF_t_Pds)]



print('Migrating Ppds data...')
print('\n')

RF_amplitude_time_Ppds = [[]]*len(Ppds_depth)
RF_amplitude_depth_Ppds = [[]]*len(Ppds_depth)
for i,j in enumerate(Ppds_depth):
    sta_t_Ppds = j
    RF_t_Ppds = camadas_terra_10_km
    RF_amplitude_time_Ppds[i] = [Ppds_time[i][sta_t_Ppds.index(l)] if l in sta_t_Ppds else -1 for k,l in enumerate(RF_t_Ppds)]
    RF_amplitude_depth_Ppds[i] = [sta_t_Ppds[sta_t_Ppds.index(l)] for k,l in enumerate(RF_t_Ppds) if l in sta_t_Ppds]

RF_amplitude_Ppds = [[]]*len(RF_amplitude_time_Ppds)

for i,j in enumerate(RF_amplitude_time_Ppds):
    sta_t_Ppds = [round(l,1) for k,l in enumerate(sta_time[i])]
    RF_t_Ppds = [round(l,1) for k,l in enumerate(j)]
    RF_amplitude_Ppds[i] = [sta_data[i][sta_t_Ppds.index(l)] if l != -1 else 0 for k,l in enumerate(RF_t_Ppds)]


print('Data stacking in each point of the filtered grid')
print('\n')

dados_grid_lat = pp_1_lat
dados_grid_lon = pp_1_long

if LINEAR_STACKING == True: 
	RF_data_raw_Pds = [[]]*len(grid_sel_x)
	RF_amplitude_depth_raw_Pds = [[]]*len(grid_sel_x)

	RF_data_raw_Ppds = [[]]*len(grid_sel_x)
	RF_amplitude_depth_raw_Ppds = [[]]*len(grid_sel_x)

	for i,j in enumerate(grid_sel_x):
		RF_data_raw_Pds[i] = [RF_amplitude_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]
		RF_amplitude_depth_raw_Pds[i] = [RF_amplitude_depth_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]

		RF_data_raw_Ppds[i] = [RF_amplitude_Ppds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]
		RF_amplitude_depth_raw_Ppds[i] = [RF_amplitude_depth_Ppds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]


	if DEPTH_ESTIMATION == True: 
		print('Depth mean and std estimation in each point of the filtered grid')
		print('\n')

		RF_DEPTH_raw_1_Pds = [[]]*len(RF_data_raw_Pds)
		RF_DEPTH_raw_1_Ppds = [[]]*len(RF_data_raw_Ppds)
		for i,j in enumerate(RF_data_raw_Pds):
			lst_DEPTH_raw = []
			for k,l in enumerate(j):
				lst_depth_amp = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if DEPTH_1-DEPTH_RANGE < c < DEPTH_1+DEPTH_RANGE]
				lst_depth_pp = [c for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if DEPTH_1-DEPTH_RANGE < c < DEPTH_1+DEPTH_RANGE]
				lst_DEPTH_raw.append(lst_depth_pp[lst_depth_amp.index(max(lst_depth_amp))])
			RF_DEPTH_raw_1_Pds[i] = lst_DEPTH_raw

		for i,j in enumerate(RF_data_raw_Ppds):
			lst_DEPTH_raw = []
			for k,l in enumerate(j):
				lst_depth_amp = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Ppds[i][k]) if DEPTH_1-DEPTH_RANGE < c < DEPTH_1+DEPTH_RANGE]
				lst_depth_pp = [c for x,c in enumerate(RF_amplitude_depth_raw_Ppds[i][k]) if DEPTH_1-DEPTH_RANGE < c < DEPTH_1+DEPTH_RANGE]
				lst_DEPTH_raw.append(lst_depth_pp[lst_depth_amp.index(max(lst_depth_amp))])
			RF_DEPTH_raw_1_Ppds[i] = lst_DEPTH_raw
			

		RF_DEPTH_mean_1_Pds = []
		RF_DEPTH_std_1_Pds = []

		RF_DEPTH_mean_1_Ppds = []
		RF_DEPTH_std_1_Ppds = []

		RF_DEPTH_mean_1_true_Pds = []
		RF_DEPTH_std_1_true_Pds = []

		RF_DEPTH_mean_1_true_Ppds = []
		RF_DEPTH_std_1_true_Ppds = []

		for i,j in enumerate(RF_DEPTH_raw_1_Pds):
				print('Bootstrap estimation of the true depths')
				print('\n')
				if len(j) > NUMBER_PP_PER_BIN:
					std_1_lst_Pds = []
					std_1_lst_Ppds = []


					std_1_lst_true_Pds = []
					std_1_lst_true_Ppds = []

					lst_1_std_true_Pds  = []
					lst_1_std_true_Ppds  = []

					Vp_Vs_ratio_depth_1 = Vs_depth_1/Vp_depth_1
					alfa = Vp_depth_1 - Vs_depth_1
					beta = Vp_depth_1 + Vs_depth_1
					gamma_vp_vs = GAMMA*Vp_Vs_ratio_depth_1

					depth_Pds = j
					depth_Ppds = RF_DEPTH_raw_1_Ppds[i]

					delta_Vp = [(alfa*beta*(depth_Ppds[_t] - _y))/(alfa*(1+gamma_vp_vs)*_y - 
								   (beta*(1-gamma_vp_vs)*depth_Ppds[_t])) for _t,_y in enumerate(depth_Pds)]

					delta_Vs = [_delta_vp * GAMMA * Vp_Vs_ratio_depth_1 for _delta_vp in delta_Vp]

					depth_true_Pds = [_y * ((Vp_depth_1-Vs_depth_1)/(Vp_depth_1*Vs_depth_1)) *
											 (((Vs_depth_1+delta_Vs[_t])*(Vp_depth_1+delta_Vp[_t]))/(Vp_depth_1+delta_Vp[_t]-Vs_depth_1-delta_Vs[_t]))
											 for _t,_y in enumerate(depth_Pds)] 

					depth_true_Ppds = [_y * ((Vp_depth_1-Vs_depth_1)/(Vp_depth_1*Vs_depth_1)) *
											 (((Vs_depth_1+delta_Vs[_t])*(Vp_depth_1+delta_Vp[_t]))/(Vp_depth_1+delta_Vp[_t]-Vs_depth_1-delta_Vs[_t]))
											 for _t,_y in enumerate(depth_Ppds)] 	
							

					RF_DEPTH_mean_1_Pds.append(np.mean(j))
					RF_DEPTH_mean_1_Ppds.append(np.mean(RF_DEPTH_raw_1_Ppds[i]))

					RF_DEPTH_mean_1_true_Pds.append(np.mean(depth_true_Pds))
					RF_DEPTH_mean_1_true_Ppds.append(np.mean(depth_true_Ppds))

					for _k in range(BOOTSTRAP_INTERATOR):
						print('Bootstrap estimation '+str(1+_k))
						BOOTSTRAP_LST_Pds = [random.choice(j) for _ in j]
						BOOTSTRAP_LST_Ppds = [random.choice(j) for _ in RF_DEPTH_raw_1_Ppds[i]]


						std_1_lst_Pds.append(np.std(BOOTSTRAP_LST_Pds))
						std_1_lst_Ppds.append(np.std(BOOTSTRAP_LST_Ppds))
						
						print('Depth - mean = '+str(np.mean(BOOTSTRAP_LST_Pds))+' -  Standard deviation = '+str(np.mean(std_1_lst_Pds))+' of Pds Conversions')
						print('Depth - mean = '+str(np.mean(BOOTSTRAP_LST_Ppds))+' -  Standard deviation = '+str(np.mean(std_1_lst_Ppds))+' of Ppds Conversions')
			
						depth_Pds_bootstrap = BOOTSTRAP_LST_Pds
						depth_Ppds_bootstrap = BOOTSTRAP_LST_Ppds

						delta_Vp_bootstrap = [(alfa*beta*(depth_Ppds_bootstrap[_t] - _y))/(alfa*(1+gamma_vp_vs)*_y - 
									   (beta*(1-gamma_vp_vs)*depth_Ppds_bootstrap[_t])) for _t,_y in enumerate(depth_Pds_bootstrap)]

						delta_Vs_bootstrap = [_delta_vp * GAMMA * Vp_Vs_ratio_depth_1 for _delta_vp in delta_Vp_bootstrap]

						depth_true_Pds_bootstrap = [_y * ((Vp_depth_1-Vs_depth_1)/(Vp_depth_1*Vs_depth_1)) *
												 (((Vs_depth_1+delta_Vs_bootstrap[_t])*(Vp_depth_1+delta_Vp_bootstrap[_t]))/
												(Vp_depth_1+delta_Vp_bootstrap[_t]-Vs_depth_1-delta_Vs[_t]))
												 for _t,_y in enumerate(depth_Pds_bootstrap)] 


						depth_true_Ppds_bootstrap = [_y * ((Vp_depth_1-Vs_depth_1)/(Vp_depth_1*Vs_depth_1)) *
												 (((Vs_depth_1+delta_Vs_bootstrap[_t])*(Vp_depth_1+delta_Vp_bootstrap[_t]))/
												(Vp_depth_1+delta_Vp_bootstrap[_t]-Vs_depth_1-delta_Vs[_t]))
												 for _t,_y in enumerate(depth_Ppds_bootstrap)] 			


						lst_1_std_true_Pds.append(np.std(depth_true_Pds_bootstrap))
						lst_1_std_true_Ppds.append(np.std(depth_true_Ppds_bootstrap))

						print('True Depth - mean = '+str(RF_DEPTH_mean_1_true_Pds[i])+' -  Standard deviation = '+str(np.mean(lst_1_std_true_Pds))+' of Pds Conversions')
						print('True Depth - mean = '+str(RF_DEPTH_mean_1_true_Ppds[i])+' -  Standard deviation = '+str(np.mean(lst_1_std_true_Ppds))+' of Ppds Conversions')

						

					RF_DEPTH_std_1_Pds.append(np.mean(std_1_lst_Pds))

					RF_DEPTH_std_1_Ppds.append(np.mean(std_1_lst_Ppds))

					RF_DEPTH_std_1_true_Pds.append(np.mean(lst_1_std_true_Pds))


					RF_DEPTH_std_1_true_Ppds.append(np.mean(lst_1_std_true_Ppds))
					print('\n')
					print('Depth - mean = '+str(RF_DEPTH_mean_1_Pds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_1_Pds[i])+' of Pds Conversions')
					print('\n')
					print('Depth - mean = '+str(RF_DEPTH_mean_1_Ppds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_1_Ppds[i])+' of Ppds Conversions')

					print('\n')
					print('True Depth - mean = '+str(RF_DEPTH_mean_1_true_Pds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_1_true_Pds[i])+' of Pds Conversions')
					print('\n')
					print('True Depth - mean = '+str(RF_DEPTH_mean_1_true_Ppds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_1_true_Ppds[i])+' of Ppds Conversions')
					print('\n')
				else:

					RF_DEPTH_std_1_Pds.append(0)
					RF_DEPTH_mean_1_Pds.append(0)

					RF_DEPTH_std_1_Ppds.append(0)
					RF_DEPTH_mean_1_Ppds.append(0)

					RF_DEPTH_mean_1_true_Pds.append(0)
					RF_DEPTH_std_1_true_Pds.append(0)

					RF_DEPTH_mean_1_true_Ppds.append(0)
					RF_DEPTH_std_1_true_Ppds.append(0)

					print('\n')
					print('Depth - mean = '+str(RF_DEPTH_mean_1_Pds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_1_Pds[i])+' of Pds Conversions')
					print('\n')
					print('Depth - mean = '+str(RF_DEPTH_mean_1_Ppds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_1_Ppds[i])+' of Ppds Conversions')

					print('\n')
					print('True Depth - mean = '+str(RF_DEPTH_mean_1_true_Pds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_1_true_Pds[i])+' of Pds Conversions')
					print('\n')
					print('True Depth - mean = '+str(RF_DEPTH_mean_1_true_Ppds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_1_true_Ppds[i])+' of Ppds Conversions')
					print('\n')





	if DEPTH_ESTIMATION == True: 
		print('Depth mean and std estimation in each point of the filtered grid')
		print('\n')

		RF_DEPTH_raw_2_Pds = [[]]*len(RF_data_raw_Pds)
		RF_DEPTH_raw_2_Ppds = [[]]*len(RF_data_raw_Ppds)
		for i,j in enumerate(RF_data_raw_Pds):
			lst_DEPTH_raw = []
			for k,l in enumerate(j):
				lst_depth_amp = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if DEPTH_2-DEPTH_RANGE < c < DEPTH_2+DEPTH_RANGE]
				lst_depth_pp = [c for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if DEPTH_2-DEPTH_RANGE < c < DEPTH_2+DEPTH_RANGE]
				lst_DEPTH_raw.append(lst_depth_pp[lst_depth_amp.index(max(lst_depth_amp))])
			RF_DEPTH_raw_2_Pds[i] = lst_DEPTH_raw

		for i,j in enumerate(RF_data_raw_Ppds):
			lst_DEPTH_raw = []
			for k,l in enumerate(j):
				lst_depth_amp = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Ppds[i][k]) if DEPTH_2-DEPTH_RANGE < c < DEPTH_2+DEPTH_RANGE]
				lst_depth_pp = [c for x,c in enumerate(RF_amplitude_depth_raw_Ppds[i][k]) if DEPTH_2-DEPTH_RANGE < c < DEPTH_2+DEPTH_RANGE]
				lst_DEPTH_raw.append(lst_depth_pp[lst_depth_amp.index(max(lst_depth_amp))])
			RF_DEPTH_raw_2_Ppds[i] = lst_DEPTH_raw
			
		RF_DEPTH_delta_Vp = []

		RF_DEPTH_mean_2_Pds = []
		RF_DEPTH_std_2_Pds = []

		RF_DEPTH_mean_2_Ppds = []
		RF_DEPTH_std_2_Ppds = []

		RF_DEPTH_mean_2_true_Pds = []
		RF_DEPTH_std_2_true_Pds = []

		RF_DEPTH_mean_2_true_Ppds = []
		RF_DEPTH_std_2_true_Ppds = []

		for i,j in enumerate(RF_DEPTH_raw_2_Pds):
				print('Bootstrap estimation of the true depths')
				print('\n')
				if len(j) > NUMBER_PP_PER_BIN:
					std_2_lst_Pds = []
					std_2_lst_Ppds = []


					std_2_lst_true_Pds = []
					std_2_lst_true_Ppds = []

					lst_2_std_true_Pds  = []
					lst_2_std_true_Ppds  = []

					Vp_Vs_ratio_depth_2 = Vs_depth_2/Vp_depth_2
					alfa = Vp_depth_2 - Vs_depth_2
					beta = Vp_depth_2 + Vs_depth_2
					gamma_vp_vs = GAMMA*Vp_Vs_ratio_depth_2

					depth_Pds = j
					depth_Ppds = RF_DEPTH_raw_2_Ppds[i]

					delta_Vp = [(alfa*beta*(depth_Ppds[_t] - _y))/(alfa*(1+gamma_vp_vs)*_y - 
								   (beta*(1-gamma_vp_vs)*depth_Ppds[_t])) for _t,_y in enumerate(depth_Pds)]
					

					delta_Vs = [_delta_vp * GAMMA * Vp_Vs_ratio_depth_2 for _delta_vp in delta_Vp]

					depth_true_Pds = [_y * ((Vp_depth_2-Vs_depth_2)/(Vp_depth_2*Vs_depth_2)) *
											 (((Vs_depth_2+delta_Vs[_t])*(Vp_depth_2+delta_Vp[_t]))/(Vp_depth_2+delta_Vp[_t]-Vs_depth_2-delta_Vs[_t]))
											 for _t,_y in enumerate(depth_Pds)] 

					depth_true_Ppds = [_y * ((Vp_depth_2-Vs_depth_2)/(Vp_depth_2*Vs_depth_2)) *
											 (((Vs_depth_2+delta_Vs[_t])*(Vp_depth_2+delta_Vp[_t]))/(Vp_depth_2+delta_Vp[_t]-Vs_depth_2-delta_Vs[_t]))
											 for _t,_y in enumerate(depth_Ppds)] 	
							
					RF_DEPTH_delta_Vp.append(np.mean(delta_Vp))

					RF_DEPTH_mean_2_Pds.append(np.mean(j))
					RF_DEPTH_mean_2_Ppds.append(np.mean(RF_DEPTH_raw_2_Ppds[i]))

					RF_DEPTH_mean_2_true_Pds.append(np.mean(depth_true_Pds))
					RF_DEPTH_mean_2_true_Ppds.append(np.mean(depth_true_Ppds))

					for _k in range(BOOTSTRAP_INTERATOR):
						print('Bootstrap estimation '+str(1+_k))
						BOOTSTRAP_LST_Pds = [random.choice(j) for _ in j]
						BOOTSTRAP_LST_Ppds = [random.choice(j) for _ in RF_DEPTH_raw_2_Ppds[i]]


						std_2_lst_Pds.append(np.std(BOOTSTRAP_LST_Pds))
						std_2_lst_Ppds.append(np.std(BOOTSTRAP_LST_Ppds))
						
						print('Depth - mean = '+str(np.mean(BOOTSTRAP_LST_Pds))+' -  Standard deviation = '+str(np.mean(std_2_lst_Pds))+' of Pds Conversions')
						print('Depth - mean = '+str(np.mean(BOOTSTRAP_LST_Ppds))+' -  Standard deviation = '+str(np.mean(std_2_lst_Ppds))+' of Ppds Conversions')
			
						depth_Pds_bootstrap = BOOTSTRAP_LST_Pds
						depth_Ppds_bootstrap = BOOTSTRAP_LST_Ppds

						delta_Vp_bootstrap = [(alfa*beta*(depth_Ppds_bootstrap[_t] - _y))/(alfa*(1+gamma_vp_vs)*_y - 
									   (beta*(1-gamma_vp_vs)*depth_Ppds_bootstrap[_t])) for _t,_y in enumerate(depth_Pds_bootstrap)]

						delta_Vs_bootstrap = [_delta_vp * GAMMA * Vp_Vs_ratio_depth_2 for _delta_vp in delta_Vp_bootstrap]

						depth_true_Pds_bootstrap = [_y * ((Vp_depth_2-Vs_depth_2)/(Vp_depth_2*Vs_depth_2)) *
												 (((Vs_depth_2+delta_Vs_bootstrap[_t])*(Vp_depth_2+delta_Vp_bootstrap[_t]))/
												(Vp_depth_2+delta_Vp_bootstrap[_t]-Vs_depth_2-delta_Vs[_t]))
												 for _t,_y in enumerate(depth_Pds_bootstrap)] 


						depth_true_Ppds_bootstrap = [_y * ((Vp_depth_2-Vs_depth_2)/(Vp_depth_2*Vs_depth_2)) *
												 (((Vs_depth_2+delta_Vs_bootstrap[_t])*(Vp_depth_2+delta_Vp_bootstrap[_t]))/
												(Vp_depth_2+delta_Vp_bootstrap[_t]-Vs_depth_2-delta_Vs[_t]))
												 for _t,_y in enumerate(depth_Ppds_bootstrap)] 			


						lst_2_std_true_Pds.append(np.std(depth_true_Pds_bootstrap))
						lst_2_std_true_Ppds.append(np.std(depth_true_Ppds_bootstrap))

						print('True Depth - mean = '+str(RF_DEPTH_mean_2_true_Pds[i])+' -  Standard deviation = '+str(np.mean(lst_2_std_true_Pds))+' of Pds Conversions')
						print('True Depth - mean = '+str(RF_DEPTH_mean_2_true_Ppds[i])+' -  Standard deviation = '+str(np.mean(lst_2_std_true_Ppds))+' of Ppds Conversions')

						

					RF_DEPTH_std_2_Pds.append(np.mean(std_2_lst_Pds))

					RF_DEPTH_std_2_Ppds.append(np.mean(std_2_lst_Ppds))

					RF_DEPTH_std_2_true_Pds.append(np.mean(lst_2_std_true_Pds))


					RF_DEPTH_std_2_true_Ppds.append(np.mean(lst_2_std_true_Ppds))
					print('\n')
					print('Depth - mean = '+str(RF_DEPTH_mean_2_Pds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_2_Pds[i])+' of Pds Conversions')
					print('\n')
					print('Depth - mean = '+str(RF_DEPTH_mean_2_Ppds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_2_Ppds[i])+' of Ppds Conversions')

					print('\n')
					print('True Depth - mean = '+str(RF_DEPTH_mean_2_true_Pds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_2_true_Pds[i])+' of Pds Conversions')
					print('\n')
					print('True Depth - mean = '+str(RF_DEPTH_mean_2_true_Ppds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_2_true_Ppds[i])+' of Ppds Conversions')
					print('\n')
				else:


					RF_DEPTH_delta_Vp.append(0)

					RF_DEPTH_std_2_Pds.append(0)
					RF_DEPTH_mean_2_Pds.append(0)

					RF_DEPTH_std_2_Ppds.append(0)
					RF_DEPTH_mean_2_Ppds.append(0)

					RF_DEPTH_mean_2_true_Pds.append(0)
					RF_DEPTH_std_2_true_Pds.append(0)

					RF_DEPTH_mean_2_true_Ppds.append(0)
					RF_DEPTH_std_2_true_Ppds.append(0)

					print('\n')
					print('Depth - mean = '+str(RF_DEPTH_mean_2_Pds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_2_Pds[i])+' of Pds Conversions')
					print('\n')
					print('Depth - mean = '+str(RF_DEPTH_mean_2_Ppds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_2_Ppds[i])+' of Ppds Conversions')

					print('\n')
					print('True Depth - mean = '+str(RF_DEPTH_mean_2_true_Pds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_2_true_Pds[i])+' of Pds Conversions')
					print('\n')
					print('True Depth - mean = '+str(RF_DEPTH_mean_2_true_Ppds[i])+' -  Standard deviation = '+str(RF_DEPTH_std_2_true_Ppds[i])+' of Ppds Conversions')
					print('\n')



	RF_stacking_Pds = []
	len_RF_stacking_Pds = []


	for i,j in enumerate(RF_data_raw_Pds):
		if len(j) > NUMBER_PP_PER_BIN:
			RF_stacking_Pds.append([sum(x)/len(j)  for x in zip(*j)])
			len_RF_stacking_Pds.append(len(j))
		else:
			RF_stacking_Pds.append([])
			len_RF_stacking_Pds.append(0),

	RF_stacking_Ppds = []
	len_RF_stacking_Ppds = []


	for i,j in enumerate(RF_data_raw_Ppds):
		if len(j) > NUMBER_PP_PER_BIN:
			RF_stacking_Ppds.append([sum(x)/len(j)  for x in zip(*j)])
			len_RF_stacking_Ppds.append(len(j))
		else:
			RF_stacking_Ppds.append([])
			len_RF_stacking_Ppds.append(0)

else: 
	pass

RF_lst_delta_Vp = []
RF_lst_DEPTH_mean_1_true_Pds = []
RF_lst_lon_mean_1_true_Pds = []
RF_lst_lat_mean_1_true_Pds = []

for _i,true_Pds_1 in enumerate(RF_DEPTH_mean_1_true_Pds):
	if true_Pds_1 > 0:
		RF_lst_delta_Vp.append(RF_DEPTH_delta_Vp[_i])		
		RF_lst_lon_mean_1_true_Pds.append(grid_sel_x[_i])
		RF_lst_lat_mean_1_true_Pds.append(grid_sel_y[_i])
		RF_lst_DEPTH_mean_1_true_Pds.append(true_Pds_1)


RF_lst_DEPTH_mean_1_true_Ppds = []
RF_lst_lon_mean_1_true_Ppds = []
RF_lst_lat_mean_1_true_Ppds = []

for _i,true_Ppds_1 in enumerate(RF_DEPTH_mean_1_true_Ppds):
	if true_Ppds_1 > 0:
		RF_lst_lon_mean_1_true_Ppds.append(grid_sel_x[_i])
		RF_lst_lat_mean_1_true_Ppds.append(grid_sel_y[_i])
		RF_lst_DEPTH_mean_1_true_Ppds.append(true_Ppds_1)

RF_lst_DEPTH_mean_2_true_Pds = []
RF_lst_lon_mean_2_true_Pds = []
RF_lst_lat_mean_2_true_Pds = []

for _i,true_Pds_2 in enumerate(RF_DEPTH_mean_2_true_Pds):
	if true_Pds_2 > 0:
		RF_lst_lon_mean_2_true_Pds.append(grid_sel_x[_i])
		RF_lst_lat_mean_2_true_Pds.append(grid_sel_y[_i])
		RF_lst_DEPTH_mean_2_true_Pds.append(true_Pds_2)


RF_lst_DEPTH_mean_2_true_Ppds = []
RF_lst_lon_mean_2_true_Ppds = []
RF_lst_lat_mean_2_true_Ppds = []

for _i,true_Ppds_2 in enumerate(RF_DEPTH_mean_2_true_Ppds):
	if true_Ppds_2 > 0:
		RF_lst_lon_mean_2_true_Ppds.append(grid_sel_x[_i])
		RF_lst_lat_mean_2_true_Ppds.append(grid_sel_y[_i])
		RF_lst_DEPTH_mean_2_true_Ppds.append(true_Ppds_2)
		

## FUnction to set midpoint in plots

class MidPointNorm(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)  
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0: 
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint



fig, axes = plt.subplots(nrows=2, ncols=2,figsize=(10,10),squeeze=False,sharex=False,sharey=False)

ax1 = axes[0, 0]
ax = axes[0, 1]
ax2 = axes[1, 0]
ax3 = axes[1, 1]

colormap = 'seismic'

m1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax1)

m1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

norm1 = MidPointNorm(midpoint=DEPTH_1)
x, y = m1(RF_lst_lon_mean_1_true_Pds,RF_lst_lat_mean_1_true_Pds)
sc = m1.scatter(x,y,40,RF_lst_DEPTH_mean_1_true_Pds,cmap=colormap,marker='o',norm=norm1)


for lon, lat in zip(sta_long,sta_lat):
    x,y = m1(lon, lat)
    msize = 10
    l1, = m1.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m1.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m1.drawcoastlines(color='k',zorder=1)
m1.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m1.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])


ax1.set_title(str(int(DEPTH_1))+' Pds conversion', y=1.08)



m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

norm = MidPointNorm(midpoint=DEPTH_1)
x, y = m(RF_lst_lon_mean_1_true_Ppds,RF_lst_lat_mean_1_true_Ppds)
sc1 = m.scatter(x,y,40,RF_lst_DEPTH_mean_1_true_Ppds,cmap=colormap,marker='o',norm=norm)



for lon, lat in zip(sta_long,sta_lat):
    x,y = m(lon, lat)
    msize = 10
    l1, = m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m.drawcoastlines(color='k',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax.set_title(str(int(DEPTH_1))+' Ppds conversion', y=1.08)





m2 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax2)

m2.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m2.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

norm2 = MidPointNorm(midpoint=DEPTH_2)
x, y = m2(RF_lst_lon_mean_2_true_Pds,RF_lst_lat_mean_2_true_Pds)
sc = m2.scatter(x,y,40,RF_lst_DEPTH_mean_2_true_Pds,cmap=colormap,marker='o',norm=norm2)


for lon, lat in zip(sta_long,sta_lat):
    x,y = m2(lon, lat)
    msize = 10
    l1, = m2.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m2.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m2.drawcoastlines(color='k',zorder=1)
m2.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m2.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax2.set_title(str(int(DEPTH_2))+' Pds conversion', y=1.08)



m3 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax3)

m3.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m3.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

norm3 = MidPointNorm(midpoint=DEPTH_2)
x, y = m3(RF_lst_lon_mean_2_true_Ppds,RF_lst_lat_mean_2_true_Ppds)
sc2 = m3.scatter(x,y,40,RF_lst_DEPTH_mean_2_true_Ppds,cmap=colormap,marker='o',norm=norm3)


for lon, lat in zip(sta_long,sta_lat):
    x,y = m3(lon, lat)
    msize = 10
    l1, = m3.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m3.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m3.drawcoastlines(color='k',zorder=1)
m3.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m3.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax3.set_title(str(int(DEPTH_2))+' Ppds conversion', y=1.08)

cbar_ax1 = fig.add_axes([0.2, 0.5, 0.6, 0.02])
fig.colorbar(sc1, cax=cbar_ax1,orientation='horizontal')


cbar_ax2 = fig.add_axes([0.2, 0.1, 0.6, 0.02])
fig.colorbar(sc2, cax=cbar_ax2,orientation='horizontal')


fig.savefig(PP_FIGURE+'TRUE_DEPTH_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)



#Delta VP of each bin


fig_delta_vp, ax_delta_vp = plt.subplots(nrows=1, ncols=1,figsize=(10,10),squeeze=False)


m_delta_vp = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL)

m_delta_vp.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_delta_vp.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

norm_delta_vp = MidPointNorm(midpoint=0)
x, y = m_delta_vp(RF_lst_lon_mean_1_true_Pds,RF_lst_lat_mean_1_true_Pds)
sc_delta_vp = m_delta_vp.scatter(x,y,40,RF_lst_delta_Vp,cmap=colormap,marker='o',norm=norm_delta_vp)


for lon, lat in zip(sta_long,sta_lat):
    x,y = m_delta_vp(lon, lat)
    msize = 10
    l1, = m_delta_vp.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_delta_vp.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_delta_vp.drawcoastlines(color='k',zorder=1)
m_delta_vp.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_delta_vp.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_delta_vp.colorbar(sc_delta_vp,orientation='horizontal')

fig_delta_vp.suptitle('Delta Vp for each bin')

#Thickness of the Mantle Transition Zone
thickness_MTZ_Pds = []
thickness_MTZ_Ppds = []

for i,j in enumerate(RF_lst_DEPTH_mean_2_true_Pds):
	thickness_MTZ_Pds.append(j-RF_lst_DEPTH_mean_1_true_Pds[i])
	thickness_MTZ_Ppds.append(RF_lst_DEPTH_mean_2_true_Ppds[i]-RF_lst_DEPTH_mean_1_true_Ppds[i])


fig_thickness, (ax_thickness1, ax_thickness2) = plt.subplots(nrows=1, ncols=2,figsize=(10,5))

m_thickness1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness1)

m_thickness1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

norm_thickness1 = MidPointNorm(midpoint=250)
x, y = m_thickness1(RF_lst_lon_mean_1_true_Pds,RF_lst_lat_mean_1_true_Pds)
sc_thickness1 = m_thickness1.scatter(x,y,40,thickness_MTZ_Pds,cmap=colormap,marker='o',norm=norm_thickness1)


for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness1(lon, lat)
    msize = 10
    l1, = m_thickness1.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness1.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness1.drawcoastlines(color='k',zorder=1)
m_thickness1.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness1.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness1,ax=ax_thickness1,orientation='horizontal')
ax_thickness1.set_title('Thickness of MTZ (Pds)',y=1.08)


m_thickness2 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness2)

m_thickness2.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness2.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

norm_thickness2 = MidPointNorm(midpoint=250)
x, y = m_thickness2(RF_lst_lon_mean_1_true_Pds,RF_lst_lat_mean_1_true_Pds)
sc_thickness2 = m_thickness2.scatter(x,y,40,thickness_MTZ_Pds,cmap=colormap,marker='o',norm=norm_thickness2)


for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness2(lon, lat)
    msize = 10
    l1, = m_thickness2.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness2.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness2.drawcoastlines(color='k',zorder=1)
m_thickness2.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness2.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness2,ax=ax_thickness2,orientation='horizontal')
ax_thickness2.set_title('Thickness of MTZ (Ppds)',y=1.08)


plt.show()
'''
#####################################################################	TERMINAR DE FAZER AS COMPRAÇÔES ENTRE TRUE E NORMAL!   #####################################################################


print('Saving Selected Piercing Points in JSON file')
print('\n')

os.makedirs(PP_SELEC_DIR,exist_ok=True)

SELECTED_BINNED_DATA_dic = {'lat':[],'lon':[],'len':[],'mean_1':[],'std_1':[],'mean_2':[],'std_2':[],'data':[]}
for i,j in enumerate(RF_stacking):
	SELECTED_BINNED_DATA_dic['lat'].append(grid_sel_y[i])
	SELECTED_BINNED_DATA_dic['lon'].append(grid_sel_x[i])
	SELECTED_BINNED_DATA_dic['len'].append(len_RF_stacking[i])
	SELECTED_BINNED_DATA_dic['mean_1'].append(RF_DEPTH_mean_1[i])
	SELECTED_BINNED_DATA_dic['std_1'].append(RF_DEPTH_std_1[i])
	SELECTED_BINNED_DATA_dic['mean_2'].append(RF_DEPTH_mean_2[i])
	SELECTED_BINNED_DATA_dic['std_2'].append(RF_DEPTH_std_2[i])
	SELECTED_BINNED_DATA_dic['data'].append(j)

with open(PP_SELEC_DIR+'SELECTED_BINNED_Ps.json', 'w') as fp:
	json.dump(SELECTED_BINNED_DATA_dic, fp)'''
