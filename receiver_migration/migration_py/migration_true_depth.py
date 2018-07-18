
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



from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,PHASES_LST,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,DEPTH_1,DEPTH_2,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP_MED,PHASES_PPvs_LST,DIST_GRID_PP,
					LINEAR_STACKING,DEPTH_ESTIMATION,DEPTH_RANGE,BOOTSTRAP_INTERATOR,BOOTSTRAP_DEPTH_ESTIMATION
				   )

print('Starting Receiver Functions migration code')
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

PdS_dist = PdS_Dic['dist']
PdS_time = PdS_Dic['time']
PdS_depth = PdS_Dic['depth']


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
			

#####################################################################	TERMINAR DE FAZER ISSO !!!!!! CANSEI!!!!!! #####################################################################
		RF_DEPTH_mean_1 = []
		RF_DEPTH_std_1 = []
		for i,j in enumerate(RF_DEPTH_raw_1_Pds):
				print('Bootstrap estimation of the true depths')
				print('\n')
				if len(j) > NUMBER_PP_PER_BIN:
					std_1_lst = []
					mean_1_lst = []
					for _k in range(BOOTSTRAP_INTERATOR):
						print('Bootstrap estimation '+str(1+_k))
						BOOTSTRAP_LST = [random.choice(j) for _ in j]
						std_1_lst.append(np.std(BOOTSTRAP_LST))
						print('Standard deviation = #'+str(std_1_lst[_k]))
					RF_DEPTH_std_1.append(np.mean(std_1_lst))
					RF_DEPTH_mean_1.append(np.mean(j))
					print('\n')
					print('Mean = '+str(RF_DEPTH_mean_1[i])+' -  Standard deviation = '+str(RF_DEPTH_std_1[i]))
				else:
					RF_DEPTH_mean_1.append([])
					RF_DEPTH_std_1.append([])

		RF_DEPTH_raw_2 = [[]]*len(RF_data_raw)
		for i,j in enumerate(RF_data_raw):
			lst_DEPTH_raw = []
			for k,l in enumerate(j):
				lst_depth_amp = [l[x] for x,c in enumerate(RF_amplitude_depth_raw[i][k]) if DEPTH_2-DEPTH_RANGE < c < DEPTH_2+DEPTH_RANGE]
				lst_depth_pp = [c for x,c in enumerate(RF_amplitude_depth_raw[i][k]) if DEPTH_2-DEPTH_RANGE < c < DEPTH_2+DEPTH_RANGE]
				lst_DEPTH_raw.append(lst_depth_pp[lst_depth_amp.index(max(lst_depth_amp))])
			RF_DEPTH_raw_2[i] = lst_DEPTH_raw
			
		RF_DEPTH_mean_2 = []
		RF_DEPTH_std_2 = []
		for i,j in enumerate(RF_DEPTH_raw_2):
			if BOOTSTRAP_DEPTH_ESTIMATION == True:
				print('Bootstrap estimation of the true depths')
				print('\n')
				if len(j) > NUMBER_PP_PER_BIN:
					std_2_lst = []
					mean_2_lst = []
					for _k in range(BOOTSTRAP_INTERATOR):
						print('Bootstrap estimation #'+str(1+_k))
						BOOTSTRAP_LST = [random.choice(j) for _ in j]
						std_2_lst.append(np.std(BOOTSTRAP_LST))
						print('Standard deviation = '+str(std_2_lst[_k]))
					RF_DEPTH_std_2.append(np.mean(std_2_lst))
					RF_DEPTH_mean_2.append(np.mean(j))
					print('\n')
					print('Mean = '+str(RF_DEPTH_mean_2[i])+' -  Standard deviation = '+str(RF_DEPTH_std_2[i]))
				else:
					RF_DEPTH_mean_2.append([])
					RF_DEPTH_std_2.append([])
			else:
				if len(j) > NUMBER_PP_PER_BIN:
					RF_DEPTH_mean_2.append(np.mean(j))
					RF_DEPTH_std_2.append(np.std(j))
				else:
					RF_DEPTH_mean_2.append([])
					RF_DEPTH_std_2.append([])

	RF_stacking = []
	len_RF_stacking = []


	for i,j in enumerate(RF_data_raw):
		if len(j) > NUMBER_PP_PER_BIN:
			RF_stacking.append([sum(x)/len(j)  for x in zip(*j)])
			len_RF_stacking.append(len(j))
		else:
			RF_stacking.append([])
			len_RF_stacking.append(0)
else: 
	pass

fig_receiver_function_per_bin=plt.figure(figsize=(20,10))

project_Lat = PROJECT_LAT
project_Lon = PROJECT_LON

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)


lats = grid_sel_y
lons = grid_sel_x
RF_number = len_RF_stacking
x, y = m(lons,lats)
sc = m.scatter(x,y,40,RF_number,cmap='viridis',marker='o')
plt.colorbar(sc,label='Total of Receiver Functions per bin')


for lon, lat in zip(sta_long,sta_lat):
    x,y = m(lon, lat)
    msize = 10
    l1, = m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m.drawcoastlines(color='k',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

plt.legend([l1,sc],['Stations','Filtered Grid Points'],scatterpoints=1, frameon=True,labelspacing=1, loc='upper right',facecolor='w')
plt.title('Receiver Functions per bin', y=1.08)
fig_receiver_function_per_bin.savefig(PP_FIGURE+'RECEIVER_FUNCTION_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

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
	json.dump(SELECTED_BINNED_DATA_dic, fp)
