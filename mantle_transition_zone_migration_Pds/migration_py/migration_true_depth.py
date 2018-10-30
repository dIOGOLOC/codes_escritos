	
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
import collections





from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP_MED,DIST_GRID_PP,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,BOOTSTRAP_DEPTH_ESTIMATION,
					GAMMA
				   )

print('Starting Receiver Functions migration code to estimate the true depths of the Earth discontinuities')
print('\n')


print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)

for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[0] == 410:
		Vp_depth_1 = j[2]
		Vs_depth_1 = j[4]
		
for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[0] == 660:
		Vp_depth_2 = j[2]
		Vs_depth_2 = j[4]
		



print('Looking for Receiver Functions data in JSON file in '+STA_DIR)
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

print('Creating the Earth layered model')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

dist_med_camada_terra = [abs(c - ((410+660)/2)) for x,c in enumerate(camadas_terra_10_km)]

DEPTH_MED = camadas_terra_10_km[dist_med_camada_terra.index(min(dist_med_camada_terra))]

print('Importing Pds piercing points to each PHASE')
print('\n')

PHASES = 'P410s','P'+"{0:.0f}".format(DEPTH_MED)+'s','P660s'

print('Importing Pds Piercing Points for '+PHASES[0])
print('\n')

filename_1 = PP_DIR+'PP_'+PHASES[0]+'_dic.json'

PP_1_dic = json.load(open(filename_1))

PP_dist_1 = []
PP_time_1 = []
PP_lat_1 = []
PP_lon_1 = []
PP_depth_1 = [] 

for i,j in enumerate(PP_1_dic):
	PP_dist_1.append(j['dist'][0])
	PP_time_1.append(j['time'][0])
	PP_lat_1.append(j['lat'][0])
	PP_lon_1.append(j['lon'][0])
	PP_depth_1.append(j['depth'][0])

print('Importing Pds Piercing Points for '+PHASES[1])
print('\n')

filename_med = PP_DIR+'PP_'+PHASES[1]+'_dic.json'

PP_med_dic = json.load(open(filename_med))

PP_dist_med = []
PP_time_med = []
PP_lat_med = []
PP_lon_med = []
PP_depth_med = [] 

for i,j in enumerate(PP_med_dic):
	PP_dist_med.append(j['dist'][0])
	PP_time_med.append(j['time'][0])
	PP_lat_med.append(j['lat'][0])
	PP_lon_med.append(j['lon'][0])
	PP_depth_med.append(j['depth'][0])

print('Importing Pds Piercing Points for '+PHASES[2])
print('\n')

filename_2 = PP_DIR+'PP_'+PHASES[2]+'_dic.json'

PP_2_dic = json.load(open(filename_2))

PP_dist_2 = []
PP_time_2 = []
PP_lat_2 = []
PP_lon_2 = []
PP_depth_2 = [] 

for i,j in enumerate(PP_2_dic):
	PP_dist_2.append(j['dist'][0])
	PP_time_2.append(j['time'][0])
	PP_lat_2.append(j['lat'][0])
	PP_lon_2.append(j['lon'][0])
	PP_depth_2.append(j['depth'][0])

print('P410s Piercing Points')
print('\n')

pp_1_lat  = [[]]*len(PP_lon_1)
pp_1_long  = [[]]*len(PP_lon_1)


for i,j in enumerate(PP_lon_1):
    for k,l in enumerate(j):
        if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_1[i][k] == 410:
                pp_1_lat[i] = PP_lat_1[i][k] 
                pp_1_long[i] = l


print('Pds Piercing Points - '+"{0:.0f}".format(DEPTH_MED))
print('\n')

pp_med_lat  = [[]]*len(PP_lon_med)
pp_med_long  = [[]]*len(PP_lon_med)


for i,j in enumerate(PP_lon_med):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_med[i][k] == DEPTH_MED:
			pp_med_lat[i] = PP_lat_med[i][k] 
			pp_med_long[i] = l

print('P660s Piercing Points')
print('\n')


pp_2_lat  = [[]]*len(PP_lon_2)
pp_2_long  = [[]]*len(PP_lon_2)


for i,j in enumerate(PP_lon_2):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_2[i][k] == 660:
			pp_2_lat[i] = PP_lat_2[i][k]
			pp_2_long[i] = l


print('Creating GRID POINTS')
print('\n')

area = (LLCRNRLON_LARGE,URCRNRLON_LARGE, LLCRNRLAT_LARGE, URCRNRLAT_LARGE)

shape = (abs(abs(URCRNRLON_LARGE) - abs(LLCRNRLON_LARGE))*GRID_PP_MULT, abs(abs(URCRNRLAT_LARGE) - abs(LLCRNRLAT_LARGE))*GRID_PP_MULT)

grdx, grdy = gridder.regular(area, shape)


print('Filtering GRID POINTS')
print('\n')


dist_pp_grid_min = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_min[i] = [np.sqrt((j - pp_1_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_1_lat)]

dist_pp_grid_med = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_med[i] = [np.sqrt((j - pp_med_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_med_lat)]
    
dist_pp_grid_max = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_max[i] = [np.sqrt((j - pp_2_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_2_lat)]

grid_sel_min = []
grid_sel_min_data = []
for i,j in enumerate(dist_pp_grid_min):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[NUMBER_PP_PER_BIN]] < DIST_GRID_PP:
        grid_sel_min.append((grdx[i],grdy[i]))


grid_sel_med = []
grid_sel_med_data = []
for i,j in enumerate(dist_pp_grid_med):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[NUMBER_PP_PER_BIN]] < DIST_GRID_PP:
        grid_sel_med.append((grdx[i],grdy[i]))

        
grid_sel_max = []
grid_sel_min_data = []

for i,j in enumerate(dist_pp_grid_max):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[NUMBER_PP_PER_BIN]] < DIST_GRID_PP:
        grid_sel_max.append((grdx[i],grdy[i]))


grid_sel = grid_sel_min+grid_sel_med+grid_sel_max


grid_selected = set(map(tuple,grid_sel))

grid_sel_x = []
grid_sel_y = []

for i,j in enumerate(grid_selected):
    grid_sel_x.append(j[0])
    grid_sel_y.append(j[1])

grdx_1 = []
grdy_1 = []
grdx_sel_diff = []
for i,j in enumerate(grdx):
    grdx_sel_diff.append((j,grdy[i]))

grid_sel_diff = []
for i,j in enumerate(grid_selected):
	grid_sel_diff.append((j[0],j[1]))

grid_xy_diff = []
for i,j in enumerate(grdx_sel_diff):
	if j not in grid_sel_diff:
		grid_xy_diff.append(j)			
	
grid_x_diff = []
grid_y_diff = []
for i,j in enumerate(grid_xy_diff):
	grid_x_diff.append(j[0])
	grid_y_diff.append(j[1])

###################################################################################################################

print('Plotting: Figure Pds Piercing Points')
print('\n')


fig_PP, ax =  plt.subplots(nrows=1, ncols=1,figsize=(20,20))

#Figure Ppds

m_PP = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_LARGE,
            llcrnrlat=LLCRNRLAT_LARGE,urcrnrlon=URCRNRLON_LARGE,urcrnrlat=URCRNRLAT_LARGE,ax=ax)

m_PP.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_PP.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_PP(lon, lat)
    msize = 10
    l1, = m_PP.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')


for lon_1_Pds, lat_1_Pds in zip(pp_1_long,pp_1_lat):
    x_1_Pds,y_1_Pds = m_PP(lon_1_Pds, lat_1_Pds)
    msize_1 = 5
    l2, = m_PP.plot(x_1_Pds, y_1_Pds, '.',markersize=msize_1,markeredgecolor='k',markerfacecolor='b',alpha=0.5)


for lon_med_Pds, lat_med_Pds in zip(pp_med_long,pp_med_lat):
    x_med_Pds,y_med_Pds = m_PP(lon_med_Pds, lat_med_Pds)
    msize_1 = 5
    l3, = m_PP.plot(x_med_Pds, y_med_Pds, '.',markersize=msize_1,markeredgecolor='k',markerfacecolor='g',alpha=0.5)


for lon_2_Pds, lat_2_Pds in zip(pp_2_long,pp_2_lat):
    x_2_Pds,y_2_Pds = m_PP(lon_2_Pds, lat_2_Pds)
    msize_2 = 5
    l4, = m_PP.plot(x_2_Pds, y_2_Pds, '.',markersize=msize_2,markeredgecolor='k',markerfacecolor='r',alpha=0.5)



for lon_sel, lat_sel in zip(grid_sel_x,grid_sel_y):
    x_sel,y_sel = m_PP(lon_sel, lat_sel)
    msize_2 = 4
    l5, = m_PP.plot(x_sel, y_sel, 's',markersize=msize_2,markeredgecolor='k',markerfacecolor='None')

for lon_sel, lat_sel in zip(grid_x_diff,grid_y_diff):
    x_grdx,y_grdy = m_PP(lon_sel, lat_sel)
    msize_2 = 2
    l6, = m_PP.plot(x_grdx, y_grdy, '.',markersize=msize_2,markeredgecolor='k',markerfacecolor='None')

m_PP.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_PP.drawcoastlines(color='k',zorder=1)
m_PP.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_PP.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax.set_title('Pds Piercing Points',ha='center',va='top',y=1.08)
ax.legend([l1,l2,l3,l4,l5,l6],['Stations','Piercing Points 410 km','Piercing Points '+"{0:.0f}".format(DEPTH_MED)+' km','Piercing Points 660 km','Selected Grid', 'Raw Grid'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w',fontsize='smaller')

plt.show()

os.makedirs(PP_FIGURE,exist_ok=True)
fig_PP.savefig(PP_FIGURE+'PP_Pds.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting: Figure Pds Piercing Points')
print('\n')


fig_PP, ax =  plt.subplots(nrows=1, ncols=1,figsize=(20,20))

#Figure Ppds

m_PP = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_LARGE,
            llcrnrlat=LLCRNRLAT_LARGE,urcrnrlon=URCRNRLON_LARGE,urcrnrlat=URCRNRLAT_LARGE,ax=ax)

m_PP.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_PP.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_PP(lon, lat)
    msize = 10
    l1, = m_PP.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

for lon_med_Pds, lat_med_Pds in zip(pp_med_long,pp_med_lat):
    x_med_Pds,y_med_Pds = m_PP(lon_med_Pds, lat_med_Pds)
    msize_1 = 5
    l3, = m_PP.plot(x_med_Pds, y_med_Pds, 'X',markersize=msize_1,markeredgecolor='k',markerfacecolor='k',alpha=0.5)


for lon_sel, lat_sel in zip(grid_sel_x,grid_sel_y):
    x_sel,y_sel = m_PP(lon_sel, lat_sel)
    msize_2 = 4
    l5, = m_PP.plot(x_sel, y_sel, 's',markersize=msize_2,markeredgecolor='k',markerfacecolor='None')

for lon_sel, lat_sel in zip(grid_x_diff,grid_y_diff):
    x_grdx,y_grdy = m_PP(lon_sel, lat_sel)
    msize_2 = 2
    l6, = m_PP.plot(x_grdx, y_grdy, '.',markersize=msize_2,markeredgecolor='k',markerfacecolor='None')

m_PP.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_PP.drawcoastlines(color='k',zorder=1)
m_PP.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_PP.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax.set_title('Pds Piercing Points',ha='center',va='top',y=1.08)
ax.legend([l1,l3,l5,l6],['Stations','Piercing Points '+"{0:.0f}".format(DEPTH_MED)+' km','Selected Grid', 'Raw Grid'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w',fontsize='smaller')

plt.show()

os.makedirs(PP_FIGURE,exist_ok=True)
fig_PP.savefig(PP_FIGURE+'PP_MED_Pds.'+EXT_FIG,dpi=DPI_FIG)

##########################################################################################################################################

print('Importing depths and times of Pds conversion  dataset')
print('\n')

filename_Pds = PdS_DIR+'Pds_dic.json'

PdS_Dic = json.load(open(filename_Pds))

Pds_dist = []
Pds_time = []
Pds_depth = [] 
Pds_number = []

for i,j in enumerate(PdS_Dic):
	Pds_dist.append(j['dist'][0])
	Pds_time.append(j['time'][0])
	Pds_depth.append(j['depth'][0])
	Pds_number.append(j['number'][0])



###################################################################################################################
print('Migrating Pds dataset')
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


###################################################################################################################

print('Selecting the migrated data per grid data')
print('\n')

dados_grid_lat = pp_med_lat
dados_grid_lon = pp_med_long

RF_data_raw_Pds = [[]]*len(grid_sel_x)
RF_amplitude_depth_raw_Pds = [[]]*len(grid_sel_x)

for i,j in enumerate(grid_sel_x):
	RF_data_raw_Pds[i] = [RF_amplitude_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]
	RF_amplitude_depth_raw_Pds[i] = [RF_amplitude_depth_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]
###################################################################################################################


print('Estimating Mean and Standard Deviation for each discontinuity')
print('\n')


if BOOTSTRAP_DEPTH_ESTIMATION == True:
	### Creating Dictionaries to allocate results ###
	def nested_dict():
		return collections.defaultdict(nested_dict)

	RF_BOOTSTRAP_ESTIMATION_Pds = nested_dict()

	for _k in range(BOOTSTRAP_INTERATOR):
		print('Bootstrap estimation '+str(_k))

		for i,j in enumerate(RF_data_raw_Pds):
			if len(j) > NUMBER_PP_PER_BIN:
				print('Grid point number: '+str(i))
				print('lat: '+str(grid_sel_y[i]))
				print('lon: '+str(grid_sel_x[i]))

				old_RF_DATA_raw_Pds_lst = np.arange(0,len(j))
				new_RANDOM_RF_DATA_raw_Pds_lst = np.random.choice(old_RF_DATA_raw_Pds_lst,size=len(old_RF_DATA_raw_Pds_lst),replace=True)
				new_RANDOM_RF_DATA_raw_Pds = [j[_t_] for _t_ in new_RANDOM_RF_DATA_raw_Pds_lst]

				#410 km

				lst_410_depth_Pds = []
				lst_410_amp_Pds = []
				for k,l in enumerate(new_RANDOM_RF_DATA_raw_Pds):
					lst_depth_amp = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
					lst_depth_pp = [c for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
					lst_410_depth_Pds.append(lst_depth_pp[lst_depth_amp.index(max(lst_depth_amp))])	
					lst_410_amp_Pds.append(lst_depth_amp.index(max(lst_depth_amp)))
				######## Estimating mean and std depth ########

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_410_DEPTH'] = lst_410_depth_Pds

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lon'] = grid_sel_x[i]
				
				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lat'] = grid_sel_y[i]

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'] = np.average(lst_410_depth_Pds,weights=[1/i  if i != 0 else 0 for i in lst_410_amp_Pds])

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_std'] = np.std(lst_410_depth_Pds)

				print('410 Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'])+' ± '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_std']))

				#660 km

				lst_660_depth_Pds = []
				lst_660_amp_Pds = []
				for k,l in enumerate(new_RANDOM_RF_DATA_raw_Pds):
					lst_depth_amp = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
					lst_depth_pp = [c for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
					lst_660_depth_Pds.append(lst_depth_pp[lst_depth_amp.index(max(lst_depth_amp))])
					lst_660_amp_Pds.append(lst_depth_amp.index(max(lst_depth_amp)))
				######## Estimating mean and std depth ########

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_660_DEPTH'] = lst_660_depth_Pds

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'] = np.average(lst_660_depth_Pds,weights=[1/i if i != 0 else 0 for i in lst_660_amp_Pds])

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_std'] = np.std(lst_660_depth_Pds)

				print('660 Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'])+' ± '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_std']))
				print('\n')


			else:

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_410_DEPTH'] = 0

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lon'] = grid_sel_x[i]
				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lat'] = grid_sel_y[i]

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'] = 0

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_std'] = 0
				
				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_660_DEPTH'] = 0

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'] = 0

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_std'] = 0
	


#############################################################################################################################################################################################
#Allocating mean and std results:


RF_lat = []
RF_lon = []

RF_DEPTH_mean_1_Pds_lst = [[]]*len(RF_data_raw_Pds)
RF_DEPTH_std_1_Pds_lst = [[]]*len(RF_data_raw_Pds)
RF_DEPTH_mean_1_Pds = []
RF_DEPTH_std_1_Pds = []

RF_DEPTH_mean_2_Pds_lst = [[]]*len(RF_data_raw_Pds)
RF_DEPTH_std_2_Pds_lst = [[]]*len(RF_data_raw_Pds)
RF_DEPTH_mean_2_Pds = []
RF_DEPTH_std_2_Pds = []


for i,j in enumerate(RF_data_raw_Pds):
	if len(j) > NUMBER_PP_PER_BIN:
		RF_lat.append(RF_BOOTSTRAP_ESTIMATION_Pds[0][i]['lat'])
		RF_lon.append(RF_BOOTSTRAP_ESTIMATION_Pds[0][i]['lon'])

		RF_DEPTH_mean_1_Pds_lst[i].append([RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'] for _k in range(BOOTSTRAP_INTERATOR)])
		RF_DEPTH_std_1_Pds_lst[i].append([RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_std'] for _k in range(BOOTSTRAP_INTERATOR)])
		#RF_DEPTH_mean_1_Pds.append(np.mean(RF_DEPTH_mean_1_Pds_lst[i]))
		RF_DEPTH_std_1_Pds.append(np.mean(RF_DEPTH_std_1_Pds_lst[i]))

		RF_DEPTH_mean_2_Pds_lst[i].append([RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'] for _k in range(BOOTSTRAP_INTERATOR)])
		RF_DEPTH_std_2_Pds_lst[i].append([RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_std'] for _k in range(BOOTSTRAP_INTERATOR)])
		#RF_DEPTH_mean_2_Pds.append(np.mean(RF_DEPTH_mean_2_Pds_lst[i]))
		RF_DEPTH_std_2_Pds.append(np.mean(RF_DEPTH_std_2_Pds_lst[i]))

#############################################################################################################################################################################################

print('Stacking Pds data')
RF_stacking_Pds = []
len_RF_stacking_Pds = []

for i,j in enumerate(RF_data_raw_Pds):
	if len(j) > NUMBER_PP_PER_BIN:
		RF_stacking_Pds.append([sum(x)/len(j)  for x in zip(*j)])
		len_RF_stacking_Pds.append(len(j))


#############################################################################################################################################################################################
print('Calculating the Mean of 410 km and 660 km in Pds stacked data')


for i,j in enumerate(RF_stacking_Pds):
	lst_depth_410_amp = [j[x] for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
	lst_depth_410_pp = [c for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
	RF_DEPTH_mean_1_Pds.append(round(lst_depth_410_pp[lst_depth_410_amp.index(max(lst_depth_410_amp))],1))

	lst_depth_660_amp = [j[x] for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
	lst_depth_660_pp = [c for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
	RF_DEPTH_mean_2_Pds.append(round(lst_depth_660_pp[lst_depth_660_amp.index(max(lst_depth_660_amp))],1))

#############################################################################################################################################################################################

print('Plotting Figure: Depth of the Mantle Transition Zone for Pds phases (410 and 660 km)')
#Figure Depth of the Mantle Transition Zone for Pds and Ppds phases for 410 and 660 km

fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(10,20),squeeze=False,sharex=False,sharey=False)

ax1 = axes[0,0]
ax2 = axes[0,1]

colormap = 'seismic'

#Figure Depth of the Mantle Transition Zone for Pds phase for 410 km

m1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax1)

m1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m1(RF_lon,RF_lat)
sc1 = m1.scatter(x,y,40,RF_DEPTH_mean_1_Pds,cmap=colormap,marker='o',edgecolors='k',vmin=410-DEPTH_RANGE, vmax=410+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m1(lon, lat)
    msize = 10
    l1, = m1.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m1.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m1.drawcoastlines(color='k',zorder=1)
m1.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m1.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])


ax1.set_title('410 km Pds', y=1.08)

#Figure Depth of the Mantle Transition Zone for Pds phase for 660 km

m2 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax2)

m2.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m2.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m2(RF_lon,RF_lat)
sc2 = m2.scatter(x,y,40,RF_DEPTH_mean_2_Pds,cmap=colormap,marker='o',edgecolors='k',vmin=660-DEPTH_RANGE, vmax=660+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m2(lon, lat)
    msize = 10
    l1, = m2.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m2.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m2.drawcoastlines(color='k',zorder=1)
m2.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m2.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax2.set_title('660 km Pds', y=1.08)

fig.colorbar(sc1, ax=ax1,orientation='horizontal')
fig.colorbar(sc2, ax=ax2,orientation='horizontal')

fig.subplots_adjust(wspace=0.25, hspace=0.25)

fig.suptitle('Depth per bin')

plt.show()

fig.savefig(PP_FIGURE+'DEPTH_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting Figure: Std (bootstraping) of the Mantle Transition Zone for Pds phases (410 and 660 km)')


fig_std, axes_std = plt.subplots(nrows=1, ncols=2,figsize=(10,20),squeeze=False,sharex=False,sharey=False)

ax1_std = axes_std[0,0]
ax2_std = axes_std[0,1]

colormap_std = 'viridis'

#Figure Depth of the Mantle Transition Zone for Pds phase for 410 km

m1_std = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax1_std)

m1_std.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m1_std.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x1_std, y1_std = m1(RF_lon,RF_lat)
sc1_std = m1_std.scatter(x1_std,y1_std,40,RF_DEPTH_std_1_Pds,cmap=colormap_std,marker='o',edgecolors='k')


for lon1_std, lat1_std in zip(sta_long,sta_lat):
    x,y = m1_std(lon1_std, lat1_std)
    msize = 10
    l1, = m1_std.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m1_std.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m1_std.drawcoastlines(color='k',zorder=1)
m1_std.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m1_std.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax1_std.set_title('410 km Pds', y=1.08)
fig_std.colorbar(sc1_std, ax=ax1_std,orientation='horizontal')

#Figure Depth of the Mantle Transition Zone for Pds phase for 660 km


m2_std = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax2_std)

m2_std.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m2_std.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x2_std, y2_std = m2_std(RF_lon,RF_lat)
sc2_std = m2_std.scatter(x2_std,y2_std,40,RF_DEPTH_std_2_Pds,cmap=colormap_std,marker='o',edgecolors='k')

for lon, lat in zip(sta_long,sta_lat):
    x,y = m2_std(lon, lat)
    msize = 10
    l1, = m2_std.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m2_std.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m2_std.drawcoastlines(color='k',zorder=1)
m2_std.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m2_std.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax2_std.set_title('660 km Pds', y=1.08)
fig_std.colorbar(sc2_std, ax=ax2_std,orientation='horizontal')

fig_std.subplots_adjust(wspace=0.25, hspace=0.25)

fig_std.suptitle('Standard Deviation (bootstraping) per bin')

plt.show()

fig_std.savefig(PP_FIGURE+'TRUE_STD_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting Figure: Thickness of the Mantle Transition Zone...')
#Figure Thickness of the Mantle Transition Zone

thickness_MTZ_Pds = []


for i,j in enumerate(RF_DEPTH_mean_2_Pds):
	thickness_MTZ_Pds.append(j-RF_DEPTH_mean_1_Pds[i])

colormap_MTZ = 'seismic_r'

fig_thickness, ax_thickness1 = plt.subplots(nrows=1, ncols=1,figsize=(10,10))

m_thickness1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness1)

m_thickness1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness1(RF_lon,RF_lat)
sc_thickness1 = m_thickness1.scatter(x,y,40,thickness_MTZ_Pds,cmap=colormap_MTZ,marker='o',edgecolors='k',vmin=250-DEPTH_RANGE, vmax=250+DEPTH_RANGE)

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

plt.show()
fig_thickness.savefig(PP_FIGURE+'THICKNESS_MTZ_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Saving Selected Piercing Points in JSON file')
print('\n')

os.makedirs(PP_SELEC_DIR,exist_ok=True)

SELECTED_BINNED_DATA_dic = {'lat':[],'lon':[],'len_Pds':[],'mean_1_Pds':[],'std_1_Pds':[],'mean_2_Pds':[],'std_2_Pds':[],'data_Pds':[],'thickness_MTZ_Pds':[]}

for i,j in enumerate(RF_stacking_Pds):

	SELECTED_BINNED_DATA_dic['lat'].append(RF_lat[i])

	SELECTED_BINNED_DATA_dic['lon'].append(RF_lon[i])

	SELECTED_BINNED_DATA_dic['len_Pds'].append(len_RF_stacking_Pds[i])

	SELECTED_BINNED_DATA_dic['mean_1_Pds'].append(float(RF_DEPTH_mean_1_Pds[i]))

	SELECTED_BINNED_DATA_dic['std_1_Pds'].append(RF_DEPTH_std_1_Pds[i])

	SELECTED_BINNED_DATA_dic['mean_2_Pds'].append(float(RF_DEPTH_mean_2_Pds[i]))

	SELECTED_BINNED_DATA_dic['std_2_Pds'].append(RF_DEPTH_std_2_Pds[i])

	SELECTED_BINNED_DATA_dic['thickness_MTZ_Pds'].append(float(thickness_MTZ_Pds[i]))
	
	SELECTED_BINNED_DATA_dic['data_Pds'].append(j)


with open(PP_SELEC_DIR+'SELECTED_BINNED_Ps.json', 'w') as fp:
	json.dump(SELECTED_BINNED_DATA_dic, fp)