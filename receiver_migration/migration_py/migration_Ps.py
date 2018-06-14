
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

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP
				   )


print('Looking for receiver functions data in JSON file in '+STA_DIR)

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

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

# Importing piercing points to P410s

filename_P410s = PP_DIR+'PP_P410s_dic.json'

PP_P410s_dic = json.load(open(filename_P410s))



PP_dist_P410s = PP_P410s_dic['dist']
PP_time_P410s = PP_P410s_dic['time']
PP_lat_P410s = PP_P410s_dic['lat']
PP_lon_P410s = PP_P410s_dic['lon']
PP_depth_P410s = PP_P410s_dic['depth']


# Importing piercing points to P660s

filename_P660s = PP_DIR+'PP_P660s_dic.json'

PP_P660s_dic = json.load(open(filename_P410s))



PP_dist_P660s = PP_P660s_dic['dist']
PP_time_P660s = PP_P660s_dic['time']
PP_lat_P660s = PP_P660s_dic['lat']
PP_lon_P660s = PP_P660s_dic['lon']
PP_depth_P660s = PP_P660s_dic['depth']

# Ray Trace PLOTS

if RAY_TRACE_PLOT == True:
	fig_ray_path, (ax1,ax2)= plt.subplots(nrows=1, ncols=2,figsize=(30, 10))

	    
	for i,j in enumerate(PP_depth_P410s):
	    ax1.plot(PP_lon_P410s[i],j,'g--',lw=0.5,alpha=0.5)
	ax1.set_title('Ps 410 km')

	for i,j in enumerate(PP_depth_P660s):
	    ax2.plot(PP_lon_P660s[i],j,'b--',lw=0.5,alpha=0.5)
	ax2.set_title('Ps 660 km')
	    
	for i,j in enumerate(camadas_terra_10_km):
	    ax1.hlines(j,LLCRNRLON_SMALL,URCRNRLON_SMALL,lw=0.5,alpha=0.5)
	    ax2.hlines(j,LLCRNRLON_SMALL,URCRNRLON_SMALL,lw=0.5,alpha=0.5)

	for i,j in enumerate(event_depth):
	    ax1.plot(sta_long[i],-20, color='w', marker='^', markersize=15,markeredgecolor='k')
	    ax2.plot(sta_long[i],-20, color='w', marker='^', markersize=15,markeredgecolor='k')

	    
	ax1.hlines(0,LLCRNRLON_SMALL,URCRNRLON_SMALL,colors='k',linestyles='solid')
	ax1.hlines(40,LLCRNRLON_SMALL,URCRNRLON_SMALL,colors='orange',linestyles='dashed',label='Moho')
	ax1.hlines(410,LLCRNRLON_SMALL,URCRNRLON_SMALL,colors='g',linestyles='dashed',label='410 km')
	ax1.hlines(660,LLCRNRLON_SMALL,URCRNRLON_SMALL,colors='b',linestyles='dashed',label='660 km')

	ax1.set_ylim(1000,-200)
	ax1.set_ylabel('Depth (km)')
	ax1.set_xlabel('Longitude')

	ax1.legend(loc=0,edgecolor='w',fancybox=True)
	ax1.set_xlim(LLCRNRLON_SMALL,URCRNRLON_SMALL)

	ax2.hlines(0,LLCRNRLON_SMALL,URCRNRLON_SMALL,colors='k',linestyles='solid')
	ax2.hlines(40,LLCRNRLON_SMALL,URCRNRLON_SMALL,colors='orange',linestyles='dashed',label='Moho')
	ax2.hlines(410,LLCRNRLON_SMALL,URCRNRLON_SMALL,colors='g',linestyles='dashed',label='410 km')
	ax2.hlines(660,LLCRNRLON_SMALL,URCRNRLON_SMALL,colors='b',linestyles='dashed',label='660 km')

	ax2.set_ylim(1000,-200)
	ax2.set_ylabel('Depth (km)')
	ax2.set_xlabel('Longitude')

	ax2.legend(loc=0,edgecolor='w',fancybox=True)
	ax2.set_xlim(max(sta_long),URCRNRLON_SMALL)

	os.makedirs(RAY_PATH_FIGURE,exist_ok=True)
	fig_ray_path.savefig(RAY_PATH_FIGURE+'ray_path.'+EXT_FIG,dpi=DPI_FIG)
else: 
	pass

if RAY_TRACE_410_660_PLOT == True:
	fig_ray_path_410_660, (ax1,ax2)= plt.subplots(nrows=1, ncols=2,figsize=(30, 10))

	   
	for i,j in enumerate(PP_depth_P410s):
	    ax1.plot(PP_time_P410s[i],j,'g--',lw=0.5,alpha=0.5)
	ax1.set_title('Ps 410 km')

	for i,j in enumerate(PP_depth_P660s):
	    ax2.plot(PP_time_P660s[i],j,'b--',lw=0.5,alpha=0.5)
	ax2.set_title('Ps 660 km')
	    
	for i,j in enumerate(camadas_terra_10_km):
	    ax1.hlines(j,0,1000,lw=0.5,alpha=0.5)
	    ax2.hlines(j,0,1000,lw=0.5,alpha=0.5)

	ax1.hlines(0,200,850,colors='k',linestyles='solid')
	ax1.hlines(40,200,850,colors='orange',linestyles='dashed',label='Moho')
	ax1.hlines(410,200,850,colors='g',linestyles='dashed',label='410 km')
	ax1.hlines(660,200,850,colors='b',linestyles='dashed',label='660 km')

	ax1.set_xlim(200,850)
	ax1.set_ylim(1000,0)
	ax1.set_ylabel('Depth (km)')
	ax1.set_xlabel('Time (s)')

	ax1.legend(loc=0,edgecolor='w',fancybox=True)

	ax2.hlines(0,200,850,colors='k',linestyles='solid')
	ax2.hlines(40,200,850,colors='orange',linestyles='dashed',label='Moho')
	ax2.hlines(410,200,850,colors='g',linestyles='dashed',label='410 km')
	ax2.hlines(660,200,850,colors='b',linestyles='dashed',label='660 km')

	ax2.set_xlim(200,850)
	ax2.set_ylim(1000,0)
	ax2.set_ylabel('Depth (km)')
	ax2.set_xlabel('Time (s)')

	ax2.legend(loc=0,edgecolor='w',fancybox=True)

	fig_ray_path_410_660.savefig(RAY_PATH_FIGURE+'ray_path_410_660.'+EXT_FIG,dpi=DPI_FIG)
else: 
	pass


# Piercing Points -410 km

pp_410_lat  = [[]]*len(PP_lon_P410s)
pp_410_long  = [[]]*len(PP_lon_P410s)


for i,j in enumerate(PP_lon_P410s):
    for k,l in enumerate(j):
        if LLCRNRLON_SMALL<= l <= URCRNRLON_SMALL and PP_depth_P410s[i][k] == 410:
                pp_410_lat[i] = PP_lat_P410s[i][k] 
                pp_410_long[i] = l


# Piercing Points -660 km

pp_660_lat  = [[]]*len(PP_lon_P660s)
pp_660_long  = [[]]*len(PP_lon_P660s)


for i,j in enumerate(PP_lon_P660s):
    for k,l in enumerate(j):
        if LLCRNRLON_SMALL <= l <= URCRNRLON_SMALL and PP_depth_P660s[i][k] == 660:
                pp_660_lat[i] = PP_lat_P660s[i][k] 
                pp_660_long[i] = l


print('Creating GRID POINTS')

area = (LLCRNRLON_SMALL,URCRNRLON_SMALL, LLCRNRLAT_SMALL, URCRNRLAT_SMALL)

shape = (abs(abs(URCRNRLON_SMALL) - abs(LLCRNRLON_SMALL))*3, abs(abs(URCRNRLAT_SMALL) - abs(LLCRNRLAT_SMALL))*3)

# First, we need to know the real data at the grid points
grdx, grdy = gridder.regular(area, shape)

fig_grid_points=plt.figure(figsize=(20,10))

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_LARGE,
            llcrnrlat=LLCRNRLAT_LARGE,urcrnrlon=URCRNRLON_LARGE,urcrnrlat=URCRNRLAT_LARGE)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m(lon, lat)
    msize = 10
    m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

    
for lon, lat in zip(pp_410_long,pp_410_lat):
    x,y = m(lon, lat)
    msize = 3
    m.plot(x, y, '+',markersize=msize,markeredgecolor='r',markerfacecolor='k')
    
for lon, lat in zip(pp_660_long,pp_660_lat):
    x,y = m(lon, lat)
    msize = 3
    m.plot(x, y, '+',markersize=msize,markeredgecolor='b',markerfacecolor='k')
    
for lon, lat in zip(grdx,grdy):
    x,y = m(lon, lat)
    msize = 5
    m.plot(x, y, '.',markersize=msize,markeredgecolor='k',markerfacecolor="None")

m.fillcontinents(color='whitesmoke',lake_color=None)
m.drawcoastlines(color='dimgray',zorder=10)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
plt.title('BINNED DATA', y=1.08)

os.makedirs(PP_FIGURE,exist_ok=True)
fig_grid_points.savefig(PP_FIGURE+'BINNED_DATA.'+EXT_FIG,dpi=DPI_FIG)


# Filtering grid points


dist_pp_grid_min = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_min[i] = [np.sqrt((j - pp_410_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_410_lat)]
    
dist_pp_grid_max = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_max[i] = [np.sqrt((j - pp_660_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_660_lat)]



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


fig_grid_points_filtered=plt.figure(figsize=(20,10))

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

for lon, lat in zip(grdx,grdy):
    x,y = m(lon, lat)
    msize = 0.5
    l5, = m.plot(x, y, 'o',markersize=msize,markeredgecolor='k',markerfacecolor="None")

for lon, lat in zip(grid_sel_x, grid_sel_y):
    x,y = m(lon, lat)
    msize = 4
    l4, = m.plot(x, y, 'o',markersize=msize,markeredgecolor='k',markerfacecolor="None")
    
for lon, lat in zip(pp_410_long,pp_410_lat):
    x,y = m(lon, lat)
    msize = 3
    l3, = m.plot(x, y, '+',markersize=msize,markeredgecolor='r',markerfacecolor='k')
    
for lon, lat in zip(pp_660_long,pp_660_lat):
    x,y = m(lon, lat)
    msize = 3
    l2, = m.plot(x, y, '+',markersize=msize,markeredgecolor='b',markerfacecolor='k')
    

for lon, lat in zip(sta_long,sta_lat):
    x,y = m(lon, lat)
    msize = 10
    l1, = m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m.fillcontinents(color='whitesmoke',lake_color=None)
m.drawcoastlines(color='dimgray',zorder=10)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

        
label=['Stations','Piercing Points 410','Piercing Points 660','Filtered Grid Points','Grid Points']
plt.legend([l1,l2,l3,l4,l5],label,scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w')

plt.title('SELECTED BINNED DATA', y=1.08)
fig_grid_points_filtered.savefig(PP_FIGURE+'SELECTED_BINNED_DATA.'+EXT_FIG,dpi=DPI_FIG)


print('Importing depths and times to the Ps conversion to each event for all stations')

filename_Pds = PdS_DIR+'Pds_dic.json'

PdS_Dic = json.load(open(filename_Pds))

P_dist = PdS_Dic['dist']
P_time = PdS_Dic['time']
P_depth = PdS_Dic['depth']


print('Migrating data...')

RF_amplitude_time = [[]]*len(P_depth)
for i,j in enumerate(P_depth):
    sta_t = j
    RF_t = camadas_terra_10_km
    RF_amplitude_time[i] = [P_time[i][sta_t.index(l)] for k,l in enumerate(RF_t)]


RF_amplitude = []
for i,j in enumerate(RF_amplitude_time):
    sta_t = [round(l,1) for k,l in enumerate(sta_time[i])]
    RF_t = [round(l,1) for k,l in enumerate(j)]
    RF_amplitude.append([sta_data[i][sta_t.index(l)] for k,l in enumerate(RF_t)])


RF_amplitude_time = [[]]*len(P_depth)
for i,j in enumerate(P_depth):
    sta_t = j
    RF_t = camadas_terra_10_km
    RF_amplitude_time[i] = [P_time[i][sta_t.index(l)] if l in sta_t else -1 for k,l in enumerate(RF_t)]


RF_amplitude = [[]]*len(RF_amplitude_time)
for i,j in enumerate(RF_amplitude_time):
    sta_t = [round(l,1) for k,l in enumerate(sta_time[i])]
    RF_t = [round(l,1) for k,l in enumerate(j)]
    RF_amplitude[i] = [sta_data[i][sta_t.index(l)] if l != -1 else 0 for k,l in enumerate(RF_t)]


# Data stacking in each point of the filtered grid

dist_between_grid_piercing_points = DIST_GRID_PP
dados_grid_lat = pp_410_lat
dados_grid_lon = pp_410_long


RF_data_raw = [[]]*len(grid_sel_x)
RF_RAY_raw = [[]]*len(grid_sel_x)

for i,j in enumerate(grid_sel_x):
    RF_data_raw[i] = [RF_amplitude[k]  for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < dist_between_grid_piercing_points]
    RF_RAY_raw[i] = [event_ray[k]  for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < dist_between_grid_piercing_points]


RF_stacking = []
len_RF_stacking = []

for i,j in enumerate(RF_data_raw):
    if len(j) > NUMBER_PP_PER_BIN:
        RF_stacking.append([sum(x)/len(j)  for x in zip(*j)])
        len_RF_stacking.append(len(j))
    else:
        RF_stacking.append([])
        len_RF_stacking.append(0)


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

#Saving Selected Piercing Points in JSON file

os.makedirs(PP_SELEC_DIR,exist_ok=True)

SELECTED_BINNED_DATA_dic = {'lat':[],'lon':[],'len':[],'data':[]}
for i,j in enumerate(RF_stacking):
	SELECTED_BINNED_DATA_dic['lat'].append(grid_sel_y[i])
	SELECTED_BINNED_DATA_dic['lon'].append(grid_sel_x[i])
	SELECTED_BINNED_DATA_dic['len'].append(len_RF_stacking[i])
	SELECTED_BINNED_DATA_dic['data'].append(j)

with open(PP_SELEC_DIR+'SELECTED_BINNED.json', 'w') as fp:
	json.dump(SELECTED_BINNED_DATA_dic, fp)
