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
from matplotlib.collections import PolyCollection
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
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,STA_DIR,
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


print('Importing Ppds piercing points to each PHASE')
print('\n')

PHASES_Ppds = 'PPv410s','PPv'+"{0:.0f}".format(DEPTH_MED)+'s','PPv660s'

print('Importing Ppds Piercing Points '+PHASES_Ppds[0])
print('\n')

filename_1_Ppds = PP_DIR+'PP_'+PHASES_Ppds[0]+'_dic.json'

PP_1_dic_Ppds = json.load(open(filename_1_Ppds))

PP_dist_1_Ppds = []
PP_time_1_Ppds = []
PP_lat_1_Ppds = []
PP_lon_1_Ppds = []
PP_depth_1_Ppds = [] 
PP_1_number = [] 
for i,j in enumerate(PP_1_dic_Ppds):
	PP_dist_1_Ppds.append(j['dist'][0])
	PP_time_1_Ppds.append(j['time'][0])
	PP_lat_1_Ppds.append(j['lat'][0])
	PP_lon_1_Ppds.append(j['lon'][0])
	PP_depth_1_Ppds.append(j['depth'][0])
	PP_1_number.append(j['number'][0])


print('Importing Ppds Piercing Points '+PHASES_Ppds[1])
print('\n')

filename_med_Ppds = PP_DIR+'PP_'+PHASES_Ppds[1]+'_dic.json'

PP_med_dic_Ppds = json.load(open(filename_med_Ppds))

PP_dist_med_Ppds = []
PP_time_med_Ppds = []
PP_lat_med_Ppds = []
PP_lon_med_Ppds = []
PP_depth_med_Ppds = [] 
PP_med_number = [] 

for i,j in enumerate(PP_med_dic_Ppds):
	PP_dist_med_Ppds.append(j['dist'][0])
	PP_time_med_Ppds.append(j['time'][0])
	PP_lat_med_Ppds.append(j['lat'][0])
	PP_lon_med_Ppds.append(j['lon'][0])
	PP_depth_med_Ppds.append(j['depth'][0])
	PP_med_number.append(j['number'][0])


print('Importing Ppds Piercing Points '+PHASES_Ppds[2])
print('\n')

filename_2_Ppds = PP_DIR+'PP_'+PHASES_Ppds[2]+'_dic.json'

PP_2_dic_Ppds = json.load(open(filename_2_Ppds))

PP_dist_2_Ppds = []
PP_time_2_Ppds = []
PP_lat_2_Ppds = []
PP_lon_2_Ppds = []
PP_depth_2_Ppds = [] 
PP_2_number = [] 

for i,j in enumerate(PP_2_dic_Ppds):
	PP_dist_2_Ppds.append(j['dist'][0])
	PP_time_2_Ppds.append(j['time'][0])
	PP_lat_2_Ppds.append(j['lat'][0])
	PP_lon_2_Ppds.append(j['lon'][0])
	PP_depth_2_Ppds.append(j['depth'][0])
	PP_2_number.append(j['number'][0])

print('PPv410s Piercing Points')
print('\n')

pp_1_lat_Ppds  = [[]]*len(PP_lon_1_Ppds)
pp_1_long_Ppds  = [[]]*len(PP_lon_1_Ppds)


for i,j in enumerate(PP_lon_1_Ppds):
    for k,l in enumerate(j):
        if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_1_Ppds[i][k] == 410:
                pp_1_lat_Ppds[i] = PP_lat_1_Ppds[i][k] 
                pp_1_long_Ppds[i] = l



print('Ppds Piercing Points - '+"{0:.0f}".format(DEPTH_MED))
print('\n')

pp_med_lat_Ppds  = [[]]*len(PP_lon_med_Ppds)
pp_med_long_Ppds  = [[]]*len(PP_lon_med_Ppds)


for i,j in enumerate(PP_lon_med_Ppds):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_med_Ppds[i][k] == DEPTH_MED:
			pp_med_lat_Ppds[i] = PP_lat_med_Ppds[i][k]
			pp_med_long_Ppds[i] = l


print('PPv660s Piercing Points')
print('\n')


pp_2_lat_Ppds  = [[]]*len(PP_lon_2_Ppds)
pp_2_long_Ppds  = [[]]*len(PP_lon_2_Ppds)


for i,j in enumerate(PP_lon_2_Ppds):
    for k,l in enumerate(j):
        if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_2_Ppds[i][k] == 660:
                pp_2_lat_Ppds[i] = PP_lat_2_Ppds[i][k] 
                pp_2_long_Ppds[i] = l


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
grid_camadas_x = []
grid_camadas_y = []
grid_camadas_z = []
for i,j in enumerate(camadas_terra_10_km):
	for k,l in enumerate(grid_sel_x):
		grid_camadas_x.append(grid_sel_x[k])
		grid_camadas_y.append(grid_sel_y[k])
		grid_camadas_z.append(-j)

###################################################################################################################
print('Plotting: Figure earth model layers')
print('\n')

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.azim = 95
ax.elev = 10
ax.dist = 6

ax.plot(sta_long,sta_lat,'^',markersize=10,markeredgecolor='k',markerfacecolor='grey')
ax.plot(grid_sel_x,grid_sel_y,'.',markersize=2,markeredgecolor='k',markerfacecolor='k')

ax.scatter3D(grid_camadas_x, grid_camadas_y, grid_camadas_z, c='k', marker='o')
ax.set_zlim(-800,0)
ax.set_xlim(LLCRNRLON_SMALL,URCRNRLON_SMALL)
ax.set_ylim(LLCRNRLAT_SMALL,URCRNRLAT_SMALL)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Depth (km)')
plt.show()

###################################################################################################################

print('Plotting: Figure 410 and 660 Pds Piercing Points')
print('\n')

PP_410_lat = []
PP_410_lon = []
PP_410_depth = []
for i,j in enumerate(PP_depth_1):
	for k,l in enumerate(j):
			PP_410_lat.append(PP_lat_1[i][k])
			PP_410_lon.append(PP_lon_1[i][k])
			PP_410_depth.append(-PP_depth_1[i][k])

PP_660_lat = []
PP_660_lon = []
PP_660_depth = []
for i,j in enumerate(PP_depth_2):
	for k,l in enumerate(j):
			PP_660_lat.append(PP_lat_2[i][k])
			PP_660_lon.append(PP_lon_2[i][k])
			PP_660_depth.append(-PP_depth_2[i][k])


fig = plt.figure()
ax = plt.axes(projection='3d')

ax.azim = 95
ax.elev = 10
ax.dist = 8

ax.plot(sta_long,sta_lat,'^',markersize=10,markeredgecolor='k',markerfacecolor='grey')

x, y = np.meshgrid(np.arange(LLCRNRLON_SMALL,URCRNRLON_SMALL,1/GRID_PP_MULT),
                      np.arange(LLCRNRLAT_SMALL,URCRNRLAT_SMALL,1/GRID_PP_MULT))

z = np.zeros_like(x)
ax.plot_surface(x, y, z,color='None',edgecolor='k')

z1660 = np.ones_like(x)
z660 = np.array([element*-660 for element in z1660])

z1410 = np.ones_like(x)
z410 = np.array([element*-410 for element in z1410])

ax.plot_surface(x, y, z660,color='None',edgecolor='g')
ax.plot_surface(x, y, z410,color='None',edgecolor='b')

ax.scatter3D(PP_410_lon,PP_410_lat, PP_410_depth, c='k',marker='.',alpha=0.1)
ax.scatter3D(PP_660_lon,PP_660_lat, PP_660_depth, c='k',marker='.',alpha=0.1)

ax.set_zlim(-800,0)
ax.set_xlim(LLCRNRLON_SMALL,URCRNRLON_SMALL)
ax.set_ylim(LLCRNRLAT_SMALL,URCRNRLAT_SMALL)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Depth (km)')
plt.show()

###################################################################################################################

print('Plotting: Figure 410 Pds Piercing Points and Model earth')
print('\n')

PP_410_lat_POINTS = []
PP_410_lon_POINTS = []
PP_410_depth_POINTS = []
for i,j in enumerate(PP_depth_1):
	for k,l in enumerate(j):
		if l == 410:
			PP_410_lat_POINTS.append(PP_lat_1[i][k])
			PP_410_lon_POINTS.append(PP_lon_1[i][k])
			PP_410_depth_POINTS.append(-PP_depth_1[i][k])

PP_660_lat_POINTS = []
PP_660_lon_POINTS = []
PP_660_depth_POINTS = []
for i,j in enumerate(PP_depth_2):
	for k,l in enumerate(j):
		if l == 660:
			PP_660_lat_POINTS.append(PP_lat_2[i][k])
			PP_660_lon_POINTS.append(PP_lon_2[i][k])
			PP_660_depth_POINTS.append(-PP_depth_2[i][k])


fig = plt.figure(figsize=(30,15))
ax = Axes3D(fig)

ax.azim = 95
ax.elev = 10
ax.dist = 8

ax.plot(grid_sel_x,grid_sel_y,'.',markersize=2,markeredgecolor='k',markerfacecolor='k')
ax.plot(sta_long,sta_lat,'^',markersize=10,markeredgecolor='k',markerfacecolor='grey')


z = np.zeros_like(x)
ax.plot_surface(x, y, z,color='None',edgecolor='k')

ax.scatter3D(PP_410_lon_POINTS,PP_410_lat_POINTS, PP_410_depth_POINTS,  c='b',marker='+')
ax.scatter3D(PP_660_lon_POINTS,PP_660_lat_POINTS, PP_660_depth_POINTS,  c='g',marker='+')

m1 = ax.scatter3D(grid_camadas_x, grid_camadas_y, grid_camadas_z,c=grid_camadas_z,edgecolor='k',cmap='bone')


ax.set_zlim(-800,0)
ax.set_xlim(LLCRNRLON_SMALL,URCRNRLON_SMALL)
ax.set_ylim(LLCRNRLAT_SMALL,URCRNRLAT_SMALL)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Depth (km)')

fig.colorbar(m1,aspect=40,shrink=0.7)

plt.show()