	
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
from shapely.geometry import Polygon, MultiPoint, Point
import shapefile





from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP_MED,DIST_GRID_PP,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,BOOTSTRAP_DEPTH_ESTIMATION,ROTATE_GRID,ROTATE_ANGLE,GAMMA
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

if ROTATE_GRID == True:

	radians = ROTATE_ANGLE*np.pi/180
	grdx_rot = []
	grdy_rot = []

	for i,j in enumerate(grdx):

		x = grdx[i]
		y = grdy[i]

		ox = URCRNRLON_LARGE
		oy = URCRNRLAT_LARGE

		qx = ox + np.cos(radians) * (x - ox) + np.sin(radians) * (y - oy)
		qy = oy + -np.sin(radians) * (x - ox) + np.cos(radians) * (y - oy)

		grdx_rot.append(float(qx))
		grdy_rot.append(float(qy))

	grdx = grdx_rot
	grdy = grdy_rot

if FILTER_BY_SHAPEFILE == True:
	polygon = shapefile.Reader(SHAPEFILE_GRID) 
	polygon = polygon.shapes()  
	shpfilePoints = []
	for shape in polygon:
		shpfilePoints = shape.points 
	polygon = shpfilePoints 
	poly = Polygon(polygon)
	
	pontos = [Point(grdx[i],grdy[i]) for i,j in enumerate(grdx)]

	inside_points = np.array(MultiPoint(pontos).intersection(poly))

	grdx = []
	grdy = []
	for i,j in enumerate(inside_points):
		grdx.append(j[0])
		grdy.append(j[1])
	

print('Filtering GRID POINTS')
print('\n')


dist_pp_grid_min = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_min[i] = [np.sqrt((j - pp_1_long_Ppds[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_1_lat_Ppds)]

dist_pp_grid_med = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_med[i] = [np.sqrt((j - pp_med_long_Ppds[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_med_lat_Ppds)]
    
dist_pp_grid_max = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_max[i] = [np.sqrt((j - pp_2_long_Ppds[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_2_lat_Ppds)]

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
print('Plotting: Figure Pds and Ppds Piercing Points')
print('\n')


fig_PP, (ax, ax1) =  plt.subplots(nrows=1, ncols=2,figsize=(20,20))

#Figure Pds

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
    l2, = m_PP.plot(x_1_Pds, y_1_Pds, '.',markersize=msize_1,markeredgecolor='k',markerfacecolor='b')


for lon_med_Pds, lat_med_Pds in zip(pp_med_long,pp_med_lat):
    x_med_Pds,y_med_Pds = m_PP(lon_med_Pds, lat_med_Pds)
    msize_1 = 5
    l3, = m_PP.plot(x_med_Pds, y_med_Pds, '.',markersize=msize_1,markeredgecolor='k',markerfacecolor='g')


for lon_2_Pds, lat_2_Pds in zip(pp_2_long,pp_2_lat):
    x_2_Pds,y_2_Pds = m_PP(lon_2_Pds, lat_2_Pds)
    msize_2 = 5
    l4, = m_PP.plot(x_2_Pds, y_2_Pds, '.',markersize=msize_2,markeredgecolor='k',markerfacecolor='r')



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

#Figure Ppds

m_PP1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_LARGE,
            llcrnrlat=LLCRNRLAT_LARGE,urcrnrlon=URCRNRLON_LARGE,urcrnrlat=URCRNRLAT_LARGE,ax=ax1)

m_PP1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_PP1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_PP1(lon, lat)
    msize = 10
    l6, = m_PP1.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')


for lon_1_Ppds, lat_1_Ppds in zip(pp_1_long_Ppds,pp_1_lat_Ppds):
    x_1_Ppds,y_1_Ppds = m_PP1(lon_1_Ppds, lat_1_Ppds)
    msize = 5
    l7, = m_PP1.plot(x_1_Ppds, y_1_Ppds, '.',markersize=msize,markeredgecolor='k',markerfacecolor='b')

for lon_med_Ppds, lat_med_Ppds in zip(pp_med_long_Ppds,pp_med_lat_Ppds):
    x_med_Ppds,y_med_Ppds = m_PP1(lon_med_Ppds, lat_med_Ppds)
    msize = 5
    l8, = m_PP1.plot(x_med_Ppds, y_med_Ppds, '.',markersize=msize,markeredgecolor='k',markerfacecolor='g')

for lon_2_Ppds, lat_2_Ppds in zip(pp_2_long_Ppds,pp_2_lat_Ppds):
    x_2_Ppds,y_2_Ppds = m_PP1(lon_2_Ppds, lat_2_Ppds)
    msize_2 = 5
    l9, = m_PP1.plot(x_2_Ppds, y_2_Ppds, '.',markersize=msize_2,markeredgecolor='k',markerfacecolor='r')

for lon_sel, lat_sel in zip(grid_sel_x,grid_sel_y):
    x_sel,y_sel = m_PP1(lon_sel, lat_sel)
    msize_2 = 4
    l10, = m_PP1.plot(x_sel, y_sel, 's',markersize=msize_2,markeredgecolor='k',markerfacecolor='None')

for lon_sel, lat_sel in zip(grid_x_diff,grid_y_diff):
    x_grdx,y_grdy = m_PP1(lon_sel, lat_sel)
    msize_2 = 2
    l11, = m_PP1.plot(x_grdx, y_grdy, '.',markersize=msize_2,markeredgecolor='k',markerfacecolor='None')

m_PP1.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_PP1.drawcoastlines(color='k',zorder=1)
m_PP1.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_PP1.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax1.set_title('Ppds Piercing Points',ha='center',va='top',y=1.08)
ax1.legend([l6,l7,l8,l9,l10,l11],['Stations','Piercing Points 410 km','Piercing Points '+"{0:.0f}".format(DEPTH_MED)+' km','Piercing Points 660 km','Selected Grid','Raw Grid'],scatterpoints=1, frameon=True,labelspacing=0.5, loc='lower right',facecolor='w',fontsize='smaller')


plt.show()

os.makedirs(PP_FIGURE,exist_ok=True)
fig_PP.savefig(PP_FIGURE+'PP_Pds_Ppds.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting: Figure Pds Average Piercing Points')
print('\n')


fig_PP, ax =  plt.subplots(nrows=1, ncols=1,figsize=(20,20))

#Figure Pds

m_PP = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_LARGE,
            llcrnrlat=LLCRNRLAT_LARGE,urcrnrlon=URCRNRLON_LARGE,urcrnrlat=URCRNRLAT_LARGE,ax=ax)

m_PP.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_PP.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_PP(lon, lat)
    msize = 10
    l1, = m_PP.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

for lon_med_Pds, lat_med_Pds in zip(pp_med_long_Ppds,pp_med_lat_Ppds):
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

###################################################################################################################

print('Plotting: Figure Ppds Average Piercing Points')
print('\n')


fig_PP, ax =  plt.subplots(nrows=1, ncols=1,figsize=(20,20))

#Figure Pds

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

ax.set_title('Ppds Piercing Points',ha='center',va='top',y=1.08)
ax.legend([l1,l3,l5,l6],['Stations','Piercing Points '+"{0:.0f}".format(DEPTH_MED)+' km','Selected Grid', 'Raw Grid'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w',fontsize='smaller')

plt.show()

os.makedirs(PP_FIGURE,exist_ok=True)
fig_PP.savefig(PP_FIGURE+'PP_MED_Ppds.'+EXT_FIG,dpi=DPI_FIG)

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



print('Importing depths and times of Ppds conversion  dataset')
print('\n')

filename_Ppds = PdS_DIR+'PPvs_dic.json'

Ppds_Dic = json.load(open(filename_Ppds))

Ppds_dist = []
Ppds_time = []
Ppds_depth = [] 
Ppds_number = []

for i,j in enumerate(Ppds_Dic):
	Ppds_dist.append(j['dist'][0])
	Ppds_time.append(j['time'][0])
	Ppds_depth.append(j['depth'][0])
	Ppds_number.append(j['number'][0])


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


print('Migrating Ppds dataset')
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




###################################################################################################################

print('Selecting the migrated data per grid data')
print('\n')

dados_grid_lat = pp_med_lat
dados_grid_lon = pp_med_long

RF_data_raw_Pds = [[]]*len(grid_sel_x)
RF_amplitude_depth_raw_Pds = [[]]*len(grid_sel_x)

RF_data_raw_Ppds = [[]]*len(grid_sel_x)
RF_amplitude_depth_raw_Ppds = [[]]*len(grid_sel_x)

for i,j in enumerate(grid_sel_x):
	RF_data_raw_Pds[i] = [RF_amplitude_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]
	RF_amplitude_depth_raw_Pds[i] = [RF_amplitude_depth_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]

	RF_data_raw_Ppds[i] = [RF_amplitude_Ppds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]
	RF_amplitude_depth_raw_Ppds[i] = [RF_amplitude_depth_Ppds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) < DIST_GRID_PP_MED]

###################################################################################################################


print('Estimating Mean and Standard Deviation for each discontinuity')
print('\n')


if BOOTSTRAP_DEPTH_ESTIMATION == True:
	### Creating Dictionaries to allocate results ###
	def nested_dict():
		return collections.defaultdict(nested_dict)

	RF_BOOTSTRAP_ESTIMATION_Pds = nested_dict()
	RF_BOOTSTRAP_ESTIMATION_Ppds = nested_dict()


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
				new_RANDOM_RF_DATA_raw_Ppds = [RF_data_raw_Ppds[i][_t_] for _t_ in new_RANDOM_RF_DATA_raw_Pds_lst]

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_DATA'] = new_RANDOM_RF_DATA_raw_Pds
				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['RF_DATA'] = new_RANDOM_RF_DATA_raw_Ppds
				


				lst_410_depth_Pds = []
				lst_410_amp_Pds = []

				lst_660_depth_Pds = []
				lst_660_amp_Pds = []
				for k,l in enumerate(new_RANDOM_RF_DATA_raw_Pds):

					#410 km

					lst_depth_amp_410_Pds = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
					lst_depth_pp_410_Pds = [c for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
					lst_410_depth_Pds.append(lst_depth_pp_410_Pds[lst_depth_amp_410_Pds.index(max(lst_depth_amp_410_Pds))])	
					lst_410_amp_Pds.append(lst_depth_amp_410_Pds.index(max(lst_depth_amp_410_Pds)))

					#660 km

					lst_depth_amp_660_Pds = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
					lst_depth_pp_660_Pds = [c for x,c in enumerate(RF_amplitude_depth_raw_Pds[i][k]) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
					lst_660_depth_Pds.append(lst_depth_pp_660_Pds[lst_depth_amp_660_Pds.index(max(lst_depth_amp_660_Pds))])
					lst_660_amp_Pds.append(lst_depth_amp_660_Pds.index(max(lst_depth_amp_660_Pds)))


				lst_410_depth_Ppds = []
				lst_410_amp_Ppds = []

				lst_660_depth_Ppds = []
				lst_660_amp_Ppds = []
				for k,l in enumerate(new_RANDOM_RF_DATA_raw_Ppds):

					#410 km

					lst_depth_amp_410_Ppds = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Ppds[i][k]) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
					lst_depth_pp_410_Ppds = [c for x,c in enumerate(RF_amplitude_depth_raw_Ppds[i][k]) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
					lst_410_depth_Ppds.append(lst_depth_pp_410_Ppds[lst_depth_amp_410_Ppds.index(max(lst_depth_amp_410_Ppds))])
					lst_410_amp_Ppds.append(lst_depth_amp_410_Ppds.index(max(lst_depth_amp_410_Ppds)))

					#660 km

					lst_depth_amp_660_Ppds = [l[x] for x,c in enumerate(RF_amplitude_depth_raw_Ppds[i][k]) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
					lst_depth_pp_660_Ppds = [c for x,c in enumerate(RF_amplitude_depth_raw_Ppds[i][k]) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
					lst_660_depth_Ppds.append(lst_depth_pp_660_Ppds[lst_depth_amp_660_Ppds.index(max(lst_depth_amp_660_Ppds))])
					lst_660_amp_Ppds.append(lst_depth_amp_660_Ppds.index(max(lst_depth_amp_660_Ppds)))

				######## Estimating mean and std depth ########

				#410 km


				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_410_DEPTH'] = lst_410_depth_Pds
				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['RF_410_DEPTH'] = lst_410_depth_Ppds

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lon'] = grid_sel_x[i]
				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lat'] = grid_sel_y[i]

				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['lon'] = grid_sel_x[i]
				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['lat'] = grid_sel_y[i]

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'] = np.average(lst_410_depth_Pds,weights=[100 if i > 0.005 else 5 for i in lst_410_amp_Pds])
				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_std'] = np.std(lst_410_depth_Pds)

				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_mean'] = np.average(lst_410_depth_Ppds,weights=[100 if i > 0.005 else 5 for i in lst_410_amp_Ppds])
				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_std'] = np.std(lst_410_depth_Ppds)

				print('410 Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'])+' ± '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_std']))
				print('410 Ppds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_mean'])+' ± '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_std']))

				######## Estimating TRUE depth ########

				Vp_Vs_ratio_depth_1 = Vs_depth_1/Vp_depth_1
				alfa = Vp_depth_1 - Vs_depth_1
				beta = Vp_depth_1 + Vs_depth_1
				gamma_vp_vs_1 = GAMMA*Vp_Vs_ratio_depth_1

				H_a_Pds = np.average(lst_410_depth_Pds,weights=[100 if i > 0.005 else 5 for i in lst_410_amp_Pds])
				H_a_Ppds = np.average(lst_410_depth_Ppds,weights=[100 if i > 0.005 else 5 for i in lst_410_amp_Ppds])

				lst_delta_Vp_410_Pds = (alfa*beta*(H_a_Ppds - H_a_Pds))/(((alfa*H_a_Pds)+(alfa*H_a_Pds*gamma_vp_vs_1)) - ((beta*H_a_Ppds)-(beta*H_a_Ppds*gamma_vp_vs_1)))

				lst_delta_Vs_410_Pds = lst_delta_Vp_410_Pds * GAMMA * Vp_Vs_ratio_depth_1

				lst_true_410_mean_Pds = H_a_Pds* ((Vp_depth_1-Vs_depth_1)/(Vp_depth_1*Vs_depth_1)) * ((Vs_depth_1+lst_delta_Vs_410_Pds)*(Vp_depth_1+lst_delta_Vp_410_Pds))/(Vp_depth_1+lst_delta_Vp_410_Pds-Vs_depth_1-lst_delta_Vs_410_Pds)

				lst_true_410_mean_Ppds = H_a_Ppds * ((Vp_depth_1+Vs_depth_1)/(Vp_depth_1*Vs_depth_1)) * ((Vs_depth_1+lst_delta_Vs_410_Pds)*(Vp_depth_1+lst_delta_Vp_410_Pds))/(Vp_depth_1+lst_delta_Vp_410_Pds+Vs_depth_1+lst_delta_Vs_410_Pds)
					
				
				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_410_mean'] = lst_delta_Vp_410_Pds

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_410_mean'] = lst_delta_Vs_410_Pds

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_410_mean'] = lst_true_410_mean_Pds

				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_410_mean'] = lst_true_410_mean_Ppds

				print('Delta Vp 410 km = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_410_mean']))
				print('Delta Vs 410 km = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_410_mean']))				
				
				print('410 km Pds True Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_410_mean']))
				print('410 km Ppds True Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_410_mean']))

				######## Estimating mean and std depth ########

				#660 km

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_660_DEPTH'] = lst_660_depth_Pds
				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['RF_660_DEPTH'] = lst_660_depth_Ppds

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'] = np.average(lst_660_depth_Pds,weights=[100 if i > 0.005 else 5 for i in lst_660_amp_Pds])
				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_std'] = np.std(lst_660_depth_Pds)

				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_mean'] = np.average(lst_660_depth_Ppds,weights=[100 if i > 0.005 else 5 for i in lst_660_amp_Ppds])
				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_std'] = np.std(lst_660_depth_Ppds)

				print('660 km Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'])+' ± '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_std']))
				print('660 km Ppds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_mean'])+' ± '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_std']))



				######## Estimating TRUE depth ########

				Vp_Vs_ratio_depth_2 = Vs_depth_2/Vp_depth_2
				alfa = Vp_depth_2 - Vs_depth_2
				beta = Vp_depth_2 + Vs_depth_2
				gamma_vp_vs_2 = GAMMA*Vp_Vs_ratio_depth_2

				H_a_Pds_660 = np.average(lst_660_depth_Pds,weights=[100 if i > 0.005 else 5 for i in lst_660_amp_Pds])
				H_a_Ppds_660 = np.average(lst_660_depth_Ppds,weights=[100 if i > 0.005 else 5 for i in lst_660_amp_Ppds])


				lst_delta_Vp_660_Pds = (alfa*beta*(H_a_Ppds_660 - H_a_Pds_660))/(((alfa*H_a_Pds_660)+(alfa*H_a_Pds_660*gamma_vp_vs_1)) - ((beta*H_a_Ppds_660)-(beta*H_a_Ppds_660*gamma_vp_vs_1)))

				lst_delta_Vs_660_Pds = lst_delta_Vp_660_Pds * GAMMA * Vp_Vs_ratio_depth_2

				lst_true_660_mean_Pds = H_a_Pds_660* ((Vp_depth_2-Vs_depth_2)/(Vp_depth_2*Vs_depth_2)) * ((Vs_depth_2+lst_delta_Vs_660_Pds)*(Vp_depth_2+lst_delta_Vp_660_Pds))/(Vp_depth_2+lst_delta_Vp_660_Pds-Vs_depth_2-lst_delta_Vs_660_Pds)

				lst_true_660_mean_Ppds = H_a_Ppds_660 * ((Vp_depth_2+Vs_depth_2)/(Vp_depth_2*Vs_depth_2)) * ((Vs_depth_2+lst_delta_Vs_660_Pds)*(Vp_depth_2+lst_delta_Vp_660_Pds))/(Vp_depth_2+lst_delta_Vp_660_Pds+Vs_depth_2+lst_delta_Vs_660_Pds)
	
				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_660_mean'] = lst_delta_Vp_660_Pds

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_660_mean'] = lst_delta_Vs_660_Pds

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_660_mean'] = lst_true_660_mean_Pds

				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_660_mean'] = lst_true_660_mean_Ppds

				print('Delta Vp 660 km = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_660_mean']))
				print('Delta Vs 660 km = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_660_mean']))

				print('660 km Pds True Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_660_mean']))
				print('660 km Ppds True Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_660_mean']))

				######## Estimating MTZ thickness ########

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean'] = np.mean(lst_660_depth_Pds) - np.mean(lst_410_depth_Pds)

				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['thickness_MTZ_mean'] = np.mean(lst_660_depth_Ppds) - np.mean(lst_410_depth_Ppds)

				print('MTZ Pds thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean']))
				print('MTZ Ppds thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['thickness_MTZ_mean']))


				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_thickness_MTZ_mean'] = np.mean(lst_true_660_mean_Pds) - np.mean(lst_true_410_mean_Pds)

				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_thickness_MTZ_mean'] = np.mean(lst_true_660_mean_Ppds) - np.mean(lst_true_410_mean_Ppds)

				print('MTZ Pds true thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_thickness_MTZ_mean']))
				print('MTZ Ppds true thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_thickness_MTZ_mean']))

				######## Estimating MTZ difference thickness ########

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ'] = (np.mean(lst_true_660_mean_Pds) - np.mean(lst_true_410_mean_Pds))  - (np.mean(lst_660_depth_Pds) - np.mean(lst_410_depth_Pds))
				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ'] = (np.mean(lst_true_660_mean_Ppds) - np.mean(lst_true_410_mean_Ppds)) - (np.mean(lst_660_depth_Ppds) - np.mean(lst_410_depth_Ppds))

				print('MTZ Pds diff thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ']))
				print('MTZ Ppds diff thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ']))
				print('\n')


#############################################################################################################################################################################################
#Allocating mean and std results:


RF_lat = []
RF_lon = []

RF_DEPTH_mean_1_Pds = []
RF_DEPTH_std_1_Pds = []
RF_DEPTH_mean_1_Ppds = []
RF_DEPTH_std_1_Ppds = []

delta_1_Vp_mean = []
delta_1_Vp_std = []
delta_1_Vs_mean = []
delta_1_Vs_std = []

RF_DEPTH_mean_1_true_Pds = []
RF_DEPTH_std_1_true_Pds = []
RF_DEPTH_mean_1_true_Ppds = []
RF_DEPTH_std_1_true_Ppds = []

RF_DEPTH_mean_2_Pds = []
RF_DEPTH_std_2_Pds = []
RF_DEPTH_mean_2_Ppds = []
RF_DEPTH_std_2_Ppds = []

delta_2_Vp_mean = []
delta_2_Vp_std = []
delta_2_Vs_mean = []
delta_2_Vs_std = []

RF_DEPTH_mean_2_true_Pds = []
RF_DEPTH_std_2_true_Pds = []
RF_DEPTH_mean_2_true_Ppds = []
RF_DEPTH_std_2_true_Ppds = []

thickness_MTZ_Pds = []
thickness_MTZ_Ppds = []

thickness_MTZ_Pds_std = []
thickness_MTZ_Ppds_std = []

true_thickness_MTZ_Pds = []
true_thickness_MTZ_Ppds = []

true_thickness_MTZ_Pds_std = []
true_thickness_MTZ_Ppds_std = []

diff_thickness_MTZ_Pds = []
diff_thickness_MTZ_Ppds = []

diff_thickness_MTZ_Pds_std = []
diff_thickness_MTZ_Ppds_std  = []

RF_BOOTSTRAP_DATA_Pds = []
RF_BOOTSTRAP_DATA_Ppds = []

RF_BOOTSTRAP_410_DEPTH_Pds = []
RF_BOOTSTRAP_410_DEPTH_Ppds = []
RF_BOOTSTRAP_660_DEPTH_Pds = []
RF_BOOTSTRAP_660_DEPTH_Ppds = []


for i,j in enumerate(RF_data_raw_Pds):
	if len(j) > NUMBER_PP_PER_BIN:
		RF_lat.append(RF_BOOTSTRAP_ESTIMATION_Pds[0][i]['lat'])
		RF_lon.append(RF_BOOTSTRAP_ESTIMATION_Pds[0][i]['lon'])

		flat_DATA_list_Pds = [val for sublist in [RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_DATA'] for _k in range(BOOTSTRAP_INTERATOR)] for val in sublist]
		flat_DATA_list_Ppds = [val for sublist in [RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['RF_DATA'] for _k in range(BOOTSTRAP_INTERATOR)] for val in sublist]
		RF_BOOTSTRAP_DATA_Pds.append(flat_DATA_list_Pds)
		RF_BOOTSTRAP_DATA_Ppds.append(flat_DATA_list_Ppds)

		
		flat_410_list_Pds = [val for sublist in [RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_410_DEPTH'] for _k in range(BOOTSTRAP_INTERATOR)] for val in sublist]
		flat_410_list_Ppds = [val for sublist in [RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['RF_410_DEPTH'] for _k in range(BOOTSTRAP_INTERATOR)] for val in sublist]
		flat_660_list_Pds = [val for sublist in [RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_660_DEPTH'] for _k in range(BOOTSTRAP_INTERATOR)] for val in sublist]
		flat_660_list_Ppds = [val for sublist in [RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['RF_660_DEPTH'] for _k in range(BOOTSTRAP_INTERATOR)] for val in sublist]
		RF_BOOTSTRAP_410_DEPTH_Pds.append(flat_410_list_Pds)
		RF_BOOTSTRAP_410_DEPTH_Ppds.append(flat_410_list_Ppds)
		RF_BOOTSTRAP_660_DEPTH_Pds.append(flat_660_list_Pds)
		RF_BOOTSTRAP_660_DEPTH_Ppds.append(flat_660_list_Ppds)

		flat_mean_1_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		flat_std_1_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_std']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_DEPTH_mean_1_Pds.append(np.mean(flat_mean_1_Pds))
		RF_DEPTH_std_1_Pds.append(np.mean(flat_std_1_Pds))

		flat_mean_1_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		flat_std_1_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_std']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_DEPTH_mean_1_Ppds.append(np.mean(flat_mean_1_Ppds))
		RF_DEPTH_std_1_Ppds.append(np.mean(flat_std_1_Ppds))

		flat_mean_2_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		flat_std_2_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_std']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_DEPTH_mean_2_Pds.append(np.mean(flat_mean_2_Pds))
		RF_DEPTH_std_2_Pds.append(np.mean(flat_std_2_Pds))

		flat_mean_2_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		flat_std_2_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_std']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_DEPTH_mean_2_Ppds.append(np.mean(flat_mean_2_Ppds))
		RF_DEPTH_std_2_Ppds.append(np.mean(flat_std_2_Ppds))


		flat_delta_1_Vp_mean = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		delta_1_Vp_mean.append(np.mean(flat_delta_1_Vp_mean))
		delta_1_Vp_std.append(np.std(flat_delta_1_Vp_mean))

		flat_delta_1_Vs_mean = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		delta_1_Vs_mean.append(np.mean(flat_delta_1_Vs_mean))
		delta_1_Vs_std.append(np.std(flat_delta_1_Vs_mean))


		flat_mean_1_true_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_410_mean']) for _k in range(BOOTSTRAP_INTERATOR)] 
		flat_mean_1_true_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]

		RF_DEPTH_mean_1_true_Pds.append(np.mean(flat_mean_1_true_Pds))
		RF_DEPTH_std_1_true_Pds.append(np.std(flat_mean_1_true_Pds))

		RF_DEPTH_mean_1_true_Ppds.append(np.mean(flat_mean_1_true_Ppds))
		RF_DEPTH_std_1_true_Ppds.append(np.std(flat_mean_1_true_Ppds))

		flat_delta_2_Vp = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		delta_2_Vp_mean.append(np.mean(flat_delta_2_Vp))
		delta_2_Vp_std.append(np.std(flat_delta_2_Vp))

		flat_delta_2_Vs = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		delta_2_Vs_mean.append(np.mean(flat_delta_2_Vs))
		delta_2_Vs_std.append(np.std(flat_delta_2_Vs))


		flat_mean_2_true_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		flat_mean_2_true_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]

		RF_DEPTH_mean_2_true_Pds.append(np.mean(flat_mean_2_true_Pds))
		RF_DEPTH_std_2_true_Pds.append(np.std(flat_mean_2_true_Pds))
		RF_DEPTH_mean_2_true_Ppds.append(np.mean(flat_mean_2_true_Ppds))
		RF_DEPTH_std_2_true_Ppds.append(np.std(flat_mean_2_true_Ppds))

		flat_thickness_MTZ_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		flat_thickness_MTZ_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['thickness_MTZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]

		thickness_MTZ_Pds.append(np.mean(flat_thickness_MTZ_Pds))
		thickness_MTZ_Pds_std.append(np.std(flat_thickness_MTZ_Pds))
		thickness_MTZ_Ppds.append(np.mean(flat_thickness_MTZ_Ppds))
		thickness_MTZ_Ppds_std.append(np.std(flat_thickness_MTZ_Ppds))

		flat_true_thickness_MTZ_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_thickness_MTZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		flat_true_thickness_MTZ_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_thickness_MTZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]

		true_thickness_MTZ_Pds.append(np.mean(flat_true_thickness_MTZ_Pds))
		true_thickness_MTZ_Pds_std.append(np.std(flat_true_thickness_MTZ_Pds))
		true_thickness_MTZ_Ppds.append(np.mean(flat_true_thickness_MTZ_Ppds))
		true_thickness_MTZ_Ppds_std.append(np.std(flat_true_thickness_MTZ_Ppds))


		flat_diff_thickness_MTZ_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ']) for _k in range(BOOTSTRAP_INTERATOR)]
		flat_diff_thickness_MTZ_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ']) for _k in range(BOOTSTRAP_INTERATOR)]

		diff_thickness_MTZ_Pds.append(np.mean(flat_diff_thickness_MTZ_Pds))
		diff_thickness_MTZ_Pds_std.append(np.std(flat_diff_thickness_MTZ_Pds))
		diff_thickness_MTZ_Ppds.append(np.mean(flat_diff_thickness_MTZ_Ppds))
		diff_thickness_MTZ_Ppds_std.append(np.std(flat_diff_thickness_MTZ_Ppds))


#############################################################################################################################################################################################

print('Stacking Pds and Ppds data')
len_RF_stacking_Pds = []
RF_stacking_Pds = []

for i,j in enumerate(RF_data_raw_Pds):
	if len(j) > NUMBER_PP_PER_BIN:
		RF_stacking_Pds.append([sum(x)/len(j)  for x in zip(*j)])
		len_RF_stacking_Pds.append(len(j))

len_RF_stacking_Ppds = []
RF_stacking_Ppds = []

for i,j in enumerate(RF_data_raw_Ppds):
	if len(j) > NUMBER_PP_PER_BIN:
		RF_stacking_Ppds.append([sum(x)/len(j)  for x in zip(*j)])
		len_RF_stacking_Ppds.append(len(j))


#############################################################################################################################################################################################

print('Plotting Figure: Depth of the Mantle Transition Zone for Pds and Ppds phases for 410 and 660 km ...')
#Figure Depth of the Mantle Transition Zone for Pds and Ppds phases for 410 and 660 km

fig, axes = plt.subplots(nrows=2, ncols=2,figsize=(10,10),squeeze=False,sharex=False,sharey=False)

ax1 = axes[0, 0]
ax = axes[0, 1]
ax2 = axes[1, 0]
ax3 = axes[1, 1]

colormap = 'seismic'

#Figure Depth of the Mantle Transition Zone for Pds phase for 410 km

m1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax1)

m1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m1(RF_lon,RF_lat)
sc1 = m1.scatter(x,y,40,RF_DEPTH_mean_1_Pds,cmap=colormap,marker='s',edgecolors='none',vmin=410-DEPTH_RANGE,vmax=410+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m1(lon, lat)
    msize = 10
    l1, = m1.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m1.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m1.drawcoastlines(color='k',zorder=1)
m1.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m1.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])


ax1.set_title('410 km Pds', y=1.08)


#Figure Depth of the Mantle Transition Zone for Ppds phase for 410 km

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m(RF_lon,RF_lat)
sc = m.scatter(x,y,40,RF_DEPTH_mean_1_Ppds,cmap=colormap,marker='s',edgecolors='none',vmin=410-DEPTH_RANGE,vmax=410+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m(lon, lat)
    msize = 10
    l1, = m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m.drawcoastlines(color='k',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax.set_title('410 km Ppds', y=1.08)


#Figure Depth of the Mantle Transition Zone for Pds phase for 660 km


m2 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax2)

m2.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m2.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m2(RF_lon,RF_lat)
sc2 = m2.scatter(x,y,40,RF_DEPTH_mean_2_Pds,cmap=colormap,marker='s',edgecolors='none',vmin=660-DEPTH_RANGE,vmax=660+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m2(lon, lat)
    msize = 10
    l1, = m2.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m2.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m2.drawcoastlines(color='k',zorder=1)
m2.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m2.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax2.set_title('660 km Pds', y=1.08)

#Figure Depth of the Mantle Transition Zone for Ppds phase for 660 km

m3 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax3)

m3.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m3.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m3(RF_lon,RF_lat)
sc3 = m3.scatter(x,y,40,RF_DEPTH_mean_2_Ppds,cmap=colormap,marker='s',edgecolors='none',vmin=660-DEPTH_RANGE,vmax=660+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m3(lon, lat)
    msize = 10
    l1, = m3.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m3.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m3.drawcoastlines(color='k',zorder=1)
m3.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m3.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax3.set_title('660 km Ppds', y=1.08)

fig.colorbar(sc, ax=ax,orientation='horizontal')
fig.colorbar(sc1, ax=ax1,orientation='horizontal')
fig.colorbar(sc2, ax=ax2,orientation='horizontal')
fig.colorbar(sc3, ax=ax3,orientation='horizontal')

fig.subplots_adjust(wspace=0.25, hspace=0.25)

fig.suptitle('Depth per bin')

plt.show()

fig.savefig(PP_FIGURE+'DEPTH_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

##################################################################################################
print('Plotting Figure: True depth of the Mantle Transition Zone for Pds and Ppds phases for 410 and 660 km ...')
#Figure True depth of the Mantle Transition Zone for Pds and Ppds phases for 410 and 660 km

fig, axes = plt.subplots(nrows=2, ncols=2,figsize=(10,10),squeeze=False,sharex=False,sharey=False)

ax1 = axes[0, 0]
ax = axes[0, 1]
ax2 = axes[1, 0]
ax3 = axes[1, 1]

colormap = 'seismic'

#Figure Depth of the Mantle Transition Zone for Pds phase for 410 km

m1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax1)

m1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m1(RF_lon,RF_lat)
sc1 = m1.scatter(x,y,40,RF_DEPTH_mean_1_true_Pds,cmap=colormap,marker='s',edgecolors='none',vmin=410-DEPTH_RANGE,vmax=410+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m1(lon, lat)
    msize = 10
    l1, = m1.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m1.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m1.drawcoastlines(color='k',zorder=1)
m1.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m1.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])


ax1.set_title('410 km Pds', y=1.08)

#Figure Depth of the Mantle Transition Zone for Ppds phase for 410 km

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m(RF_lon,RF_lat)
sc = m.scatter(x,y,40,RF_DEPTH_mean_1_true_Ppds,cmap=colormap,marker='s',edgecolors='none',vmin=410-DEPTH_RANGE,vmax=410+DEPTH_RANGE)


for lon, lat in zip(sta_long,sta_lat):
    x,y = m(lon, lat)
    msize = 10
    l1, = m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m.drawcoastlines(color='k',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax.set_title('410 km Ppds', y=1.08)


#Figure Depth of the Mantle Transition Zone for Pds phase for 660 km


m2 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax2)

m2.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m2.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m2(RF_lon,RF_lat)
sc2 = m2.scatter(x,y,40,RF_DEPTH_mean_2_true_Pds,cmap=colormap,marker='s',edgecolors='none',vmin=660-DEPTH_RANGE,vmax=660+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m2(lon, lat)
    msize = 10
    l1, = m2.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m2.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m2.drawcoastlines(color='k',zorder=1)
m2.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m2.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax2.set_title('660 km Pds', y=1.08)

#Figure Depth of the Mantle Transition Zone for Ppds phase for 660 km

m3 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax3)

m3.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m3.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m3(RF_lon,RF_lat)
sc3 = m3.scatter(x,y,40,RF_DEPTH_mean_2_true_Ppds,cmap=colormap,marker='s',edgecolors='none',vmin=660-DEPTH_RANGE,vmax=660+DEPTH_RANGE)


for lon, lat in zip(sta_long,sta_lat):
    x,y = m3(lon, lat)
    msize = 10
    l1, = m3.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m3.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m3.drawcoastlines(color='k',zorder=1)
m3.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m3.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax3.set_title('660 km Ppds', y=1.08)

#cbar_ax1 = fig.add_axes([0.2, 0.5, 0.6, 0.02])
fig.colorbar(sc, ax=ax,orientation='horizontal')
fig.colorbar(sc1, ax=ax1,orientation='horizontal')
fig.colorbar(sc2, ax=ax2,orientation='horizontal')

#cbar_ax2 = fig.add_axes([0.2, 0.1, 0.6, 0.02])
fig.colorbar(sc3, ax=ax3,orientation='horizontal')

fig.subplots_adjust(wspace=0.25, hspace=0.25)

fig.suptitle('True Depth per bin')

plt.show()

fig.savefig(PP_FIGURE+'TRUE_DEPTH_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting Figure: Std (bootstraping) of the Mantle Transition Zone for Pds and Ppds phases for 410 and 660 km ...')

#Figure True std (bootstraping) of the Mantle Transition Zone for Pds and Ppds phases for 410 and 660 km

fig_std, axes_std = plt.subplots(nrows=2, ncols=2,figsize=(10,10),squeeze=False,sharex=False,sharey=False)

ax1_std = axes_std[0, 0]
ax_std = axes_std[0, 1]
ax2_std = axes_std[1, 0]
ax3_std = axes_std[1, 1]

colormap_std = 'viridis'

#Figure Depth of the Mantle Transition Zone for Pds phase for 410 km

m1_std = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax1_std)

m1_std.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m1_std.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x1_std, y1_std = m1(RF_lon,RF_lat)
sc1_std = m1_std.scatter(x1_std,y1_std,40,RF_DEPTH_std_1_Pds,cmap=colormap_std,marker='s',edgecolors='none')


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

#Figure Depth of the Mantle Transition Zone for Ppds phase for 410 km

m_std = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_std)

m_std.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_std.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x_std, y_std = m_std(RF_lon,RF_lat)
sc_std = m_std.scatter(x_std,y_std,40,RF_DEPTH_std_1_Ppds,cmap=colormap_std,marker='s',edgecolors='none')

for lon_std, lat_std in zip(sta_long,sta_lat):
    x,y = m_std(lon_std, lat_std)
    msize = 10
    l1, = m_std.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_std.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_std.drawcoastlines(color='k',zorder=1)
m_std.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_std.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax_std.set_title('410 km Ppds', y=1.08)
fig_std.colorbar(sc_std, ax=ax_std,orientation='horizontal')

#Figure Depth of the Mantle Transition Zone for Pds phase for 660 km


m2_std = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax2_std)

m2_std.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m2_std.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x2_std, y2_std = m2_std(RF_lon,RF_lat)
sc2_std = m2_std.scatter(x2_std,y2_std,40,RF_DEPTH_std_2_Pds,cmap=colormap_std,marker='s',edgecolors='none')

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

#Figure Depth of the Mantle Transition Zone for Ppds phase for 660 km

m3_std = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax3_std)

m3_std.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m3_std.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x3_std, y3_std = m3(RF_lon,RF_lat)
sc3_std = m3_std.scatter(x3_std,y3_std,40,RF_DEPTH_std_2_Ppds,cmap=colormap_std,marker='s',edgecolors='none')


for lon3_std, lat3_std in zip(sta_long,sta_lat):
    x,y = m3_std(lon3_std, lat3_std)
    msize = 10
    l1, = m3_std.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m3_std.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m3_std.drawcoastlines(color='k',zorder=1)
m3_std.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m3_std.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

ax3_std.set_title('660 km Ppds', y=1.08)
fig_std.colorbar(sc3_std, ax=ax3_std,orientation='horizontal')


fig_std.subplots_adjust(wspace=0.25, hspace=0.25)

fig_std.suptitle('Standard Deviation (bootstraping) per bin')

plt.show()

fig_std.savefig(PP_FIGURE+'TRUE_STD_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting Figure: Delta Vp of each bin...')
#Figure Delta Vp of each bin

colormap_delta_vp = 'viridis'

fig_delta_vp, (ax1_delta_vp, ax2_delta_vp) =  plt.subplots(nrows=1, ncols=2,figsize=(10,5))

m1_delta_vp = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax1_delta_vp)

m1_delta_vp.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m1_delta_vp.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m1_delta_vp(RF_lon,RF_lat)
sc1_delta_vp = m1_delta_vp.scatter(x,y,40,delta_1_Vp_mean,cmap=colormap_delta_vp,marker='s',edgecolors='none')

for lon, lat in zip(sta_long,sta_lat):
    x,y = m1_delta_vp(lon, lat)
    msize = 10
    l1, = m1_delta_vp.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m1_delta_vp.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m1_delta_vp.drawcoastlines(color='k',zorder=1)
m1_delta_vp.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m1_delta_vp.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_delta_vp.colorbar(sc1_delta_vp,ax=ax1_delta_vp,orientation='horizontal')
ax1_delta_vp.set_title('Delta Vp - 410 km',y=1.08)

m_delta_vp = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax2_delta_vp)

m_delta_vp.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_delta_vp.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_delta_vp(RF_lon,RF_lat)
sc_delta_vp = m_delta_vp.scatter(x,y,40,delta_2_Vp_mean,cmap=colormap_delta_vp,marker='s',edgecolors='none')

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_delta_vp(lon, lat)
    msize = 10
    l1, = m_delta_vp.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_delta_vp.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_delta_vp.drawcoastlines(color='k',zorder=1)
m_delta_vp.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_delta_vp.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_delta_vp.colorbar(sc_delta_vp,ax=ax2_delta_vp,orientation='horizontal')
ax2_delta_vp.set_title('Delta Vp - 660 km',y=1.08)

fig_delta_vp.suptitle(r'$\delta V_{p}$ for each bin'+'\n'+
r'$\delta V_{p} = \frac{\alpha . \beta . (H_{A}^{(Ppds)} - H_{A}^{(Pds)})}{\alpha . (1 + \frac{\gamma.V_{s0}}{V_{p0}}) . H_{A}^{(Pds)} - \beta . (1 - \frac{\gamma.V_{s0}}{V_{p0}}). H_{A}^{(Ppds)}}$',
ha='center',va='top')

plt.show()

fig_delta_vp.savefig(PP_FIGURE+'DELTA_VP_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting Figure: Thickness of the Mantle Transition Zone...')
#Figure Thickness of the Mantle Transition Zone

colormap_MTZ = 'seismic_r'

fig_thickness, (ax_thickness1, ax_thickness2) = plt.subplots(nrows=1, ncols=2,figsize=(10,5))

m_thickness1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness1)

m_thickness1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness1(RF_lon,RF_lat)
sc_thickness1 = m_thickness1.scatter(x,y,40,thickness_MTZ_Pds,cmap=colormap_MTZ,marker='s',edgecolors='none',vmin=250-DEPTH_RANGE,vmax=250+DEPTH_RANGE)

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

x, y = m_thickness2(RF_lon,RF_lat)
sc_thickness2 = m_thickness2.scatter(x,y,40,thickness_MTZ_Ppds,cmap=colormap_MTZ,marker='s',edgecolors='none',vmin=250-DEPTH_RANGE,vmax=250+DEPTH_RANGE)

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
fig_thickness.savefig(PP_FIGURE+'THICKNESS_MTZ_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting Figure: True Thickness of the Mantle Transition Zone...')
#Figure Thickness of the Mantle Transition Zone

colormap_MTZ = 'seismic_r'

fig_thickness, (ax_thickness1, ax_thickness2) = plt.subplots(nrows=1, ncols=2,figsize=(10,5))

m_thickness1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness1)

m_thickness1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness1(RF_lon,RF_lat)
sc_thickness1 = m_thickness1.scatter(x,y,40,true_thickness_MTZ_Pds,cmap=colormap_MTZ,marker='s',edgecolors='none',vmin=250-DEPTH_RANGE,vmax=250+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness1(lon, lat)
    msize = 10
    l1, = m_thickness1.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness1.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness1.drawcoastlines(color='k',zorder=1)
m_thickness1.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness1.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness1,ax=ax_thickness1,orientation='horizontal')
ax_thickness1.set_title('True Thickness of MTZ (Pds)',y=1.08)


m_thickness2 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness2)

m_thickness2.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness2.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness2(RF_lon,RF_lat)
sc_thickness2 = m_thickness2.scatter(x,y,40,true_thickness_MTZ_Ppds,cmap=colormap_MTZ,marker='s',edgecolors='none',vmin=250-DEPTH_RANGE,vmax=250+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness2(lon, lat)
    msize = 10
    l1, = m_thickness2.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness2.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness2.drawcoastlines(color='k',zorder=1)
m_thickness2.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness2.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness2,ax=ax_thickness2,orientation='horizontal')
ax_thickness2.set_title('True Thickness of MTZ (Ppds)',y=1.08)

plt.show()
fig_thickness.savefig(PP_FIGURE+'TRUE_THICKNESS_MTZ_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting Figure: Difference between True Thickness and Apparent Thickness of the Mantle Transition Zone...')

colormap_MTZ = 'viridis'

fig_thickness, (ax_thickness1, ax_thickness2) = plt.subplots(nrows=1, ncols=2,figsize=(10,5))

m_thickness1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness1)

m_thickness1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness1(RF_lon,RF_lat)
sc_thickness1 = m_thickness1.scatter(x,y,40,diff_thickness_MTZ_Pds,cmap=colormap_MTZ,marker='s',edgecolors='none')

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness1(lon, lat)
    msize = 10
    l1, = m_thickness1.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness1.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness1.drawcoastlines(color='k',zorder=1)
m_thickness1.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness1.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness1,ax=ax_thickness1,orientation='horizontal')
ax_thickness1.set_title('Difference Thickness of MTZ (Pds)',y=1.08)


m_thickness2 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness2)

m_thickness2.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness2.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness2(RF_lon,RF_lat)
sc_thickness2 = m_thickness2.scatter(x,y,40,diff_thickness_MTZ_Ppds,cmap=colormap_MTZ,marker='s',edgecolors='none')

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness2(lon, lat)
    msize = 10
    l1, = m_thickness2.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness2.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness2.drawcoastlines(color='k',zorder=1)
m_thickness2.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness2.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness2,ax=ax_thickness2,orientation='horizontal')
ax_thickness2.set_title('Difference Thickness of MTZ (Ppds)',y=1.08)

plt.show()
fig_thickness.savefig(PP_FIGURE+'DIFFERENCE_BETWEEN_TRUE_THICKNESS_APPARENT_THICKNESS_MTZ_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################


print('Plotting Figure: Thickness, True Thickness and Difference Thickness of the Mantle Transition Zone...')

colormap_MTZ = 'seismic_r'

fig_thickness, ax_thickness = plt.subplots(nrows=3, ncols=2,figsize=(10,40))

ax_thickness1 = ax_thickness[0,0]
ax_thickness2 = ax_thickness[0,1]
ax_thickness3 = ax_thickness[1,0]
ax_thickness4 = ax_thickness[1,1]
ax_thickness5 = ax_thickness[2,0]
ax_thickness6 = ax_thickness[2,1]

m_thickness1 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness1)

m_thickness1.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness1.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness1(RF_lon,RF_lat)
sc_thickness1 = m_thickness1.scatter(x,y,40,thickness_MTZ_Pds,cmap=colormap_MTZ,marker='s',edgecolors='none',vmin=250-DEPTH_RANGE,vmax=250+DEPTH_RANGE)

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

########

m_thickness2 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness2)

m_thickness2.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness2.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness2(RF_lon,RF_lat)
sc_thickness2 = m_thickness2.scatter(x,y,40,thickness_MTZ_Ppds,cmap=colormap_MTZ,marker='s',edgecolors='none',vmin=250-DEPTH_RANGE,vmax=250+DEPTH_RANGE)

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

########

m_thickness3 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness3)

m_thickness3.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness3.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness3(RF_lon,RF_lat)
sc_thickness3 = m_thickness3.scatter(x,y,40,true_thickness_MTZ_Pds,cmap=colormap_MTZ,marker='s',edgecolors='none',vmin=250-DEPTH_RANGE,vmax=250+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness3(lon, lat)
    msize = 10
    l1, = m_thickness3.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness3.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness3.drawcoastlines(color='k',zorder=1)
m_thickness3.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness3.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness3,ax=ax_thickness3,orientation='horizontal')
ax_thickness3.set_title('True Thickness of MTZ (Pds)',y=1.08)

########

m_thickness4 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness4)

m_thickness4.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness4.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness4(RF_lon,RF_lat)
sc_thickness4 = m_thickness4.scatter(x,y,40,true_thickness_MTZ_Ppds,cmap=colormap_MTZ,marker='s',edgecolors='none',vmin=250-DEPTH_RANGE,vmax=250+DEPTH_RANGE)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness4(lon, lat)
    msize = 10
    l1, = m_thickness4.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness4.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness4.drawcoastlines(color='k',zorder=1)
m_thickness4.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness4.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness4,ax=ax_thickness4,orientation='horizontal')
ax_thickness4.set_title('True Thickness of MTZ (Ppds)',y=1.08)

########

m_thickness5 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness5)

m_thickness5.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness5.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness5(RF_lon,RF_lat)
sc_thickness5 = m_thickness5.scatter(x,y,40,diff_thickness_MTZ_Pds,cmap=colormap_MTZ,marker='s',edgecolors='none')

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness5(lon, lat)
    msize = 10
    l1, = m_thickness5.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness5.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness5.drawcoastlines(color='k',zorder=1)
m_thickness5.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness5.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness5,ax=ax_thickness5,orientation='horizontal')
ax_thickness5.set_title('Difference Thickness of MTZ (Pds)',y=1.08)

########

m_thickness6 = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax_thickness6)

m_thickness6.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m_thickness6.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

x, y = m_thickness6(RF_lon,RF_lat)
sc_thickness6 = m_thickness6.scatter(x,y,40,diff_thickness_MTZ_Ppds,cmap=colormap_MTZ,marker='s',edgecolors='none')

for lon, lat in zip(sta_long,sta_lat):
    x,y = m_thickness6(lon, lat)
    msize = 10
    l1, = m_thickness6.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m_thickness6.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m_thickness6.drawcoastlines(color='k',zorder=1)
m_thickness6.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m_thickness6.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
fig_thickness.colorbar(sc_thickness6,ax=ax_thickness6,orientation='horizontal')
ax_thickness6.set_title('Difference Thickness of MTZ (Ppds)',y=1.08)


plt.show()
fig_thickness.subplots_adjust(wspace=0.25, hspace=0.4)
fig_thickness.savefig(PP_FIGURE+'MOSAIC_DIFFERENCE_BETWEEN_TRUE_THICKNESS_APPARENT_THICKNESS_MTZ_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)


print('Saving Selected Piercing Points in JSON file')
print('\n')

os.makedirs(PP_SELEC_DIR,exist_ok=True)

SELECTED_BINNED_DATA_dic = {'lat':[],'lon':[],'len_Pds':[],'len_Ppds':[],'true_mean_1_Pds':[],'true_std_1_Pds':[],'true_mean_2_Pds':[],'true_std_2_Pds':[],
'true_mean_1_Ppds':[],'true_std_1_Ppds':[],'true_mean_2_Ppds':[],'true_std_2_Ppds':[],'mean_1_Pds':[],'std_1_Pds':[],'mean_2_Pds':[],'std_2_Pds':[],'mean_1_Ppds':[],'std_1_Ppds':[],'mean_2_Ppds':[],
'std_2_Ppds':[],'delta_1_Vp_mean':[],'delta_1_Vp_std':[],'delta_2_Vp_mean':[],'delta_2_Vp_std':[],'delta_1_Vs_mean':[],'delta_1_Vs_std':[],'delta_2_Vs_mean':[],'delta_2_Vs_std':[],
'mtz_thickness_Pds':[],'mtz_thickness_Pds_std':[],'true_thickness_MTZ_Pds':[],'true_thickness_MTZ_Pds_std':[],'true_thickness_MTZ_Ppds':[],'true_thickness_MTZ_Ppds_std':[],
'mtz_thickness_Ppds':[],'mtz_thickness_Ppds_std':[],'difference_thickness_MTZ_Pds':[],'difference_thickness_MTZ_Pds_std':[],'difference_thickness_MTZ_Ppds':[],'difference_thickness_MTZ_Ppds_std':[],
'data_Pds':[],'data_Ppds':[],'data_BOOTSTRAP_Pds':[],'data_BOOTSTRAP_Ppds':[],'RF_BOOTSTRAP_410_DEPTH_Pds':[],'RF_BOOTSTRAP_410_DEPTH_Ppds':[],'RF_BOOTSTRAP_660_DEPTH_Pds':[],'RF_BOOTSTRAP_660_DEPTH_Ppds':[],}
for i,j in enumerate(RF_BOOTSTRAP_DATA_Pds):
	SELECTED_BINNED_DATA_dic['data_BOOTSTRAP_Pds'].append(j)
	SELECTED_BINNED_DATA_dic['data_BOOTSTRAP_Ppds'].append(RF_BOOTSTRAP_DATA_Ppds[i])

	SELECTED_BINNED_DATA_dic['data_Pds'].append(RF_stacking_Pds[i])
	SELECTED_BINNED_DATA_dic['data_Ppds'].append(RF_stacking_Ppds[i])

	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_410_DEPTH_Pds'].append(RF_BOOTSTRAP_410_DEPTH_Pds[i])
	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_410_DEPTH_Ppds'].append(RF_BOOTSTRAP_410_DEPTH_Ppds[i])
	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_660_DEPTH_Pds'].append(RF_BOOTSTRAP_660_DEPTH_Pds[i])
	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_660_DEPTH_Ppds'].append(RF_BOOTSTRAP_660_DEPTH_Ppds[i])

	SELECTED_BINNED_DATA_dic['lat'].append(RF_lat[i])
	SELECTED_BINNED_DATA_dic['lon'].append(RF_lon[i])

	SELECTED_BINNED_DATA_dic['len_Pds'].append(len_RF_stacking_Pds[i])
	SELECTED_BINNED_DATA_dic['len_Ppds'].append(len_RF_stacking_Ppds[i])

	SELECTED_BINNED_DATA_dic['true_mean_1_Pds'].append(float(RF_DEPTH_mean_1_true_Pds[i]))
	SELECTED_BINNED_DATA_dic['true_std_1_Pds'].append(float(RF_DEPTH_std_1_true_Pds[i]))
	SELECTED_BINNED_DATA_dic['true_mean_2_Pds'].append(float(RF_DEPTH_mean_2_true_Pds[i]))
	SELECTED_BINNED_DATA_dic['true_std_2_Pds'].append(float(RF_DEPTH_std_2_true_Pds[i]))

	SELECTED_BINNED_DATA_dic['mean_1_Pds'].append(float(RF_DEPTH_mean_1_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_1_Pds'].append(float(RF_DEPTH_std_1_Pds[i]))
	SELECTED_BINNED_DATA_dic['mean_2_Pds'].append(float(RF_DEPTH_mean_2_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_2_Pds'].append(float(RF_DEPTH_std_2_Pds[i]))

	SELECTED_BINNED_DATA_dic['true_mean_1_Ppds'].append(float(RF_DEPTH_mean_1_true_Ppds[i]))
	SELECTED_BINNED_DATA_dic['true_std_1_Ppds'].append(float(RF_DEPTH_std_1_true_Ppds[i]))
	SELECTED_BINNED_DATA_dic['true_mean_2_Ppds'].append(float(RF_DEPTH_mean_2_true_Ppds[i]))
	SELECTED_BINNED_DATA_dic['true_std_2_Ppds'].append(float(RF_DEPTH_std_2_true_Ppds[i]))

	SELECTED_BINNED_DATA_dic['mean_1_Ppds'].append(float(RF_DEPTH_mean_1_Ppds[i]))
	SELECTED_BINNED_DATA_dic['std_1_Ppds'].append(float(RF_DEPTH_std_1_Ppds[i]))
	SELECTED_BINNED_DATA_dic['mean_2_Ppds'].append(float(RF_DEPTH_mean_2_Ppds[i]))
	SELECTED_BINNED_DATA_dic['std_2_Ppds'].append(float(RF_DEPTH_std_2_Ppds[i]))

	SELECTED_BINNED_DATA_dic['delta_1_Vp_mean'].append(float(delta_1_Vp_mean[i]))
	SELECTED_BINNED_DATA_dic['delta_1_Vp_std'].append(float(delta_1_Vp_std[i]))
	SELECTED_BINNED_DATA_dic['delta_1_Vs_mean'].append(float(delta_1_Vs_mean[i]))
	SELECTED_BINNED_DATA_dic['delta_1_Vs_std'].append(float(delta_1_Vs_std[i]))

	SELECTED_BINNED_DATA_dic['delta_2_Vp_mean'].append(float(delta_2_Vp_mean[i]))
	SELECTED_BINNED_DATA_dic['delta_2_Vp_std'].append(float(delta_2_Vp_std[i]))
	SELECTED_BINNED_DATA_dic['delta_2_Vs_mean'].append(float(delta_2_Vs_mean[i]))
	SELECTED_BINNED_DATA_dic['delta_2_Vs_std'].append(float(delta_2_Vs_std[i]))

	SELECTED_BINNED_DATA_dic['mtz_thickness_Pds'].append(float(thickness_MTZ_Pds[i]))
	SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std'].append(float(thickness_MTZ_Pds_std[i]))
	
	SELECTED_BINNED_DATA_dic['mtz_thickness_Ppds'].append(float(thickness_MTZ_Ppds[i]))
	SELECTED_BINNED_DATA_dic['mtz_thickness_Ppds_std'].append(float(thickness_MTZ_Ppds_std[i]))

	SELECTED_BINNED_DATA_dic['true_thickness_MTZ_Pds'].append(float(true_thickness_MTZ_Pds[i]))
	SELECTED_BINNED_DATA_dic['true_thickness_MTZ_Pds_std'].append(float(true_thickness_MTZ_Pds_std[i]))

	SELECTED_BINNED_DATA_dic['true_thickness_MTZ_Ppds'].append(float(true_thickness_MTZ_Ppds[i]))
	SELECTED_BINNED_DATA_dic['true_thickness_MTZ_Ppds_std'].append(float(true_thickness_MTZ_Ppds_std[i]))

	SELECTED_BINNED_DATA_dic['difference_thickness_MTZ_Pds'].append(float(diff_thickness_MTZ_Pds[i]))
	SELECTED_BINNED_DATA_dic['difference_thickness_MTZ_Pds_std'].append(float(diff_thickness_MTZ_Pds_std[i]))

	SELECTED_BINNED_DATA_dic['difference_thickness_MTZ_Ppds'].append(float(diff_thickness_MTZ_Ppds[i]))
	SELECTED_BINNED_DATA_dic['difference_thickness_MTZ_Ppds_std'].append(float(diff_thickness_MTZ_Ppds_std[i]))


with open(PP_SELEC_DIR+'SELECTED_BINNED_Ps.json', 'w') as fp:
	json.dump(SELECTED_BINNED_DATA_dic, fp)
