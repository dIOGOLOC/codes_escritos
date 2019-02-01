# coding: utf-8


import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import copy
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature
import shapefile
from fatiando import gridder, utils
from scipy.stats import mode
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import json
import random
from matplotlib.colors import Normalize
from numpy import ma
from matplotlib import cbook
import collections
from matplotlib.collections import PatchCollection

from shapely.geometry import Polygon, MultiPoint, Point
from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle,Rectangle
import math




from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					NUMBER_PP_PER_BIN,MIN_AMP_PDS_PPDS,DEPTH_TARGET,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,
					LLCRNRLON_SMALL,URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,
					PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,BOUNDARY_1_SHP,BOUNDARY_2_SHP,
					EXT_FIG,DPI_FIG,FRESNEL_ZONE_RADIUS,DIST_GRID_PP,NUMBER_STA_PER_BIN,OUTPUT_DIR,MIN_AMP_GOOD,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,BOOTSTRAP_DEPTH_ESTIMATION,GAMMA,COLORMAP_STD,COLORMAP_VEL
				   )

print('Starting Receiver Functions migration code to estimate the true depths of the Earth discontinuities')
print('\n')

print('Importing depths and times of Pds conversion  dataset')
print('\n')

print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)
print('\n')


for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 410:
		Vp_depth_1 = j[3]
		Vs_depth_1 = j[5]

for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 520:
		Vp_depth_520 = j[3]
		Vs_depth_520 = j[5]
		
for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 660:
		Vp_depth_2 = j[3]
		Vs_depth_2 = j[5]

print('410 km earth model Vp : '+str(Vp_depth_1))
print('410 km earth model Vs : '+str(Vs_depth_1))
print('520 km earth model Vp : '+str(Vp_depth_520))
print('520 km earth model Vs : '+str(Vs_depth_520))
print('660 km earth model Vp : '+str(Vp_depth_2))
print('660 km earth model Vs : '+str(Vs_depth_2))
print('\n')

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Stations'+'/'

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

print('Importing Pds piercing points to each PHASE')
print('\n')

PHASES = 'P410s','P'+str(DEPTH_TARGET)+'s','P660s'

print('Importing Pds Piercing Points for '+PHASES[0])
print('\n')

PP_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Piercing_Points'+'/'


filename_1 = PP_DIR+'PP_'+PHASES[0]+'_dic.json'

PP_1_dic = json.load(open(filename_1))


PP_time_1 = PP_1_dic['time']
PP_lat_1 = PP_1_dic['lat']
PP_lon_1 = PP_1_dic['lon']
PP_depth_1 = PP_1_dic['depth']


print('Importing Pds Piercing Points for '+PHASES[1])
print('\n')

filename_med = PP_DIR+'PP_'+PHASES[1]+'_dic.json'

PP_med_dic = json.load(open(filename_med))

PP_time_med = PP_med_dic['time']
PP_lat_med = PP_med_dic['lat']
PP_lon_med = PP_med_dic['lon']
PP_depth_med = PP_med_dic['depth']


print('Importing Pds Piercing Points for '+PHASES[2])
print('\n')

filename_2 = PP_DIR+'PP_'+PHASES[2]+'_dic.json'

PP_2_dic = json.load(open(filename_2))

PP_time_2 = PP_2_dic['time']
PP_lat_2 = PP_2_dic['lat']
PP_lon_2 = PP_2_dic['lon']
PP_depth_2 = PP_2_dic['depth'] 

print('P410s Piercing Points')
print('\n')

pp_1_lat  = [[]]*len(PP_lon_1)
pp_1_long  = [[]]*len(PP_lon_1)


for i,j in enumerate(PP_lon_1):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_1[i][k] == 410:
				pp_1_lat[i] = PP_lat_1[i][k]
				pp_1_long[i] = l

pp_1_lat = [i for i in pp_1_lat if type(i) == float ]
pp_1_long = [i for i in pp_1_long if type(i) == float ]

print('P'+str(DEPTH_TARGET)+'s Piercing Points')
print('\n')

pp_med_lat  = [[]]*len(PP_lon_med)
pp_med_long  = [[]]*len(PP_lon_med)


for i,j in enumerate(PP_lon_med):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_med[i][k] == DEPTH_TARGET:
				pp_med_lat[i] = PP_lat_med[i][k]
				pp_med_long[i] = l

pp_med_lat = [i for i in pp_med_lat if type(i) == float ]
pp_med_long = [i for i in pp_med_long if type(i) == float ]

print('P660s Piercing Points')
print('\n')

pp_2_lat  = [[]]*len(PP_lon_2)
pp_2_long  = [[]]*len(PP_lon_2)


for i,j in enumerate(PP_lon_2):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_2[i][k] == 660:
			if PP_lon_2[i][k]  != [] and PP_lat_2[i][k] != []:
				pp_2_lat[i] = PP_lat_2[i][k]
				pp_2_long[i] = l

pp_2_lat = [i for i in pp_2_lat if type(i) == float ]
pp_2_long = [i for i in pp_2_long if type(i) == float ]

print('Importing Ppds piercing points to each PHASE')
print('\n')

PHASES_Ppds = 'PPv410s','PPv'+str(DEPTH_TARGET)+'s','PPv660s'

print('Importing Ppds Piercing Points '+PHASES_Ppds[0])
print('\n')

filename_1_Ppds = PP_DIR+'PP_'+PHASES_Ppds[0]+'_dic.json'

PP_1_dic_Ppds = json.load(open(filename_1_Ppds))

PP_time_1_Ppds = PP_1_dic_Ppds['time']
PP_lat_1_Ppds = PP_1_dic_Ppds['lat']
PP_lon_1_Ppds = PP_1_dic_Ppds['lon']
PP_depth_1_Ppds = PP_1_dic_Ppds['depth'] 

print('Importing Ppds Piercing Points '+PHASES_Ppds[1])
print('\n')

filename_med_Ppds = PP_DIR+'PP_'+PHASES_Ppds[1]+'_dic.json'

PP_med_dic_Ppds = json.load(open(filename_med_Ppds))

PP_time_med_Ppds = PP_med_dic_Ppds['time']
PP_lat_med_Ppds = PP_med_dic_Ppds['lat']
PP_lon_med_Ppds = PP_med_dic_Ppds['lon']
PP_depth_med_Ppds = PP_med_dic_Ppds['depth']


print('Importing Ppds Piercing Points '+PHASES_Ppds[2])
print('\n')

filename_2_Ppds = PP_DIR+'PP_'+PHASES_Ppds[2]+'_dic.json'

PP_2_dic_Ppds = json.load(open(filename_2_Ppds))

PP_time_2_Ppds =  PP_2_dic_Ppds['time']
PP_lat_2_Ppds = PP_2_dic_Ppds['lat']
PP_lon_2_Ppds = PP_2_dic_Ppds['lon']
PP_depth_2_Ppds = PP_2_dic_Ppds['depth']

print('PPv410s Piercing Points')
print('\n')

pp_1_lat_Ppds  = [[]]*len(PP_lon_1_Ppds)
pp_1_long_Ppds  = [[]]*len(PP_lon_1_Ppds)


for i,j in enumerate(PP_lon_1_Ppds):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_1_Ppds[i][k] == 410:
				pp_1_lat_Ppds[i] = PP_lat_1_Ppds[i][k]
				pp_1_long_Ppds[i] = l

pp_1_lat_Ppds = [i for i in pp_1_lat_Ppds if type(i) == float ]
pp_1_long_Ppds = [i for i in pp_1_long_Ppds if type(i) == float ]


print('PPv'+str(DEPTH_TARGET)+'s Piercing Points')
print('\n')

pp_med_lat_Ppds  = [[]]*len(PP_lon_med_Ppds)
pp_med_long_Ppds  = [[]]*len(PP_lon_med_Ppds)


for i,j in enumerate(PP_lon_med_Ppds):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_med_Ppds[i][k] == DEPTH_TARGET:
				pp_med_lat_Ppds[i] = PP_lat_med_Ppds[i][k]
				pp_med_long_Ppds[i] = l

pp_med_lat_Ppds = [i for i in pp_med_lat_Ppds if type(i) == float ]
pp_med_long_Ppds = [i for i in pp_med_long_Ppds if type(i) == float ]

print('PPv660s Piercing Points')
print('\n')


pp_2_lat_Ppds  = [[]]*len(PP_lon_2_Ppds)
pp_2_long_Ppds  = [[]]*len(PP_lon_2_Ppds)
pp_2_depth_Ppds  = [[]]*len(PP_lon_2_Ppds)


for i,j in enumerate(PP_lon_2_Ppds):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_2_Ppds[i][k] == 660:
				pp_2_lat_Ppds[i] = PP_lat_2_Ppds[i][k]
				pp_2_long_Ppds[i] = l

pp_2_lat_Ppds = [i for i in pp_2_lat_Ppds if type(i) == float ]
pp_2_long_Ppds = [i for i in pp_2_long_Ppds if type(i) == float ]


print('Creating GRID POINTS')
print('\n')

area = (LLCRNRLON_SMALL,URCRNRLON_SMALL, LLCRNRLAT_SMALL, URCRNRLAT_SMALL)

shape = (abs(abs(URCRNRLON_SMALL) - abs(LLCRNRLON_SMALL))*GRID_PP_MULT, abs(abs(URCRNRLAT_SMALL) - abs(LLCRNRLAT_SMALL))*GRID_PP_MULT)

grdx, grdy = gridder.regular(area, shape)

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
    if vect_j[indices[NUMBER_PP_PER_BIN]] < FRESNEL_ZONE_RADIUS:
        grid_sel_min.append((grdx[i],grdy[i]))


grid_sel_med = []
grid_sel_med_data = []
for i,j in enumerate(dist_pp_grid_med):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[NUMBER_PP_PER_BIN]] < FRESNEL_ZONE_RADIUS:
        grid_sel_med.append((grdx[i],grdy[i]))

        
grid_sel_max = []
grid_sel_min_data = []

for i,j in enumerate(dist_pp_grid_max):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[NUMBER_PP_PER_BIN]] < FRESNEL_ZONE_RADIUS:
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


fig_PP_Pds_Ppds, (ax, ax1) = plt.subplots(ncols=2, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(20,10),sharey=True)

#Figure Pds

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
l1, = ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.Geodetic())
l2, = ax.plot(pp_1_long,pp_1_lat, '.',markersize=5,markeredgecolor='k',markerfacecolor='b',transform=ccrs.Geodetic())
l3, = ax.plot(pp_med_long,pp_med_lat, '.',markersize=5,markeredgecolor='k',markerfacecolor='g',transform=ccrs.Geodetic())
l4, = ax.plot(pp_2_long,pp_2_lat, '.',markersize=5,markeredgecolor='k',markerfacecolor='r',transform=ccrs.Geodetic())

for i,j in enumerate(grdx):
	circulo = Circle(radius=DIST_GRID_PP,xy=(grdx[i], grdy[i]),color='None', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
	ax.add_patch(circulo)
ax.set_title('Pds Piercing Points',ha='center',va='top',y=1.08)
ax.legend([l1,l2,l3,l4,circulo],['Stations','Piercing Points 410 km','Piercing Points 530 km','Piercing Points 660 km','Selected Grid'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w',fontsize='smaller')

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax.gridlines(draw_labels=True)


##############################################################################################

ax1.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
l1, = ax1.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.Geodetic())
l2, = ax1.plot(pp_1_long_Ppds,pp_1_lat_Ppds, '.',markersize=5,markeredgecolor='k',markerfacecolor='b',transform=ccrs.Geodetic())
l3, = ax1.plot(pp_med_long_Ppds,pp_med_lat_Ppds, '.',markersize=5,markeredgecolor='k',markerfacecolor='g',transform=ccrs.Geodetic())
l4, = ax1.plot(pp_2_long_Ppds,pp_2_lat_Ppds, '.',markersize=5,markeredgecolor='k',markerfacecolor='r',transform=ccrs.Geodetic())
for i,j in enumerate(grdx):
	circulo = Circle(radius=DIST_GRID_PP,xy=(grdx[i], grdy[i]),color='None', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
	ax1.add_patch(circulo)
ax1.set_title('Ppds Piercing Points',ha='center',va='top',y=1.08)
ax1.legend([l1,l2,l3,l4,circulo],['Stations','Piercing Points 410 km','Piercing Points '+str(DEPTH_TARGET)+' km','Piercing Points 660 km','Selected Grid'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w',fontsize='smaller')

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax1.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax1.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax1.gridlines(draw_labels=True)

print('GRID Check!')
#plt.show()

PP_FIGURE = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Figures'+'/'

RESULTS_FOLDER = PP_FIGURE+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)

print('Folder to save Figures files:')
print(RESULTS_FOLDER)
print('\n')

fig_PP_Pds_Ppds.savefig(RESULTS_FOLDER+'PP_Pds_Ppds.'+EXT_FIG,dpi=DPI_FIG)


###################################################################################################################

print('Plotting: Figure Pds Depth Target Piercing Points')
print('\n')


fig_PP, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#Figure Pds

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
l1, = ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.Geodetic())
l3, = ax.plot(pp_med_long_Ppds,pp_med_lat_Ppds, 'X',markersize=5,markeredgecolor='k',markerfacecolor='k',alpha=0.5,transform=ccrs.Geodetic())

for i,j in enumerate(grdx):
	circulo = Circle(radius=DIST_GRID_PP,xy=(grdx[i], grdy[i]),color='None', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
	ax.add_patch(circulo)
for i,j in enumerate(pp_med_long_Ppds):
	circulo_fresnel = Circle(radius=FRESNEL_ZONE_RADIUS, xy=(pp_med_long_Ppds[i],pp_med_lat_Ppds[i]), color='gray',linewidth=0,alpha=0.2,transform=ccrs.Geodetic(),zorder=1)
	ax.add_patch(circulo_fresnel)

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

ax.set_title('Pds Piercing Points',ha='center',va='top',y=1.08)
ax.legend([l1,l3,circulo,circulo_fresnel],['Stations','Piercing Points '+str(DEPTH_TARGET)+' km','Selected Grid','Piercing Points Fresnel Zone'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w',fontsize='smaller')
ax.gridlines(draw_labels=True)

#plt.show()

fig_PP.savefig(RESULTS_FOLDER+'PP_MED_Pds.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################

print('Plotting: Figure Ppds Depth Target Piercing Points')
print('\n')
fig_PP, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#Figure Ppds

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
l1, = ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.Geodetic())
l3, = ax.plot(pp_med_long_Ppds,pp_med_lat_Ppds, 'X',markersize=5,markeredgecolor='k',markerfacecolor='k',alpha=0.5,transform=ccrs.Geodetic())
for i,j in enumerate(grdx):
	circulo = Circle(radius=DIST_GRID_PP,xy=(grdx[i], grdy[i]),color='None', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
	ax.add_patch(circulo)
for i,j in enumerate(pp_med_long_Ppds):
	circulo_fresnel = Circle(radius=FRESNEL_ZONE_RADIUS, xy=(pp_med_long_Ppds[i],pp_med_lat_Ppds[i]), color='gray',linewidth=0,alpha=0.2,transform=ccrs.Geodetic(),zorder=1)
	ax.add_patch(circulo_fresnel)

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax.gridlines(draw_labels=True)

ax.set_title('Ppds Piercing Points',ha='center',va='top',y=1.08)
ax.legend([l1,l3,circulo,circulo_fresnel],['Stations','Piercing Points '+str(DEPTH_TARGET)+' km','Selected Grid','Piercing Points Fresnel Zone'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w',fontsize='smaller')

#plt.show()

fig_PP.savefig(RESULTS_FOLDER+'PP_MED_Ppds.'+EXT_FIG,dpi=DPI_FIG)

##########################################################################################################################################

print('Importing depths and times of Pds conversion  dataset')
print('\n')

PdS_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Phases'+'/'

filename_Pds = PdS_DIR+'Pds_dic.json'

PdS_Dic = json.load(open(filename_Pds))

depth_str_to_float_Pds = []
for i,j in enumerate(PdS_Dic['depth']):
	depth_str_to_float_Pds.append([float(l) for l in j])

Pds_time = PdS_Dic['time']
Pds_st_lat = PdS_Dic['st_lat']
Pds_st_lon = PdS_Dic['st_long']
Pds_ev_lat = PdS_Dic['ev_lat']
Pds_ev_lon = PdS_Dic['ev_long']
Pds_depth = depth_str_to_float_Pds

print('Importing depths and times of Ppds conversion  dataset')
print('\n')

filename_Ppds = PdS_DIR+'PPvs_dic.json'

Ppds_Dic = json.load(open(filename_Ppds))

depth_str_to_float_Ppds = []
for i,j in enumerate(Ppds_Dic['depth']):
	depth_str_to_float_Ppds.append([float(l) for l in j])

Ppds_time = Ppds_Dic['time']
Ppds_st_lat = Ppds_Dic['st_lat']
Ppds_st_lon = Ppds_Dic['st_long']
Ppds_ev_lat = Ppds_Dic['ev_lat']
Ppds_ev_lon = Ppds_Dic['ev_long']
Ppds_depth = depth_str_to_float_Ppds


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

##################################################################################################################

print('Filtering migrated data per grid data')
print('\n')

dados_grid_lat = pp_med_lat_Ppds
dados_grid_lon = pp_med_long_Ppds

RF_data_raw_Pds = [[]]*len(grid_sel_x)
RF_amplitude_depth_raw_Pds = [[]]*len(grid_sel_x)

RF_data_raw_Ppds = [[]]*len(grid_sel_x)
RF_amplitude_depth_raw_Ppds = [[]]*len(grid_sel_x)

RF_STA_number_raw = [[]]*len(grid_sel_x)

for i,j in enumerate(grid_sel_x):
	RF_data_raw_Pds[i] = [RF_amplitude_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) <= FRESNEL_ZONE_RADIUS]
	RF_amplitude_depth_raw_Pds[i] = [RF_amplitude_depth_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) <= FRESNEL_ZONE_RADIUS]

	RF_data_raw_Ppds[i] = [RF_amplitude_Ppds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) <= FRESNEL_ZONE_RADIUS]
	RF_amplitude_depth_raw_Ppds[i] = [RF_amplitude_depth_Ppds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) <= FRESNEL_ZONE_RADIUS]

	RF_STA_number_raw[i] = len(set([(sta_long[k],sta_lat[k]) for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) <= FRESNEL_ZONE_RADIUS]))

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
			if len(j) >= NUMBER_PP_PER_BIN and RF_STA_number_raw[i] >= NUMBER_STA_PER_BIN:

				print('Grid point number: '+str(i))
				print('lat: '+str(grid_sel_y[i]))
				print('lon: '+str(grid_sel_x[i]))

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lon'] = grid_sel_x[i]
				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lat'] = grid_sel_y[i]

				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['lon'] = grid_sel_x[i]
				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['lat'] = grid_sel_y[i]

				old_RF_DATA_raw_Pds_lst = np.arange(0,len(j))
				new_RANDOM_RF_DATA_raw_Pds_lst = np.random.choice(old_RF_DATA_raw_Pds_lst,size=len(old_RF_DATA_raw_Pds_lst),replace=True)
				
				new_RANDOM_RF_DATA_raw_Pds = [j[_t_] for _t_ in new_RANDOM_RF_DATA_raw_Pds_lst]
				new_RANDOM_RF_DATA_raw_Ppds = [RF_data_raw_Ppds[i][_t_] for _t_ in new_RANDOM_RF_DATA_raw_Pds_lst]

				RF_STACKING_BOOTSTRAP_Pds = [sum(i)/len(new_RANDOM_RF_DATA_raw_Pds) for i in zip(*new_RANDOM_RF_DATA_raw_Pds)]
				RF_STACKING_BOOTSTRAP_Ppds = [sum(i)/len(new_RANDOM_RF_DATA_raw_Ppds) for i in zip(*new_RANDOM_RF_DATA_raw_Ppds)]

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_DATA'] = RF_STACKING_BOOTSTRAP_Pds
				RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['RF_DATA'] = RF_STACKING_BOOTSTRAP_Ppds

				######## Estimating 410 km apparent depth ########

				#410 km Pds

				lst_depth_amp_410_Pds = [RF_STACKING_BOOTSTRAP_Pds[x] for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
				lst_depth_pp_410_Pds = [c for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
				lst_410_depth_Pds = lst_depth_pp_410_Pds[lst_depth_amp_410_Pds.index(max(lst_depth_amp_410_Pds))]
				lst_410_amp_Pds = lst_depth_amp_410_Pds.index(max(lst_depth_amp_410_Pds))

				#410 km Ppds

				lst_depth_amp_410_Ppds = [RF_STACKING_BOOTSTRAP_Ppds[x] for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
				lst_depth_pp_410_Ppds = [c for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
				lst_410_depth_Ppds = lst_depth_pp_410_Ppds[lst_depth_amp_410_Ppds.index(max(lst_depth_amp_410_Ppds))]
				lst_410_amp_Ppds = lst_depth_amp_410_Ppds.index(max(lst_depth_amp_410_Ppds))

				######## Estimating 520 km apparent depth ########

				#520 km Pds

				lst_depth_amp_520_Pds = [RF_STACKING_BOOTSTRAP_Pds[x] for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]
				lst_depth_pp_520_Pds = [c for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]
				lst_520_depth_Pds = lst_depth_pp_520_Pds[lst_depth_amp_520_Pds.index(max(lst_depth_amp_520_Pds))]
				lst_520_amp_Pds = lst_depth_amp_520_Pds.index(max(lst_depth_amp_520_Pds))

				#520 km Ppds

				lst_depth_amp_520_Ppds = [RF_STACKING_BOOTSTRAP_Ppds[x] for x,c in enumerate(camadas_terra_10_km) if 520-DEPTH_RANGE <= c <= 520+DEPTH_RANGE]
				lst_depth_pp_520_Ppds = [c for x,c in enumerate(camadas_terra_10_km) if 520-DEPTH_RANGE <= c <= 520+DEPTH_RANGE]
				lst_520_depth_Ppds = lst_depth_pp_520_Ppds[lst_depth_amp_520_Ppds.index(max(lst_depth_amp_520_Ppds))]
				lst_520_amp_Ppds = lst_depth_amp_520_Ppds.index(max(lst_depth_amp_520_Ppds))


				######## Estimating 660 km apparent depth ########

				#660 km Pds

				lst_depth_amp_660_Pds = [RF_STACKING_BOOTSTRAP_Pds[x] for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
				lst_depth_pp_660_Pds = [c for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
				lst_660_depth_Pds = lst_depth_pp_660_Pds[lst_depth_amp_660_Pds.index(max(lst_depth_amp_660_Pds))]
				lst_660_amp_Pds = lst_depth_amp_660_Pds.index(max(lst_depth_amp_660_Pds))
				
				#660 km Ppds

				lst_depth_amp_660_Ppds = [RF_STACKING_BOOTSTRAP_Ppds[x] for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
				lst_depth_pp_660_Ppds = [c for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
				lst_660_depth_Ppds = lst_depth_pp_660_Ppds[lst_depth_amp_660_Ppds.index(max(lst_depth_amp_660_Ppds))]
				lst_660_amp_Ppds = lst_depth_amp_660_Ppds.index(max(lst_depth_amp_660_Ppds))


				##############################

				if lst_410_amp_Pds >= MIN_AMP_PDS_PPDS:

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'] = lst_410_depth_Pds
				
				else: 

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'] = np.nan

				if lst_410_amp_Ppds >= MIN_AMP_PDS_PPDS:

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_mean'] = lst_410_depth_Ppds
				
				else:

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_mean'] = np.nan
				
				print('410 Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean']))
				print('410 Ppds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_mean']))

				if lst_520_amp_Pds >= MIN_AMP_PDS_PPDS:

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['520_mean'] = lst_520_depth_Pds
				
				else: 

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['520_mean'] = np.nan

				if lst_520_amp_Ppds >= MIN_AMP_PDS_PPDS:

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['520_mean'] = lst_520_depth_Ppds
				
				else:

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['520_mean'] = np.nan
				
				print('520 Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['520_mean']))
				print('520 Ppds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['520_mean']))

				if lst_660_depth_Pds >= MIN_AMP_PDS_PPDS:

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'] = lst_660_depth_Pds
			
				else:

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'] = np.nan

				if lst_660_depth_Ppds >= MIN_AMP_PDS_PPDS:

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_mean'] = lst_660_depth_Ppds

				else:

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_mean'] = np.nan

				print('660 km Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean']))
				print('660 km Ppds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_mean']))


				######## Estimating MTZ thickness ########

				if  lst_410_amp_Pds >= MIN_AMP_PDS_PPDS  and lst_660_depth_Pds >= MIN_AMP_PDS_PPDS:

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean'] = lst_660_depth_Pds - lst_410_depth_Pds

				else: 

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean'] = np.nan


				if  lst_410_amp_Ppds >= MIN_AMP_PDS_PPDS and lst_660_depth_Ppds >= MIN_AMP_PDS_PPDS:

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['thickness_MTZ_mean'] = lst_660_depth_Ppds - lst_410_depth_Ppds
				
				else: 

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['thickness_MTZ_mean'] = np.nan

				print('MTZ Pds thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean']))
				print('MTZ Ppds thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['thickness_MTZ_mean']))

				######## Estimating TRUE depth ########

				if  lst_410_amp_Pds >= MIN_AMP_PDS_PPDS  and lst_660_depth_Pds >= MIN_AMP_PDS_PPDS and lst_410_amp_Ppds >= MIN_AMP_PDS_PPDS and lst_660_depth_Ppds >= MIN_AMP_PDS_PPDS:

					Vp_Vs_ratio_depth_1 = Vs_depth_1/Vp_depth_1
					alfa = Vp_depth_1 - Vs_depth_1
					beta = Vp_depth_1 + Vs_depth_1
					gamma_vp_vs_1 = GAMMA*Vp_Vs_ratio_depth_1

					H_a_Pds = lst_410_depth_Pds
					H_a_Ppds = lst_410_depth_Ppds

					lst_delta_Vp_410_Pds = (alfa*beta*(H_a_Ppds - H_a_Pds)) * ((alfa*(1+gamma_vp_vs_1)*H_a_Pds) - (beta*(1-gamma_vp_vs_1)*H_a_Ppds))**(-1)

					lst_delta_Vs_410_Pds = lst_delta_Vp_410_Pds * GAMMA * Vp_Vs_ratio_depth_1
				
					lst_true_410_mean_Pds = (H_a_Pds) * ((Vp_depth_1-Vs_depth_1)/(Vp_depth_1*Vs_depth_1)) * (((Vs_depth_1+lst_delta_Vs_410_Pds)*(Vp_depth_1+lst_delta_Vp_410_Pds))/(Vp_depth_1+lst_delta_Vp_410_Pds-Vs_depth_1-lst_delta_Vs_410_Pds))

					lst_true_410_mean_Ppds = (H_a_Ppds) * ((Vp_depth_1+Vs_depth_1)/(Vp_depth_1*Vs_depth_1)) * (((Vs_depth_1+lst_delta_Vs_410_Pds)*(Vp_depth_1+lst_delta_Vp_410_Pds))/(Vp_depth_1+lst_delta_Vp_410_Pds+Vs_depth_1+lst_delta_Vs_410_Pds))
							
						
					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_410_mean'] = lst_delta_Vp_410_Pds

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_410_mean'] = lst_delta_Vs_410_Pds

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_410_mean'] = lst_true_410_mean_Pds

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_410_mean'] = lst_true_410_mean_Ppds

					print('Equation 10 (Gao & Liu,2013) is true? = '+str(round(H_a_Pds/H_a_Ppds,2) == round(((alfa+lst_delta_Vp_410_Pds-lst_delta_Vs_410_Pds)/(beta+lst_delta_Vp_410_Pds+lst_delta_Vs_410_Pds))*(beta/alfa),2)))

					print('Delta Vp 410 km = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_410_mean']))
					print('Delta Vs 410 km = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_410_mean']))				
						
					print('410 km Pds True Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_410_mean']))
					print('410 km Ppds True Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_410_mean']))
		
					######## Estimating TRUE depth ########

					Vp_Vs_ratio_depth_2 = Vs_depth_2/Vp_depth_2
					alfa = Vp_depth_2 - Vs_depth_2
					beta = Vp_depth_2 + Vs_depth_2
					gamma_vp_vs_2 = GAMMA*Vp_Vs_ratio_depth_2

					H_a_Pds_660 = lst_660_depth_Pds
					H_a_Ppds_660 = lst_660_depth_Ppds

					lst_delta_Vp_660_Pds = (alfa*beta*(H_a_Ppds_660 - H_a_Pds_660)) / (alfa*(1+gamma_vp_vs_2)*H_a_Pds_660 - beta*(1-gamma_vp_vs_2)*H_a_Ppds_660)

					lst_delta_Vs_660_Pds = lst_delta_Vp_660_Pds * GAMMA * Vp_Vs_ratio_depth_2

					lst_true_660_mean_Pds = (H_a_Pds_660) * ((Vp_depth_2-Vs_depth_2)/(Vp_depth_2*Vs_depth_2)) * (((Vs_depth_2+lst_delta_Vs_660_Pds)*(Vp_depth_2+lst_delta_Vp_660_Pds))/(Vp_depth_2+lst_delta_Vp_660_Pds-Vs_depth_2-lst_delta_Vs_660_Pds))

					lst_true_660_mean_Ppds = (H_a_Ppds_660) * ((Vp_depth_2+Vs_depth_2)/(Vp_depth_2*Vs_depth_2)) * (((Vs_depth_2+lst_delta_Vs_660_Pds)*(Vp_depth_2+lst_delta_Vp_660_Pds))/(Vp_depth_2+lst_delta_Vp_660_Pds+Vs_depth_2+lst_delta_Vs_660_Pds))
			
					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_660_mean'] = lst_delta_Vp_660_Pds

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_660_mean'] = lst_delta_Vs_660_Pds

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_660_mean'] = lst_true_660_mean_Pds

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_660_mean'] = lst_true_660_mean_Ppds

					print('Equation 10 (Gao & Liu,2013) is true? = '+str(round(H_a_Pds_660/H_a_Ppds_660,2) == round(((alfa+lst_delta_Vp_660_Pds-lst_delta_Vs_660_Pds)/(beta+lst_delta_Vp_660_Pds+lst_delta_Vs_660_Pds))*(beta/alfa),2)))

					print('Delta Vp 660 km = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_660_mean']))
					print('Delta Vs 660 km = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_660_mean']))

					print('660 km Pds True Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_660_mean']))
					print('660 km Ppds True Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_660_mean']))

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_thickness_MTZ_mean'] = lst_true_660_mean_Pds - lst_true_410_mean_Pds

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_thickness_MTZ_mean'] = lst_true_660_mean_Ppds - lst_true_410_mean_Ppds

					print('MTZ Pds true thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_thickness_MTZ_mean']))
					print('MTZ Ppds true thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_thickness_MTZ_mean']))

					######## Estimating MTZ difference True and Apparent thickness ########

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ'] = (lst_true_660_mean_Pds - lst_true_410_mean_Pds)  - (lst_660_depth_Pds - lst_410_depth_Pds)
					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ'] = (lst_true_660_mean_Ppds - lst_true_410_mean_Ppds) - (lst_660_depth_Ppds - lst_410_depth_Ppds)

					print('MTZ Pds diff True and Apparent thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ']))
					print('MTZ Ppds diff True and Apparent thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ']))

					######## Estimating MTZ difference True and Model thickness ########

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ_model'] = (lst_true_660_mean_Pds - lst_true_410_mean_Pds)  - 250
					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ_model'] = (lst_true_660_mean_Ppds - lst_true_410_mean_Ppds) - 250

					print('MTZ Pds diff True and Model thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ_model']))
					print('MTZ Ppds diff True and Model thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ_model']))
					print('\n')
				
				else:
											
					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_410_mean'] = np.nan

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_410_mean'] = np.nan

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_410_mean'] = np.nan

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_410_mean'] = np.nan

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_660_mean'] = np.nan

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_660_mean'] = np.nan

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_660_mean'] = np.nan

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_660_mean'] = np.nan

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_thickness_MTZ_mean'] = np.nan

					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_thickness_MTZ_mean'] = np.nan

				
					######## Estimating MTZ difference True and Apparent thickness ########

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ'] = np.nan
					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ'] = np.nan

					######## Estimating MTZ difference True and Model thickness ########

					RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ_model'] = np.nan
					RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ_model'] = np.nan
					print('\n')

###################################################################################################################

print('Allocating results and Stacking Pds and Ppds data')

RF_lat = []
RF_lon = []

RF_lat_true = []
RF_lon_true = []

RF_DEPTH_mean_520_Pds = []
RF_DEPTH_std_520_Pds = []
RF_DEPTH_mean_520_Ppds = []
RF_DEPTH_std_520_Ppds = []

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

difference_thickness_MTZ_model_Pds = []
difference_thickness_MTZ_model_Pds_std = []

difference_thickness_MTZ_model_Ppds = []
difference_thickness_MTZ_model_Ppds_std  = []

RF_BOOTSTRAP_DATA_Pds = []
RF_BOOTSTRAP_DATA_Ppds = []

RF_BOOTSTRAP_DEPTH_mean_1_Pds = []
RF_BOOTSTRAP_DEPTH_mean_1_Ppds = []

RF_BOOTSTRAP_DEPTH_mean_520_Pds = []
RF_BOOTSTRAP_DEPTH_mean_520_Ppds = []

RF_BOOTSTRAP_DEPTH_mean_2_Pds = []
RF_BOOTSTRAP_DEPTH_mean_2_Ppds = []

len_RF_stacking_Pds = []
RF_stacking_Pds = []

len_RF_stacking_Ppds = []
RF_stacking_Ppds = []


for i,j in enumerate(RF_data_raw_Pds):
	if len(j) >= NUMBER_PP_PER_BIN and RF_STA_number_raw[i] >= NUMBER_STA_PER_BIN:
		stacking_Pds_data = [sum(x)/len(j)  for x in zip(*j)]
		RF_stacking_Pds.append(stacking_Pds_data)
		len_RF_stacking_Pds.append(len(j))

		stacking_Ppds_data = [sum(x)/len(RF_data_raw_Ppds[i])  for x in zip(*RF_data_raw_Ppds[i])]
		RF_stacking_Ppds.append(stacking_Ppds_data)
		len_RF_stacking_Ppds.append(len(RF_data_raw_Ppds[i]))

		RF_lat.append(grid_sel_y[i])
		RF_lon.append(grid_sel_x[i])

		flat_DATA_list_Pds = [RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_DATA'] for _k in range(BOOTSTRAP_INTERATOR)]
		flat_DATA_list_Ppds = [RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['RF_DATA'] for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DATA_Pds.append(flat_DATA_list_Pds)
		RF_BOOTSTRAP_DATA_Ppds.append(flat_DATA_list_Ppds)

		#Analysing stacked data amplitude in d410Pds

		flat_mean_1_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_1_Pds.append(flat_mean_1_Pds)

		lst_stacking_data_410_Pds = [stacking_Pds_data[x] for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
		d410Pds_candidate = [abs(np.nanmean(flat_mean_1_Pds) - c) for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
		
		amp_d410Pds = lst_stacking_data_410_Pds[d410Pds_candidate.index(min(d410Pds_candidate))]

		
		if  amp_d410Pds >= MIN_AMP_GOOD:
			
			RF_DEPTH_mean_1_Pds.append(np.nanmean(flat_mean_1_Pds))
			RF_DEPTH_std_1_Pds.append(np.nanstd(flat_mean_1_Pds))

		else: 

			RF_DEPTH_mean_1_Pds.append(np.nan)
			RF_DEPTH_std_1_Pds.append(np.nan)


		#Analysing stacked data amplitude in d410Ppds

		flat_mean_1_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_1_Ppds.append(flat_mean_1_Ppds)		

		lst_stacking_data_410_Ppds = [stacking_Ppds_data[x] for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
		d410Ppds_candidate = [abs(np.nanmean(flat_mean_1_Ppds) - c) for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]

		amp_d410Ppds = lst_stacking_data_410_Ppds[d410Ppds_candidate.index(min(d410Ppds_candidate))]

		if  amp_d410Ppds >= MIN_AMP_GOOD:

			RF_DEPTH_mean_1_Ppds.append(np.nanmean(flat_mean_1_Ppds))
			RF_DEPTH_std_1_Ppds.append(np.nanstd(flat_mean_1_Ppds))

		else: 

			RF_DEPTH_mean_1_Ppds.append(np.nan)
			RF_DEPTH_std_1_Ppds.append(np.nan)


		#Analysing stacked data amplitude in d520Pds

		flat_mean_520_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['520_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_520_Pds.append(flat_mean_520_Pds)

		lst_stacking_data_520_Pds = [stacking_Pds_data[x] for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]
		d520Pds_candidate = [abs(np.nanmean(flat_mean_520_Pds) - c) for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]
		
		amp_d520Pds = lst_stacking_data_520_Pds[d520Pds_candidate.index(min(d520Pds_candidate))]

		
		if  amp_d520Pds >= MIN_AMP_GOOD:
			
			RF_DEPTH_mean_520_Pds.append(np.nanmean(flat_mean_520_Pds))
			RF_DEPTH_std_520_Pds.append(np.nanstd(flat_mean_520_Pds))

		else: 

			RF_DEPTH_mean_520_Pds.append(np.nan)
			RF_DEPTH_std_520_Pds.append(np.nan)


		#Analysing stacked data amplitude in d520Ppds

		flat_mean_520_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['520_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_520_Ppds.append(flat_mean_520_Ppds)		

		lst_stacking_data_520_Ppds = [stacking_Ppds_data[x] for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]
		d520Ppds_candidate = [abs(np.nanmean(flat_mean_520_Ppds) - c) for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]

		amp_d520Ppds = lst_stacking_data_520_Ppds[d520Ppds_candidate.index(min(d520Ppds_candidate))]

		if  amp_d520Ppds >= MIN_AMP_GOOD:

			RF_DEPTH_mean_520_Ppds.append(np.nanmean(flat_mean_520_Ppds))
			RF_DEPTH_std_520_Ppds.append(np.nanstd(flat_mean_520_Ppds))

		else: 

			RF_DEPTH_mean_520_Ppds.append(np.nan)
			RF_DEPTH_std_520_Ppds.append(np.nan)

		#Analysing stacked data amplitude in d660Pds

		flat_mean_2_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_2_Pds.append(flat_mean_2_Pds)

		lst_stacking_data_660_Pds = [stacking_Pds_data[x] for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
		d660Pds_candidate = [abs(np.nanmean(flat_mean_2_Pds) - c) for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
		
		amp_d660Pds = lst_stacking_data_660_Pds[d660Pds_candidate.index(min(d660Pds_candidate))]


		if  amp_d660Pds >= MIN_AMP_GOOD:

			counts = mode(flat_mean_2_Pds).count
			RF_DEPTH_mean_2_Pds.append(np.nanmean(flat_mean_2_Pds))
			RF_DEPTH_std_2_Pds.append(np.nanstd(flat_mean_2_Pds))

		else: 

			RF_DEPTH_mean_2_Pds.append(np.nan)
			RF_DEPTH_std_2_Pds.append(np.nan)

		#Analysing stacked data amplitude in d660Ppds

		flat_mean_2_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_2_Ppds.append(flat_mean_2_Ppds)

		lst_stacking_data_660_Ppds = [stacking_Ppds_data[x] for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
		d660Ppds_candidate = [abs(np.nanmean(flat_mean_2_Ppds) - c) for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
		
		amp_d660Ppds = lst_stacking_data_660_Pds[d660Ppds_candidate.index(min(d660Ppds_candidate))]


		if  amp_d660Ppds >= MIN_AMP_GOOD:

			counts = mode(flat_mean_2_Ppds).count
			RF_DEPTH_mean_2_Ppds.append(np.nanmean(flat_mean_2_Ppds))
			RF_DEPTH_std_2_Ppds.append(np.nanstd(flat_mean_2_Ppds))

		else: 

			RF_DEPTH_mean_2_Ppds.append(np.nan)
			RF_DEPTH_std_2_Ppds.append(np.nan)


		#Analysing stacked data amplitude to calculate MTZ THICKNESS Pds


		if  amp_d410Pds >= MIN_AMP_GOOD and amp_d660Pds >= MIN_AMP_GOOD:

			flat_thickness_MTZ_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			thickness_MTZ_Pds.append(np.nanmean(flat_thickness_MTZ_Pds))
			thickness_MTZ_Pds_std.append(np.nanstd(flat_thickness_MTZ_Pds))

		else:

			thickness_MTZ_Pds.append(np.nan)
			thickness_MTZ_Pds_std.append(np.nan)

		#Analysing stacked data amplitude to calculate MTZ THICKNESS Ppds

		if  amp_d410Ppds >= MIN_AMP_GOOD and amp_d660Ppds >= MIN_AMP_GOOD:

			flat_thickness_MTZ_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['thickness_MTZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			thickness_MTZ_Ppds.append(np.nanmean(flat_thickness_MTZ_Ppds))
			thickness_MTZ_Ppds_std.append(np.nanstd(flat_thickness_MTZ_Ppds))

		else:
			
			thickness_MTZ_Ppds.append(np.nan)
			thickness_MTZ_Ppds_std.append(np.nan)


		#Analysing stacked data amplitude to calculate TRUE MTZ THICKNESS


		if amp_d410Pds >= MIN_AMP_GOOD and amp_d660Pds >= MIN_AMP_GOOD and amp_d410Ppds >= MIN_AMP_GOOD and amp_d660Ppds >= MIN_AMP_GOOD:
			RF_lat_true.append(grid_sel_y[i])
			RF_lon_true.append(grid_sel_x[i])

			flat_delta_1_Vp_mean = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			delta_1_Vp_mean.append(np.nanmean(flat_delta_1_Vp_mean))
			delta_1_Vp_std.append(np.nanstd(flat_delta_1_Vp_mean))

			flat_delta_1_Vs_mean = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			delta_1_Vs_mean.append(np.nanmean(flat_delta_1_Vs_mean))
			delta_1_Vs_std.append(np.nanstd(flat_delta_1_Vs_mean))

			flat_mean_1_true_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			RF_DEPTH_mean_1_true_Pds.append(np.nanmean(flat_mean_1_true_Pds))
			RF_DEPTH_std_1_true_Pds.append(np.nanstd(flat_mean_1_true_Pds))

			flat_mean_1_true_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			RF_DEPTH_mean_1_true_Ppds.append(np.nanmean(flat_mean_1_true_Ppds))
			RF_DEPTH_std_1_true_Ppds.append(np.nanstd(flat_mean_1_true_Ppds))

			flat_delta_2_Vp = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vp_660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			delta_2_Vp_mean.append(np.nanmean(flat_delta_2_Vp))
			delta_2_Vp_std.append(np.nanstd(flat_delta_2_Vp))

			flat_delta_2_Vs = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['delta_Vs_660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			delta_2_Vs_mean.append(np.nanmean(flat_delta_2_Vs))
			delta_2_Vs_std.append(np.nanstd(flat_delta_2_Vs))

			flat_mean_2_true_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			RF_DEPTH_mean_2_true_Pds.append(np.nanmean(flat_mean_2_true_Pds))
			RF_DEPTH_std_2_true_Pds.append(np.nanstd(flat_mean_2_true_Pds))

			flat_mean_2_true_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			RF_DEPTH_mean_2_true_Ppds.append(np.nanmean(flat_mean_2_true_Ppds))
			RF_DEPTH_std_2_true_Ppds.append(np.nanstd(flat_mean_2_true_Ppds))

			flat_true_thickness_MTZ_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['true_thickness_MTZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			true_thickness_MTZ_Pds.append(np.nanmean(flat_true_thickness_MTZ_Pds))
			true_thickness_MTZ_Pds_std.append(np.nanstd(flat_true_thickness_MTZ_Pds))

			flat_true_thickness_MTZ_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['true_thickness_MTZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			true_thickness_MTZ_Ppds.append(np.nanmean(flat_true_thickness_MTZ_Ppds))
			true_thickness_MTZ_Ppds_std.append(np.nanstd(flat_true_thickness_MTZ_Ppds))

			flat_diff_thickness_MTZ_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ']) for _k in range(BOOTSTRAP_INTERATOR)]
			diff_thickness_MTZ_Pds.append(np.nanmean(flat_diff_thickness_MTZ_Pds))
			diff_thickness_MTZ_Pds_std.append(np.nanstd(flat_diff_thickness_MTZ_Pds))

			flat_diff_thickness_MTZ_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ']) for _k in range(BOOTSTRAP_INTERATOR)]
			diff_thickness_MTZ_Ppds.append(np.nanmean(flat_diff_thickness_MTZ_Ppds))
			diff_thickness_MTZ_Ppds_std.append(np.nanstd(flat_diff_thickness_MTZ_Ppds))

			flat_difference_thickness_MTZ_model_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['difference_thickness_MTZ_model']) for _k in range(BOOTSTRAP_INTERATOR)]
			difference_thickness_MTZ_model_Pds.append(np.nanmean(flat_difference_thickness_MTZ_model_Pds))
			difference_thickness_MTZ_model_Pds_std.append(np.nanstd(flat_difference_thickness_MTZ_model_Pds))

			flat_difference_thickness_MTZ_model_Ppds = [float(RF_BOOTSTRAP_ESTIMATION_Ppds[_k][i]['difference_thickness_MTZ_model']) for _k in range(BOOTSTRAP_INTERATOR)]
			difference_thickness_MTZ_model_Ppds.append(np.nanmean(flat_difference_thickness_MTZ_model_Ppds))
			difference_thickness_MTZ_model_Ppds_std.append(np.nanstd(flat_difference_thickness_MTZ_model_Ppds))
			
		else: 

			RF_lat_true.append(np.nan)
			RF_lon_true.append(np.nan)

			delta_1_Vp_mean.append(np.nan)
			delta_1_Vp_std.append(np.nan)
			delta_1_Vs_mean.append(np.nan)
			delta_1_Vs_std.append(np.nan)

			RF_DEPTH_mean_1_true_Pds.append(np.nan)
			RF_DEPTH_std_1_true_Pds.append(np.nan)
			RF_DEPTH_mean_1_true_Ppds.append(np.nan)
			RF_DEPTH_std_1_true_Ppds.append(np.nan)

			delta_2_Vp_mean.append(np.nan)
			delta_2_Vp_std.append(np.nan)
			delta_2_Vs_mean.append(np.nan)
			delta_2_Vs_std.append(np.nan)

			RF_DEPTH_mean_2_true_Pds.append(np.nan)
			RF_DEPTH_std_2_true_Pds.append(np.nan)
			RF_DEPTH_mean_2_true_Ppds.append(np.nan)
			RF_DEPTH_std_2_true_Ppds.append(np.nan)

			true_thickness_MTZ_Pds.append(np.nan)
			true_thickness_MTZ_Pds_std.append(np.nan)
			true_thickness_MTZ_Ppds.append(np.nan)
			true_thickness_MTZ_Ppds_std.append(np.nan)

			diff_thickness_MTZ_Pds.append(np.nan)
			diff_thickness_MTZ_Pds_std.append(np.nan)
			diff_thickness_MTZ_Ppds.append(np.nan)
			diff_thickness_MTZ_Ppds_std.append(np.nan)

			difference_thickness_MTZ_model_Pds.append(np.nan)
			difference_thickness_MTZ_model_Pds_std.append(np.nan)
			difference_thickness_MTZ_model_Ppds.append(np.nan)
			difference_thickness_MTZ_model_Ppds_std.append(np.nan)

#############################################################################################################################################################################################

print('Plotting: Figure Final Grid and Ppds Average Piercing Points')
print('\n')

fig_PP, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#Figure

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
l1, = ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.Geodetic())
l3, = ax.plot(pp_med_long_Ppds,pp_med_lat_Ppds, 'X',markersize=5,markeredgecolor='k',markerfacecolor='k',alpha=0.5,transform=ccrs.Geodetic())
for i,j in enumerate(RF_lon):
	circulo = Circle(radius=DIST_GRID_PP,xy=(RF_lon[i], RF_lat[i]),color='None', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
	ax.add_patch(circulo)
for i,j in enumerate(pp_med_long_Ppds):
	circulo_fresnel = Circle(radius=FRESNEL_ZONE_RADIUS, xy=(pp_med_long_Ppds[i],pp_med_lat_Ppds[i]), color='gray',linewidth=0,alpha=0.2,transform=ccrs.Geodetic(),zorder=1)
	ax.add_patch(circulo_fresnel)

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax.gridlines(draw_labels=True)

#ax.set_title('Figure: Final Grid and Ppds Average Piercing Points',ha='center',va='top',y=1.08)
ax.legend([l1,l3,circulo,circulo_fresnel],['Stations','Piercing Points '+str(DEPTH_TARGET),'Selected Grid','Piercing Points Fresnel Zone'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w',fontsize='smaller')

#plt.show()

fig_PP.savefig(RESULTS_FOLDER+'PP_FINAL_GRID.'+EXT_FIG,dpi=DPI_FIG)

#############################################################################################################################################################################################

print('Saving Selected Piercing Points in JSON file')
print('\n')

PP_SELEC_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'SELECTED_BINNED_DATA'+'/'

os.makedirs(PP_SELEC_DIR,exist_ok=True)

SELECTED_BINNED_DATA_dic = {
	'lat':[],'lon':[],
	'lat_true':[],'lon_true':[],
	'len_Pds':[],'len_Ppds':[],
	'true_mean_1_Pds':[],'true_std_1_Pds':[],
	'true_mean_2_Pds':[],'true_std_2_Pds':[],
	'true_mean_1_Ppds':[],'true_std_1_Ppds':[],
	'true_mean_2_Ppds':[],'true_std_2_Ppds':[],
	'mean_1_Pds':[],'std_1_Pds':[],
	'mean_520_Pds':[],'std_520_Pds':[],
	'mean_2_Pds':[],'std_2_Pds':[],
	'mean_1_Ppds':[],'std_1_Ppds':[],
	'mean_520_Ppds':[],'std_520_Ppds':[],
	'mean_2_Ppds':[],'std_2_Ppds':[],
	'delta_1_Vp_mean':[],'delta_1_Vp_std':[],
	'delta_2_Vp_mean':[],'delta_2_Vp_std':[],
	'delta_1_Vs_mean':[],'delta_1_Vs_std':[],
	'delta_2_Vs_mean':[],'delta_2_Vs_std':[],
	'mtz_thickness_Pds':[],'mtz_thickness_Pds_std':[],
	'mtz_thickness_Ppds':[],'mtz_thickness_Ppds_std':[],
	'true_thickness_MTZ_Pds':[],'true_thickness_MTZ_Pds_std':[],
	'true_thickness_MTZ_Ppds':[],'true_thickness_MTZ_Ppds_std':[],
	'difference_thickness_MTZ_Pds':[],'difference_thickness_MTZ_Pds_std':[],
	'difference_thickness_MTZ_Ppds':[],'difference_thickness_MTZ_Ppds_std':[],
	'data_Pds':[],'data_Ppds':[],
	'data_BOOTSTRAP_Pds':[],'data_BOOTSTRAP_Ppds':[],
	'RF_BOOTSTRAP_DEPTH_mean_1_Pds':[],'RF_BOOTSTRAP_DEPTH_mean_1_Ppds':[],
	'RF_BOOTSTRAP_DEPTH_mean_520_Pds':[],'RF_BOOTSTRAP_DEPTH_mean_520_Ppds':[],
	'RF_BOOTSTRAP_DEPTH_mean_2_Pds':[],'RF_BOOTSTRAP_DEPTH_mean_2_Ppds':[],
	'difference_thickness_MTZ_model_Pds':[],'difference_thickness_MTZ_model_Pds_std':[],
	'difference_thickness_MTZ_model_Ppds':[],'difference_thickness_MTZ_model_Ppds_std':[]
	}

for i,j in enumerate(RF_BOOTSTRAP_DATA_Pds):
	SELECTED_BINNED_DATA_dic['data_BOOTSTRAP_Pds'].append(j)
	SELECTED_BINNED_DATA_dic['data_BOOTSTRAP_Ppds'].append(RF_BOOTSTRAP_DATA_Ppds[i])

	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_1_Pds'].append(RF_BOOTSTRAP_DEPTH_mean_1_Pds[i])
	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_1_Ppds'].append(RF_BOOTSTRAP_DEPTH_mean_1_Ppds[i])
	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_520_Pds'].append(RF_BOOTSTRAP_DEPTH_mean_520_Pds[i])
	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_520_Ppds'].append(RF_BOOTSTRAP_DEPTH_mean_520_Ppds[i])
	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_2_Pds'].append(RF_BOOTSTRAP_DEPTH_mean_2_Pds[i])
	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_2_Ppds'].append(RF_BOOTSTRAP_DEPTH_mean_2_Ppds[i])

	SELECTED_BINNED_DATA_dic['data_Pds'].append(RF_stacking_Pds[i])
	SELECTED_BINNED_DATA_dic['data_Ppds'].append(RF_stacking_Ppds[i])

	SELECTED_BINNED_DATA_dic['lat'].append(float("%.4f" % round(RF_lat[i],4)))
	SELECTED_BINNED_DATA_dic['lon'].append(float("%.4f" % round(RF_lon[i],4)))

	SELECTED_BINNED_DATA_dic['lat_true'].append(float("%.4f" % round(RF_lat_true[i],4)))
	SELECTED_BINNED_DATA_dic['lon_true'].append(float("%.4f" % round(RF_lon_true[i],4)))

	SELECTED_BINNED_DATA_dic['len_Pds'].append(len_RF_stacking_Pds[i])
	SELECTED_BINNED_DATA_dic['len_Ppds'].append(len_RF_stacking_Ppds[i])

	SELECTED_BINNED_DATA_dic['true_mean_1_Pds'].append(float(RF_DEPTH_mean_1_true_Pds[i]))
	SELECTED_BINNED_DATA_dic['true_std_1_Pds'].append(float(RF_DEPTH_std_1_true_Pds[i]))
	SELECTED_BINNED_DATA_dic['true_mean_2_Pds'].append(float(RF_DEPTH_mean_2_true_Pds[i]))
	SELECTED_BINNED_DATA_dic['true_std_2_Pds'].append(float(RF_DEPTH_std_2_true_Pds[i]))

	SELECTED_BINNED_DATA_dic['mean_1_Pds'].append(float(RF_DEPTH_mean_1_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_1_Pds'].append(float(RF_DEPTH_std_1_Pds[i]))
	SELECTED_BINNED_DATA_dic['mean_520_Pds'].append(float(RF_DEPTH_mean_520_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_520_Pds'].append(float(RF_DEPTH_std_520_Pds[i]))
	SELECTED_BINNED_DATA_dic['mean_2_Pds'].append(float(RF_DEPTH_mean_2_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_2_Pds'].append(float(RF_DEPTH_std_2_Pds[i]))

	SELECTED_BINNED_DATA_dic['true_mean_1_Ppds'].append(float(RF_DEPTH_mean_1_true_Ppds[i]))
	SELECTED_BINNED_DATA_dic['true_std_1_Ppds'].append(float(RF_DEPTH_std_1_true_Ppds[i]))
	SELECTED_BINNED_DATA_dic['true_mean_2_Ppds'].append(float(RF_DEPTH_mean_2_true_Ppds[i]))
	SELECTED_BINNED_DATA_dic['true_std_2_Ppds'].append(float(RF_DEPTH_std_2_true_Ppds[i]))

	SELECTED_BINNED_DATA_dic['mean_1_Ppds'].append(float(RF_DEPTH_mean_1_Ppds[i]))
	SELECTED_BINNED_DATA_dic['std_1_Ppds'].append(float(RF_DEPTH_std_1_Ppds[i]))
	SELECTED_BINNED_DATA_dic['mean_520_Ppds'].append(float(RF_DEPTH_mean_520_Ppds[i]))
	SELECTED_BINNED_DATA_dic['std_520_Ppds'].append(float(RF_DEPTH_std_520_Ppds[i]))
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

	SELECTED_BINNED_DATA_dic['difference_thickness_MTZ_model_Pds'].append(float(difference_thickness_MTZ_model_Pds[i]))
	SELECTED_BINNED_DATA_dic['difference_thickness_MTZ_model_Pds_std'].append(float(difference_thickness_MTZ_model_Pds_std[i]))

	SELECTED_BINNED_DATA_dic['difference_thickness_MTZ_model_Ppds'].append(float(difference_thickness_MTZ_model_Ppds[i]))
	SELECTED_BINNED_DATA_dic['difference_thickness_MTZ_model_Ppds_std'].append(float(difference_thickness_MTZ_model_Ppds_std[i]))

RESULTS_FOLDER_BINS = PP_SELEC_DIR+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER_BINS,exist_ok=True)
with open(RESULTS_FOLDER_BINS+'SELECTED_BINNED.json', 'w') as fp:
	json.dump(SELECTED_BINNED_DATA_dic, fp)