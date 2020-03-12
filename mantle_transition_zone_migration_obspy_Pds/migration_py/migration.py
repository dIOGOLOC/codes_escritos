# coding: utf-8


import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees, degrees2kilometers
import copy
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature
import shapefile
from scipy.stats import mode
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import json
import random
from matplotlib.colors import Normalize
from numpy import ma
from matplotlib import cbook
import collections
from matplotlib.collections import PatchCollection
import verde as vd
from shapely.geometry import Polygon, MultiPoint, Point
from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle,Rectangle
import math




from parameters_py.mgconfig import (
					RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					NUMBER_PP_PER_BIN,DEPTH_TARGET,RF_FREQUENCY,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,
					LLCRNRLON_SMALL,URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,
					PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,BOUNDARY_1_SHP,BOUNDARY_2_SHP,
					EXT_FIG,DPI_FIG,DIST_GRID_PP,NUMBER_STA_PER_BIN,OUTPUT_DIR,CONFIDENCE_BOUND,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,COLORMAP_STD,COLORMAP_VEL
								   )


print('Starting Receiver Functions migration code to estimate the depths of the Mantle discontinuities')
print('\n')

print('Importing depths and times of Pds conversion  dataset')
print('\n')

print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)


for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 410:
		Vp_depth_1 = j[3]
		Vs_depth_1 = j[5]

print('410 km earth model Vp : '+str(Vp_depth_1))
print('410 km earth model Vs : '+str(Vs_depth_1))
print('----')
#====================================================

for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 520:
		Vp_depth_520 = j[3]
		Vs_depth_520 = j[5]

print('520 km earth model Vp : '+str(Vp_depth_520))
print('520 km earth model Vs : '+str(Vs_depth_520))
print('----')
#====================================================
		
for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 660:
		Vp_depth_2 = j[3]
		Vs_depth_2 = j[5]

print('660 km earth model Vp : '+str(Vp_depth_2))
print('660 km earth model Vs : '+str(Vs_depth_2))
print('----')

#====================================================
print('Target Depth = '+str(DEPTH_TARGET))

for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == DEPTH_TARGET:
		Vp_depth_DEPTH_TARGET = j[3]
		Vs_depth_DEPTH_TARGET = j[5]

print('Target Depth earth model Vp : '+str(Vp_depth_DEPTH_TARGET))
print('Target Depth earth model Vs : '+str(Vs_depth_DEPTH_TARGET))
#====================================================

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
pp_time_DEPTH_TARGET  = [[]]*len(PP_lon_med)

for i,j in enumerate(PP_lon_med):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_med[i][k] == DEPTH_TARGET:
				pp_med_lat[i] = PP_lat_med[i][k]
				pp_med_long[i] = l
				pp_time_DEPTH_TARGET[i] = PP_time_med[i][-1] - PP_time_med[i][k]

pp_med_lat = [i for i in pp_med_lat if type(i) == float ]
pp_med_long = [i for i in pp_med_long if type(i) == float ]
pp_time_DEPTH_TARGET = [i for i in pp_time_DEPTH_TARGET if type(i) == float ]

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


print('Calculating MEAN FIRST FRESNEL ZONE RADIUS')
print('DEPTH TARGET MEAN TIME: '+str(np.mean(pp_time_DEPTH_TARGET)))

FRESNEL_ZONE_RADIUS_km = (Vp_depth_DEPTH_TARGET/2)* np.sqrt(np.mean(pp_time_DEPTH_TARGET) / RF_FREQUENCY)
FRESNEL_ZONE_RADIUS = kilometer2degrees(FRESNEL_ZONE_RADIUS_km)
print('MEAN FIRST FRESNEL ZONE RADIUS : '+str(FRESNEL_ZONE_RADIUS))
print('\n')

print('Creating GRID POINTS')
print('\n')

spacing = GRID_PP_MULT
west, east, south, north = LLCRNRLON_SMALL, URCRNRLON_SMALL, LLCRNRLAT_SMALL, URCRNRLAT_SMALL 
region = (west, east, south, north) 
rows, cols = vd.grid_coordinates(region=region, spacing=spacing)

grdx = rows.ravel()
grdy = cols.ravel()

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

dist_pp_grid_med = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_med[i] = [np.sqrt((j - pp_med_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_med_lat)]
    dist_pp_grid_med[i] = [np.sqrt((j - pp_med_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_med_lat)]


grid_sel_med = []
for i,j in enumerate(dist_pp_grid_med):
	vect_j = np.array(j)
	indices = vect_j.argsort()
	if vect_j[indices[NUMBER_PP_PER_BIN]] <= FRESNEL_ZONE_RADIUS:
		grid_sel_med.append((grdx[i],grdy[i]))

grid_sel = grid_sel_med

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
lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()
###################################################################################################################

print('Plotting: Figure Pds Piercing Points')
print('\n')


fig_PP_Pds, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10),sharey=True)

#Figure Pds

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
l2, = ax.plot(pp_1_long,pp_1_lat, '.',markersize=13,markeredgecolor='k',markerfacecolor='b',transform=ccrs.Geodetic())
l3, = ax.plot(pp_med_long,pp_med_lat, '.',markersize=10,markeredgecolor='k',markerfacecolor='g',transform=ccrs.Geodetic())
l4, = ax.plot(pp_2_long,pp_2_lat, '.',markersize=10,markeredgecolor='k',markerfacecolor='r',transform=ccrs.Geodetic())
l1, = ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.Geodetic())

legend = ax.legend([l1,l2,l3,l4],['Stations','Piercing Points 410 km','Piercing Points '+str(DEPTH_TARGET)+' km','Piercing Points 660 km'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)

PP_FIGURE = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Figures'+'/'

RESULTS_FOLDER = PP_FIGURE+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)

print('Folder to save Figures files:')
print(RESULTS_FOLDER)
print('\n')

fig_PP_Pds.savefig(RESULTS_FOLDER+'PP_Pds_no_GRID.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################
print('Plotting: Figure Pds Piercing Points')
print('\n')


fig_PP_Pds, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10),sharey=True)

#Figure Pds

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
l2, = ax.plot(pp_1_long,pp_1_lat, '.',markersize=13,markeredgecolor='k',markerfacecolor='b',transform=ccrs.Geodetic())
l3, = ax.plot(pp_med_long,pp_med_lat, '.',markersize=10,markeredgecolor='k',markerfacecolor='g',transform=ccrs.Geodetic())
l4, = ax.plot(pp_2_long,pp_2_lat, '.',markersize=10,markeredgecolor='k',markerfacecolor='r',transform=ccrs.Geodetic())
l1, = ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.Geodetic())

for i,j in enumerate(grdx):
	circulo = Circle(radius=DIST_GRID_PP,xy=(grdx[i], grdy[i]),color='None', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
	ax.add_patch(circulo)
legend = ax.legend([l1,l2,l3,l4,circulo],['Stations','Piercing Points 410 km','Piercing Points '+str(DEPTH_TARGET)+' km','Piercing Points 660 km','Selected Grid'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)

fig_PP_Pds.savefig(RESULTS_FOLDER+'PP_Pds.'+EXT_FIG,dpi=DPI_FIG)


###################################################################################################################

print('Plotting: Figure Pds Depth Target Piercing Points')
print('\n')


fig_PP, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#Figure Pds

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
l1, = ax.plot(sta_long,sta_lat, '^',markersize=13,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.Geodetic())
l3, = ax.plot(pp_med_long,pp_med_lat, 'X',markersize=10,markeredgecolor='k',markerfacecolor='r',alpha=0.5,transform=ccrs.Geodetic())

for i,j in enumerate(grdx):
	circulo = Circle(radius=DIST_GRID_PP,xy=(grdx[i], grdy[i]),color='None', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
	ax.add_patch(circulo)
for i,j in enumerate(pp_med_long):
	circulo_fresnel = Circle(radius=FRESNEL_ZONE_RADIUS, xy=(pp_med_long[i],pp_med_lat[i]), color='orange',linewidth=0,alpha=0.2,transform=ccrs.Geodetic(),zorder=1)
	ax.add_patch(circulo_fresnel)

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)

legend = ax.legend([l1,l3,circulo,circulo_fresnel],['Stations','Piercing Points '+str(DEPTH_TARGET)+' km','Selected Grid','Piercing Points Fresnel Zone'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)
fig_PP.savefig(RESULTS_FOLDER+'PP_MED_Pds.'+EXT_FIG,dpi=DPI_FIG)


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

##################################################################################################################

print('Filtering migrated data per grid data raw')
print('\n')

dados_grid_lat = pp_med_lat
dados_grid_lon = pp_med_long

number_PP_per_bin_raw = [[]]*len(grdx)

for i,j in enumerate(grdx):
	number_PP_per_bin_raw[i] = [RF_amplitude_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grdy[i] - l)**2) <= FRESNEL_ZONE_RADIUS]


#############################################################################################################################################################################################

print('Plotting: Figure Piercing Points per bin')
print('\n')

fig_PP, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#Figure

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

number_RF_per_bin = [len(i) for i in number_PP_per_bin_raw]

colormap = plt.get_cmap("magma")


norm_number = mpl.colors.Normalize(vmin=min(number_RF_per_bin),vmax=max(number_RF_per_bin),clip=True)
colors_number = colormap(norm_number(np.array(number_RF_per_bin,dtype='float64')))

for i,j in enumerate(number_RF_per_bin):
	if j > 9:
		circulo = Circle(radius=DIST_GRID_PP,xy=(grdx[i], grdy[i]),color=colors_number[i], ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(circulo)


reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)

sm_number = plt.cm.ScalarMappable(cmap=colormap,norm=norm_number)
sm_number._A = []

fig_PP.colorbar(sm_number,ax=ax,orientation='horizontal',shrink=0.8,label='Number of Piercing Points per Bin')

fig_PP.savefig(RESULTS_FOLDER+'NUMBER_PP_PER_BIN.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################


print('Filtering migrated data per selected grid data')
print('\n')

dados_grid_lat = pp_med_lat
dados_grid_lon = pp_med_long

RF_data_raw_Pds = [[]]*len(grid_sel_x)
RF_amplitude_depth_raw_Pds = [[]]*len(grid_sel_x)

RF_STA_number_raw = [[]]*len(grid_sel_x)

for i,j in enumerate(grid_sel_x):
	RF_data_raw_Pds[i] = [RF_amplitude_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) <= FRESNEL_ZONE_RADIUS]
	RF_amplitude_depth_raw_Pds[i] = [RF_amplitude_depth_Pds[k] for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) <= FRESNEL_ZONE_RADIUS]

	RF_STA_number_raw[i] = len(set([(sta_long[k],sta_lat[k]) for k,l in enumerate(dados_grid_lat) if np.sqrt((j - dados_grid_lon[k])**2 + (grid_sel_y[i] - l)**2) <= FRESNEL_ZONE_RADIUS]))

###################################################################################################################

print('Estimating Mean and Standard Deviation for each discontinuity')
print('\n')


### Creating Dictionaries to allocate results ###
def nested_dict():
	return collections.defaultdict(nested_dict)

RF_BOOTSTRAP_ESTIMATION_Pds = nested_dict()

for _k in range(BOOTSTRAP_INTERATOR):
	print('Bootstrap estimation '+str(_k+1))

	for i,j in enumerate(RF_data_raw_Pds):
		if len(j) >= NUMBER_PP_PER_BIN and RF_STA_number_raw[i] >= NUMBER_STA_PER_BIN:

			print('Grid point number: '+str(i))
			print('lat: '+str(grid_sel_y[i]))
			print('lon: '+str(grid_sel_x[i]))

			RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lon'] = grid_sel_x[i]
			RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['lat'] = grid_sel_y[i]

			old_RF_DATA_raw_Pds_lst = np.arange(0,len(j))
			new_RANDOM_RF_DATA_raw_Pds_lst = np.random.choice(old_RF_DATA_raw_Pds_lst,size=len(old_RF_DATA_raw_Pds_lst),replace=True)
				
			new_RANDOM_RF_DATA_raw_Pds = [j[_t_] for _t_ in new_RANDOM_RF_DATA_raw_Pds_lst]

			RF_STACKING_BOOTSTRAP_Pds = [sum(i)/len(new_RANDOM_RF_DATA_raw_Pds) for i in zip(*new_RANDOM_RF_DATA_raw_Pds)]

			RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_DATA'] = RF_STACKING_BOOTSTRAP_Pds

			######## Estimating LVZ atop the 410-km discontinuity  ########

			#LVZ atop 410 km

			lst_depth_amp_LVZ_Pds = [RF_STACKING_BOOTSTRAP_Pds[x] for x,c in enumerate(camadas_terra_10_km) if 350-(DEPTH_RANGE*2) <= c <= 350+(DEPTH_RANGE*2)]
			lst_depth_pp_LVZ_Pds = [c for x,c in enumerate(camadas_terra_10_km) if 350-(DEPTH_RANGE*2) <= c <= 350+(DEPTH_RANGE*2)]
			lst_LVZ_depth_Pds = lst_depth_pp_LVZ_Pds[lst_depth_amp_LVZ_Pds.index(min(lst_depth_amp_LVZ_Pds))]
			lst_LVZ_amp_Pds = lst_depth_amp_LVZ_Pds.index(min(lst_depth_amp_LVZ_Pds))

			######## Estimating 410 km apparent depth ########

			#410 km Pds

			lst_depth_amp_410_Pds = [RF_STACKING_BOOTSTRAP_Pds[x] for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
			lst_depth_pp_410_Pds = [c for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
			lst_410_depth_Pds = lst_depth_pp_410_Pds[lst_depth_amp_410_Pds.index(max(lst_depth_amp_410_Pds))]
			lst_410_amp_Pds = lst_depth_amp_410_Pds.index(max(lst_depth_amp_410_Pds))

			######## Estimating 520 km apparent depth ########

			#520 km Pds

			lst_depth_amp_520_Pds = [RF_STACKING_BOOTSTRAP_Pds[x] for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]
			lst_depth_pp_520_Pds = [c for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]
			lst_520_depth_Pds = lst_depth_pp_520_Pds[lst_depth_amp_520_Pds.index(max(lst_depth_amp_520_Pds))]
			lst_520_amp_Pds = lst_depth_amp_520_Pds.index(max(lst_depth_amp_520_Pds))

			######## Estimating 660 km apparent depth ########

			#660 km Pds

			lst_depth_amp_660_Pds = [RF_STACKING_BOOTSTRAP_Pds[x] for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
			lst_depth_pp_660_Pds = [c for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
			lst_660_depth_Pds = lst_depth_pp_660_Pds[lst_depth_amp_660_Pds.index(max(lst_depth_amp_660_Pds))]
			lst_660_amp_Pds = lst_depth_amp_660_Pds.index(max(lst_depth_amp_660_Pds))


			######## Estimating LVZ below the 660-km discontinuity  ########

			#LVZ below the 660-km

			lst_depth_amp_LVZ_700_Pds = [RF_STACKING_BOOTSTRAP_Pds[x] for x,c in enumerate(camadas_terra_10_km) if 700-(DEPTH_RANGE*2) <= c <= 700+(DEPTH_RANGE*2)]
			lst_depth_pp_LVZ_700_Pds = [c for x,c in enumerate(camadas_terra_10_km) if 700-(DEPTH_RANGE*2) <= c <= 700+(DEPTH_RANGE*2)]
			lst_LVZ_700_depth_Pds = lst_depth_pp_LVZ_700_Pds[lst_depth_amp_LVZ_700_Pds.index(min(lst_depth_amp_LVZ_700_Pds))]
			lst_LVZ_700_amp_Pds = lst_depth_amp_LVZ_700_Pds.index(min(lst_depth_amp_LVZ_700_Pds))
				
			######################################################################################################################################################

			if abs(lst_LVZ_amp_Pds) > 0:

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['LVZ_mean'] = lst_LVZ_depth_Pds
				
			else: 

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['LVZ_mean'] = np.nan
				
			print('LVZ Pds = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['LVZ_mean']))


			if lst_410_amp_Pds > 0:

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'] = lst_410_depth_Pds
				
			else: 

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean'] = np.nan
				
			print('410 Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean']))


			if lst_520_amp_Pds > 0:

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['520_mean'] = lst_520_depth_Pds
				
			else: 

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['520_mean'] = np.nan
				
			print('520 Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['520_mean']))

			if lst_660_depth_Pds > 0:

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'] = lst_660_depth_Pds
			
			else:

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean'] = np.nan

			print('660 km Pds Depth = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean']))


			if abs(lst_LVZ_700_amp_Pds) > 0:

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['LVZ_700_mean'] = lst_LVZ_700_depth_Pds
				
			else: 

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['LVZ_700_mean'] = np.nan
				
			print('LVZ 700 Pds = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['LVZ_700_mean']))


			######## Estimating MTZ thickness and difference between MTZ and Model thickness ########

			if  lst_410_amp_Pds > 0  and lst_660_depth_Pds > 0:

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean'] = lst_660_depth_Pds - lst_410_depth_Pds

			else: 

				RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean'] = np.nan

			print('MTZ Pds thickness = '+str(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean']))
			print('\n')
				
###################################################################################################################

print('Allocating results and Stacking Pds data')

RF_lat = []
RF_lon = []

RF_DEPTH_mean_LVZ_Pds = []
RF_DEPTH_std_LVZ_Pds = []

RF_DEPTH_mean_520_Pds = []
RF_DEPTH_std_520_Pds = []

RF_DEPTH_mean_1_Pds = []
RF_DEPTH_std_1_Pds = []

RF_DEPTH_mean_2_Pds = []
RF_DEPTH_std_2_Pds = []

RF_DEPTH_mean_LVZ_700_Pds = []
RF_DEPTH_std_LVZ_700_Pds = []

thickness_MTZ_Pds = []
thickness_MTZ_Pds_std = []

RF_BOOTSTRAP_DATA_Pds = []
RF_BOOTSTRAP_DATA_Pds_std = []

RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds = []

RF_BOOTSTRAP_DEPTH_mean_1_Pds = []

RF_BOOTSTRAP_DEPTH_mean_520_Pds = []

RF_BOOTSTRAP_DEPTH_mean_2_Pds = []

RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds = []

len_RF_stacking_Pds = []
RF_stacking_Pds = []

for i,j in enumerate(RF_data_raw_Pds):
	if len(j) >= NUMBER_PP_PER_BIN and RF_STA_number_raw[i] >= NUMBER_STA_PER_BIN:
		
		stacking_Pds_data = [sum(x)/len(j)  for x in zip(*j)]
		RF_stacking_Pds.append(stacking_Pds_data)
		len_RF_stacking_Pds.append(len(j))

		RF_lat.append(grid_sel_y[i])
		RF_lon.append(grid_sel_x[i])

		#Analysing BOOTSTRAP data amplitude

		flat_DATA_list_Pds = [RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['RF_DATA'] for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DATA_Pds.append(flat_DATA_list_Pds)
		
		BOOTSTRAP_DATA_Pds_per_depth = [[]]*len(flat_DATA_list_Pds[0])
		for u,w in enumerate(flat_DATA_list_Pds):
			for z,q in enumerate(w):
				BOOTSTRAP_DATA_Pds_per_depth[z].append(q)

		BOOTSTRAP_DATA_Pds_std = [np.std(i)*CONFIDENCE_BOUND for i in BOOTSTRAP_DATA_Pds_per_depth]
		RF_BOOTSTRAP_DATA_Pds_std.append(BOOTSTRAP_DATA_Pds_std)

		#Analysing stacked data amplitude in LVZ atop 410 km

		flat_mean_LVZ_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['LVZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds.append(flat_mean_LVZ_Pds)

		lst_stacking_data_LVZ_Pds = [stacking_Pds_data[x] for x,c in enumerate(camadas_terra_10_km) if 400-(DEPTH_RANGE*2) <= c <= 400-(DEPTH_RANGE*2)]
		LVZ_candidate = [abs(np.nanmean(flat_mean_LVZ_Pds) - c) for x,c in enumerate(camadas_terra_10_km) if 400-(DEPTH_RANGE*2) <= c <= 400-(DEPTH_RANGE*2)]
		BOOTSTRAP_DATA_LVZ_std_lst = [BOOTSTRAP_DATA_Pds_std[x] for x,c in enumerate(camadas_terra_10_km) if 400-(DEPTH_RANGE*2) <= c <= 400-(DEPTH_RANGE*2)]

		amp_LVZ = lst_stacking_data_LVZ_Pds[LVZ_candidate.index(min(LVZ_candidate))]

		BOOTSTRAP_DATA_LVZ_std_amp = BOOTSTRAP_DATA_LVZ_std_lst[LVZ_candidate.index(min(LVZ_candidate))]
		
		if  abs(amp_LVZ)-(BOOTSTRAP_DATA_LVZ_std_amp*CONFIDENCE_BOUND) > 0:
			
			RF_DEPTH_mean_LVZ_Pds.append(np.nanmean(flat_mean_LVZ_Pds))
			RF_DEPTH_std_LVZ_Pds.append(np.nanstd(flat_mean_LVZ_Pds)*CONFIDENCE_BOUND)

		else: 

			RF_DEPTH_mean_LVZ_Pds.append(np.nan)
			RF_DEPTH_std_LVZ_Pds.append(np.nan)

		
		#Analysing stacked data amplitude in d410Pds

		flat_mean_1_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['410_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_1_Pds.append(flat_mean_1_Pds)

		lst_stacking_data_410_Pds = [stacking_Pds_data[x] for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]
		d410Pds_candidate = [abs(np.nanmean(flat_mean_1_Pds) - c) for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]	
		BOOTSTRAP_DATA_410_std_lst = [BOOTSTRAP_DATA_Pds_std[x] for x,c in enumerate(camadas_terra_10_km) if 410-DEPTH_RANGE <= c <= 410+DEPTH_RANGE]

		amp_d410Pds = lst_stacking_data_410_Pds[d410Pds_candidate.index(min(d410Pds_candidate))]
		BOOTSTRAP_DATA_410_std_amp = BOOTSTRAP_DATA_410_std_lst[d410Pds_candidate.index(min(d410Pds_candidate))]

		
		if  amp_d410Pds-(BOOTSTRAP_DATA_410_std_amp*CONFIDENCE_BOUND) > 0:
				
			RF_DEPTH_mean_1_Pds.append(np.nanmean(flat_mean_1_Pds))
			RF_DEPTH_std_1_Pds.append(np.nanstd(flat_mean_1_Pds)*CONFIDENCE_BOUND)

		else: 

			RF_DEPTH_mean_1_Pds.append(np.nan)
			RF_DEPTH_std_1_Pds.append(np.nan)


		#Analysing stacked data amplitude in d520Pds

		flat_mean_520_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['520_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_520_Pds.append(flat_mean_520_Pds)

		lst_stacking_data_520_Pds = [stacking_Pds_data[x] for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]
		d520Pds_candidate = [abs(np.nanmean(flat_mean_520_Pds) - c) for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]
		BOOTSTRAP_DATA_520_std_lst = [BOOTSTRAP_DATA_Pds_std[x] for x,c in enumerate(camadas_terra_10_km) if 520-(DEPTH_RANGE*2) <= c <= 520+(DEPTH_RANGE*2)]

		amp_d520Pds = lst_stacking_data_520_Pds[d520Pds_candidate.index(min(d520Pds_candidate))]
		BOOTSTRAP_DATA_520_std_amp = BOOTSTRAP_DATA_520_std_lst[d520Pds_candidate.index(min(d520Pds_candidate))]

		
		if  amp_d520Pds-(BOOTSTRAP_DATA_520_std_amp*CONFIDENCE_BOUND) > 0:
					
			RF_DEPTH_mean_520_Pds.append(np.nanmean(flat_mean_520_Pds))
			RF_DEPTH_std_520_Pds.append(np.nanstd(flat_mean_520_Pds)*CONFIDENCE_BOUND)

		else: 

			RF_DEPTH_mean_520_Pds.append(np.nan)
			RF_DEPTH_std_520_Pds.append(np.nan)


		#Analysing stacked data amplitude in d660Pds

		flat_mean_2_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['660_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_2_Pds.append(flat_mean_2_Pds)

		lst_stacking_data_660_Pds = [stacking_Pds_data[x] for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
		d660Pds_candidate = [abs(np.nanmean(flat_mean_2_Pds) - c) for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]
		BOOTSTRAP_DATA_Pds_660_lst = [BOOTSTRAP_DATA_Pds_std[x] for x,c in enumerate(camadas_terra_10_km) if 660-DEPTH_RANGE <= c <= 660+DEPTH_RANGE]

		amp_d660Pds = lst_stacking_data_660_Pds[d660Pds_candidate.index(min(d660Pds_candidate))]
		BOOTSTRAP_DATA_660_std_amp = BOOTSTRAP_DATA_Pds_660_lst[d660Pds_candidate.index(min(d660Pds_candidate))]


		if  amp_d660Pds-(BOOTSTRAP_DATA_660_std_amp*CONFIDENCE_BOUND) > 0:

			RF_DEPTH_mean_2_Pds.append(np.nanmean(flat_mean_2_Pds))
			RF_DEPTH_std_2_Pds.append(np.nanstd(flat_mean_2_Pds)*CONFIDENCE_BOUND)

		else: 

			RF_DEPTH_mean_2_Pds.append(np.nan)
			RF_DEPTH_std_2_Pds.append(np.nan)

		#Analysing stacked data amplitude in LVZ below 700 km

		flat_mean_LVZ_700_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['LVZ_700_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
		RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds.append(flat_mean_LVZ_700_Pds)

		lst_stacking_data_LVZ_700_Pds = [stacking_Pds_data[x] for x,c in enumerate(camadas_terra_10_km) if 700-(DEPTH_RANGE*2) <= c <= 700-(DEPTH_RANGE*2)]
		LVZ_700_candidate = [abs(np.nanmean(flat_mean_LVZ_700_Pds) - c) for x,c in enumerate(camadas_terra_10_km) if 700-(DEPTH_RANGE*2) <= c <= 700-(DEPTH_RANGE*2)]
		BOOTSTRAP_DATA_LVZ_700_std_lst = [BOOTSTRAP_DATA_Pds_std[x] for x,c in enumerate(camadas_terra_10_km) if 700-(DEPTH_RANGE*2) <= c <= 700-(DEPTH_RANGE*2)]

		amp_LVZ_700 = lst_stacking_data_LVZ_700_Pds[LVZ_700_candidate.index(min(LVZ_700_candidate))]

		BOOTSTRAP_DATA_LVZ_700_std_amp = BOOTSTRAP_DATA_LVZ_700_std_lst[LVZ_700_candidate.index(min(LVZ_700_candidate))]
		
		if  abs(amp_LVZ_700)-(BOOTSTRAP_DATA_LVZ_700_std_amp*CONFIDENCE_BOUND) > 0:
			
			RF_DEPTH_mean_LVZ_700_Pds.append(np.nanmean(flat_mean_LVZ_700_Pds))
			RF_DEPTH_std_LVZ_700_Pds.append(np.nanstd(flat_mean_LVZ_700_Pds)*CONFIDENCE_BOUND)

		else: 

			RF_DEPTH_mean_LVZ_700_Pds.append(np.nan)
			RF_DEPTH_std_LVZ_700_Pds.append(np.nan)


		#Analysing stacked data amplitude to calculate MTZ THICKNESS Pds


		if  amp_d410Pds-(BOOTSTRAP_DATA_410_std_amp*CONFIDENCE_BOUND) > 0 and amp_d660Pds-(BOOTSTRAP_DATA_660_std_amp*CONFIDENCE_BOUND) > 0:

			flat_thickness_MTZ_Pds = [float(RF_BOOTSTRAP_ESTIMATION_Pds[_k][i]['thickness_MTZ_mean']) for _k in range(BOOTSTRAP_INTERATOR)]
			thickness_MTZ_Pds.append(np.nanmean(flat_thickness_MTZ_Pds))
			thickness_MTZ_Pds_std.append(np.nanstd(flat_thickness_MTZ_Pds)*CONFIDENCE_BOUND)

		else:

			thickness_MTZ_Pds.append(np.nan)
			thickness_MTZ_Pds_std.append(np.nan)

#############################################################################################################################################################################################

print('Plotting: Figure Final Grid and Piercing Points First Fresnel Zone')
print('\n')

fig_PP, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#Figure

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
l1, = ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.Geodetic())
l3, = ax.plot(pp_med_long,pp_med_lat, 'X',markersize=5,markeredgecolor='k',markerfacecolor='r',alpha=0.5,transform=ccrs.Geodetic())
for i,j in enumerate(RF_lon):
	circulo = Circle(radius=DIST_GRID_PP,xy=(RF_lon[i], RF_lat[i]),color='None', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
	ax.add_patch(circulo)
for i,j in enumerate(pp_med_long):
	circulo_fresnel = Circle(radius=FRESNEL_ZONE_RADIUS, xy=(pp_med_long[i],pp_med_lat[i]), color='orange',linewidth=0,alpha=0.2,transform=ccrs.Geodetic(),zorder=1)
	ax.add_patch(circulo_fresnel)

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)

legend = ax.legend([l1,l3,circulo,circulo_fresnel],['Stations','Piercing Points '+str(DEPTH_TARGET),'Selected Grid','Piercing Points Fresnel Zone'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)

fig_PP.savefig(RESULTS_FOLDER+'PP_FINAL_GRID.'+EXT_FIG,dpi=DPI_FIG)

#############################################################################################################################################################################################

print('Saving Selected Piercing Points in JSON file')
print('\n')

PP_SELEC_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'SELECTED_BINNED_DATA'+'/'

os.makedirs(PP_SELEC_DIR,exist_ok=True)

SELECTED_BINNED_DATA_dic = {
	'lat':[],'lon':[],
	'len_Pds':[],
	'mean_1_Pds':[],'std_1_Pds':[],
	'mean_LVZ_Pds':[],'std_LVZ_Pds':[],
	'mean_520_Pds':[],'std_520_Pds':[],
	'mean_2_Pds':[],'std_2_Pds':[],
	'mean_LVZ_700_Pds':[],'std_LVZ_700_Pds':[],

	'mtz_thickness_Pds':[],'mtz_thickness_Pds_std':[],

	'data_Pds':[],'data_BOOTSTRAP_Pds':[],'data_Pds_std':[],

	'RF_BOOTSTRAP_DEPTH_mean_1_Pds':[],
	'RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds':[],
	'RF_BOOTSTRAP_DEPTH_mean_520_Pds':[],
	'RF_BOOTSTRAP_DEPTH_mean_2_Pds':[],
	'RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds':[]


	}

for i,j in enumerate(RF_BOOTSTRAP_DATA_Pds):
	SELECTED_BINNED_DATA_dic['data_BOOTSTRAP_Pds'].append(j)

	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds'].append(RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds[i])

	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_1_Pds'].append(RF_BOOTSTRAP_DEPTH_mean_1_Pds[i])

	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_520_Pds'].append(RF_BOOTSTRAP_DEPTH_mean_520_Pds[i])

	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_2_Pds'].append(RF_BOOTSTRAP_DEPTH_mean_2_Pds[i])

	SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds'].append(RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds[i])

	SELECTED_BINNED_DATA_dic['data_Pds'].append(RF_stacking_Pds[i])
	SELECTED_BINNED_DATA_dic['data_Pds_std'].append(RF_BOOTSTRAP_DATA_Pds_std[i])


	SELECTED_BINNED_DATA_dic['lat'].append(float("%.4f" % round(RF_lat[i],4)))
	SELECTED_BINNED_DATA_dic['lon'].append(float("%.4f" % round(RF_lon[i],4)))

	SELECTED_BINNED_DATA_dic['len_Pds'].append(len_RF_stacking_Pds[i])


	SELECTED_BINNED_DATA_dic['mean_LVZ_Pds'].append(float(RF_DEPTH_mean_LVZ_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_LVZ_Pds'].append(float(RF_DEPTH_std_LVZ_Pds[i]))

	SELECTED_BINNED_DATA_dic['mean_1_Pds'].append(float(RF_DEPTH_mean_1_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_1_Pds'].append(float(RF_DEPTH_std_1_Pds[i]))

	SELECTED_BINNED_DATA_dic['mean_520_Pds'].append(float(RF_DEPTH_mean_520_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_520_Pds'].append(float(RF_DEPTH_std_520_Pds[i]))
	
	SELECTED_BINNED_DATA_dic['mean_2_Pds'].append(float(RF_DEPTH_mean_2_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_2_Pds'].append(float(RF_DEPTH_std_2_Pds[i]))

	SELECTED_BINNED_DATA_dic['mean_LVZ_700_Pds'].append(float(RF_DEPTH_mean_LVZ_700_Pds[i]))
	SELECTED_BINNED_DATA_dic['std_LVZ_700_Pds'].append(float(RF_DEPTH_std_LVZ_700_Pds[i]))

	SELECTED_BINNED_DATA_dic['mtz_thickness_Pds'].append(float(thickness_MTZ_Pds[i]))
	SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std'].append(float(thickness_MTZ_Pds_std[i]))
	
RESULTS_FOLDER_BINS = PP_SELEC_DIR+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER_BINS,exist_ok=True)
with open(RESULTS_FOLDER_BINS+'SELECTED_BINNED.json', 'w') as fp:
	json.dump(SELECTED_BINNED_DATA_dic, fp)