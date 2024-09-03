# coding: utf-8

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.legend_handler import HandlerPatch
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees, degrees2kilometers
import copy
import matplotlib
from matplotlib.cm import get_cmap
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature
import scipy.io
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
import json
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import interpolate
from shapely.geometry import Polygon, MultiPoint, Point, LinearRing
import shapefile
from matplotlib.colors import Normalize
from matplotlib.patches import Circle,Rectangle,Ellipse
import matplotlib.patches as mpatches

import math






from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					OUTPUT_DIR,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,
					EXT_FIG,DPI_FIG,DIST_GRID_PP,NUMBER_STA_PER_BIN,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,COLORMAP_STD,COLORMAP_VEL,DEPTH_MOHO,DEPTH_LAB
				   )


print('Starting Final Plot CODE')
print('\n')

print('Calculating earth model layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_MOHO_'+str(DEPTH_MOHO)+'_DEPTH_LAB_'+str(DEPTH_LAB)+'/'+'Stations'+'/'

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

print('Importing selected binned data')
print('\n')

PP_SELEC_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_MOHO_'+str(DEPTH_MOHO)+'_DEPTH_LAB_'+str(DEPTH_LAB)+'/'+'SELECTED_BINNED_DATA'+'/'

RESULTS_FOLDER_BINS = PP_SELEC_DIR+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
filename = RESULTS_FOLDER_BINS+'SELECTED_BINNED.json'

SELECTED_BINNED_DATA_dic = json.load(open(filename))

lats = SELECTED_BINNED_DATA_dic['lat']
lons = SELECTED_BINNED_DATA_dic['lon']

RF_number = SELECTED_BINNED_DATA_dic['len_Pds']

RF_stacking_Pds = SELECTED_BINNED_DATA_dic['data_Pds']

RF_DEPTH_mean_MOHO_Pds = SELECTED_BINNED_DATA_dic['mean_MOHO_Pds']
RF_DEPTH_std_MOHO_Pds = SELECTED_BINNED_DATA_dic['std_MOHO_Pds']

RF_DEPTH_mean_LAB_Pds = SELECTED_BINNED_DATA_dic['mean_LAB_Pds']
RF_DEPTH_std_LAB_Pds = SELECTED_BINNED_DATA_dic['std_LAB_Pds']

PP_FIGURE = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_MOHO_'+str(DEPTH_MOHO)+'_DEPTH_LAB_'+str(DEPTH_LAB)+'/'+'Figures'+'/'

RESULTS_FOLDER = PP_FIGURE+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)


###################################################################################################################

print('Total of bins: '+str(len(lats)))
print('\n')

print('Number of bins - MOHO: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_MOHO_Pds))))
print('Max - MOHO: '+str(np.nanmax(RF_DEPTH_mean_MOHO_Pds)))
print('Min - MOHO: '+str(np.nanmin(RF_DEPTH_mean_MOHO_Pds)))
print(r'Bins - MOHO: '+str(np.nanmean(RF_DEPTH_mean_MOHO_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_MOHO_Pds)))
print('\n')

print('Number of bins - LAB: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_LAB_Pds))))
print('Max - LAB: '+str(np.nanmax(RF_DEPTH_mean_LAB_Pds)))
print('Min - LAB: '+str(np.nanmin(RF_DEPTH_mean_LAB_Pds)))
print(r'Bins - LAB: '+str(np.nanmean(RF_DEPTH_mean_LAB_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_LAB_Pds)))
print('\n')

print('\n')

###################################################################################################################

#Color Maps

colormap = plt.get_cmap(COLORMAP_VEL)


colormap_std = plt.get_cmap(COLORMAP_STD)

colormap_segmentation = INTER_DEPTH/100

n=41
x = 0.5
upper = plt.cm.seismic(np.linspace(0, x, n)[::-1])
white = plt.cm.seismic(np.ones(18)*0.5)
lower= plt.cm.seismic(np.linspace(1-x, 1, n)[::-1])
colors = np.vstack((lower, white,upper))
tmap = matplotlib.colors.LinearSegmentedColormap.from_list('map_white', colors)

#############################################################################################################################################################################################

print('Plotting Figure: Thickness of MOHO')

fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
ax.gridlines(draw_labels=True)

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

#bounds = np.arange(30, 50+colormap_segmentation, colormap_segmentation)
#norm_MOHO = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)
norm_MOHO = Normalize(vmin=20,vmax=50)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_MOHO_Pds[i]) == False:
		if RF_DEPTH_std_MOHO_Pds[i] < INTER_DEPTH:
			circulo_MOHO = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_MOHO(RF_DEPTH_mean_MOHO_Pds[i])),  ec='k',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_MOHO)
		else: 
			circulo_MOHO = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_MOHO(RF_DEPTH_mean_MOHO_Pds[i])),  ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_MOHO)
	else:
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

#______________________________________________________________________

sm_MOHO = plt.cm.ScalarMappable(cmap=tmap,norm=norm_MOHO)
sm_MOHO._A = []
cbar_MOHO = fig.colorbar(sm_MOHO,ax=ax,orientation='horizontal',shrink=0.8,label='MOHO depth (km)')

cbar_MOHO.set_ticks(np.arange(20, 50+INTER_DEPTH, INTER_DEPTH))
cbar_MOHO.set_ticklabels(np.arange(20, 50+INTER_DEPTH, INTER_DEPTH))

fig.savefig(RESULTS_FOLDER+'THICKNESS_MOHO.'+EXT_FIG,dpi=DPI_FIG)


#############################################################################################################################################################################################


print('Plotting Figure: Thickness of LAB')

fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)
ax.gridlines(draw_labels=True)

#bounds = np.arange(100, 200+colormap_segmentation, colormap_segmentation)
#norm_LAB = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)
norm_LAB = Normalize(vmin=50,vmax=200)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_LAB_Pds[i]) == False:
		if RF_DEPTH_std_LAB_Pds[i] < INTER_DEPTH:
			circulo_LAB = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_LAB(RF_DEPTH_mean_LAB_Pds[i])), ec='k',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_LAB)
		else:
			circulo_LAB = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_LAB(RF_DEPTH_mean_LAB_Pds[i])), ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_LAB)			
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

#______________________________________________________________________

sm_LAB = plt.cm.ScalarMappable(cmap=tmap,norm=norm_LAB)
sm_LAB._A = []
cbar_LAB = fig.colorbar(sm_LAB,ax=ax,orientation='horizontal',shrink=0.8,label='LAB depth (km)')

cbar_LAB.set_ticks(np.arange(50, 200+INTER_DEPTH, INTER_DEPTH))
cbar_LAB.set_ticklabels(np.arange(50, 200+INTER_DEPTH, INTER_DEPTH))

fig.savefig(RESULTS_FOLDER+'THICKNESS_LAB.'+EXT_FIG,dpi=DPI_FIG)

print('Ending Final Plot CODE')