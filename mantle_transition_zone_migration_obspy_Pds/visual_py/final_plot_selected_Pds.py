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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator, FixedLocator
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import pyarrow.feather as feather
import pandas as pd
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
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,
					EXT_FIG,DPI_FIG,DIST_GRID_PP,NUMBER_STA_PER_BIN,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,COLORMAP_STD,COLORMAP_VEL,DEPTH_TARGET
				   )


print('Starting Final Plot CODE')
print('\n')

print('Calculating earth model layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Stations'+'/'

print('Looking for Receiver Functions data in FEATHER file in '+STA_DIR)
print('\n')

filename_STA = STA_DIR+'sta_dic.feather'

sta_dic = pd.read_feather(filename_STA)  

event_depth = sta_dic['event_depth'].tolist()
event_lat = sta_dic['event_lat'].tolist()
event_long = sta_dic['event_long'].tolist()
event_dist = sta_dic['event_dist'].tolist()
event_gcarc = sta_dic['event_gcarc'].tolist()
event_sta = sta_dic['event_sta'].tolist()
event_ray = sta_dic['event_ray'].tolist()
sta_lat = sta_dic['sta_lat'].tolist()
sta_long = sta_dic['sta_long'].tolist()
sta_data = sta_dic['sta_data'].tolist()
sta_time = sta_dic['sta_time'].tolist()

print('Importing selected binned data')
print('\n')

PP_SELEC_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'SELECTED_BINNED_DATA'+'/'


RESULTS_FOLDER_BINS = PP_SELEC_DIR+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
filename = RESULTS_FOLDER_BINS+'SELECTED_BINNED.feather'

SELECTED_BINNED_DATA_dic = pd.read_feather(filename)  

lats = SELECTED_BINNED_DATA_dic['lat']
lons = SELECTED_BINNED_DATA_dic['lon']

RF_number = SELECTED_BINNED_DATA_dic['len_Pds']

RF_stacking_Pds = SELECTED_BINNED_DATA_dic['data_Pds']

RF_DEPTH_mean_LVZ_Pds = SELECTED_BINNED_DATA_dic['mean_LVZ_Pds']
RF_DEPTH_std_LVZ_Pds = SELECTED_BINNED_DATA_dic['std_LVZ_Pds']

RF_DEPTH_mean_1_Pds = SELECTED_BINNED_DATA_dic['mean_1_Pds']
RF_DEPTH_std_1_Pds = SELECTED_BINNED_DATA_dic['std_1_Pds']

RF_DEPTH_mean_520_Pds = SELECTED_BINNED_DATA_dic['mean_520_Pds']
RF_DEPTH_std_520_Pds = SELECTED_BINNED_DATA_dic['std_520_Pds']

RF_DEPTH_mean_2_Pds = SELECTED_BINNED_DATA_dic['mean_2_Pds']
RF_DEPTH_std_2_Pds = SELECTED_BINNED_DATA_dic['std_2_Pds']

RF_DEPTH_mean_LVZ_700_Pds = SELECTED_BINNED_DATA_dic['mean_LVZ_700_Pds']
RF_DEPTH_std_LVZ_700_Pds = SELECTED_BINNED_DATA_dic['std_LVZ_700_Pds']

RF_DEPTH_mtz_thickness_Pds = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds']
RF_DEPTH_mtz_thickness_Pds_std = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std']

PP_FIGURE = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Figures'+'/'

RESULTS_FOLDER = PP_FIGURE+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)


###################################################################################################################

print('Total of bins: '+str(len(RF_DEPTH_mean_1_Pds)))
print('\n')

print('Number of bins - 350 LVZ: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_LVZ_Pds))))
print('Max - 350 LVZ: '+str(np.nanmax(RF_DEPTH_mean_LVZ_Pds)))
print('Min - 350 LVZ: '+str(np.nanmin(RF_DEPTH_mean_LVZ_Pds)))
print(r'Bins - 350 LVZ depth: '+str(np.nanmean(RF_DEPTH_mean_LVZ_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_LVZ_Pds)))
print('\n')

print('Number of bins - 410 km depth: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_1_Pds))))
print('Max - 410 km depth: '+str(np.nanmax(RF_DEPTH_mean_1_Pds)))
print('Min - 410 km depth: '+str(np.nanmin(RF_DEPTH_mean_1_Pds)))
print(r'Bins - 410 km depth: '+str(np.nanmean(RF_DEPTH_mean_1_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_1_Pds)))
print('\n')

print('Number of bins - 520 km depth: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_520_Pds))))
print('Max - 520 km depth: '+str(np.nanmax(RF_DEPTH_mean_520_Pds)))
print('Min - 520 km depth: '+str(np.nanmin(RF_DEPTH_mean_520_Pds)))
print(r'Bins - 520 km depth: '+str(np.nanmean(RF_DEPTH_mean_520_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_520_Pds)))
print('\n')

print('Number of bins - 660 km depth: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_2_Pds))))
print('Max - 660 km depth: '+str(np.nanmax(RF_DEPTH_mean_2_Pds)))
print('Min - 660 km depth: '+str(np.nanmin(RF_DEPTH_mean_2_Pds)))
print(r'Bins - 660 km depth: '+str(np.nanmean(RF_DEPTH_mean_2_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_2_Pds)))
print('\n')

print('Number of bins - 700 LVZ: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_LVZ_700_Pds))))
print('Max - 700 LVZ: '+str(np.nanmax(RF_DEPTH_mean_LVZ_700_Pds)))
print('Min - 700 LVZ: '+str(np.nanmin(RF_DEPTH_mean_LVZ_700_Pds)))
print(r'Bins - 700 LVZ depth: '+str(np.nanmean(RF_DEPTH_mean_LVZ_700_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_LVZ_700_Pds)))
print('\n')

print('Number of bins - MTZ Thickness: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mtz_thickness_Pds))))
print('Max - MTZ Thickness: '+str(np.nanmax(RF_DEPTH_mtz_thickness_Pds)))
print('Min - MTZ Thickness: '+str(np.nanmin(RF_DEPTH_mtz_thickness_Pds)))
print(r'Bins - MTZ Thickness: '+str(np.nanmean(RF_DEPTH_mtz_thickness_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mtz_thickness_Pds)))
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

lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()

#############################################################################################################################################################################################
'''
print('Plotting Figure: Apparent Depth estimates')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,10))
	
idx_Pds = np.argsort(RF_DEPTH_mtz_thickness_Pds)

#bounds = np.arange(200, 300+colormap_segmentation, colormap_segmentation)
#norm_MTZ = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)
norm_MTZ = Normalize(vmin=200,vmax=300)
colors_MTZ_Pds = tmap(norm_MTZ(np.array(RF_DEPTH_mtz_thickness_Pds,dtype='float64')))

for j,i in enumerate(idx_Pds):
	if math.isnan(RF_DEPTH_mtz_thickness_Pds[i]) == False:
		ax.errorbar(RF_DEPTH_mean_2_Pds[i],RF_DEPTH_mean_1_Pds[i], yerr=RF_DEPTH_std_1_Pds[i],xerr=RF_DEPTH_std_2_Pds[i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
		ax.scatter(RF_DEPTH_mean_2_Pds[i],RF_DEPTH_mean_1_Pds[i],color=colors_MTZ_Pds[i],zorder=10)

ax.text(612,456,'N = '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mtz_thickness_Pds))),fontsize=12, fontweight='bold',bbox={'facecolor':'white','edgecolor':'none','pad':1})

ax.yaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

ax.set_ylim(360,460)
ax.set_xlim(610,710)
ax.set_ylabel('Depth d410 (km)')
ax.set_xlabel('Depth d660 (km)')
ax.tick_params(labeltop=True, labelright=True,labelsize=13)


sm_MTZ = plt.cm.ScalarMappable(cmap=tmap,norm=norm_MTZ)
sm_MTZ._A = []

cbar = fig.colorbar(sm_MTZ,ax=ax,orientation='horizontal',shrink=0.8,label='MTZ Thickness')

cbar.set_ticks(np.arange(200, 300+INTER_DEPTH, INTER_DEPTH))
cbar.set_ticklabels(np.arange(200, 300+INTER_DEPTH, INTER_DEPTH))

fig.savefig(RESULTS_FOLDER+'APPARENT_DEPTH_PLOT.'+EXT_FIG,dpi=DPI_FIG)

'''
#############################################################################################################################################################################################

print('Plotting Figure: Apparent Depth of 410 km and 660 km')

fig, axes = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(20,10),sharey=False)

ax = axes[0]
ax2 = axes[1]

#410 km

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

#bounds = np.arange(360, 460+colormap_segmentation,colormap_segmentation)
#norm_410 = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)
norm_410 = Normalize(vmin=360,vmax=460)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_1_Pds[i]) == False:
		if RF_DEPTH_std_1_Pds[i] < 10:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_410(RF_DEPTH_mean_1_Pds[i])), ec='k',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
		else:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_410(RF_DEPTH_mean_1_Pds[i])), ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree(),label='Stations')

ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)


#660 km

ax2.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax2.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax2.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

#bounds = np.arange(610, 710+colormap_segmentation, colormap_segmentation)
#norm_660 = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)
norm_660 = Normalize(vmin=610,vmax=710)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_2_Pds[i]) == False:
		if RF_DEPTH_std_1_Pds[i] < 10:
			circulo_660 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_660(RF_DEPTH_mean_2_Pds[i])),  ec='k',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax2.add_patch(circulo_660)
		else:
			circulo_660 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_660(RF_DEPTH_mean_2_Pds[i])),  ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax2.add_patch(circulo_660)

	else: 
		pass

ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree(),label='Stations')
ax2.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax2.grid(color='grey', linestyle='--', linewidth=1)
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.tick_params(labeltop=True, labelright=True,labelsize=15)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=tmap,norm=norm_410)
sm_410._A = []
cbar_410 = fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8,label='410 km Depth')

cbar_410.set_ticks(np.arange(360, 460+INTER_DEPTH, INTER_DEPTH))
cbar_410.set_ticklabels(np.arange(360, 460+INTER_DEPTH, INTER_DEPTH))

sm_660 = plt.cm.ScalarMappable(cmap=tmap,norm=norm_660)
sm_660._A = []

cbar_660 = fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8,label='660 km Depth')

cbar_660.set_ticks(np.arange(610, 710+INTER_DEPTH, INTER_DEPTH))
cbar_660.set_ticklabels(np.arange(610, 710+INTER_DEPTH, INTER_DEPTH))


fig.savefig(RESULTS_FOLDER+'Apparent_depth_Pds.'+EXT_FIG,dpi=DPI_FIG)

#############################################################################################################################################################################################

print('Plotting Figure: Thickness of the Mantle Transition Zone')

fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#Pds phase

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

#bounds = np.arange(200, 300+colormap_segmentation, colormap_segmentation)
#norm_MTZ = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)
norm_MTZ = Normalize(vmin=200,vmax=300)
colors_MTZ_Pds = tmap(norm_MTZ(np.array(RF_DEPTH_mtz_thickness_Pds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mtz_thickness_Pds[i]) == False:
		if RF_DEPTH_mtz_thickness_Pds_std[i] < 10:
			circulo_MTZ = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_MTZ(RF_DEPTH_mtz_thickness_Pds[i])),  ec='k',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_MTZ)
		else: 
			circulo_MTZ = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_MTZ(RF_DEPTH_mtz_thickness_Pds[i])),  ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_MTZ)
	else:
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())
ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)

#______________________________________________________________________

sm_MTZ = plt.cm.ScalarMappable(cmap=tmap,norm=norm_MTZ)
sm_MTZ._A = []
cbar_MTZ = fig.colorbar(sm_MTZ,ax=ax,orientation='horizontal',shrink=0.8,label='Thickness of MTZ')

cbar_MTZ.set_ticks(np.arange(200, 300+INTER_DEPTH, INTER_DEPTH))
cbar_MTZ.set_ticklabels(np.arange(200, 300+INTER_DEPTH, INTER_DEPTH))

fig.savefig(RESULTS_FOLDER+'THICKNESS_MTZ.'+EXT_FIG,dpi=DPI_FIG)


#############################################################################################################################################################################################


print('Plotting Figure: 520 Pds Depth')

fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#520 km Pds

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

#bounds = np.arange(470, 570+colormap_segmentation, colormap_segmentation)
#norm_410 = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)
norm_520 = Normalize(vmin=470,vmax=570)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_520_Pds[i]) == False:
		if RF_DEPTH_std_520_Pds[i] < 10:
			circulo_520 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_520(RF_DEPTH_mean_520_Pds[i])), ec='k',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_520)
		else:
			circulo_520 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_520(RF_DEPTH_mean_520_Pds[i])), ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_520)			
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())
ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)

#______________________________________________________________________

sm_520 = plt.cm.ScalarMappable(cmap=tmap,norm=norm_520)
sm_520._A = []
cbar_520 = fig.colorbar(sm_520,ax=ax,orientation='horizontal',shrink=0.8,label='520 km Depth')

cbar_520.set_ticks(np.arange(470, 570+INTER_DEPTH, INTER_DEPTH))
cbar_520.set_ticklabels(np.arange(470, 570+INTER_DEPTH, INTER_DEPTH))

fig.savefig(RESULTS_FOLDER+'520_APPARENT_DEPTH.'+EXT_FIG,dpi=DPI_FIG)

########################################################################################################################################################################


print('Plotting Figure: LVZ atop 410 km')

fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#LVZ

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

#bounds = np.arange(300, 400+colormap_segmentation, colormap_segmentation)
#norm_410 = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)
norm_410 = Normalize(vmin=300,vmax=400)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_LVZ_Pds[i]) == False:
		if RF_DEPTH_std_LVZ_Pds[i] < 10:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_410(RF_DEPTH_mean_LVZ_Pds[i])), ec='k',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
		else:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_410(RF_DEPTH_mean_LVZ_Pds[i])), ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())
ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=tmap,norm=norm_410)
sm_410._A = []
cbar_410 = fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8,label='LVZ atop 410 km')

cbar_410.set_ticks(np.arange(300, 400+INTER_DEPTH, INTER_DEPTH))
cbar_410.set_ticklabels(np.arange(300, 400+INTER_DEPTH, INTER_DEPTH))

fig.savefig(RESULTS_FOLDER+'LVZ_ATOP_410_KM.'+EXT_FIG,dpi=DPI_FIG)

########################################################################################################################################################################

print('Plotting Figure: LVZ 700 km')

fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#LVZ

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

#bounds = np.arange(650, 750+colormap_segmentation, colormap_segmentation)
norm_700 = Normalize(vmin=650,vmax=750)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_LVZ_700_Pds[i]) == False:
		if RF_DEPTH_std_LVZ_700_Pds[i] < 10:
			circulo_700 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_700(RF_DEPTH_mean_LVZ_700_Pds[i])), ec='k',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_700)
		else:
			circulo_700 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=tmap(norm_700(RF_DEPTH_mean_LVZ_700_Pds[i])), ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_700)
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())
ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+0,3), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+0,3), crs=ccrs.PlateCarree())
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labeltop=True, labelright=True,labelsize=15)

#______________________________________________________________________

sm_700 = plt.cm.ScalarMappable(cmap=tmap,norm=norm_700)
sm_700._A = []
cbar_700 = fig.colorbar(sm_700,ax=ax,orientation='horizontal',shrink=0.8,label='LVZ 700 km')

cbar_700.set_ticks(np.arange(650, 750+INTER_DEPTH, INTER_DEPTH))
cbar_700.set_ticklabels(np.arange(650, 750+INTER_DEPTH, INTER_DEPTH))

fig.savefig(RESULTS_FOLDER+'LVZ_700_KM.'+EXT_FIG,dpi=DPI_FIG)

print('Ending Final Plot CODE')