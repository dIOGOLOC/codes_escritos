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
from fatiando import gridder, utils
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

PP_SELEC_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'SELECTED_BINNED_DATA'+'/'

RESULTS_FOLDER_BINS = PP_SELEC_DIR+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
filename = RESULTS_FOLDER_BINS+'SELECTED_BINNED.json'

SELECTED_BINNED_DATA_dic = json.load(open(filename))

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

RF_DEPTH_mtz_thickness_Pds = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds']
RF_DEPTH_mtz_thickness_Pds_std = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std']

PP_FIGURE = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Figures'+'/'

RESULTS_FOLDER = PP_FIGURE+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)


###################################################################################################################

print('Total of bins: '+str(len(RF_DEPTH_mean_1_Pds)))
print('\n')

print('Pds Phases')
print('Number of bins - LVZ: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_LVZ_Pds))))
print(r'Bins - 410 km depth: '+str(np.nanmean(RF_DEPTH_mean_LVZ_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_LVZ_Pds)))

print('Number of bins - 410 km depth: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_1_Pds))))
print(r'Bins - 410 km depth: '+str(np.nanmean(RF_DEPTH_mean_1_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_1_Pds)))

print('Number of bins - 520 km depth: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_520_Pds))))
print(r'Bins - 520 km depth: '+str(np.nanmean(RF_DEPTH_mean_520_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_520_Pds)))

print('Number of bins - 660 km depth: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mean_2_Pds))))
print(r'Bins - 660 km depth: '+str(np.nanmean(RF_DEPTH_mean_2_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mean_2_Pds)))

print('Number of bins - MTZ Thickness: '+str(np.count_nonzero(~np.isnan(RF_DEPTH_mtz_thickness_Pds))))
print(r'Bins - MTZ Thickness: '+str(np.nanmean(RF_DEPTH_mtz_thickness_Pds))+' ± '+str(np.nanstd(RF_DEPTH_mtz_thickness_Pds)))
print('\n')

###################################################################################################################

#Color Maps

colormap = plt.get_cmap(COLORMAP_VEL)


colormap_std = plt.get_cmap(COLORMAP_STD)


#############################################################################################################################################################################################

print('Plotting Figure: Apparent Depth estimates (Pds phase)')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,10))
	
idx_Pds = np.argsort(RF_DEPTH_mtz_thickness_Pds)

norm_410 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)
colors_410_Pds = colormap(norm_410(np.array(RF_DEPTH_mtz_thickness_Pds,dtype='float64')))

for j,i in enumerate(idx_Pds):
	if math.isnan(RF_DEPTH_mtz_thickness_Pds[i]) == False:
		ax.errorbar(RF_DEPTH_mean_2_Pds[i],RF_DEPTH_mean_1_Pds[i], yerr=RF_DEPTH_std_1_Pds[i],xerr=RF_DEPTH_std_2_Pds[i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
		ax.scatter(RF_DEPTH_mean_2_Pds[i],RF_DEPTH_mean_1_Pds[i],color=colors_410_Pds[i],zorder=10)

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


sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []

fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8,label='MTZ Thickness Pds')


#plt.show()
fig.savefig(RESULTS_FOLDER+'APPARENT_DEPTH_PLOT.'+EXT_FIG,dpi=DPI_FIG)


#############################################################################################################################################################################################

print('Plotting Figure: Apparent Depth of 410 km and 660 km (Pds phase)')

fig, axes = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(20,10),sharey=True)

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
ax.gridlines(draw_labels=True)

norm_410 = mpl.colors.Normalize(vmin=360,vmax=460,clip=True)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_1_Pds[i]) == False:
		if RF_DEPTH_std_1_Pds[i] < 10:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_410(RF_DEPTH_mean_1_Pds[i])), ec='None',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
		else:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_410(RF_DEPTH_mean_1_Pds[i])), ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('410 km Pds', y=1.08)


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
ax2.gridlines(draw_labels=True)

norm_660 = mpl.colors.Normalize(vmin=610,vmax=710,clip=True)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_2_Pds[i]) == False:
		if RF_DEPTH_std_1_Pds[i] < 10:
			circulo_660 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_660(RF_DEPTH_mean_2_Pds[i])),  ec='None',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax2.add_patch(circulo_660)
		else:
			circulo_660 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_660(RF_DEPTH_mean_2_Pds[i])),  ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax2.add_patch(circulo_660)

	else: 
		pass

ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax2.set_title('660 km Pds', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

sm_660 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_660)
sm_660._A = []

fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8)

fig.savefig(RESULTS_FOLDER+'Apparent_depth_Pds.'+EXT_FIG,dpi=DPI_FIG)

#############################################################################################################################################################################################

print('Plotting Figure: Thickness of the Mantle Transition Zone (Pds Phase)')

fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,10))

#Pds phase

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

norm_410 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mtz_thickness_Pds[i]) == False:
		if RF_DEPTH_mtz_thickness_Pds_std[i] < 10:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_410(RF_DEPTH_mtz_thickness_Pds[i])),  ec='None',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
		else: 
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_410(RF_DEPTH_mtz_thickness_Pds[i])),  ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
	else:
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('Thickness of MTZ (Pds)', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

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
ax.gridlines(draw_labels=True)

norm_410 = mpl.colors.Normalize(vmin=470,vmax=570,clip=True)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_520_Pds[i]) == False:
		if RF_DEPTH_std_520_Pds[i] < 10:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_410(RF_DEPTH_mean_520_Pds[i])), ec='None',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
		else:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_410(RF_DEPTH_mean_520_Pds[i])), ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)			
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('520 km Pds', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)


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
ax.gridlines(draw_labels=True)

norm_410 = mpl.colors.Normalize(vmin=300,vmax=400,clip=True)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_LVZ_Pds[i]) == False:
		if RF_DEPTH_std_LVZ_Pds[i] < 10:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_410(RF_DEPTH_mean_LVZ_Pds[i])), ec='None',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
		else:
			circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i], lats[i]),color=colormap(norm_410(RF_DEPTH_mean_LVZ_Pds[i])), ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
			ax.add_patch(circulo_410)
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('LVZ atop 410 km', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)


fig.savefig(RESULTS_FOLDER+'LVZ_ATOP_410_KM.'+EXT_FIG,dpi=DPI_FIG)

print('Ending Final Plot CODE')