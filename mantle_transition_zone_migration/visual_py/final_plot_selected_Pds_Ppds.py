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
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature
from fatiando import gridder, utils
import scipy.io
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import json
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import interpolate
from shapely.geometry import Polygon, MultiPoint, Point, LinearRing
import shapefile
from matplotlib.colors import Normalize
from matplotlib.patches import Circle,Rectangle
import math






from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,MIN_AMP_PDS_PPDS,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,
					PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP_MED,DIST_GRID_PP,NUMBER_STA_PER_BIN,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,BOOTSTRAP_DEPTH_ESTIMATION,GAMMA,COLORMAP_STD,COLORMAP_VEL
				   )


print('Starting Final Plot CODE')
print('\n')

print('Calculating earth model layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)


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

filename = PP_SELEC_DIR+'SELECTED_BINNED_Ps.json'

SELECTED_BINNED_DATA_dic = json.load(open(filename))

lats = SELECTED_BINNED_DATA_dic['lat']
lons = SELECTED_BINNED_DATA_dic['lon']
RF_number = SELECTED_BINNED_DATA_dic['len_Pds']

RF_stacking_Pds = SELECTED_BINNED_DATA_dic['data_Pds']
RF_stacking_Ppds = SELECTED_BINNED_DATA_dic['data_Ppds']

RF_DEPTH_mean_1_Pds = SELECTED_BINNED_DATA_dic['mean_1_Pds']
RF_DEPTH_std_1_Pds = SELECTED_BINNED_DATA_dic['std_1_Pds']
RF_DEPTH_mean_2_Pds = SELECTED_BINNED_DATA_dic['mean_2_Pds']
RF_DEPTH_std_2_Pds = SELECTED_BINNED_DATA_dic['std_2_Pds']

RF_DEPTH_mean_1_Ppds = SELECTED_BINNED_DATA_dic['mean_1_Ppds']
RF_DEPTH_std_1_Ppds = SELECTED_BINNED_DATA_dic['std_1_Ppds']
RF_DEPTH_mean_2_Ppds = SELECTED_BINNED_DATA_dic['mean_2_Ppds']
RF_DEPTH_std_2_Ppds = SELECTED_BINNED_DATA_dic['std_2_Ppds']

RF_DEPTH_mtz_thickness_Pds = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds']
RF_DEPTH_mtz_thickness_Ppds = SELECTED_BINNED_DATA_dic['mtz_thickness_Ppds']

RF_DEPTH_true_thickness_MTZ_Pds = SELECTED_BINNED_DATA_dic['true_thickness_MTZ_Pds']
RF_DEPTH_true_thickness_MTZ_Ppds = SELECTED_BINNED_DATA_dic['true_thickness_MTZ_Ppds']

RF_DEPTH_mean_1_true_Pds = SELECTED_BINNED_DATA_dic['true_mean_1_Pds']
RF_DEPTH_std_1_true_Pds = SELECTED_BINNED_DATA_dic['true_std_1_Pds']

RF_DEPTH_mean_2_true_Pds = SELECTED_BINNED_DATA_dic['true_mean_2_Pds']
RF_DEPTH_std_2_true_Pds = SELECTED_BINNED_DATA_dic['true_std_2_Pds']

RF_DEPTH_mean_1_true_Ppds = SELECTED_BINNED_DATA_dic['true_mean_1_Ppds']
RF_DEPTH_std_1_true_Ppds = SELECTED_BINNED_DATA_dic['true_std_1_Ppds']

RF_DEPTH_mean_2_true_Ppds = SELECTED_BINNED_DATA_dic['true_mean_2_Ppds']
RF_DEPTH_std_2_true_Ppds = SELECTED_BINNED_DATA_dic['true_std_2_Ppds']

RF_delta_1_Vp_mean = SELECTED_BINNED_DATA_dic['delta_1_Vp_mean']
RF_delta_1_Vp_std = SELECTED_BINNED_DATA_dic['delta_1_Vp_std']
RF_delta_1_Vs_mean = SELECTED_BINNED_DATA_dic['delta_1_Vs_mean']
RF_delta_1_Vs_std = SELECTED_BINNED_DATA_dic['delta_1_Vs_std']

RF_delta_2_Vp_mean = SELECTED_BINNED_DATA_dic['delta_2_Vp_mean']
RF_delta_2_Vp_std = SELECTED_BINNED_DATA_dic['delta_2_Vp_std']
RF_delta_2_Vs_mean = SELECTED_BINNED_DATA_dic['delta_2_Vs_mean']
RF_delta_2_Vs_std = SELECTED_BINNED_DATA_dic['delta_2_Vs_std']

RESULTS_FOLDER = PP_FIGURE+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)


###################################################################################################################

#Color Maps

colormap = plt.get_cmap(COLORMAP_VEL)


colormap_std = plt.get_cmap(COLORMAP_STD)


#############################################################################################################################################################################################

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,20))

majorLocatorY = MultipleLocator(20)
minorLocatorY = MultipleLocator(5)

	
idx_Pds = np.argsort(RF_DEPTH_true_thickness_MTZ_Pds)

norm_410 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)
colors_410_Pds = colormap(norm_410(np.array(RF_DEPTH_true_thickness_MTZ_Pds,dtype='float64')))

for j,i in enumerate(idx_Pds):
	if math.isnan(RF_DEPTH_true_thickness_MTZ_Pds[i]) == False:
		ax.errorbar(RF_DEPTH_mean_2_true_Pds[i],RF_DEPTH_mean_1_true_Pds[i], yerr=RF_DEPTH_std_1_true_Pds[i],xerr=RF_DEPTH_std_2_true_Pds[i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
		ax.scatter(RF_DEPTH_mean_2_true_Pds[i],RF_DEPTH_mean_1_true_Pds[i],color=colors_410_Pds[i],zorder=10)

ax.yaxis.set_major_locator(majorLocatorY)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
ax.yaxis.set_ticks_position('both')

ax.set_ylim(360,460)
ax.set_xlim(610,710)
ax.set_ylabel('True Depth d410 (km)')
ax.set_xlabel('True Depth d660 (km)')

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8,label='MTZ True Thickness Pds')

plt.show()
fig.savefig(RESULTS_FOLDER+'TRUE_DEPTH_PLOT.'+EXT_FIG,dpi=DPI_FIG)


#############################################################################################################################################################################################

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20,10),sharey=True)

ax = axes[0]
ax1 = axes[1]


majorLocatorY = MultipleLocator(20)
minorLocatorY = MultipleLocator(5)

	
idx_Pds = np.argsort(RF_DEPTH_mtz_thickness_Pds)
idx_Ppds = np.argsort(RF_DEPTH_mtz_thickness_Ppds)

norm_410 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)
colors_410_Pds = colormap(norm_410(np.array(RF_DEPTH_mtz_thickness_Pds,dtype='float64')))
colors_410_Ppds = colormap(norm_410(np.array(RF_DEPTH_mtz_thickness_Ppds,dtype='float64')))

for j,i in enumerate(idx_Pds):
	if math.isnan(RF_DEPTH_mtz_thickness_Pds[i]) == False:
		ax.errorbar(RF_DEPTH_mean_2_Pds[i],RF_DEPTH_mean_1_Pds[i], yerr=RF_DEPTH_std_1_Pds[i],xerr=RF_DEPTH_std_2_Pds[i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
		ax.scatter(RF_DEPTH_mean_2_Pds[i],RF_DEPTH_mean_1_Pds[i],color=colors_410_Pds[i],zorder=10)

ax.yaxis.set_major_locator(majorLocatorY)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
ax.yaxis.set_ticks_position('both')

ax.set_ylim(360,460)
ax.set_xlim(610,710)
ax.set_ylabel('Depth d410 (km)')
ax.set_xlabel('Depth d660 (km)')
ax.set_title('Pds Phases')


for j,i in enumerate(idx_Ppds):
	if math.isnan(RF_DEPTH_mtz_thickness_Ppds[i]) == False:
		ax1.errorbar(RF_DEPTH_mean_2_Ppds[i],RF_DEPTH_mean_1_Ppds[i], yerr=RF_DEPTH_std_1_Ppds[i],xerr=RF_DEPTH_std_2_Ppds[i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
		ax1.scatter(RF_DEPTH_mean_2_Ppds[i],RF_DEPTH_mean_1_Ppds[i],color=colors_410_Ppds[i],zorder=10)

ax1.yaxis.set_major_locator(majorLocatorY)
ax1.yaxis.set_minor_locator(minorLocatorY)
ax1.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
ax1.yaxis.set_ticks_position('both')

ax1.tick_params(labelright=True,labelleft=False)
ax1.yaxis.set_label_position("right")

ax1.set_ylim(360,460)
ax1.set_xlim(610,710)
ax1.set_ylabel('Depth d410 (km)')
ax1.set_xlabel('Depth d660 (km)')
ax1.set_title('Ppds Phases')

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8,label='MTZ Thickness Pds')
fig.colorbar(sm_410,ax=ax1,orientation='horizontal',shrink=0.8,label='MTZ Thickness Ppds')

plt.show()
fig.savefig(RESULTS_FOLDER+'APPARENT_DEPTH_PDS_PPDS_PLOT.'+EXT_FIG,dpi=DPI_FIG)


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
colors_410 = colormap(norm_410(np.array(RF_DEPTH_mean_1_Pds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_1_Pds[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)
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

norm_660 = mpl.colors.Normalize(vmin=560,vmax=760,clip=True)
colors_660 = colormap(norm_660(np.array(RF_DEPTH_mean_2_Pds,dtype='float64')))


for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_2_Pds[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
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

print('Plotting Figure: Apparent Depth of 410 km and 660 km (Ppds phase)')

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
colors_410 = colormap(norm_410(np.array(RF_DEPTH_mean_1_Ppds,dtype='float64')))
for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_1_Ppds[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('410 km Ppds', y=1.08)


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

norm_660 = mpl.colors.Normalize(vmin=560,vmax=760,clip=True)
colors_660 = colormap(norm_660(np.array(RF_DEPTH_mean_2_Ppds,dtype='float64')))


for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_2_Ppds[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
	else:
		pass

ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax2.set_title('660 km Ppds', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

sm_660 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_660)
sm_660._A = []
fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8)

fig.savefig(RESULTS_FOLDER+'Apparent_depth_Ppds.'+EXT_FIG,dpi=DPI_FIG)

#############################################################################################################################################################################################

print('Plotting Figure: True Depth of 410 km and 660 km (Pds phases)')

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
colors_410 = colormap(norm_410(np.array(RF_DEPTH_mean_1_true_Pds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_1_true_Pds[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)
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

norm_660 = mpl.colors.Normalize(vmin=560,vmax=760,clip=True)
colors_660 = colormap(norm_660(np.array(RF_DEPTH_mean_2_true_Pds,dtype='float64')))


for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_2_true_Pds[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
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

fig.savefig(RESULTS_FOLDER+'TRUE_DEPTH_Pds.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################
print('Plotting Figure: True Depth of 410 km and 660 km (Ppds phases)')

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
colors_410 = colormap(norm_410(np.array(RF_DEPTH_mean_1_true_Ppds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_1_true_Ppds[i]) == False:

		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)

	else:
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('410 km Ppds', y=1.08)


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

norm_660 = mpl.colors.Normalize(vmin=560,vmax=760,clip=True)
colors_660 = colormap(norm_660(np.array(RF_DEPTH_mean_2_true_Ppds,dtype='float64')))


for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mean_1_true_Ppds[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
	else:
		pass

ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax2.set_title('660 km Ppds', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

sm_660 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_660)
sm_660._A = []
fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8)

fig.savefig(RESULTS_FOLDER+'TRUE_DEPTH_Ppds.'+EXT_FIG,dpi=DPI_FIG)

###################################################################################################################
print('Plotting Figure: Uncertainty (1 sigma) of 410 km and 660 km (Pds phases)')

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

norm_410 = mpl.colors.Normalize(vmin=0,vmax=50,clip=True)
colors_410 = colormap_std(norm_410(np.array(RF_DEPTH_std_1_Pds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_std_1_Pds[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)
	else:
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title(r'Uncertainty (1$\sigma$) of  410 km Pds', y=1.08)


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

norm_660 = mpl.colors.Normalize(vmin=0,vmax=50,clip=True)
colors_660 = colormap_std(norm_660(np.array(RF_DEPTH_std_2_Pds,dtype='float64')))


for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_std_2_Pds[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
	else: 
		pass

ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax2.set_title(r'Uncertainty (1$\sigma$) of  660 km Pds', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap_std,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

sm_660 = plt.cm.ScalarMappable(cmap=colormap_std,norm=norm_660)
sm_660._A = []
fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8)

fig.savefig(RESULTS_FOLDER+'Uncertainty_DEPTH_Pds.'+EXT_FIG,dpi=DPI_FIG)

#######################################################################################################################################
print('Plotting Figure: Uncertainty (1 sigma) of 410 km and 660 km (Ppds phases)')

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

norm_410 = mpl.colors.Normalize(vmin=0,vmax=50,clip=True)
colors_410 = colormap_std(norm_410(np.array(RF_DEPTH_std_1_Ppds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_std_1_Ppds[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)
	else:
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title(r'Uncertainty (1$\sigma$) of  410 km Pds', y=1.08)


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

norm_660 = mpl.colors.Normalize(vmin=0,vmax=50,clip=True)
colors_660 = colormap_std(norm_660(np.array(RF_DEPTH_std_2_Ppds,dtype='float64')))


for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_std_2_Ppds[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
	else:
		pass
ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax2.set_title(r'Uncertainty (1$\sigma$) of  660 km Pds', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap_std,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

sm_660 = plt.cm.ScalarMappable(cmap=colormap_std,norm=norm_660)
sm_660._A = []
fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8)

fig.savefig(RESULTS_FOLDER+'Uncertainty_DEPTH_Ppds.'+EXT_FIG,dpi=DPI_FIG)

#######################################################################################################################################

print('Plotting Figure: Delta Vp (410 km/660 km) ')


fig, axes = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(20,10),sharey=True)

ax = axes[0]
ax2 = axes[1]

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

norm_410 = mpl.colors.Normalize(vmin=-1,vmax=1,clip=True)
colors_410 = colormap(norm_410(np.array(RF_delta_1_Vp_mean,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_delta_1_Vp_mean[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)
	else: 
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('Delta Vp - 410 km', y=1.08)


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

norm_660 = mpl.colors.Normalize(vmin=-1,vmax=1,clip=True)
colors_660 = colormap(norm_660(np.array(RF_delta_2_Vp_mean,dtype='float64')))


for i,j in enumerate(lons):
	if math.isnan(RF_delta_2_Vp_mean[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
	else:
		pass

ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax2.set_title('Delta Vp - 660 km', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

sm_660 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_660)
sm_660._A = []
fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8)

fig.savefig(RESULTS_FOLDER+'DELTA_VP.'+EXT_FIG,dpi=DPI_FIG)

#######################################################################################################################################

print('Plotting Figure: Delta Vs (410 km/660 km) ')


fig, axes = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(20,10),sharey=True)

ax = axes[0]
ax2 = axes[1]

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

norm_410 = mpl.colors.Normalize(vmin=-1,vmax=1,clip=True)
colors_410 = colormap(norm_410(np.array(RF_delta_1_Vs_mean,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_delta_1_Vs_mean[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)
	else: 
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('Delta Vs - 410 km', y=1.08)


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

norm_660 = mpl.colors.Normalize(vmin=-1,vmax=1,clip=True)
colors_660 = colormap(norm_660(np.array(RF_delta_2_Vs_mean,dtype='float64')))


for i,j in enumerate(lons):
	if math.isnan(RF_delta_2_Vs_mean[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
	else:
		pass

ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax2.set_title('Delta Vs - 660 km', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

sm_660 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_660)
sm_660._A = []
fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8)

fig.savefig(RESULTS_FOLDER+'DELTA_VS.'+EXT_FIG,dpi=DPI_FIG)


#######################################################################################################################################

print('Plotting Figure: Thickness of the Mantle Transition Zone (Pds and Ppds Phases)')

fig, axes = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(20,10),sharey=True)

ax = axes[0]
ax2 = axes[1]

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
ax.gridlines(draw_labels=True)

norm_410 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)
colors_410 = colormap(norm_410(np.array(RF_DEPTH_mtz_thickness_Pds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mtz_thickness_Pds[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)
	else:
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('Thickness of MTZ (Pds)', y=1.08)


#Ppds phase

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

norm_660 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)
colors_660 = colormap(norm_660(np.array(RF_DEPTH_mtz_thickness_Ppds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mtz_thickness_Ppds[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
	else:
		pass

ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax2.set_title('Thickness of MTZ (Ppds)', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

sm_660 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_660)
sm_660._A = []
fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8)

fig.savefig(RESULTS_FOLDER+'THICKNESS_MTZ.'+EXT_FIG,dpi=DPI_FIG)

##############################################################################################

print('Plotting Figure: True Thickness of the Mantle Transition Zone')

fig, axes = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(20,10),sharey=True)

ax = axes[0]
ax2 = axes[1]

#Pds Phase


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

norm_410 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)
colors_410 = colormap(norm_410(np.array(RF_DEPTH_true_thickness_MTZ_Pds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_true_thickness_MTZ_Pds[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax.add_patch(retangulo_410)
	else:
		pass


ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax.set_title('True Thickness of MTZ (Pds)', y=1.08)

#Ppds Phase

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

norm_660 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)
colors_660 = colormap(norm_660(np.array(RF_DEPTH_true_thickness_MTZ_Ppds,dtype='float64')))

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_true_thickness_MTZ_Ppds[i]) == False:
		retangulo_660 = Rectangle(xy=(lons[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2, lats[i]-(DIST_GRID_PP_MED/(GRID_PP_MULT/2))/2),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_660[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
		ax2.add_patch(retangulo_660)
	else:
		pass

ax2.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

ax2.set_title('True Thickness of MTZ (Ppds)', y=1.08)

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

sm_660 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_660)
sm_660._A = []
fig.colorbar(sm_660,ax=ax2,orientation='horizontal',shrink=0.8)

fig.savefig(RESULTS_FOLDER+'TRUE_THICKNESS_MTZ.'+EXT_FIG,dpi=DPI_FIG)

########################################################################################################################################################################

print('Ending Final Plot CODE')