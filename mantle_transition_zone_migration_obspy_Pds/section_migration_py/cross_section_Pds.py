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
import scipy.io
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator, FixedLocator
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import json
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import interpolate
from shapely.geometry import Polygon, MultiPoint, Point, LineString
import shapefile
from matplotlib.colors import Normalize
from matplotlib.patches import Circle,Rectangle
import math
import verde as vd






from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,TECTO_SHP,COLORMAP_VEL,COLORMAP_STD,OUTPUT_DIR,
					EXT_FIG,DPI_FIG,DIST_GRID_PP,NUMBER_STA_PER_BIN,NUMBER_PP_PER_BIN,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,CROSS_SECTION_AXIS,DEPTH_TARGET
				   )

print('Starting Cross section CODE')
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

#Lat/Lon lists:

lats = SELECTED_BINNED_DATA_dic['lat']
lons = SELECTED_BINNED_DATA_dic['lon']

#Receiver Functions:

RF_number = SELECTED_BINNED_DATA_dic['len_Pds']

RF_stacking_Pds = SELECTED_BINNED_DATA_dic['data_Pds']

#Estimates 350 LVZ:

RF_DEPTH_mean_LVZ_Pds = SELECTED_BINNED_DATA_dic['mean_LVZ_Pds']
RF_DEPTH_std_LVZ_Pds = SELECTED_BINNED_DATA_dic['std_LVZ_Pds']


#Estimates P410s:

RF_DEPTH_mean_1_Pds = SELECTED_BINNED_DATA_dic['mean_1_Pds']
RF_DEPTH_std_1_Pds = SELECTED_BINNED_DATA_dic['std_1_Pds']

#Estimates P520s:

RF_DEPTH_mean_520_Pds = SELECTED_BINNED_DATA_dic['mean_520_Pds']
RF_DEPTH_std_520_Pds = SELECTED_BINNED_DATA_dic['std_520_Pds']

#Estimates P660s:

RF_DEPTH_mean_2_Pds = SELECTED_BINNED_DATA_dic['mean_2_Pds']
RF_DEPTH_std_2_Pds = SELECTED_BINNED_DATA_dic['std_2_Pds']

#Estimates 700 LVZ:

RF_DEPTH_mean_LVZ_700_Pds = SELECTED_BINNED_DATA_dic['mean_LVZ_700_Pds']
RF_DEPTH_std_LVZ_700_Pds = SELECTED_BINNED_DATA_dic['std_LVZ_700_Pds']

#Estimates MTZ Pds:

RF_DEPTH_mtz_thickness_Pds = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds']
RF_DEPTH_mtz_thickness_Pds_std = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std']

#############################################################################################################################3


# create the grid coordinates

spacing = GRID_PP_MULT
region = (LLCRNRLON_SMALL, URCRNRLON_SMALL, LLCRNRLAT_SMALL, URCRNRLAT_SMALL)

rows, cols = vd.grid_coordinates(region=region, spacing=spacing)

PP_FIGURE = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Figures'+'/'
RESULTS_FOLDER = PP_FIGURE+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)


print('Allocating cross-section data')
print('\n')

if CROSS_SECTION_AXIS == 'x':

	#Profile lat/lon

	AB_lon = [[]]*len(rows[:,0])
	AB_lat = [[]]*len(rows[:,0])

	#Receiver Functions:

	RF_data_profile_Pds = [[]]*len(rows[:,0])

	#Estimates 350 LVZ:

	RF_DEPTH_mean_LVZ_profile_Pds = [[]]*len(rows[:,0])
	RF_DEPTH_std_LVZ_profile_Pds = [[]]*len(rows[:,0])

	#Estimates P410s:

	RF_DEPTH_mean_1_profile_Pds = [[]]*len(rows[:,0])
	RF_DEPTH_std_1_profile_Pds = [[]]*len(rows[:,0])

	#Estimates P520s:

	RF_DEPTH_mean_520_profile_Pds = [[]]*len(rows[:,0])
	RF_DEPTH_std_520_profile_Pds = [[]]*len(rows[:,0])

	#Estimates P660s:

	RF_DEPTH_mean_2_profile_Pds = [[]]*len(rows[:,0])
	RF_DEPTH_std_2_profile_Pds = [[]]*len(rows[:,0])

	#Estimates 700 LVZ:

	RF_DEPTH_mean_LVZ_700_profile_Pds = [[]]*len(rows[:,0])
	RF_DEPTH_std_LVZ_700_profile_Pds = [[]]*len(rows[:,0])

	#Estimates MTZ Pds:

	RF_DEPTH_mtz_thickness_profile_Pds = [[]]*len(rows[:,0]) 
	RF_DEPTH_mtz_thickness_profile_Pds_std = [[]]*len(rows[:,0]) 

	for i,j in enumerate(rows[:,0]):

		lat_lon = [(lons[k],lats[k]) for k,l in enumerate(lats)]
		grid_column = [(float("%.4f" % round(rows[i,:][k],4)),float("%.4f" % round(cols[i,:][k],4))) for k,l in enumerate(rows[i,:])]

		#Profile lat/lon

		AB_lon[i] = [rows[i,:][k] for k,l in enumerate(rows[i,:])]
		AB_lat[i] = [cols[i,:][k] for k,l in enumerate(rows[i,:])]


		#Receiver Functions:

		RF_data_profile_Pds[i] = [RF_stacking_Pds[lat_lon.index(l)] if l in lat_lon else np.zeros_like(RF_stacking_Pds[k]) for k,l in enumerate(grid_column)]

		#Estimates 350 LVZ:

		RF_DEPTH_mean_LVZ_profile_Pds[i] = [RF_DEPTH_mean_LVZ_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_LVZ_profile_Pds[i] = [RF_DEPTH_std_LVZ_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]


		#Estimates P410s:

		RF_DEPTH_mean_1_profile_Pds[i] = [RF_DEPTH_mean_1_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_1_profile_Pds[i] = [RF_DEPTH_std_1_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates P520s:

		RF_DEPTH_mean_520_profile_Pds[i] = [RF_DEPTH_mean_520_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_520_profile_Pds[i] = [RF_DEPTH_std_520_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates P660s:

		RF_DEPTH_mean_2_profile_Pds[i] = [RF_DEPTH_mean_2_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_2_profile_Pds[i] = [RF_DEPTH_std_2_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates 700 LVZ:

		RF_DEPTH_mean_LVZ_700_profile_Pds[i] = [RF_DEPTH_mean_LVZ_700_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_LVZ_700_profile_Pds[i] = [RF_DEPTH_std_LVZ_700_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates MTZ Pds:

		RF_DEPTH_mtz_thickness_profile_Pds[i] = [RF_DEPTH_mtz_thickness_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_mtz_thickness_profile_Pds_std[i] = [RF_DEPTH_mtz_thickness_Pds_std[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

else:

	#Profile lat/lon

	AB_lon = [[]]*len(rows[0,:])
	AB_lat = [[]]*len(rows[0,:])

	#Receiver Functions:

	RF_data_profile_Pds = [[]]*len(rows[0,:])

	#Estimates 350 LVZ:

	RF_DEPTH_mean_LVZ_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_LVZ_profile_Pds = [[]]*len(rows[0,:])


	#Estimates P410s:

	RF_DEPTH_mean_1_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_1_profile_Pds = [[]]*len(rows[0,:])


	#Estimates P520s:

	RF_DEPTH_mean_520_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_520_profile_Pds = [[]]*len(rows[0,:])

	#Estimates P660s:

	RF_DEPTH_mean_2_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_2_profile_Pds = [[]]*len(rows[0,:])

	#Estimates 700 LVZ:

	RF_DEPTH_mean_LVZ_700_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_LVZ_700_profile_Pds = [[]]*len(rows[0,:])

	#Estimates MTZ Pds:

	RF_DEPTH_mtz_thickness_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_mtz_thickness_profile_Pds_std = [[]]*len(rows[0,:])

	for i,j in enumerate(rows[0,:]):
		
		lat_lon = [(lons[k],lats[k]) for k,l in enumerate(lats)]
		grid_column = [(float("%.4f" % round(rows[:,i][k],4)),float("%.4f" % round(cols[:,i][k],4))) for k,l in enumerate(rows[:,i])]

		#Profile lat/lon

		AB_lon[i] = [rows[:,i][k] for k,l in enumerate(rows[:,i])]
		AB_lat[i] = [cols[:,i][k] for k,l in enumerate(rows[:,i])]

		#Receiver Functions:

		RF_data_profile_Pds[i] = [RF_stacking_Pds[lat_lon.index(l)] if l in lat_lon else np.zeros_like(RF_stacking_Pds[k]) for k,l in enumerate(grid_column)]

		#Estimates 350 LVZ:

		RF_DEPTH_mean_LVZ_profile_Pds[i] = [RF_DEPTH_mean_LVZ_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_LVZ_profile_Pds[i] = [RF_DEPTH_std_LVZ_Pds[lat_lon.index(l)]  if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates P410s:

		RF_DEPTH_mean_1_profile_Pds[i] = [RF_DEPTH_mean_1_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_1_profile_Pds[i] = [RF_DEPTH_std_1_Pds[lat_lon.index(l)]  if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates P520s:

		RF_DEPTH_mean_520_profile_Pds[i] = [RF_DEPTH_mean_520_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_520_profile_Pds[i] = [RF_DEPTH_std_520_Pds[lat_lon.index(l)]  if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates P660s:

		RF_DEPTH_mean_2_profile_Pds[i] = [RF_DEPTH_mean_2_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_2_profile_Pds[i] = [RF_DEPTH_std_2_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates 700 LVZ:

		RF_DEPTH_mean_LVZ_700_profile_Pds[i] = [RF_DEPTH_mean_LVZ_700_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_LVZ_700_profile_Pds[i] = [RF_DEPTH_std_LVZ_700_Pds[lat_lon.index(l)]  if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates MTZ Pds:

		RF_DEPTH_mtz_thickness_profile_Pds[i] = [RF_DEPTH_mtz_thickness_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_mtz_thickness_profile_Pds_std[i] = [RF_DEPTH_mtz_thickness_Pds_std[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

	
print('Plotting cross-sections according to the '+CROSS_SECTION_AXIS+' direction')

for i,j in enumerate(RF_data_profile_Pds):

	fig = plt.figure(figsize=(30, 10))

	fig.suptitle('Cross section for Pds')

	gs = gridspec.GridSpec(4, 3)
	gs.update(wspace=0.2, hspace=0.5)

	MTZ_thickness = fig.add_subplot(gs[1,1:])
	diff_MTZ_thickness = fig.add_subplot(gs[0,1:],sharex=MTZ_thickness)

	pefil_pds = fig.add_subplot(gs[2:4,1:],sharex=MTZ_thickness)

	#_____________________________________________

	map_MTZ_thickness =  fig.add_subplot(gs[0:2,0],projection=ccrs.PlateCarree())

	#_____________________________________________

	apparent_410 = fig.add_subplot(gs[2,0])
	apparent_660 = fig.add_subplot(gs[3,0],sharex=apparent_410)


	#######################################################################

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

	map_MTZ_thickness.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
	map_MTZ_thickness.yaxis.set_ticks_position('both')
	map_MTZ_thickness.xaxis.set_ticks_position('both')

	map_MTZ_thickness.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE,3), crs=ccrs.PlateCarree())
	map_MTZ_thickness.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE,3), crs=ccrs.PlateCarree())
	map_MTZ_thickness.tick_params(labelbottom=False,labeltop=True,labelleft=True,labelright=True, labelsize=15)

	map_MTZ_thickness.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

	reader_1_SHP = Reader(BOUNDARY_1_SHP)
	shape_1_SHP = list(reader_1_SHP.geometries())
	plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
	map_MTZ_thickness.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

	reader_2_SHP = Reader(BOUNDARY_2_SHP)
	shape_2_SHP = list(reader_2_SHP.geometries())
	plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
	map_MTZ_thickness.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

	bounds = np.arange(200, 300+colormap_segmentation, colormap_segmentation)
	norm_map_MTZ_thickness = Normalize(vmin=200,vmax=300)

	for t,y in enumerate(lons):
		if math.isnan(RF_DEPTH_mtz_thickness_Pds[t]) == False:
			if RF_DEPTH_mtz_thickness_Pds_std[t] < 10:
				circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[t], lats[t]),color=tmap(norm_map_MTZ_thickness(RF_DEPTH_mtz_thickness_Pds[t])),ec='k',linewidth=0.5,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
				map_MTZ_thickness.add_patch(circulo_410)
			else:
				circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[t], lats[t]),color=tmap(norm_map_MTZ_thickness(RF_DEPTH_mtz_thickness_Pds[t])),  ec='k',linewidth=0.5,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
				map_MTZ_thickness.add_patch(circulo_410)
		else:
			pass

	
	for x,c in enumerate(AB_lon[i]):
		circulo_410_profile = Circle(radius=DIST_GRID_PP,xy=(AB_lon[i][x], AB_lat[i][x]),fc='None',ec='k',transform=ccrs.Geodetic(),zorder=10)
		map_MTZ_thickness.add_patch(circulo_410_profile)

	map_MTZ_thickness.set_title('MTZ Thickness', y=1.1, fontsize=20)

	sm_map_MTZ_thickness = plt.cm.ScalarMappable(cmap=tmap,norm=norm_map_MTZ_thickness)
	sm_map_MTZ_thickness._A = []
	cbar = fig.colorbar(sm_map_MTZ_thickness,ax=map_MTZ_thickness,orientation='vertical',shrink=0.9,pad=0.1)

	cbar.set_ticks(np.arange(200, 300+INTER_DEPTH, INTER_DEPTH*2))
	cbar.set_ticklabels(np.arange(200, 300+INTER_DEPTH, INTER_DEPTH*2))
	cbar.ax.tick_params(labelsize=15)


	#### Profile  ####

	factor_Pds = 300

	majorLocatorY = MultipleLocator(100)
	minorLocatorY = MultipleLocator(10)

	
	for _i, _j in enumerate(RF_data_profile_Pds[i]):
		RF_data_factor_Pds = [_i/factor_Pds+l for k, l in enumerate(_j)]
		pefil_pds.plot(RF_data_factor_Pds,camadas_terra_10_km,'k',linewidth=1)

		pefil_pds.plot(_i/factor_Pds,RF_DEPTH_mean_LVZ_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')
		pefil_pds.plot(_i/factor_Pds,RF_DEPTH_mean_1_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')
		pefil_pds.plot(_i/factor_Pds,RF_DEPTH_mean_520_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')
		pefil_pds.plot(_i/factor_Pds,RF_DEPTH_mean_2_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')
		pefil_pds.plot(_i/factor_Pds,RF_DEPTH_mean_LVZ_700_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')

		pefil_pds.errorbar(_i/factor_Pds,RF_DEPTH_mean_LVZ_profile_Pds[i][_i], yerr=RF_DEPTH_std_LVZ_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_pds.errorbar(_i/factor_Pds,RF_DEPTH_mean_1_profile_Pds[i][_i], yerr=RF_DEPTH_std_1_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_pds.errorbar(_i/factor_Pds,RF_DEPTH_mean_520_profile_Pds[i][_i], yerr=RF_DEPTH_std_520_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_pds.errorbar(_i/factor_Pds,RF_DEPTH_mean_2_profile_Pds[i][_i], yerr=RF_DEPTH_std_2_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_pds.errorbar(_i/factor_Pds,RF_DEPTH_mean_LVZ_700_profile_Pds[i][_i], yerr=RF_DEPTH_std_LVZ_700_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)

		pefil_pds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)>=_i/factor_Pds, facecolor='black',alpha=0.3, interpolate=True)
		pefil_pds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)<=_i/factor_Pds, facecolor='gray',alpha=0.3, interpolate=True)
		
		pefil_pds.set_ylim(MAX_DEPTH,MIN_DEPTH)
		pefil_pds.yaxis.set_ticks_position('both')
		pefil_pds.yaxis.set_major_locator(majorLocatorY)
		pefil_pds.yaxis.set_minor_locator(minorLocatorY)

		pefil_pds.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		pefil_pds.set_ylabel('Depth (km)', fontsize=18)
		pefil_pds.yaxis.set_label_position("right")
		pefil_pds.tick_params(labelright=True,labelsize=15)
		pefil_pds.set_ylim(800,300)

		if CROSS_SECTION_AXIS == 'x':
			pefil_pds.set_xticks(np.linspace(pefil_pds.get_xlim()[0],pefil_pds.get_xlim()[1],10))
			#pefil_pds.set_xticklabels(np.linspace(LLCRNRLON_SMALL,URCRNRLON_SMALL,10))
			pefil_pds.set_xlabel('Longitude ($^\circ$)',labelpad=30,fontsize=20)
			pefil_pds.set_title('Latitude= '+str(AB_lat[i][0])+'$^\circ$', fontsize=20)

		else:
			pefil_pds.set_xticks(np.linspace(pefil_pds.get_xlim()[0],pefil_pds.get_xlim()[1],10))
			#pefil_pds.set_xticklabels(np.linspace(LLCRNRLAT_SMALL,URCRNRLAT_SMALL,10))
			pefil_pds.set_xlabel('Latitude ($^\circ$)',labelpad=30,fontsize=20)
			pefil_pds.set_title('Longitude = '+str(AB_lon[i][0])+'$^\circ$', fontsize=20)


	#### Figure Apparent  410 km Pds  ####

	for _i, _j in enumerate(RF_data_profile_Pds[i]):
		apparent_410.plot(_i,RF_DEPTH_mean_1_profile_Pds[i][_i]-410,marker='o',markerfacecolor='none',markeredgecolor='dimgray')

		apparent_410.errorbar(_i,RF_DEPTH_mean_1_profile_Pds[i][_i]-410, yerr=RF_DEPTH_std_1_profile_Pds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)

		apparent_410.set_title('diff 410 km', fontsize=20)
		apparent_410.yaxis.set_ticks_position('both')
		apparent_410.yaxis.set_major_locator(MultipleLocator(25))
		apparent_410.yaxis.set_minor_locator(MultipleLocator(10))
		apparent_410.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_410.tick_params(labelleft=True,labelright=True,labelsize=15)
		apparent_410.yaxis.set_label_position("right")
		apparent_410.set_ylim(50,-50)
		apparent_410.set_xticks([])



	#### Figure Apparent  660 km Pds  ####


	for _i, _j in enumerate(RF_DEPTH_mean_2_profile_Pds[i]):
		apparent_660.plot(_i,RF_DEPTH_mean_2_profile_Pds[i][_i]-660,marker='o',markerfacecolor='none',markeredgecolor='dimgray')

		apparent_660.errorbar(_i,RF_DEPTH_mean_2_profile_Pds[i][_i]-660, yerr=RF_DEPTH_std_2_profile_Pds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)
		
		apparent_660.set_ylim(50,-50)
		apparent_660.set_title('diff 660 km ', fontsize=20)
		apparent_660.yaxis.set_ticks_position('both')
		apparent_660.yaxis.set_major_locator(MultipleLocator(25))
		apparent_660.yaxis.set_minor_locator(MultipleLocator(10))
		apparent_660.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_660.tick_params(labelleft=True,labelright=True,labelsize=15)
		apparent_660.yaxis.set_label_position("right")
		apparent_660.set_xticks([])


	#### Figure MTZ Apparent thickness  ####

	for _i, _j in enumerate(RF_DEPTH_mtz_thickness_profile_Pds[i]):
		MTZ_thickness.plot(_i/factor_Pds,RF_DEPTH_mtz_thickness_profile_Pds[i][_i]-250,marker='o',markerfacecolor='none',markeredgecolor='gray')

	for _i, _j in enumerate(RF_DEPTH_mtz_thickness_profile_Pds[i]):
		MTZ_thickness.errorbar(_i/factor_Pds,RF_DEPTH_mtz_thickness_profile_Pds[i][_i]-250, yerr=RF_DEPTH_mtz_thickness_profile_Pds_std[i][_i], ecolor='gray',elinewidth=1,capsize=2,capthick=1)
	
	MTZ_thickness.set_ylim(50,-50)
	MTZ_thickness.yaxis.set_label_position("right")
	MTZ_thickness.set_title('diff MTZ Thickness', fontsize=20)
	MTZ_thickness.yaxis.set_ticks_position('both')
	MTZ_thickness.yaxis.set_major_locator(MultipleLocator(25))
	MTZ_thickness.yaxis.set_minor_locator(MultipleLocator(10))
	MTZ_thickness.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
	MTZ_thickness.tick_params(labelleft=True,labelright=True,labelbottom=False,labelsize=15)


	#### Figure diff MTZ Apparent thickness  ####

	for _i, _j in enumerate(RF_DEPTH_mean_LVZ_profile_Pds[i]):
		diff_MTZ_thickness.plot(_i/factor_Pds,RF_DEPTH_mean_LVZ_profile_Pds[i][_i]-350,marker='o',markerfacecolor='none',markeredgecolor='gray')

	for _i, _j in enumerate(RF_DEPTH_mean_LVZ_profile_Pds[i]):
		diff_MTZ_thickness.errorbar(_i/factor_Pds,RF_DEPTH_mean_LVZ_profile_Pds[i][_i]-350, yerr=RF_DEPTH_std_LVZ_profile_Pds[i][_i], ecolor='gray',elinewidth=1,capsize=2,capthick=1)
	
	diff_MTZ_thickness.set_ylim(50,-50)
	diff_MTZ_thickness.yaxis.set_label_position("right")
	diff_MTZ_thickness.set_title('diff LVZ', fontsize=20)
	diff_MTZ_thickness.yaxis.set_ticks_position('both')
	diff_MTZ_thickness.yaxis.set_major_locator(MultipleLocator(25))
	diff_MTZ_thickness.yaxis.set_minor_locator(MultipleLocator(10))
	diff_MTZ_thickness.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
	diff_MTZ_thickness.tick_params(labelleft=True,labelright=True,labelbottom=False,labelsize=15)


	fig.savefig(RESULTS_FOLDER+'SELECTED_BINNED_DATA_'+CROSS_SECTION_AXIS+'_CROSS_SECTION_Pds_PROFILE_'+str(i+1)+'.'+EXT_FIG,dpi=DPI_FIG)
print('Ending the Cross-section CODE')