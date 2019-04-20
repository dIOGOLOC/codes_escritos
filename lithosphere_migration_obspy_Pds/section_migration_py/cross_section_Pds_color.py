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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes







from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_SPACE,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,TECTO_SHP,COLORMAP_VEL,COLORMAP_STD,OUTPUT_DIR,
					EXT_FIG,DPI_FIG,DIST_GRID_PP,NUMBER_STA_PER_BIN,NUMBER_PP_PER_BIN,VMIN,VMAX,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,CROSS_SECTION_AXIS,DEPTH_MOHO,DEPTH_LAB
				   )


print('Starting Cross section CODE')
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

#Lat/Lon lists:

lats = SELECTED_BINNED_DATA_dic['lat']
lons = SELECTED_BINNED_DATA_dic['lon']

#Receiver Functions:

RF_number = SELECTED_BINNED_DATA_dic['len_Pds']

RF_stacking_Pds = SELECTED_BINNED_DATA_dic['data_Pds']

#Estimates MOHO:

RF_DEPTH_mean_MOHO_Pds = SELECTED_BINNED_DATA_dic['mean_MOHO_Pds']
RF_DEPTH_std_MOHO_Pds = SELECTED_BINNED_DATA_dic['std_MOHO_Pds']

#Estimates LAB:

RF_DEPTH_mean_LAB_Pds = SELECTED_BINNED_DATA_dic['mean_LAB_Pds']
RF_DEPTH_std_LAB_Pds = SELECTED_BINNED_DATA_dic['std_LAB_Pds']

#############################################################################################################################3

region = (LLCRNRLON_SMALL,URCRNRLON_SMALL, LLCRNRLAT_SMALL, URCRNRLAT_SMALL)

grdx, grdy = vd.grid_coordinates(region=region,spacing=(GRID_PP_SPACE,GRID_PP_SPACE))

rows = grdx
cols = grdy

PP_FIGURE = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_MOHO_'+str(DEPTH_MOHO)+'_DEPTH_LAB_'+str(DEPTH_LAB)+'/'+'Figures'+'/'

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

	#Estimates MOHO:

	RF_DEPTH_mean_MOHO_profile_Pds = [[]]*len(rows[:,0])
	RF_DEPTH_std_MOHO_profile_Pds = [[]]*len(rows[:,0])

	#Estimates LAB:

	RF_DEPTH_mean_LAB_profile_Pds = [[]]*len(rows[:,0])
	RF_DEPTH_std_LAB_profile_Pds = [[]]*len(rows[:,0])

	for i,j in enumerate(rows[:,0]):

		lat_lon = [(lons[k],lats[k]) for k,l in enumerate(lats)]
		grid_column = [(float("%.4f" % round(rows[i,:][k],4)),float("%.4f" % round(cols[i,:][k],4))) for k,l in enumerate(rows[i,:])]

		#Profile lat/lon

		AB_lon[i] = [rows[i,:][k] for k,l in enumerate(rows[i,:])]
		AB_lat[i] = [cols[i,:][k] for k,l in enumerate(rows[i,:])]


		#Receiver Functions:

		RF_data_profile_Pds[i] = [RF_stacking_Pds[lat_lon.index(l)] if l in lat_lon else np.zeros_like(RF_stacking_Pds[k]) for k,l in enumerate(grid_column)]

		#Estimates MOHO:

		RF_DEPTH_mean_MOHO_profile_Pds[i] = [RF_DEPTH_mean_MOHO_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_MOHO_profile_Pds[i] = [RF_DEPTH_std_MOHO_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates LAB:

		RF_DEPTH_mean_LAB_profile_Pds[i] = [RF_DEPTH_mean_LAB_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_LAB_profile_Pds[i] = [RF_DEPTH_std_LAB_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

else:

	#Profile lat/lon

	AB_lon = [[]]*len(rows[0,:])
	AB_lat = [[]]*len(rows[0,:])

	#Receiver Functions:

	RF_data_profile_Pds = [[]]*len(rows[0,:])

	#Estimates MOHO:

	RF_DEPTH_mean_MOHO_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_MOHO_profile_Pds = [[]]*len(rows[0,:])

	#Estimates LAB:

	RF_DEPTH_mean_LAB_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_LAB_profile_Pds = [[]]*len(rows[0,:])


	for i,j in enumerate(rows[0,:]):
		
		lat_lon = [(lons[k],lats[k]) for k,l in enumerate(lats)]
		grid_column = [(float("%.4f" % round(rows[:,i][k],4)),float("%.4f" % round(cols[:,i][k],4))) for k,l in enumerate(rows[:,i])]

		#Profile lat/lon

		AB_lon[i] = [rows[:,i][k] for k,l in enumerate(rows[:,i])]
		AB_lat[i] = [cols[:,i][k] for k,l in enumerate(rows[:,i])]

		#Receiver Functions:

		RF_data_profile_Pds[i] = [RF_stacking_Pds[lat_lon.index(l)] if l in lat_lon else np.zeros_like(RF_stacking_Pds[k]) for k,l in enumerate(grid_column)]

		#Estimates MOHO:

		RF_DEPTH_mean_MOHO_profile_Pds[i] = [RF_DEPTH_mean_MOHO_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_MOHO_profile_Pds[i] = [RF_DEPTH_std_MOHO_Pds[lat_lon.index(l)]  if l in lat_lon else np.nan for k,l in enumerate(grid_column)]


		#Estimates LAB:

		RF_DEPTH_mean_LAB_profile_Pds[i] = [RF_DEPTH_mean_LAB_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_LAB_profile_Pds[i] = [RF_DEPTH_std_LAB_Pds[lat_lon.index(l)]  if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

print('Plotting cross-sections according to the '+CROSS_SECTION_AXIS+' direction')
for i,j in enumerate(RF_data_profile_Pds):

	#Cross section figure

	fig = plt.figure(figsize=(30, 10))

	fig.suptitle('Cross section for Pds')

	gs = gridspec.GridSpec(4, 3)
	gs.update(wspace=0.2, hspace=0.5)

	pefil_pds = fig.add_subplot(gs[0:4,1:])

	#_____________________________________________

	map_MOHO_thickness =  fig.add_subplot(gs[0:2,0],projection=ccrs.PlateCarree())

	#_____________________________________________

	apparent_MOHO = fig.add_subplot(gs[2,0])
	apparent_LAB = fig.add_subplot(gs[3,0],sharex=apparent_MOHO)


	#######################################################################

	colormap = plt.get_cmap(COLORMAP_VEL)

	colormap_std = plt.get_cmap(COLORMAP_STD)

	colormap_segmentation = INTER_DEPTH/10
	n=41
	x = 0.5
	upper = plt.cm.seismic(np.linspace(0, x, n)[::-1])
	white = plt.cm.seismic(np.ones(18)*0.5)
	lower= plt.cm.seismic(np.linspace(1-x, 1, n)[::-1])
	colors = np.vstack((lower, white,upper))
	tmap = matplotlib.colors.LinearSegmentedColormap.from_list('map_white', colors)

	map_MOHO_thickness.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
	map_MOHO_thickness.yaxis.set_ticks_position('both')
	map_MOHO_thickness.xaxis.set_ticks_position('both')

	map_MOHO_thickness.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE,4), crs=ccrs.PlateCarree())
	map_MOHO_thickness.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE,4), crs=ccrs.PlateCarree())
	map_MOHO_thickness.tick_params(labelbottom=False,labeltop=True,labelleft=True,labelright=True)

	map_MOHO_thickness.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

	reader_1_SHP = Reader(BOUNDARY_1_SHP)
	shape_1_SHP = list(reader_1_SHP.geometries())
	plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
	map_MOHO_thickness.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

	reader_2_SHP = Reader(BOUNDARY_2_SHP)
	shape_2_SHP = list(reader_2_SHP.geometries())
	plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
	map_MOHO_thickness.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

	bounds = np.arange(30, 50+colormap_segmentation, colormap_segmentation)
	norm_map_MOHO_thickness = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)

	for t,y in enumerate(lons):
		if math.isnan(RF_DEPTH_mean_MOHO_Pds[t]) == False:
			if RF_DEPTH_std_MOHO_Pds[t] < INTER_DEPTH:
				circulo_MOHO = Circle(radius=DIST_GRID_PP,xy=(lons[t], lats[t]),color=tmap(norm_map_MOHO_thickness(RF_DEPTH_mean_MOHO_Pds[t])),ec='k',linewidth=0.5,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
				map_MOHO_thickness.add_patch(circulo_MOHO)
			else:
				circulo_MOHO = Circle(radius=DIST_GRID_PP,xy=(lons[t], lats[t]),color=tmap(norm_map_MOHO_thickness(RF_DEPTH_mean_MOHO_Pds[t])),  ec='k',linewidth=0.5,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
				map_MOHO_thickness.add_patch(circulo_MOHO)
		else:
			pass

	for x,c in enumerate(AB_lon[i]):
		circulo_MOHO_profile = Circle(radius=DIST_GRID_PP,xy=(AB_lon[i][x], AB_lat[i][x]),fc='None',ec='k',transform=ccrs.Geodetic(),zorder=10)
		map_MOHO_thickness.add_patch(circulo_MOHO_profile)

	map_MOHO_thickness.set_title('MOHO Thickness', y=1.1)

	sm_map_MOHO_thickness = plt.cm.ScalarMappable(cmap=tmap,norm=norm_map_MOHO_thickness)
	sm_map_MOHO_thickness._A = []
	cbar = fig.colorbar(sm_map_MOHO_thickness,ax=map_MOHO_thickness,orientation='vertical',shrink=0.9,pad=0.1,label='MOHO Thickness (km)')

	cbar.set_ticks(np.arange(30, 50+INTER_DEPTH, INTER_DEPTH))
	cbar.set_ticklabels(np.arange(30, 50+INTER_DEPTH, INTER_DEPTH))

	#### Profile  ####
	
	majorLocatorY = MultipleLocator(10)
	minorLocatorY = MultipleLocator(5)

	grid_Pds = np.array(RF_data_profile_Pds[i])
	extent_Pds = [0,len(RF_data_profile_Pds[i]),MAX_DEPTH,MIN_DEPTH]


	im = pefil_pds.imshow(grid_Pds.T,extent=extent_Pds,interpolation='hamming', cmap='seismic',vmin=VMIN,vmax=VMAX)
	pefil_pds.set_aspect('auto')
	axins = inset_axes(pefil_pds,
                   width="10%",  # width = 10% of parent_bbox width
                   height="5%",  # height : 50%
                   loc='upper left',
                   bbox_to_anchor=(0.85, 0.07, 1, 1),
                   bbox_transform=pefil_pds.transAxes,
                   borderpad=0,
                   )
	plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')
	for _i, _j in enumerate(RF_data_profile_Pds[i]):
		pefil_pds.plot(_i,RF_DEPTH_mean_MOHO_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')
		pefil_pds.plot(_i,RF_DEPTH_mean_LAB_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')

		pefil_pds.errorbar(_i,RF_DEPTH_mean_MOHO_profile_Pds[i][_i], yerr=RF_DEPTH_std_MOHO_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_pds.errorbar(_i,RF_DEPTH_mean_LAB_profile_Pds[i][_i], yerr=RF_DEPTH_std_LAB_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)

		if CROSS_SECTION_AXIS == 'y':
			pefil_pds.text(_i,820,"{0:.1f}".format(AB_lon[i][_i]),rotation=-45,fontsize=10)
			pefil_pds.set_xlabel('Longitude ($^\circ$)',labelpad=30)
		else:
			pefil_pds.text(_i,820,"{0:.1f}".format(AB_lat[i][_i]),rotation=-45,fontsize=10)
			pefil_pds.set_xlabel('Latitude ($^\circ$)',labelpad=30)

	pefil_pds.yaxis.set_ticks_position('both')
	pefil_pds.yaxis.set_major_locator(majorLocatorY)
	pefil_pds.yaxis.set_minor_locator(minorLocatorY)
	pefil_pds.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
	pefil_pds.set_title('Cross-section - Pds')
	pefil_pds.set_xticks([])
	pefil_pds.set_ylabel('Depth (km)')
	pefil_pds.yaxis.set_label_position("right")
	pefil_pds.tick_params(labelright=True)

	#### Figure Apparent  MOHO  ####

	for _i, _j in enumerate(RF_data_profile_Pds[i]):
		apparent_MOHO.plot(_i,RF_DEPTH_mean_MOHO_profile_Pds[i][_i],marker='o',markerfacecolor='none',markeredgecolor='dimgray')

		apparent_MOHO.errorbar(_i,RF_DEPTH_mean_MOHO_profile_Pds[i][_i], yerr=RF_DEPTH_std_MOHO_profile_Pds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)

		apparent_MOHO.set_title('MOHO')
		apparent_MOHO.yaxis.set_ticks_position('both')
		apparent_MOHO.yaxis.set_major_locator(MultipleLocator(10))
		apparent_MOHO.yaxis.set_minor_locator(MultipleLocator(5))
		apparent_MOHO.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_MOHO.tick_params(labelleft=True,labelright=True)
		apparent_MOHO.yaxis.set_label_position("right")
		apparent_MOHO.set_xticks([])
		apparent_MOHO.set_ylim(60,20)


	#### Figure Apparent  LAB  ####


	for _i, _j in enumerate(RF_data_profile_Pds[i]):
		apparent_LAB.plot(_i,RF_DEPTH_mean_LAB_profile_Pds[i][_i],marker='o',markerfacecolor='none',markeredgecolor='dimgray')

		apparent_LAB.errorbar(_i,RF_DEPTH_mean_LAB_profile_Pds[i][_i], yerr=RF_DEPTH_std_LAB_profile_Pds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)
		
		apparent_LAB.set_title('LAB ')
		apparent_LAB.yaxis.set_ticks_position('both')
		apparent_LAB.yaxis.set_major_locator(MultipleLocator(10))
		apparent_LAB.yaxis.set_minor_locator(MultipleLocator(5))
		apparent_LAB.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_LAB.tick_params(labelleft=True,labelright=True)
		apparent_LAB.yaxis.set_label_position("right")
		apparent_LAB.set_xticks([])
		apparent_LAB.set_ylim(200,150)


	fig.savefig(RESULTS_FOLDER+'SELECTED_BINNED_DATA_'+CROSS_SECTION_AXIS+'_CROSS_SECTION_Pds_Ppds_PROFILE_'+str(i+1)+'_COLOR.'+EXT_FIG,dpi=DPI_FIG)
print('Ending the Cross-section CODE')