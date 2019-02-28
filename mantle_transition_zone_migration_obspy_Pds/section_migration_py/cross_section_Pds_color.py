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
from shapely.geometry import Polygon, MultiPoint, Point, LineString
import shapefile
from matplotlib.colors import Normalize
from matplotlib.patches import Circle,Rectangle
import math
from mpl_toolkits.axes_grid1.inset_locator import inset_axes







from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,TECTO_SHP,COLORMAP_VEL,COLORMAP_STD,OUTPUT_DIR,
					EXT_FIG,DPI_FIG,DIST_GRID_PP,NUMBER_STA_PER_BIN,NUMBER_PP_PER_BIN,VMIN,VMAX,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,CROSS_SECTION_AXIS,DEPTH_TARGET,
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

#Estimates LVZ:

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

#Estimates MTZ Pds:

RF_DEPTH_mtz_thickness_Pds = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds']
RF_DEPTH_mtz_thickness_Pds_std = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std']

#############################################################################################################################3


area = (LLCRNRLON_SMALL,URCRNRLON_SMALL, LLCRNRLAT_SMALL, URCRNRLAT_SMALL)

shape = (int(abs(abs(URCRNRLON_SMALL) - abs(LLCRNRLON_SMALL))*GRID_PP_MULT),int(abs(abs(URCRNRLAT_SMALL) - abs(LLCRNRLAT_SMALL))*GRID_PP_MULT))

grdx, grdy = gridder.regular(area, shape)

shape_new = (int(len(set(grdx))),int(len(set(grdy))))

rows = np.array(grdx).reshape(shape_new)
cols = np.array(grdy).reshape(shape_new)

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

	#Estimates LVZ:

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

		#Estimates LVZ:

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

		#Estimates MTZ Pds:

		RF_DEPTH_mtz_thickness_profile_Pds[i] = [RF_DEPTH_mtz_thickness_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_mtz_thickness_profile_Pds_std[i] = [RF_DEPTH_mtz_thickness_Pds_std[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

else:

	#Profile lat/lon

	AB_lon = [[]]*len(rows[0,:])
	AB_lat = [[]]*len(rows[0,:])

	#Receiver Functions:

	RF_data_profile_Pds = [[]]*len(rows[0,:])

	#Estimates P410s:

	RF_DEPTH_mean_1_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_1_profile_Pds = [[]]*len(rows[0,:])

	#Estimates LVZ:

	RF_DEPTH_mean_LVZ_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_LVZ_profile_Pds = [[]]*len(rows[0,:])

	#Estimates P520s:

	RF_DEPTH_mean_520_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_520_profile_Pds = [[]]*len(rows[0,:])

	#Estimates P660s:

	RF_DEPTH_mean_2_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_2_profile_Pds = [[]]*len(rows[0,:])

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

		#Estimates LVZ:

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

		#Estimates MTZ Pds:

		RF_DEPTH_mtz_thickness_profile_Pds[i] = [RF_DEPTH_mtz_thickness_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_mtz_thickness_profile_Pds_std[i] = [RF_DEPTH_mtz_thickness_Pds_std[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

print('Plotting cross-sections according to the '+CROSS_SECTION_AXIS+' direction')
for i,j in enumerate(RF_data_profile_Pds):

	#Cross section figure

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
	apparent_660 = fig.add_subplot(gs[3,0])


	#######################################################################

	colormap = plt.get_cmap(COLORMAP_VEL)

	colormap_std = plt.get_cmap(COLORMAP_STD)


	map_MTZ_thickness.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
	map_MTZ_thickness.yaxis.set_ticks_position('both')
	map_MTZ_thickness.xaxis.set_ticks_position('both')

	map_MTZ_thickness.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE,4), crs=ccrs.PlateCarree())
	map_MTZ_thickness.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE,4), crs=ccrs.PlateCarree())
	map_MTZ_thickness.tick_params(labelbottom=False,labeltop=True,labelleft=True,labelright=True)

	map_MTZ_thickness.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

	reader_1_SHP = Reader(BOUNDARY_1_SHP)
	shape_1_SHP = list(reader_1_SHP.geometries())
	plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
	map_MTZ_thickness.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

	reader_2_SHP = Reader(BOUNDARY_2_SHP)
	shape_2_SHP = list(reader_2_SHP.geometries())
	plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
	map_MTZ_thickness.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

	norm_map_MTZ_thickness = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)

	for t,y in enumerate(lons):
		if math.isnan(RF_DEPTH_mtz_thickness_Pds[t]) == False:
			if RF_DEPTH_mtz_thickness_Pds_std[t] < 10:
				circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[t], lats[t]),color=colormap(norm_map_MTZ_thickness(RF_DEPTH_mtz_thickness_Pds[t])),ec='None',linewidth=1,linestyle='-',transform=ccrs.Geodetic(),zorder=2)
				map_MTZ_thickness.add_patch(circulo_410)
			else:
				circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[t], lats[t]),color=colormap(norm_map_MTZ_thickness(RF_DEPTH_mtz_thickness_Pds[t])),  ec='k',linewidth=1,linestyle=':',transform=ccrs.Geodetic(),zorder=2)
				map_MTZ_thickness.add_patch(circulo_410)
		else:
			pass

	for x,c in enumerate(AB_lon[i]):
		circulo_410_profile = Circle(radius=DIST_GRID_PP,xy=(AB_lon[i][x], AB_lat[i][x]),fc='None',ec='k',transform=ccrs.Geodetic(),zorder=10)
		map_MTZ_thickness.add_patch(circulo_410_profile)

	map_MTZ_thickness.set_title('MTZ Thickness', y=1.1)

	sm_map_MTZ_thickness = plt.cm.ScalarMappable(cmap=colormap,norm=norm_map_MTZ_thickness)
	sm_map_MTZ_thickness._A = []
	fig.colorbar(sm_map_MTZ_thickness,ax=map_MTZ_thickness,orientation='vertical',shrink=0.9,pad=0.1,label='Thickness (km)')

	#### Figure Pds  ####
	
	majorLocatorY = MultipleLocator(50)
	minorLocatorY = MultipleLocator(10)

	grid_Pds = np.array(RF_data_profile_Pds[i])
	extent_Pds = [0,len(RF_data_profile_Pds[i]),800,300]


	im = pefil_pds.imshow(grid_Pds.T,extent=extent_Pds,interpolation='bicubic', cmap='seismic',vmin=VMIN,vmax=VMAX)
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
		pefil_pds.plot(_i,RF_DEPTH_mean_LVZ_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')
		pefil_pds.plot(_i,RF_DEPTH_mean_1_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')
		pefil_pds.plot(_i,RF_DEPTH_mean_520_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')
		pefil_pds.plot(_i,RF_DEPTH_mean_2_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')

		pefil_pds.errorbar(_i,RF_DEPTH_mean_LVZ_profile_Pds[i][_i], yerr=RF_DEPTH_std_LVZ_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_pds.errorbar(_i,RF_DEPTH_mean_1_profile_Pds[i][_i], yerr=RF_DEPTH_std_1_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_pds.errorbar(_i,RF_DEPTH_mean_520_profile_Pds[i][_i], yerr=RF_DEPTH_std_520_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_pds.errorbar(_i,RF_DEPTH_mean_2_profile_Pds[i][_i], yerr=RF_DEPTH_std_2_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)

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

	#### Figure Apparent  410 km Pds  ####

	for _i, _j in enumerate(RF_data_profile_Pds[i]):
		apparent_410.plot(_i,RF_DEPTH_mean_1_profile_Pds[i][_i]-410,marker='o',markerfacecolor='none',markeredgecolor='dimgray')

		apparent_410.errorbar(_i,RF_DEPTH_mean_1_profile_Pds[i][_i]-410, yerr=RF_DEPTH_std_1_profile_Pds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)

		apparent_410.set_title('diff 410 km')
		apparent_410.yaxis.set_ticks_position('both')
		apparent_410.yaxis.set_major_locator(MultipleLocator(25))
		apparent_410.yaxis.set_minor_locator(MultipleLocator(10))
		apparent_410.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_410.tick_params(labelleft=True,labelright=True)
		apparent_410.yaxis.set_label_position("right")
		apparent_410.set_xticks([])
		apparent_410.set_ylim(-50,50)


	#### Figure Apparent  660 km Pds  ####


	for _i, _j in enumerate(RF_DEPTH_mean_2_profile_Pds[i]):
		apparent_660.plot(_i,RF_DEPTH_mean_2_profile_Pds[i][_i]-660,marker='o',markerfacecolor='none',markeredgecolor='dimgray')

		apparent_660.errorbar(_i,RF_DEPTH_mean_2_profile_Pds[i][_i]-660, yerr=RF_DEPTH_std_2_profile_Pds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)
		
		apparent_660.set_ylim(-50,50)
		apparent_660.set_title('diff 660 km ')
		apparent_660.yaxis.set_ticks_position('both')
		apparent_660.yaxis.set_major_locator(MultipleLocator(25))
		apparent_660.yaxis.set_minor_locator(MultipleLocator(10))
		apparent_660.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_660.tick_params(labelleft=True,labelright=True)
		apparent_660.yaxis.set_label_position("right")
		apparent_660.set_xticks([])


	#### Figure MTZ Apparent thickness  ####

	for _i, _j in enumerate(RF_DEPTH_mtz_thickness_profile_Pds[i]):
		MTZ_thickness.plot(_i,RF_DEPTH_mtz_thickness_profile_Pds[i][_i]-250,marker='o',markerfacecolor='none',markeredgecolor='gray')

	for _i, _j in enumerate(RF_DEPTH_mtz_thickness_profile_Pds[i]):
		MTZ_thickness.errorbar(_i,RF_DEPTH_mtz_thickness_profile_Pds[i][_i]-250, yerr=RF_DEPTH_mtz_thickness_profile_Pds_std[i][_i], ecolor='gray',elinewidth=1,capsize=2,capthick=1)
	
	MTZ_thickness.set_ylim(-50,50)
	MTZ_thickness.yaxis.set_label_position("right")
	MTZ_thickness.set_title('diff MTZ Thickness')
	MTZ_thickness.yaxis.set_ticks_position('both')
	MTZ_thickness.yaxis.set_major_locator(MultipleLocator(25))
	MTZ_thickness.yaxis.set_minor_locator(MultipleLocator(10))
	MTZ_thickness.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
	MTZ_thickness.tick_params(labelleft=True,labelright=True)


	#### Figure diff MTZ Apparent thickness  ####

	for _i, _j in enumerate(RF_DEPTH_mean_LVZ_profile_Pds[i]):
		diff_MTZ_thickness.plot(_i,RF_DEPTH_mean_LVZ_profile_Pds[i][_i]-350,marker='o',markerfacecolor='none',markeredgecolor='gray')

	for _i, _j in enumerate(RF_DEPTH_mean_LVZ_profile_Pds[i]):
		diff_MTZ_thickness.errorbar(_i,RF_DEPTH_mean_LVZ_profile_Pds[i][_i]-350, yerr=RF_DEPTH_std_LVZ_profile_Pds[i][_i], ecolor='gray',elinewidth=1,capsize=2,capthick=1)
	
	diff_MTZ_thickness.set_ylim(-100,100)
	diff_MTZ_thickness.yaxis.set_label_position("right")
	diff_MTZ_thickness.set_title('diff LVZ')
	diff_MTZ_thickness.yaxis.set_ticks_position('both')
	diff_MTZ_thickness.yaxis.set_major_locator(MultipleLocator(50))
	diff_MTZ_thickness.yaxis.set_minor_locator(MultipleLocator(20))
	diff_MTZ_thickness.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
	diff_MTZ_thickness.tick_params(labelleft=True,labelright=True)

	fig.savefig(RESULTS_FOLDER+'SELECTED_BINNED_DATA_'+CROSS_SECTION_AXIS+'_CROSS_SECTION_Pds_Ppds_PROFILE_'+str(i+1)+'_COLOR.'+EXT_FIG,dpi=DPI_FIG)
print('Ending the Cross-section CODE')