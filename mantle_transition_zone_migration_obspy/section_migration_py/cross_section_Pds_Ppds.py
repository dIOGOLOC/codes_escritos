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






from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,TECTO_SHP,COLORMAP_VEL,COLORMAP_STD,OUTPUT_DIR,
					EXT_FIG,DPI_FIG,FRESNEL_ZONE_RADIUS,DIST_GRID_PP,NUMBER_STA_PER_BIN,NUMBER_PP_PER_BIN,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,BOOTSTRAP_DEPTH_ESTIMATION,GAMMA,CROSS_SECTION_AXIS
				   )


print('Starting Cross section CODE')
print('\n')

print('Calculating earth model layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'/'+'Stations'+'/'

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

PP_SELEC_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'/'+'SELECTED_BINNED_DATA'+'/'


RESULTS_FOLDER_BINS = PP_SELEC_DIR+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
filename = RESULTS_FOLDER_BINS+'SELECTED_BINNED.json'

SELECTED_BINNED_DATA_dic = json.load(open(filename))

#Lat/Lon lists:

lats = SELECTED_BINNED_DATA_dic['lat']
lons = SELECTED_BINNED_DATA_dic['lon']

lats_true = SELECTED_BINNED_DATA_dic['lat_true']
lons_true = SELECTED_BINNED_DATA_dic['lon_true']

#Receiver Functions:

RF_number = SELECTED_BINNED_DATA_dic['len_Pds']

RF_stacking_Pds = SELECTED_BINNED_DATA_dic['data_Pds']
RF_stacking_Ppds = SELECTED_BINNED_DATA_dic['data_Ppds']

#Estimates P410s:

RF_DEPTH_mean_1_Pds = SELECTED_BINNED_DATA_dic['mean_1_Pds']
RF_DEPTH_std_1_Pds = SELECTED_BINNED_DATA_dic['std_1_Pds']

#Estimates P660s:

RF_DEPTH_mean_2_Pds = SELECTED_BINNED_DATA_dic['mean_2_Pds']
RF_DEPTH_std_2_Pds = SELECTED_BINNED_DATA_dic['std_2_Pds']

#Estimates Pp410s:

RF_DEPTH_mean_1_Ppds = SELECTED_BINNED_DATA_dic['mean_1_Ppds']
RF_DEPTH_std_1_Ppds = SELECTED_BINNED_DATA_dic['std_1_Ppds']

#Estimates Pp660s:

RF_DEPTH_mean_2_Ppds = SELECTED_BINNED_DATA_dic['mean_2_Ppds']
RF_DEPTH_std_2_Ppds = SELECTED_BINNED_DATA_dic['std_2_Ppds']

#Estimates MTZ Pds:

RF_DEPTH_mtz_thickness_Pds = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds']
RF_DEPTH_mtz_thickness_Pds_std = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std']

#Estimates MTZ Ppds:

RF_DEPTH_mtz_thickness_Ppds = SELECTED_BINNED_DATA_dic['mtz_thickness_Ppds']
RF_DEPTH_mtz_thickness_Ppds_std = SELECTED_BINNED_DATA_dic['mtz_thickness_Ppds_std']

#Estimates True MTZ:

RF_DEPTH_true_thickness_MTZ = SELECTED_BINNED_DATA_dic['true_thickness_MTZ_Pds']
RF_DEPTH_true_thickness_MTZ_std = SELECTED_BINNED_DATA_dic['true_thickness_MTZ_Pds_std']

#Estimates True Depth 410:

RF_DEPTH_mean_1_true = SELECTED_BINNED_DATA_dic['true_mean_1_Pds']
RF_DEPTH_std_1_true = SELECTED_BINNED_DATA_dic['true_std_1_Pds']

#Estimates True Depth 660:

RF_DEPTH_mean_2_true = SELECTED_BINNED_DATA_dic['true_mean_2_Pds']
RF_DEPTH_std_2_true = SELECTED_BINNED_DATA_dic['true_std_2_Pds']

#############################################################################################################################3


area = (LLCRNRLON_SMALL,URCRNRLON_SMALL, LLCRNRLAT_SMALL, URCRNRLAT_SMALL)

shape = (int(abs(abs(URCRNRLON_SMALL) - abs(LLCRNRLON_SMALL))*GRID_PP_MULT),int(abs(abs(URCRNRLAT_SMALL) - abs(LLCRNRLAT_SMALL))*GRID_PP_MULT))

grdx, grdy = gridder.regular(area, shape)

shape_new = (int(len(set(grdx))),int(len(set(grdy))))

rows = np.array(grdx).reshape(shape_new)
cols = np.array(grdy).reshape(shape_new)

PP_FIGURE = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'/'+'Figures'+'/'
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
	RF_data_profile_Ppds = [[]]*len(rows[:,0])

	#Estimates P410s:

	RF_DEPTH_mean_1_profile_Pds = [[]]*len(rows[:,0])
	RF_DEPTH_std_1_profile_Pds = [[]]*len(rows[:,0])

	#Estimates P660s:

	RF_DEPTH_mean_2_profile_Pds = [[]]*len(rows[:,0])
	RF_DEPTH_std_2_profile_Pds = [[]]*len(rows[:,0])

	#Estimates Pp410s:

	RF_DEPTH_mean_1_profile_Ppds = [[]]*len(rows[:,0])
	RF_DEPTH_std_1_profile_Ppds = [[]]*len(rows[:,0])

	#Estimates Pp660s:

	RF_DEPTH_mean_2_profile_Ppds = [[]]*len(rows[:,0])
	RF_DEPTH_std_2_profile_Ppds = [[]]*len(rows[:,0])


	#Estimates MTZ Pds:

	RF_DEPTH_mtz_thickness_profile_Pds = [[]]*len(rows[:,0]) 
	RF_DEPTH_mtz_thickness_profile_Pds_std = [[]]*len(rows[:,0]) 

	#Estimates MTZ Ppds:

	RF_DEPTH_mtz_thickness_profile_Ppds = [[]]*len(rows[:,0])
	RF_DEPTH_mtz_thickness_profile_Ppds_std = [[]]*len(rows[:,0])

	#Estimates True Depth 410:

	RF_DEPTH_mean_1_true_profile = [[]]*len(rows[:,0])
	RF_DEPTH_std_1_true_profile = [[]]*len(rows[:,0])

	#Estimates True Depth 660:

	RF_DEPTH_mean_2_true_profile = [[]]*len(rows[:,0])
	RF_DEPTH_std_2_true_profile = [[]]*len(rows[:,0])

	#Estimates True MTZ:

	RF_DEPTH_true_thickness_MTZ_profile = [[]]*len(rows[:,0]) 
	RF_DEPTH_true_thickness_MTZ_profile_std = [[]]*len(rows[:,0]) 

	for i,j in enumerate(rows[:,0]):

		lat_lon = [(lons[k],lats[k]) for k,l in enumerate(lats)]
		lat_lon_true = [(lons_true[k],lats_true[k]) for k,l in enumerate(lat_true)]
		grid_column = [(float("%.4f" % round(rows[i,:][k],4)),float("%.4f" % round(cols[i,:][k],4))) for k,l in enumerate(rows[i,:])]

		#Profile lat/lon

		AB_lon[i] = [rows[i,:][k] for k,l in enumerate(rows[i,:])]
		AB_lat[i] = [cols[i,:][k] for k,l in enumerate(rows[i,:])]


		#Receiver Functions:

		RF_data_profile_Pds[i] = [RF_stacking_Pds[lat_lon.index(l)] if l in lat_lon else np.zeros_like(RF_stacking_Pds[k]) for k,l in enumerate(grid_column)]
		RF_data_profile_Ppds[i] = [RF_stacking_Ppds[lat_lon.index(l)] if l in lat_lon else np.zeros_like(RF_stacking_Ppds[k]) for k,l in enumerate(grid_column)]

		#Estimates P410s:

		RF_DEPTH_mean_1_profile_Pds[i] = [RF_DEPTH_mean_1_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_1_profile_Pds[i] = [RF_DEPTH_std_1_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates P660s:

		RF_DEPTH_mean_2_profile_Pds[i] = [RF_DEPTH_mean_2_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_2_profile_Pds[i] = [RF_DEPTH_std_2_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates Pp410s:

		RF_DEPTH_mean_1_profile_Ppds[i] = [RF_DEPTH_mean_1_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_1_profile_Ppds[i] = [RF_DEPTH_std_1_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates Pp660s:

		RF_DEPTH_mean_2_profile_Ppds[i] = [RF_DEPTH_mean_2_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_2_profile_Ppds[i] = [RF_DEPTH_std_2_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates MTZ Pds:

		RF_DEPTH_mtz_thickness_profile_Pds[i] = [RF_DEPTH_mtz_thickness_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_mtz_thickness_profile_Pds_std[i] = [RF_DEPTH_mtz_thickness_Pds_std[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates MTZ Ppds:

		RF_DEPTH_mtz_thickness_profile_Ppds[i] = [RF_DEPTH_mtz_thickness_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_mtz_thickness_profile_Ppds_std[i] = [RF_DEPTH_mtz_thickness_Ppds_std[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates True Depth 410:

		RF_DEPTH_mean_1_true_profile[i] = [RF_DEPTH_mean_1_true[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_1_true_profile[i] = [RF_DEPTH_std_1_true[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates True Depth 660:

		RF_DEPTH_mean_2_true_profile[i] = [RF_DEPTH_mean_2_true[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_2_true_profile[i] = [RF_DEPTH_std_2_true[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates True MTZ:

		RF_DEPTH_true_thickness_MTZ_profile[i] = [RF_DEPTH_true_thickness_MTZ[lat_lon.index(l)] if l in lat_lon_true else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_true_thickness_MTZ_profile_std[i] = [RF_DEPTH_true_thickness_MTZ_std[lat_lon.index(l)] if l in lat_lon_true else np.nan for k,l in enumerate(grid_column)]


else:

	#Profile lat/lon

	AB_lon = [[]]*len(rows[0,:])
	AB_lat = [[]]*len(rows[0,:])


	#Receiver Functions:

	RF_data_profile_Pds = [[]]*len(rows[0,:])
	RF_data_profile_Ppds = [[]]*len(rows[0,:])

	#Estimates P410s:

	RF_DEPTH_mean_1_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_1_profile_Pds = [[]]*len(rows[0,:])

	#Estimates P660s:

	RF_DEPTH_mean_2_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_std_2_profile_Pds = [[]]*len(rows[0,:])

	#Estimates Pp410s:

	RF_DEPTH_mean_1_profile_Ppds = [[]]*len(rows[0,:])
	RF_DEPTH_std_1_profile_Ppds = [[]]*len(rows[0,:])

	#Estimates Pp660s:

	RF_DEPTH_mean_2_profile_Ppds = [[]]*len(rows[0,:])
	RF_DEPTH_std_2_profile_Ppds = [[]]*len(rows[0,:])


	#Estimates MTZ Pds:

	RF_DEPTH_mtz_thickness_profile_Pds = [[]]*len(rows[0,:])
	RF_DEPTH_mtz_thickness_profile_Pds_std = [[]]*len(rows[0,:])

	#Estimates MTZ Ppds:

	RF_DEPTH_mtz_thickness_profile_Ppds = [[]]*len(rows[0,:])
	RF_DEPTH_mtz_thickness_profile_Ppds_std = [[]]*len(rows[0,:])

	#Estimates True Depth 410:

	RF_DEPTH_mean_1_true_profile = [[]]*len(rows[0,:])
	RF_DEPTH_std_1_true_profile = [[]]*len(rows[0,:])

	#Estimates True Depth 660:

	RF_DEPTH_mean_2_true_profile = [[]]*len(rows[0,:])
	RF_DEPTH_std_2_true_profile = [[]]*len(rows[0,:])

	#Estimates True MTZ:

	RF_DEPTH_true_thickness_MTZ_profile = [[]]*len(rows[0,:])
	RF_DEPTH_true_thickness_MTZ_profile_std = [[]]*len(rows[0,:])

	for i,j in enumerate(rows[0,:]):
		
		lat_lon = [(lons[k],lats[k]) for k,l in enumerate(lats)]
		grid_column = [(float("%.4f" % round(rows[:,i][k],4)),float("%.4f" % round(cols[:,i][k],4))) for k,l in enumerate(rows[:,i])]

		#Profile lat/lon

		AB_lon[i] = [rows[:,i][k] for k,l in enumerate(rows[:,i])]
		AB_lat[i] = [cols[:,i][k] for k,l in enumerate(rows[:,i])]

		#Receiver Functions:

		RF_data_profile_Pds[i] = [RF_stacking_Pds[lat_lon.index(l)] if l in lat_lon else np.zeros_like(RF_stacking_Pds[k]) for k,l in enumerate(grid_column)]
		RF_data_profile_Ppds[i] = [RF_stacking_Ppds[lat_lon.index(l)] if l in lat_lon else np.zeros_like(RF_stacking_Ppds[k]) for k,l in enumerate(grid_column)]

		#Estimates P410s:

		RF_DEPTH_mean_1_profile_Pds[i] = [RF_DEPTH_mean_1_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_1_profile_Pds[i] = [RF_DEPTH_std_1_Pds[lat_lon.index(l)]  if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates P660s:

		RF_DEPTH_mean_2_profile_Pds[i] = [RF_DEPTH_mean_2_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_2_profile_Pds[i] = [RF_DEPTH_std_2_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates Pp410s:

		RF_DEPTH_mean_1_profile_Ppds[i] = [RF_DEPTH_mean_1_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_1_profile_Ppds[i] = [RF_DEPTH_std_1_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates Pp660s:

		RF_DEPTH_mean_2_profile_Ppds[i] = [RF_DEPTH_mean_2_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_2_profile_Ppds[i] = [RF_DEPTH_std_2_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates MTZ Pds:

		RF_DEPTH_mtz_thickness_profile_Pds[i] = [RF_DEPTH_mtz_thickness_Pds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_mtz_thickness_profile_Pds_std[i] = [RF_DEPTH_mtz_thickness_Pds_std[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates MTZ Ppds:

		RF_DEPTH_mtz_thickness_profile_Ppds[i] = [RF_DEPTH_mtz_thickness_Ppds[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_mtz_thickness_profile_Ppds_std[i] = [RF_DEPTH_mtz_thickness_Ppds_std[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]


		#Estimates True Depth 410:

		RF_DEPTH_mean_1_true_profile[i] = [RF_DEPTH_mean_1_true[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_1_true_profile[i] = [RF_DEPTH_std_1_true[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates True Depth 660:

		RF_DEPTH_mean_2_true_profile[i] = [RF_DEPTH_mean_2_true[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_std_2_true_profile[i] = [RF_DEPTH_std_2_true[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

		#Estimates True MTZ:

		RF_DEPTH_true_thickness_MTZ_profile[i] = [RF_DEPTH_true_thickness_MTZ[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]
		RF_DEPTH_true_thickness_MTZ_profile_std[i] = [RF_DEPTH_true_thickness_MTZ_std[lat_lon.index(l)] if l in lat_lon else np.nan for k,l in enumerate(grid_column)]

	
print('Plotting cross-sections according to the '+CROSS_SECTION_AXIS+' direction')
for i,j in enumerate(RF_data_profile_Pds):

	#Cross section figure

	fig = plt.figure(figsize=(25, 10))

	fig.suptitle('Cross section for Pds and Ppds')

	gs = gridspec.GridSpec(6, 3)
	gs.update(wspace=0.2, hspace=0.5)

	MTZ_thickness = fig.add_subplot(gs[:2,1:])

	pefil_pds = fig.add_subplot(gs[2:4,1:],sharex=MTZ_thickness)
	pefil_ppds = fig.add_subplot(gs[4:6,1:],sharex=MTZ_thickness)

	#_____________________________________________

	apparent_410 = fig.add_subplot(gs[2,0])
	apparent_660 = fig.add_subplot(gs[3,0],sharex=apparent_410)

	apparent_410_ppds = fig.add_subplot(gs[4, 0],sharex=apparent_410)
	apparent_660_ppds = fig.add_subplot(gs[5, 0],sharex=apparent_410)

	#_____________________________________________

	map_MTZ_thickness =  fig.add_subplot(gs[0:2,0], projection=ccrs.PlateCarree())

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

	reader_3_SHP = Reader(TECTO_SHP)
	shape_3_SHP = list(reader_3_SHP.geometries())
	plot_shape_3_SHP = cfeature.ShapelyFeature(shape_3_SHP, ccrs.PlateCarree())
	map_MTZ_thickness.add_feature(plot_shape_3_SHP, facecolor='none', edgecolor='k',linewidth=2)

	norm_map_MTZ_thickness = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)

	for t,y in enumerate(lons_true):
		if math.isnan(RF_DEPTH_true_thickness_MTZ[t]) == False:
			circulo_410 = Circle(radius=DIST_GRID_PP*(1-(RF_DEPTH_true_thickness_MTZ_std[t]/50)),xy=(lons_true[t], lats_true[t]),color=colormap(norm_map_MTZ_thickness(RF_DEPTH_true_thickness_MTZ[t])), ec='None',transform=ccrs.Geodetic(),zorder=2)
			map_MTZ_thickness.add_patch(circulo_410)
		else:
			pass

	for x,c in enumerate(AB_lon[i]):
		circulo_410_profile = Circle(radius=DIST_GRID_PP,xy=(AB_lon[i][x], AB_lat[i][x]),fc='k',alpha=0.5,ec='None',transform=ccrs.Geodetic(),zorder=10)
		map_MTZ_thickness.add_patch(circulo_410_profile)

	map_MTZ_thickness.set_title('MTZ True Thickness', y=1.1)

	sm_map_MTZ_thickness = plt.cm.ScalarMappable(cmap=colormap,norm=norm_map_MTZ_thickness)
	sm_map_MTZ_thickness._A = []
	fig.colorbar(sm_map_MTZ_thickness,ax=map_MTZ_thickness,orientation='vertical',shrink=0.9,pad=0.1,label='Thickness (km)')

	#### Figure Pds  ####


	factor_Pds = 150

	majorLocatorY = MultipleLocator(50)
	minorLocatorY = MultipleLocator(10)

	
	for _i, _j in enumerate(RF_data_profile_Pds[i]):
		RF_data_factor_Pds = [_i/factor_Pds+l for k, l in enumerate(_j)]
		pefil_pds.plot(RF_data_factor_Pds,camadas_terra_10_km,'k',linewidth=1)
		pefil_pds.yaxis.set_major_locator(majorLocatorY)
		pefil_pds.yaxis.set_minor_locator(minorLocatorY)
		pefil_pds.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

		pefil_pds.plot(_i/factor_Pds,RF_DEPTH_mean_1_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')
		pefil_pds.plot(_i/factor_Pds,RF_DEPTH_mean_2_profile_Pds[i][_i],'ok',ms=3,markerfacecolor='none')

		pefil_pds.errorbar(_i/factor_Pds,RF_DEPTH_mean_1_profile_Pds[i][_i], yerr=RF_DEPTH_std_1_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_pds.errorbar(_i/factor_Pds,RF_DEPTH_mean_2_profile_Pds[i][_i], yerr=RF_DEPTH_std_2_profile_Pds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)

		pefil_pds.yaxis.set_ticks_position('both')

		pefil_pds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)>=_i/factor_Pds, facecolor='black',alpha=0.3, interpolate=True)
		pefil_pds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)<=_i/factor_Pds, facecolor='gray',alpha=0.3, interpolate=True)
		
		pefil_pds.set_title('Cross-section - Pds')
		pefil_pds.set_xticks([])
		pefil_pds.set_ylabel('Depth (km)')
		pefil_pds.yaxis.set_label_position("right")
		pefil_pds.tick_params(labelright=True)
		pefil_pds.set_ylim(MAX_DEPTH,MIN_DEPTH)

	#### Figure Ppds  ####


	factor_Ppds = 150

	for _i, _j in enumerate(RF_data_profile_Ppds[i]):
		RF_data_factor_Ppds = [_i/factor_Ppds+l for k, l in enumerate(_j)]
		pefil_ppds.plot(RF_data_factor_Ppds,camadas_terra_10_km,'k',linewidth=1)
		pefil_ppds.yaxis.set_major_locator(majorLocatorY)
		pefil_ppds.yaxis.set_minor_locator(minorLocatorY)
		pefil_ppds.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

		pefil_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[i][_i],marker='p',ms=3,markerfacecolor='none',markeredgecolor='k')
		pefil_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[i][_i],marker='p',ms=3,markerfacecolor='none',markeredgecolor='k')

		pefil_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[i][_i], yerr=RF_DEPTH_std_1_profile_Ppds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
		pefil_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[i][_i], yerr=RF_DEPTH_std_2_profile_Ppds[i][_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)

		pefil_ppds.yaxis.set_ticks_position('both')
		pefil_ppds.yaxis.set_label_position("right")

		pefil_ppds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)>=_i/factor_Ppds, facecolor='black',alpha=0.3, interpolate=True)
		pefil_ppds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)<=_i/factor_Ppds, facecolor='gray',alpha=0.3, interpolate=True)
		
		pefil_ppds.set_xticks([])
		pefil_ppds.set_title('Cross-section - Ppds')
		pefil_ppds.set_ylabel('Depth (km)')
		pefil_ppds.tick_params(labelleft=True,labelright=True)
		pefil_ppds.set_ylim(MAX_DEPTH,MIN_DEPTH)
		if CROSS_SECTION_AXIS == 'y':
			pefil_ppds.text(_i/factor_Ppds,820,"{0:.1f}".format(AB_lon[i][_i]),rotation=-45,fontsize=10)
			pefil_ppds.set_xlabel('Longitude ($^\circ$)',labelpad=30)
		else:
			pefil_ppds.text(_i/factor_Ppds,820,"{0:.1f}".format(AB_lat[i][_i]),rotation=-45,fontsize=10)
			pefil_ppds.set_xlabel('Latitude ($^\circ$)',labelpad=30)


	#### Figure True and Apparent  410 km Pds  ####


	for _i, _j in enumerate(RF_data_profile_Ppds[i]):
		apparent_410.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Pds[i][_i],marker='o',markerfacecolor='none',markeredgecolor='dimgray')
		apparent_410.plot(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile[i][_i],marker='s',markerfacecolor='none',markeredgecolor='k')

		apparent_410.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Pds[i][_i], yerr=RF_DEPTH_std_1_profile_Pds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)
		apparent_410.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile[i][_i], yerr=RF_DEPTH_std_1_true_profile[i][_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)

		apparent_410.set_title('410 km Pds')
		apparent_410.set_ylabel('Depth (km)')
		apparent_410.yaxis.set_ticks_position('both')
		apparent_410.yaxis.set_major_locator(MultipleLocator(50))
		apparent_410.yaxis.set_minor_locator(MultipleLocator(10))
		apparent_410.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_410.tick_params(labelleft=True,labelright=False)
		apparent_410.set_xticks([])
		apparent_410.set_ylim(500,300)


	#### Figure True and Apparent  660 km Pds  ####


	for _i, _j in enumerate(RF_DEPTH_mean_2_profile_Pds[i]):
		apparent_660.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Pds[i][_i],marker='o',markerfacecolor='none',markeredgecolor='dimgray')
		apparent_660.plot(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile[i][_i],marker='s',markerfacecolor='none',markeredgecolor='k')

		apparent_660.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Pds[i][_i], yerr=RF_DEPTH_std_2_profile_Pds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)
		apparent_660.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile[i][_i], yerr=RF_DEPTH_std_2_true_profile[i][_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)

		apparent_660.set_title('660 km Pds')
		apparent_660.set_ylim(800,600)
		apparent_660.set_ylabel('Depth (km)')
		apparent_660.yaxis.set_ticks_position('both')
		apparent_660.yaxis.set_major_locator(MultipleLocator(50))
		apparent_660.yaxis.set_minor_locator(MultipleLocator(10))
		apparent_660.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_660.tick_params(labelleft=True,labelright=False)
		apparent_660.set_xticks([])

	#### Figure True and Apparent  410 km Ppds  ####


	for _i, _j in enumerate(RF_DEPTH_mean_2_profile_Pds[i]):
		apparent_410_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[i][_i],marker='p',markerfacecolor='none',markeredgecolor='dimgray')
		apparent_410_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile[i][_i],marker='s',markerfacecolor='none',markeredgecolor='k')

		apparent_410_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[i][_i], yerr=RF_DEPTH_std_1_profile_Ppds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)
		apparent_410_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile[i][_i], yerr=RF_DEPTH_std_1_true_profile[i][_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)

		apparent_410_ppds.set_title('410 km Ppds')
		apparent_410_ppds.set_ylim(500,300)
		apparent_410_ppds.set_ylabel('Depth (km)')
		apparent_410_ppds.yaxis.set_ticks_position('both')
		apparent_410_ppds.yaxis.set_major_locator(MultipleLocator(50))
		apparent_410_ppds.yaxis.set_minor_locator(MultipleLocator(10))
		apparent_410_ppds.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_410_ppds.tick_params(labelleft=True,labelright=False)
		apparent_410_ppds.set_xticks([])


	#### Figure True and Apparent  660 km Ppds  ####


	for _i, _j in enumerate(RF_DEPTH_mean_2_profile_Pds[i]):
		apparent_660_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[i][_i],marker='p',markerfacecolor='none',markeredgecolor='dimgray')
		apparent_660_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile[i][_i],marker='s',markerfacecolor='none',markeredgecolor='k')

		apparent_660_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[i][_i], yerr=RF_DEPTH_std_2_profile_Ppds[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)
		apparent_660_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile[i][_i], yerr=RF_DEPTH_std_2_true_profile[i][_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)

		apparent_660_ppds.set_title('660 km Ppds')
		apparent_660_ppds.set_ylim(800,600)
		apparent_660_ppds.set_ylabel('Depth (km)')
		apparent_660_ppds.yaxis.set_ticks_position('both')
		apparent_660_ppds.yaxis.set_major_locator(MultipleLocator(50))
		apparent_660_ppds.yaxis.set_minor_locator(MultipleLocator(10))
		apparent_660_ppds.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
		apparent_660_ppds.tick_params(labelleft=True,labelright=False)
		apparent_660_ppds.set_xticks([])
		if CROSS_SECTION_AXIS == 'y':
			apparent_660_ppds.text(_i/factor_Ppds,820,"{0:.1f}".format(AB_lon[i][_i]),rotation=-45,fontsize=10)
			apparent_660_ppds.set_xlabel('Longitude ($^\circ$)',labelpad=30)
		else:
			apparent_660_ppds.text(_i/factor_Ppds,820,"{0:.1f}".format(AB_lat[i][_i]),rotation=-45,fontsize=10)
			apparent_660_ppds.set_xlabel('Latitude ($^\circ$)',labelpad=30)


		#### Figure MTZ True and Apparent thickness  ####

	for _i, _j in enumerate(RF_DEPTH_mtz_thickness_profile_Pds[i]):
		MTZ_thickness.plot(_i/factor_Ppds,RF_DEPTH_mtz_thickness_profile_Pds[i][_i],marker='o',markerfacecolor='none',markeredgecolor='gray')
	for _i, _j in enumerate(RF_DEPTH_mtz_thickness_profile_Ppds[i]):
		MTZ_thickness.plot(_i/factor_Ppds,RF_DEPTH_mtz_thickness_profile_Ppds[i][_i],marker='p',markerfacecolor='none',markeredgecolor='dimgray')

	for _i, _j in enumerate(RF_DEPTH_true_thickness_MTZ_profile[i]):
		MTZ_thickness.plot(_i/factor_Ppds,RF_DEPTH_true_thickness_MTZ_profile[i][_i],marker='s',markerfacecolor='none',markeredgecolor='k')

	for _i, _j in enumerate(RF_DEPTH_mtz_thickness_profile_Pds[i]):
		MTZ_thickness.errorbar(_i/factor_Ppds,RF_DEPTH_mtz_thickness_profile_Pds[i][_i], yerr=RF_DEPTH_mtz_thickness_profile_Pds_std[i][_i], ecolor='gray',elinewidth=1,capsize=2,capthick=1)

	for _i, _j in enumerate(RF_DEPTH_mtz_thickness_profile_Ppds[i]):
		MTZ_thickness.errorbar(_i/factor_Ppds,RF_DEPTH_mtz_thickness_profile_Ppds[i][_i], yerr=RF_DEPTH_mtz_thickness_profile_Ppds_std[i][_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)

	for _i, _j in enumerate(RF_DEPTH_true_thickness_MTZ_profile[i]):
		MTZ_thickness.errorbar(_i/factor_Ppds,RF_DEPTH_true_thickness_MTZ_profile[i][_i], yerr=RF_DEPTH_true_thickness_MTZ_profile_std[i][_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	
	MTZ_thickness.axhline(y=250,linewidth=0.25,linestyle='--', color='gray')
	MTZ_thickness.set_ylim(300,200)
	MTZ_thickness.set_ylabel('Depth (km)')
	MTZ_thickness.yaxis.set_label_position("right")
	MTZ_thickness.set_title('MTZ Apparent and True Thickness')
	MTZ_thickness.yaxis.set_ticks_position('both')
	MTZ_thickness.yaxis.set_major_locator(MultipleLocator(20))
	MTZ_thickness.yaxis.set_minor_locator(MultipleLocator(10))
	MTZ_thickness.grid(True,which='major',color='gray',linewidth=1,linestyle='--')
	MTZ_thickness.tick_params(labelleft=True,labelright=True)

	fig.savefig(RESULTS_FOLDER+'SELECTED_BINNED_DATA_'+CROSS_SECTION_AXIS+'_CROSS_SECTION_Pds_Ppds_PROFILE_'+str(i+1)+'.'+EXT_FIG,dpi=DPI_FIG)
print('Ending the Cross-section CODE')