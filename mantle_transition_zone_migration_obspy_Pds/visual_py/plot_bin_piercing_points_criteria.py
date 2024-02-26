# coding: utf-8

import matplotlib.pyplot as plt
import matplotlib as mpl

import numpy as np
import obspy
import os
import glob
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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
import pandas as pd
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import interpolate
from shapely.geometry import Polygon, MultiPoint, Point, LinearRing
import shapefile
from matplotlib.colors import Normalize
from matplotlib.patches import Circle,Rectangle
import math
import operator
from matplotlib.transforms import offset_copy






from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					NUMBER_PP_PER_BIN,LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,OUTPUT_DIR,
					EXT_FIG,DPI_FIG,DIST_GRID_PP,NUMBER_STA_PER_BIN,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,COLORMAP_STD,COLORMAP_VEL,DEPTH_TARGET
				   )


print('Starting Check resolution criteria CODE')
print('\n')

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Stations'+'/'

print('Looking for receiver functions data in JSON file in '+STA_DIR)
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

print('Looking for selected binned data')
print('\n')

PP_SELEC_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'SELECTED_BINNED_DATA'+'/'


lst_json_file = []
for root, dirs, files in os.walk(PP_SELEC_DIR):
    for file in files:
        if file.endswith(".feather"):
             lst_json_file.append(os.path.join(root, file))

lst_json_file_PP = [int(i.split('NUMBER_PP_PER_BIN_')[1].split('_NUMBER_STA_PER_BIN_')[0]) for i in lst_json_file ]
lst_json_file_STA =  [int(i.split('_NUMBER_STA_PER_BIN_')[1].split('/SELECTED_BINNED.feather')[0]) for i in lst_json_file ]

ind = np.lexsort((lst_json_file_PP,lst_json_file_STA))

PP_index =  [lst_json_file_PP[i] for i in ind]
STA_index =  [lst_json_file_STA[i] for i in ind]

sort_lst_json = [lst_json_file[i] for i in ind]

print('Importing selected binned data')
print('\n')

lats = []
lons =  []

RF_DEPTH_mean_1_Pds =  []
RF_DEPTH_std_1_Pds =  []

RF_DEPTH_mean_2_Pds =  []
RF_DEPTH_std_2_Pds =  []

RF_DEPTH_mtz_thickness_Pds =  []
RF_DEPTH_mtz_thickness_Pds_std =  []

for i,j in enumerate(sort_lst_json):

	SELECTED_BINNED_DATA_dic = pd.read_feather(j)

	lats.append(SELECTED_BINNED_DATA_dic['lat'].tolist())
	lons.append(SELECTED_BINNED_DATA_dic['lon'].tolist())

	RF_DEPTH_mean_1_Pds.append(SELECTED_BINNED_DATA_dic['mean_1_Pds'].tolist())
	RF_DEPTH_std_1_Pds.append(SELECTED_BINNED_DATA_dic['std_1_Pds'].tolist())

	RF_DEPTH_mean_2_Pds.append(SELECTED_BINNED_DATA_dic['mean_2_Pds'].tolist())
	RF_DEPTH_std_2_Pds.append(SELECTED_BINNED_DATA_dic['std_2_Pds'].tolist())

	RF_DEPTH_mtz_thickness_Pds.append(SELECTED_BINNED_DATA_dic['mtz_thickness_Pds'].tolist())
	RF_DEPTH_mtz_thickness_Pds_std.append(SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std'].tolist())


RESULTS_FOLDER = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Figures'+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)

###################################################################################################################

#Color Maps

colormap = plt.get_cmap(COLORMAP_VEL)


colormap_std = plt.get_cmap(COLORMAP_STD)


#############################################################################################################################################################################################


mosaic_lst_MTZ = RF_DEPTH_mtz_thickness_Pds
mosaic_lst_MTZ_name = 'MTZ_thickness_Pds'
mosaic_lst_MTZ_label = 'MTZ thickness Pds (km)'

mosaic_lst_660 = RF_DEPTH_mean_2_Pds
mosaic_lst_660_name = 'DEPTH_660_Pds'
mosaic_lst_660_label = '660 depth Pds (km)'

mosaic_lst_410 = RF_DEPTH_mean_1_Pds
mosaic_lst_410_name = 'DEPTH_410_Pds'
mosaic_lst_410_label = '410 depth Pds (km)'


#############################################################################################################################################################################################

def plot_mosaic_MTZ(mosaic_lst,mosaic_lst_name,mosaic_lst_label):

		fig, axes = plt.subplots(nrows=len(set(lst_json_file_PP)), ncols=len(set(lst_json_file_STA)), gridspec_kw = {'wspace':0, 'hspace':0}, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(1+len(set(lst_json_file_PP))*2,len(set(lst_json_file_STA))*2),sharex='col', sharey='row',constrained_layout=True)

		for k,ax in zip(range(len(mosaic_lst)), axes.flat):

			ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

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
			ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=0.1, color='k', alpha=0.5, linestyle='dotted')
			
			geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
			text_transform = offset_copy(geodetic_transform, units='dots', x=0)

			ax.text(-40.5, -1.8,'p = '+str(PP_index[k]),fontsize=7, verticalalignment='center', 
						horizontalalignment='right', transform=text_transform,bbox=dict(facecolor='white',edgecolor='white',pad=0.1))

			ax.text(-40.5, -0.6,'s = '+str(STA_index[k]),fontsize=7, verticalalignment='center', 
						horizontalalignment='right', transform=text_transform,bbox=dict(facecolor='white',edgecolor='white',pad=0.1))

			n=30
			x = 0.5
			lower = colormap(np.linspace(0, x, n))
			white = colormap(np.ones(80-2*n)*x)
			upper = colormap(np.linspace(1-x, 1, n))
			colors = np.vstack((lower, white, upper))
			tmap = mpl.colors.LinearSegmentedColormap.from_list('map_white', colors)

			norm_map_MTZ = Normalize(vmin=200,vmax=300)

			for i,j in enumerate(lons[k]):
				if math.isnan(mosaic_lst[k][i]) == False:
					retangulo_660 = Circle(radius=DIST_GRID_PP,xy=(lons[k][i], lats[k][i]),color=tmap(norm_map_MTZ(mosaic_lst[k][i])), ec='k',linewidth=0.2,transform=ccrs.Geodetic(),zorder=2)
					ax.add_patch(retangulo_660)
				else: 
					pass
		#______________________________________________________________________

		#sm_MTZ = plt.cm.ScalarMappable(cmap=tmap,norm=norm_map_MTZ)
		#sm_MTZ._A = []

		#fig.colorbar(sm_MTZ,  ax=axes.flat, orientation='horizontal',shrink=0.5,fraction=0.05,label=mosaic_lst_label)

		fig.savefig(RESULTS_FOLDER+mosaic_lst_name+'_mosaic.'+EXT_FIG,dpi=DPI_FIG)


def plot_mosaic_660(mosaic_lst,mosaic_lst_name,mosaic_lst_label):

		fig, axes = plt.subplots(nrows=len(set(lst_json_file_PP)), ncols=len(set(lst_json_file_STA)), gridspec_kw = {'wspace':0, 'hspace':0}, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(1+len(set(lst_json_file_PP))*2,len(set(lst_json_file_STA))*2),sharex='col', sharey='row',constrained_layout=True)

		for k,ax in zip(range(len(mosaic_lst)), axes.flat):

			ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

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
			ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=0.1, color='k', alpha=0.5, linestyle='dotted')
			
			geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
			text_transform = offset_copy(geodetic_transform, units='dots', x=0)

			ax.text(-40.5, -1.8,'p = '+str(PP_index[k]),fontsize=7, verticalalignment='center', 
						horizontalalignment='right', transform=text_transform,bbox=dict(facecolor='white',edgecolor='white',pad=0.1))

			ax.text(-40.5, -0.6,'s = '+str(STA_index[k]),fontsize=7, verticalalignment='center', 
						horizontalalignment='right', transform=text_transform,bbox=dict(facecolor='white',edgecolor='white',pad=0.1))


			n=30
			x = 0.5
			lower = colormap(np.linspace(0, x, n))
			white = colormap(np.ones(80-2*n)*x)
			upper = colormap(np.linspace(1-x, 1, n))
			colors = np.vstack((lower, white, upper))
			tmap = mpl.colors.LinearSegmentedColormap.from_list('map_white', colors)

			norm_660 = Normalize(vmin=610,vmax=710)

			for i,j in enumerate(lons[k]):
				if math.isnan(mosaic_lst[k][i]) == False:
					retangulo_660 = Circle(radius=DIST_GRID_PP,xy=(lons[k][i],lats[k][i]),color=tmap(norm_660(mosaic_lst[k][i])), ec='k',linewidth=0.2,transform=ccrs.Geodetic(),zorder=2)
					ax.add_patch(retangulo_660)
				else: 
					pass
		#______________________________________________________________________

		#sm_660 = plt.cm.ScalarMappable(cmap=tmap,norm=norm_660)
		#sm_660._A = []
		
		#fig.colorbar(sm_660,  ax=axes.flat, orientation='horizontal',shrink=0.5,fraction=0.05,label=mosaic_lst_label)

		fig.savefig(RESULTS_FOLDER+mosaic_lst_name+'_mosaic.'+EXT_FIG,dpi=DPI_FIG)


def plot_mosaic_410(mosaic_lst,mosaic_lst_name,mosaic_lst_label):

		fig, axes = plt.subplots(nrows=len(set(lst_json_file_PP)), ncols=len(set(lst_json_file_STA)), gridspec_kw = {'wspace':0, 'hspace':0}, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(1+len(set(lst_json_file_PP))*2,len(set(lst_json_file_STA))*2),sharex='col', sharey='row',constrained_layout=True)

		for k,ax in zip(range(len(mosaic_lst)), axes.flat):

			ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])

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
			ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=0.1, color='k', alpha=0.5, linestyle='dotted')
			
			geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
			text_transform = offset_copy(geodetic_transform, units='dots', x=0)

			ax.text(-40.5, -1.8,'p = '+str(PP_index[k]),fontsize=7, verticalalignment='center', 
						horizontalalignment='right', transform=text_transform,bbox=dict(facecolor='white',edgecolor='white',pad=0.1))

			ax.text(-40.5, -0.6,'s = '+str(STA_index[k]),fontsize=7, verticalalignment='center', 
						horizontalalignment='right', transform=text_transform,bbox=dict(facecolor='white',edgecolor='white',pad=0.1))

			n=30
			x = 0.5
			lower = colormap(np.linspace(0, x, n))
			white = colormap(np.ones(80-2*n)*x)
			upper = colormap(np.linspace(1-x, 1, n))
			colors = np.vstack((lower, white, upper))
			tmap = mpl.colors.LinearSegmentedColormap.from_list('map_white', colors)

			norm_410 = Normalize(vmin=360,vmax=460)

			for i,j in enumerate(lons[k]):
				if math.isnan(mosaic_lst[k][i]) == False:
					retangulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[k][i], lats[k][i]),color=tmap(norm_410(mosaic_lst[k][i])), ec='k',linewidth=0.2,transform=ccrs.Geodetic(),zorder=2)
					ax.add_patch(retangulo_410)
				else: 
					pass
		#______________________________________________________________________

		#sm_410 = plt.cm.ScalarMappable(cmap=tmap,norm=norm_410)
		#sm_410._A = []

		#fig.colorbar(sm_410,  ax=axes.flat, orientation='horizontal',shrink=0.5,fraction=0.05,label=mosaic_lst_label)


		fig.savefig(RESULTS_FOLDER+mosaic_lst_name+'_mosaic.'+EXT_FIG,dpi=DPI_FIG)

########################################################################################################################################################################
print('Plotting Figure: Mosaic of each estimates')

plot_mosaic_MTZ(mosaic_lst_MTZ,mosaic_lst_MTZ_name,mosaic_lst_MTZ_label)
plot_mosaic_660(mosaic_lst_660,mosaic_lst_660_name,mosaic_lst_660_label)
plot_mosaic_410(mosaic_lst_410,mosaic_lst_410_name,mosaic_lst_410_label)

print('Ending Final Plot CODE')