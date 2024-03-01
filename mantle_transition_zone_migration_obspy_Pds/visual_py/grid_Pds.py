# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import copy
import matplotlib
from matplotlib.cm import get_cmap
from mpl_toolkits.mplot3d import Axes3D
import shapefile
import scipy.io
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import json
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import interpolate
from mayavi import mlab

from parameters_py.mgconfig import (
					RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,
					NUMBER_PP_PER_BIN,GRID_PP_MULT,OUTPUT_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,NUMBER_STA_PER_BIN,			
					EXT_FIG,DPI_FIG,DEPTH_RANGE,DEPTH_TARGET
				   )


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


#Lat/Lon lists:

lats = SELECTED_BINNED_DATA_dic['lat'].tolist()
lons = SELECTED_BINNED_DATA_dic['lon'].tolist()


#Receiver Functions:

RF_number = SELECTED_BINNED_DATA_dic['len_Pds'].tolist()
RF_stacking_Pds = SELECTED_BINNED_DATA_dic['data_Pds'].tolist()

#Estimates 350 LVZ:

RF_DEPTH_mean_LVZ_Pds = SELECTED_BINNED_DATA_dic['mean_LVZ_Pds'].tolist()
RF_DEPTH_std_LVZ_Pds = SELECTED_BINNED_DATA_dic['std_LVZ_Pds'].tolist()


#Estimates P410s:

RF_DEPTH_mean_1_Pds = SELECTED_BINNED_DATA_dic['mean_1_Pds'].tolist()
RF_DEPTH_std_1_Pds = SELECTED_BINNED_DATA_dic['std_1_Pds'].tolist()

#Estimates P520s:

RF_DEPTH_mean_520_Pds = SELECTED_BINNED_DATA_dic['mean_520_Pds'].tolist()
RF_DEPTH_std_520_Pds = SELECTED_BINNED_DATA_dic['std_520_Pds'].tolist()

#Estimates P660s:

RF_DEPTH_mean_2_Pds = SELECTED_BINNED_DATA_dic['mean_2_Pds'].tolist()
RF_DEPTH_std_2_Pds = SELECTED_BINNED_DATA_dic['std_2_Pds'].tolist()

#Estimates 700 LVZ:

RF_DEPTH_mean_LVZ_700_Pds = SELECTED_BINNED_DATA_dic['mean_LVZ_700_Pds'].tolist()
RF_DEPTH_std_LVZ_700_Pds = SELECTED_BINNED_DATA_dic['std_LVZ_700_Pds'].tolist()

#Estimates MTZ Pds:

RF_DEPTH_mtz_thickness_Pds = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds'].tolist()
RF_DEPTH_mtz_thickness_Pds_std = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std'].tolist()

#############################################################################################################################

grid_camadas_x = []
grid_camadas_y = []
grid_camadas_z = []
grid_Pds_data = []

grid_camadas_x_410 = []
grid_camadas_y_410 = []
grid_camadas_z_410 = []
grid_Pds_data_410 = []

grid_camadas_x_660 = []
grid_camadas_y_660 = []
grid_camadas_z_660 = []
grid_Pds_data_660 = []


grid_PP_x_660 = []
grid_PP_y_660 = []
grid_PP_z_660 = []
grid_PP_data_660 = []


for i,j in enumerate(camadas_terra_10_km):
	if j in range(380, 450) or j in range(630, 700):
		for k,l in enumerate(lons):
			grid_camadas_x.append(lons[k])
			grid_camadas_y.append(lats[k])
			grid_camadas_z.append(-j)
			grid_Pds_data.append(RF_stacking_Pds[k][i])
	else:
		for k,l in enumerate(lons):
			grid_camadas_x.append(lons[k])
			grid_camadas_y.append(lats[k])
			grid_camadas_z.append(-j)
			grid_Pds_data.append(RF_stacking_Pds[k][i])
			#grid_Pds_data.append(RF_stacking_Pds[k][i]/100)


for i,j in enumerate(camadas_terra_10_km):
	if j == 410:
		for k,l in enumerate(lons):
			grid_camadas_x_410.append(lons[k])
			grid_camadas_y_410.append(lats[k])
			grid_camadas_z_410.append(-j)
			grid_Pds_data_410.append(RF_stacking_Pds[k][i])


for i,j in enumerate(camadas_terra_10_km):
	if j == 660:
		for k,l in enumerate(lons):
			grid_camadas_x_660.append(lons[k])
			grid_camadas_y_660.append(lats[k])
			grid_camadas_z_660.append(-j)
			grid_Pds_data_660.append(RF_stacking_Pds[k][i])

#for i,j in enumerate(pp_2_lat):
	#if j != []:
		#grid_PP_x_660.append(pp_2_long[i])
		#grid_PP_y_660.append(pp_2_lat[i])
		#grid_PP_z_660.append(-660)


################################################################################
print('Plotting: Figure 3D mesh with Mayavi')
print('\n')


xi = np.linspace(LLCRNRLON_SMALL, URCRNRLON_SMALL, abs(abs(URCRNRLON_SMALL) - abs(LLCRNRLON_SMALL))*GRID_PP_MULT) 
yi = np.linspace(LLCRNRLAT_SMALL, URCRNRLAT_SMALL, abs(abs(URCRNRLAT_SMALL) - abs(LLCRNRLAT_SMALL))*GRID_PP_MULT)
zi = np.linspace(-800,-300,50)


gridx, gridy, gridz = np.meshgrid(xi,yi,zi)

volume = interpolate.griddata((np.array(grid_camadas_x),np.array(grid_camadas_y),np.array(grid_camadas_z)), np.array(grid_Pds_data), (gridx, gridy,gridz),method='linear')


mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))

#mlab.pipeline.volume(mlab.pipeline.scalar_field(volume), vmin=-0.003, vmax=0.02)

#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(volume),plane_orientation='z_axes',slice_index=20,colormap="seismic")

#vol = mlab.points3d(np.array(grid_camadas_x),np.array(grid_camadas_y),np.array(grid_camadas_z), np.array(grid_Pds_data),colormap="seismic",scale_factor=0.2)
vol = mlab.contour3d(gridx, gridy, gridz, volume,colormap="seismic")
mlab.colorbar(vol)
mlab.xlabel('Longitude (degree)')
mlab.ylabel('Latitude (degree)')
mlab.zlabel('Depth (km)')
mlab.outline()
mlab.show()

print('Ending the grid 3D CODE')