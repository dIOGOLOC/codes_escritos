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
from fatiando import gridder, utils
import scipy.io
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import json
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import interpolate
from mayavi import mlab
from tvtk.api import tvtk







from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,GRID_PP_MULT,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,SECTION_NUM,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP_MED,DIST_GRID_PP,DEPTH_RANGE
				   )


print('Starting Cross section CODE')
print('\n')

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

print('Creating my earth model')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

dist_med_camada_terra = [abs(c - ((410+660)/2)) for x,c in enumerate(camadas_terra_10_km)]

DEPTH_MED = camadas_terra_10_km[dist_med_camada_terra.index(min(dist_med_camada_terra))]

print('Importing Pds piercing points to each PHASE')
print('\n')

PHASES = 'P410s','P'+"{0:.0f}".format(DEPTH_MED)+'s','P660s'

print('Importing Pds Piercing Points for '+PHASES[0])
print('\n')

filename_1 = PP_DIR+'PP_'+PHASES[0]+'_dic.json'

PP_1_dic = json.load(open(filename_1))

PP_dist_1 = []
PP_time_1 = []
PP_lat_1 = []
PP_lon_1 = []
PP_depth_1 = [] 

for i,j in enumerate(PP_1_dic):
	PP_dist_1.append(j['dist'][0])
	PP_time_1.append(j['time'][0])
	PP_lat_1.append(j['lat'][0])
	PP_lon_1.append(j['lon'][0])
	PP_depth_1.append(j['depth'][0])

print('Importing Pds Piercing Points for '+PHASES[1])
print('\n')

filename_med = PP_DIR+'PP_'+PHASES[1]+'_dic.json'

PP_med_dic = json.load(open(filename_med))

PP_dist_med = []
PP_time_med = []
PP_lat_med = []
PP_lon_med = []
PP_depth_med = [] 

for i,j in enumerate(PP_med_dic):
	PP_dist_med.append(j['dist'][0])
	PP_time_med.append(j['time'][0])
	PP_lat_med.append(j['lat'][0])
	PP_lon_med.append(j['lon'][0])
	PP_depth_med.append(j['depth'][0])

print('Importing Pds Piercing Points for '+PHASES[2])
print('\n')

filename_2 = PP_DIR+'PP_'+PHASES[2]+'_dic.json'

PP_2_dic = json.load(open(filename_2))

PP_dist_2 = []
PP_time_2 = []
PP_lat_2 = []
PP_lon_2 = []
PP_depth_2 = [] 

for i,j in enumerate(PP_2_dic):
	PP_dist_2.append(j['dist'][0])
	PP_time_2.append(j['time'][0])
	PP_lat_2.append(j['lat'][0])
	PP_lon_2.append(j['lon'][0])
	PP_depth_2.append(j['depth'][0])

print('P410s Piercing Points')
print('\n')

pp_1_lat  = [[]]*len(PP_lon_1)
pp_1_long  = [[]]*len(PP_lon_1)


for i,j in enumerate(PP_lon_1):
    for k,l in enumerate(j):
        if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_1[i][k] == 410:
                pp_1_lat[i] = PP_lat_1[i][k] 
                pp_1_long[i] = l


print('Pds Piercing Points - '+"{0:.0f}".format(DEPTH_MED))
print('\n')

pp_med_lat  = [[]]*len(PP_lon_med)
pp_med_long  = [[]]*len(PP_lon_med)


for i,j in enumerate(PP_lon_med):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_med[i][k] == DEPTH_MED:
			pp_med_lat[i] = PP_lat_med[i][k] 
			pp_med_long[i] = l

print('P660s Piercing Points')
print('\n')


pp_2_lat  = [[]]*len(PP_lon_2)
pp_2_long  = [[]]*len(PP_lon_2)


for i,j in enumerate(PP_lon_2):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_2[i][k] == 660:
			pp_2_lat[i] = PP_lat_2[i][k]
			pp_2_long[i] = l


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
			#grid_Pds_data.append(RF_stacking_Ppds[k][i])
			grid_Pds_data.append(RF_stacking_Pds[k][i]/100)


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

for i,j in enumerate(pp_2_lat):
	if j != []:
		grid_PP_x_660.append(pp_2_long[i])
		grid_PP_y_660.append(pp_2_lat[i])
		grid_PP_z_660.append(-660)


################################################################################
print('Plotting: Figure 3D mesh with Mayavi')
print('\n')


xi = np.linspace(LLCRNRLON_SMALL, URCRNRLON_SMALL, abs(abs(URCRNRLON_SMALL) - abs(LLCRNRLON_SMALL))*GRID_PP_MULT) 
yi = np.linspace(LLCRNRLAT_SMALL, URCRNRLAT_SMALL, abs(abs(URCRNRLAT_SMALL) - abs(LLCRNRLAT_SMALL))*GRID_PP_MULT)
zi = np.linspace(-800, -300, 50)

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