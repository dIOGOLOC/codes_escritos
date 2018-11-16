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
from mpl_toolkits.basemap import Basemap
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
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,GRID_PP_MULT,ROTATE_ANGLE,
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

print('Creating my earth model')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

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


################################################################################
print('Plotting: Figure 3D mesh with Mayavi')
print('\n')

xi = np.linspace(LLCRNRLON_SMALL, URCRNRLON_SMALL, 100)
yi = np.linspace(LLCRNRLAT_SMALL, URCRNRLAT_SMALL, 100)
zi = np.linspace(-800, -300, 100)

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

'''
print('Plotting: Figure 3D Grid with Mayavi')
print('\n')

mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))

# Visualize the points
pts_410 = mlab.points3d(grid_camadas_x_410, grid_camadas_y_410, grid_camadas_z_410,grid_Pds_data_410, scale_mode='none', scale_factor=0.2,colormap="seismic")
pts_660 = mlab.points3d(grid_camadas_x_660, grid_camadas_y_660, grid_camadas_z_660,grid_Pds_data_660, scale_mode='none', scale_factor=0.2,colormap="seismic")

# Create and visualize the mesh
mesh_410 = mlab.pipeline.delaunay2d(pts_410)
surf_410 = mlab.pipeline.surface(mesh_410,colormap="seismic")

mesh_660 = mlab.pipeline.delaunay2d(pts_660)
surf_660 = mlab.pipeline.surface(mesh_660,colormap="seismic")

mlab.show()

print('Ending the grid 3D CODE')
'''