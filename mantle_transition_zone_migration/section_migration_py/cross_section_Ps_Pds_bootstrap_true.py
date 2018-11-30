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
from matplotlib.patches import Circle,Rectangle
import math



from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,GRID_PP_MULT,COLORMAP_STD,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,					
					PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP_MED,DIST_GRID_PP,DEPTH_RANGE,COLORMAP_VEL,
					NUMBER_PP_PER_BIN,NUMBER_STA_PER_BIN
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

RF_stacking_Pds_BOOTSTRAP = SELECTED_BINNED_DATA_dic['data_BOOTSTRAP_Pds']
RF_stacking_Ppds_BOOTSTRAP = SELECTED_BINNED_DATA_dic['data_BOOTSTRAP_Ppds']

RF_BOOTSTRAP_DEPTH_mean_1_Pds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_1_Pds']
RF_BOOTSTRAP_DEPTH_mean_1_Ppds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_1_Ppds']

RF_BOOTSTRAP_DEPTH_mean_2_Pds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_2_Pds']
RF_BOOTSTRAP_DEPTH_mean_2_Ppds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_2_Ppds']

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

print('Calculating earth model layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

print('Plotting ...')
print('\n')

#Color Maps

colormap = plt.get_cmap(COLORMAP_VEL)

fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=PROJECT_LON, globe=None)},figsize=(10,5))

ax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

ax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE,4), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE,4), crs=ccrs.PlateCarree())
ax.tick_params(labelbottom='off',labeltop='on',labelleft='on',labelright='on')

ax.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

norm_410 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)
colors_410 = colormap(norm_410(RF_DEPTH_true_thickness_MTZ_Pds))

for i,j in enumerate(lons):
	retangulo = Rectangle(xy=(lons[i] - DIST_GRID_PP_MED/(GRID_PP_MULT/2), lats[i] - DIST_GRID_PP_MED/(GRID_PP_MULT/2)),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color='None', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=2)
	ax.add_patch(retangulo)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_true_thickness_MTZ_Pds[i]) == False:
		retangulo_410 = Rectangle(xy=(lons[i] - DIST_GRID_PP_MED/(GRID_PP_MULT/2), lats[i] - DIST_GRID_PP_MED/(GRID_PP_MULT/2)),width=DIST_GRID_PP_MED/(GRID_PP_MULT/2), height=DIST_GRID_PP_MED/(GRID_PP_MULT/2),color=colors_410[i], ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=3)
		ax.add_patch(retangulo_410)
	else:
		pass

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)


plt.title('Pick four points for bootstrapping cross-section and them close the windows', y=1.08)


lon_click = []
lat_click = []

def onclick(event):
	dist = []

	for i,j in enumerate(lons):
		latx,laty = ax.projection.transform_point(lons[i],lats[i],src_crs=ccrs.PlateCarree())
		dist.append(np.sqrt((event.xdata-latx)**2+(event.ydata-laty)**2))
	ind= np.argmin(dist)
	print('Cross section lon and lat selected')
	print('lon=%f, lat=%f' %(lons[ind],lats[ind]))
	lon_click.append(lons[ind])
	lat_click.append(lats[ind])
	ax.plot(lons[ind],lats[ind],'Xr',ms=10, transform=ccrs.PlateCarree(),zorder=20)
	plt.draw()
	if len(lon_click) == 4:
		fig.canvas.mpl_disconnect(cid)

	return lon_click,lat_click

cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()

RESULTS_FOLDER = PP_FIGURE+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)
fig.savefig(RESULTS_FOLDER+'SELECTED_BINNED_DATA_CROSS_SECTION_Pds_Ppds_bootstrap_points.'+EXT_FIG,dpi=DPI_FIG)

print('Calculating the distance between selected points and selected grid')
print('\n')

RF_data_profile_Pds = []
RF_DEPTH_mean_1_profile_Pds = []
RF_DEPTH_std_1_profile_Pds = []
RF_DEPTH_mean_2_profile_Pds = []
RF_DEPTH_std_2_profile_Pds = []
RF_DEPTH_mean_1_true_profile_Pds = []
RF_DEPTH_std_1_true_profile_Pds = []
RF_DEPTH_mean_2_true_profile_Pds = []
RF_DEPTH_std_2_true_profile_Pds = []
RF_DEPTH_mtz_thickness_profile_Pds = [] 
RF_DEPTH_true_thickness_MTZ_profile_Pds = [] 

RF_delta_1_Vp_mean_profile = []
RF_delta_1_Vp_std_profile = []
RF_delta_1_Vs_mean_profile = []
RF_delta_1_Vs_std_profile = []

RF_delta_2_Vp_mean_profile = []
RF_delta_2_Vp_std_profile = []
RF_delta_2_Vs_mean_profile = []
RF_delta_2_Vs_std_profile = []


RF_data_profile_Ppds = []
RF_DEPTH_mean_1_profile_Ppds = []
RF_DEPTH_std_1_profile_Ppds = []
RF_DEPTH_mean_2_profile_Ppds = []
RF_DEPTH_std_2_profile_Ppds = []
RF_DEPTH_mean_1_true_profile_Ppds = []
RF_DEPTH_std_1_true_profile_Ppds = []
RF_DEPTH_mean_2_true_profile_Ppds = []
RF_DEPTH_std_2_true_profile_Ppds = []
RF_DEPTH_mtz_thickness_profile_Ppds = [] 
RF_DEPTH_true_thickness_MTZ_profile_Ppds = [] 

RF_stacking_Pds_BOOTSTRAP_profile = []
RF_stacking_Ppds_BOOTSTRAP_profile = []

RF_lat_profile = []
RF_lon_profile = []

RF_BOOTSTRAP_DEPTH_mean_1_Pds_profile = []
RF_BOOTSTRAP_DEPTH_mean_1_Ppds_profile = []

RF_BOOTSTRAP_DEPTH_mean_2_Pds_profile = []
RF_BOOTSTRAP_DEPTH_mean_2_Ppds_profile = []


for i,j in enumerate(lon_click):
		dist = [np.sqrt((j - lons[k])**2 + (lat_click[i] - l)**2)  for k,l in enumerate(lats)]

		RF_lat_profile.append(lats[dist.index(min(dist))])
		RF_lon_profile.append(lons[dist.index(min(dist))])
				
		RF_BOOTSTRAP_DEPTH_mean_1_Pds_profile.append(RF_BOOTSTRAP_DEPTH_mean_1_Pds[dist.index(min(dist))])
		RF_BOOTSTRAP_DEPTH_mean_1_Ppds_profile.append(RF_BOOTSTRAP_DEPTH_mean_1_Ppds[dist.index(min(dist))])

		RF_BOOTSTRAP_DEPTH_mean_2_Pds_profile.append(RF_BOOTSTRAP_DEPTH_mean_2_Pds[dist.index(min(dist))])
		RF_BOOTSTRAP_DEPTH_mean_2_Ppds_profile.append(RF_BOOTSTRAP_DEPTH_mean_2_Ppds[dist.index(min(dist))])

		RF_data_profile_Pds.append(RF_stacking_Pds[dist.index(min(dist))])
		RF_DEPTH_mean_1_profile_Pds.append(RF_DEPTH_mean_1_Pds[dist.index(min(dist))])
		RF_DEPTH_std_1_profile_Pds.append(RF_DEPTH_std_1_Pds[dist.index(min(dist))])
		RF_DEPTH_mean_2_profile_Pds.append(RF_DEPTH_mean_2_Pds[dist.index(min(dist))])
		RF_DEPTH_std_2_profile_Pds.append(RF_DEPTH_std_2_Pds[dist.index(min(dist))])
		RF_DEPTH_mean_1_true_profile_Pds.append(RF_DEPTH_mean_1_true_Pds[dist.index(min(dist))])
		RF_DEPTH_std_1_true_profile_Pds.append(RF_DEPTH_std_1_true_Pds[dist.index(min(dist))])
		RF_DEPTH_mean_2_true_profile_Pds.append(RF_DEPTH_mean_2_true_Pds[dist.index(min(dist))])
		RF_DEPTH_std_2_true_profile_Pds.append(RF_DEPTH_std_2_true_Pds[dist.index(min(dist))])
		RF_DEPTH_mtz_thickness_profile_Pds.append(RF_DEPTH_mtz_thickness_Pds[dist.index(min(dist))])
		RF_DEPTH_true_thickness_MTZ_profile_Pds.append(RF_DEPTH_true_thickness_MTZ_Pds[dist.index(min(dist))]) 

		RF_delta_1_Vp_mean_profile.append(RF_delta_1_Vp_mean[dist.index(min(dist))])
		RF_delta_1_Vp_std_profile.append(RF_delta_1_Vp_std[dist.index(min(dist))])
		RF_delta_1_Vs_mean_profile.append(RF_delta_1_Vs_mean[dist.index(min(dist))])
		RF_delta_1_Vs_std_profile.append(RF_delta_1_Vs_std[dist.index(min(dist))])

		RF_delta_2_Vp_mean_profile.append(RF_delta_2_Vp_mean[dist.index(min(dist))])
		RF_delta_2_Vp_std_profile.append(RF_delta_2_Vp_std[dist.index(min(dist))])
		RF_delta_2_Vs_mean_profile.append(RF_delta_2_Vs_mean[dist.index(min(dist))])
		RF_delta_2_Vs_std_profile.append(RF_delta_2_Vs_std[dist.index(min(dist))])



		RF_data_profile_Ppds.append(RF_stacking_Ppds[dist.index(min(dist))])
		RF_DEPTH_mean_1_profile_Ppds.append(RF_DEPTH_mean_1_Ppds[dist.index(min(dist))])
		RF_DEPTH_std_1_profile_Ppds.append(RF_DEPTH_std_1_Ppds[dist.index(min(dist))])
		RF_DEPTH_mean_2_profile_Ppds.append(RF_DEPTH_mean_2_Ppds[dist.index(min(dist))])
		RF_DEPTH_std_2_profile_Ppds.append(RF_DEPTH_std_2_Ppds[dist.index(min(dist))])
		RF_DEPTH_mean_1_true_profile_Ppds.append(RF_DEPTH_mean_1_true_Ppds[dist.index(min(dist))])
		RF_DEPTH_std_1_true_profile_Ppds.append(RF_DEPTH_std_1_true_Ppds[dist.index(min(dist))])
		RF_DEPTH_mean_2_true_profile_Ppds.append(RF_DEPTH_mean_2_true_Ppds[dist.index(min(dist))])
		RF_DEPTH_std_2_true_profile_Ppds.append(RF_DEPTH_std_2_true_Ppds[dist.index(min(dist))])
		RF_DEPTH_mtz_thickness_profile_Ppds.append(RF_DEPTH_mtz_thickness_Ppds[dist.index(min(dist))])
		RF_DEPTH_true_thickness_MTZ_profile_Ppds.append(RF_DEPTH_true_thickness_MTZ_Ppds[dist.index(min(dist))]) 

		RF_stacking_Pds_BOOTSTRAP_profile.append(RF_stacking_Pds_BOOTSTRAP[dist.index(min(dist))]) 
		RF_stacking_Ppds_BOOTSTRAP_profile.append(RF_stacking_Ppds_BOOTSTRAP[dist.index(min(dist))]) 

print('Plotting the Final Figure')

#Cross section figure

fig = plt.figure(figsize=(30, 15))

fig.suptitle('Pds and Ppds Bootstrapping')


gs = gridspec.GridSpec(6, 8)
gs.update(wspace=0.5, hspace=0.75)


#Figure Pds
for _i, _j in enumerate(RF_data_profile_Pds):

		pds_grid = fig.add_subplot(gs[0:3, _i*2:_i*2+1])
		ppds_grid = fig.add_subplot(gs[3:6, _i*2:_i*2+1])

		pds_grid_410_660 = fig.add_subplot(gs[0:3, _i*2+1])

		ppds_grid_410_660 = fig.add_subplot(gs[3:6,_i*2+1])

		

		factor_Pds = 1

		majorLocatorY = MultipleLocator(50)
		minorLocatorY = MultipleLocator(25)

	
		x_data_Pds= []
		for x,c in enumerate(RF_stacking_Pds_BOOTSTRAP_profile[_i]):
			RF_data_factor_Pds_bootstrap = [_i/factor_Pds+l for k, l in enumerate(c)]
			x_data_Pds.append(RF_data_factor_Pds_bootstrap)
			pds_grid.plot(RF_data_factor_Pds_bootstrap,camadas_terra_10_km,'silver',linewidth=0.1, zorder=10)

		min_x = [min(a) for a in zip(*x_data_Pds)]
		max_x = [max(a) for a in zip(*x_data_Pds)]
		pds_grid.fill_betweenx(y=camadas_terra_10_km,x1=min_x, x2=max_x, facecolor='whitesmoke',alpha=0.8, interpolate=True, zorder=5)

		pds_grid.text(min(min_x),RF_DEPTH_mean_1_true_profile_Pds[_i],str(round(RF_DEPTH_mean_1_true_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_std_1_true_profile_Pds[_i])),zorder=40, weight = 'bold',fontsize='x-small')
		pds_grid.text(min(min_x),RF_DEPTH_mean_2_true_profile_Pds[_i],str(round(RF_DEPTH_mean_2_true_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_std_2_true_profile_Pds[_i])),zorder=41, weight = 'bold',fontsize='x-small')


		RF_data_factor_Pds = [_i/factor_Pds+l for k, l in enumerate(_j)]
		pds_grid.plot(RF_data_factor_Pds,camadas_terra_10_km,'k',linewidth=2, zorder=30)

		pds_grid.yaxis.set_ticks_position('both')
		pds_grid.yaxis.set_major_locator(majorLocatorY)
		pds_grid.grid(True,which='major',linestyle='--')

		pds_grid.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)>=_i/factor_Pds, facecolor='dimgrey',interpolate=True, zorder=19)
		pds_grid.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)<=_i/factor_Pds, facecolor='lightgrey', interpolate=True, zorder=20)
		pds_grid.set_xticks([])

	
		pds_grid.set_title('Lat = '+str(round(RF_lat_profile[_i],1))+' - Lon = '+str(round(RF_lon_profile[_i],1)))
		pds_grid.set_ylim(800,300)

		if _i == 0:
			pds_grid.set_ylabel('Depth (km)')
			pds_grid.yaxis.set_label_position("left")
		if _i != 0:
			pds_grid.axes.axes.yaxis.set_ticklabels([])

	#Figure Ppds


		factor_Ppds = 1

		x_data_Ppds = []
		for x,c in enumerate(RF_stacking_Ppds_BOOTSTRAP_profile[_i]):
			RF_data_factor_Ppds_bootstrap = [_i/factor_Ppds+l for k, l in enumerate(c)]
			x_data_Ppds.append(RF_data_factor_Ppds_bootstrap)
			ppds_grid.plot(RF_data_factor_Ppds_bootstrap,camadas_terra_10_km,'silver',linewidth=0.1, zorder=10)

		min_x = [min(a) for a in zip(*x_data_Ppds)]
		max_x = [max(a) for a in zip(*x_data_Ppds)]
		ppds_grid.fill_betweenx(y=camadas_terra_10_km,x1=min_x, x2=max_x, facecolor='whitesmoke',alpha=0.8, interpolate=True, zorder=5)

		RF_data_factor_Ppds = [_i/factor_Ppds+l for k, l in enumerate(RF_data_profile_Ppds[_i])]
		ppds_grid.plot(RF_data_factor_Ppds,camadas_terra_10_km,'k',linewidth=2, zorder=30)

		ppds_grid.text(min(min_x),RF_DEPTH_mean_1_true_profile_Ppds[_i],str(round(RF_DEPTH_mean_1_true_profile_Ppds[_i]))+'±'+str(round(RF_DEPTH_std_1_true_profile_Ppds[_i])),zorder=40, weight = 'bold',fontsize='x-small')
		ppds_grid.text(min(min_x),RF_DEPTH_mean_2_true_profile_Ppds[_i],str(round(RF_DEPTH_mean_2_true_profile_Ppds[_i]))+'±'+str(round(RF_DEPTH_std_2_true_profile_Ppds[_i])),zorder=41, weight = 'bold',fontsize='x-small')


		ppds_grid.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)>=_i/factor_Ppds, facecolor='dimgrey',interpolate=True, zorder=19)
		ppds_grid.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)<=_i/factor_Ppds, facecolor='lightgrey',interpolate=True, zorder=20)
		
		ppds_grid.yaxis.set_major_locator(majorLocatorY)
		ppds_grid.grid(True,which='major',linestyle='--')

		ppds_grid.yaxis.set_ticks_position('both')

		#ppds_grid.set_title('Lat = '+str(round(RF_lat_profile[_i],1))+' - Lon = '+str(round(RF_lon_profile[_i],1)))
		ppds_grid.set_xticks([])
		ppds_grid.set_ylim(800,300)

		if _i == 0:
			ppds_grid.set_ylabel('Depth (km)')
			ppds_grid.yaxis.set_label_position("left")

		if _i != 0:
			ppds_grid.axes.axes.yaxis.set_ticklabels([])


		#### Plot Depth 410 Pds ####

		pds_grid_410_660.hist(RF_BOOTSTRAP_DEPTH_mean_1_Pds_profile[_i],bins=5,orientation='horizontal',color='k')

		#### Plot Depth 660 Pds ####

		pds_grid_410_660.hist(RF_BOOTSTRAP_DEPTH_mean_2_Pds_profile[_i],bins=5,orientation='horizontal',color='k')

		pds_grid_410_660.yaxis.set_ticks_position('both')
		pds_grid_410_660.yaxis.set_ticks_position('both')
		pds_grid_410_660.yaxis.set_major_locator(MultipleLocator(50))
		pds_grid_410_660.grid(True,which='major',linestyle='--')
		pds_grid_410_660.set_xlim(0,100)
		pds_grid_410_660.set_ylim(800,300)
		pds_grid_410_660.axes.axes.xaxis.set_ticklabels([])


		if _i != 3:
			pds_grid_410_660.axes.axes.yaxis.set_ticklabels([])

		if _i == 3:
			pds_grid_410_660.set_ylabel('Depth (km)')
			pds_grid_410_660.yaxis.set_label_position("right")
			pds_grid_410_660.tick_params(labelright=True,labelleft=False)


		#### Plot Depth 410 Ppds ####

		ppds_grid_410_660.hist(RF_BOOTSTRAP_DEPTH_mean_1_Ppds_profile[_i],bins=5,orientation='horizontal',color='k')

		#### Plot Depth 660 Ppds ####

		ppds_grid_410_660.hist(RF_BOOTSTRAP_DEPTH_mean_2_Ppds_profile[_i],bins=5,orientation='horizontal',color='k')
		ppds_grid_410_660.yaxis.set_ticks_position('both')
		ppds_grid_410_660.set_xlim(0,100)
		ppds_grid_410_660.yaxis.set_ticks_position('both')
		ppds_grid_410_660.yaxis.set_major_locator(MultipleLocator(50))
		ppds_grid_410_660.grid(True,which='major',linestyle='--')
		ppds_grid_410_660.set_xlabel('Population')
		ppds_grid_410_660.set_ylim(800,300)


		if _i != 3:
			ppds_grid_410_660.axes.axes.yaxis.set_ticklabels([])

		if _i == 3:
			ppds_grid_410_660.set_ylabel('Depth (km)')
			ppds_grid_410_660.yaxis.set_label_position("right")
			ppds_grid_410_660.tick_params(labelright=True,labelleft=False)
		

plt.show()

fig.savefig(RESULTS_FOLDER+'SELECTED_BINNED_DATA_CROSS_SECTION_Pds_Ppds_bootstrap.'+EXT_FIG,dpi=DPI_FIG)

print('Ending the Cross section CODE')