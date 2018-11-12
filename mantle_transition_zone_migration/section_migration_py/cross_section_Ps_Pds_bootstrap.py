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
from mpl_toolkits.basemap import Basemap
import shapefile
from fatiando import gridder, utils
import scipy.io
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import json
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable



from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,SECTION_NUM,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP_MED,DIST_GRID_PP,DEPTH_RANGE
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

print('Plotting Total of Receiver Functions per bin')
print('\n')

fig_receiver_function_per_bin,ax = plt.subplots(1,1,figsize=(10,5))

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)



x, y = m(lons,lats)
sc = m.scatter(x,y,40,RF_number,cmap='copper',marker='s',edgecolors='none')
plt.colorbar(sc,label='Total of Receiver Functions per bin')



m.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m.drawcoastlines(color='k',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

plt.title('Pick two points for cross-section and them close the windows', y=1.08)


lon_click = []
lat_click = []
def onclick(event):
	dist = []
	for i,j in enumerate(lons):
		lon_x,lat_y = m(j,lats[i])
		dist.append(np.sqrt((event.xdata-lon_x)**2+(event.ydata-lat_y)**2))
	ind= np.argmin(dist)
	print('Cross section lon and lat selected')
	print('lon=%f, lat=%f' %(lons[ind],lats[ind]))
	lon_click.append(lons[ind])
	lat_click.append(lats[ind])
	ax.plot(event.xdata,event.ydata,'xr',ms=10)
	plt.draw()
	if len(lon_click) == 4:
		fig_receiver_function_per_bin.canvas.mpl_disconnect(cid)

	return lon_click,lat_click

cid = fig_receiver_function_per_bin.canvas.mpl_connect('button_press_event', onclick)

plt.show()

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
RF_DEPTH_mean_2_true_profile_Ppds = []
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
		RF_DEPTH_mean_2_true_profile_Ppds.append(RF_DEPTH_mean_2_true_Ppds[dist.index(min(dist))])
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

		pds_grid.text(min(min_x),RF_DEPTH_mean_1_profile_Pds[_i],str(round(RF_DEPTH_mean_1_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_std_1_profile_Pds[_i])),zorder=40, weight = 'bold',fontsize='x-small')
		pds_grid.text(min(min_x),RF_DEPTH_mean_2_profile_Pds[_i],str(round(RF_DEPTH_mean_2_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_std_2_profile_Pds[_i])),zorder=41, weight = 'bold',fontsize='x-small')


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

		ppds_grid.text(min(min_x),RF_DEPTH_mean_1_profile_Ppds[_i],str(round(RF_DEPTH_mean_1_profile_Ppds[_i]))+'±'+str(round(RF_DEPTH_std_1_profile_Ppds[_i])),zorder=40, weight = 'bold',fontsize='x-small')
		ppds_grid.text(min(min_x),RF_DEPTH_mean_2_profile_Ppds[_i],str(round(RF_DEPTH_mean_2_profile_Ppds[_i]))+'±'+str(round(RF_DEPTH_std_2_profile_Ppds[_i])),zorder=41, weight = 'bold',fontsize='x-small')


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
		pds_grid_410_660.set_xlim(0,60)
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
		ppds_grid_410_660.set_xlim(0,60)
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

fig.savefig(PP_FIGURE+'SELECTED_BINNED_DATA_CROSS_SECTION_Pds_Ppds_bootstrap.'+EXT_FIG,dpi=DPI_FIG)

print('Ending the Cross section CODE')