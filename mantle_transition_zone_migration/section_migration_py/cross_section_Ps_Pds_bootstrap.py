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

RF_BOOTSTRAP_410_DEPTH_Pds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_410_DEPTH_Pds']
RF_BOOTSTRAP_410_DEPTH_Ppds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_410_DEPTH_Ppds']
RF_BOOTSTRAP_660_DEPTH_Pds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_660_DEPTH_Pds']
RF_BOOTSTRAP_660_DEPTH_Ppds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_660_DEPTH_Ppds']

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
	if len(lon_click) == 2:
		fig_receiver_function_per_bin.canvas.mpl_disconnect(cid)

	return lon_click,lat_click

cid = fig_receiver_function_per_bin.canvas.mpl_connect('button_press_event', onclick)

plt.show()

print('Calculating cross section coordinates')
print('\n')

AB_lon_line = lon_click
AB_lat_line = lat_click



if abs(abs(AB_lon_line[1]) - abs(AB_lon_line[0])) > abs(abs(AB_lat_line[1]) - abs(AB_lat_line[0])):
	AB_lon = np.linspace(AB_lon_line[0],AB_lon_line[1],abs(abs(AB_lon_line[1])-abs(AB_lon_line[0]))*SECTION_NUM)
	AB_lat = np.linspace(AB_lat_line[0],AB_lat_line[1],abs(abs(AB_lon_line[1])-abs(AB_lon_line[0]))*SECTION_NUM)
else:
	AB_lon = np.linspace(AB_lon_line[0],AB_lon_line[1],abs(abs(AB_lat_line[1])-abs(AB_lat_line[0]))*SECTION_NUM)
	AB_lat = np.linspace(AB_lat_line[0],AB_lat_line[1],abs(abs(AB_lat_line[1])-abs(AB_lat_line[0]))*SECTION_NUM)


print('Calculating the distance between cross section and selected grid')
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

RF_BOOTSTRAP_410_DEPTH_Pds_profile = []
RF_BOOTSTRAP_410_DEPTH_Ppds_profile = []
RF_BOOTSTRAP_660_DEPTH_Pds_profile = []
RF_BOOTSTRAP_660_DEPTH_Ppds_profile = []

profile_lon = []

for i,j in enumerate(AB_lon):
	dist = [np.sqrt((j - lons[k])**2 + (AB_lat[i] - l)**2)  for k,l in enumerate(lats)]
	if min(dist) < DIST_GRID_PP/SECTION_NUM:
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

		RF_BOOTSTRAP_410_DEPTH_Pds_profile.append(RF_BOOTSTRAP_410_DEPTH_Pds[dist.index(min(dist))]) 
		RF_BOOTSTRAP_410_DEPTH_Ppds_profile.append(RF_BOOTSTRAP_410_DEPTH_Ppds[dist.index(min(dist))]) 
		RF_BOOTSTRAP_660_DEPTH_Pds_profile.append(RF_BOOTSTRAP_660_DEPTH_Pds[dist.index(min(dist))]) 
		RF_BOOTSTRAP_660_DEPTH_Ppds_profile.append(RF_BOOTSTRAP_660_DEPTH_Ppds[dist.index(min(dist))]) 


	else:
		RF_data_profile_Pds.append(np.zeros_like(RF_stacking_Pds[dist.index(min(dist))]))
		RF_DEPTH_mean_1_profile_Pds.append(np.zeros_like(RF_DEPTH_mean_1_Pds[dist.index(min(dist))]))
		RF_DEPTH_std_1_profile_Pds.append(np.zeros_like(RF_DEPTH_std_1_Pds[dist.index(min(dist))]))
		RF_DEPTH_mean_2_profile_Pds.append(np.zeros_like(RF_DEPTH_mean_2_Pds[dist.index(min(dist))]))
		RF_DEPTH_std_2_profile_Pds.append(np.zeros_like(RF_DEPTH_std_2_Pds[dist.index(min(dist))]))
		RF_DEPTH_mean_1_true_profile_Pds.append(np.zeros_like(RF_DEPTH_mean_1_true_Pds[dist.index(min(dist))]))
		RF_DEPTH_mean_2_true_profile_Pds.append(np.zeros_like(RF_DEPTH_mean_2_true_Pds[dist.index(min(dist))]))
		RF_DEPTH_mtz_thickness_profile_Pds.append(np.zeros_like(RF_DEPTH_mtz_thickness_Pds[dist.index(min(dist))]))
		RF_DEPTH_true_thickness_MTZ_profile_Pds.append(np.zeros_like(RF_DEPTH_true_thickness_MTZ_Pds[dist.index(min(dist))]))


		RF_delta_1_Vp_mean_profile.append(np.zeros_like(RF_delta_1_Vp_mean[dist.index(min(dist))]))
		RF_delta_1_Vp_std_profile.append(np.zeros_like(RF_delta_1_Vp_std[dist.index(min(dist))]))
		RF_delta_1_Vs_mean_profile.append(np.zeros_like(RF_delta_1_Vs_mean[dist.index(min(dist))]))
		RF_delta_1_Vs_std_profile.append(np.zeros_like(RF_delta_1_Vs_std[dist.index(min(dist))]))

		RF_delta_2_Vp_mean_profile.append(np.zeros_like(RF_delta_2_Vp_mean[dist.index(min(dist))]))
		RF_delta_2_Vp_std_profile.append(np.zeros_like(RF_delta_2_Vp_std[dist.index(min(dist))]))
		RF_delta_2_Vs_mean_profile.append(np.zeros_like(RF_delta_2_Vs_mean[dist.index(min(dist))]))
		RF_delta_2_Vs_std_profile.append(np.zeros_like(RF_delta_2_Vs_std[dist.index(min(dist))]))

		RF_data_profile_Ppds.append(np.zeros_like(RF_stacking_Ppds[dist.index(min(dist))]))
		RF_DEPTH_mean_1_profile_Ppds.append(np.zeros_like(RF_DEPTH_mean_1_Ppds[dist.index(min(dist))]))
		RF_DEPTH_std_1_profile_Ppds.append(np.zeros_like(RF_DEPTH_std_1_Ppds[dist.index(min(dist))]))
		RF_DEPTH_mean_2_profile_Ppds.append(np.zeros_like(RF_DEPTH_mean_2_Ppds[dist.index(min(dist))]))
		RF_DEPTH_std_2_profile_Ppds.append(np.zeros_like(RF_DEPTH_std_2_Ppds[dist.index(min(dist))]))
		RF_DEPTH_mean_1_true_profile_Ppds.append(np.zeros_like(RF_DEPTH_mean_1_true_Ppds[dist.index(min(dist))]))
		RF_DEPTH_mean_2_true_profile_Ppds.append(np.zeros_like(RF_DEPTH_mean_2_true_Ppds[dist.index(min(dist))]))
		RF_DEPTH_mtz_thickness_profile_Ppds.append(np.zeros_like(RF_DEPTH_mtz_thickness_Ppds[dist.index(min(dist))]))
		RF_DEPTH_true_thickness_MTZ_profile_Ppds.append(np.zeros_like(RF_DEPTH_true_thickness_MTZ_Ppds[dist.index(min(dist))]))

		RF_stacking_Pds_BOOTSTRAP_profile.append(np.zeros_like(RF_stacking_Pds_BOOTSTRAP[dist.index(min(dist))]))
		RF_stacking_Ppds_BOOTSTRAP_profile.append(np.zeros_like(RF_stacking_Ppds_BOOTSTRAP[dist.index(min(dist))]))

		RF_BOOTSTRAP_410_DEPTH_Pds_profile.append(np.zeros_like(RF_BOOTSTRAP_410_DEPTH_Pds[dist.index(min(dist))])) 
		RF_BOOTSTRAP_410_DEPTH_Ppds_profile.append(np.zeros_like(RF_BOOTSTRAP_410_DEPTH_Ppds[dist.index(min(dist))])) 
		RF_BOOTSTRAP_660_DEPTH_Pds_profile.append(np.zeros_like(RF_BOOTSTRAP_660_DEPTH_Pds[dist.index(min(dist))]))  
		RF_BOOTSTRAP_660_DEPTH_Ppds_profile.append(np.zeros_like(RF_BOOTSTRAP_660_DEPTH_Ppds[dist.index(min(dist))]))  

print('Plotting the Final Figure')

#Cross section figure

fig = plt.figure(figsize=(30, 15))

fig.suptitle('Cross section for Pds and Ppds')


gs = gridspec.GridSpec(4, 6)
gs.update(wspace=0.5, hspace=0.25)

ax3 = fig.add_subplot(gs[0,0:2])
ax4= fig.add_subplot(gs[1,0:2],sharex=ax3)
ax5 = fig.add_subplot(gs[2,0:2],sharex=ax3)
ax6 = fig.add_subplot(gs[3,0:2],sharex=ax3)

ax01 = fig.add_subplot(gs[0, 2:],sharey=ax3)
ax02 = fig.add_subplot(gs[1, 2:],sharex=ax01,sharey=ax4)

ax11 = fig.add_subplot(gs[2, 2:],sharex=ax01,sharey=ax5)
ax12 = fig.add_subplot(gs[3, 2:],sharex=ax01,sharey=ax6)




#Migration figure


factor_Pds = 25

majorLocatorY = MultipleLocator(20)
minorLocatorY = MultipleLocator(10)


for _i, _j in enumerate(RF_data_profile_Pds):
	x_data_410_Pds= []
	for x,c in enumerate(RF_stacking_Pds_BOOTSTRAP_profile[_i]):
		RF_data_factor_Pds_bootstrap = [_i/factor_Pds+l for k, l in enumerate(c)]
		x_data_410_Pds.append(RF_data_factor_Pds_bootstrap)
		ax01.plot(RF_data_factor_Pds_bootstrap,camadas_terra_10_km,'gainsboro',linewidth=0.1,alpha=0.2, zorder=10)

	min_x = [min(a) for a in zip(*x_data_410_Pds)]
	max_x = [max(a) for a in zip(*x_data_410_Pds)]
	ax01.fill_betweenx(y=camadas_terra_10_km,x1=min_x, x2=max_x, facecolor='whitesmoke',alpha=1, interpolate=True, zorder=5)

	RF_data_factor_Pds = [_i/factor_Pds+l for k, l in enumerate(_j)]
	ax01.plot(RF_data_factor_Pds,camadas_terra_10_km,'k',linewidth=2, zorder=30)

	ax01.yaxis.set_ticks_position('both')
	ax01.yaxis.set_major_locator(majorLocatorY)
	ax01.grid(True,which='major',linestyle='--')

	ax01.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)>=_i/factor_Pds, facecolor='dimgrey',interpolate=True, zorder=19)
	ax01.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)<=_i/factor_Pds, facecolor='lightgrey', interpolate=True, zorder=20)
	ax01.set_yticks([])
	ax01.tick_params(labelright=False,labelleft=True)


	
	ax01.set_title('410 km - Pds')
	ax01.set_ylim(410+(DEPTH_RANGE*2)+INTER_DEPTH,410-(DEPTH_RANGE*2)-INTER_DEPTH)

for _i, _j in enumerate(RF_data_profile_Pds):
	x_data_660_Pds = []
	for x,c in enumerate(RF_stacking_Pds_BOOTSTRAP_profile[_i]):
		RF_data_factor_Pds_bootstrap = [_i/factor_Pds+l for k, l in enumerate(c)]
		x_data_660_Pds.append(RF_data_factor_Pds_bootstrap)
		ax02.plot(RF_data_factor_Pds_bootstrap,camadas_terra_10_km,'gainsboro',linewidth=0.1,alpha=0.2, zorder=10)

	min_x = [min(a) for a in zip(*x_data_660_Pds)]
	max_x = [max(a) for a in zip(*x_data_660_Pds)]
	ax02.fill_betweenx(y=camadas_terra_10_km,x1=min_x, x2=max_x, facecolor='whitesmoke',alpha=1, interpolate=True, zorder=5)

	RF_data_factor_Pds = [_i/factor_Pds+l for k, l in enumerate(_j)]
	ax02.plot(RF_data_factor_Pds,camadas_terra_10_km,'k',linewidth=2, zorder=30)

	ax02.yaxis.set_major_locator(majorLocatorY)
	ax02.grid(True,which='major',linestyle='--')

	ax02.yaxis.set_ticks_position('both')

	ax02.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)>=_i/factor_Pds, facecolor='dimgrey',interpolate=True, zorder=19)
	ax02.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)<=_i/factor_Pds, facecolor='lightgrey',interpolate=True, zorder=20)
	ax02.set_title('660 km - Pds')
	ax02.set_ylim(660+(DEPTH_RANGE*2)+INTER_DEPTH,660-(DEPTH_RANGE*2)-INTER_DEPTH)
	ax02.tick_params(labelright=False,labelleft=True)


#Migration figure


factor_Ppds = 25

for _i, _j in enumerate(RF_data_profile_Ppds):
	x_data_410_Ppds = []
	for x,c in enumerate(RF_stacking_Ppds_BOOTSTRAP_profile[_i]):
		RF_data_factor_Ppds_bootstrap = [_i/factor_Ppds+l for k, l in enumerate(c)]
		x_data_410_Ppds.append(RF_data_factor_Ppds_bootstrap)
		ax11.plot(RF_data_factor_Ppds_bootstrap,camadas_terra_10_km,'gainsboro',linewidth=0.1,alpha=0.2, zorder=10)		

	min_x = [min(a) for a in zip(*x_data_410_Ppds)]
	max_x = [max(a) for a in zip(*x_data_410_Ppds)]
	ax11.fill_betweenx(y=camadas_terra_10_km,x1=min_x, x2=max_x,  facecolor='whitesmoke',alpha=1, interpolate=True, zorder=5)

	RF_data_factor_Ppds = [_i/factor_Ppds+l for k, l in enumerate(_j)]
	ax11.plot(RF_data_factor_Ppds,camadas_terra_10_km,'k',linewidth=2, zorder=30)

	ax11.yaxis.set_major_locator(majorLocatorY)
	ax11.yaxis.set_minor_locator(minorLocatorY)
	ax11.grid(True,which='major',linestyle='--')

	ax11.yaxis.set_ticks_position('both')

	ax11.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)>=_i/factor_Ppds, facecolor='dimgrey',interpolate=True, zorder=19)
	ax11.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)<=_i/factor_Ppds, facecolor='lightgrey',interpolate=True, zorder=20)
	
	ax11.set_title('410 km - Ppds')
	ax11.tick_params(labelright=False,labelleft=True)

	ax11.set_ylim(410+(DEPTH_RANGE*2)+INTER_DEPTH,410-(DEPTH_RANGE*2)-INTER_DEPTH)

for _i, _j in enumerate(RF_data_profile_Ppds):
	x_data_660_Ppds = []
	for x,c in enumerate(RF_stacking_Ppds_BOOTSTRAP_profile[_i]):
		RF_data_factor_Ppds_bootstrap = [_i/factor_Ppds+l for k, l in enumerate(c)]
		x_data_660_Ppds.append(RF_data_factor_Ppds_bootstrap)
		ax12.plot(RF_data_factor_Ppds_bootstrap,camadas_terra_10_km,'gainsboro',linewidth=0.1,alpha=0.2, zorder=10)


	min_x = [min(a) for a in zip(*x_data_660_Ppds)]
	max_x = [max(a) for a in zip(*x_data_660_Ppds)]
	ax12.fill_betweenx(y=camadas_terra_10_km,x1=min_x, x2=max_x, facecolor='whitesmoke',alpha=1, interpolate=True, zorder=5)

	RF_data_factor_Ppds = [_i/factor_Ppds+l for k, l in enumerate(_j)]
	ax12.plot(RF_data_factor_Ppds,camadas_terra_10_km,'k',linewidth=2, zorder=30)

	ax12.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)>=_i/factor_Ppds, facecolor='dimgrey',interpolate=True, zorder=19)
	ax12.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)<=_i/factor_Ppds, facecolor='lightgrey',interpolate=True, zorder=20)
	
	ax12.yaxis.set_major_locator(majorLocatorY)
	ax12.grid(True,which='major',linestyle='--')

	ax12.yaxis.set_ticks_position('both')

	ax12.set_title('660 km - Ppds')
	ax12.set_xticks([])
	ax12.tick_params(labelright=False,labelleft=True)


	ax12.set_ylim(660+(DEPTH_RANGE*2)+INTER_DEPTH,660-(DEPTH_RANGE*2)-INTER_DEPTH)
	if abs(abs(AB_lon_line[1]) - abs(AB_lon_line[0])) > abs(abs(AB_lat_line[1]) - abs(AB_lat_line[0])):
		ax12.text(_i/factor_Ppds*0.95,660+(DEPTH_RANGE*2)+2*INTER_DEPTH,"{0:.1f}".format(AB_lon[_i]),rotation=-45,fontsize=10)
		ax12.set_xlabel('Longitude ($^\circ$)',labelpad=30)
	else:
		ax12.text(_i/factor_Ppds*0.95,660+(DEPTH_RANGE*2)+2*INTER_DEPTH,"{0:.1f}".format(AB_lat[_i]),rotation=-45,fontsize=10)
		ax12.set_xlabel('Latitude ($^\circ$)',labelpad=30)

#### Plot Depth ####

for _i, _j in enumerate(RF_data_profile_Ppds):
	ax3.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Pds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')

	for x,c in enumerate(RF_BOOTSTRAP_410_DEPTH_Pds_profile[_i]):
		ax3.plot(_i/factor_Pds,c,'o',markeredgecolor='gray',ms=3,markerfacecolor='none',alpha=0.2)

	ax3.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Pds[_i], yerr=RF_DEPTH_std_1_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)

	ax3.set_title('410 km - Pds')
	ax3.set_ylim(410+(DEPTH_RANGE*2)+INTER_DEPTH,410-(DEPTH_RANGE*2)-INTER_DEPTH)
	ax3.set_ylabel('Depth (km)')
	ax3.yaxis.set_ticks_position('both')
	ax3.yaxis.set_major_locator(MultipleLocator(20))
	ax3.grid(True,which='major',linestyle='--')


for _i, _j in enumerate(RF_data_profile_Ppds):
	ax4.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Pds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')

	for x,c in enumerate(RF_BOOTSTRAP_660_DEPTH_Pds_profile[_i]):
		ax4.plot(_i/factor_Pds,c,'o',markeredgecolor='gray',ms=3,markerfacecolor='none',alpha=0.2)


	ax4.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Pds[_i], yerr=RF_DEPTH_std_2_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)

	ax4.set_title('660 km - Pds')
	ax4.set_ylim(660+(DEPTH_RANGE*2)+INTER_DEPTH,660-(DEPTH_RANGE*2)-INTER_DEPTH)
	ax4.set_ylabel('Depth (km)')
	ax4.yaxis.set_ticks_position('both')
	ax4.yaxis.set_major_locator(MultipleLocator(20))
	ax4.grid(True,which='major',linestyle='--')

for _i, _j in enumerate(RF_data_profile_Ppds):
	ax5.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')

	for x,c in enumerate(RF_BOOTSTRAP_410_DEPTH_Pds_profile[_i]):
		ax5.plot(_i/factor_Pds,c,'o',markeredgecolor='gray',ms=3,markerfacecolor='none',alpha=0.2)

	ax5.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i], yerr=RF_DEPTH_std_1_profile_Ppds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)

	ax5.set_title('410 km - Ppds')
	ax5.set_ylim(410+(DEPTH_RANGE*2)+INTER_DEPTH,410-(DEPTH_RANGE*2)-INTER_DEPTH)
	ax5.set_ylabel('Depth (km)')
	ax5.yaxis.set_ticks_position('both')
	ax5.yaxis.set_major_locator(MultipleLocator(20))
	ax5.yaxis.set_minor_locator(MultipleLocator(10))
	ax5.grid(True,which='major',linestyle='--')

for _i, _j in enumerate(RF_data_profile_Ppds):
	ax6.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')


	for x,c in enumerate(RF_BOOTSTRAP_660_DEPTH_Ppds_profile[_i]):
		ax6.plot(_i/factor_Ppds,c,'o',markeredgecolor='gray',ms=3,markerfacecolor='none',alpha=0.2)


	ax6.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i], yerr=RF_DEPTH_std_2_profile_Ppds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)

	ax6.set_title('660 km - Ppds')
	ax6.set_ylim(660+(DEPTH_RANGE*2)+INTER_DEPTH,660-(DEPTH_RANGE*2)-INTER_DEPTH)
	ax6.set_ylabel('Depth (km)')
	ax6.yaxis.set_ticks_position('both')
	ax6.yaxis.set_major_locator(MultipleLocator(20))
	ax6.yaxis.set_minor_locator(MultipleLocator(10))
	ax6.grid(True,which='major',linestyle='--')
	ax6.set_xticks([])

	if abs(abs(AB_lon_line[1]) - abs(AB_lon_line[0])) > abs(abs(AB_lat_line[1]) - abs(AB_lat_line[0])):
		ax6.text(_i/factor_Ppds*0.95,660+(DEPTH_RANGE*2)+2*INTER_DEPTH,"{0:.1f}".format(AB_lon[_i]),rotation=-45,fontsize=10)
		ax6.set_xlabel('Longitude ($^\circ$)',labelpad=30)
	else:
		ax6.text(_i/factor_Ppds*0.95,660+(DEPTH_RANGE*2)+INTER_DEPTH+10,"{0:.1f}".format(AB_lat[_i]),rotation=-45,fontsize=10)
		ax6.set_xlabel('Latitude ($^\circ$)',labelpad=30)

plt.show()

fig.savefig(PP_FIGURE+'SELECTED_BINNED_DATA_CROSS_SECTION_Pds_Ppds_bootstrap.'+EXT_FIG,dpi=DPI_FIG)

print('Ending the Cross section CODE')