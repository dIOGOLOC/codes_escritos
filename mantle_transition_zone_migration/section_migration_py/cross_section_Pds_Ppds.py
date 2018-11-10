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
sc = m.scatter(x,y,40,RF_number,cmap='magma',marker='s',edgecolors='none')
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

print('Plotting the Final Figure')

#Cross section figure

fig = plt.figure(figsize=(20, 30))

fig.suptitle('Cross section for Pds and Ppds')

gs = gridspec.GridSpec(5, 3)
gs.update(wspace=0.2, hspace=0.5)


MTZ_thickness = fig.add_subplot(gs[0,1:])

pefil_pds = fig.add_subplot(gs[1:3,1:],sharex=MTZ_thickness)
pefil_ppds = fig.add_subplot(gs[3:5,1:],sharex=MTZ_thickness)


P_anomaly = fig.add_subplot(gs[0,0])

apparent_410 = fig.add_subplot(gs[1,0],sharex=P_anomaly)
apparent_660 = fig.add_subplot(gs[2,0],sharex=P_anomaly)

apparent_410_ppds = fig.add_subplot(gs[3, 0],sharex=P_anomaly)
apparent_660_ppds = fig.add_subplot(gs[4, 0],sharex=P_anomaly)




#Migration figure


#### Figure Pds  ####


factor_Pds = 150

majorLocatorY = MultipleLocator(100)
minorLocatorY = MultipleLocator(50)


for _i, _j in enumerate(RF_data_profile_Pds):
	RF_data_factor_Pds = [_i/factor_Pds+l for k, l in enumerate(_j)]
	pefil_pds.plot(RF_data_factor_Pds,camadas_terra_10_km,'k',linewidth=0.5)
	pefil_pds.yaxis.set_major_locator(majorLocatorY)
	pefil_pds.yaxis.set_minor_locator(minorLocatorY)
	pefil_pds.grid(True,which='minor',linestyle='--')
	pefil_pds.grid(True,which='major',color='k',linewidth=1)

	pefil_pds.plot(_i/factor_Pds,RF_DEPTH_mean_1_profile_Pds[_i],'ok',ms=3,markerfacecolor='none')
	pefil_pds.plot(_i/factor_Pds,RF_DEPTH_mean_2_profile_Pds[_i],'ok',ms=3,markerfacecolor='none')

	pefil_pds.errorbar(_i/factor_Pds,RF_DEPTH_mean_1_profile_Pds[_i], yerr=RF_DEPTH_std_1_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
	pefil_pds.errorbar(_i/factor_Pds,RF_DEPTH_mean_2_profile_Pds[_i], yerr=RF_DEPTH_std_2_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)

	pefil_pds.yaxis.set_ticks_position('both')

	pefil_pds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)>=_i/factor_Pds, facecolor='black',alpha=0.5, interpolate=True)
	pefil_pds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)<=_i/factor_Pds, facecolor='gray',alpha=0.5, interpolate=True)
	pefil_pds.set_title('Cross-section - Pds')
	pefil_pds.set_xticks([])
	pefil_pds.set_ylabel('Depth (km)')
	pefil_pds.yaxis.set_label_position("right")
	pefil_pds.tick_params(labelright=True)
	pefil_pds.set_ylim(MAX_DEPTH,MIN_DEPTH)
'''
if abs(abs(AB_lon_line[1]) - abs(AB_lon_line[0])) > abs(abs(AB_lat_line[1]) - abs(AB_lat_line[0])):
		pefil_pds.text(_i/factor_Pds*0.95,820,"{0:.1f}".format(AB_lon[_i]),rotation=-45,fontsize=10)
		pefil_pds.set_xlabel('Longitude ($^\circ$)',labelpad=30)
else:
		pefil_pds.text(_i/factor_Pds*0.95,820,"{0:.1f}".format(AB_lat[_i]),rotation=-45,fontsize=10)
		pefil_pds.set_xlabel('Latitude ($^\circ$)',labelpad=30)
'''

#### Figure Ppds  ####


factor_Ppds = 150

for _i, _j in enumerate(RF_data_profile_Ppds):
	RF_data_factor_Ppds = [_i/factor_Ppds+l for k, l in enumerate(_j)]
	pefil_ppds.plot(RF_data_factor_Ppds,camadas_terra_10_km,'k',linewidth=0.5)
	pefil_ppds.yaxis.set_major_locator(majorLocatorY)
	pefil_ppds.yaxis.set_minor_locator(minorLocatorY)
	pefil_ppds.grid(True,which='minor',linestyle='--')
	pefil_ppds.grid(True,which='major',color='k',linewidth=1)

	pefil_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i],'ok',ms=3,markerfacecolor='none')
	pefil_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i],'ok',ms=3,markerfacecolor='none')

	pefil_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i], yerr=RF_DEPTH_std_1_profile_Ppds[_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
	pefil_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i], yerr=RF_DEPTH_std_2_profile_Ppds[_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)

	pefil_ppds.yaxis.set_ticks_position('both')
	pefil_ppds.yaxis.set_label_position("right")

	pefil_ppds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)>=_i/factor_Ppds, facecolor='black',alpha=0.5, interpolate=True)
	pefil_ppds.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)<=_i/factor_Ppds, facecolor='gray',alpha=0.5, interpolate=True)
	pefil_ppds.set_xticks([])
	pefil_ppds.set_title('Cross-section - Ppds')
	pefil_ppds.set_ylabel('Depth (km)')
	pefil_ppds.tick_params(labelleft=True,labelright=True)
	pefil_ppds.set_ylim(MAX_DEPTH,MIN_DEPTH)
	if abs(abs(AB_lon_line[1]) - abs(AB_lon_line[0])) > abs(abs(AB_lat_line[1]) - abs(AB_lat_line[0])):
		pefil_ppds.text(_i/factor_Ppds*0.95,820,"{0:.1f}".format(AB_lon[_i]),rotation=-45,fontsize=10)
		pefil_ppds.set_xlabel('Longitude ($^\circ$)',labelpad=30)
	else:
		pefil_ppds.text(_i/factor_Ppds*0.95,820,"{0:.1f}".format(AB_lat[_i]),rotation=-45,fontsize=10)
		pefil_ppds.set_xlabel('Latitude ($^\circ$)',labelpad=30)


#### Figure True and Apparent  410 km Pds  ####


for _i, _j in enumerate(RF_data_profile_Ppds):
	apparent_410.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Pds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')
	apparent_410.plot(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile_Pds[_i],marker='s',markerfacecolor='none',markeredgecolor='dimgray')

	apparent_410.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Pds[_i], yerr=RF_DEPTH_std_1_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	apparent_410.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile_Pds[_i], yerr=RF_DEPTH_std_1_true_profile_Pds[_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)

	apparent_410.set_title('410 km Pds')
	apparent_410.set_ylabel('Depth (km)')
	apparent_410.yaxis.set_ticks_position('both')
	apparent_410.yaxis.set_major_locator(MultipleLocator(100))
	apparent_410.yaxis.set_minor_locator(MultipleLocator(50))
	apparent_410.grid(True,which='minor',linestyle='--')
	apparent_410.grid(True,which='major',color='k',linewidth=1)
	apparent_410.tick_params(labelleft=True,labelright=False)
	apparent_410.set_xticks([])
	apparent_410.set_ylim(500,300)


#### Figure True and Apparent  660 km Pds  ####


for _i, _j in enumerate(RF_data_profile_Ppds):
	apparent_660.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Pds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')
	apparent_660.plot(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile_Pds[_i],marker='s',markerfacecolor='none',markeredgecolor='dimgray')

	apparent_660.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Pds[_i], yerr=RF_DEPTH_std_2_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	apparent_660.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile_Pds[_i], yerr=RF_DEPTH_std_2_true_profile_Pds[_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)

	apparent_660.set_title('660 km Pds')
	apparent_660.set_ylim(800,600)
	apparent_660.set_ylabel('Depth (km)')
	apparent_660.yaxis.set_ticks_position('both')
	apparent_660.yaxis.set_major_locator(MultipleLocator(100))
	apparent_660.yaxis.set_minor_locator(MultipleLocator(50))
	apparent_660.grid(True,which='minor',linestyle='--')
	apparent_660.grid(True,which='major',color='k',linewidth=1)
	apparent_660.tick_params(labelleft=True,labelright=False)
	apparent_660.set_xticks([])

#### Figure True and Apparent  410 km Ppds  ####


for _i, _j in enumerate(RF_data_profile_Ppds):
	apparent_410_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')
	apparent_410_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile_Pds[_i],marker='s',markerfacecolor='none',markeredgecolor='dimgray')

	apparent_410_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i], yerr=RF_DEPTH_std_1_profile_Ppds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	apparent_410_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile_Pds[_i], yerr=RF_DEPTH_std_1_true_profile_Pds[_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)

	apparent_410_ppds.set_title('410 km Ppds')
	apparent_410_ppds.set_ylim(500,300)
	apparent_410_ppds.set_ylabel('Depth (km)')
	apparent_410_ppds.yaxis.set_ticks_position('both')
	apparent_410_ppds.yaxis.set_major_locator(MultipleLocator(100))
	apparent_410_ppds.yaxis.set_minor_locator(MultipleLocator(50))
	apparent_410_ppds.grid(True,which='minor',linestyle='--')
	apparent_410_ppds.grid(True,which='major',color='k',linewidth=1)
	apparent_410_ppds.tick_params(labelleft=True,labelright=False)
	apparent_410_ppds.set_xticks([])


#### Figure True and Apparent  660 km Ppds  ####


for _i, _j in enumerate(RF_data_profile_Ppds):
	apparent_660_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')
	apparent_660_ppds.plot(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile_Pds[_i],marker='s',markerfacecolor='none',markeredgecolor='dimgray')

	apparent_660_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i], yerr=RF_DEPTH_std_2_profile_Ppds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	apparent_660_ppds.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile_Pds[_i], yerr=RF_DEPTH_std_2_true_profile_Pds[_i], ecolor='dimgray',elinewidth=1,capsize=2,capthick=1)

	apparent_660_ppds.set_title('660 km Ppds')
	apparent_660_ppds.set_ylim(800,600)
	apparent_660_ppds.set_ylabel('Depth (km)')
	apparent_660_ppds.yaxis.set_ticks_position('both')
	apparent_660_ppds.yaxis.set_major_locator(MultipleLocator(100))
	apparent_660_ppds.yaxis.set_minor_locator(MultipleLocator(50))
	apparent_660_ppds.grid(True,which='minor',linestyle='--')
	apparent_660_ppds.grid(True,which='major',color='k',linewidth=1)
	apparent_660_ppds.tick_params(labelleft=True,labelright=False)
	apparent_660_ppds.set_xticks([])
	if abs(abs(AB_lon_line[1]) - abs(AB_lon_line[0])) > abs(abs(AB_lat_line[1]) - abs(AB_lat_line[0])):
		apparent_660_ppds.text(_i/factor_Pds*0.95,820,"{0:.1f}".format(AB_lon[_i]),rotation=-45,fontsize=10)
		apparent_660_ppds.set_xlabel('Longitude ($^\circ$)',labelpad=30)
	else:
		apparent_660_ppds.text(_i/factor_Pds*0.95,820,"{0:.1f}".format(AB_lat[_i]),rotation=-45,fontsize=10)
		apparent_660_ppds.set_xlabel('Latitude ($^\circ$)',labelpad=30)


	#### Figure MTZ True and Apparent thickness  ####

for _i, _j in enumerate(RF_data_profile_Ppds):
	MTZ_thickness.axhline(y=250,linewidth=0.5, color='k')

	MTZ_thickness.plot(_i/factor_Ppds,RF_DEPTH_mtz_thickness_profile_Pds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')
	MTZ_thickness.plot(_i/factor_Ppds,RF_DEPTH_true_thickness_MTZ_profile_Pds[_i],marker='s',markerfacecolor='none',markeredgecolor='gray')

	MTZ_thickness.errorbar(_i/factor_Ppds,RF_DEPTH_mtz_thickness_profile_Pds[_i], yerr=RF_DEPTH_std_2_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	MTZ_thickness.errorbar(_i/factor_Ppds,RF_DEPTH_true_thickness_MTZ_profile_Pds[_i], yerr=RF_DEPTH_std_2_profile_Ppds[_i], ecolor='gray',elinewidth=1,capsize=2,capthick=1)
	MTZ_thickness.set_ylim(350,150)
	MTZ_thickness.set_ylabel('Depth (km)')
	MTZ_thickness.yaxis.set_label_position("right")
	MTZ_thickness.set_title('MTZ Apparent and True Thickness')
	MTZ_thickness.yaxis.set_ticks_position('both')
	MTZ_thickness.yaxis.set_major_locator(MultipleLocator(100))
	MTZ_thickness.yaxis.set_minor_locator(MultipleLocator(50))
	MTZ_thickness.grid(True,which='major',color='k',linewidth=1)
	MTZ_thickness.grid(True,which='minor',linestyle='--')
	MTZ_thickness.tick_params(labelleft=True,labelright=True)


#### Figure 410 km and 660 km P-velocity anomaly  ####


for _i, _j in enumerate(RF_data_profile_Ppds):
	P_anomaly.axhline(y=0,linewidth=0.5, color='k')

	P_anomaly.plot(_i/factor_Ppds,RF_delta_1_Vp_mean_profile[_i],marker='^',markerfacecolor='none',markeredgecolor='k')
	P_anomaly.plot(_i/factor_Ppds,RF_delta_2_Vp_mean_profile[_i],marker='v',markerfacecolor='none',markeredgecolor='k')

	P_anomaly.errorbar(_i/factor_Ppds,RF_delta_1_Vp_mean_profile[_i], yerr=RF_delta_1_Vp_std_profile[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	P_anomaly.errorbar(_i/factor_Ppds,RF_delta_2_Vp_mean_profile[_i], yerr=RF_delta_2_Vp_std_profile[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	P_anomaly.set_ylim(-0.6,0.6)
	P_anomaly.set_title('P-velocity anomaly (%)')
	P_anomaly.yaxis.set_ticks_position('both')
	P_anomaly.yaxis.set_major_locator(MultipleLocator(0.4))
	P_anomaly.yaxis.set_minor_locator(MultipleLocator(0.2))
	P_anomaly.grid(True,which='major',color='k',linewidth=1)
	P_anomaly.grid(True,which='minor',linestyle='--')
	P_anomaly.tick_params(labelleft=True,labelright=False)
	P_anomaly.set_xticks([])

plt.show()

fig.savefig(PP_FIGURE+'SELECTED_BINNED_DATA_CROSS_SECTION_Pds_Ppds.'+EXT_FIG,dpi=DPI_FIG)

print('Ending the Cross section CODE')