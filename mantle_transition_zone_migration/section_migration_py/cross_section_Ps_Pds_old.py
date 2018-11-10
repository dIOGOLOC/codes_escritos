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
sc = m.scatter(x,y,40,RF_number,cmap='gist_heat',marker='s',edgecolors='none')
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

fig = plt.figure(figsize=(20, 15))

fig.suptitle('Cross section for Pds and Ppds')


gs = gridspec.GridSpec(4, 8)
gs.update(wspace=1, hspace=0.5)

ax = fig.add_subplot(gs[0:2, 0:2])
ax7 = fig.add_subplot(gs[2:4, 0:2])
ax1 = fig.add_subplot(gs[0:2, 4:])
ax2 = fig.add_subplot(gs[2:4, 4:],sharex=ax1)


ax3 = fig.add_subplot(gs[0,2:4])
ax4= fig.add_subplot(gs[1,2:4],sharex=ax3)
ax5 = fig.add_subplot(gs[2,2:4],sharex=ax3)
ax6 = fig.add_subplot(gs[3,2:4],sharex=ax3)

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)


l2, = m.drawgreatcircle(AB_lon_line[0],AB_lat_line[0],AB_lon_line[1],AB_lat_line[1],linewidth=2,color='k')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
x, y = m(lons,lats)
sc = m.scatter(x,y,40,RF_delta_1_Vp_mean,cmap='RdGy',marker='s',edgecolors='none',vmin=-1,vmax=1)
fig.colorbar(sc,cax=cax)

for lon, lat in zip(sta_long,sta_lat):
	x,y = m(lon, lat)
	msize = 5
	l1, = m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m.fillcontinents(color='whitesmoke',lake_color='#46bcec',zorder=0)
m.drawcoastlines(color='dimgray',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

#ax.legend([l1,sc,l2],['Stations','Grid Points','Cross section'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w')
ax.set_title('P-velocity anomaly (%)', y=1.08)

# figure 2 

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL,ax=ax7)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)


l2, = m.drawgreatcircle(AB_lon_line[0],AB_lat_line[0],AB_lon_line[1],AB_lat_line[1],linewidth=2,color='k')

divider = make_axes_locatable(ax7)
cax = divider.append_axes("right", size="5%", pad=0.05)
x, y = m(lons,lats)
sc = m.scatter(x,y,40,RF_DEPTH_true_thickness_MTZ_Pds,cmap='RdGy',marker='s',edgecolors='none',vmin=250-(DEPTH_RANGE*2)-INTER_DEPTH,vmax=250+(DEPTH_RANGE*2)+INTER_DEPTH)
fig.colorbar(sc,cax=cax)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m(lon, lat)
    msize = 5
    l1, = m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m.fillcontinents(color='whitesmoke',lake_color='#46bcec',zorder=0)
m.drawcoastlines(color='dimgray',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])
ax7.set_title('MTZ True Thickness', y=1.08)



#Migration figure


#### Figure Pds  ####


factor_Pds = 300

majorLocatorY = MultipleLocator(50)
minorLocatorY = MultipleLocator(10)


for _i, _j in enumerate(RF_data_profile_Pds):
	RF_data_factor_Pds = [_i/factor_Pds+l for k, l in enumerate(_j)]
	ax1.plot(RF_data_factor_Pds,camadas_terra_10_km,'k',linewidth=0.5)
	ax1.yaxis.set_major_locator(majorLocatorY)
	ax1.yaxis.set_minor_locator(minorLocatorY)
	ax1.grid(True,which='minor',linestyle='--')
	ax1.grid(True,which='major',color='k',linewidth=1)

	ax1.plot(_i/factor_Pds,RF_DEPTH_mean_1_profile_Pds[_i],'ok',ms=3,markerfacecolor='none')
	ax1.plot(_i/factor_Pds,RF_DEPTH_mean_2_profile_Pds[_i],'ok',ms=3,markerfacecolor='none')

	ax1.errorbar(_i/factor_Pds,RF_DEPTH_mean_1_profile_Pds[_i], yerr=RF_DEPTH_std_1_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
	ax1.errorbar(_i/factor_Pds,RF_DEPTH_mean_2_profile_Pds[_i], yerr=RF_DEPTH_std_2_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)

	ax1.yaxis.set_ticks_position('both')

	ax1.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)>=_i/factor_Pds, facecolor='black',alpha=0.5, interpolate=True)
	ax1.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,_i/factor_Pds,where=np.array(RF_data_factor_Pds)<=_i/factor_Pds, facecolor='gray',alpha=0.5, interpolate=True)
	ax1.set_title('Depth - Pds')
	ax1.set_xticks([])
	ax1.set_ylabel('Depth (km)')
	ax1.tick_params(labelright=True)
	ax1.set_ylim(MAX_DEPTH,MIN_DEPTH,MAX_DEPTH)


#Migration figure

#### Figure Ppds  ####


factor_Ppds = 300

for _i, _j in enumerate(RF_data_profile_Ppds):
	RF_data_factor_Ppds = [_i/factor_Ppds+l for k, l in enumerate(_j)]
	ax2.plot(RF_data_factor_Ppds,camadas_terra_10_km,'k',linewidth=0.5)
	ax2.yaxis.set_major_locator(majorLocatorY)
	ax2.yaxis.set_minor_locator(minorLocatorY)
	ax2.grid(True,which='minor',linestyle='--')
	ax2.grid(True,which='major',color='k',linewidth=1)

	ax2.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i],'ok',ms=3,markerfacecolor='none')
	ax2.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i],'ok',ms=3,markerfacecolor='none')

	ax2.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i], yerr=RF_DEPTH_std_1_profile_Ppds[_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)
	ax2.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i], yerr=RF_DEPTH_std_2_profile_Ppds[_i], ecolor='k',elinewidth=1,capsize=1,capthick=1)

	ax2.yaxis.set_ticks_position('both')

	ax2.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)>=_i/factor_Ppds, facecolor='black',alpha=0.5, interpolate=True)
	ax2.fill_betweenx(camadas_terra_10_km,RF_data_factor_Ppds,_i/factor_Ppds,where=np.array(RF_data_factor_Ppds)<=_i/factor_Ppds, facecolor='gray',alpha=0.5, interpolate=True)
	ax2.set_xticks([])
	ax2.set_title('Depth - Ppds')
	ax2.set_ylabel('Depth (km)')

	ax2.tick_params(labelright=True)
	ax2.set_ylim(MAX_DEPTH,MIN_DEPTH,MAX_DEPTH)
	if abs(abs(AB_lon_line[1]) - abs(AB_lon_line[0])) > abs(abs(AB_lat_line[1]) - abs(AB_lat_line[0])):
		ax2.text(_i/factor_Ppds*0.95,820,"{0:.1f}".format(AB_lon[_i]),rotation=-45,fontsize=10)
		ax2.set_xlabel('Longitude ($^\circ$)',labelpad=30)
	else:
		ax2.text(_i/factor_Ppds*0.95,820,"{0:.1f}".format(AB_lat[_i]),rotation=-45,fontsize=10)
		ax2.set_xlabel('Latitude ($^\circ$)',labelpad=30)


#### Figure True and Apparent  410 km Pds and Ppds  ####


for _i, _j in enumerate(RF_data_profile_Ppds):
	ax3.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Pds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')
	ax3.plot(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i],marker='o',markerfacecolor='none',markeredgecolor='gray')
	ax3.plot(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile_Pds[_i],marker='s',markerfacecolor='none',markeredgecolor='red')

	ax3.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Pds[_i], yerr=RF_DEPTH_std_1_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	ax3.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_profile_Ppds[_i], yerr=RF_DEPTH_std_1_profile_Ppds[_i], ecolor='gray',elinewidth=1,capsize=2,capthick=1)
	ax3.errorbar(_i/factor_Ppds,RF_DEPTH_mean_1_true_profile_Pds[_i], yerr=RF_DEPTH_std_1_true_profile_Pds[_i], ecolor='red',elinewidth=1,capsize=2,capthick=1)

	ax3.set_title('410 km Apparent and True Depth')
	ax3.set_ylim(410+(DEPTH_RANGE*2)+INTER_DEPTH,410-(DEPTH_RANGE*2)-INTER_DEPTH)
	ax3.set_ylabel('Depth (km)')
	ax3.yaxis.set_ticks_position('both')
	ax3.yaxis.set_major_locator(MultipleLocator(50))
	ax3.yaxis.set_minor_locator(MultipleLocator(10))
	ax3.grid(True,which='minor',linestyle='--')
	ax3.grid(True,which='major',color='k',linewidth=1)
	ax3.tick_params(labelright=True)
	ax3.set_xticks([])


#### Figure True and Apparent  660 km Pds and Ppds  ####


for _i, _j in enumerate(RF_data_profile_Ppds):
	ax4.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Pds[_i],marker='^',markerfacecolor='none',markeredgecolor='k')
	ax4.plot(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i],marker='^',markerfacecolor='none',markeredgecolor='gray')
	ax4.plot(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile_Pds[_i],marker='s',markerfacecolor='none',markeredgecolor='red')

	ax4.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Pds[_i], yerr=RF_DEPTH_std_2_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	ax4.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_profile_Ppds[_i], yerr=RF_DEPTH_std_2_profile_Ppds[_i], ecolor='gray',elinewidth=1,capsize=2,capthick=1)
	ax4.errorbar(_i/factor_Ppds,RF_DEPTH_mean_2_true_profile_Pds[_i], yerr=RF_DEPTH_std_2_true_profile_Pds[_i], ecolor='red',elinewidth=1,capsize=2,capthick=1)

	ax4.set_title('660 km Apparent and True Depth')
	ax4.set_ylim(660+(DEPTH_RANGE*2)+INTER_DEPTH,660-(DEPTH_RANGE*2)-INTER_DEPTH)
	ax4.set_ylabel('Depth (km)')
	ax4.yaxis.set_ticks_position('both')
	ax4.yaxis.set_major_locator(MultipleLocator(50))
	ax4.yaxis.set_minor_locator(MultipleLocator(10))
	ax4.grid(True,which='minor',linestyle='--')
	ax4.grid(True,which='major',color='k',linewidth=1)
	ax4.tick_params(labelright=True)
	ax4.set_xticks([])

#### Figure 410 km and 660 km P-velocity anomaly  ####


for _i, _j in enumerate(RF_data_profile_Ppds):
	ax5.plot(_i/factor_Ppds,RF_delta_1_Vp_mean_profile[_i],marker='^',markerfacecolor='none',markeredgecolor='k')
	ax5.plot(_i/factor_Ppds,RF_delta_2_Vp_mean_profile[_i],marker='o',markerfacecolor='none',markeredgecolor='k')

	ax5.errorbar(_i/factor_Ppds,RF_delta_1_Vp_mean_profile[_i], yerr=RF_delta_1_Vp_std_profile[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	ax5.errorbar(_i/factor_Ppds,RF_delta_2_Vp_mean_profile[_i], yerr=RF_delta_2_Vp_std_profile[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	ax5.set_ylim(-1,1)
	ax5.set_title('P-velocity anomaly (%)')
	ax5.yaxis.set_ticks_position('both')
	ax5.yaxis.set_major_locator(MultipleLocator(0.4))
	ax5.yaxis.set_minor_locator(MultipleLocator(0.2))
	ax5.grid(True,which='major',color='k',linewidth=1)
	ax5.grid(True,which='minor',linestyle='--')
	ax5.tick_params(labelright=True)
	ax5.set_xticks([])

	#### Figure MTZ True and Apparent thickness  ####


for _i, _j in enumerate(RF_data_profile_Ppds):
	ax6.plot(_i/factor_Ppds,RF_DEPTH_mtz_thickness_profile_Pds[_i],marker='o',markerfacecolor='none',markeredgecolor='k')
	ax6.plot(_i/factor_Ppds,RF_DEPTH_true_thickness_MTZ_profile_Pds[_i],marker='s',markerfacecolor='none',markeredgecolor='r')

	ax6.errorbar(_i/factor_Ppds,RF_DEPTH_mtz_thickness_profile_Pds[_i], yerr=RF_DEPTH_std_2_profile_Pds[_i], ecolor='k',elinewidth=1,capsize=2,capthick=1)
	ax6.errorbar(_i/factor_Ppds,RF_DEPTH_true_thickness_MTZ_profile_Pds[_i], yerr=RF_DEPTH_std_2_profile_Ppds[_i], ecolor='r',elinewidth=1,capsize=2,capthick=1)
	ax6.set_ylim(250+(DEPTH_RANGE*2)+INTER_DEPTH,250-(DEPTH_RANGE*2)-INTER_DEPTH)
	ax6.set_ylabel('Depth (km)')
	ax6.set_title('MTZ Apparent and True Thickness')
	ax6.yaxis.set_ticks_position('both')
	ax6.yaxis.set_major_locator(MultipleLocator(50))
	ax6.yaxis.set_minor_locator(MultipleLocator(10))
	ax6.grid(True,which='major',color='k',linewidth=1)
	ax6.grid(True,which='minor',linestyle='--')
	ax6.tick_params(labelright=True)

	if abs(abs(AB_lon_line[1]) - abs(AB_lon_line[0])) > abs(abs(AB_lat_line[1]) - abs(AB_lat_line[0])):
		ax6.text(_i/factor_Ppds*0.95,250+(DEPTH_RANGE*2)+(2*INTER_DEPTH),"{0:.1f}".format(AB_lon[_i]),rotation=-45,fontsize=10)
		ax6.set_xlabel('Longitude ($^\circ$)',labelpad=30)
	else:
		ax6.text(_i/factor_Ppds*0.95,250+(DEPTH_RANGE*2)+(2*INTER_DEPTH),"{0:.1f}".format(AB_lat[_i]),rotation=-45,fontsize=10)
		ax6.set_xlabel('Latitude ($^\circ$)',labelpad=30)




plt.show()

fig.savefig(PP_FIGURE+'SELECTED_BINNED_DATA_CROSS_SECTION_Pds_Ppds.'+EXT_FIG,dpi=DPI_FIG)

print('Ending the Cross section CODE')