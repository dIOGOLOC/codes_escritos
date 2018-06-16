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


from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,STA_DIR,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP
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

filename = PP_SELEC_DIR+'SELECTED_BINNED.json'

SELECTED_BINNED_DATA_dic = json.load(open(filename))



lats = SELECTED_BINNED_DATA_dic['lat']
lons = SELECTED_BINNED_DATA_dic['lon']
RF_number = SELECTED_BINNED_DATA_dic['len']
RF_stacking = SELECTED_BINNED_DATA_dic['data']


print('Calculating earth model layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

print('Plotting Total of Receiver Functions per bin')
print('\n')

fig_receiver_function_per_bin,ax = plt.subplots(1,1,figsize=(10,5))

project_Lat = PROJECT_LAT
project_Lon = PROJECT_LON

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)



x, y = m(lons,lats)
sc = m.scatter(x,y,40,RF_number,cmap='viridis',marker='o')
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
	AB_lon = np.linspace(AB_lon_line[0],AB_lon_line[1],abs(abs(AB_lon_line[1])-abs(AB_lon_line[0]))*4)
	AB_lat = np.linspace(AB_lat_line[0],AB_lat_line[1],abs(abs(AB_lon_line[1])-abs(AB_lon_line[0]))*4)
else:
	AB_lon = np.linspace(AB_lon_line[0],AB_lon_line[1],abs(abs(AB_lat_line[1])-abs(AB_lat_line[0]))*4)
	AB_lat = np.linspace(AB_lat_line[0],AB_lat_line[1],abs(abs(AB_lat_line[1])-abs(AB_lat_line[0]))*4)


fig_receiver_function_profile = plt.figure(figsize=(10,5))

project_Lat = PROJECT_LAT
project_Lon = PROJECT_LON

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)

#xA, yA = m(AB_lon_line[0],AB_lat_line[0])
#plt.text(xA-100000, yA, 'A',fontsize=20)

#xB, yB = m(AB_lon_line[1],AB_lat_line[1])
#plt.text(xB+40000, yB, 'B',fontsize=20)

m.drawgreatcircle(AB_lon_line[0],AB_lat_line[0],AB_lon_line[1],AB_lat_line[1],linewidth=4,color='r')


x, y = m(lons,lats)
sc = m.scatter(x,y,40,RF_number,cmap='viridis',marker='o')
plt.colorbar(sc,label='Total of Receiver Functions per bin')



m.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m.drawcoastlines(color='k',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

plt.title('Receiver Functions per bin', y=1.08)
plt.show()


print('Calculating the distance between cross section and selected grid')
print('\n')

RF_data_profile = []

for i,j in enumerate(AB_lon):
	dist = [np.sqrt((j - lons[k])**2 + (AB_lat[i] - l)**2)  for k,l in enumerate(lats)]
	RF_data_profile.append(RF_stacking[dist.index(min(dist))])


RF_data_profile_stacking = []
for i,j in enumerate(RF_data_profile):
    if len(j) > 10:
        RF_data_profile_stacking.append(j)

    #else:
        #RF_data_profile_stacking.append(np.zeros_like(camadas_terra_10_km))



print('Plotting the Final Figure')

#Cross section figure

fig, (ax1,ax) = plt.subplots(1,2, figsize=(20, 10))


project_Lat = PROJECT_LAT
project_Lon = PROJECT_LON

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)


#xA, yA = m(AB_lon_line[0],AB_lat_line[0])
#plt.text(xA-100000, yA, 'A',fontsize=20)

#xB, yB = m(AB_lon_line[1],AB_lat_line[1])
#plt.text(xB+40000, yB, 'B',fontsize=20)

l2, = m.drawgreatcircle(AB_lon_line[0],AB_lat_line[0],AB_lon_line[1],AB_lat_line[1],linewidth=5,color='r')

x, y = m(lons,lats)
sc = m.scatter(x,y,40,RF_number,cmap='viridis',marker='o')
plt.colorbar(sc,label='Total of Receiver Functions per bin')

for lon, lat in zip(sta_long,sta_lat):
    x,y = m(lon, lat)
    msize = 10
    l1, = m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='grey')

m.fillcontinents(color='whitesmoke',lake_color='#46bcec',zorder=0)
m.drawcoastlines(color='dimgray',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

plt.legend([l1,sc,l2],['Stations','Grid Points','Cross section'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w')
plt.title('BINNED DATA and Cross section', y=1.08)


#Migration figure


factor = 100

majorLocatorX = MultipleLocator(10)
majorLocatorY = MultipleLocator(50)
minorLocatorY = MultipleLocator(10)
minorLocatorX = MultipleLocator(5)


for _i, _j in enumerate(RF_data_profile_stacking):
	RF_data_factor = [_i/factor+l for k, l in enumerate(_j)]
	#ax1.text(_i/factor,820,"{0:.1f}".format(AB_lon[_i]),rotation=45)
	ax1.plot(RF_data_factor,camadas_terra_10_km,'k',linewidth=0.5)
	ax1.yaxis.set_major_locator(majorLocatorY)
	#ax1.xaxis.set_minor_locator(minorLocatorX)
	ax1.yaxis.set_minor_locator(minorLocatorY)
	ax1.grid(True,which='minor',linestyle='--')
	ax1.grid(True,which='major',color='k',linewidth=1)

	ax1.yaxis.set_ticks_position('both')

	ax1.fill_betweenx(camadas_terra_10_km,RF_data_factor,_i/factor,where=np.array(RF_data_factor)>=_i/factor, facecolor='black',alpha=0.8, interpolate=True)
	ax1.fill_betweenx(camadas_terra_10_km,RF_data_factor,_i/factor,where=np.array(RF_data_factor)<=_i/factor, facecolor='gray',alpha=0.6, interpolate=True)
	ax1.set_title('Receiver Function - a = 0.5')
	ax1.set_xticks([])
	ax1.set_ylabel('Depth (km)')
	#ax1.set_xlabel('Longitude')
	ax1.tick_params(labelright=True)
	ax1.set_ylim(MAX_DEPTH,MIN_DEPTH,MAX_DEPTH)
plt.show()
fig.savefig(PP_FIGURE+'SELECTED_BINNED_DATA_CROSS_SECTION.'+EXT_FIG,dpi=DPI_FIG)

print('Ending the Cross section CODE')
