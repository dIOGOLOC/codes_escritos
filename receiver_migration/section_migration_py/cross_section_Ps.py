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
					RF_DIR,RF_EXT,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,
					DIST_T_DIR,DEPTH_T_DIR,TIME_T_DIR,DIST_PP_DIR,DEPTH_PP_DIR,TIME_PP_DIR,LAT_PP_DIR,
					LON_PP_DIR,PP_SELEC_DIR,DIST_GRID_PP,NUMBER_PP_PER_BIN,RAY_TRACE_PLOT,RAY_TRACE_410_660_PLOT,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG
				   )
filename = PP_SELEC_DIR+'SELECTED_BINNED.json'

print(filename)
SELECTED_BINNED_DATA_dic = json.load(open(filename))


#definition of the coor catcher


fig_receiver_function_per_bin=plt.figure(figsize=(10,5))

project_Lat = PROJECT_LAT
project_Lon = PROJECT_LON

m = Basemap(resolution='l',projection='merc',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,llcrnrlon=LLCRNRLON_SMALL,
            llcrnrlat=LLCRNRLAT_SMALL,urcrnrlon=URCRNRLON_SMALL,urcrnrlat=URCRNRLAT_SMALL)

m.readshapefile(BOUNDARY_1_SHP,name=BOUNDARY_1_SHP_NAME,linewidth=3)
m.readshapefile(BOUNDARY_2_SHP,name=BOUNDARY_2_SHP_NAME,linewidth=0.7)


lats = SELECTED_BINNED_DATA_dic['lat']
lons = SELECTED_BINNED_DATA_dic['lon']
RF_number = SELECTED_BINNED_DATA_dic['len']
x, y = m(lons,lats)
sc = m.scatter(x,y,40,RF_number,cmap='viridis',marker='o')
plt.colorbar(sc,label='Total of Receiver Functions per bin')



m.fillcontinents(color='whitesmoke',lake_color=None,zorder=2,alpha=0.1)
m.drawcoastlines(color='k',zorder=1)
m.drawmeridians(np.arange(0, 360, 5),color='lightgrey',labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 90, 5),color='lightgrey',labels=[True,True,True,True])

plt.title('Receiver Functions per bin', y=1.08)

ax1 = fig_receiver_function_per_bin.add_subplot(111)
ax1.set_title('custom picker for line data')
line, = ax1.plot(rand(100), rand(100), 'o', picker=line_picker)
fig.canvas.mpl_connect('pick_event', onpick2)

plt.show()


'''

# Creating cross-section

# A-B cross-section

AB_lon_line = [-50,-41]
AB_lat_line = [-6,-4]

#AB_lon_line = [-50,-41]
#AB_lat_line = [-7.5,-5.5]

AB_lon = np.linspace(AB_lon_line[0],AB_lon_line[1],abs(abs(AB_lon_line[1])-abs(AB_lon_line[0]))*4)
AB_lat = np.linspace(AB_lat_line[0],AB_lat_line[1],abs(abs(AB_lon_line[1])-abs(AB_lon_line[0]))*4)

fig=plt.figure(figsize=(20,10))

project_Lat = -5
project_Lon = -45

m = Basemap(resolution='l',projection='merc',lat_0=project_Lat, lon_0=project_Lon,llcrnrlon=-51.,
            llcrnrlat=-12.,urcrnrlon=-40.,urcrnrlat=-2)

m.readshapefile('/home/diogo/dados_doutorado/parnaiba_basin/SIG_dados/shapes/bacia_parnaiba/bacia_parnaiba',name='bacia',linewidth=3)
m.readshapefile('/home/diogo/dados_doutorado/parnaiba_basin/SIG_dados/shapes/Estados/Brasil',name='estados',linewidth=0.7)


xA, yA = m(AB_lon_line[0],AB_lat_line[0])
plt.text(xA-100000, yA, 'A',fontsize=20)

xB, yB = m(AB_lon_line[1],AB_lat_line[1])
plt.text(xB+40000, yB, 'B',fontsize=20)

m.drawgreatcircle(AB_lon_line[0],AB_lat_line[0],AB_lon_line[1],AB_lat_line[1],linewidth=5,color='r')

for lon, lat in zip(AB_lon,AB_lat):
    x,y = m(lon, lat)
    msize = 10
    #l1, = m.plot(x, y, 'o',markersize=msize,markeredgecolor='k',markerfacecolor='y')

lats = grid_sel_y
lons = grid_sel_x
RF_number = len_RF_stacking
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

plt.legend([l1,sc],['Stations','Grid Points'],scatterpoints=1, frameon=True,labelspacing=1, loc='upper right',facecolor='w')
plt.title('BIN DATA and Profile', y=1.08)


RF_data_profile = []
RF_RAY_profile = []

for i,j in enumerate(AB_lon):
    dist = [np.sqrt((j - grid_sel_x[k])**2 + (AB_lat[i] - l)**2)  for k,l in enumerate(grid_sel_y)]
    RF_data_profile.append(RF_stacking[dist.index(min(dist))])
    RF_RAY_profile.append(RF_ray_stacking[dist.index(min(dist))])


# In[46]:


RF_data_profile_stacking = []
RF_RAY_profile_stacking = []
for i,j in enumerate(RF_data_profile):
    if len(j) > 10:
        RF_data_profile_stacking.append(j)
        RF_RAY_profile_stacking.append(RF_RAY_profile[i])

    #else:
        #RF_data_profile_stacking.append(np.zeros_like(camadas_terra_10_km))
        #RF_RAY_profile_stacking.append(0)


# In[47]:


factor = 80

fig, ax = plt.subplots(1, 1, figsize=(15, 5))
#fig, ax = plt.subplots(1, 1)

majorLocatorX = MultipleLocator(10)
majorLocatorY = MultipleLocator(50)
minorLocatorY = MultipleLocator(10)
minorLocatorX = MultipleLocator(5)


for i, j in enumerate(RF_data_profile_stacking): 
    #plt.plot(time_real-10,i/factor+j,'k',linewidth=0.5)
    #ax.text(i/factor,-25,"{0:.2f}".format(RF_RAY_profile_stacking[i]),rotation=45)
    ax.plot(i/factor+j,camadas_terra_10_km,'k',linewidth=0.5)
    ax.yaxis.set_major_locator(majorLocatorY)
    ax.xaxis.set_minor_locator(minorLocatorX)
    ax.yaxis.set_minor_locator(minorLocatorY)
    ax.grid(True,which='minor',linestyle='--')
    ax.grid(True,which='major',color='k',linewidth=1)

    ax.yaxis.set_ticks_position('both')

    ax.fill_betweenx(camadas_terra_10_km,i/factor+j,i/factor,where=(i/factor+j)>=i/factor, facecolor='black',alpha=0.6, interpolate=True)
    ax.fill_betweenx(camadas_terra_10_km,i/factor+j,i/factor,where=(i/factor+j)<=i/factor, facecolor='gray',alpha=0.6, interpolate=True)
    ax.set_title('Funcoes do Receptor - a = 0.5')
    ax.set_xticks([])
    ax.set_ylabel('Depth (km)')
    
    ax.tick_params(labelright=True)
    ax.set_ylim(800,300)


# # Figura completa

# In[48]:


fig, (ax1,ax) = plt.subplots(1,2, figsize=(20, 10))

project_Lat = -5
project_Lon = -45

m = Basemap(resolution='l',projection='merc',lat_0=project_Lat, lon_0=project_Lon,llcrnrlon=-51.,
            llcrnrlat=-12.,urcrnrlon=-40.,urcrnrlat=-2,ax=ax)

m.readshapefile('/home/diogo/dados_doutorado/parnaiba_basin/SIG_dados/shapes/bacia_parnaiba/bacia_parnaiba',name='bacia',linewidth=3)
m.readshapefile('/home/diogo/dados_doutorado/parnaiba_basin/SIG_dados/shapes/Estados/Brasil',name='estados',linewidth=0.7)


xA, yA = m(AB_lon_line[0],AB_lat_line[0])
plt.text(xA-100000, yA, 'A',fontsize=20)

xB, yB = m(AB_lon_line[1],AB_lat_line[1])
plt.text(xB+40000, yB, 'B',fontsize=20)

m.drawgreatcircle(AB_lon_line[0],AB_lat_line[0],AB_lon_line[1],AB_lat_line[1],linewidth=5,color='r')

for lon, lat in zip(AB_lon,AB_lat):
    x,y = m(lon, lat)
    msize = 10
    #l1, = m.plot(x, y, 'o',markersize=msize,markeredgecolor='k',markerfacecolor='y')

lats = grid_sel_y
lons = grid_sel_x
RF_number = len_RF_stacking
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

plt.legend([l1,sc],['Stations','Grid Points'],scatterpoints=1, frameon=True,labelspacing=1, loc='lower right',facecolor='w')
plt.title('BIN DATA and Profile', y=1.08)

factor = 100

majorLocatorX = MultipleLocator(10)
majorLocatorY = MultipleLocator(50)
minorLocatorY = MultipleLocator(10)
minorLocatorX = MultipleLocator(5)


for i, j in enumerate(RF_data_profile_stacking): 
    #plt.plot(time_real-10,i/factor+j,'k',linewidth=0.5)
    #ax.text(i/factor,-25,"{0:.2f}".format(RF_RAY_profile_stacking[i]),rotation=45)
    ax1.plot(i/factor+j,camadas_terra_10_km,'k',linewidth=0.5)
    ax1.yaxis.set_major_locator(majorLocatorY)
    ax1.xaxis.set_minor_locator(minorLocatorX)
    ax1.yaxis.set_minor_locator(minorLocatorY)
    ax1.grid(True,which='minor',linestyle='--')
    ax1.grid(True,which='major',color='k',linewidth=1)

    ax1.yaxis.set_ticks_position('both')

    ax1.fill_betweenx(camadas_terra_10_km,i/factor+j,i/factor,where=(i/factor+j)>=i/factor, facecolor='black',alpha=0.6, interpolate=True)
    ax1.fill_betweenx(camadas_terra_10_km,i/factor+j,i/factor,where=(i/factor+j)<=i/factor, facecolor='gray',alpha=0.6, interpolate=True)
    ax1.set_title('Receiver Function - a = 0.5')
    ax1.set_xticks([])
    ax1.set_ylabel('Depth (km)')
    
    ax1.tick_params(labelright=True)
    ax1.set_ylim(800,300)
'''
