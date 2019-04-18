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
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,
					NUMBER_PP_PER_BIN,GRID_PP_MULT,COLORMAP_STD,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,OUTPUT_DIR,		
					EXT_FIG,DPI_FIG,DIST_GRID_PP,DEPTH_RANGE,COLORMAP_VEL,
					NUMBER_PP_PER_BIN,NUMBER_STA_PER_BIN,DEPTH_TARGET
				   )


print('Starting Cross section CODE')
print('\n')

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Stations'+'/'


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

PP_SELEC_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'SELECTED_BINNED_DATA'+'/'

RESULTS_FOLDER_BINS = PP_SELEC_DIR+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
filename = RESULTS_FOLDER_BINS+'SELECTED_BINNED.json'

SELECTED_BINNED_DATA_dic = json.load(open(filename))

lats = SELECTED_BINNED_DATA_dic['lat']
lons = SELECTED_BINNED_DATA_dic['lon']

RF_number = SELECTED_BINNED_DATA_dic['len_Pds']

RF_stacking_Pds = SELECTED_BINNED_DATA_dic['data_Pds']
RF_stacking_Pds_std = SELECTED_BINNED_DATA_dic['data_Pds_std']

RF_stacking_Pds_BOOTSTRAP = SELECTED_BINNED_DATA_dic['data_BOOTSTRAP_Pds']

RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds']

RF_BOOTSTRAP_DEPTH_mean_1_Pds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_1_Pds']

RF_BOOTSTRAP_DEPTH_mean_520_Pds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_520_Pds']

RF_BOOTSTRAP_DEPTH_mean_2_Pds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_2_Pds']

RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds = SELECTED_BINNED_DATA_dic['RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds']


RF_DEPTH_mean_LVZ_Pds = SELECTED_BINNED_DATA_dic['mean_LVZ_Pds']
RF_DEPTH_std_LVZ_Pds = SELECTED_BINNED_DATA_dic['std_LVZ_Pds']

RF_DEPTH_mean_1_Pds = SELECTED_BINNED_DATA_dic['mean_1_Pds']
RF_DEPTH_std_1_Pds = SELECTED_BINNED_DATA_dic['std_1_Pds']

RF_DEPTH_mean_520_Pds = SELECTED_BINNED_DATA_dic['mean_520_Pds']
RF_DEPTH_std_520_Pds = SELECTED_BINNED_DATA_dic['std_520_Pds']

RF_DEPTH_mean_2_Pds = SELECTED_BINNED_DATA_dic['mean_2_Pds']
RF_DEPTH_std_2_Pds = SELECTED_BINNED_DATA_dic['std_2_Pds']

RF_DEPTH_mean_LVZ_700_Pds = SELECTED_BINNED_DATA_dic['mean_LVZ_700_Pds']
RF_DEPTH_std_LVZ_700_Pds = SELECTED_BINNED_DATA_dic['std_LVZ_700_Pds']

RF_DEPTH_mtz_thickness_Pds = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds']
RF_DEPTH_mtz_thickness_Pds_std = SELECTED_BINNED_DATA_dic['mtz_thickness_Pds_std']

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
ax.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True)

ax.grid(True,which='major',color='gray',linewidth=1,linestyle='None')

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
ax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

norm_410 = mpl.colors.Normalize(vmin=200,vmax=300,clip=True)

for i,j in enumerate(lons):
#	if math.isnan(RF_DEPTH_mtz_thickness_Pds[i]) == False and math.isnan(RF_DEPTH_mean_1_Pds[i]) == False  and math.isnan(RF_DEPTH_mean_1_Pds[i]) == False and  math.isnan(RF_DEPTH_mean_LVZ_Pds[i]) == False and math.isnan(RF_DEPTH_mean_520_Pds[i]) == False:
	if math.isnan(RF_DEPTH_mtz_thickness_Pds[i]) == False and math.isnan(RF_DEPTH_mean_1_Pds[i]) == False  and math.isnan(RF_DEPTH_mean_1_Pds[i]) == False and  math.isnan(RF_DEPTH_mean_LVZ_Pds[i]) == False:
		circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i],lats[i]),color=colormap(norm_410(RF_DEPTH_mtz_thickness_Pds[i])), ec='k',transform=ccrs.Geodetic(),zorder=3)
		ax.add_patch(circulo_410)
		circulo_410.pickable()
		circulo_410.set_picker(True)

ax.plot(sta_long,sta_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())

#______________________________________________________________________

sm_410 = plt.cm.ScalarMappable(cmap=colormap,norm=norm_410)
sm_410._A = []
fig.colorbar(sm_410,ax=ax,orientation='horizontal',shrink=0.8)

plt.title('Pick four points for bootstrapping cross-section and them close the windows', y=1.08)


lon_click = []
lat_click = []

def onpick1(event):
	if isinstance(event.artist,Circle):
		patch = event.artist

		patch_lon = float("%.4f" % round(patch.center[0],4))
		patch_lat = float("%.4f" % round(patch.center[1],4))

		print('Return lon/lat of the rectangle selected (left,bottom)')
		print('lon='+str(patch_lon), 'lat='+str(patch_lat))

		lon_click.append(patch_lon)
		lat_click.append(patch_lat)
		retangulo_PICK = Circle(radius=DIST_GRID_PP,xy=(patch_lon, patch_lat),color='k', ec='None',linewidth=1,transform=ccrs.Geodetic(),zorder=5)
		ax.add_patch(retangulo_PICK)
	plt.draw()
	if len(lon_click) == 4:
		fig.canvas.mpl_disconnect(cid)

	return lon_click,lat_click
    

cid = fig.canvas.mpl_connect('pick_event', onpick1)

plt.show()

PP_FIGURE = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Figures'+'/'

RESULTS_FOLDER = PP_FIGURE+'/'+'RESULTS_NUMBER_PP_PER_BIN_'+str(NUMBER_PP_PER_BIN)+'_NUMBER_STA_PER_BIN_'+str(NUMBER_STA_PER_BIN)+'/'
os.makedirs(RESULTS_FOLDER,exist_ok=True)
fig.savefig(RESULTS_FOLDER+'SELECTED_BINNED_DATA_CROSS_SECTION_Pds_Ppds_bootstrap_points.'+EXT_FIG,dpi=DPI_FIG)

print('Allocating points')
print('\n')

#Profile lat/lon

RF_lat_profile = []
RF_lon_profile = []

#Profile Data

RF_data_profile_Pds = []
RF_data_profile_Pds_std = []

#350 LVZ

RF_DEPTH_mean_LVZ_profile_Pds = []
RF_DEPTH_std_LVZ_profile_Pds = []

#P410s

RF_DEPTH_mean_1_profile_Pds = []
RF_DEPTH_std_1_profile_Pds = []

#P520s

RF_DEPTH_mean_520_profile_Pds = []
RF_DEPTH_std_520_profile_Pds = []

#P660s

RF_DEPTH_mean_2_profile_Pds = []
RF_DEPTH_std_2_profile_Pds = []

#700 LVZ

RF_DEPTH_mean_LVZ_700_profile_Pds = []
RF_DEPTH_std_LVZ_700_profile_Pds = []

#MTZ Pds

RF_DEPTH_mtz_thickness_profile_Pds = [] 
RF_DEPTH_mtz_thickness_profile_Pds_std = []

#Bootstrap Data Receiver Functions Pds

RF_stacking_Pds_BOOTSTRAP_profile = []

#Bootstrap Data Mean 350 LVZ

RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds_profile = []

#Bootstrap Data Mean P410s

RF_BOOTSTRAP_DEPTH_mean_1_Pds_profile = []

#Bootstrap Data Mean P520s

RF_BOOTSTRAP_DEPTH_mean_520_Pds_profile = []

#Bootstrap Data Mean P660s

RF_BOOTSTRAP_DEPTH_mean_2_Pds_profile = []

#Bootstrap Data Mean 700 LVZ

RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds_profile = []

lat_lon = [(lons[k],lats[k]) for k,l in enumerate(lats)]

for i,j in enumerate(lon_click):

	idx = lat_lon.index((lon_click[i],lat_click[i]))

	#Profile lat/lon

	RF_lat_profile.append(lats[idx])
	RF_lon_profile.append(lons[idx])
	
	#Profile Data

	RF_data_profile_Pds.append(RF_stacking_Pds[idx])
	RF_data_profile_Pds_std.append(RF_stacking_Pds_std[idx])

	#350 LVZ
	
	RF_DEPTH_mean_LVZ_profile_Pds.append(RF_DEPTH_mean_LVZ_Pds[idx])
	RF_DEPTH_std_LVZ_profile_Pds.append(RF_DEPTH_std_LVZ_Pds[idx])

	#P410s
	
	RF_DEPTH_mean_1_profile_Pds.append(RF_DEPTH_mean_1_Pds[idx])
	RF_DEPTH_std_1_profile_Pds.append(RF_DEPTH_std_1_Pds[idx])

	#P520s
	
	RF_DEPTH_mean_520_profile_Pds.append(RF_DEPTH_mean_520_Pds[idx])
	RF_DEPTH_std_520_profile_Pds.append(RF_DEPTH_std_520_Pds[idx])

	#P660s
	
	RF_DEPTH_mean_2_profile_Pds.append(RF_DEPTH_mean_2_Pds[idx])
	RF_DEPTH_std_2_profile_Pds.append(RF_DEPTH_std_2_Pds[idx])

	#700 LVZ
	
	RF_DEPTH_mean_LVZ_700_profile_Pds.append(RF_DEPTH_mean_LVZ_700_Pds[idx])
	RF_DEPTH_std_LVZ_700_profile_Pds.append(RF_DEPTH_std_LVZ_700_Pds[idx])
	
	#MTZ Pds

	RF_DEPTH_mtz_thickness_profile_Pds.append(RF_DEPTH_mtz_thickness_Pds[idx])
	RF_DEPTH_mtz_thickness_profile_Pds_std.append(RF_DEPTH_mtz_thickness_Pds_std[idx])

	#Bootstrap Data Receiver Functions Pds
	
	RF_stacking_Pds_BOOTSTRAP_profile.append(RF_stacking_Pds_BOOTSTRAP[idx])
 
	#Bootstrap Data Mean 350 LVZ
	
	RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds_profile.append(RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds[idx])

	#Bootstrap Data Mean P410s
	
	RF_BOOTSTRAP_DEPTH_mean_1_Pds_profile.append(RF_BOOTSTRAP_DEPTH_mean_1_Pds[idx])

	#Bootstrap Data Mean P520s
	
	RF_BOOTSTRAP_DEPTH_mean_520_Pds_profile.append(RF_BOOTSTRAP_DEPTH_mean_520_Pds[idx])

	#Bootstrap Data Mean P660s

	RF_BOOTSTRAP_DEPTH_mean_2_Pds_profile.append(RF_BOOTSTRAP_DEPTH_mean_2_Pds[idx])

	#Bootstrap Data Mean 700 LVZ
	
	RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds_profile.append(RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds[idx])


print('Plotting the Final Figure')

#Cross section figure

fig = plt.figure(figsize=(40, 10))

fig.suptitle('Pds Bootstrapping points')


gs = gridspec.GridSpec(3, 8)
gs.update(wspace=0.5, hspace=0.75)


#Inset image
newax = fig.add_axes([0.8, 0.8, 0.15,0.15],anchor='NE', zorder=30,projection=ccrs.Mercator(central_longitude=PROJECT_LON, globe=None))

newax.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
newax.yaxis.set_ticks_position('both')
newax.xaxis.set_ticks_position('both')

newax.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE,4), crs=ccrs.PlateCarree())
newax.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE,4), crs=ccrs.PlateCarree())
newax.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True)

newax.grid(True,which='major',color='gray',linewidth=1,linestyle='None')

reader_1_SHP = Reader(BOUNDARY_1_SHP)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
newax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=3)

reader_2_SHP = Reader(BOUNDARY_2_SHP)
shape_2_SHP = list(reader_2_SHP.geometries())
plot_shape_2_SHP = cfeature.ShapelyFeature(shape_2_SHP, ccrs.PlateCarree())
newax.add_feature(plot_shape_2_SHP, facecolor='none', edgecolor='k',linewidth=1)

n=40
x = 0.5
upper = plt.cm.seismic(np.linspace(0, x, n)[::-1])
white = plt.cm.seismic(np.ones(20)*0.5)
lower= plt.cm.seismic(np.linspace(1-x, 1, n)[::-1])
colors = np.vstack((lower, white,upper))
tmap = matplotlib.colors.LinearSegmentedColormap.from_list('map_white', colors)

colormap_segmentation = INTER_DEPTH/10

bounds = np.arange(200, 300+colormap_segmentation, colormap_segmentation)
norm_410 = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=colormap.N)

for i,j in enumerate(lons):
	if math.isnan(RF_DEPTH_mtz_thickness_Pds[i]) == False:
		circulo_410 = Circle(radius=DIST_GRID_PP,xy=(lons[i],lats[i]),color=tmap(norm_410(RF_DEPTH_mtz_thickness_Pds[i])), ec='None',transform=ccrs.Geodetic(),zorder=3)
		newax.add_patch(circulo_410)
		circulo_410.pickable()
		circulo_410.set_picker(True)

for i,j in enumerate(lon_click):
	retangulo_PICK = Circle(radius=DIST_GRID_PP*2,xy=(lon_click[i],lat_click[i]),color='r', ec='k',linewidth=1,transform=ccrs.Geodetic(),zorder=5)
	newax.add_patch(retangulo_PICK)

newax.plot(sta_long,sta_lat, '^',markersize=5,markeredgecolor='k',markerfacecolor='grey',transform=ccrs.PlateCarree())
newax.axis('off')


#Figure Pds
for _i, _j in enumerate(RF_data_profile_Pds):
		pds_grid = fig.add_subplot(gs[0:3, _i*2:_i*2+1])

		pds_grid_410_660 = fig.add_subplot(gs[0:3, _i*2+1])

		factor_Pds = 1

		majorLocatorY = MultipleLocator(50)
		minorLocatorY = MultipleLocator(10)

	
		x_data_Pds= []
		for x,c in enumerate(RF_stacking_Pds_BOOTSTRAP_profile[_i]):
			RF_data_factor_Pds_bootstrap = [l for k, l in enumerate(c)]
			x_data_Pds.append(RF_data_factor_Pds_bootstrap)
			pds_grid.plot(RF_data_factor_Pds_bootstrap,camadas_terra_10_km,'silver',linewidth=0.1, zorder=10)

		min_x = [min(a) for a in zip(*x_data_Pds)]
		max_x = [max(a) for a in zip(*x_data_Pds)]
		pds_grid.fill_betweenx(y=camadas_terra_10_km,x1=min_x, x2=max_x, facecolor='whitesmoke',alpha=0.3, interpolate=True, zorder=5)
		
		if math.isnan(RF_DEPTH_mean_1_profile_Pds[i]) == False:
			pds_grid.text(-0.0095,RF_DEPTH_mean_1_profile_Pds[_i],str(round(RF_DEPTH_mean_1_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_std_1_profile_Pds[_i])),zorder=40,fontsize=9, fontweight='bold',ha='left',bbox={'facecolor':'white','edgecolor':'none','pad':1})
		
		if math.isnan(RF_DEPTH_mean_520_profile_Pds[i]) == False and math.isnan(RF_DEPTH_std_520_profile_Pds[i]) == False:
			pds_grid.text(-0.0095,RF_DEPTH_mean_520_profile_Pds[_i],str(round(RF_DEPTH_mean_520_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_std_520_profile_Pds[_i])),zorder=41,fontsize=9, fontweight='bold',ha='left',bbox={'facecolor':'white','edgecolor':'none','pad':1})

		if math.isnan(RF_DEPTH_mean_2_profile_Pds[i]) == False:
			pds_grid.text(-0.0095,RF_DEPTH_mean_2_profile_Pds[_i],str(round(RF_DEPTH_mean_2_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_std_2_profile_Pds[_i])),zorder=42,fontsize=9, fontweight='bold',ha='left',bbox={'facecolor':'white','edgecolor':'none','pad':1})

		if math.isnan(RF_DEPTH_mean_LVZ_profile_Pds[i]) == False and math.isnan(RF_DEPTH_std_LVZ_profile_Pds[i]) == False:
			pds_grid.text(0.001,RF_DEPTH_mean_LVZ_profile_Pds[_i],str(round(RF_DEPTH_mean_LVZ_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_std_LVZ_profile_Pds[_i])),zorder=42,fontsize=9, fontweight='bold',ha='left',bbox={'facecolor':'white','edgecolor':'none','pad':1})

		if math.isnan(RF_DEPTH_mean_LVZ_700_profile_Pds[i]) == False and math.isnan(RF_DEPTH_std_LVZ_700_profile_Pds[i]) == False:
			pds_grid.text(0.001,RF_DEPTH_mean_LVZ_700_profile_Pds[_i],str(round(RF_DEPTH_mean_LVZ_700_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_std_LVZ_700_profile_Pds[_i])),zorder=42,fontsize=9, fontweight='bold',ha='left',bbox={'facecolor':'white','edgecolor':'none','pad':1})

		RF_data_factor_Pds = [l for k, l in enumerate(_j)]
		pds_grid.plot(RF_data_factor_Pds,camadas_terra_10_km,'k',linewidth=2, zorder=30)

		std_mais = []
		std_menos = []
		for w,t in enumerate(RF_data_factor_Pds):
			std_mais.append(t+RF_data_profile_Pds_std[i][w])
			std_menos.append(t-RF_data_profile_Pds_std[i][w])

		pds_grid.plot(std_mais,camadas_terra_10_km,color='gray',linewidth=1,linestyle='dashed', zorder=30)
		pds_grid.plot(std_menos,camadas_terra_10_km,color='gray',linewidth=1,linestyle='dashed', zorder=30)

		pds_grid.yaxis.set_ticks_position('both')
		pds_grid.yaxis.set_major_locator(majorLocatorY)
		pds_grid.yaxis.set_minor_locator(minorLocatorY)
		pds_grid.grid(True,which='major',linestyle='None')

		pds_grid.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,0,where=np.array(RF_data_factor_Pds)>=0,alpha=0.5, facecolor='dimgrey',interpolate=True, zorder=19)
		pds_grid.fill_betweenx(camadas_terra_10_km,RF_data_factor_Pds,0,where=np.array(RF_data_factor_Pds)<=0,alpha=0.5, facecolor='lightgrey', interpolate=True, zorder=20)
		pds_grid.set_xlim(-0.01,0.01)
		pds_grid.xaxis.label.set_size(24)
		pds_grid.yaxis.label.set_size(24)

	
		pds_grid.set_title('Lat = '+str(round(RF_lat_profile[_i],1))+' - Lon = '+str(round(RF_lon_profile[_i],1)))
		pds_grid.set_ylim(799,301)

		if _i == 0:
			pds_grid.set_ylabel('Depth (km)',fontsize=15)
			pds_grid.yaxis.set_label_position("left")
		if _i != 0:
			pds_grid.axes.axes.yaxis.set_ticklabels([])

		#### Plot 350 LVZ ####

		pds_grid_410_660.hist(RF_BOOTSTRAP_DEPTH_mean_LVZ_Pds_profile[_i],bins=10,orientation='horizontal',color='k')

		#### Plot Depth 410 Pds ####

		pds_grid_410_660.hist(RF_BOOTSTRAP_DEPTH_mean_1_Pds_profile[_i],bins=10,orientation='horizontal',color='k')

		#### Plot Depth 520 Pds ####

		pds_grid_410_660.hist(RF_BOOTSTRAP_DEPTH_mean_520_Pds_profile[_i],bins=10,orientation='horizontal',color='k')

		#### Plot Depth 660 Pds ####

		pds_grid_410_660.hist(RF_BOOTSTRAP_DEPTH_mean_2_Pds_profile[_i],bins=10,orientation='horizontal',color='k')

		#### Plot 700 LVZ ####

		pds_grid_410_660.hist(RF_BOOTSTRAP_DEPTH_mean_LVZ_700_Pds_profile[_i],bins=10,orientation='horizontal',color='k')

		pds_grid_410_660.yaxis.set_ticks_position('both')
		pds_grid_410_660.yaxis.set_ticks_position('both')
		pds_grid_410_660.yaxis.set_major_locator(majorLocatorY)
		pds_grid_410_660.yaxis.set_minor_locator(minorLocatorY)
		pds_grid_410_660.grid(True,which='major',linestyle='None')
		pds_grid_410_660.set_xlabel('Population')
		pds_grid_410_660.set_xlim(0,100)
		pds_grid_410_660.set_ylim(799,301)


		if _i != 3:
			pds_grid_410_660.axes.axes.yaxis.set_ticklabels([])

		if _i == 3:
			pds_grid_410_660.set_ylabel('Depth (km)',fontsize=15)
			pds_grid_410_660.yaxis.set_label_position("right")
			pds_grid_410_660.tick_params(labelright=True,labelleft=False,labelsize=10)


		pds_grid_410_660.text(5,790,' MTZ = '+str(round(RF_DEPTH_mtz_thickness_profile_Pds[_i]))+'±'+str(round(RF_DEPTH_mtz_thickness_profile_Pds_std[_i])),zorder=40,fontsize=8, fontweight='bold',bbox={'facecolor':'white','edgecolor':'none','pad':1})		

plt.show()

fig.savefig(RESULTS_FOLDER+'SELECTED_BINNED_DATA_CROSS_SECTION_Pds_Ppds_bootstrap.'+EXT_FIG,dpi=DPI_FIG)

print('Ending the Cross section CODE')