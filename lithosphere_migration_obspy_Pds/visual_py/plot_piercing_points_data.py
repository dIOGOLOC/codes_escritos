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
import shapefile
from fatiando import gridder, utils
import scipy.io
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import json
import random
from matplotlib.colors import Normalize
from numpy import ma
from matplotlib import cbook
import collections
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.gridspec as gridspec
from shapely.geometry import Polygon, MultiPoint, Point
import shapefile
from mpl_toolkits.axes_grid1.inset_locator import inset_axes








from parameters_py.mgconfig import (
					RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					NUMBER_PP_PER_BIN,LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG,NUMBER_STA_PER_BIN,OUTPUT_DIR,RF_FREQUENCY,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,COLORMAP_STD,COLORMAP_VEL,DEPTH_TARGET
				   )

print('Starting Receiver Functions migration code to estimate the depths of the Mantle discontinuities')
print('\n')

print('Importing depths and times of Pds conversion  dataset')
print('\n')

print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)


for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 410:
		Vp_depth_1 = j[3]
		Vs_depth_1 = j[5]

print('410 km earth model Vp : '+str(Vp_depth_1))
print('410 km earth model Vs : '+str(Vs_depth_1))
print('----')
#====================================================

for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 520:
		Vp_depth_520 = j[3]
		Vs_depth_520 = j[5]

print('520 km earth model Vp : '+str(Vp_depth_520))
print('520 km earth model Vs : '+str(Vs_depth_520))
print('----')
#====================================================
		
for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 660:
		Vp_depth_2 = j[3]
		Vs_depth_2 = j[5]

print('660 km earth model Vp : '+str(Vp_depth_2))
print('660 km earth model Vs : '+str(Vs_depth_2))
print('----')

#====================================================
print('Target Depth = '+str(DEPTH_TARGET))

for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == DEPTH_TARGET:
		Vp_depth_DEPTH_TARGET = j[3]
		Vs_depth_DEPTH_TARGET = j[5]

print('Target Depth earth model Vp : '+str(Vp_depth_DEPTH_TARGET))
print('Target Depth earth model Vs : '+str(Vs_depth_DEPTH_TARGET))
#====================================================

print('\n')

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Stations'+'/'

print('Looking for Receiver Functions data in JSON file in '+STA_DIR)
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

print('Creating the Earth layered model')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

print('Importing Pds piercing points to each PHASE')
print('\n')

PHASES = 'P410s','P'+str(DEPTH_TARGET)+'s','P660s'

print('Importing Pds Piercing Points for '+PHASES[0])
print('\n')

PP_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Piercing_Points'+'/'


filename_1 = PP_DIR+'PP_'+PHASES[0]+'_dic.json'

PP_1_dic = json.load(open(filename_1))


PP_time_1 = PP_1_dic['time']
PP_lat_1 = PP_1_dic['lat']
PP_lon_1 = PP_1_dic['lon']
PP_depth_1 = PP_1_dic['depth']


print('Importing Pds Piercing Points for '+PHASES[1])
print('\n')

filename_med = PP_DIR+'PP_'+PHASES[1]+'_dic.json'

PP_med_dic = json.load(open(filename_med))

PP_time_med = PP_med_dic['time']
PP_lat_med = PP_med_dic['lat']
PP_lon_med = PP_med_dic['lon']
PP_depth_med = PP_med_dic['depth']


print('Importing Pds Piercing Points for '+PHASES[2])
print('\n')

filename_2 = PP_DIR+'PP_'+PHASES[2]+'_dic.json'

PP_2_dic = json.load(open(filename_2))

PP_time_2 = PP_2_dic['time']
PP_lat_2 = PP_2_dic['lat']
PP_lon_2 = PP_2_dic['lon']
PP_depth_2 = PP_2_dic['depth'] 

print('P410s Piercing Points')
print('\n')

pp_1_lat  = [[]]*len(PP_lon_1)
pp_1_long  = [[]]*len(PP_lon_1)


for i,j in enumerate(PP_lon_1):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_1[i][k] == 410:
				pp_1_lat[i] = PP_lat_1[i][k]
				pp_1_long[i] = l

pp_1_lat = [i for i in pp_1_lat if type(i) == float ]
pp_1_long = [i for i in pp_1_long if type(i) == float ]

print('P'+str(DEPTH_TARGET)+'s Piercing Points')
print('\n')

pp_med_lat  = [[]]*len(PP_lon_med)
pp_med_long  = [[]]*len(PP_lon_med)
pp_time_DEPTH_TARGET  = [[]]*len(PP_lon_med)

for i,j in enumerate(PP_lon_med):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_med[i][k] == DEPTH_TARGET:
				pp_med_lat[i] = PP_lat_med[i][k]
				pp_med_long[i] = l
				pp_time_DEPTH_TARGET[i] = PP_time_med[i][-1] - PP_time_med[i][k]

pp_med_lat = [i for i in pp_med_lat if type(i) == float ]
pp_med_long = [i for i in pp_med_long if type(i) == float ]
pp_time_DEPTH_TARGET = [i for i in pp_time_DEPTH_TARGET if type(i) == float ]

print('P660s Piercing Points')
print('\n')

pp_2_lat  = [[]]*len(PP_lon_2)
pp_2_long  = [[]]*len(PP_lon_2)


for i,j in enumerate(PP_lon_2):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_2[i][k] == 660:
			if PP_lon_2[i][k]  != [] and PP_lat_2[i][k] != []:
				pp_2_lat[i] = PP_lat_2[i][k]
				pp_2_long[i] = l

pp_2_lat = [i for i in pp_2_lat if type(i) == float ]
pp_2_long = [i for i in pp_2_long if type(i) == float ]


print('Calculating MEAN FIRST FRESNEL ZONE RADIUS')
print('DEPTH TARGET MEAN TIME: '+str(np.mean(pp_time_DEPTH_TARGET)))

FRESNEL_ZONE_RADIUS_km = (Vp_depth_DEPTH_TARGET/2)* np.sqrt(np.mean(pp_time_DEPTH_TARGET) / RF_FREQUENCY)
FRESNEL_ZONE_RADIUS = kilometer2degrees(FRESNEL_ZONE_RADIUS_km)
print('MEAN FIRST FRESNEL ZONE RADIUS : '+str(FRESNEL_ZONE_RADIUS))
print('\n')

print('Creating GRID POINTS')
print('\n')

area = (LLCRNRLON_SMALL,URCRNRLON_SMALL, LLCRNRLAT_SMALL, URCRNRLAT_SMALL)

shape = (abs(abs(URCRNRLON_SMALL) - abs(LLCRNRLON_SMALL))*GRID_PP_MULT, abs(abs(URCRNRLAT_SMALL) - abs(LLCRNRLAT_SMALL))*GRID_PP_MULT)

grdx, grdy = gridder.regular(area, shape)

if FILTER_BY_SHAPEFILE == True:
	polygon = shapefile.Reader(SHAPEFILE_GRID) 
	polygon = polygon.shapes()  
	shpfilePoints = []
	for shape in polygon:
		shpfilePoints = shape.points 
	polygon = shpfilePoints 
	poly = Polygon(polygon)
	
	pontos = [Point(grdx[i],grdy[i]) for i,j in enumerate(grdx)]
	inside_points = np.array(MultiPoint(pontos).intersection(poly))
	
	grdx = []
	grdy = []
	for i,j in enumerate(inside_points):
		grdx.append(j[0])
		grdy.append(j[1])


print('Filtering GRID POINTS')
print('\n')

dist_pp_grid_med = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_med[i] = [np.sqrt((j - pp_med_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_med_lat)]
    dist_pp_grid_med[i] = [np.sqrt((j - pp_med_long[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_med_lat)]


grid_sel_med = []
for i,j in enumerate(dist_pp_grid_med):
	vect_j = np.array(j)
	indices = vect_j.argsort()
	if vect_j[indices[NUMBER_PP_PER_BIN]] <= FRESNEL_ZONE_RADIUS:
		grid_sel_med.append((grdx[i],grdy[i]))

grid_sel = grid_sel_med

grid_selected = set(map(tuple,grid_sel))

grid_sel_x = []
grid_sel_y = []

for i,j in enumerate(grid_selected):
    grid_sel_x.append(j[0])
    grid_sel_y.append(j[1])

grdx_1 = []
grdy_1 = []
grdx_sel_diff = []
for i,j in enumerate(grdx):
    grdx_sel_diff.append((j,grdy[i]))

grid_sel_diff = []
for i,j in enumerate(grid_selected):
	grid_sel_diff.append((j[0],j[1]))

grid_xy_diff = []
for i,j in enumerate(grdx_sel_diff):
	if j not in grid_sel_diff:
		grid_xy_diff.append(j)			
	
grid_x_diff = []
grid_y_diff = []
for i,j in enumerate(grid_xy_diff):
	grid_x_diff.append(j[0])
	grid_y_diff.append(j[1])

##########################################################################################################################################

print('Importing depths and times of Pds conversion  dataset')
print('\n')

PdS_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Phases'+'/'

filename_Pds = PdS_DIR+'Pds_dic.json'

PdS_Dic = json.load(open(filename_Pds))

depth_str_to_float_Pds = []
for i,j in enumerate(PdS_Dic['depth']):
	depth_str_to_float_Pds.append([float(l) for l in j])

Pds_time = PdS_Dic['time']
Pds_st_lat = PdS_Dic['st_lat']
Pds_st_lon = PdS_Dic['st_long']
Pds_ev_lat = PdS_Dic['ev_lat']
Pds_ev_lon = PdS_Dic['ev_long']
Pds_depth = depth_str_to_float_Pds

###################################################################################################################
'''
print('Plotting: Figure Radial receiver functions stacked in ray parameter')
print('\n')

time_PP_wave = []
time_P410s_wave = []
time_P660s_wave = []

for i,j in enumerate(event_depth):
	model_RF = TauPyModel(model="iasp91")
	arrivalsP = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["P"])
	arrP = arrivalsP[0]

	arrivalsP410s = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["P410s"])
	arrP410s = arrivalsP410s[0]
     
	arrivalsP660s = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["P660s"])
	arrP660s = arrivalsP660s[0]
    
	arrivalsPP = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["PP"])
	arrPP = arrivalsPP[0]
	
	time_PP_wave.append(arrPP.time - arrP.time)

	time_P410s_wave.append(arrP410s.time - arrP.time)
    
	time_P660s_wave.append(arrP660s.time - arrP.time)
    

RF_orglisl = np.argsort(event_ray)[::-1] 
#RF_orglisl = np.argsort(event_gcarc) 

time = sta_time[0]

FR = []
GC = []
RP = [] 
time_PP_wave_corrected = []
time_P410s_wave_corrected = []
time_P660s_wave_corrected = []
for i,j in enumerate(RF_orglisl):
	FR.append(sta_data[j])
	RP.append(event_ray[j])
	GC.append(event_gcarc[j])
	time_PP_wave_corrected.append(time_PP_wave[j])
	time_P410s_wave_corrected.append(time_P410s_wave[j])
	time_P660s_wave_corrected.append(time_P660s_wave[j])

Z = np.array(FR)

########################################

fig, ax = plt.subplots(1, 1, figsize=(5, 20))

majorLocatorX = MultipleLocator(100)
majorLocatorY = MultipleLocator(100)
minorLocatorY = MultipleLocator(20)
minorLocatorX = MultipleLocator(20)


#Sismograma sem filtro PP
v=0.01
im = ax.imshow(Z.T, interpolation='bicubic', cmap=cm.viridis,
                origin='upper', aspect='auto',
                vmax=v, vmin=-v)

number_RF = np.arange(0,len(RP))

number_RF = number_RF[::20]

for i,j in enumerate(time_PP_wave_corrected[::20]):
    ax.plot(number_RF[i],(j-30)*10,marker="o",color='k',markersize=5)
    ax.plot(number_RF[i],j*10,marker="o",color='r',markersize=5)
    ax.plot(number_RF[i],(j+30)*10,marker="o",color='k',markersize=5)
    

axins = inset_axes(ax,
                   width="20%",  # width = 10% of parent_bbox width
                   height="2%",  # height : 50%
                   loc='upper left',
                   bbox_to_anchor=(0.75, 0.03, 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')


ax.set_ylim(2500,0)
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('RF Data')
ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position('both')
ax.tick_params(labelleft=True,labelright=True)

ax.set_yticklabels(["{0:.0f}".format(time[i]) for i in np.arange(-100,len(time),100)])
ax.set_xticklabels(["{0:.1f}".format(RP[i]*100) for i in np.arange(0,720,100)])
plt.show()

########################################

fig, ax = plt.subplots(1, 1, figsize=(5, 20))

majorLocatorX = MultipleLocator(100)
majorLocatorY = MultipleLocator(100)
minorLocatorY = MultipleLocator(20)
minorLocatorX = MultipleLocator(20)


#Sismograma sem filtro PP
v=0.01
im = ax.imshow(Z.T, interpolation='bicubic', cmap=cm.viridis,
                origin='upper', aspect='auto',
                vmax=v, vmin=-v)

for i,j in enumerate(time_P410s_wave_corrected[::20]):
    ax.plot(number_RF[i],j*10,marker="_",color='k',markersize=5)


for i,j in enumerate(time_P660s_wave_corrected[::20]):
    ax.plot(number_RF[i],j*10,marker="_",color='k',markersize=5)    


axins = inset_axes(ax,
                   width="20%",  # width = 10% of parent_bbox width
                   height="2%",  # height : 50%
                   loc='upper left',
                   bbox_to_anchor=(0.75, 0.03, 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

ax.set_ylim(2500,0)
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('RF Data')

ax.yaxis.set_label_position("right")
ax.yaxis.set_ticks_position('both')
ax.tick_params(labelleft=True,labelright=True)

ax.set_yticklabels(["{0:.0f}".format(time[i]) for i in np.arange(-100,len(time),100)])
ax.set_xticklabels(["{0:.1f}".format(RP[i]*100) for i in np.arange(0,len(RP),100)])
plt.show()

###################################################################################################################

print('Plotting: Figure Distribution and usual statistics about events used in the study')
print('\n')

fig = plt.figure(figsize=(5, 10))
gs = gridspec.GridSpec(5,2)
gs.update(wspace=0.5,hspace=1)

ax1 = fig.add_subplot(gs[:3,:])

ax2 = fig.add_subplot(gs[3, 0])
ax3 = fig.add_subplot(gs[3, 1])
ax4 = fig.add_subplot(gs[4, 0])
ax5 = fig.add_subplot(gs[4, 1])



m = Basemap(resolution='l',projection='ortho',lat_0=PROJECT_LAT, lon_0=PROJECT_LON,ax=ax1)

for lon, lat in zip(sta_long,sta_lat):
    x,y = m(lon, lat)
    msize = 6
    l1, = m.plot(x, y, '^',markersize=msize,markeredgecolor='k',markerfacecolor='w')

for lon, lat in zip(evlo_0,evla_0):
	x,y = m(lon, lat)
	msize = 4
	m.plot(x, y, '*',markersize=msize,markeredgecolor='grey',markerfacecolor='whitesmoke')

for lon, lat in zip(event_lon,event_lat):
    x,y = m(lon, lat)
    msize = 6
    m.plot(x, y, '*',markersize=msize,markeredgecolor='k',markerfacecolor='dimgrey')
    
#m.tissot(PROJECT_LON, PROJECT_LAT, 30,100,zorder=10,edgecolor='dimgray',linewidth=1,facecolor='none')
#m.tissot(PROJECT_LON, PROJECT_LAT, 90,100,zorder=10,edgecolor='dimgray',linewidth=1,facecolor='none')


m.fillcontinents(color='whitesmoke',lake_color=None)
m.drawcoastlines(color='dimgray',zorder=10)
m.drawmeridians(np.arange(0, 360, 20),color='lightgrey')
m.drawparallels(np.arange(-90, 90, 10),color='lightgrey')

#############
ax2.hist(event_gcarc,bins=50,orientation='vertical',color='k')
#ax2.set_yticklabels([])
ax2.set_title("Distance (Degrees)")
ax2.set_xlim(0,100)

#############
ax3.hist(event_ray,bins=50,orientation='vertical',color='k')
#ax3.set_yticklabels([])
ax3.set_title("Ray Parameter (s/rad)")


#############
ax4.hist(event_mag,bins=50,orientation='vertical',color='k') 
#ax4.set_yticklabels([])
ax4.set_title("Magnitude (mb)")
ax4.set_xlim(0,10)

#############
ax5.hist(event_depth,bins=50,orientation='vertical',color='k')
#ax5.set_yticklabels([])
ax5.set_title("Depth (km)")


plt.show()
###################################################################################################################
print('Plotting: Figure earth model layers')

grid_camadas_x = grdx
grid_camadas_y = grdy
grid_camadas_z = camadas_terra_10_km


print('\n')

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.azim = 95
ax.elev = 10
ax.dist = 6

ax.plot(sta_long,sta_lat,'^',markersize=10,markeredgecolor='k',markerfacecolor='grey')
ax.plot(grid_sel_x,grid_sel_y,'.',markersize=2,markeredgecolor='k',markerfacecolor='k')

ax.scatter3D(grid_camadas_x, grid_camadas_y, grid_camadas_z, c='k', marker='o')
ax.set_zlim(-800,0)
ax.set_xlim(LLCRNRLON_SMALL,URCRNRLON_SMALL)
ax.set_ylim(LLCRNRLAT_SMALL,URCRNRLAT_SMALL)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Depth (km)')
plt.show()
'''

###################################################################################################################

print('Plotting: Figure 410 and 660 Pds Piercing Points')
print('\n')

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.azim = 95
ax.elev = 10
ax.dist = 8

ax.plot(sta_long,sta_lat,'^',markersize=10,markeredgecolor='k',markerfacecolor='grey')

x, y = np.meshgrid(np.arange(LLCRNRLON_SMALL,URCRNRLON_SMALL,1/GRID_PP_MULT),
                      np.arange(LLCRNRLAT_SMALL,URCRNRLAT_SMALL,1/GRID_PP_MULT))

z = np.zeros_like(x)
ax.plot_surface(x, y, z,color='None',edgecolor='k')

z1660 = np.ones_like(x)
z660 = np.array([element*-660 for element in z1660])

z1410 = np.ones_like(x)
z410 = np.array([element*-410 for element in z1410])

ax.plot_surface(x, y, z660,color='None',edgecolor='g')
ax.plot_surface(x, y, z410,color='None',edgecolor='b')

intp = cbook.simple_linear_interpolation

ax.plot3D(intp(np.array(PP_410_lon),50),intp(np.array(PP_410_lat), 50), intp(np.array(PP_410_depth), 50), c='r',alpha=0.3)
ax.plot3D(intp(np.array(PP_660_lon),50),intp(np.array(PP_660_lat), 50), intp(np.array(PP_660_depth), 50), c='g',alpha=0.3)

ax.scatter3D(PP_410_lon_POINTS,PP_410_lat_POINTS, PP_410_depth_POINTS,  c='k',marker='X',s=50)
ax.scatter3D(PP_660_lon_POINTS,PP_660_lat_POINTS, PP_660_depth_POINTS,  c='k',marker='X',s=50)

ax.set_zlim(-800,0)
ax.set_xlim(LLCRNRLON_SMALL,URCRNRLON_SMALL)
ax.set_ylim(LLCRNRLAT_SMALL,URCRNRLAT_SMALL)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Depth (km)')
plt.show()

###################################################################################################################

print('Plotting: Figure Middle Pds Piercing Points and Model earth')
print('\n')



fig = plt.figure(figsize=(30,15))
ax = Axes3D(fig)

ax.azim = 95
ax.elev = 10
ax.dist = 8

ax.plot(grid_sel_x,grid_sel_y,'.',markersize=2,markeredgecolor='k',markerfacecolor='k')
ax.plot(sta_long,sta_lat,'^',markersize=10,markeredgecolor='k',markerfacecolor='grey')


z = np.zeros_like(x)
ax.plot_surface(x, y, z,color='None',edgecolor='k')

ax.scatter3D(PP_410_lon_POINTS,PP_410_lat_POINTS, PP_410_depth_POINTS,  c='b',marker='+')
ax.scatter3D(MED_PP_lon_POINTS,MED_PP_lat_POINTS, MED_PP_depth_POINTS,  c='r',marker='+',s=50)
ax.scatter3D(PP_660_lon_POINTS,PP_660_lat_POINTS, PP_660_depth_POINTS,  c='g',marker='+')

m1 = ax.scatter3D(grid_camadas_x, grid_camadas_y, grid_camadas_z,c=grid_camadas_z,edgecolor='k',cmap='bone')


ax.set_zlim(-800,0)
ax.set_xlim(LLCRNRLON_SMALL,URCRNRLON_SMALL)
ax.set_ylim(LLCRNRLAT_SMALL,URCRNRLAT_SMALL)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Depth (km)')

fig.colorbar(m1,aspect=40,shrink=0.7)

plt.show()

###################################################################################################################

print('Plotting: Figure 410 Pds Piercing Points and Model earth')
print('\n')


fig = plt.figure(figsize=(30,15))
ax = Axes3D(fig)

ax.azim = 95
ax.elev = 10
ax.dist = 8

ax.plot(grid_sel_x,grid_sel_y,'.',markersize=2,markeredgecolor='k',markerfacecolor='k')
ax.plot(sta_long,sta_lat,'^',markersize=10,markeredgecolor='k',markerfacecolor='grey')


z = np.zeros_like(x)
#ax.plot_surface(x, y, z,color='None',edgecolor='k')

ax.scatter3D(MED_PP_lon_POINTS,MED_PP_lat_POINTS, MED_PP_depth_POINTS,  c='r',marker='X',s=50)
ax.plot3D(intp(np.array(MED_PP_lon),100),intp(np.array(MED_PP_lat), 100), intp(np.array(MED_PP_depth), 100), c='k',alpha=0.3)


m1 = ax.scatter3D(grid_camadas_x, grid_camadas_y, grid_camadas_z,c=grid_camadas_z,edgecolor='k',cmap='bone')


ax.set_zlim(-800,0)
ax.set_xlim(LLCRNRLON_SMALL,URCRNRLON_SMALL)
ax.set_ylim(LLCRNRLAT_SMALL,URCRNRLAT_SMALL)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Depth (km)')

fig.colorbar(m1,aspect=40,shrink=0.7)

plt.show()