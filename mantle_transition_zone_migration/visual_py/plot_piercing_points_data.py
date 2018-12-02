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







from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,STA_DIR,MIN_AMP_PDS_PPDS,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,
					PP_FIGURE,EXT_FIG,DPI_FIG,DIST_GRID_PP_MED,DIST_GRID_PP,NUMBER_STA_PER_BIN,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,BOOTSTRAP_DEPTH_ESTIMATION,GAMMA,COLORMAP_STD,COLORMAP_VEL
				   )

print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)
print('\n')


for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 410:
		Vp_depth_1 = j[3]
		Vs_depth_1 = j[5]
		
for i,j in enumerate(model_10_km.model.s_mod.v_mod.layers):
	if j[1] == 660:
		Vp_depth_2 = j[3]
		Vs_depth_2 = j[5]

print('410 km earth model Vp : '+str(Vp_depth_1))
print('410 km earth model Vs : '+str(Vs_depth_1))
print('660 km earth model Vp : '+str(Vp_depth_2))
print('660 km earth model Vs : '+str(Vs_depth_2))
print('\n')


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

dist_med_camada_terra = [abs(c - ((410+660)/2)) for x,c in enumerate(camadas_terra_10_km)]

DEPTH_MED = camadas_terra_10_km[dist_med_camada_terra.index(min(dist_med_camada_terra))]

print('Importing Pds piercing points to each PHASE')
print('\n')

PHASES = 'P410s','P'+"{0:.0f}".format(DEPTH_MED)+'s','P660s'

print('Importing Pds Piercing Points for '+PHASES[0])
print('\n')

filename_1 = PP_DIR+'PP_'+PHASES[0]+'_dic.json'

PP_1_dic = json.load(open(filename_1))

PP_dist_1 = []
PP_time_1 = []
PP_lat_1 = []
PP_lon_1 = []
PP_depth_1 = [] 

for i,j in enumerate(PP_1_dic):
	PP_dist_1.append(j['dist'][0])
	PP_time_1.append(j['time'][0])
	PP_lat_1.append(j['lat'][0])
	PP_lon_1.append(j['lon'][0])
	PP_depth_1.append(j['depth'][0])

print('Importing Pds Piercing Points for '+PHASES[1])
print('\n')

filename_med = PP_DIR+'PP_'+PHASES[1]+'_dic.json'

PP_med_dic = json.load(open(filename_med))

PP_dist_med = []
PP_time_med = []
PP_lat_med = []
PP_lon_med = []
PP_depth_med = [] 

for i,j in enumerate(PP_med_dic):
	PP_dist_med.append(j['dist'][0])
	PP_time_med.append(j['time'][0])
	PP_lat_med.append(j['lat'][0])
	PP_lon_med.append(j['lon'][0])
	PP_depth_med.append(j['depth'][0])

print('Importing Pds Piercing Points for '+PHASES[2])
print('\n')

filename_2 = PP_DIR+'PP_'+PHASES[2]+'_dic.json'

PP_2_dic = json.load(open(filename_2))

PP_dist_2 = []
PP_time_2 = []
PP_lat_2 = []
PP_lon_2 = []
PP_depth_2 = [] 

for i,j in enumerate(PP_2_dic):
	PP_dist_2.append(j['dist'][0])
	PP_time_2.append(j['time'][0])
	PP_lat_2.append(j['lat'][0])
	PP_lon_2.append(j['lon'][0])
	PP_depth_2.append(j['depth'][0])

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

print('Pds Piercing Points - '+"{0:.0f}".format(DEPTH_MED))
print('\n')

pp_med_lat  = [[]]*len(PP_lon_med)
pp_med_long  = [[]]*len(PP_lon_med)


for i,j in enumerate(PP_lon_med):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_med[i][k] == DEPTH_MED:
				pp_med_lat[i] = PP_lat_med[i][k]
				pp_med_long[i] = l

pp_med_lat = [i for i in pp_med_lat if type(i) == float ]
pp_med_long = [i for i in pp_med_long if type(i) == float ]

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

print('Importing Ppds piercing points to each PHASE')
print('\n')

PHASES_Ppds = 'PPv410s','PPv'+"{0:.0f}".format(DEPTH_MED)+'s','PPv660s'

print('Importing Ppds Piercing Points '+PHASES_Ppds[0])
print('\n')

filename_1_Ppds = PP_DIR+'PP_'+PHASES_Ppds[0]+'_dic.json'

PP_1_dic_Ppds = json.load(open(filename_1_Ppds))

PP_dist_1_Ppds = []
PP_time_1_Ppds = []
PP_lat_1_Ppds = []
PP_lon_1_Ppds = []
PP_depth_1_Ppds = [] 
PP_1_number = [] 
for i,j in enumerate(PP_1_dic_Ppds):
	PP_dist_1_Ppds.append(j['dist'][0])
	PP_time_1_Ppds.append(j['time'][0])
	PP_lat_1_Ppds.append(j['lat'][0])
	PP_lon_1_Ppds.append(j['lon'][0])
	PP_depth_1_Ppds.append(j['depth'][0])
	PP_1_number.append(j['number'][0])


print('Importing Ppds Piercing Points '+PHASES_Ppds[1])
print('\n')

filename_med_Ppds = PP_DIR+'PP_'+PHASES_Ppds[1]+'_dic.json'

PP_med_dic_Ppds = json.load(open(filename_med_Ppds))

PP_dist_med_Ppds = []
PP_time_med_Ppds = []
PP_lat_med_Ppds = []
PP_lon_med_Ppds = []
PP_depth_med_Ppds = [] 
PP_med_number = [] 

for i,j in enumerate(PP_med_dic_Ppds):
	PP_dist_med_Ppds.append(j['dist'][0])
	PP_time_med_Ppds.append(j['time'][0])
	PP_lat_med_Ppds.append(j['lat'][0])
	PP_lon_med_Ppds.append(j['lon'][0])
	PP_depth_med_Ppds.append(j['depth'][0])
	PP_med_number.append(j['number'][0])


print('Importing Ppds Piercing Points '+PHASES_Ppds[2])
print('\n')

filename_2_Ppds = PP_DIR+'PP_'+PHASES_Ppds[2]+'_dic.json'

PP_2_dic_Ppds = json.load(open(filename_2_Ppds))

PP_dist_2_Ppds = []
PP_time_2_Ppds = []
PP_lat_2_Ppds = []
PP_lon_2_Ppds = []
PP_depth_2_Ppds = [] 
PP_2_number = [] 

for i,j in enumerate(PP_2_dic_Ppds):
	PP_dist_2_Ppds.append(j['dist'][0])
	PP_time_2_Ppds.append(j['time'][0])
	PP_lat_2_Ppds.append(j['lat'][0])
	PP_lon_2_Ppds.append(j['lon'][0])
	PP_depth_2_Ppds.append(j['depth'][0])
	PP_2_number.append(j['number'][0])

print('PPv410s Piercing Points')
print('\n')

pp_1_lat_Ppds  = [[]]*len(PP_lon_1_Ppds)
pp_1_long_Ppds  = [[]]*len(PP_lon_1_Ppds)


for i,j in enumerate(PP_lon_1_Ppds):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_1_Ppds[i][k] == 410:
				pp_1_lat_Ppds[i] = PP_lat_1_Ppds[i][k]
				pp_1_long_Ppds[i] = l

pp_1_lat_Ppds = [i for i in pp_1_lat_Ppds if type(i) == float ]
pp_1_long_Ppds = [i for i in pp_1_long_Ppds if type(i) == float ]


print('Ppds Piercing Points - '+"{0:.0f}".format(DEPTH_MED))
print('\n')

pp_med_lat_Ppds  = [[]]*len(PP_lon_med_Ppds)
pp_med_long_Ppds  = [[]]*len(PP_lon_med_Ppds)


for i,j in enumerate(PP_lon_med_Ppds):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE<= l <= URCRNRLON_LARGE and PP_depth_med_Ppds[i][k] == DEPTH_MED:
				pp_med_lat_Ppds[i] = PP_lat_med_Ppds[i][k]
				pp_med_long_Ppds[i] = l

pp_med_lat_Ppds = [i for i in pp_med_lat_Ppds if type(i) == float ]
pp_med_long_Ppds = [i for i in pp_med_long_Ppds if type(i) == float ]

print('PPv660s Piercing Points')
print('\n')


pp_2_lat_Ppds  = [[]]*len(PP_lon_2_Ppds)
pp_2_long_Ppds  = [[]]*len(PP_lon_2_Ppds)
pp_2_depth_Ppds  = [[]]*len(PP_lon_2_Ppds)


for i,j in enumerate(PP_lon_2_Ppds):
	for k,l in enumerate(j):
		if LLCRNRLON_LARGE <= l <= URCRNRLON_LARGE and PP_depth_2_Ppds[i][k] == 660:
				pp_2_lat_Ppds[i] = PP_lat_2_Ppds[i][k]
				pp_2_long_Ppds[i] = l

pp_2_lat_Ppds = [i for i in pp_2_lat_Ppds if type(i) == float ]
pp_2_long_Ppds = [i for i in pp_2_long_Ppds if type(i) == float ]



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


dist_pp_grid_min = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_min[i] = [np.sqrt((j - pp_1_long_Ppds[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_1_lat_Ppds)]

dist_pp_grid_med = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_med[i] = [np.sqrt((j - pp_med_long_Ppds[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_med_lat_Ppds)]
    
dist_pp_grid_max = [[]]*len(grdx)
for i,j in enumerate(grdx):
    dist_pp_grid_max[i] = [np.sqrt((j - pp_2_long_Ppds[k])**2 + (grdy[i] - l)**2) for k,l in enumerate(pp_2_lat_Ppds)]


grid_sel_min = []
grid_sel_min_data = []
for i,j in enumerate(dist_pp_grid_min):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[NUMBER_PP_PER_BIN]] < DIST_GRID_PP:
        grid_sel_min.append((grdx[i],grdy[i]))


grid_sel_med = []
grid_sel_med_data = []
for i,j in enumerate(dist_pp_grid_med):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[NUMBER_PP_PER_BIN]] < DIST_GRID_PP:
        grid_sel_med.append((grdx[i],grdy[i]))

        
grid_sel_max = []
grid_sel_min_data = []

for i,j in enumerate(dist_pp_grid_max):
    vect_j = np.array(j) 
    indices = vect_j.argsort()
    if vect_j[indices[NUMBER_PP_PER_BIN]] < DIST_GRID_PP:
        grid_sel_max.append((grdx[i],grdy[i]))


grid_sel = grid_sel_min+grid_sel_med+grid_sel_max


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

###################################################################################################################

print('Plotting: Figure Radial receiver functions stacked in ray parameter')
print('\n')

time_PP_wave = []
time_P410s_wave = []
time_Pp410s_wave = []
time_P660s_wave = []
time_Pp660s_wave = []
for i,j in enumerate(event_depth):
	model_RF = TauPyModel(model="iasp91")
	arrivalsP = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["P"])
	arrP = arrivalsP[0]

	arrivalsP410s = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["P410s"])
	arrP410s = arrivalsP410s[0]
    
	arrivalsPp410s = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["PPv410s"])
	arrPp410s = arrivalsPp410s[0]
    
	arrivalsP660s = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["P660s"])
	arrP660s = arrivalsP660s[0]
    
	arrivalsPp660s = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["PPv660s"])
	arrPp660s = arrivalsPp660s[0]

	arrivalsPP = model_RF.get_travel_times(source_depth_in_km=event_depth[i], distance_in_degree=event_gcarc[i], phase_list=["PP"])
	arrPP = arrivalsPP[0]
	
	time_PP_wave.append(arrPP.time - arrP.time)

	time_P410s_wave.append(arrP410s.time - arrP.time)
    
	time_Pp410s_wave.append(arrPp410s.time - arrP.time)
    
	time_P660s_wave.append(arrP660s.time - arrP.time)
    
	time_Pp660s_wave.append(arrPp660s.time - arrP.time)

RF_orglisl = np.argsort(event_ray)[::-1] 

time = sta_time[0]

FR = []
GCARC = []
RP = [] 
time_PP_wave_corrected = []
time_P410s_wave_corrected = []
time_Pp410s_wave_corrected = []
time_P660s_wave_corrected = []
time_Pp660s_wave_corrected = []
for i,j in enumerate(RF_orglisl):
	FR.append(sta_data[j])
	RP.append(event_ray[j])
	time_PP_wave_corrected.append(time_PP_wave[j])
	time_P410s_wave_corrected.append(time_P410s_wave[j])
	time_Pp410s_wave_corrected.append(time_Pp410s_wave[j])
	time_P660s_wave_corrected.append(time_P660s_wave[j])
	time_Pp660s_wave_corrected.append(time_Pp660s_wave[j])

Z = np.array(FR)

########################################

fig, ax = plt.subplots(1, 1, figsize=(5, 20))

majorLocatorX = MultipleLocator(300)
majorLocatorY = MultipleLocator(100)
minorLocatorY = MultipleLocator(20)
minorLocatorX = MultipleLocator(10)


#Sismograma sem filtro PP
v=0.008
im = ax.imshow(Z.T, interpolation='bicubic', cmap=cm.viridis,
                origin='upper', aspect='auto',
                vmax=v, vmin=-v)
for i,j in enumerate(time_PP_wave_corrected):
    ax.plot(i,(j-30)*10,'.k',markersize=2)
    ax.plot(i,j*10,'^r',markersize=2)
    ax.plot(i,(j+30)*10,'.k',markersize=2)
    
ax.set_ylim(2500,0)
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('RF Data')

ax.grid(True)
ax.set_yticklabels(["{0:.0f}".format(time[i]) for i in np.arange(-100,len(time),100)])
ax.set_xticklabels(["{0:.1f}".format(RP[i]*100) for i in np.arange(0,len(RP),100)])
plt.show()

########################################

fig, ax = plt.subplots(1, 1, figsize=(5, 20))

majorLocatorX = MultipleLocator(300)
majorLocatorY = MultipleLocator(100)
minorLocatorY = MultipleLocator(20)
minorLocatorX = MultipleLocator(10)


#Sismograma sem filtro PP
v=0.008
im = ax.imshow(Z.T, interpolation='bicubic', cmap=cm.viridis,
                origin='upper', aspect='auto',
                vmax=v, vmin=-v)

for i,j in enumerate(time_P410s_wave_corrected):
    ax.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)

for i,j in enumerate(time_Pp410s_wave_corrected):
    ax.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)
    
for i,j in enumerate(time_P660s_wave_corrected):
    ax.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)
    
for i,j in enumerate(time_Pp660s_wave_corrected):
    ax.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)
    
    
ax.set_ylim(2500,0)
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('RF Data')

ax.grid(True)
ax.set_yticklabels(["{0:.0f}".format(time[i]) for i in np.arange(-100,len(time),100)])
ax.set_xticklabels(["{0:.1f}".format(RP[i]*100) for i in np.arange(0,len(RP),100)])
plt.show()
'''
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

###################################################################################################################

print('Plotting: Figure 410 and 660 Pds Piercing Points')
print('\n')

#Pds 410

PP_410_lat = []
PP_410_lon = []
PP_410_depth = []
for i,j in enumerate(PP_depth_1):
	for k,l in enumerate(j):
		if LLCRNRLAT_SMALL <= PP_lat_1[i][k] <= URCRNRLAT_SMALL and LLCRNRLON_SMALL <= PP_lon_1[i][k] <= URCRNRLON_SMALL:
			PP_410_lat.append(PP_lat_1[i][k])
			PP_410_lon.append(PP_lon_1[i][k])
			PP_410_depth.append(-PP_depth_1[i][k])

#Pds 660

PP_660_lat = []
PP_660_lon = []
PP_660_depth = []
for i,j in enumerate(PP_depth_2):
	for k,l in enumerate(j):
		if LLCRNRLAT_SMALL <= PP_lat_2[i][k] <= URCRNRLAT_SMALL and LLCRNRLON_SMALL <= PP_lon_2[i][k] <= URCRNRLON_SMALL:
			PP_660_lat.append(PP_lat_2[i][k])
			PP_660_lon.append(PP_lon_2[i][k])
			PP_660_depth.append(-PP_depth_2[i][k])

#Pds Middle Layer

MED_PP_lat = []
MED_PP_lon = []
MED_PP_depth = []
for i,j in enumerate(PP_depth_2):
	for k,l in enumerate(j):
		if LLCRNRLAT_SMALL <= PP_lat_med[i][k] <= URCRNRLAT_SMALL and LLCRNRLON_SMALL <= PP_lon_med[i][k] <= URCRNRLON_SMALL:
			MED_PP_lat.append(PP_lat_med[i][k])
			MED_PP_lon.append(PP_lon_med[i][k])
			MED_PP_depth.append(-PP_depth_med[i][k])


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
'''