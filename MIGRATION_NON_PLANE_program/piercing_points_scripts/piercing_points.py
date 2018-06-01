
# coding: utf-8

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees


print('Looking for receiver functions files in '+RF_DIR)
ev_list = []
ev_listS = []
for root, dirs, files in os.walk(RF_DIR):
    for datafile in files:
        if datafile.endswith(RF_EXT):
            ev_list.append(os.path.join(root, datafile))
ev_listS = sorted(ev_list)


ev = obspy.Stream()
for i,j in enumerate(ev_listS):
    ev += obspy.read(j)


event_DD = []
event_MM = []
event_YYYY = []
event_hh = []
event_mm = []
event_julday = []
event_depth = []
event_lat = []
event_long = []
event_dist = []
event_gcarc = []
event_sta = []
event_channel = []
event_ray = []
sta_lat = []
sta_long = []
sta_channel = []
sta_data = []
sta_time = []
event_starttime = []
event_endtime = []


for i,j in enumerate(ev):
    if j.stats.sac.gcarc > 30:
            event_time = (j.stats.starttime)
            event_starttime.append(j.stats.starttime)
            event_endtime.append(j.stats.endtime)
            event_DD.append("{0:02.0f}".format(event_time.day))
            event_MM.append("{0:02.0f}".format(event_time.month))
            event_YYYY.append(event_time.year)
            event_hh.append("{0:02.0f}".format(event_time.hour))
            event_mm.append("{0:02.0f}".format(event_time.minute))
            event_julday.append(event_time.julday)
            #event_depth.append(j.stats.sac.evdp)
            event_depth.append(j.stats.sac.evdp/1000) #para os dados sintéticos
            event_lat.append(j.stats.sac.evla)
            event_long.append(j.stats.sac.evlo)
            event_dist.append(j.stats.sac.dist)
            event_gcarc.append(j.stats.sac.gcarc)
            event_sta.append(j.stats.station)
            event_ray.append(j.stats.sac.user8)
            sta_lat.append(j.stats.sac.stla)
            sta_long.append(j.stats.sac.stlo)
            sta_data.append(j.data[100:2700])
            sta_time.append(j.times()[100:2700]-10)


print('Importing earth model from obspy.taup.TauPyModel')
model_THICKNESS_km = TauPyModel(model=MODEL_DIR_NPZ)


# # Calculando os caminhos dos raios para os eventos das estações BP
fase = PHASES_PIERCING_POINTS_PS


print('Phase '+fase+' calculation')
arrivalsP410s = []
for i,j in enumerate(event_depth):
	print('Event '+str(i)+' of '+str(len(event_depth))+' events')
	arrivalsP410s.append(model_THICKNESS_km.get_pierce_points_geo(
                                    source_depth_in_km=j, source_latitude_in_deg=event_lat[i], 
                                    source_longitude_in_deg=event_long[i], receiver_latitude_in_deg=sta_lat[i], 
                                    receiver_longitude_in_deg=sta_long[i], phase_list=[fase]))


dist_P410s = [[]]*len(arrivalsP410s)
depth_P410s = [[]]*len(arrivalsP410s)
time_P410s = [[]]*len(arrivalsP410s)
lat_P410s = [[]]*len(arrivalsP410s)
lon_P410s = [[]]*len(arrivalsP410s)

for i,j in enumerate(arrivalsP410s):
	time_P410s[i] = [l.pierce['time'] for k,l in enumerate(j)][0]
	depth_P410s[i] = [l.pierce['depth'] for k,l in enumerate(j)][0]
	dist_P410s[i] = [l.pierce['dist'] for k,l in enumerate(j)][0]
	lat_P410s[i] = [l.pierce['lat'] for k,l in enumerate(j)][0]
	lon_P410s[i] = [l.pierce['lon'] for k,l in enumerate(j)][0]

print('Saving .TXT files')
#distance
P410s_dist_txt = open(DIST_PP_DIR+fase+'_dist.txt', 'w')

for i,j in enumerate(dist_P410s):
    P410s_dist_txt.write(str(list(j))+'\n')
P410s_dist_txt.close()


#depth
depth_P410s_txt = open(DEPTH_PP_DIR+fase+'_depth.txt', 'w')

for i,j in enumerate(depth_P410s):
    depth_P410s_txt.write(str(list(j))+'\n')
depth_P410s_txt.close()


#time
time_P410s_txt = open(TIME_PP_DIR+fase+'_time.txt', 'w')

for i,j in enumerate(time_P410s):
    time_P410s_txt.write(str(list(j))+'\n')
time_P410s_txt.close()


#latitude
lat_P410s_txt = open(LAT_PP_DIR+fase+'_lat.txt', 'w')

for i,j in enumerate(lat_P410s):
    lat_P410s_txt.write(str(list(j))+'\n')
lat_P410s_txt.close()


#longitude
lon_P410s_txt = open(LON_PP_DIR+fase+'_lon.txt', 'w')

for i,j in enumerate(lon_P410s):
    lon_P410s_txt.write(str(list(j))+'\n')
lon_P410s_txt.close()

print("Piercing Points estimated!")

