# coding: utf-8

import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees

from variable_scripts.mgconfig import (
					RF_DIR,RF_EXT,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,
					DIST_T_DIR,DEPTH_T_DIR,TIME_T_DIR
				      )

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
print('Importing earth model from : '+MODEL_FILE_NPZ)
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)

print('P-wave time calculation')
arrivalsP = []
for i,j in enumerate(event_depth):
	print('Event '+str(i)+' of '+str(len(event_depth))+' events')
	arrivalsP.append(model_10_km.get_travel_times_geo(
                                    source_depth_in_km=j, source_latitude_in_deg=event_lat[i], 
                                    source_longitude_in_deg=event_long[i], receiver_latitude_in_deg=sta_lat[i], 
                                    receiver_longitude_in_deg=sta_long[i], phase_list='P'))

print('Ps conversion time calculation:')
phase_lst = ['P'+str(i)+'s' for i in range(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)]
print(phase_lst)



arrivals = []
for i,j in enumerate(event_depth):
	print('Event '+str(i)+' of '+str(len(event_depth))+' events')
	arrivals.append(model_10_km.get_travel_times_geo(
		                            source_depth_in_km=j, source_latitude_in_deg=event_lat[i], 
		                            source_longitude_in_deg=event_long[i], receiver_latitude_in_deg=sta_lat[i], 
		                            receiver_longitude_in_deg=sta_long[i], phase_list=phase_lst))




dist_event = [[]]*len(arrivals)
depth = [[]]*len(arrivals)
time = [[]]*len(arrivals)

for i,j in enumerate(arrivals):
            time[i] = [l.time - arrivalsP[i][0].time for k,l in enumerate(j)]
            depth[i] = [float(l.name[1:-1]) for k,l in enumerate(j)]
            dist_event[i] = [l.distance for k,l in enumerate(j)]


print('Saving .TXT files')
os.makedirs(DIST_T_DIR,exist_ok=True)
dist_event_txt = open(DIST_T_DIR+'Ps_dist_event.txt', 'w')

for i,j in enumerate(dist_event):
    dist_event_txt.write(str(list(j))+'\n')
dist_event_txt.close()


os.makedirs(DEPTH_T_DIR,exist_ok=True)
depth_txt = open(DEPTH_T_DIR+'Ps_depth.txt', 'w')

for i,j in enumerate(depth):
    depth_txt.write(str(list(j))+'\n')
depth_txt.close()

os.makedirs(TIME_T_DIR,exist_ok=True)
time_txt = open(TIME_T_DIR+'Ps_time.txt', 'w')

for i,j in enumerate(time):
    time_txt.write(str(list(j))+'\n')
time_txt.close()

print("Pds estimated!")

