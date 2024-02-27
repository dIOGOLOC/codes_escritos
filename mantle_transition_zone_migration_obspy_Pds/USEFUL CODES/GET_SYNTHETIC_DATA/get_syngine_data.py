import matplotlib.pyplot as plt
import numpy as np
import obspy
import os
import glob
from obspy.clients.syngine import Client as SyngineClient
from obspy.taup import TauPyModel
from tqdm import tqdm
import sys

# ========== 
# INPUT DATA
# ========== 

INPUT_DATA = '/home/sysop/dados_posdoc/MTZ_2024/PRF_SEISPY_DATA_NO_FILTER_PP/'

OUTPUT_DATA = '/home/sysop/dados_posdoc/MTZ_2024/SYNTHETIC_DATA_MTZ/'

# -----------
# GET FR DATA
# -----------

PRF_folders = sorted(glob.glob(INPUT_DATA+'*/*_P_R.sac'))

# ---------------
# READING FR DATA
# ---------------

ev = obspy.Stream()
for i,j in enumerate(PRF_folders):
    ev += obspy.read(j,headonly=True)

# ------------------
# ALLOCATING FR DATA
# ------------------
    
eventid_lst = []
lat_lst = []
long_lst = []
sta_lst = []

for i,j in tqdm(enumerate(ev),total=len(ev),desc='Getting the GCMT event name'):
    event_time = obspy.UTCDateTime(year=j.stats.sac.nzyear, julday=j.stats.sac.nzjday)
    event_DD = '{:02}'.format(event_time.day)
    event_MM = '{:02}'.format(event_time.month)
    event_YYYY = j.stats.sac.nzyear
    event_hh = '{:02}'.format(j.stats.sac.nzhour)
    event_mm = '{:02}'.format(j.stats.sac.nzmin)
    event_julday = j.stats.sac.nzjday
    event_depth = j.stats.sac.evdp
    event_lat = j.stats.sac.evla
    event_long = j.stats.sac.evlo
    event_dist = j.stats.sac.dist
    event_gcarc = j.stats.sac.gcarc
    event_sta = j.stats.station
    sta_lat = j.stats.sac.stla
    sta_long = j.stats.sac.stlo

    lat_lst.append(sta_lat)
    long_lst.append(sta_long)
    sta_lst.append(event_sta)

    eventid_lst.append('GCMT:C'+str(event_YYYY)+str(event_MM)+str(event_DD)+str(event_hh)+str(event_mm)+'A')

# --------------------------
# DOWNLOADING SYNTHETIC DATA
# --------------------------
    
print(eventid_lst[-10:-1])

c_s = SyngineClient()
model = "iasp91_2s"

for k,l in tqdm(enumerate(eventid_lst),total=len(eventid_lst),desc='Downloading Synthetic data'):
    try:
        st_synth = obspy.Stream(c_s.get_waveforms(model=model,receiverlatitude=lat_lst[k],receiverlongitude=long_lst[k],networkcode='BP',stationcode=sta_lst[k],
                                                eventid=eventid_lst[k],dt="0.1",units="velocity",starttime="P-10", endtime="P+260",format='saczip'))
    


        for r,t in enumerate(st_synth):
            evdp = t.stats.sac.evdp
            gcarc = t.stats.sac.gcarc

            model_time = TauPyModel(model="iasp91")
            arrivals = model_time.get_travel_times(source_depth_in_km=evdp, distance_in_degree=gcarc, phase_list=["P"])
            arr = arrivals[0]
            j = t.stats.starttime
            t.stats.sac.user0 = arr.ray_param/6371
            folder_loc_string = OUTPUT_DATA+t.stats.station+'/'+str(j.year)+'/'+str("{0:0=3d}".format(j.julday))+'/'+str(j.year)+'.'+str(j.julday)+'.'+str(j.hour)+'.'+str(j.minute)+'.'+str(j.second)+'.'+str(j.microsecond)
            os.makedirs(folder_loc_string,exist_ok=True)
            t.write(folder_loc_string+'/SYN.'+t.stats.network+'.'+t.stats.station+'.'+str(j.year)+'.'+str(j.hour)+'.'+str(j.minute)+'.'+str(j.second)+'.'+str(j.microsecond)+'.'+t.stats.channel[-1],format='SAC')

    except:
        print('Problema no Evento: '+eventid_lst[k])