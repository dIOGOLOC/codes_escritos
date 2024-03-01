import matplotlib.pyplot as plt
import numpy as np
import obspy
import os
import glob
from obspy.clients.syngine import Client as SyngineClient
from obspy.taup import TauPyModel
from tqdm import tqdm
import sys
from multiprocessing import Pool

# ========== 
# INPUT DATA
# ========== 

INPUT_DATA = '/home/sysop/dados_posdoc/MTZ_2024/PRF_selected_YES_PP_FILTER_POST/'

OUTPUT_DATA = '/home/sysop/dados_posdoc/DATA_MTZ/MTZ_2024/DATA_2024_SYNTHETIC/'

# -----------
# GET FR DATA
# -----------

PRF_files = sorted(glob.glob(INPUT_DATA+'*/*_P_R.sac'))

# ---------------
# READING FR DATA
# ---------------

def read_down_save_syn(RF_file):
    #reading the data
    RF_ = obspy.read(RF_file,headonly=True)[0]

    #adjusting for a error time
    model_time = TauPyModel(model="iasp91")
    arrivals = model_time.get_travel_times(source_depth_in_km=RF_.stats.sac.evdp, distance_in_degree= RF_.stats.sac.gcarc, phase_list=["P"])
    arr = arrivals[0]
    event_time_o = arr.time-10
    event_time = RF_.stats.starttime - event_time_o

    event_DD = '{:02}'.format(event_time.day)
    event_MM = '{:02}'.format(event_time.month)
    event_YYYY = event_time.year
    event_hh = '{:02}'.format(event_time.hour)
    event_mm = '{:02}'.format(event_time.minute)
    event_sta = RF_.stats.station
    event_ray = arr.ray_param/6371
    sta_lat = RF_.stats.sac.stla
    sta_long = RF_.stats.sac.stlo

    eventid = 'GCMT:C'+str(event_YYYY)+str(event_MM)+str(event_DD)+str(event_hh)+str(event_mm)+'A'

    # --------------------------
    # DOWNLOADING SYNTHETIC DATA
    # --------------------------
    
    c_s = SyngineClient()
    model = "iasp91_2s"

    try:

        st_synth = obspy.Stream(c_s.get_waveforms(model=model,receiverlatitude=sta_lat,receiverlongitude=sta_long,networkcode='BP',stationcode=event_sta,
                                                eventid=eventid,dt="0.1",units="velocity",starttime="P-10", endtime="P+260",format='saczip'))
    
        for r,t in enumerate(st_synth):
            evdp = t.stats.sac.evdp
            gcarc = t.stats.sac.gcarc

            j = t.stats.starttime
            
            t.stats.sac.user0 = event_ray

            folder_loc_string = OUTPUT_DATA+t.stats.station+'/'+str(j.year)+'/'+str("{0:0=3d}".format(j.julday))+'/'+str(j.year)+'.'+str(j.julday)+'.'+str(j.hour)+'.'+str(j.minute)+'.'+str(j.second)+'.'+str(j.microsecond)
            os.makedirs(folder_loc_string,exist_ok=True)
            t.write(folder_loc_string+'/SYN.'+t.stats.network+'.'+t.stats.station+'.'+str(j.year)+'.'+str(j.hour)+'.'+str(j.minute)+'.'+str(j.second)+'.'+str(j.microsecond)+'.'+t.stats.channel[-1],format='SAC')

    except:
        print('Problema no Evento: '+eventid)
# -------------------------------------------------------------------------------------------
        

# ===================
# Download event data
# ===================

print('Download event data for each station via syngine')
print('\n')

with Pool(processes=20) as p:
	max_ = len(PRF_files)
	with tqdm(total=max_,desc='GCMT synthetic event download') as pbar:
		for i, _ in enumerate(p.imap_unordered(read_down_save_syn,PRF_files)):
			pbar.update()

print('Finished!')