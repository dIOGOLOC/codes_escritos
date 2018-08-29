import obspy 
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, Stream,read_inventory
import os
import glob
from obspy.io.xseed import Parser
from obspy import read_inventory
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx

from parameters_py.config import (
					DIR_DATA,SOURCE,NETWORK_CODE,NETWORK_DESCRIPTION,START_DATE,SAMPLING_RATE,LOCATION,
                    OUTPUT_XML_FILE_DIR,OUTPUT_JSON_FILE_DIR,OUTPUT_FIGURE_DIR
                    
				   )

print('Importing XML file')

inv = read_inventory(OUTPUT_XML_FILE_DIR+NETWORK_CODE+".xml")
print(inv)


def calc_PSD(data,sta_name):
    st = Stream()
    for i,j in enumerate(data):
        st += read(j)
    st.merge()
    for k,l in enumerate(st):
        l.stats.station = sta_name
        l.stats.network = NETWORK_CODE
        l.stats.location = LOCATION
     
        sta_channel = l.stats.channel
        if sta_channel == 'HH1':
            l.stats.channel = 'HHN'
            print('Calculating PPSD: '+sta_name+' station / channel: '+sta_channel)

            try:
                ppsd = PPSD(l.stats, metadata=inv)
                ppsd.add(st) 
                print(ppsd.id)
                ppsd.plot(show_coverage=False,show=True,cmap=pqlx,filename=OUTPUT_FIGURE_DIR+'/'+sta_name+'_PPSD_pqlx'+sta_channel+'.pdf')
            except:
                pass

        elif sta_channel == 'HH1j':
            l.stats.channel = 'HHN'
            print('Calculating PPSD: '+sta_name+' station / channel: '+sta_channel)
            try:
                ppsd = PPSD(l.stats, metadata=inv)
                ppsd.add(st) 
                print(ppsd.id)
                ppsd.plot(show_coverage=False,show=True,cmap=pqlx,filename=OUTPUT_FIGURE_DIR+'/'+sta_name+'_PPSD_pqlx'+sta_channel+'.pdf')
            except:
                pass
        
        elif sta_channel == 'HH2':
            l.stats.channel = 'HHE'
            print('Calculating PPSD: '+sta_name+' station / channel: '+sta_channel)
            try:
                ppsd = PPSD(l.stats, metadata=inv)
                ppsd.add(st) 
                print(ppsd.id)
                ppsd.plot(show_coverage=False,show=True,cmap=pqlx,filename=OUTPUT_FIGURE_DIR+'/'+sta_name+'_PPSD_pqlx'+sta_channel+'.pdf')
            except:
                pass
        else:
            print('Calculating PPSD: '+sta_name+' station / channel: '+sta_channel)
            try:
                ppsd = PPSD(l.stats, metadata=inv)
                ppsd.add(st) 
                print(ppsd.id)
                ppsd.plot(show_coverage=False,show=True,cmap=pqlx,filename=OUTPUT_FIGURE_DIR+'/'+sta_name+'_PPSD_pqlx'+sta_channel+'.pdf')
            except:
                pass