'''
Script to copy raw data before merge
'''

#Useful modules

import os
from obspy import read
import shutil

def plot_station_raw_RF(RF_data,STA):
	plt.figure(figsize = (15,25))
	for j, i in enumerate(RF_data): 
	    plt.plot(0,j,'k',linewidth=0.7)
	    plt.title('Funcoes do Receptor '+STA)
	    #plt.xlim(-4,100)
