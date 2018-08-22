'''
Script to copy raw data before merge
'''

#Useful modules

import os
from obspy import read
import shutil
import matplotlib.pyplot as plt

from parameters_py.config import (
					CUT_BEFORE_P,CUT_AFTER_P,SAMPLING_RATE
				   )


def plot_station_raw_RF(RF_data,STA):
	plt.figure(figsize = (15,25))
	for i, j in enumerate(RF_data): 
		plt.plot(j,'k',linewidth=0.5,alpha=0.5)
		plt.text(-0.1,0,str(len(RF_data)))
		plt.title('Funcoes do Receptor '+STA)
		plt.xlim(0,700)
	plt.show()
