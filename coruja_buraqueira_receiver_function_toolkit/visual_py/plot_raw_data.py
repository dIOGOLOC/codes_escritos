'''
Script to copy raw data before stack
'''

#Useful modules

import os
from obspy import read
import shutil
import matplotlib.pyplot as plt

from parameters_py.config import (
					CUT_BEFORE_P,CUT_AFTER_P,SAMPLING_RATE,CODA_TRACE_CHECK,X_LIM_MIN,X_LIM_MAX
				   )


def plot_station_raw_RF(RF_data,RF_data_time,STA):
	RF_stack_data = [sum(i)/len(RF_data) for i in zip(*RF_data)]
	min_y = [min(a) for a in zip(*RF_data)]
	max_y = [max(a) for a in zip(*RF_data)]
	plt.figure(figsize = (30,10))
	for i, j in enumerate(RF_data): 
		plt.fill_between(RF_data_time,min_y,max_y, facecolor='grey',alpha=0.5,label='Max/Min RF')
		plt.plot(RF_data_time,RF_stack_data,'k',linewidth=2,label='RF stack')
		plt.text(X_LIM_MIN,max(max_y),'N = '+str(len(RF_data)))
		plt.title('Receiver Functions - '+STA)
		plt.xlim(X_LIM_MIN,X_LIM_MAX)

	plt.show()

def plot_station_raw_RF_TRACE(RF_data,RF_data_time,STA):
	RF_stack_data = [sum(i)/len(RF_data) for i in zip(*RF_data)]
	plt.figure(figsize = (30,10))
	for i, j in enumerate(RF_data): 
		plt.plot(RF_data_time,j,'gray',linewidth=0.5,label='RF data')
		plt.plot(RF_data_time,RF_stack_data,'k',linewidth=2,label='RF stack')
		plt.text(-5,0.20,'N = '+str(len(RF_data)))
		plt.title('Receiver Functions - '+STA)
		plt.xlim(X_LIM_MIN,X_LIM_MAX)

	plt.show()

