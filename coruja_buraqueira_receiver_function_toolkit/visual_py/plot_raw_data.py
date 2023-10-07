'''
Script to copy raw data before stack
'''

#Useful modules

import os
from obspy import read
import shutil
import matplotlib.pyplot as plt
import pandas as pd

from parameters_py.config import (
					CUT_BEFORE_P,CUT_AFTER_P,SAMPLING_RATE,CODA_TRACE_CHECK,X_LIM_MIN,X_LIM_MAX
				   )


def plot_station_raw_RF(RF_data_df):
	RF_stack_data = RF_data_df['dataR'].mean()
	RF_stack_data_time = RF_data_df['dataR_time'].mean()
	min_y = [min(a) for a in zip(*RF_data_df['dataR'])]
	max_y = [max(a) for a in zip(*RF_data_df['dataR'])]
	plt.figure(figsize = (30,10))
	for i, row in RF_data_df.iterrows():
		plt.fill_between(RF_stack_data_time,min_y,max_y, facecolor='grey',alpha=0.5,label='Max/Min RF')
		plt.plot(RF_stack_data_time,RF_stack_data,'k',linewidth=2,label='RF stack')
		plt.text(X_LIM_MIN,max(max_y),'N = '+str(len(RF_data_df['dataR'])))
		plt.title('Receiver Functions - '+row['kstnm'])
		plt.xlim(X_LIM_MIN,X_LIM_MAX)

	plt.show()

def plot_station_raw_RF_TRACE(RF_all_data_df):
	
	fig, (ax1,ax2) = plt.subplots(2, 1, layout='constrained',sharex=True)
	# ----------------------------------
	# Axis 1 

	RF_stack_all_data = RF_all_data_df['dataR'].mean()
	RF_stack_all_data_time = RF_all_data_df['dataR_time'].mean()

	for i, row in RF_all_data_df.iterrows():
		ax1.plot(row['dataR_time'],row['dataR'],'gray',linewidth=0.5,label='RF data')
		ax1.text(-5,0.20,'N = '+str(len(RF_all_data_df)))
		ax1.set_title('Receiver Functions - '+row['kstnm'])
		ax1.set_xlim(X_LIM_MIN,X_LIM_MAX)
	ax1.plot(RF_stack_all_data_time,RF_stack_all_data,'k',linewidth=2,label='RF stack')

	# ----------------------------------
	# Axis 2

	RF_data_df = RF_all_data_df[RF_all_data_df['selection'] == True]
	RF_stack_data = RF_data_df['dataR'].mean()
	RF_stack_data_time = RF_data_df['dataR_time'].mean()

	for i, row in RF_data_df.iterrows():
		ax2.plot(row['dataR_time'],row['dataR'],'gray',linewidth=0.5,label='RF data')
		ax2.text(-5,0.20,'N = '+str(len(RF_data_df)))
		ax2.set_title('Receiver Functions - '+row['kstnm'])
		ax2.set_xlim(X_LIM_MIN,X_LIM_MAX)
	ax2.plot(RF_stack_data_time,RF_stack_data,'k',linewidth=2,label='RF stack')

	plt.show()


def plot_station_quadrant_RF_TRACE(RF_all_data_df):

	RF_data_df = RF_all_data_df[RF_all_data_df['selection'] == True]

	
	# Especifique o número desejado de grupos (quadrantes)
	num_grupos = 5

	# Use pd.cut para dividir as linhas em grupos com base na coluna 'azimute'
	RF_data_df['grupo'] = pd.cut(RF_data_df['baz'], bins=num_grupos, labels=[i for i in range(num_grupos)])

	fig, axs = plt.subplots(num_grupos,1,figsize=(5,20), layout='constrained',sharey=True,sharex=True)

	for g in range(num_grupos):
		RF_data_df_grupo = RF_data_df[RF_data_df['grupo'] == g]

		# ----------------------------------
		# Axis

		RF_stack_data = RF_data_df_grupo['dataR'].mean()
		RF_stack_data_time = RF_data_df_grupo['dataR_time'].mean()

		for i, row in RF_data_df_grupo.iterrows():
			axs[g].plot(row['dataR_time'],row['dataR'],'gray',linewidth=0.5,label='RF data')
			axs[g].text(0.8,0.95,'N = '+str(len(RF_data_df_grupo)),transform=axs[g].transAxes)
			axs[g].set_title('Grupo: '+str(g)+'('+'~baz='+str(round(RF_data_df_grupo['baz'].mean()))+'°)')
			axs[g].set_xlim(X_LIM_MIN,X_LIM_MAX)
		axs[g].plot(RF_stack_data_time,RF_stack_data,'k',linewidth=2,label='RF stack')

	plt.show()