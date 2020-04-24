'''
Script to plot PPSD data based in 
https://docs.obspy.org/tutorial/code_snippets/probabilistic_power_spectral_density.html
'''

from obspy.imaging.cm import pqlx
import matplotlib.pyplot as plt
import pandas as pd
import obspy
import os
import glob
import json
import numpy as np
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib as mpl
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from obspy.signal import PPSD
from obspy.clients.arclink.client import Client
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,INITIAL_DATE,FINAL_DATE,TIME_OF_WEEKDAY_DAY,TIME_OF_WEEKDAY_START_HOUR,
					TIME_OF_WEEKDAY_FINAL_HOUR,PERIOD_PSD,AMP_PSD_MIN,AMP_PSD_MAX
				   )

# ==================================
# Function to plot TOTAL PPSD DATA
# ==================================

def plot_PPSD_TOTAL_data(date_lst):
    os.chdir(date_lst)
    files = sorted(glob.glob('*.npz'))
    ppsd = PPSD.load_npz(files[0],allow_pickle=True)

    [ppsd.add_npz(i,allow_pickle=True) for i in files[1:]]
    os.makedirs(OUTPUT_FIGURE_DIR+'TOTAL/'+ppsd.station+'/',exist_ok=True)
    ppsd.plot(cmap=pqlx,filename=OUTPUT_FIGURE_DIR+'TOTAL/'+ppsd.station+'/'+ppsd.station+'.'+ppsd.channel+'.pdf')

def plot_PPSD_WINDOWED_data(date_lst):
    os.chdir(date_lst)
    files = sorted(glob.glob('*.npz'))
    ppsd = PPSD.load_npz(files[0],allow_pickle=True)

    [ppsd.add_npz(i,allow_pickle=True) for i in files[1:]]
    ppsd.calculate_histogram(starttime=obspy.UTCDateTime(INITIAL_DATE),endtime=obspy.UTCDateTime(FINAL_DATE),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])    
    folder_output = OUTPUT_FIGURE_DIR+'WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd.station+'/'
    os.makedirs(folder_output,exist_ok=True)
    ppsd.plot(cmap=pqlx,filename=folder_output+ppsd.station+'.'+ppsd.channel+'.pdf')


def plot_PPSD_by_period(directory_data):
	
	data_lista = []

	print('Looking for data in the directory = '+directory_data)

	for root, dirs, files in os.walk(directory_data):
		for name in files:
			data_lista.append(os.path.join(root, name))

	data_lista = sorted(data_lista)
    
	dataframe_lista = []
	#create a empty dataframe with pandas
	for i,j in enumerate(data_lista):
		print("Extracting data from PPSD header "+str(i+1)+" of "+str(len(data_lista)))
                
		#Reading header from data
		ppsd = PPSD.load_npz(j,allow_pickle=True)
                
		#----------------------------
		#Dataframe starting

		network = ppsd.network
		station = ppsd.station
		channel = ppsd.channel

		flat_time_lst = ppsd.times_processed

		DATETIME = str(ppsd.times_data[0][0].year)+','+str(ppsd.times_data[0][0].month)+','+str(ppsd.times_data[0][0].day)

		#Contador da lista de horas
		time_flat_time_lst = [[]]*24
		for g,h in enumerate(np.arange(24)):
			lst_time = []
			for x,c in enumerate(flat_time_lst):		
				if c.hour == h:   		
					lst_time.append(ppsd.extract_psd_values(PERIOD_PSD)[0][x])
			time_flat_time_lst[g] = lst_time

		AMPLITUDE_HOUR = [[]]*24		
		for q,w in enumerate(time_flat_time_lst):
			AMPLITUDE_HOUR[q] = np.mean(w)
          

		dataframe_lista.append(pd.DataFrame([[network],[station],[channel],[DATETIME],[AMPLITUDE_HOUR]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','AMPLITUDE_HOUR']).T)
		print(pd.DataFrame([[network],[station],[channel],[DATETIME],[AMPLITUDE_HOUR]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','AMPLITUDE_HOUR']).T)
		print()
		print()

		#Dataframe ending
		#----------------------------


	df = pd.concat(dataframe_lista, ignore_index=True)

	#Sorting according to station
	station_lista = list(set(df['STATION']))

	for i,j in enumerate(station_lista):
		df_sta = df[df['STATION'] == j]

		channel_lista = list(set(df_sta['CHANNEL']))
		channel_lista = sorted(channel_lista)
		
		# ==========================================================
		# Calculating datetime between INITIAL_DATE and  FINAL_DATE		
		# ==========================================================

		datatime_initial = datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day)

		datatime_final = datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day)

		datetime_lista = np.arange(datatime_initial, datatime_final, datetime.timedelta(days=1)).astype(datetime.datetime)


		xlim_initial = mdates.date2num(datatime_initial)
		xlim_final = mdates.date2num(datatime_final)
            
		#----------------------------
		#Function to check if the dates in data set are inside the period chosen (INITIAL_DATE to FINAL_DATE)

		def check_datetime_in_period(datetime_lst,df_DATETIME,df_AMPLITUDE_HOUR):

			array_to_plot_by_xlim = []
			for x,c in enumerate(datetime_lst):
				lista_temp = []
				for t,y in enumerate(df_DATETIME):
					if datetime.datetime(obspy.UTCDateTime(y).year,obspy.UTCDateTime(y).month,obspy.UTCDateTime(y).day) == c:
						lista_temp.append(df_AMPLITUDE_HOUR[df_DATETIME[df_DATETIME == y].index[0]])
				array_to_plot_by_xlim.append(lista_temp)

			data_x_axis = []
			for x,c in enumerate(array_to_plot_by_xlim):
				if c != []:
					data_x_axis.append(c[0][::-1])
				else:
					data_x_axis.append(np.full_like(np.arange(24),np.nan,dtype=np.double))

			data_x_axis = np.array(data_x_axis).T

			return data_x_axis

        # ====================================
        # Function to plot DATA availability
        # ====================================
		
		#x axis parameters

		days1 = DayLocator(interval=1)   # every day
		days5 = DayLocator(interval=int(len(datetime_lista)*5/100))   # every day
		months = MonthLocator()  # every month
		yearsFmt = DateFormatter('%Y-%m-%d')
        
		days1.MAXTICKS = 10000


		#Matplotlib parameters
		fig, ax = plt.subplots(nrows=len(channel_lista), ncols=1,sharex=True,sharey=True,figsize=(40,15))
		fig.suptitle(j,fontsize=25,y=0.9)
		for k,l in enumerate(channel_lista):

			df_ch = df_sta[df_sta['CHANNEL'] == l]

			data_x_axis = check_datetime_in_period(datetime_lista,df_ch['DATETIME'],df_ch['AMPLITUDE_HOUR'])

			im = ax[k].imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.viridis,interpolation=None,vmin=AMP_PSD_MIN,vmax=AMP_PSD_MAX)
			ax[k].set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
			ax[k].yaxis.set_major_locator(MultipleLocator(4))
			ax[k].yaxis.set_minor_locator(MultipleLocator(1))
			ax[k].xaxis.set_major_locator(days5)
			ax[k].xaxis.set_major_formatter(yearsFmt)
			ax[k].xaxis.set_minor_locator(days1)
			ax[k].tick_params(which='minor', length=4)
			ax[k].tick_params(which='major', length=10)
			ax[k].set_ylim(0,24)
			ax[k].set_yticklabels([' ',' ', '1h', '5h', '9h', '13h', '17h', '21h',' ']) #y axis according to Brazil UTC-3
			ax[k].set_ylabel(l,fontsize=15)
			ax[k].grid(b=True, which='major', color='k', linestyle='-')
			ax[k].grid(b=True, which='minor', color='k', linestyle='-')

        
		plt.setp(ax[k].xaxis.get_majorticklabels(), fontsize=10, rotation=30)
		ax[-1].set_xlabel('Time', fontsize=20)
        
        #criando a localização da barra de cores:
		axins = inset_axes(ax[0],
							width="10%",  # width = 10% of parent_bbox width
							height="5%",  # height : 50%
							loc='upper left',
							bbox_to_anchor=(0.85, 0.1, 1, 1),
							bbox_transform=ax[0].transAxes,
							borderpad=0,
                           )
		cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[AMP_PSD_MIN,np.mean([AMP_PSD_MIN,AMP_PSD_MAX]),AMP_PSD_MAX],label='Amplitude '+r'$[m^2/s^4/Hz][dB]$')
        
		os.makedirs(OUTPUT_FIGURE_DIR,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+j+'_'+'PSD_BY_PERIOD_'+str(PERIOD_PSD)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).year)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).month)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).day)+'_'+str(obspy.UTCDateTime(FINAL_DATE).year)+'_'+str(obspy.UTCDateTime(FINAL_DATE).month)+'_'+str(obspy.UTCDateTime(FINAL_DATE).day)+'.pdf',dpi=500)