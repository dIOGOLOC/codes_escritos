#!/usr/bin/python -u
'''
--------------------------------------------------------------------------------
       Function to plot local events and their Wadati diagram
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 02/2022


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
This code will plot the local dataset according to a given list of an events.


Inputs:
Event traces (format: SAC)


Outputs:
Figures (PNG)


Examples of Usage (in command line):
   >> python plot_LOCAL_EVENT_WADATI.py

'''

import time
import numpy as np
import os
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool
from glob import glob

from obspy import read,read_inventory, UTCDateTime, Stream
from obspy.geodetics import gps2dist_azimuth

from matplotlib.dates import YearLocator, MonthLocator, DayLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy

import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter,LatitudeLocator,LongitudeLocator

from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression,HuberRegressor,TheilSenRegressor
from sklearn.metrics import mean_squared_error

# ------------------------------------------------------------------------------
from visual_py.event_plot import plot_event_data,plot_map_event_data,plot_map_event_data_hydrophone

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,OUTPUT_EV_DIR,DIR_STATUS,NUM_PROCESS,OUTPUT_PSD_DIR,XML_FILE,
					BOUNDARY_STATES_SHP
					)
# ------------------------------------------------------------------------------

# =======================
# Retrieving MARKES files
# =======================

print('\n')
print('Retrieving MARKERS for each event')
print('\n')


markes_file = sorted(glob(OUTPUT_EV_DIR+'Local/NETWORK_MSEED_FILES/**/*_markers'))

event_dic_lst = []
for mid,mfile in enumerate(markes_file):

	# splitting subdir/basename
	subdir, filename = os.path.split(mfile)
	event_name = subdir.split('/')[-1]
	#--------------------------------------

	dic_lst = []

	picks_markers = np.genfromtxt(mfile,dtype='str')

	for pick in picks_markers:

		day = pick[1]
		time = pick[2]
		id = pick[4]
		phase = pick[8]
		datetime_UTC = UTCDateTime(day+'T'+time)
		network = pick[4].split('.')[0]
		station = pick[4].split('.')[1]
		channel = pick[4].split('..')[1]

		inv = read_inventory(XML_FILE+'.'.join([network,station,'xml']))
		coordinates_lst = inv[0][0]

		dic = {
				'datetime':datetime_UTC.datetime,
				'id':id,
				'network':network,
				'station':station,
				'channel':channel,
				'phase':phase,
				'latitude': coordinates_lst.latitude,
				'longitude': coordinates_lst.longitude,
				'event_name': event_name,
				}
		dic_lst.append(dic)

	event_dic_lst.append(dic_lst)

# ------------------------------------------------------------------------------
event_dic_time_lst = []
for d in tqdm(event_dic_lst,desc='Dic loop'):
	df = pd.DataFrame(d)
	station_lst = sorted(list(set(df['station'].values)))

	dic_sta_lst = []
	for sta in station_lst:

		df_sta = df[df['station'] == sta]
		df_P = df_sta[(df_sta['phase'] == 'P') | (df_sta['phase'] == 'P')]
		df_S = df_sta[(df_sta['phase'] == 'S') | (df_sta['phase'] == 'S')]

		time_P = (df_P['datetime'] - df_P['datetime'].min()).mean() + df_P['datetime'].min()
		time_S = (df_S['datetime'] - df_S['datetime'].min()).mean() + df_S['datetime'].min()

		S_P_time = (time_S - time_P).total_seconds()

		datetime_event_year = df_P['datetime'].dt.strftime('%Y').values[0]
		datetime_event_julday = df_P['datetime'].dt.strftime('%j').values[0]

		# ------------------------------------------
		# Retrieving mseed files
		sta_ev_name = df_sta['event_name'].values[0]
		station_file_HHZ = read(OUTPUT_EV_DIR+'Local/'+df_sta['network'].values[0]+'/'+sta+'/'+datetime_event_year+'/'+datetime_event_julday+'/'+sta_ev_name+'/*'+sta_ev_name+'.Z')
		station_file_HHZ.detrend("linear")
		station_file_HHZ.detrend("demean")
		station_file_HHZ.taper(max_percentage=0.05, type="hann")
		pre_filt = [0.001, 0.005, 45., 50.]
		inv = read_inventory(XML_FILE+'.'.join([df_sta['network'].values[0],sta,'xml']))
		station_file_HHZ.remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)
		station_file_HHZ.filter('bandpass',freqmin=4, freqmax=16)

		station_file_HHN = read(OUTPUT_EV_DIR+'Local/'+df_sta['network'].values[0]+'/'+sta+'/'+datetime_event_year+'/'+datetime_event_julday+'/'+sta_ev_name+'/*'+sta_ev_name+'.N')
		station_file_HHN.detrend("linear")
		station_file_HHN.detrend("demean")
		station_file_HHN.taper(max_percentage=0.05, type="hann")
		pre_filt = [0.001, 0.005, 45., 50.]
		inv = read_inventory(XML_FILE+'.'.join([df_sta['network'].values[0],sta,'xml']))
		station_file_HHN.remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)
		station_file_HHN.filter('bandpass',freqmin=4, freqmax=16)

		station_file_HHE = read(OUTPUT_EV_DIR+'Local/'+df_sta['network'].values[0]+'/'+sta+'/'+datetime_event_year+'/'+datetime_event_julday+'/'+sta_ev_name+'/*'+sta_ev_name+'.E')
		station_file_HHE.detrend("linear")
		station_file_HHE.detrend("demean")
		station_file_HHE.taper(max_percentage=0.05, type="hann")
		pre_filt = [0.001, 0.005, 45., 50.]
		inv = read_inventory(XML_FILE+'.'.join([df_sta['network'].values[0],sta,'xml']))
		station_file_HHE.remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)
		station_file_HHE.filter('bandpass',freqmin=4, freqmax=16)

		# ------------------------------------------
		evla = station_file_HHE[0].stats.sac.evla
		evlo = station_file_HHE[0].stats.sac.evlo
		# ------------------------------------------

		sta_latitude = df_sta['latitude'].values[0]
		sta_longitude = df_sta['longitude'].values[0]

		dic_sta = {
					'station':sta,
					'latitude':sta_latitude,
					'longitude':sta_longitude,
					'time_P':UTCDateTime(time_P),
					'time_S':UTCDateTime(time_S),
					'S-P':S_P_time,
					'HHZ':station_file_HHZ,
					'HHN':station_file_HHN,
					'HHE':station_file_HHE,
					'evla':evla,
					'evlo':evlo,
					'event_name':sta_ev_name
					}

		dic_sta_lst.append(dic_sta)
	event_dic_time_lst.append(pd.DataFrame(dic_sta_lst))

# ------------------------------------------------------------------------------

for times_dictionary in event_dic_time_lst:

	evla = times_dictionary['evla'].to_numpy()[0]
	evlo = times_dictionary['evlo'].to_numpy()[0]
	sta_ev_name = times_dictionary['event_name'].to_numpy()[0]
	# --------------------------------------------------------------------------
	# Simple Linear Regression With scikit-learn

	#Provide data:
	X_datetime = np.array([t_value.datetime for t_value in times_dictionary['time_P'].to_numpy()])
	X_min = X_datetime.min()

	X_date = np.array([(t_value.datetime-X_min).total_seconds() for t_value in times_dictionary['time_P'].to_numpy()])
	X = np.array([(t_value.datetime-X_min).total_seconds() for t_value in times_dictionary['time_P'].to_numpy()]).reshape(-1, 1)
	Y = times_dictionary['S-P'].to_numpy().reshape(-1, 1)

	#Create a model and fit it:
	model = TheilSenRegressor(n_jobs=12,max_iter=500)

	#It’s time to start using the model. First, you need to call .fit() on model:
	model.fit(X, Y)

	#Get results:
	#coefficient of determination
	r_sq = round(model.score(X, Y),3)
	#intercept:
	intercept = model.intercept_
	#slope:
	slope = model.coef_[0]
	Vp_Vs = round(slope+1,2)

	#Predict response (y_pred = model.intercept_ + model.coef_ * x):
	x_0 = intercept/slope

	print('--------------------------------------------')
	print('Event time origin:',UTCDateTime(X_min)+x_0)
	print('--------------------------------------------')

	y_pred = model.predict(X)

	# --------------------------------------------------------------------------
	# for far field the hypocentre distance can be estimated as:
	#D = 8 (Ts − Tp)
	# for near field the hypocentre distance can be estimated as:
	#D = 10 (Ts − Tp)
	hypocentre_dist = []
	for sp_id,sp_time in enumerate(times_dictionary['S-P'].to_numpy()):
		if 'OBS' in times_dictionary['station']:
			hypocentre_dist.append(10*sp_time)
		else:
			hypocentre_dist.append(8*sp_time)

	# -----------------------------------------------------------------------
	# MAP
	# -----------------------------------------------------------------------

	fig = plt.figure(figsize=(7, 7))
	gs = gridspec.GridSpec(nrows=1, ncols=1)

	#-------------------------------------------
	crs = ccrs.PlateCarree(central_longitude=-40)
	map_loc = fig.add_subplot(gs[0],projection=crs)
	map_loc.set_title('Data: '+UTCDateTime(X_min).strftime('%d/%m/%Y - %H:%M:%S'),fontsize=20)

	LLCRNRLON_LARGE = -50
	URCRNRLON_LARGE = -35
	LLCRNRLAT_LARGE = -30
	URCRNRLAT_LARGE = -15

	map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE],crs=ccrs.PlateCarree())
	# Create a Stamen Terrain instance.
	stamen_terrain = cimgt.Stamen('terrain-background')

	# Add the Stamen data at zoom level 8.
	map_loc.add_image(stamen_terrain, 5)

	# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
	states_provinces = cfeature.NaturalEarthFeature(
													category='cultural',
													name='admin_1_states_provinces_lines',
													scale='50m',
													facecolor='none'
													)

	map_loc.add_feature(cfeature.LAND,)
	map_loc.add_feature(cfeature.COASTLINE)
	map_loc.add_feature(states_provinces, edgecolor='k',linewidth=0.5)

	gl = map_loc.gridlines(color='gray',linewidth=0.5,linestyle='--',draw_labels=True)

	gl.xlocator = LongitudeLocator(4)
	gl.ylocator = LatitudeLocator(4)
	gl.xformatter = LongitudeFormatter()
	gl.yformatter = LatitudeFormatter()

	gl.xlabel_style = {'size': 15, 'color': 'gray'}
	gl.xlabel_style = {'color': 'black', 'weight': 'bold'}
	gl.ylabel_style = {'size': 15, 'color': 'gray'}
	gl.ylabel_style = {'color': 'black', 'weight': 'bold'}

	map_loc.yaxis.set_ticks_position('both')
	map_loc.xaxis.set_ticks_position('both')
	map_loc.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True)
	# Use the cartopy interface to create a matplotlib transform object
	# for the Geodetic coordinate system. We will use this along with
	# matplotlib's offset_copy function to define a coordinate system which
	# translates the text by 25 pixels to the left.
	geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
	text_transform = offset_copy(geodetic_transform, units='dots', y=17,x=33)
	text_transform_mag = offset_copy(geodetic_transform, units='dots', y=-25,x=20)

	#Ploting circles to find the event location
	for i,j in enumerate(times_dictionary['longitude'].to_numpy()):
		map_loc.tissot(rad_km=hypocentre_dist[i], lons=times_dictionary['longitude'].to_numpy()[i],lats= times_dictionary['latitude'].to_numpy()[i], n_samples=50,color='white',edgecolors='k',alpha=0.25,zorder=31)

	#Plotting station
	map_loc.scatter(times_dictionary['longitude'].to_numpy(), times_dictionary['latitude'].to_numpy(), marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
	map_loc.scatter(evlo, evla, marker='*',s=200,c='r',edgecolors='k', transform=ccrs.PlateCarree())

	#Saving figure
	folder_name = OUTPUT_FIGURE_DIR+'EVENTS/Local/Map/'
	os.makedirs(folder_name,exist_ok=True)
	fig.savefig(folder_name+'MAP_Event_circle_'+sta_ev_name+'.png')

	# -------------------------------------------
	# WADATI DIAGRAM
	# -------------------------------------------
	fig = plt.figure(figsize=(10, 10))
	gs = gridspec.GridSpec(nrows=1, ncols=1)
	ax1 = fig.add_subplot(gs[0])

	# --------------------------------------------------------------------------
	for i,j in enumerate(times_dictionary['time_P'].to_numpy()):
		ax1.scatter(X_date[i],times_dictionary['S-P'].to_numpy()[i],marker='+',s=50,c='k')

	ax1.set_ylabel('Tempo S-P (s)',fontsize=14)
	ax1.set_xlabel('Tempo de chegada da onda P (s) após '+X_min.strftime('%H:%M:%S'),fontsize=14)
	ax1.set_title('Diagrama de Wadati - data: '+times_dictionary['time_P'].to_numpy()[0].strftime('%d/%m/%Y'),fontsize=14)

	ax1.plot(X,y_pred,c='gray',ls='--',lw=1,alpha=0.7,zorder=-1)
	ax1.text(0.1, 0.9,'Vp/Vs:'+str(Vp_Vs)+'\n'+'\n'+'R²:'+str(r_sq),horizontalalignment='center',verticalalignment='center',weight='bold',bbox=dict(facecolor='white',edgecolor='k', alpha=0.5, boxstyle='round'),transform=ax1.transAxes)

	# format the ticks
	ax1.xaxis.set_major_locator(MultipleLocator(20))
	ax1.xaxis.set_minor_locator(MultipleLocator(1))

	folder_name = OUTPUT_FIGURE_DIR+'EVENTS/Local/Wadati/'
	os.makedirs(folder_name,exist_ok=True)
	fig.savefig(folder_name+'Wadati_Event_circle_'+sta_ev_name+'.png')

	#-------------------------------------------
	# SESIMOGRAMS and MARKERS
	# -------------------------------------------

	fig = plt.figure(figsize=(20, 10))
	gs = gridspec.GridSpec(nrows=1, ncols=3)

	axs0 = fig.add_subplot(gs[0])

	minutes = mdates.MinuteLocator(1)   # every year
	seconds = mdates.SecondLocator(1)  # every second
	years_fmt = mdates.DateFormatter('%H:%M:%S')

	inset_size = 1/len(times_dictionary['HHZ'].to_numpy())
	axis_size = np.linspace(0, 1-inset_size, num=len(times_dictionary['HHZ'].to_numpy()))
	org_lst = np.argsort(times_dictionary['S-P'].to_numpy())

	# this is an inset axes over the main axes
	for idx,i in enumerate(org_lst):

		inset_ax = axs0.inset_axes([0,axis_size[idx], 1, inset_size])
		data_to_plot = times_dictionary['HHZ'].to_numpy()[i][0].trim(starttime=times_dictionary['time_P'].to_numpy()[i]-5, endtime=times_dictionary['time_S'].to_numpy()[i]+30)
		data_to_plot.taper(max_percentage=0.1, type="hann")
		data_to_plot.filter('bandpass',freqmin=4, freqmax=16)
		data_y = data_to_plot.data
		data_x_utc = data_to_plot.times('utcdatetime')-times_dictionary['time_P'].to_numpy()[i]
		data_x = data_x_utc


		#plotting data
		inset_ax.plot(data_x,data_y,c='k',ls='-',lw=1)

		inset_ax.axvline(x=(times_dictionary['time_S'].to_numpy()[i]-times_dictionary['time_P'].to_numpy()[i]),ymin=0, ymax=1,c='r',lw=3,ls='--')
		inset_ax.axvline(x=0,ymin=0, ymax=1,c='k',lw=1,ls='--')

		inset_ax.set_xlim(-3,80)
		inset_ax.set_ylabel(times_dictionary['station'].to_numpy()[i],fontsize=15)
		inset_ax.set_yticks([])
		inset_ax.set_xticks([])

		#------------------------------------------------------------
		epi_dist, az, baz  = gps2dist_azimuth(times_dictionary['evla'].to_numpy()[i], times_dictionary['evlo'].to_numpy()[i], times_dictionary['latitude'].to_numpy()[i], times_dictionary['longitude'].to_numpy()[i], a=6378137.0, f=0.0033528106647474805)
		epi_dist = epi_dist / 1000
		#------------------------------------------------------------

		inset_ax.text(0.9, 0.85, str(round(epi_dist))+' km',horizontalalignment='center',fontsize=15,verticalalignment='center',bbox=dict(facecolor='none', alpha=0.6,edgecolor='w'),transform=inset_ax.transAxes)

	axs0.axvline(x=0,ymin=0, ymax=1,c='k',lw=1,ls='--')
	axs0.set_xlim(-3,80)
	axs0.set_yticks([])
	axs0.set_title('HHZ (filter: 4-16 Hz)',fontsize=20)
	axs0.set_xlabel('Tempo após a onda P (s)',fontsize=20)

    # format the ticks
	axs0.xaxis.set_major_locator(MultipleLocator(20))
	axs0.xaxis.set_minor_locator(MultipleLocator(5))
	axs0.xaxis.set_tick_params(labelsize=15)

	# -------------------------------------------

	axs1 = fig.add_subplot(gs[1])

	minutes = mdates.MinuteLocator(1)   # every year
	seconds = mdates.SecondLocator(1)  # every second
	years_fmt = mdates.DateFormatter('%H:%M:%S')

	inset_size = 1/len(times_dictionary['HHN'].to_numpy())
	axis_size = np.linspace(0, 1-inset_size, num=len(times_dictionary['HHN'].to_numpy()))
	org_lst = np.argsort(times_dictionary['S-P'].to_numpy())

	# this is an inset axes over the main axes
	for idx,i in enumerate(org_lst):

		inset_ax = axs1.inset_axes([0,axis_size[idx], 1, inset_size])
		data_to_plot = times_dictionary['HHN'].to_numpy()[i][0].trim(starttime=times_dictionary['time_P'].to_numpy()[i]-5, endtime=times_dictionary['time_S'].to_numpy()[i]+30)
		data_to_plot.taper(max_percentage=0.1, type="hann")
		data_to_plot.filter('bandpass',freqmin=4, freqmax=16)
		data_y = data_to_plot.data
		data_x_utc = data_to_plot.times('utcdatetime')-times_dictionary['time_P'].to_numpy()[i]
		data_x = data_x_utc

		inset_ax.plot(data_x,data_y,c='k',ls='-',lw=1)

		#plotting data
		inset_ax.plot(data_x,data_y,c='k',ls='-',lw=1)
		inset_ax.axvline(x=(times_dictionary['time_S'].to_numpy()[i]-times_dictionary['time_P'].to_numpy()[i]),ymin=0, ymax=1,c='r',lw=3,ls='--')
		inset_ax.axvline(x=0,ymin=0, ymax=1,c='k',lw=1,ls='--')

		inset_ax.set_xlim(-3,80)
		inset_ax.set_ylabel(times_dictionary['station'].to_numpy()[i],fontsize=15)
		inset_ax.set_yticks([])
		inset_ax.set_xticks([])

		#------------------------------------------------------------
		epi_dist, az, baz  = gps2dist_azimuth(times_dictionary['evla'].to_numpy()[i], times_dictionary['evlo'].to_numpy()[i], times_dictionary['latitude'].to_numpy()[i], times_dictionary['longitude'].to_numpy()[i], a=6378137.0, f=0.0033528106647474805)
		epi_dist = epi_dist / 1000
		#------------------------------------------------------------

		inset_ax.text(0.9, 0.85, str(round(epi_dist))+' km',horizontalalignment='center',fontsize=15,verticalalignment='center',bbox=dict(facecolor='none', alpha=0.6,edgecolor='w'),transform=inset_ax.transAxes)


	axs1.axvline(x=0,ymin=0, ymax=1,c='k',lw=1,ls='--')
	axs1.set_xlim(-3,80)
	axs1.set_yticks([])
	axs1.set_title('HHN (filter: 4-16 Hz)',fontsize=20)
	axs1.set_xlabel('Tempo após a onda P (s)',fontsize=20)

    # format the ticks
	axs1.xaxis.set_major_locator(MultipleLocator(20))
	axs1.xaxis.set_minor_locator(MultipleLocator(5))
	axs1.xaxis.set_tick_params(labelsize=15)

	# -------------------------------------------

	axs2 = fig.add_subplot(gs[2])

	minutes = mdates.MinuteLocator(1)   # every year
	seconds = mdates.SecondLocator(1)  # every second
	years_fmt = mdates.DateFormatter('%H:%M:%S')

	inset_size = 1/len(times_dictionary['HHE'].to_numpy())
	axis_size = np.linspace(0, 1-inset_size, num=len(times_dictionary['HHE'].to_numpy()))
	org_lst = np.argsort(times_dictionary['S-P'].to_numpy())

	# this is an inset axes over the main axes
	for idx,i in enumerate(org_lst):

		inset_ax = axs2.inset_axes([0,axis_size[idx], 1, inset_size])
		data_to_plot = times_dictionary['HHE'].to_numpy()[i][0].trim(starttime=times_dictionary['time_P'].to_numpy()[i]-5, endtime=times_dictionary['time_S'].to_numpy()[i]+30)
		data_to_plot.taper(max_percentage=0.1, type="hann")
		data_to_plot.filter('bandpass',freqmin=4, freqmax=16)
		data_y = data_to_plot.data
		data_x_utc = data_to_plot.times('utcdatetime')-times_dictionary['time_P'].to_numpy()[i]
		data_x = data_x_utc

		inset_ax.plot(data_x,data_y,c='k',ls='-',lw=1)

		#plotting data
		inset_ax.plot(data_x,data_y,c='k',ls='-',lw=1)
		inset_ax.axvline(x=(times_dictionary['time_S'].to_numpy()[i]-times_dictionary['time_P'].to_numpy()[i]),ymin=0, ymax=1,c='r',lw=3,ls='--')
		inset_ax.axvline(x=0,ymin=0, ymax=1,c='k',lw=1,ls='--')

		inset_ax.set_xlim(-3,80)
		inset_ax.set_ylabel(times_dictionary['station'].to_numpy()[i],fontsize=15)
		inset_ax.set_yticks([])
		inset_ax.set_xticks([])

		#------------------------------------------------------------
		epi_dist, az, baz  = gps2dist_azimuth(times_dictionary['evla'].to_numpy()[i], times_dictionary['evlo'].to_numpy()[i], times_dictionary['latitude'].to_numpy()[i], times_dictionary['longitude'].to_numpy()[i], a=6378137.0, f=0.0033528106647474805)
		epi_dist = epi_dist / 1000
		#------------------------------------------------------------

		inset_ax.text(0.9, 0.85, str(round(epi_dist))+' km',horizontalalignment='center',fontsize=15,verticalalignment='center',bbox=dict(facecolor='none', alpha=0.6,edgecolor='w'),transform=inset_ax.transAxes)

	axs2.axvline(x=0,ymin=0, ymax=1,c='k',lw=1,ls='--')
	axs2.set_xlim(-3,80)
	axs2.set_yticks([])
	axs2.set_title('HHE (filter: 4-16 Hz)',fontsize=20)
	axs2.set_xlabel('Tempo após a onda P (s)',fontsize=20)

	# format the ticks
	axs2.xaxis.set_major_locator(MultipleLocator(20))
	axs2.xaxis.set_minor_locator(MultipleLocator(5))
	axs2.xaxis.set_tick_params(labelsize=15)

	folder_name = OUTPUT_FIGURE_DIR+'EVENTS/Local/Seismogram/'
	os.makedirs(folder_name,exist_ok=True)
	fig.savefig(folder_name+'Seismogram_Event_'+sta_ev_name+'.png')
