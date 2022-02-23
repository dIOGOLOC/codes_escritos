'''
--------------------------------------------------------------------------------
 Function to trim/plot local the dataset according to local events time
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 02/2022


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
This code will trim and plot the local datase according to a given an event time
and a list of stations.

More information in:
https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.trim.html


Inputs:
JSON file with event description:
    ev_timeUTC: event time in UTC (str)
    ev_year: year of the event
    ev_month: month of the event
    ev_day: day of the event
    ev_julday: julian day of the event
    ev_hour: hour of the event
    ev_minute: minute of the event
    ev_second: second of the event
    ev_microsecond: microsecond of the event
    evla: latitude of the event
    evlo: longitude of the event
    evdp: depth of the event
    mag: magnitude of the event

'''


import os
import glob
import obspy as op
from obspy import read,read_inventory, UTCDateTime, Stream
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.io.sac.sactrace import SACTrace
from obspy.signal.trigger import classic_sta_lta, trigger_onset, coincidence_trigger,recursive_sta_lta,plot_trigger

import numpy as np
import pandas as pd
import geopandas as geopd
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.transforms import offset_copy

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter,LatitudeLocator,LongitudeLocator

from parameters_py.config import (
					DIR_DATA,TAUPY_MODEL,EV_GCARC_MIN,EV_GCARC_MAX,CUT_BEFORE_P,CUT_AFTER_P,XML_FILE,
					OUTPUT_EV_DIR,OUTPUT_FIGURE_DIR,BOUNDARY_STATES_SHP,LOCAL_EV_DISTANCE_MIN
							   )
# ===============================
# Function to cut and plot event:
# ===============================

def cut_data_by_event(knetwk,kstnm,stla,stlo,ev_timeUTC,ev_julday,ev_year,ev_month,ev_day,ev_hour,ev_minute,ev_second,ev_microsecond,ev_lat,ev_long,ev_depth,ev_mag):
	#Calculating distance, azimuth and backazimuth
	dist,az,baz = op.geodetics.gps2dist_azimuth(ev_lat,ev_long,stla,stlo)
	gcarc = op.geodetics.kilometer2degrees(dist/1000)
	if EV_GCARC_MIN <= gcarc <= EV_GCARC_MAX:

		#Calculating ray parameter
		model = TauPyModel(model=TAUPY_MODEL)
		arrivals = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=gcarc, phase_list=["P"])
		arr = arrivals[0]

		#Reference time
		starttime = (op.UTCDateTime(ev_timeUTC)+arr.time)-CUT_BEFORE_P
		endtime = (op.UTCDateTime(ev_timeUTC)+arr.time)+CUT_AFTER_P

		########################################################################################################################################################
		#STREAM

		#-----------------------------------
		#Component E

		if os.path.isfile(DIR_DATA+knetwk+'/'+kstnm+'/HHE.D'+'/'+knetwk+'.'+kstnm+'..HHE.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)) == True:

			#Creating Event Directory
			event_directory = OUTPUT_EV_DIR+'Regional/'+knetwk+'/'+kstnm+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]
			os.makedirs(event_directory, exist_ok=True)

			stE = read(DIR_DATA+knetwk+'/'+kstnm+'/HHE.D'+'/'+knetwk+'.'+kstnm+'..HHE.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
			stE.trim(starttime,endtime)

			headerHHE = {
						'kstnm': kstnm, 'kcmpnm': 'HHE','knetwk':knetwk,
						'stla': float(stla), 'stlo': float(stlo),
						'evdp': float(ev_depth), 'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag),
						'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute),'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]), 'nzsec': int(starttime.second), 'nzyear': int(starttime.year),
						'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(dist/1000), 'gcarc': float(gcarc),'az': float(az), 'baz': float(baz),
						'o':float(CUT_BEFORE_P),
						'delta':stE[0].stats.delta
						}

			sacHHE = SACTrace(data=stE[0].data, **headerHHE)
			sacHHE.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.E')

		#-----------------------------------
		#Component N
		if os.path.isfile(DIR_DATA+knetwk+'/'+kstnm+'/HHN.D'+'/'+knetwk+'.'+kstnm+'..HHN.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)) == True:

			#Creating Event Directory
			event_directory = OUTPUT_EV_DIR+'Regional/'+knetwk+'/'+kstnm+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]
			os.makedirs(event_directory, exist_ok=True)

			stN = read(DIR_DATA+knetwk+'/'+kstnm+'/HHN.D'+'/'+knetwk+'.'+kstnm+'..HHN.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
			stN.trim(starttime,endtime)

			headerHHY = {
						'kstnm': kstnm, 'kcmpnm': 'HHN','knetwk':knetwk,
						'stla': float(stla), 'stlo': float(stlo),
						'evdp': float(ev_depth), 'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag),
						'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),
						'cmpaz': 0.0, 'cmpinc': 90.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
						'o':float(CUT_BEFORE_P),
						'delta':stN[0].stats.delta
						}

			sacHHY = SACTrace(data=stN[0].data, **headerHHY)
			sacHHY.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.N')

		#-----------------------------------
		#Component Z

		if os.path.isfile(DIR_DATA+knetwk+'/'+kstnm+'/HHZ.D'+'/'+knetwk+'.'+kstnm+'..HHZ.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)) == True:

			#Creating Event Directory
			event_directory = OUTPUT_EV_DIR+'Regional/'+knetwk+'/'+kstnm+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]
			os.makedirs(event_directory, exist_ok=True)

			stZ = read(DIR_DATA+knetwk+'/'+kstnm+'/HHZ.D'+'/'+knetwk+'.'+kstnm+'..HHZ.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
			stZ.trim(starttime,endtime)

			headerHHZ = {
						'kstnm': kstnm, 'kcmpnm': 'HHZ','knetwk':knetwk,
						'stla': float(stla), 'stlo': float(stlo),
						'evdp': float(ev_depth), 'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag),
						'nzhour': int(starttime.hour),'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),
						'cmpaz': 0.0, 'cmpinc': 0.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
						'o':float(CUT_BEFORE_P),
						'delta':stZ[0].stats.delta
						}

			sacHHZ = SACTrace(data=stZ[0].data, **headerHHZ)
			sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.Z')

		#-----------------------------------
		#Component X
		if os.path.isfile(DIR_DATA+knetwk+'/'+kstnm+'/HHX.D'+'/'+knetwk+'.'+kstnm+'..HHX.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)) == True:

			#Creating Event Directory
			event_directory = OUTPUT_EV_DIR+'Regional/'+knetwk+'/'+kstnm+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]
			os.makedirs(event_directory, exist_ok=True)

			stX = read(DIR_DATA+knetwk+'/'+kstnm+'/HHX.D'+'/'+knetwk+'.'+kstnm+'..HHX.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
			stX.trim(starttime,endtime)

			headerHHX = {
						'kstnm': kstnm, 'kcmpnm': 'HHX','knetwk':knetwk,
						'stla': float(stla), 'stlo': float(stlo),
						'evdp': float(ev_depth), 'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag),
						'nzhour': int(starttime.hour),'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),
						'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
						'o':float(CUT_BEFORE_P),
						'delta':stX[0].stats.delta
						}

			sacHHX = SACTrace(data=stX[0].data, **headerHHX)
			sacHHX.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.X')


def plot_event_data(event_lst, inv, event_name,folder_name,lf,hf):
	if len(event_lst) > 1:
		st = Stream()
		for i in event_lst:
			st += read(i)

		st.detrend("linear")
		st.detrend("demean")
		st.taper(max_percentage=0.05, type="hann")

		if st[0].stats.channel != 'HHX':
			for tr in st:
				pre_filt = [0.001, 0.005, 45., 50.]
				tr.remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)

		fig, axes = plt.subplots(len(st),2, sharex=True,figsize=(20, 15))

		cols = ['Raw Data', 'Filterd Data ('+str(lf)+' Hz to '+str(hf)+' Hz)']
		cols_br = ['Dado bruto', 'Dado filtrado ('+str(lf)+' Hz to '+str(hf)+' Hz)']

		for ax, col in zip(axes[0], cols):
			ax.set_title(col)
			if LABEL_LANG == 'br':
				ax.set_title(cols_br)

			else:
				ax.set_title(col)

		for i,j in enumerate(st):
			axes[i,0].plot(j.times(),j.data,'k')
			axes[i,0].set_xlim(5,100)

		if LABEL_LANG == 'br':
			axes[i,0].set_xlabel('Tempos após P (s)')

		else:
			axes[i,0].set_xlabel('Time after P (s)')

		for i,j in enumerate(st):
			j.filter('bandpass',freqmin=lf, freqmax=hf)
			axes[i,1].set_xlim(5,100)
			axes[i,1].plot(j.times(),j.data,'k')
			axes[i,1].text(100.5,0,j.stats.station)

		if LABEL_LANG == 'br':
			axes[i,1].set_xlabel('Tempos após P (s)')
			fig.suptitle('Evento: '+event_name)

		else:
			axes[i,1].set_xlabel('Time after P (s)')
			fig.suptitle('Event: '+event_name)

		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'Event - '+event_name+'.png')
		plt.tight_layout()

	else:
		st = read(event_lst[0])
		st.detrend("linear")
		st.detrend("demean")
		st.taper(max_percentage=0.05, type="hann")

		if st[0].stats.channel != 'HHX':
			for tr in st:
				pre_filt = [0.001, 0.005, 45., 50.]
				tr.remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)

		fig, axes = plt.subplots(1,2, sharex=True,figsize=(20, 15))

		cols = ['Raw Data', 'Filterd Data ('+str(lf)+' Hz to '+str(hf)+' Hz)']
		cols_br = ['Dado bruto', 'Dado filtrado ('+str(lf)+' Hz to '+str(hf)+' Hz)']

		for i,j in enumerate(st):
			axes[0].plot(j.times(),j.data,'k')
			axes[0].set_xlim(5,100)

		if LABEL_LANG == 'br':
			axes[0].set_xlabel('Tempos após P (s)')
			axes[0].set_title(cols_br[0])

		else:
			axes[0].set_xlabel('Time after P (s)')
			axes[0].set_title(cols[0])

		for i,j in enumerate(st):
			j.filter('bandpass',freqmin=lf, freqmax=hf)
			axes[1].set_xlim(5,100)
			axes[1].plot(j.times(),j.data,'k')
			axes[1].text(100.5,0,j.stats.station)

		if LABEL_LANG == 'br':
			axes[1].set_xlabel('Tempos após P (s)')
			axes[1].set_title(cols_br[1])
			fig.suptitle('Evento: '+event_name)

		else:
			axes[1].set_xlabel('Time after P (s)')
			axes[1].set_title(cols[1])
			fig.suptitle('Event: '+event_name)

		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'Event - '+event_name+'.png')
		plt.tight_layout()

def plot_map_event_data(event_lstZ,event_lstN,event_lstE,event_name,folder_name,lf,hf):
		stZ = Stream()
		for i in event_lstZ:
			temp_z = read(i)
			if round(temp_z[0].stats.sac.dist) < LOCAL_EV_DISTANCE_MIN:
				stZ += temp_z

		stZ.detrend("linear")
		stZ.detrend("demean")
		stZ.taper(max_percentage=0.05, type="hann")
		stZ.filter('bandpass',freqmin=lf, freqmax=hf)

		stN = Stream()
		for i in event_lstN:
			temp_N = read(i)
			if round(temp_N[0].stats.sac.dist) < LOCAL_EV_DISTANCE_MIN:
				stN += read(i)

		stN.detrend("linear")
		stN.detrend("demean")
		stN.taper(max_percentage=0.05, type="hann")
		stN.filter('bandpass',freqmin=lf, freqmax=hf)

		stE = Stream()
		for i in event_lstE:
			temp_E = read(i)
			if round(temp_E[0].stats.sac.dist) < LOCAL_EV_DISTANCE_MIN:
				stE += read(i)

		stE.detrend("linear")
		stE.detrend("demean")
		stE.taper(max_percentage=0.05, type="hann")
		stE.filter('bandpass',freqmin=lf, freqmax=hf)

		# -----------------------------------------------------------------------
		# MAP
		# -----------------------------------------------------------------------

		fig = plt.figure(figsize=(10, 10))
		event_date = event_name.split('.')
		gs = gridspec.GridSpec(nrows=1, ncols=1)
		#-------------------------------------------
		crs = ccrs.Orthographic(central_longitude=-40, central_latitude=-20)
		map_loc = fig.add_subplot(gs[0],projection=crs)

		LLCRNRLON_LARGE = -52
		URCRNRLON_LARGE = -28
		LLCRNRLAT_LARGE = -30
		URCRNRLAT_LARGE = -12

		# Create a Stamen Terrain instance.
		stamen_terrain = cimgt.StamenTerrain()

		# Add the Stamen data at zoom level 8.
		map_loc.add_image(stamen_terrain, 10)

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

		#Plotting event
		map_loc.plot(stZ[0].stats.sac.evlo, stZ[0].stats.sac.evla, marker='*', color='red', markersize=20,transform=ccrs.PlateCarree())
		map_loc.text(stZ[0].stats.sac.evlo, stZ[0].stats.sac.evla, str(stZ[0].stats.sac.mag),fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform_mag)

		#Plotting stations
		stlo = []
		stla = []
		station_name = []
		gcarc_lst = []
		for i,j in enumerate(stZ):
			    stlo.append(j.stats.sac.stlo)
			    stla.append(j.stats.sac.stla)
			    station_name.append(j.stats.station)
			    gcarc_lst.append(j.stats.sac.gcarc)

		gcarc_sort = np.argsort(gcarc_lst)

		map_loc.scatter(stlo, stla, marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

		for i,j in enumerate(station_name):
			map_loc.text(stlo[i], stla[i], station_name[i],fontsize=12,verticalalignment='center', horizontalalignment='right',transform=text_transform,bbox=dict(facecolor='white',edgecolor='none', alpha=0.5, boxstyle='round'))

		#Saving figure
		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'MAP_Event - '+event_name+'.png')

		# -----------------------------------------------------------------------
		# Waveforms
		# -----------------------------------------------------------------------

		fig, axes = plt.subplots(nrows=len(stZ), ncols=3, sharex='col',figsize=(20, 20))
		event_date = event_name.split('.')
		fig.suptitle('Dia do Evento - '+UTCDateTime(year=int(event_date[0]),julday=int(event_date[1])).strftime('%d/%m/%Y')+' - Magnitude:'+str(stZ[0].stats.sac.mag),fontsize=20)

		#-------------------------------------------

		minutes = mdates.MinuteLocator(3)   # every minute
		#seconds = mdates.SecondLocator(1)  # every second
		seconds = mdates.MinuteLocator(1)  # every second
		years_fmt = mdates.DateFormatter('%H-%M-%S')

		for i,j in enumerate(gcarc_sort):
				ax = axes[i,0]
				ax.plot(stZ[j].times("matplotlib"),stZ[j].data,color='k')

				ax.set_xlim((stZ[j].times("utcdatetime")[0]-100).matplotlib_date,(stZ[j].times("utcdatetime")[0]+300).matplotlib_date)
				ax.set_yticks([])

				if LABEL_LANG == 'br':
					ax.set_title(stZ[j].id+' - dist='+str(int(round(stZ[j].stats.sac.dist)))+' km',fontsize=15)

				else:
					ax.set_title(stZ[j].id+' - dist='+str(int(round(stZ[j].stats.sac.dist)))+' km',fontsize=15)

				# format the ticks
				ax.xaxis.set_major_locator(minutes)
				ax.xaxis.set_major_formatter(years_fmt)
				ax.xaxis.set_minor_locator(seconds)

				#--------------------------------------------------------------
				ax = axes[i,1]
				ax.plot(stN[j].times('matplotlib'),stN[j].data,color='k')

				ax.set_xlim((stN[j].times("utcdatetime")[0]-100).matplotlib_date,(stN[j].times("utcdatetime")[0]+300).matplotlib_date)
				ax.set_yticks([])
				ax.set_title(stN[j].id+' - dist='+str(int(round(stN[j].stats.sac.dist)))+' km',fontsize=15)

				# format the ticks
				ax.xaxis.set_major_locator(minutes)
				ax.xaxis.set_major_formatter(years_fmt)
				ax.xaxis.set_minor_locator(seconds)

				#--------------------------------------------------------------
				ax = axes[i,2]
				ax.plot(stE[j].times('matplotlib'),stE[j].data,color='k')

				ax.set_xlim((stE[j].times("utcdatetime")[0]-100).matplotlib_date,(stE[j].times("utcdatetime")[0]+300).matplotlib_date)
				ax.set_yticks([])
				ax.set_title(stE[j].id+' - dist='+str(int(round(stE[j].stats.sac.dist)))+' km',fontsize=15)

				# format the ticks
				ax.xaxis.set_major_locator(minutes)
				ax.xaxis.set_major_formatter(years_fmt)
				ax.xaxis.set_minor_locator(seconds)

		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'Stations_Event_'+event_name+'_mag_'+str(stZ[0].stats.sac.mag)+'.png')


		# -----------------------------------------------------------------------
		# LTA/STA
		# -----------------------------------------------------------------------

		stalen = 4
		ltalen = 10
		trig_on = 1.5
		trig_off = 0.5

		npts = stZ[0].stats.npts
		dt = stZ[0].stats.delta
		df = stZ[0].stats.sampling_rate

		# -----------------------------------------------------------------------

		fig, axes = plt.subplots(nrows=len(stZ), ncols=3, sharex='col',figsize=(20, 20))
		event_date = event_name.split('.')
		if LABEL_LANG == 'br':
			fig.suptitle('Dia do Evento: '+UTCDateTime(year=int(event_date[0]),julday=int(event_date[1])).strftime('%d/%m/%Y')+' - Magnitude:'+str(stZ[0].stats.sac.mag),fontsize=20)

		else:
			fig.suptitle('Event day: '+UTCDateTime(year=int(event_date[0]),julday=int(event_date[1])).strftime('%d/%m/%Y')+' - Magnitude:'+str(stZ[0].stats.sac.mag),fontsize=20)

		#-------------------------------------------

		minutes = mdates.MinuteLocator(3)   # every minute
		#seconds = mdates.SecondLocator(1)  # every second
		seconds = mdates.MinuteLocator(1)  # every second
		years_fmt = mdates.DateFormatter('%H-%M-%S')
		for i,j in enumerate(gcarc_sort):
				ax = axes[i,0]
				# format the ticks
				ax.xaxis.set_major_locator(minutes)
				ax.xaxis.set_major_formatter(years_fmt)
				ax.xaxis.set_minor_locator(seconds)

				cft_Z = classic_sta_lta(stZ[j].data, int(stalen * df), int(ltalen * df))
				on_off_Z = np.array(trigger_onset(cft_Z, trig_on, trig_off))

				#plot_trigger(stZ[j], cft_Z, trig_on, trig_off, show=True)

				on_Z = float((on_off_Z[:,0]/df)[0])
				off_Z = float((on_off_Z[:,1]/df)[0])

				time_on_Z = stZ[j].stats.starttime+on_Z
				time_off_Z = stZ[j].stats.starttime+off_Z

				ax.plot(stZ[j].times('matplotlib'), cft_Z, 'k')
				ax.axhline(trig_on, color='red', lw=1, ls='--')
				ax.axhline(trig_off, color='blue', lw=1, ls='--')

				#--------------------------------------------------------------
				ax = axes[i,1]
				# format the ticks
				ax.xaxis.set_major_locator(minutes)
				ax.xaxis.set_major_formatter(years_fmt)
				ax.xaxis.set_minor_locator(seconds)

				cft_N = classic_sta_lta(stN[j].data, int(stalen * df), int(ltalen * df))
				on_off_N= np.array(trigger_onset(cft_N, trig_on, trig_off))

				on_N = float((on_off_N[:,0]/df)[0])
				off_N = float((on_off_N[:,1]/df)[0])

				time_on_N = stN[j].stats.starttime+on_N
				time_off_N = stN[j].stats.starttime+off_N

				ax.plot(stN[j].times('matplotlib'), cft_N, 'k')
				ax.axhline(trig_on, color='red', lw=1, ls='--')
				ax.axhline(trig_off, color='blue', lw=1, ls='--')

				#--------------------------------------------------------------
				ax = axes[i,2]
				# format the ticks
				ax.xaxis.set_major_locator(minutes)
				ax.xaxis.set_major_formatter(years_fmt)
				ax.xaxis.set_minor_locator(seconds)

				cft_E = classic_sta_lta(stE[j].data, int(stalen * df), int(ltalen * df))
				on_off_E = np.array(trigger_onset(cft_E, trig_on, trig_off))

				on_E = float((on_off_E[:,0]/df)[0])
				off_E = float((on_off_E[:,1]/df)[0])

				time_on_E = stE[j].stats.starttime+on_E
				time_off_E = stE[j].stats.starttime+off_E

				ax.plot(stE[j].times('matplotlib'), cft_E, 'k')
				ax.axhline(trig_on, color='red', lw=1, ls='--')
				ax.axhline(trig_off, color='blue', lw=1, ls='--')

		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'Stations_Event_LTA_STA_'+event_name+'_mag_'+str(stZ[0].stats.sac.mag)+'.png')

def plot_map_event_data_hydrophone(event_lstX,event_name,folder_name,lf,hf):
	if len(event_lstX) > 1:
		stZ = Stream()
		for i in event_lstX:
			temp_z = read(i)
			stZ += temp_z

		stZ.detrend("linear")
		stZ.detrend("demean")
		stZ.taper(max_percentage=0.05, type="hann")
		stZ.filter('bandpass',freqmin=lf, freqmax=hf)

		#-------------------------------------------

		fig = plt.figure(figsize=(10, 10))
		event_date = event_name.split('.')
		gs = gridspec.GridSpec(nrows=1, ncols=1)
		#-------------------------------------------
		crs = ccrs.Orthographic(central_longitude=-40, central_latitude=-20)
		map_loc = fig.add_subplot(gs[0],projection=crs)

		LLCRNRLON_LARGE = -52
		URCRNRLON_LARGE = -28
		LLCRNRLAT_LARGE = -30
		URCRNRLAT_LARGE = -12

		# Create a Stamen Terrain instance.
		stamen_terrain = cimgt.Stamen('terrain-background')

		# Add the Stamen data at zoom level 8.
		map_loc.add_image(stamen_terrain, 10)

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
		map_loc.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True, labelsize=15)


		# Use the cartopy interface to create a matplotlib transform object
		# for the Geodetic coordinate system. We will use this along with
		# matplotlib's offset_copy function to define a coordinate system which
		# translates the text by 25 pixels to the left.
		geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
		text_transform = offset_copy(geodetic_transform, units='dots', y=17,x=33)
		text_transform_mag = offset_copy(geodetic_transform, units='dots', y=-25,x=20)

		stlo = []
		stla = []
		station_name = []
		gcarc_lst = []
		for i,j in enumerate(stZ):
			    stlo.append(j.stats.sac.stlo)
			    stla.append(j.stats.sac.stla)
			    station_name.append(j.stats.station)
			    gcarc_lst.append(j.stats.sac.gcarc)

		gcarc_sort = np.argsort(gcarc_lst)

		map_loc.scatter(stlo, stla, marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

		for i,j in enumerate(station_name):
			map_loc.text(stlo[i], stla[i], station_name[i],fontsize=12,verticalalignment='center', horizontalalignment='right',transform=text_transform,bbox=dict(facecolor='white',edgecolor='none', alpha=0.5, boxstyle='round'))

		map_loc.plot(stZ[0].stats.sac.evlo, stZ[0].stats.sac.evla, marker='*', color='red', markersize=20,transform=ccrs.PlateCarree())
		map_loc.text(stZ[0].stats.sac.evlo, stZ[0].stats.sac.evla, str(stZ[0].stats.sac.mag),fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform_mag)

		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'MAP_Event_hydrophone_'+event_name+'.png')

		# -----------------------------------------------------------------------
		# Waveforms
		# -----------------------------------------------------------------------

		fig, axes = plt.subplots(nrows=len(stZ), ncols=1, sharex='col',figsize=(20, 20))
		event_date = event_name.split('.')

		if LABEL_LANG == 'br':
			fig.suptitle('Dia do Evento - '+UTCDateTime(year=int(event_date[0]),julday=int(event_date[1])).strftime('%d/%m/%Y')+' - Magnitude:'+str(stZ[0].stats.sac.mag),fontsize=20)

		else:
			fig.suptitle('Event day: '+UTCDateTime(year=int(event_date[0]),julday=int(event_date[1])).strftime('%d/%m/%Y')+' - Magnitude:'+str(stZ[0].stats.sac.mag),fontsize=20)


		#-------------------------------------------

		minutes = mdates.MinuteLocator(2)   # every year
		#seconds = mdates.SecondLocator(60)  # every second
		seconds = mdates.MinuteLocator(1)  # every second
		years_fmt = mdates.DateFormatter('%H-%M-%S')

		for i,j in enumerate(gcarc_sort):
			ax = axes[i]
			ax.plot(stZ[j].times("matplotlib"),stZ[j].data,color='k')

			ax.set_xlim((stZ[j].times("utcdatetime")[0]-50).matplotlib_date,(stZ[j].times("utcdatetime")[0]+200).matplotlib_date)
			ax.set_yticks([])
			ax.set_title(stZ[j].id+' - dist='+str(int(round(stZ[j].stats.sac.dist)))+' km',fontsize=15)

			# format the ticks
			ax.xaxis.set_major_locator(minutes)
			ax.xaxis.set_major_formatter(years_fmt)
			ax.xaxis.set_minor_locator(seconds)

		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'Stations_Event_hydrophone_'+event_name+'.png')

	#-------------------------------------------

	if len(event_lstX) == 1:
		stZ = Stream()
		stZ += read(event_lstX[0])
		stZ.detrend("linear")
		stZ.detrend("demean")
		stZ.taper(max_percentage=0.05, type="hann")
		stZ.filter('bandpass',freqmin=lf, freqmax=hf)


		#-------------------------------------------

		fig = plt.figure(figsize=(10, 10))
		event_date = event_name.split('.')
		gs = gridspec.GridSpec(nrows=1, ncols=1)
		#-------------------------------------------
		crs = ccrs.Orthographic(central_longitude=-40, central_latitude=-20)
		map_loc = fig.add_subplot(gs[0],projection=crs)

		LLCRNRLON_LARGE = -52
		URCRNRLON_LARGE = -28
		LLCRNRLAT_LARGE = -30
		URCRNRLAT_LARGE = -12

		# Create a Stamen Terrain instance.
		stamen_terrain = cimgt.StamenTerrain()

		# Add the Stamen data at zoom level 8.
		map_loc.add_image(stamen_terrain, 10)

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

		gl.xlabel_style = {'size': 20, 'color': 'gray'}
		gl.xlabel_style = {'color': 'black', 'weight': 'bold'}
		gl.ylabel_style = {'size': 20, 'color': 'gray'}
		gl.ylabel_style = {'color': 'black', 'weight': 'bold'}

		map_loc.yaxis.set_ticks_position('both')
		map_loc.xaxis.set_ticks_position('both')
		map_loc.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True, labelsize=15)


		# Use the cartopy interface to create a matplotlib transform object
		# for the Geodetic coordinate system. We will use this along with
		# matplotlib's offset_copy function to define a coordinate system which
		# translates the text by 25 pixels to the left.
		geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
		text_transform = offset_copy(geodetic_transform, units='dots', y=17,x=33)
		text_transform_mag = offset_copy(geodetic_transform, units='dots', y=-25,x=20)

		stlo = stZ[0].stats.sac.stlo
		stla = stZ[0].stats.sac.stla
		station_name = stZ[0].stats.station
		gcarc_lst = stZ[0].stats.sac.gcarc

		map_loc.scatter(stlo, stla, marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

		map_loc.text(stlo, stla, station_name,fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform,bbox=dict(facecolor='white',edgecolor='none', alpha=0.5, boxstyle='round'))
		map_loc.plot(stZ[0].stats.sac.evlo, stZ[0].stats.sac.evla, marker='*', color='red', markersize=15,transform=ccrs.Geodetic())
		map_loc.text(stZ[0].stats.sac.evlo, stZ[0].stats.sac.evla, str(stZ[0].stats.sac.mag),fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform_mag)

		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'MAP_Event_hydrophone_'+event_name+'.png')

		#-------------------------------------------

		minutes = mdates.MinuteLocator()   # every year
		seconds = mdates.SecondLocator()  # every second
		years_fmt = mdates.DateFormatter('%H-%M-%S')

		ax = fig.add_subplot(gs[1])
		ax.plot(stZ[0].times("matplotlib"),stZ[0].data,color='k')

		ax.set_xlim((stZ[0].times("utcdatetime")[0]-10).matplotlib_date,(stZ[0].times("utcdatetime")[0]+200).matplotlib_date)
		ax.set_yticks([])
		ax.set_title(stZ[0].id+' - dist='+str(int(round(stZ[0].stats.sac.dist)))+' km',fontsize=15)

		# format the ticks
		ax.xaxis.set_major_locator(minutes)
		ax.xaxis.set_major_formatter(years_fmt)
		ax.xaxis.set_minor_locator(seconds)

		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'Stations_Event_hydrophone_'+event_name+'.png')
