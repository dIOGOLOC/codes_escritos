'''
Script to download and trim local event data 
(https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.trim.html#obspy.core.stream.Stream.trim)
'''


import os
import glob
import obspy as op
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.io.sac.sactrace import SACTrace
from obspy.clients.arclink.client import Client


from parameters_py.config import (
					TAUPY_MODEL,USER,HOST,PORT,
					CUT_BEFORE_P_LOCAL,CUT_AFTER_P_LOCAL,
					DIR_DATA,OUTPUT_EV_DIR,OUTPUT_FIGURE_DIR
							   )

#Function to cut and plot event:
			
def cut_download_data_by_event(knetwk,kstnm,stla,stlo,evla,evlo,evdp,evmag,ev_timeUTC):
	client = Client(user=USER,host=HOST, port=PORT)
	#Calculating distance, azimuth and backazimuth
	dist,az,baz = op.geodetics.gps2dist_azimuth(evla,evlo,stla,stlo)
	gcarc = op.geodetics.kilometer2degrees(dist/1000)

	
	#Calculating ray parameter
	model = TauPyModel(model=TAUPY_MODEL)
	arrivals = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=gcarc, phase_list=["P"])
	arr = arrivals[0]

	#Reference time
	starttime = (op.UTCDateTime(ev_timeUTC)+arr.time)-CUT_BEFORE_P_LOCAL
	endtime = (op.UTCDateTime(ev_timeUTC)+arr.time)+CUT_AFTER_P_LOCAL	

	########################################################################################################################################################
	#STREAM 

	#Creating Event Directory 
	event_directory = OUTPUT_EV_DIR+'Local/'+knetwk+'/'+kstnm+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]
	os.makedirs(event_directory, exist_ok=True)
	print('Event Directory - ',event_directory)
			
	#-----------------------------------
	#Component E
	stE = client.get_waveforms(knetwk,kstnm,'00','HHE',starttime,endtime)
	headerHHX = {
				'kstnm': kstnm, 'kcmpnm': 'HHE', 
				'stla': float(stla), 'stlo': float(stlo), 
				'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag), 
				'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute),'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]), 'nzsec': int(starttime.second), 'nzyear': int(starttime.year),
				'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(dist/1000), 'gcarc': float(gcarc),'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),
				'o':float(CUT_BEFORE_P_LOCAL),
				'delta':stE[0].stats.delta
				}

	sacHHX = SACTrace(data=stE[0].data, **headerHHX)
	sacHHX.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.E')
								
	#-----------------------------------
	#Component N
				
	stN = client.get_waveforms(knetwk,kstnm,'00','HHN',starttime,endtime)

	headerHHY = {
				'kstnm': kstnm, 'kcmpnm': 'HHN',
				'stla': float(stla), 'stlo': float(stlo),
				'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag), 
				'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),
				'cmpaz': 0.0, 'cmpinc': 90.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),
				'o':float(CUT_BEFORE_P_LOCAL),
				'delta':stN[0].stats.delta
				}

	sacHHY = SACTrace(data=stN[0].data, **headerHHY)
	sacHHY.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.N')
				
	#-----------------------------------
	#Component Z
	stZ = client.get_waveforms(knetwk,kstnm,'00','HHZ',starttime,endtime)

	headerHHZ = {
				'kstnm': kstnm, 'kcmpnm': 'HHZ', 
				'stla': float(stla), 'stlo': float(stlo), 
				'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag), 
				'nzhour': int(starttime.hour),'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),	
				'cmpaz': 0.0, 'cmpinc': 0.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),
				'o':float(CUT_BEFORE_P_LOCAL),
				'delta':stZ[0].stats.delta
						}

	sacHHZ = SACTrace(data=stZ[0].data, **headerHHZ)
	sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.Z')

def cut_data_by_local_event(knetwk,kstnm,stla,stlo,evla,evlo,evdp,evmag,ev_timeUTC):
	data_sta = DIR_DATA+knetwk+'/'+kstnm+'/'
	print(data_sta)
	if os.path.isdir(data_sta) == True:

		#Calculating distance, azimuth and backazimuth
		dist,az,baz = op.geodetics.gps2dist_azimuth(evla,evlo,stla,stlo)
		gcarc = op.geodetics.kilometer2degrees(dist/1000)
		
		#Reference time
		starttime = op.UTCDateTime(ev_timeUTC)-CUT_BEFORE_P_LOCAL
		endtime = op.UTCDateTime(ev_timeUTC)+CUT_AFTER_P_LOCAL	

		########################################################################################################################################################
		#STREAM 

		#-----------------------------------
		#Component E

		if os.path.isfile(data_sta+'HHE.D'+'/'+knetwk+'.'+kstnm+'..HHE.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)) == True:
		
			stE = op.read(data_sta+'HHE.D'+'/'+knetwk+'.'+kstnm+'..HHE.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
			stE.trim(starttime,endtime)		

			#Creating Event Directory 
			event_directory = OUTPUT_EV_DIR+'Local/'+knetwk+'/'+kstnm+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]
			os.makedirs(event_directory, exist_ok=True)
			print('Event Directory - ',event_directory)

			headerHHE = {
						'kstnm': kstnm, 'kcmpnm': 'HHE', 
						'stla': float(stla), 'stlo': float(stlo), 
						'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag), 
						'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute),'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]), 'nzsec': int(starttime.second), 'nzyear': int(starttime.year),
						'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(dist/1000), 'gcarc': float(gcarc),'az': float(az), 'baz': float(baz),
						'o':float(CUT_BEFORE_P_LOCAL),
						'delta':stE[0].stats.delta
						}

			sacHHE = SACTrace(data=stE[0].data, **headerHHE)
			sacHHE.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.E')
									
			#-----------------------------------
			#Component N
						
			stN = op.read(data_sta+'HHN.D'+'/'+knetwk+'.'+kstnm+'..HHN.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
			stN.trim(starttime,endtime)		

			headerHHY = {
						'kstnm': kstnm, 'kcmpnm': 'HHN',
						'stla': float(stla), 'stlo': float(stlo),
						'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag), 
						'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),
						'cmpaz': 0.0, 'cmpinc': 90.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
						'o':float(CUT_BEFORE_P_LOCAL),
						'delta':stN[0].stats.delta
						}

			sacHHY = SACTrace(data=stN[0].data, **headerHHY)
			sacHHY.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.N')
						
			#-----------------------------------
			#Component Z
			stZ = op.read(data_sta+'HHZ.D'+'/'+knetwk+'.'+kstnm+'..HHZ.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
			stZ.trim(starttime,endtime)	
			
			headerHHZ = {
						'kstnm': kstnm, 'kcmpnm': 'HHZ', 
						'stla': float(stla), 'stlo': float(stlo), 
						'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag), 
						'nzhour': int(starttime.hour),'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),	
						'cmpaz': 0.0, 'cmpinc': 0.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
						'o':float(CUT_BEFORE_P_LOCAL),
						'delta':stZ[0].stats.delta
						}

			sacHHZ = SACTrace(data=stZ[0].data, **headerHHZ)
			sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.Z')

			#-----------------------------------
			#Component X
			stX = op.read(data_sta+'HHX.D'+'/'+knetwk+'.'+kstnm+'..HHX.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
			stX.trim(starttime,endtime)	
			
			headerHHX = {
						'kstnm': kstnm, 'kcmpnm': 'HHX', 
						'stla': float(stla), 'stlo': float(stlo), 
						'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag), 
						'nzhour': int(starttime.hour),'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),	
						'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
						'o':float(CUT_BEFORE_P_LOCAL),
						'delta':stX[0].stats.delta
						}

			sacHHX = SACTrace(data=stX[0].data, **headerHHX)
			sacHHX.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.X')

def plot_event_data(direc):
	os.chdir(direc)
	event_date = direc.split('/')[-1]
	stZ = op.read('*HHZ*')
	
	fig, axes = plt.subplots(len(stZ),2, sharex=True,figsize=(20, 15))

	cols = ['Raw Data', 'Filterd Data (2 Hz to 10 Hz)']
	
	for ax, col in zip(axes[0], cols):
		ax.set_title(col)

	for i,j in enumerate(stZ):
		axes[i,0].plot(j.times(),j.data,'k')
		axes[i,0].set_xlim(5,45)
	axes[i,0].set_xlabel('Time after P (s)')

	for i,j in enumerate(stZ):
		j.filter('bandpass',freqmin=2.0, freqmax=10.0)
		axes[i,1].set_xlim(5,45)
		if j.stats.station in ['9FE7','9FF5','9FF9','A031','A062','9F74','9FF3','9FF6','A002','A037','A06A']:
			axes[i,1].plot(j.times(),[i*(-1) for i in j.data],'k')
		else:
			axes[i,1].plot(j.times(),j.data,'k')
		axes[i,1].text(45.5,0,j.stats.station)
		#axes[i,1].set_ylim(-0.000001,0.000001)
	axes[i,1].set_xlabel('Time after P (s)')
	fig.suptitle('Event - '+event_date)
	os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/',exist_ok=True)
	fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/Event - '+event_date+'.pdf')
	plt.tight_layout()
	plt.show()

def plot_event_dataset(direc):
	os.chdir(direc)
	event_date = direc.split('/')[-1]
	stZ = op.read('*HHZ*')
	
	fig, axes = plt.subplots(1,1, sharex=True,figsize=(20, 15))

	for i,j in enumerate(stZ):
		j.filter('bandpass',freqmin=2.0, freqmax=10.0)
		axes.set_xlim(5,45)
		if j.stats.station in ['9FE7','9FF5','9FF9','A031','A062','9F74','9FF3','9FF6','A002','A037','A06A']:
			axes.plot(j.times(),[i*(-1) for i in j.data],alpha=0.5,label=j.stats.station)
		else:
			axes.plot(j.times(),j.data,alpha=0.1,label=j.stats.station)
		axes.text(45.5,0,j.stats.station)
		axes.set_ylim(-0.000001,0.000001)
	axes.set_xlabel('Time after P (s)')
	axes.legend(loc=0)
	fig.suptitle('Event - '+event_date)
	os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/',exist_ok=True)
	fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/Event - '+event_date+'_DATASET.pdf')
	plt.tight_layout()
	plt.show()