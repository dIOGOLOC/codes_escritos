'''
Script to trim event data of large earthquakes 
(https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.trim.html#obspy.core.stream.Stream.trim)
and plot a mosaic of event raw and filtered data.
'''


import os
import glob
import obspy as op
from obspy import read,read_inventory, UTCDateTime, Stream
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.io.sac.sactrace import SACTrace


from parameters_py.config import (
					DIR_DATA,TAUPY_MODEL,EV_GCARC_MIN,EV_GCARC_MAX,CUT_BEFORE_P,CUT_AFTER_P,XML_FILE,
					OUTPUT_EV_DIR,OUTPUT_FIGURE_DIR
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

		for tr in st:
			tr.remove_response(inventory=inv,output="VEL",water_level=60)
		
		fig, axes = plt.subplots(len(st),2, sharex=True,figsize=(20, 15))

		cols = ['Raw Data', 'Filterd Data ('+str(lf)+' Hz to '+str(hf)+' Hz)']
		
		for ax, col in zip(axes[0], cols):
			ax.set_title(col)

		for i,j in enumerate(st):
			axes[i,0].plot(j.times(),j.data,'k')
			axes[i,0].set_xlim(5,100)
		axes[i,0].set_xlabel('Time after P (s)')

		for i,j in enumerate(st):
			j.filter('bandpass',freqmin=lf, freqmax=hf)
			axes[i,1].set_xlim(5,100)
			axes[i,1].plot(j.times(),j.data,'k')
			axes[i,1].text(100.5,0,j.stats.station)
		axes[i,1].set_xlabel('Time after P (s)')
		fig.suptitle('Event - '+event_name)
		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'Event - '+event_name+'.pdf')
		plt.tight_layout()

	else:
		st = read(event_lst[0])
		st.detrend("linear")
		st.detrend("demean")
		st.taper(max_percentage=0.05, type="hann")

		for tr in st:
			tr.remove_response(inventory=inv,output="VEL",water_level=60)
		
		fig, axes = plt.subplots(1,2, sharex=True,figsize=(20, 15))

		cols = ['Raw Data', 'Filterd Data ('+str(lf)+' Hz to '+str(hf)+' Hz)']
		
		for i,j in enumerate(st):
			axes[0].plot(j.times(),j.data,'k')
			axes[0].set_xlim(5,100)
		axes[0].set_xlabel('Time after P (s)')
		axes[0].set_title(cols[0])

		for i,j in enumerate(st):
			j.filter('bandpass',freqmin=lf, freqmax=hf)
			axes[1].set_xlim(5,100)
			axes[1].plot(j.times(),j.data,'k')
			axes[1].text(100.5,0,j.stats.station)
		axes[1].set_xlabel('Time after P (s)')
		axes[1].set_title(cols[1])

		fig.suptitle('Event - '+event_name)
		os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name,exist_ok=True)
		fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/'+folder_name+'Event - '+event_name+'.pdf')
		plt.tight_layout()