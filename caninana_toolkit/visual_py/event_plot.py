'''
Script to trim event data of large earthquakes 
(https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.trim.html#obspy.core.stream.Stream.trim)
and plot a mosaic of event raw and filtered data.
'''


import os
import glob
import obspy as op
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt


from parameters_py.config import (
					DIR_DATA,TAUPY_MODEL,EV_GCARC_MIN,EV_GCARC_MAX,CUT_BEFORE_P,CUT_AFTER_P,OUTPUT_XML_FILE_DIR,
					NETWORK_CODE,LOCATION,OUTPUT_EV_DIR,OUTPUT_FIGURE_DIR
							   )

print('Importing XML file')

inv = op.read_inventory(OUTPUT_XML_FILE_DIR+NETWORK_CODE+".xml")
print(inv)

#Function to cut and plot event:
			
def cut_data_by_event(kstnm,stla,stlo,ev_timeUTC,ev_julday,ev_year,ev_month,ev_day,ev_hour,ev_minute,ev_second,ev_microsecond,ev_lat,ev_long,ev_depth,ev_mag):
	if os.path.isdir(DIR_DATA+kstnm+'/'+ev_year+'/'+ev_julday) == True:
			#Calculating distance, azimuth and backazimuth
			dist,az,baz = op.geodetics.gps2dist_azimuth(ev_lat,ev_long,stla,stlo)
			gcarc = op.geodetics.kilometer2degrees(dist/1000)
			if EV_GCARC_MIN <= gcarc <= EV_GCARC_MAX:
				print('Current Directory = '+DIR_DATA+kstnm+'/'+ev_year+'/'+ev_julday)	
				os.chdir(DIR_DATA+kstnm+'/'+ev_year+'/'+ev_julday)

				#Calculating ray parameter
				model = TauPyModel(model=TAUPY_MODEL)
				arrivals = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=gcarc, phase_list=["P"])
				arr = arrivals[0]

				#Reference time
				starttime = (op.UTCDateTime(ev_timeUTC)+arr.time)-CUT_BEFORE_P
				endtime = (op.UTCDateTime(ev_timeUTC)+arr.time)+CUT_AFTER_P	

				#########################################################################################################################################################
				#STREAM 
				
				HH = op.read('*')
				HH.merge()
				for i,j in enumerate(HH):
					j.stats.station = kstnm
					j.stats.network = NETWORK_CODE
					j.stats.location = LOCATION
					if j.stats.channel in ['HH1','HH1j','HHY']:
						j.stats.channel = 'HHN'
					elif j.stats.channel in ['HH2','HHX']:
						j.stats.channel = 'HHE'
					elif j.stats.channel in ['HHZ']:
						j.stats.channel = 'HHZ'

				HH.trim(starttime,endtime)
				HH.remove_response(inventory=inv,output='DISP')
				output_folder = OUTPUT_EV_DIR+str(int(starttime.year))+'.'+str(int(starttime.julday))+'.'+str(int(starttime.hour))+'.'+str(int(starttime.minute))+'.'+str(int(starttime.second))
				os.makedirs(output_folder,exist_ok=True)
				for tr in HH: 
					tr.write(output_folder+'/'+tr.id+'.'+str(int(starttime.year))+'.'+str(int(starttime.julday))+'.'+str(int(starttime.hour))+'.'+str(int(starttime.minute))+'.'+str(int(starttime.second))+".MSEED") 
			
	return 'No data for this event = '+DIR_DATA+kstnm+'/'+ev_year+'/'+ev_julday

def plot_event_data(direc):
	os.chdir(direc)
	event_date = direc.split('/')[-1]
	stZ = op.read('*HHZ*')
	fig, axes = plt.subplots(len(stZ),2, sharex=True,figsize=(20, 15))

	cols = ['Raw Data', 'Filterd Data (1 Hz to 20 Hz)']
	
	for ax, col in zip(axes[0], cols):
		ax.set_title(col)

	for i,j in enumerate(stZ):
		axes[i,0].plot(j.times(),j.data,'k')
		axes[i,0].set_xlim(5,45)
	axes[i,0].set_xlabel('Time after P (s)')

	for i,j in enumerate(stZ):
		j.filter('bandpass',freqmin=1.0, freqmax=20.0)
		axes[i,1].set_xlim(5,45)
		axes[i,1].plot(j.times(),j.data,'k')
		axes[i,1].text(45.5,0,j.stats.station)
		axes[i,1].set_ylim(-0.000001,0.000001)
	axes[i,1].set_xlabel('Time after P (s)')
	fig.suptitle('Event - '+event_date)
	os.makedirs(OUTPUT_FIGURE_DIR+'EVENTS/',exist_ok=True)
	fig.savefig(OUTPUT_FIGURE_DIR+'EVENTS/Event - '+event_date+'.pdf')
	plt.tight_layout()
	plt.show()

	'''
	stZ.filter('bandpass',freqmin=1.0, freqmax=20.0)
	stZ.plot()

	fig = plt.figure()
	stN = op.read('*HHN*')
	stN.filter('bandpass',freqmin=1.0, freqmax=20.0)
	stN.plot()

	fig = plt.figure()
	stE = op.read('*HHE*')
	stE.filter('bandpass',freqmin=1.0, freqmax=20.0)
	stE.plot()
	'''