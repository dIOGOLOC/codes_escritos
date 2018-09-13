'''
Function to trim raw sac files
'''


import os
import obspy as op
from obspy.taup import TauPyModel

from parameters_py.config import (
					DIR_SAC,DIR_EVENT,STA_CSV_FILE,OUTPUT_JSON_FILE_DIR,TAUPY_MODEL,
					EV_GCARC_MIN,EV_GCARC_MAX,EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,KCMPNM_N,
					KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,NAME_SUFFIX_Z,MP_PROCESSES
				   )
#Function to cut each file in each folder from every event:

def cut_data_by_event(kstnm,stla,stlo,ev_timeUTC,ev_julday,ev_year,ev_month,ev_day,ev_hour,ev_minute,ev_second,ev_microsecond,ev_lat,ev_long,ev_depth,ev_mag):
	if os.path.isdir(DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday) == True:
		try:
			#Calculating distance, azimuth and backazimuth
			dist,az,baz = op.geodetics.gps2dist_azimuth(ev_lat,ev_long,stla,stlo)
			gcarc = op.geodetics.kilometer2degrees(dist/1000)
			if EV_GCARC_MIN <= gcarc <= EV_GCARC_MAX and EV_MAGNITUDE_MB <= ev_mag:
				print('Current Directory = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday)		
				os.chdir(DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday)
				#Calculating ray parameter
				model = TauPyModel(model=TAUPY_MODEL)
				arrivals = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=gcarc, phase_list=["P"])
				arr = arrivals[0]
				#Reference time
				starttime = (op.UTCDateTime(ev_timeUTC)+arr.time)-CUT_BEFORE_P
				endtime = (op.UTCDateTime(ev_timeUTC)+arr.time)+CUT_AFTER_P	

				#Creating Event Directory 
				event_directory = DIR_EVENT+kstnm+'/'+'{:04}'.format(starttime.year)+'/'+'{:03}'.format(starttime.julday)+'/'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]			
				os.makedirs(event_directory, exist_ok=True)
				#########################################################################################################################################################
				#Component E
				HHX = op.read('*.'+NAME_SUFFIX_E)
				print(HHX)
				HHX.trim(starttime,endtime)
				
				headerHHX = {
						'kstnm': kstnm, 'kcmpnm': KCMPNM_E, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
						'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(starttime.hour),
						'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute),'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]), 'nzsec': int(starttime.second), 
						'nzyear': int(starttime.year),'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(dist/1000), 'gcarc': float(gcarc), 
						'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'o':float(CUT_BEFORE_P)
						}

				sacHHX = op.io.sac.sactrace.SACTrace(data=HHX[0].data, **headerHHX)
				sacHHX.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]+'.'+NAME_SUFFIX_E)
				
				#########################################################################################################################################################
				#Component N
				HHY = op.read('*.'+NAME_SUFFIX_N)
				print(HHY)
				HHY.trim(starttime,endtime)

				headerHHY = {
						'kstnm': kstnm, 'kcmpnm': KCMPNM_N, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
						'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(starttime.hour),
						'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),
						'nzyear': int(starttime.year),	'cmpaz': 0.0, 'cmpinc': 90.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 
						'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'o':float(CUT_BEFORE_P)
						}

				sacHHY = op.io.sac.sactrace.SACTrace(data=HHY[0].data, **headerHHY)
				sacHHY.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:03}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]+'.'+NAME_SUFFIX_N)


				#########################################################################################################################################################
				#Component Z

				HHZ = op.read('*.'+NAME_SUFFIX_Z)
				print(HHZ)
				HHZ.trim(starttime,endtime)

				headerHHZ = {
						'kstnm': kstnm, 'kcmpnm': KCMPNM_Z, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
						'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': float(ev_hour),
						'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),
						'nzyear': int(starttime.year),	'cmpaz': 0.0, 'cmpinc': 0.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 
						'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'o':float(CUT_BEFORE_P)
						}

				sacHHZ = op.io.sac.sactrace.SACTrace(data=HHZ[0].data, **headerHHZ)
				sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]+'.'+NAME_SUFFIX_Z)
				
			return 'Data OK with this event = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday
		except: 
			return 'Data error with this event = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday
		
	return 'No data for this event = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday
