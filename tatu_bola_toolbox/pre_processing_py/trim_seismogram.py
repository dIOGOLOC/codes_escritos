'''
Function to trim raw sac files
'''


import os
import obspy as op
from obspy.taup import TauPyModel

#Function to cut each file in each folder from every event:

from parameters_py.config import (
					DIR_SAC,DIR_EVENT,NEIC_CSV_FILE,STA_CSV_FILE,OUTPUT_JSON_FILE_DIR,TAUPY_MODEL,
					EV_GCARC_MIN,EV_GCARC_MAX,EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,KCMPNM_N,
					KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,NAME_SUFFIX_Z,MP_PROCESSES,
					DETREND,DETREND_TYPE,TAPER,TAPER_TYPE,TAPER_MAX_PERCENTAGE,LOWPASS_FREQ,LOWPASS_CORNER,
					LOWPASS_ZEROPHASE,HIGHPASS_FREQ,HIGHPASS_CORNER,HIGHPASS_ZEROPHASE,RMEAN,RMEAN_TYPE,
					INTERPOLATE,SAMPLING_RATE,INTERPOLATE_TYPE
				   )

def cut_data_by_event(kstnm,stla,stlo,ev_timeUTC,ev_julday,ev_year,ev_month,ev_day,ev_hour,ev_minute,ev_second,ev_microsecond,ev_lat,ev_long,ev_depth,ev_mag):
	if os.path.isdir(DIR_SAC+kstnm+'/'+ev_year+'/'+ev_month+'/'+ev_day) == True:
		try:
			#Calculating distance, azimuth and backazimuth
			dist,az,baz = op.geodetics.gps2dist_azimuth(ev_lat,ev_long,stla,stlo)
			gcarc = op.geodetics.kilometer2degrees(dist/1000)
			if EV_GCARC_MIN <= gcarc <= EV_GCARC_MAX and EV_MAGNITUDE_MB <= ev_mag:
				print('Current Directory = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_month+'/'+ev_day)		
				os.chdir(DIR_SAC+kstnm+'/'+ev_year+'/'+ev_month+'/'+ev_day)

				#Calculating ray parameter
				model = TauPyModel(model=TAUPY_MODEL)
				arrivals = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=gcarc, phase_list=["P"])
				arr = arrivals[0]

				#Reference time
				starttime = (op.UTCDateTime(ev_timeUTC)+arr.time)-CUT_BEFORE_P
				endtime = (op.UTCDateTime(ev_timeUTC)+arr.time)+CUT_AFTER_P	

				#Creating Event Directory 
				event_directory = DIR_EVENT+kstnm+'/'+'{:04}'.format(starttime.year)+'/'+'{:03}'.format(starttime.julday)+'/'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:04}'.format(starttime.microsecond)			
				os.makedirs(event_directory, exist_ok=True)
				#########################################################################################################################################################
				#Component E
				HHX = op.read('*'+KCMPNM_E+'*')
				print(HHX)
				HHX.trim(starttime,endtime)
				
				if RMEAN == True:
					HHX.detrend(type=RMEAN_TYPE)  
				else:
					pass

				if DETREND == True:
					HHX.detrend(type=DETREND_TYPE) 
				else:
					pass

				if TAPER == True:
					HHX.taper(type=TAPER_TYPE,max_percentage=TAPER_MAX_PERCENTAGE) 
				else:
					pass

				if FILTER == True:
					HHX.filter('lowpass',freq=LOWPASS_FREQ,corners=LOWPASS_CORNER,zerophase=LOWPASS_ZEROPHASE) 
					HHX.filter('highpass',freq=HIGHPASS_FREQ,corners=HIGHPASS_CORNER,zerophase=HIGHPASS_ZEROPHASE)
				else:
					pass

				if INTERPOLATE == True:
					HHX.interpolate(sampling_rate=SAMPLING_RATE, method=INTERPOLATE_TYPE)
				else:
					pass
				
	
				headerHHX = {
						'kstnm': kstnm, 'kcmpnm': KCMPNM_E, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
						'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(starttime.hour),
						'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute),'nzmsec': int(starttime.microsecond), 'nzsec': int(starttime.second), 
						'nzyear': int(starttime.year),'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(dist/1000), 'gcarc': float(gcarc), 
						'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'o':float(CUT_BEFORE_P)
						}

				sacHHX = op.io.sac.sactrace.SACTrace(data=HHX[0].data, **headerHHX)
				sacHHX.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:04}'.format(starttime.microsecond)+'.'+NAME_SUFFIX_E)

				#########################################################################################################################################################
				#Component N
				HHY = op.read('*'+KCMPNM_N+'*')
				print(HHY)
				HHY.trim(starttime,endtime)

				if RMEAN == True:
					HHY.detrend(type=RMEAN_TYPE)  
				else:
					pass

				if DETREND == True:
					HHY.detrend(type=DETREND_TYPE) 
				else:
					pass

				if TAPER == True:
					HHY.taper(type=TAPER_TYPE,max_percentage=TAPER_MAX_PERCENTAGE) 
				else:
					pass

				if FILTER == True:
					HHY.filter('lowpass',freq=LOWPASS_FREQ,corners=LOWPASS_CORNER,zerophase=LOWPASS_ZEROPHASE) 
					HHY.filter('highpass',freq=HIGHPASS_FREQ,corners=HIGHPASS_CORNER,zerophase=HIGHPASS_ZEROPHASE)
				else:
					pass

				if INTERPOLATE == True:
					HHY.interpolate(sampling_rate=SAMPLING_RATE, method=INTERPOLATE_TYPE)
				else:
					pass


				headerHHY = {
						'kstnm': kstnm, 'kcmpnm': KCMPNM_N, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
						'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(starttime.hour),
						'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int(starttime.microsecond),'nzsec': int(starttime.second),
						'nzyear': int(starttime.year),	'cmpaz': 0.0, 'cmpinc': 90.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 
						'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'o':float(CUT_BEFORE_P)
						}

				sacHHY = op.io.sac.sactrace.SACTrace(data=HHY[0].data, **headerHHY)
				sacHHY.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:03}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:04}'.format(starttime.microsecond)+'.'+NAME_SUFFIX_N)


				#########################################################################################################################################################
				#Component Z

				HHZ = op.read('*'+KCMPNM_Z+'*')
				print(HHZ)
				HHZ.trim(starttime,endtime)

				if RMEAN == True:
					HHZ.detrend(type=RMEAN_TYPE) 
				else:
					pass
				
				if DETREND == True:
					HHZ.detrend(type=DETREND_TYPE) 
				else:
					pass

				if TAPER == True:
					HHZ.taper(type=TAPER_TYPE,max_percentage=TAPER_MAX_PERCENTAGE) 
				else:
					pass

				if FILTER == True:
					HHZ.filter('lowpass',freq=LOWPASS_FREQ,corners=LOWPASS_CORNER,zerophase=LOWPASS_ZEROPHASE) 
					HHZ.filter('highpass',freq=HIGHPASS_FREQ,corners=HIGHPASS_CORNER,zerophase=HIGHPASS_ZEROPHASE)
				else:
					pass

				if INTERPOLATE == True:
					HHZ.interpolate(sampling_rate=SAMPLING_RATE, method=INTERPOLATE_TYPE)
				else:
					pass

				headerHHZ = {
						'kstnm': kstnm, 'kcmpnm': KCMPNM_Z, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
						'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(ev_hour),
						'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int(starttime.microsecond),'nzsec': int(starttime.second),
						'nzyear': int(starttime.year),	'cmpaz': 0.0, 'cmpinc': 0.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 
						'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'o':float(CUT_BEFORE_P)
						}

				sacHHZ = op.io.sac.sactrace.SACTrace(data=HHZ[0].data, **headerHHZ)
				sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:04}'.format(starttime.microsecond)+'.'+NAME_SUFFIX_Z)
			return 'Data OK with this event = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_month+'/'+ev_day
		except: 
			return 'Data error with this event = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_month+'/'+ev_day
	return 'No data for this event = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_month+'/'+ev_day
