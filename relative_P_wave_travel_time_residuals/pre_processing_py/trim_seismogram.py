'''
Function to trim raw sac files
'''


import os
from obspy.taup import TauPyModel


sta_dic = {'kstnm':['BPPF','STSN','GRJU','GENI','PRDT','TRZN','BUCO','STSR','BCDO'],
	'stla':[-6.2271,-6.0787,-5.8308,-5.4612,-5.3241,-5.1056,-5.1586,-5.2889,-5.4517],
	'stlo':[-47.2518,-46.5986,-46.0882,-45.5344,-44.3974,-42.6344,-43.2010,-43.8063,-45.0203]}

#Function to cut each file in each folder from every event:

def cut_data_by_event(kstnm,stla,stlo,ev_timeUTC,ev_julday,ev_year,ev_month,ev_day,ev_hour,ev_minute,ev_second,ev_microsecond,ev_lat,ev_long,ev_depth,ev_mag)
	if os.path.isdir(DIR_SAC+kstnm+'/'+ev_year+'/'+ev_month+'/'+ev_day) == True:
		#Calculating distance, azimuth and backazimuth
		dist,az,baz = op.geodetics.gps2dist_azimuth(ev_lat,ev_long,stla,stlo)
		gcarc = op.geodetics.kilometer2degrees(dist/1000)		
		if   EV_GCARC_MIN < gcarc < EV_GCARC_MAX and ev_mag < EV_MAGNITUDE_MB:
			print('Current Directory = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_month+'/'+ev_day)		
			os.chdir(DIR_SAC+kstnm+'/'+ev_year+'/'+ev_month+'/'+ev_day)
			event_directory = DIR_EVENT+kstnm+'/'+ev_year+'/'+ev_julday+'/'+ev_year+'.'+ev_month+'.'+ev_day+'.'+ev_hour+'.'+ev_minute+'.'+ev_second+'.'+ev_microsecond			
			os.makedirs(event_directory, exist_ok=True)	

			#Calculating ray parameter
			model = TauPyModel(model=TAUPY_MODEL)
			arrivals = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=gcarc, phase_list=["P"])
			arr = arrivals[0]

			#########################################################################################################################################################
			#Component E

			HHX = op.read('*'+KCMPNM_E+'*')
			HHX.trim((ev_timeUTC+arr.time)-CUT_BEFORE_P,(ev_timeUTC+arr.time)+CUT_AFTER_P)
			headerHHX = {
					'kstnm': kstnm, 'kcmpnm': KCMPNM_E, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
					'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(ev_hour),
					'nzjday': int(ev_julday), 'nzmin': int(ev_minute),'nzmsec': int(ev_microsecond), 'nzsec': int(ev_second), 
					'nzyear': int(ev_year),'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(dist/1000), 'gcarc': float(gcarc), 
					'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'t1':float(ev_timeUTC+arr.time)
					}

			sacHHX = op.io.sac.sactrace.SACTrace(data=HHX[0].data, **headerHHX)
			sacHHX.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+ev_year+'.'+ev_julday+'.'+ev_hour+'.'+ev_minute+'.'+ev_second+'.'+ev_microsecond+'.'+NAME_SUFFIX_E)

			#########################################################################################################################################################
			#Component N
			HHY = op.read('*'+KCMPNM_N+'*')
			HHY.trim((ev_timeUTC+arr.time)-CUT_BEFORE_P,(ev_timeUTC[i]+arr.time)+CUT_AFTER_P)
			headerHHY = {
					'kstnm': kstnm, 'kcmpnm': KCMPNM_N, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
					'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(ev_hour),
					'nzjday': int(ev_julday), 'nzmin': int(ev_minute), 'nzmsec': int(ev_microsecond),'nzsec': int(ev_second),
					'nzyear': int(ev_year),	'cmpaz': 0.0, 'cmpinc': 90.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 
					'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'t1':float(ev_timeUTC+arr.time)
					}

			sacHHY = op.io.sac.sactrace.SACTrace(data=HHY[0].data, **headerHHY)
			sacHHY.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+ev_year+'.'+ev_julday+'.'+ev_hour+'.'+ev_minute+'.'+ev_second+'.'+ev_microsecond+'.'+NAME_SUFFIX_N)


			#########################################################################################################################################################
			#Component Z

			HHZ = op.read('*'+KCMPNM_Z+'*')
			HHZ.trim((ev_timeUTC+arr.time)-CUT_BEFORE_P,(ev_timeUTC[i]+arr.time)+CUT_AFTER_P)
			headerHHZ = {
					'kstnm': kstnm, 'kcmpnm': KCMPNM_Z, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
					'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(ev_hour),
					'nzjday': int(ev_julday), 'nzmin': int(ev_minute), 'nzmsec': int(ev_microsecond),'nzsec': int(ev_second),
					'nzyear': int(ev_year),	'cmpaz': 0.0, 'cmpinc': 0.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 
					'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'t1':float(ev_timeUTC+arr.time)
					}

			sacHHZ = op.io.sac.sactrace.SACTrace(data=HHZ[0].data, **headerHHZ)
			sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+ev_year+'.'+ev_julday+'.'+ev_hour+'.'+ev_minute+'.'+ev_second+'.'+ev_microsecond+'.'+NAME_SUFFIX_Z)
