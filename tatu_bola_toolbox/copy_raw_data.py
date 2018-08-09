#!/usr/bin/python -u
"""
This script copy raw data.
"""
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool


# =====================================
# Importing copy data script_py 
# =====================================

from pre_processing_py.copy_convert_data import copy_convert_raw_data

# ==================================================
#  Importing some parameters from configuration file 
# ==================================================

from parameters_py.config import ( 
					DIR_RAW_DATA,DIR_SAC,KCMPNM_N,KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,
					NAME_SUFFIX_Z,MP_PROCESSES,OUTPUT_JSON_FILE_DIR
				   )



# ==================================
#  Function to call copy data script
# ==================================

def parallel_copy_data(folder_to_send,name_raw_data,name_file):
			
	result = copy_convert_raw_data(FOLDER_OUTPUT=folder_to_send,DATA_RAW=name_raw_data,FILE_NAME=name_file)

	return print(result)


# ============================================
#  Importing station dictionary from JSON file 
# ============================================

print('\n')
print('Looking for STATIONS data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'STA_dic.json'

sta_dic = json.load(open(filename_STA))

kstnm = sta_dic['kstnm']
stla = sta_dic['stla']
stlo = sta_dic['stlo']

for i in kstnm:
	print('Station = '+i)
print('\n')

# ==============================
#  Creating stations Input lists
# ==============================

datafile_lst1 = [] 
datafile_lst2 = [] 
datafile_lst3 = [] 
for root, dirs, files in os.walk(DIR_RAW_DATA):
	for datafile in files:
		if datafile.endswith('.1.m'):
			datafile_lst1.append(os.path.join(root, datafile))
		elif datafile.endswith('.2.m'):
			datafile_lst2.append(os.path.join(root, datafile))
		elif datafile.endswith('.3.m'):
			datafile_lst3.append(os.path.join(root, datafile))

datafile_lst1S = sorted(datafile_lst1)
datafile_lst2S = sorted(datafile_lst2)
datafile_lst3S = sorted(datafile_lst3)

print('========================= Number of files of each component ========================= ')
print('X component = '+str(len(datafile_lst1S)))
print('Y component = '+str(len(datafile_lst2S)))
print('Z component = '+str(len(datafile_lst3S)))
print('\n')


print('Creating stations input lists')
print('\n')

data_raw_list1 = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	data_raw_list1[i] = [l for l in datafile_lst1S if j in l]


data_raw_list2 = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	data_raw_list2[i] = [l for l in datafile_lst2S if j in l]

data_raw_list3 = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	data_raw_list3[i] = [l for l in datafile_lst3S if j in l]

input_list1 = [[]]*len(kstnm)
input_list2 = [[]]*len(kstnm)
input_list3 = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	date_file1 = [obspy.UTCDateTime(year=int(l.split('/')[-1].split('.')[0]), julday=int(l.split('/')[-1].split('.')[1])) for l in data_raw_list1[i]]
	file_name1 = [knetwk+'.'+j+'.'+KCMPNM_Z+'.'+'{:02}'.format(l.year)+'.'+'{:02}'.format(l.month)+'.'+'{:02}'.format(l.day)+'.sac' for k,l in enumerate(date_file1)]
	output_folder_sac1 = [DIR_SAC+j+'/'+'{:02}'.format(l.year)+'/'+'{:02}'.format(l.month)+'/'+'{:02}'.format(l.day)+'/' for k,l in enumerate(date_file1)]
	
	input_list1[i] = [output_folder_sac1,data_raw_list1[i],file_name1]
# ===========================
	date_file2 = [obspy.UTCDateTime(year=int(l.split('/')[-1].split('.')[0]), julday=int(l.split('/')[-1].split('.')[1])) for l in data_raw_list2[i]]
	file_name2 = [knetwk+'.'+j+'.'+KCMPNM_N+'.'+'{:02}'.format(l.year)+'.'+'{:02}'.format(l.month)+'.'+'{:02}'.format(l.day)+'.sac' for k,l in enumerate(date_file2)]
	output_folder_sac2 = [DIR_SAC+j+'/'+'{:02}'.format(l.year)+'/'+'{:02}'.format(l.month)+'/'+'{:02}'.format(l.day)+'/' for k,l in enumerate(date_file2)]
	
	input_list2[i] = [output_folder_sac2,data_raw_list2[i],file_name2]
# ===========================
	date_file3 = [obspy.UTCDateTime(year=int(l.split('/')[-1].split('.')[0]), julday=int(l.split('/')[-1].split('.')[1])) for l in data_raw_list3[i]]
	file_name3 = [knetwk+'.'+j+'.'+KCMPNM_E+'.'+'{:02}'.format(l.year)+'.'+'{:02}'.format(l.month)+'.'+'{:02}'.format(l.day)+'.sac' for k,l in enumerate(date_file3)]
	output_folder_sac3 = [DIR_SAC+j+'/'+'{:02}'.format(l.year)+'/'+'{:02}'.format(l.month)+'/'+'{:02}'.format(l.day)+'/' for k,l in enumerate(date_file3)]
	
	input_list3[i] = [output_folder_sac3,data_raw_list3[i],file_name3]


# ===========================
#  Copying Data Z Component
# ===========================

print('Copying, converting and saving Z Component for each station')
print('\n')

pool_1 = Pool(MP_PROCESSES)
for i,j in enumerate(input_list1):
	pool_1.starmap(parallel_copy_data, j)
pool_1.close()


# ===========================
#  Copying Data N Component
# ===========================

print('Copying, converting and saving N Component for each station')
print('\n')

pool_2 = Pool(MP_PROCESSES)
for i,j in enumerate(input_list2):
	pool_2.starmap(parallel_copy_data, j)
pool_2.close()

# ===========================
#  Copying Data E Component
# ===========================

print('Copying, converting and saving E Component for each station')
print('\n')

pool_3 = Pool(MP_PROCESSES)
for i,j in enumerate(input_list3):
	pool_3.starmap(parallel_copy_data, j)
pool_3.close()


print('Cutting finished!')
