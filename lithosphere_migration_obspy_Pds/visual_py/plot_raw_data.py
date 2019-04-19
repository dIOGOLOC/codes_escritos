'''
Script to copy raw data before stack
'''

#Useful modules

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import copy
import matplotlib
from matplotlib.cm import get_cmap
from mpl_toolkits.mplot3d import Axes3D
import shapefile
import scipy.io
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import json
import random
from matplotlib.colors import Normalize
from numpy import ma
from matplotlib import cbook
import collections


from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,GRID_PP_MULT,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG,DIST_GRID_PP,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,OUTPUT_DIR,DEPTH_TARGET
				   )

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Stations'+'/'

print('Looking for Receiver Functions data in JSON file in '+STA_DIR)
print('\n')
filename_STA = STA_DIR+'sta_dic.json'

sta_dic = json.load(open(filename_STA))

event_depth = sta_dic['event_depth']
event_lat = sta_dic['event_lat']
event_long = sta_dic['event_long']
event_dist = sta_dic['event_dist']
event_gcarc = sta_dic['event_gcarc']
event_sta = sta_dic['event_sta']
event_ray = sta_dic['event_ray']
sta_lat = sta_dic['sta_lat']
sta_long = sta_dic['sta_long']
sta_data = sta_dic['sta_data']
sta_time = sta_dic['sta_time']

RF_stack_data = [sum(i)/len(sta_data) for i in zip(*sta_data)]

plt.figure(figsize = (30,10))
for i, j in enumerate(sta_data): 
	plt.plot(sta_time[i],j,'gray',linewidth=0.5,label='RF data')
	plt.plot(sta_time[0],RF_stack_data,'k',linewidth=2,label='RF stack')
	plt.axhline(y=0, xmin=sta_time[0][0], xmax=sta_time[0][-1],c='k')
	plt.text(50,max(RF_stack_data),'N = '+str(len(sta_data)))
	plt.title('Receiver Functions')
	plt.xlim(0,25)
	#plt.ylim(-	0.1,0.1)
plt.show()

