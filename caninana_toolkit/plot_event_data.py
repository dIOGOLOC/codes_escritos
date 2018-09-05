#!/usr/bin/python -u
"""
This script reads sac data from a set of stations and trim the data.
"""
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool


# =====================================
# Importing trim data script_py 
# =====================================

from visual_py.event_plot import plot_event_data

# ==================================================
#  Importing some parameters from configuration file 
# ==================================================

from parameters_py.config import (
										DIR_DATA,TAUPY_MODEL,EV_GCARC_MIN,EV_GCARC_MAX,OUTPUT_JSON_FILE_DIR,MP_PROCESSES,
										OUTPUT_EV_DIR
				   )


# ============
#  event data
# ============

datafile_lst = [] 
for root, dirs, files in os.walk(OUTPUT_EV_DIR):
	for directories in dirs:
		datafile_name = os.path.join(root, directories)
		datafile_lst.append(datafile_name)
datafile_lstS = sorted(datafile_lst)

for i,j in enumerate(datafile_lstS):
	plot_event_data(j)

print('Plotting finished!')
