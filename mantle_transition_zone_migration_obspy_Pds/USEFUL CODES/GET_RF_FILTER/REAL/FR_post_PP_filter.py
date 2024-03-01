import os
import glob
import obspy as op
import json
from tqdm import tqdm
from multiprocessing import Pool

from obspy import read,read_inventory, UTCDateTime, Stream
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.io.sac.sactrace import SACTrace
from obspy.signal.trigger import classic_sta_lta, trigger_onset, coincidence_trigger,recursive_sta_lta,plot_trigger

import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.transforms import offset_copy

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter,LatitudeLocator,LongitudeLocator
import shutil

# =====
# paths
# =====

# directory of raw files
INPUT_DIR = '/home/sysop/dados_posdoc/MTZ_2024/PRF_selected_NO_PP_FILTER/'

#Directory to save EVENT data
OUTPUT_EV_DIR = '/home/sysop/dados_posdoc/MTZ_2024/DATA_MTZ/PRF_selected_YES_PP_FILTER_POST/'

# =====
# event
# =====

#Taup_time model to calculate travel times
TAUPY_MODEL = 'iasp91'

#Minimum event distance
EV_GCARC_MIN = 30

# ==================================================
#  Importing Local Event dictionary from JSON file
# ==================================================

print('\n')
print('Looking for RF data')
print('\n')

RF = glob.glob(INPUT_DIR+'*/*P_R*')

# ====================
# saving PRFs selected
# ====================

for i in RF:
    a = read(i,headonly=True)
    if a[0].stats.sac.gcarc > 40:
        RF_sel_directory = OUTPUT_EV_DIR+'/'.join(i.split('/')[-2:-1])+'/'

        file_name = i.split('/')[-1].split('_')[0]

        os.makedirs(RF_sel_directory,exist_ok=True)
        print(RF_sel_directory+file_name+'_P_R.sac')
        shutil.copy2(i,RF_sel_directory+file_name+'_P_R.sac')
