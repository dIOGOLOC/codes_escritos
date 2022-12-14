'''
--------------------------------------------------------------------------------
  Function calculate local magnitude
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2022

'''

import os
import glob
import numpy as np
from tarfile import StreamError
import obspy as op
from obspy import Stream, Inventory
import pandas as pd
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.io.sac.sactrace import SACTrace

DIR_DATA = '/home/diogoloc/dados_posdoc/ON02_analysis/'

# --------- #
# Functions #
# --------- #


#Creating Event Directory
event_mseed = "/home/diogoloc/dados_posdoc/ON02_analysis/Events_Snuffler/2022.318.21.28.00.00/2022.318.21.28.00.00.mseed"

st = op.read(event_mseed)
st_sac = op.read("/home/diogoloc/dados_posdoc/ON02_analysis/Events_data/2022/318/2022.318.21.28.00.00/ON.ON02.2022.318.21.28.00.00.Z")
epi_dist = st_sac[0].stats.sac.dist
print(epi_dist)
on_st = st.select(station="ON02")
on_st.plot()

# Reading events files

paz_wa = {'sensitivity': 2800, 'zeros': [0j], 'gain': 1,
          'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]}

on_st.simulate(paz_remove=None, paz_simulate=paz_wa, water_level=10)
on_st.filter('bandpass', freqmin=15,freqmax=35, zerophase=True)
on_st.plot()


tr_n = on_st.select(component="N")[0]
ampl_n = max(abs(tr_n.data))

tr_e = on_st.select(component="E")[0]
ampl_e = max(abs(tr_e.data))


ampl = max(ampl_n, ampl_e)

a = 0.018
b = 2.17
ml = np.log10(ampl * 1000) + a * epi_dist + b

print(ml)
