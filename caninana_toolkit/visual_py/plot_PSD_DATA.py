import os
import glob
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx

from parameters_py.config import (
					OUTPUT_FIGURE_DIR
				   )

# ===========================
# Function to plot PPSD DATA
# ===========================

def plot_PPSD_TOTAL_data(date_lst):
    os.chdir(date_lst)
    files = sorted(glob.glob('*.npz'))
    ppsd = PPSD.load_npz(files[0])
    [ppsd.add_npz(i) for i in files[1:]]
    os.makedirs(OUTPUT_FIGURE_DIR+'TOTAL/'+ppsd.station+'/',exist_ok=True)
    ppsd.plot(cmap=pqlx,filename=OUTPUT_FIGURE_DIR+'TOTAL/'+ppsd.station+'/'+ppsd.network+'.'+ppsd.station+'.'+ppsd.channel+'.'+str(ppsd.times_processed[0].year)+'.pdf')