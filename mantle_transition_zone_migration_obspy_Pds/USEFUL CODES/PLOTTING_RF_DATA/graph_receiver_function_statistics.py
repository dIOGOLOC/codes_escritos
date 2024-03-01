import numpy as np
import os
import obspy 
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.gridspec as gridspec
import scipy
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from obspy.taup import TauPyModel
from tqdm import tqdm
import pandas as pd
import glob
import multiprocessing
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# - - - - - - - - - - - - - - - - - - - - - -
# Finding RF files with and without PP filer
# - - - - - - - - - - - - - - - - - - - - - -

list_RF = sorted(glob.glob('/home/sysop/dados_posdoc/MTZ_2024/PRF_selected_YES_PP_FILTER_POST/*/*_P_R.sac'))

# -------------------------------------------------------------------------------------------------------------------------------------------------------

def storage_data(RF_file):
    #Reading RF data
    b = obspy.read(RF_file)

    #Allocating RF data to plot
    RF_dic = {
				'evla':b[0].stats.sac.evla,
				'evlo':b[0].stats.sac.evlo,
    			'evdp':b[0].stats.sac.evdp,
				'stla':b[0].stats.sac.stla,
				'stlo':b[0].stats.sac.stlo,
				'ray':b[0].stats.sac.user0,
				'gcarc':b[0].stats.sac.gcarc,
				'mag':b[0].stats.sac.mag,
				'year':b[0].stats.sac.nzyear            
			    }

    return RF_dic
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------

def calcular_em_paralelo(lst, num_processos):
    with multiprocessing.Pool(processes=num_processos) as pool:
        # Use tqdm para criar uma barra de progresso
        resultados = list(tqdm(pool.imap(storage_data, lst), total=len(lst), desc="Reading data"))

    return resultados

# -------------------------------------------------------------------------------------------------------------------------------------------------------

# Defina o número desejado de processos
num_processos = 20

# Calcular em paralelo com um número específico de processos e barra de progresso
resultado_final = calcular_em_paralelo(list_RF, num_processos)
RF_df = pd.DataFrame.from_dict(resultado_final)
RF_df_gcarc = RF_df.sort_values("gcarc")

###################################################################################################################

print('Plotting: Figure Distribution and usual statistics about events used in the study')
print('\n')


# Lat and Lon of large maps and small maps
# center (lat/lon) of the map

PROJECT_LAT = -5
PROJECT_LON = -46

#lower and upper corner (lat/lon) of the large map

LLCRNRLON_LARGE=-52
LLCRNRLAT_LARGE=-13
URCRNRLON_LARGE=-38
URCRNRLAT_LARGE=1

#lower and upper corner (lat/lon) of the small map

LLCRNRLON_SMALL=-49
LLCRNRLAT_SMALL=-12
URCRNRLON_SMALL=-40
URCRNRLAT_SMALL=-1

fig = plt.figure(figsize=(5, 10))
gs = gridspec.GridSpec(5,2)
gs.update(wspace=0.5,hspace=1)

ax1 = fig.add_subplot(gs[:3,:], projection=ccrs.Orthographic(central_longitude=PROJECT_LON, central_latitude=PROJECT_LAT,))

ax2 = fig.add_subplot(gs[3, 0])
ax3 = fig.add_subplot(gs[3, 1])
ax4 = fig.add_subplot(gs[4, 0])
ax5 = fig.add_subplot(gs[4, 1])

for lon, lat in zip(RF_df_gcarc['stlo'].values,RF_df_gcarc['stla'].values):
	ax1.plot(lon, lat, '^',markersize=6,markeredgecolor='k',markerfacecolor='w', transform=ccrs.Geodetic())

for lon, lat in zip(RF_df_gcarc['evlo'].values,RF_df_gcarc['evla'].values):
	ax1.plot(lon, lat, '*',markersize=6,markeredgecolor='k',markerfacecolor='dimgrey', transform=ccrs.Geodetic())

# draw coastlines and borders
ax1.add_feature(cfeature.LAND)
ax1.add_feature(cfeature.COASTLINE)
ax1.add_feature(cfeature.BORDERS, lw=0.5)
  
# draw meridians and parallels
gl = ax1.gridlines(color='k', linestyle=(0, (1, 1)), 
                  xlocs=range(0, 390, 20),
                  ylocs=[-80, -60, -30, 0, 30, 60, 80])

#############
ax2.hist(RF_df_gcarc['gcarc'].values,bins=20,orientation='vertical',color='k')
#ax2.set_yticklabels([])
ax2.set_title("Distance (Degrees)")
ax2.set_xlim(0,100)

#############
ax3.hist(RF_df_gcarc['ray'].values,bins=20,orientation='vertical',color='k')
#ax3.set_yticklabels([])
ax3.set_title("Ray Parameter (s/rad)")


#############
ax4.hist(RF_df_gcarc['mag'].values,bins=20,orientation='vertical',color='k') 
#ax4.set_yticklabels([])
ax4.set_title("Magnitude (mb)")
ax4.set_xlim(0,10)

#############
ax5.hist(RF_df_gcarc['year'].values,bins=20,orientation='vertical',color='k')
#ax5.set_yticklabels([])
ax5.set_title("Period (year)")
ax5.set_xlim(2012,2025)

plt.show()