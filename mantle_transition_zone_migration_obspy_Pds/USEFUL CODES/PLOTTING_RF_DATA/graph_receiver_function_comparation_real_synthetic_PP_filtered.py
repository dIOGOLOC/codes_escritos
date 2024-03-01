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

ev_listr_NO_FILTER = sorted(glob.glob('/home/sysop/dados_posdoc/MTZ_2024/PRF_selected_YES_PP_FILTER_POST/*/*_P_R.sac'))
ev_listr_YES_FILTER = sorted(glob.glob('/home/sysop/dados_posdoc/MTZ_2024/PRF_selected_YES_PP_FILTER_POST_SYNTHETIC/*/*_P_R.sac'))

# -------------------------------------------------------------------------------------------------------------------------------------------------------

def calc_P_PP_times(RF_file):
    #Reading RF data
    b = obspy.read(RF_file)

    #Calculating the P and PP times:    
    model = TauPyModel(model="iasp91")
    arrivalsP = model.get_travel_times(source_depth_in_km=b[0].stats.sac.evdp/1000, distance_in_degree=b[0].stats.sac.gcarc, phase_list=["P"])
    arrP = arrivalsP[0]
    
    arrivalsP410 = model.get_travel_times(source_depth_in_km=b[0].stats.sac.evdp/1000, distance_in_degree=b[0].stats.sac.gcarc, phase_list=["P410s"])
    arrP410 = arrivalsP410[0]

    arrivalsP660 = model.get_travel_times(source_depth_in_km=b[0].stats.sac.evdp/1000, distance_in_degree=b[0].stats.sac.gcarc, phase_list=["P660s"])
    arrP660 = arrivalsP660[0]

    arrivalsPP = model.get_travel_times(source_depth_in_km=b[0].stats.sac.evdp/1000, distance_in_degree=b[0].stats.sac.gcarc, phase_list=["PP"])
    arrPP = arrivalsPP[0]

    #Allocating RF data to plot
    RF_dic = {
				'data':b[0].data,
				'time':b[0].times(),
				'ray':b[0].stats.sac.user0,
				'gcarc':b[0].stats.sac.gcarc,
				'time_PP':arrPP.time - arrP.time,
				'time_P410':arrP410.time - arrP.time,
				'time_P660':arrP660.time - arrP.time                
			    }

    return RF_dic
    
# -------------------------------------------------------------------------------------------------------------------------------------------------------

def calcular_em_paralelo(lst, num_processos):
    with multiprocessing.Pool(processes=num_processos) as pool:
        # Use tqdm para criar uma barra de progresso
        resultados = list(tqdm(pool.imap(calc_P_PP_times, lst), total=len(lst), desc="Calculating"))

    return resultados

# -------------------------------------------------------------------------------------------------------------------------------------------------------

# Defina o número desejado de processos
num_processos = 20

# Calcular em paralelo com um número específico de processos e barra de progresso
resultado_final_YES_FILTER = calcular_em_paralelo(ev_listr_YES_FILTER, num_processos)
RF_df_YES_FILTER = pd.DataFrame.from_dict(resultado_final_YES_FILTER)
RF_df_ray_YES_FILTER = RF_df_YES_FILTER.sort_values("gcarc")

resultado_final_NO_FILTER = calcular_em_paralelo(ev_listr_NO_FILTER, num_processos)
RF_df_NO_FILTER = pd.DataFrame.from_dict(resultado_final_NO_FILTER)
RF_df_ray_NO_FILTER = RF_df_NO_FILTER.sort_values("gcarc")

RF_NO_FILTER = np.array(RF_df_ray_NO_FILTER['data'].tolist()).T
RF_YES_FILTER = np.array(RF_df_ray_YES_FILTER['data'].tolist()).T

# -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, (ax,ax1) = plt.subplots(1, 2, figsize=(10, 10))

majorLocatorY = MultipleLocator(20)
majorLocatorX = MultipleLocator(10)
minorLocatorY = MultipleLocator(10)
minorLocatorX = MultipleLocator(5)

#Sismograma sem filtro PP
v=0.01
  
im = ax.imshow(RF_NO_FILTER,extent=[RF_df_ray_NO_FILTER['gcarc'].tolist()[0],RF_df_ray_NO_FILTER['gcarc'].tolist()[-1],160,-10] ,interpolation='spline16', cmap=cm.viridis,
               origin='upper', aspect='auto',vmax=v, vmin=-v)

for i in range(len(RF_df_ray_NO_FILTER['time_PP'])):
    ax.plot(RF_df_ray_NO_FILTER['gcarc'][i],RF_df_ray_NO_FILTER['time_PP'][i]-10,'.k',markersize=1)
    ax.plot(RF_df_ray_NO_FILTER['gcarc'][i],RF_df_ray_NO_FILTER['time_PP'][i],'.r',markersize=1)
    ax.plot(RF_df_ray_NO_FILTER['gcarc'][i],RF_df_ray_NO_FILTER['time_PP'][i]+10,'.k',markersize=1)
    
ax.set_ylim(160,0)
ax.set_xlim(RF_df_ray_NO_FILTER['gcarc'].tolist()[0],RF_df_ray_NO_FILTER['gcarc'].tolist()[-1])
#ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
#ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('REAL DATA (n='+str(len(RF_df_ray_NO_FILTER['time_PP']))+')')
ax.grid(True)

#Sismograma com filtro PP

im = ax1.imshow(RF_YES_FILTER,extent=[RF_df_ray_YES_FILTER['gcarc'].tolist()[0],RF_df_ray_YES_FILTER['gcarc'].tolist()[-1],160,-10] ,interpolation='spline16', cmap=cm.viridis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

for i in range(len(RF_df_ray_YES_FILTER['time_PP'])):
    ax1.plot(RF_df_ray_YES_FILTER['gcarc'][i],RF_df_ray_YES_FILTER['time_PP'][i]-10,'.k',markersize=1)
    ax1.plot(RF_df_ray_YES_FILTER['gcarc'][i],RF_df_ray_YES_FILTER['time_PP'][i],'^r',markersize=1)
    ax1.plot(RF_df_ray_YES_FILTER['gcarc'][i],RF_df_ray_YES_FILTER['time_PP'][i]+10,'.k',markersize=1)
               
ax1.set_ylim(160,0)
ax1.set_xlim(RF_df_ray_YES_FILTER['gcarc'].tolist()[0],RF_df_ray_YES_FILTER['gcarc'].tolist()[-1])
ax1.xaxis.set_major_locator(majorLocatorX)
ax1.yaxis.set_major_locator(majorLocatorY)
ax1.xaxis.set_minor_locator(minorLocatorX)
ax1.yaxis.set_minor_locator(minorLocatorY)
ax1.set_ylabel('Time after P (s)')
ax1.set_xlabel('Slowness')
ax1.set_title('SYNTHETIC DATA (n='+str(len(RF_df_ray_YES_FILTER['time_PP']))+')')
ax1.grid(True)

plt.show()

# -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, (ax,ax1) = plt.subplots(1, 2, figsize=(10, 10))

majorLocatorY = MultipleLocator(20)
majorLocatorX = MultipleLocator(10)
minorLocatorY = MultipleLocator(10)
minorLocatorX = MultipleLocator(5)

#Sismograma sem filtro PP
v=0.01
  
im = ax.imshow(RF_NO_FILTER,extent=[RF_df_ray_NO_FILTER['gcarc'].tolist()[0],RF_df_ray_NO_FILTER['gcarc'].tolist()[-1],160,-10] ,interpolation='spline16', cmap=cm.viridis,
               origin='upper', aspect='auto',vmax=v, vmin=-v)

for i in range(len(RF_df_ray_NO_FILTER['time_PP'])):
    ax.plot(RF_df_ray_NO_FILTER['gcarc'][i],RF_df_ray_NO_FILTER['time_P410'][i],'.k',markersize=0.5,alpha=0.75)
    ax.plot(RF_df_ray_NO_FILTER['gcarc'][i],RF_df_ray_NO_FILTER['time_P660'][i],'.k',markersize=0.5,alpha=0.75)

print((RF_df_ray_NO_FILTER['time_P410'].values[-1]/10))

ax.text(RF_df_ray_NO_FILTER['gcarc'].values[-1]+(RF_df_ray_NO_FILTER['gcarc'].values[-1]/95),RF_df_ray_NO_FILTER['time_P410'].values[-1]+(RF_df_ray_NO_FILTER['time_P410'].values[-1]/95),'P410s')
ax.text(RF_df_ray_NO_FILTER['gcarc'].values[-1]+(RF_df_ray_NO_FILTER['gcarc'].values[-1]/95),RF_df_ray_NO_FILTER['time_P660'].values[-1]+(RF_df_ray_NO_FILTER['time_P410'].values[-1]/95),'P660s')
  
ax.set_ylim(160,0)
ax.set_xlim(RF_df_ray_NO_FILTER['gcarc'].tolist()[0],RF_df_ray_NO_FILTER['gcarc'].tolist()[-1])
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Distance (º)')
ax.set_title('REAL DATA (n='+str(len(RF_df_ray_NO_FILTER['time_PP']))+')')
ax.grid(True)

#Sismograma com filtro PP

im = ax1.imshow(RF_YES_FILTER,extent=[RF_df_ray_YES_FILTER['gcarc'].tolist()[0],RF_df_ray_YES_FILTER['gcarc'].tolist()[-1],160,-10] ,interpolation='spline16', cmap=cm.viridis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

for i in range(len(RF_df_ray_YES_FILTER['time_PP'])):
    ax1.plot(RF_df_ray_YES_FILTER['gcarc'][i],RF_df_ray_YES_FILTER['time_P410'][i],'.k',markersize=0.5,alpha=0.75)
    ax1.plot(RF_df_ray_YES_FILTER['gcarc'][i],RF_df_ray_YES_FILTER['time_P660'][i],'.k',markersize=0.5,alpha=0.75)

ax1.text(RF_df_ray_YES_FILTER['gcarc'].values[-1]+(RF_df_ray_YES_FILTER['gcarc'].values[-1]/95),RF_df_ray_YES_FILTER['time_P410'].values[-1]+(RF_df_ray_YES_FILTER['time_P410'].values[-1]/95),'P410s')
ax1.text(RF_df_ray_YES_FILTER['gcarc'].values[-1]+(RF_df_ray_YES_FILTER['gcarc'].values[-1]/95),RF_df_ray_YES_FILTER['time_P660'].values[-1]+(RF_df_ray_YES_FILTER['time_P410'].values[-1]/95),'P660s')
               
ax1.set_ylim(160,0)
ax1.set_xlim(RF_df_ray_YES_FILTER['gcarc'].tolist()[0],RF_df_ray_YES_FILTER['gcarc'].tolist()[-1])
ax1.xaxis.set_major_locator(majorLocatorX)
ax1.yaxis.set_major_locator(majorLocatorY)
ax1.xaxis.set_minor_locator(minorLocatorX)
ax1.yaxis.set_minor_locator(minorLocatorY)
ax1.set_ylabel('Time after P (s)')
ax1.set_xlabel('Distance (º)')
ax1.set_title('SYNTHETIC DATA (n='+str(len(RF_df_ray_YES_FILTER['time_PP']))+')')
ax1.grid(True)

plt.show()