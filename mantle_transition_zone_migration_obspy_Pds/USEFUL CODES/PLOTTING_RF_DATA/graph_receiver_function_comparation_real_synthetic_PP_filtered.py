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

ev_listr_real = []
for root, dirs, files in os.walk('/home/sysop/dados_posdoc/MTZ_2024/PRF_SEISPY_DATA_NO_FILTER_PP/'):
    for datafile in files:
        if datafile.endswith('_P_R.sac'):
            ev_listr_real.append(os.path.join(root, datafile))

# - - - - - - - - - - - - - - - - - - - - - -
# Finding RF files with and without PP filer
# - - - - - - - - - - - - - - - - - - - - - -

ev_listr_NO_FILTER = sorted(glob.glob('/home/sysop/dados_posdoc/MTZ_2024/PRF_selected_NO_PP_FILTER/*/*_P_R.sac'))[:10]
ev_listr_YES_FILTER = sorted(glob.glob('/home/sysop/dados_posdoc/MTZ_2024/PRF_selected_YES_PP_FILTER/*/*_P_R.sac'))[:10]

# -------------------------------------------------------------------------------------------------------------------------------------------------------

def calc_P_PP_times(RF_file):
    #Reading RF data
    b = obspy.read(RF_file)

    #Calculating the P and PP times:    
    model = TauPyModel(model="iasp91")
    arrivalsP = model.get_travel_times(source_depth_in_km=b[0].stats.sac.evdp, distance_in_degree=b[0].stats.sac.gcarc, phase_list=["P"])
    arrP = arrivalsP[0]
    
    arrivalsP410 = model.get_travel_times(source_depth_in_km=b[0].stats.sac.evdp, distance_in_degree=b[0].stats.sac.gcarc, phase_list=["P410s"])
    arrP410 = arrivalsP410[0]

    arrivalsP660 = model.get_travel_times(source_depth_in_km=b[0].stats.sac.evdp, distance_in_degree=b[0].stats.sac.gcarc, phase_list=["P660s"])
    arrP660 = arrivalsP660[0]

    arrivalsPP = model.get_travel_times(source_depth_in_km=b[0].stats.sac.evdp, distance_in_degree=b[0].stats.sac.gcarc, phase_list=["PP"])
    arrPP = arrivalsPP[0]

    #Allocating RF data to plot
    RF_dic = {
				'data':b[0].data,
				'time':b[0].times(),
				'ray':b[0].stats.sac.user0,
				'gcarc':b[0].stats.sac.gcarc,
				'time_PP':arrPP.time - arrP.time + 10,
				'time_P410':arrP410.time - arrP.time + 10,
				'time_P660':arrP660.time - arrP.time + 10                
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
RF_df_ray_YES_FILTER = RF_df_YES_FILTER.sort_values("ray")

resultado_final_NO_FILTER = calcular_em_paralelo(ev_listr_NO_FILTER, num_processos)
RF_df_NO_FILTER = pd.DataFrame.from_dict(resultado_final_NO_FILTER)
RF_df_ray_NO_FILTER = RF_df_NO_FILTER.sort_values("ray")

# -------------------------------------------------------------------------------------------------------------------------------------------------------
fig, (ax,ax1) = plt.subplots(1, 2, figsize=(10, 10))

majorLocatorX = MultipleLocator(300)
majorLocatorY = MultipleLocator(0.01)
minorLocatorY = MultipleLocator(20)
minorLocatorX = MultipleLocator(0.001)


#Sismograma sem filtro PP
v=0.01
Z_real = RF_df_ray_NO_FILTER['data'].values
for i in Z_real:
    print(type(i))
#im = ax.imshow(Z_real,extent=[0,160,RF_df_ray_NO_FILTER['ray'].tolist()[0],RF_df_ray_NO_FILTER['ray'].tolist()[-1]] ,interpolation='bicubic', cmap=cm.cividis,origin='upper', aspect='auto',vmax=v, vmin=-v)
#im = ax.imshow(Z_real,interpolation='bicubic', cmap=cm.cividis,origin='upper', aspect='auto',vmax=v, vmin=-v)
#plt.show()

'''
for i,j in enumerate(RF_df_ray_NO_FILTER['time_PP']):
    ax.plot(i,(j-30)*10,'.k',markersize=2)
    ax.plot(i,j*10,'^r',markersize=2)
    ax.plot(i,(j+30)*10,'.k',markersize=2)
    
ax.set_ylim(0,160)
ax.set_xlim(RF_df_ray_NO_FILTER['ray'].tolist()[0],RF_df_ray_NO_FILTER['ray'].tolist()[-1])
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('YES PP PHASE')

ax.grid(True)

#Sismograma com filtro PP
Z_synth = RF_df_ray_YES_FILTER['data'].values

im = ax.imshow(Z_synth,extent=[0,160,RF_df_ray_YES_FILTER['ray'].tolist()[0],RF_df_ray_YES_FILTER['ray'].tolist()[-1]] ,interpolation='bicubic', cmap=cm.cividis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

for i,j in enumerate(RF_df_ray_NO_FILTER['ray']):
    ax1.plot(i,(j-30)*10,'.k',markersize=2)
    ax1.plot(i,(j)*10,'^r',markersize=2)
    ax1.plot(i,(j+30)*10,'.k',markersize=2)
               
ax1.set_ylim(0,160)
ax1.set_xlim(RF_df_ray_NO_FILTER['ray'].tolist()[0],RF_df_ray_NO_FILTER['ray'].tolist()[-1])
ax1.set_ylabel('Time after P (s)')
ax1.set_xlabel('Slowness')
ax1.set_title('NO PP PHASE')

ax1.xaxis.set_major_locator(majorLocatorX)
ax1.yaxis.set_major_locator(majorLocatorY)
ax1.xaxis.set_minor_locator(minorLocatorX)
ax1.yaxis.set_minor_locator(minorLocatorY)


ax1.grid(True)

plt.show()


# -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, (ax,ax1) = plt.subplots(1, 2, figsize=(10, 10))

#Sismograma sem filtro PP
v=0.01
im = ax.imshow(Z_real.T,extent=[0,160,RP_real[0],RP_real[-1]], interpolation='bicubic', cmap=cm.cividis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

for i,j in enumerate(time_P410_wave_corrected_real):
    ax.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)
    
for i,j in enumerate(time_P660_wave_corrected_real):
    ax.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)
    
ax.set_ylim(0,160)
ax.set_xlim(RP_real[0],RP_real[-1])
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('YES PP PHASE')

ax.grid(True)

#Sismograma com filtro PP

im = ax1.imshow(Z_synth.T,extent=[0,160,RP_synth[0],RP_synth[-1]],interpolation='bicubic', cmap=cm.viridis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

for i,j in enumerate(time_P410_wave_corrected_synth):
    ax1.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)
    
for i,j in enumerate(time_P660_wave_corrected_synth):
    ax1.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)
               
ax.set_ylim(0,160)
ax.set_xlim(RP_synth[0],RP_synth[-1])
ax1.set_ylabel('Time after P (s)')
ax1.set_xlabel('Slowness')
ax1.set_title('NO PP PHASE')

ax1.xaxis.set_major_locator(majorLocatorX)
ax1.yaxis.set_major_locator(majorLocatorY)
ax1.xaxis.set_minor_locator(minorLocatorX)
ax1.yaxis.set_minor_locator(minorLocatorY)


ax1.grid(True)
plt.show()


# -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, (ax,ax1) = plt.subplots(1, 2, figsize=(10, 10))

#Sismograma sem filtro PP
v=0.01
im = ax.imshow(Z_real.T,extent=[0,160,RP_real[0],RP_real[-1]], interpolation='bicubic', cmap=cm.cividis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

ax.set_ylim(0,160)
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('YES PP PHASE')

ax.grid(True)

#Sismograma com filtro PP

im = ax1.imshow(Z_synth.T,extent=[0,160,RP_synth[0],RP_synth[-1]],interpolation='bicubic', cmap=cm.viridis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

ax1.set_ylim(0,160)
ax1.set_ylabel('Time after P (s)')
ax1.set_xlabel('Slowness')
ax1.set_title('NO PP PHASE')

ax1.xaxis.set_major_locator(majorLocatorX)
ax1.yaxis.set_major_locator(majorLocatorY)
ax1.xaxis.set_minor_locator(minorLocatorX)
ax1.yaxis.set_minor_locator(minorLocatorY)


ax1.grid(True)

# -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, (ax,ax1) = plt.subplots(1, 2, figsize=(10, 10))

#Sismograma sem filtro PP
im = ax.imshow(Z_real.T,extent=[0,160,RP_real[0],RP_real[-1]], interpolation='bicubic', cmap=cm.cividis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

for i,j in enumerate(time_P410_wave_corrected_real):
    ax.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)

    
for i,j in enumerate(time_P660_wave_corrected_real):
    ax.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)
    
    
ax.set_ylim(0,160)
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('Real Data')

ax.grid(True)
ax.set_yticklabels(["{0:.0f}".format(time_real[i]) for i in np.arange(-100,len(time_real),100)])
ax.set_xticklabels(["{0:.1f}".format(RP_real[i]*100) for i in np.arange(0,len(RP_real),100)])
#ax.set_xticklabels(["{0:.1f}".format(GCARC_real[i]) for i in np.arange(0,len(GCARC_real),100)])

#Sismograma com filtro PP

im = ax1.imshow(Z_synth.T,extent=[0,160,RP_synth[0],RP_synth[-1]],interpolation='bicubic', cmap=cm.viridis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

for i,j in enumerate(time_P410_wave_corrected_synth):
    ax1.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)
    
for i,j in enumerate(time_P660_wave_corrected_synth):
    ax1.plot(i,j*10,'.k',markersize=0.5,alpha=0.75)

               
ax1.set_ylim(0,160)
ax1.set_ylabel('Time after P (s)')
ax1.set_xlabel('Slowness')
ax1.set_title('Synthetic Data')

ax1.xaxis.set_major_locator(majorLocatorX)
ax1.yaxis.set_major_locator(majorLocatorY)
ax1.xaxis.set_minor_locator(minorLocatorX)
ax1.yaxis.set_minor_locator(minorLocatorY)


ax1.grid(True)

# -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, (ax,ax1) = plt.subplots(1, 2, figsize=(10, 10))

#Sismograma sem filtro PP
v=0.008
im = ax.imshow(Z_real.T, interpolation='none', cmap=cm.viridis,
                origin='upper', aspect='auto',
                vmax=v, vmin=-v)
                #vmax=abs(Z.min()), vmin=Z.min())


    
ax.set_ylim(0,160)
ax.xaxis.set_major_locator(majorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_minor_locator(minorLocatorY)
ax.set_ylabel('Time after P (s)')
ax.set_xlabel('Slowness')
ax.set_title('Real Data')

ax.grid(True)

#Sismograma com filtro PP

im = ax1.imshow(Z_synth.T,extent=[0,160,RP_synth[0],RP_synth[-1]],interpolation='bicubic', cmap=cm.viridis,
                origin='upper', aspect='auto',vmax=v, vmin=-v)

ax1.set_ylim(0,160)
ax1.set_ylabel('Time after P (s)')
ax1.set_xlabel('Slowness')
ax1.set_title('Synthetic Data')

ax1.xaxis.set_major_locator(majorLocatorX)
ax1.yaxis.set_major_locator(majorLocatorY)
ax1.xaxis.set_minor_locator(minorLocatorX)
ax1.yaxis.set_minor_locator(minorLocatorY)


ax1.grid(True)

plt.show()
'''