# coding: utf-8


import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import obspy
import os
import glob
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees, degrees2kilometers
from geopy.distance import geodesic
import copy
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature
import shapefile
from scipy.stats import mode
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import random
from matplotlib.colors import Normalize
from numpy import ma
from matplotlib import cbook
import collections
from matplotlib.collections import PatchCollection
import verde as vd
from shapely.geometry import Polygon, MultiPoint, Point
from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle,Rectangle
import math
from tqdm import tqdm
import pyarrow.feather as feather
import pandas as pd
from sklearn import preprocessing

from parameters_py.mgconfig import (
					RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,SHAPEFILE_GRID,FILTER_BY_SHAPEFILE,
					NUMBER_PP_PER_BIN,RF_FREQUENCY,GRID_PP_MULT,CROSS_SECTION_AXIS,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,
					LLCRNRLON_SMALL,URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,
					PROJECT_LAT,PROJECT_LON,BOUNDARY_1_SHP,BOUNDARY_2_SHP,Ps_OR_Sp_PHASE,
					EXT_FIG,DPI_FIG,DIST_GRID_PP,NUMBER_STA_PER_BIN,OUTPUT_DIR,CONFIDENCE_BOUND,
					DEPTH_RANGE,BOOTSTRAP_INTERATOR,COLORMAP_STD,COLORMAP_VEL,DEPTH_MOHO,DEPTH_LAB,VMIN,VMAX
								   )


print('Starting Receiver Functions migration code to estimate the depths of the Mantle discontinuities')
print('\n')

print('Importing depths and times of '+Ps_OR_Sp_PHASE+' conversion  dataset')
print('\n')

print('Importing earth model from obspy.taup.TauPyModel')
print('Importing earth model from : '+MODEL_FILE_NPZ)
model_10_km = TauPyModel(model=MODEL_FILE_NPZ)

# =============================================
# Importing station dictionary from feather file 
# =============================================

print('\n')
print('Looking for receiver functions data in FEATHER file in '+OUTPUT_DIR)
print('\n')

filename_STA = glob.glob(OUTPUT_DIR+'*/*/sta_dic.feather')[0]
sta_dic = pd.read_feather(filename_STA)  

event_depth = sta_dic['event_depth'].tolist()
event_lat = sta_dic['event_lat'].tolist()
event_long = sta_dic['event_long'].tolist()
event_dist = sta_dic['event_dist'].tolist()
event_gcarc = sta_dic['event_gcarc'].tolist()
event_ray = sta_dic['event_ray'].tolist()
sta_lat = sta_dic['sta_lat'].tolist()
sta_long = sta_dic['sta_long'].tolist()
sta_data = sta_dic['sta_data'].tolist()
sta_time = sta_dic['sta_time'].tolist()
sta_baz = sta_dic['sta_baz'].tolist()
sta_network = sta_dic['sta_knetwk'][0]
sta_station = sta_dic['sta_kstnm'][0]

print('Plotting: Receiver Functions Figure')
print('\n')

ray_lst = event_ray

# Organazing according to ray parameter:
indices = np.argsort(ray_lst)

# Criação da figura e do primeiro eixo (eixo esquerdo)
fig, (ax1,ax2) = plt.subplots(2,1,height_ratios=[1,10],figsize = (10,15))

stack_amp = []
for i,j in enumerate(indices): 
    SRF_st_data = sta_data[j]
    SRF_st_time = sta_time[j]
    
    SRF_st_data = (preprocessing.normalize([SRF_st_data])[0])

    stack_amp.append(SRF_st_data)

    # ======================================================================================================================================================================
    
    factor_amp = 50

    ax2.plot(SRF_st_time,i/factor_amp+SRF_st_data,'k',linewidth=0.5)
    ax2.fill_between(SRF_st_time,i/factor_amp+SRF_st_data,i/factor_amp,where=(i/factor_amp+SRF_st_data)>=i/factor_amp, facecolor='red',alpha=0.6, interpolate=True)
    ax2.fill_between(SRF_st_time,i/factor_amp+SRF_st_data,i/factor_amp,where=(i/factor_amp+SRF_st_data)<=i/factor_amp, facecolor='blue',alpha=0.6, interpolate=True)
    ax2.text(30.5,i/factor_amp,"{0:.3f}".format(ray_lst[j]),fontsize=10)
    ax2.grid(which='major')
    ax2.xaxis.set_major_locator(MultipleLocator(5))
    ax2.xaxis.set_minor_locator(MultipleLocator(1))
    ax2.set_yticks([])
    ax2.set_xlim(0,30)

    # ======================================================================================================================================================================
stack_amp_mean = np.array(stack_amp).sum(axis=0)/len(stack_amp)

min_y = [min(a) for a in zip(*stack_amp)]
max_y = [max(a) for a in zip(*stack_amp)]

ax1.fill_between(SRF_st_time,min_y,max_y,alpha=0.01, facecolor='grey',interpolate=True)
ax1.plot(SRF_st_time,stack_amp_mean,'k',linewidth=0.5)
ax1.fill_between(SRF_st_time,stack_amp_mean,0,where=(stack_amp_mean)>=0, facecolor='red',alpha=0.6, interpolate=True)
ax1.fill_between(SRF_st_time,stack_amp_mean,0,where=(stack_amp_mean)<=0, facecolor='blue',alpha=0.6, interpolate=True)
ax1.grid(which='major')
ax1.xaxis.set_major_locator(MultipleLocator(5))
ax1.xaxis.set_minor_locator(MultipleLocator(1))
ax1.set_yticks([])
ax1.set_xlim(0,30)

OUTPUT_FIGURE_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_MOHO_'+str(DEPTH_MOHO)+'_DEPTH_LAB_'+str(DEPTH_LAB)+'/'+'Figures'+'/'
os.makedirs(OUTPUT_FIGURE_DIR,exist_ok=True)
fig.savefig(OUTPUT_FIGURE_DIR+'Receiver_Functions_'+Ps_OR_Sp_PHASE+'_'+sta_network+'_'+sta_station+'.png')

print('Creating the Earth layered model')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

print('Importing '+Ps_OR_Sp_PHASE+' piercing points')
print('\n')

PP_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_MOHO_'+str(DEPTH_MOHO)+'_DEPTH_LAB_'+str(DEPTH_LAB)+'/'+'Piercing_Points'+'/'

filename_1 = PP_DIR+'PP_dic.feather'

PP_dic = pd.read_feather(filename_1)

PP_time = PP_dic['time'].tolist()
PP_lat = PP_dic['lat'].tolist()
PP_lon = PP_dic['lon'].tolist()
PP_depth = PP_dic['depth'].tolist()


print('Piercing Points')
print('\n')

print('Migrating '+Ps_OR_Sp_PHASE+' dataset')
print('\n')

RF_amplitude_Pds_Sdp = [[]]*len(sta_time)
RF_amplitude_Pds_Sdp_dep = [[]]*len(sta_time)

for i,j in tqdm(enumerate(sta_time),desc='Migrating RF',total=len(sta_time)):
    
	PP_t_Pds = [round(l,1) for k,l in enumerate(PP_time[i])]
	RF_t_Pds = [round(l,1) for k,l in enumerate(sta_time[i])]

	RF_amplitude_Pds_Sdp[i] = [sta_data[i][np.argmin(np.abs(PP_t_Pds[k] - RF_t_Pds))] for k,l in enumerate(PP_t_Pds)]

print('\n')
print('Estimating Mean and Standard Deviation for RF (BOOTSTRAP)')
print('\n')

BOOTSTRAP_RF_DATA_Pds_Sdp = []
BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_dep = []
BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_amp = []
BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep = []
BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_amp = []

for _k in tqdm(range(BOOTSTRAP_INTERATOR),total=BOOTSTRAP_INTERATOR,desc='Bootstrap estimation'):
	RF_DATA_lst = np.arange(0,len(RF_amplitude_Pds_Sdp))
	
	new_RANDOM_RF_DATA = np.random.choice(RF_DATA_lst,size=len(RF_amplitude_Pds_Sdp),replace=True)
				
	new_RANDOM_RF_DATA_Pds_Sdp = [RF_amplitude_Pds_Sdp[_t_] for _t_ in new_RANDOM_RF_DATA]

	RF_STACKING_BOOTSTRAP_Pds_Sdp = (np.array(new_RANDOM_RF_DATA_Pds_Sdp).sum(axis=0)/len(new_RANDOM_RF_DATA_Pds_Sdp)).tolist()
	
	BOOTSTRAP_RF_DATA_Pds_Sdp.append(RF_STACKING_BOOTSTRAP_Pds_Sdp)

	# ---------------------------------------------------------------------------------------------

	######## Estimating MOHO depth ########

	lst_depth_amp_MOHO_Pds_Sdp = [RF_STACKING_BOOTSTRAP_Pds_Sdp[x] for x,c in enumerate(PP_depth[0]) if DEPTH_MOHO-DEPTH_RANGE <= c <= DEPTH_MOHO+DEPTH_RANGE]
	lst_depth_pp_MOHO_Pds_Sdp = [c for x,c in enumerate(PP_depth[0]) if DEPTH_MOHO-DEPTH_RANGE <= c <= DEPTH_MOHO+DEPTH_RANGE]
	lst_MOHO_depth_Pds_Sdp = lst_depth_pp_MOHO_Pds_Sdp[lst_depth_amp_MOHO_Pds_Sdp.index(max(lst_depth_amp_MOHO_Pds_Sdp))]
	lst_MOHO_amp_Pds_Sdp = lst_depth_amp_MOHO_Pds_Sdp[lst_depth_amp_MOHO_Pds_Sdp.index(max(lst_depth_amp_MOHO_Pds_Sdp))]
	BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_dep.append(lst_MOHO_depth_Pds_Sdp)
	BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_amp.append(lst_MOHO_amp_Pds_Sdp)

	######## Estimating LAB depth ########

	lst_depth_amp_LAB_Pds_Sdp = [RF_STACKING_BOOTSTRAP_Pds_Sdp[x] for x,c in enumerate(PP_depth[0]) if DEPTH_LAB-2*DEPTH_RANGE <= c <= DEPTH_LAB+2*DEPTH_RANGE]
	lst_depth_pp_LAB_Pds_Sdp = [c for x,c in enumerate(PP_depth[0]) if DEPTH_LAB-2*DEPTH_RANGE <= c <= DEPTH_LAB+2*DEPTH_RANGE]
	lst_LAB_depth_Pds_Sdp = lst_depth_pp_LAB_Pds_Sdp[lst_depth_amp_LAB_Pds_Sdp.index(min(lst_depth_amp_LAB_Pds_Sdp))]

	lst_LAB_amp_Pds_Sdp = lst_depth_amp_LAB_Pds_Sdp[lst_depth_amp_LAB_Pds_Sdp.index(min(lst_depth_amp_LAB_Pds_Sdp))]
	BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep.append(lst_LAB_depth_Pds_Sdp)
	BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_amp.append(lst_LAB_amp_Pds_Sdp)

# ----- # 

# Analysing BOOTSTRAP results:

# Waveforms
BOOTSTRAP_RF_DATA_STD = np.std(BOOTSTRAP_RF_DATA_Pds_Sdp,axis=0)*CONFIDENCE_BOUND
BOOTSTRAP_RF_DATA_MEAN = np.mean(BOOTSTRAP_RF_DATA_Pds_Sdp,axis=0)

# MOHO estimates
BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_dep_STD = np.std(BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_dep)*CONFIDENCE_BOUND
BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_dep_MEAN = np.mean(BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_dep)
BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_amp_STD = np.std(BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_amp)*CONFIDENCE_BOUND
BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_amp_MEAN = np.mean(BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_amp)

# LAB estimates
BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep_STD = np.std(BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep)*CONFIDENCE_BOUND
BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep_MEAN = np.mean(BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep)
BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_amp_STD = np.std(BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_amp)*CONFIDENCE_BOUND
BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_amp_MEAN = np.mean(BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_amp)

print('Plotting: Final Figure')
print('\n')

# Definindo o tamanho global da fonte
plt.rcParams.update({'font.size': 15})  # Define o tamanho da fonte para 15

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10, 10))  # Ajusta o tamanho para 8x6 polegadas

extent_Pds = [event_gcarc[0],event_gcarc[-1],MAX_DEPTH,MIN_DEPTH]

dados = np.array(BOOTSTRAP_RF_DATA_MEAN)
dados_std = np.array(BOOTSTRAP_RF_DATA_STD)
dados_time = PP_time[0]
dados_depth = PP_depth[0]

sta_data_mean = np.array(sta_data).sum(axis=0)/len(sta_data)

for i in range(len(sta_data)):
	ax1.plot(sta_data[i],sta_time[i],ls='-',lw=0.5,alpha=0.5,c='grey')
ax1.plot(sta_data_mean,sta_time[0],ls='-',lw=2,c='k')
ax1.set_ylim(30,0)
ax1.set_xlim(-0.5,0.5)
ax1.grid(which='both',ls='--',alpha=0.5)
ax1.yaxis.set_major_locator(MultipleLocator(5))
ax1.yaxis.set_minor_locator(MultipleLocator(1))

# Configurações adicionais
ax1.set_ylabel('Time (s)', fontweight='bold')
ax1.set_aspect('auto')
ax1.tick_params(axis='both', direction='inout', which='both', width=2)

for spine in ax1.spines.values():
    spine.set_linewidth(2) 

ax2.plot(dados.T,dados_depth.T,ls='-',lw=0.5,c='k')
ax2.fill_betweenx(dados_depth.T,0,dados.T,where=(dados.T)>=0, facecolor='red',alpha=0.9, interpolate=False)
ax2.fill_betweenx(dados_depth.T,0,dados.T,where=(dados.T)<=0, facecolor='blue',alpha=0.9, interpolate=False)

ax2.plot(dados.T+dados_std.T,dados_depth.T,ls='--',lw=0.5,c='grey')
ax2.plot(dados.T-dados_std.T,dados_depth.T,ls='--',lw=0.5,c='grey')
ax2.fill_betweenx(dados_depth.T,dados.T-dados_std.T,dados.T+dados_std.T,where=(dados.T+dados_std.T)>=dados.T-dados_std.T, facecolor='gray',alpha=0.2, interpolate=False)

ax2.errorbar(BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_amp_MEAN,BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_dep_MEAN,yerr=BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_dep_STD,c='k', fmt=".", capsize=3, capthick=1)
ax2.errorbar(BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_amp_MEAN,BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep_MEAN,yerr=BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep_STD,c='k', fmt=".", capsize=3, capthick=1)

xdata, ydata = BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_amp_MEAN,BOOTSTRAP_MOHO_RF_DATA_Pds_Sdp_dep_MEAN
xdisplay, ydisplay = ax2.transData.transform((xdata, ydata))

bbox = dict(boxstyle="round", fc="0.8")
arrowprops = dict(
    arrowstyle="->")

ax2.annotate(
    f'{ydata:.1f} $\pm$ {BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep_STD:.2f} km',
    (xdata, ydata),
	fontsize=12,
    xytext=(30, 0), textcoords='offset points',
    bbox=bbox, arrowprops=arrowprops)

xdata, ydata = BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_amp_MEAN,BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep_MEAN
xdisplay, ydisplay = ax2.transData.transform((xdata, ydata))

bbox = dict(boxstyle="round", fc="0.8")
arrowprops = dict(
    arrowstyle="->")

ax2.annotate(
    f'{ydata:.1f} $\pm$ {BOOTSTRAP_LAB_RF_DATA_Pds_Sdp_dep_STD:.2f} km',
    (xdata, ydata),
	fontsize=12,
    xytext=(30, 0), textcoords='offset points',
    bbox=bbox, arrowprops=arrowprops)



ax2.set_xlim(-0.5,0.5)
ax2.set_ylim(200,0)
ax2.grid(which='both',ls='--',alpha=0.5)
ax2.yaxis.set_major_locator(MultipleLocator(50))
ax2.yaxis.set_minor_locator(MultipleLocator(10))

# Configurações adicionais

if Ps_OR_Sp_PHASE == 'Ps':
	ax1.set_xlabel('Amplitude (% P-wave)', fontweight='bold')
	ax2.set_xlabel('Amplitude (% P-wave)', fontweight='bold')
	fig.suptitle(sta_network+'.'+sta_station+'\n P-Receiver Function', fontweight='bold')

if Ps_OR_Sp_PHASE == 'Sp':
	ax1.set_xlabel('Amplitude (% S-wave)', fontweight='bold')
	ax2.set_xlabel('Amplitude (% S-wave)', fontweight='bold')
	fig.suptitle(sta_network+'.'+sta_station+'\n S-Receiver Function', fontweight='bold')

ax2.set_ylabel('Depth (km)', fontweight='bold')
ax2.set_aspect('auto')
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
ax2.tick_params(axis='both', direction='inout', which='both', width=2)

for spine in ax2.spines.values():
    spine.set_linewidth(2) 

fig.savefig(OUTPUT_FIGURE_DIR+'Receiver_Functions_migration_'+Ps_OR_Sp_PHASE+'_'+sta_network+'_'+sta_station+'.png')

# --------