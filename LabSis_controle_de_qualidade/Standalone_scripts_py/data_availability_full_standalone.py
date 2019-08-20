'''
Script to get information about the header of the raw data
(https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html#obspy.core.stream.read)
and plot a mosaic the the Data availability.
'''

import matplotlib.pyplot as plt
import obspy
import os
import glob
import json
import numpy as np
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib as mpl
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from obspy.signal import PPSD
import collections
from obspy.clients.arclink.client import Client
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


NETWORK = 'NB'

STATION = 'NBPB'

CHANNEL = 'HHZ'

INITIAL_DATE = '2019,1,1'

FINAL_DATE = '2019,7,1'

#-----------------------------------------------------------------------
DIR_DATA = '/home/sysop/seiscomp3/var/lib/archive/2019/NB/NBPB/HHZ.D/'


OUTPUT_FIGURE_DIR = '/home/sysop/Documentos/teste_lab_nbpb/FIGURES/'
#-----------------------------------------------------------------------


 

# ==============================
# Function to store data in dic
# ==============================

def nested_dict():
    return collections.defaultdict(nested_dict)

# ==============================

print('Processing folder: '+DIR_DATA)
os.chdir(DIR_DATA)
files = sorted(glob.glob('*'))

time_dic = nested_dict()

for i, j in enumerate(files):
    print("Extracting data from header "+str(i+1)+" of "+str(len(files)))
    print("File = "+DIR_DATA+j)

    #Lendo o cabeçalho da forma de onda diária do arclink
    st = obspy.read(DIR_DATA+j,headonly=True)


    #----------------------------
    #Começo do dicionário
    time_dic['network'] = st[0].stats.network
    time_dic['station'] = st[0].stats.station
    time_dic['channel'] = st[0].stats.channel
    
    time_lst = []
    for t,trace in enumerate(st):
        starttime = trace.stats.starttime
        endtime = trace.stats.endtime
        time_lst.append(np.arange(starttime,endtime,60))

    flat_time_lst = [item for sublist in time_lst for item in sublist]
    
    time_dic[files[i]]['DATETIME'] = str(st[0].stats.starttime.year)+','+str(st[0].stats.starttime.month)+','+str(st[0].stats.starttime.day)


    #Contador da lista de horas
    time_flat_time_lst = [[]]*24
    for g,h in enumerate(np.arange(24)):
        lst_time = []
        for x,c in enumerate(flat_time_lst):
            if c.hour == h:   
                lst_time.append(c.hour)
        time_flat_time_lst[g] = lst_time

    for q,w in enumerate(np.arange(24)):
        time_dic[files[i]]['HOUR'][str(w)] = len(time_flat_time_lst[int(w)])


    #Fim do dicionário
    #----------------------------

# ======================================================================
# Alocando os valores do contador de horas na variável time_plot_lst 
# ======================================================================

time_plot_lst = [[]]*len(files)
for i,j in enumerate(files):
    a = []
    for q,w in enumerate(np.arange(24)):
        a.append(time_dic[j]['HOUR'][str(w)])
    time_plot_lst[i] = a

# =========================================================================
# Calculando as datas entre o período inicial e final para plotar no eixo x 
# =========================================================================

datatime_initial = datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day)

datatime_final = datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day)

datetime_lista = np.arange(datatime_initial, datatime_final, datetime.timedelta(days=1)).astype(datetime.datetime)


xlim_initial = mdates.date2num(datatime_initial)
xlim_final = mdates.date2num(datatime_final)

#----------------------------
#Checando se as datas do banco de dados batem com o período escolhido (INITIAL_DATE e FINAL_DATE)

array_to_plot_by_xlim = []
for x,c in enumerate(datetime_lista):
    lista_temp = []
    for i,j in enumerate(files):
        if datetime.datetime(obspy.UTCDateTime(time_dic[j]['DATETIME']).year,obspy.UTCDateTime(time_dic[j]['DATETIME']).month,obspy.UTCDateTime(time_dic[j]['DATETIME']).day) == c:
            lista_temp.append(time_plot_lst[i])
    array_to_plot_by_xlim.append(lista_temp)


data_x_axis = []
for x,c in enumerate(array_to_plot_by_xlim):
    if c != []:
        data_x_axis.append(c[0][::-1])
    else:
        data_x_axis.append(np.zeros_like(np.arange(24)))

data_x_axis = np.array(data_x_axis).T


#----------------------------

# ====================================
# Function to plot DATA availability
# ====================================

#x axis parameters
days1 = DayLocator(interval=1)   # every day
days5 = DayLocator(interval=int(len(datetime_lista)*5/100))   # every day
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

#Matplotlib parameters
fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(20,10))
im = ax.imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.Greens,interpolation=None, vmin=0, vmax=60)
ax.set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
ax.yaxis.set_major_locator(MultipleLocator(4))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(days5)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days1)
ax.tick_params(which='minor', length=4)
ax.tick_params(which='major', length=10)
ax.set_ylim(0,24)
ax.set_ylabel('Hora do Dia')
ax.set_xlabel('Tempo')
ax.grid(b=True, which='major', color='k', linestyle='-')
ax.grid(b=True, which='minor', color='k', linestyle='-')
plt.setp(ax.xaxis.get_majorticklabels(), fontsize=10, rotation=30)
ax.set_title(time_dic['network']+'.'+time_dic['station']+'.'+time_dic['channel'], fontsize=20	,y=1.05)

#criando a localização da barra de cores:
axins = inset_axes(ax,
                   width="10%",  # width = 10% of parent_bbox width
                   height="5%",  # height : 50%
                   loc='upper left',
                   bbox_to_anchor=(0.85, 0.07, 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[0,30,60],label='Quantidade de arquivos por hora')
cbar.ax.set_xticklabels(['0%','50%','100%'])
os.makedirs(OUTPUT_FIGURE_DIR,exist_ok=True)
fig.savefig(OUTPUT_FIGURE_DIR+time_dic['network']+'_'+time_dic['station']+'_'+time_dic['channel']+'_'+'DATA.pdf',dpi=300)
plt.show()
