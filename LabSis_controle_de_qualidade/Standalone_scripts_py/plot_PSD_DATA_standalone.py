'''
Script to plot PPSD data based in 
https://docs.obspy.org/tutorial/code_snippets/probabilistic_power_spectral_density.html
'''

import os
import glob
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx
from obspy import UTCDateTime

#Diretório com os arquivos .npz
DIR_DATA = '/home/sysop/Documentos/teste_lab_nbpb/PPSD/2019/NBPB/HHE.PPSD/'

#Diretório para salvar as figuras
OUTPUT_FIGURE_DIR = '/home/sysop/Documentos/teste_lab_nbpb/TESTE_STANDALONE/FIGURES/'

#Initial date of the Data
INITIAL_DATE = '2019,1,1'

#Final date of the data
FINAL_DATE = '2019,7,1'

#Percentage fo the days to process and plot the PPSD?
DAY_PERCENTAGE = 2

#Restricts the data that is included in the stack by time of day and weekday. 
#Monday is 1, Sunday is 7, -1 for any day of week. 
#For example, using time_of_weekday=[(-1, 22, 24)] 
#only individual spectra that have a starttime in between 10pm and 12am are used in the stack for all days of week
#time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
TIME_OF_WEEKDAY_DAY = -1
TIME_OF_WEEKDAY_START_HOUR = 11
TIME_OF_WEEKDAY_FINAL_HOUR = 15

# ==================================
# Ploting TOTAL PPSD DATA
# ==================================
os.chdir(DIR_DATA)
files = sorted(glob.glob('*.npz'))
ppsd = PPSD.load_npz(files[0])

[ppsd.add_npz(i) for i in files[1:]]


ppsd.calculate_histogram(starttime=UTCDateTime(INITIAL_DATE),endtime=UTCDateTime(FINAL_DATE),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])    
folder_output = OUTPUT_FIGURE_DIR+'WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd.station+'/'
os.makedirs(folder_output,exist_ok=True)
ppsd.plot(cmap=pqlx,filename=folder_output+ppsd.network+'.'+ppsd.station+'.'+ppsd.channel+'.'+str(ppsd.times_processed[0].year)+'.pdf')
