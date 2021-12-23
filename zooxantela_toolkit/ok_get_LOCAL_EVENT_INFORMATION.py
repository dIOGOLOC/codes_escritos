#!/usr/bin/python -u
'''
--------------------------------------------------------------------------------
                   Collecting information of a local events
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2021


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
Given a starttime and endtime, this code will return a JSON file with a list of
events downloaded from Preliminary Seismic Boletim by Centro de Sismologia da
USP that are located inside a determined area given by a shapefile.

More information in:
http://moho.iag.usp.br/eq/bulletin/

Dataset in:
http://moho.iag.usp.br/boletim/boletim_txt/boletim2000.txt
and
http://moho.iag.usp.br/boletim/boletim_txt/boletim2001p.txt

Inputs:
LOCAL_CSV_FILE

An example of LOCAL_CSV_FILE download in:
http://moho.iag.usp.br/boletim/boletim_txt/boletim2001p.txt is shown bellow:

 YEAR MMDD HHMMSS  LAT. LONG.  ERR ST DEPTH MAG. T CAT Io  AREA LOCALITY   COMMENTS
 2001 0107 035015  -17.70 -44.70  10 MG   0.  3.4  1  I   -       Pirapora     (UnB)
 2001 0123 092131  -05.28 -39.42  50 CE   0.  3.3  1  I   -       Quixeramobim (IAG,UFRN)
 2001 0221 152021  -11.28 -74.51  10 PU  33.  5.5  0  I   2       Central Peru (IRIS)AC-IIMM
 2001 0226 204200  -04.41 -38.29  05 CE   0.  3.7  1  I   -       Cascavel     (IAG,UnB)


Data description:
    YEAR: year of the event
    MMDD: month and day of the event
    HHMMSS: hour,minute second of the event
    LAT: latitude of the event
    LONG: longitude of the event
    DEPTH: depth of the event
    MAG: magnitude of the event
    ...


Outputs:
JSON file with event description:
    ev_timeUTC: event time in UTC (str)
    ev_year: year of the event
    ev_month: month of the event
    ev_day: day of the event
    ev_julday: julian day of the event
    ev_hour: hour of the event
    ev_minute: minute of the event
    ev_second: second of the event
    ev_microsecond: microsecond of the event
    evla: latitude of the event
    evlo: longitude of the event
    evdp: depth of the event
    mag: magnitude of the event

Examples of Usage (in command line):
   >> python get_LOCAL_EVENT_INFORMATION.py

'''
#Import scripts to get information:

from get_information_py import get_local_events_information
