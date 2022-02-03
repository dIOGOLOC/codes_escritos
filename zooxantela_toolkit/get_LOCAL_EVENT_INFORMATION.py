#!/usr/bin/python -u
'''
--------------------------------------------------------------------------------
            Function to collect information of a local events
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 02/2022


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
Given a starttime and endtime, this code will return a JSON file with a list of
events downloaded from Preliminary Seismic Boletim by Centro de Sismologia da
USP that are located inside a determined area given by a shapefile.

More information in:
http://moho.iag.usp.br/eq/latest

Inputs:
LOCAL_CSV_FILE


An example of LOCAL_CSV_FILE download in:
http://moho.iag.usp.br/eq/latest

evid;origin;longitude;latitude;depth;magnitude;magnitudet;region;author;mode
usp2020mums;2020-06-30T23:35:40.31Z;-66.721;-23.861;203.5;4.1;mb;"Jujuy Province, Argentina";jroberto;M
usp2020muij;2020-06-30T21:21:44.662Z;-69.448;-25.673;10.0;4.7;mb;"Northern Chile";cleusa;M
usp2020mtqe;2020-06-30T12:10:54.829Z;-66.554;-23.669;208.1;4.5;mb;"Jujuy Province, Argentina";cleusa;M
usp2020mtgy;2020-06-30T07:31:46.66Z;-66.742;-23.558;231.0;3.9;mb;"Jujuy Province, Argentina";cleusa;M
usp2020mskk;2020-06-29T20:07:40.43Z;-40.434;-3.896;0.0;1.9;MLv;"Groairas/CE";jroberto;M
usp2020mrlp;2020-06-29T07:36:22.4Z;-40.199;-3.409;0.0;1.4;mR;"Santana do Acaraú/CE";jroberto;M
usp2020mrkj;2020-06-29T06:59:14.529Z;-40.278;-3.350;0.0;1.3;mR;"Santana do Acaraú/CE";jroberto;M
usp2020mrhl;2020-06-29T05:29:47.46Z;-63.906;-18.941;78.3;3.9;mb;"Central Bolivia";jroberto;M
usp2020mrgr;2020-06-29T05:06:31.226Z;-63.787;-18.674;10.0;4.6;mb;"Central Bolivia";jroberto;M
usp2020mrbx;2020-06-29T02:42:26.453Z;-71.789;-15.723;124.6;4.6;mb;"Southern Peru";cleusa;M

Data description:
    evid: event name
    origin: year-month-dayThour:minute:second of the event
    longitude: latitude of the event
    latitude: longitude of the event
    depth: depth of the event
    magnitude: magnitude of the event

    see http://moho.iag.usp.br/eq/bulletin and http://moho.iag.usp.br/eq/latest
    for more informations
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
