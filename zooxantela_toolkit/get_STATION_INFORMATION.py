#!/usr/bin/python -u
"""
--------------------------------------------------------------------------------
                Collecting information of a selected group of stations
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2021


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
Given a CSV file in a specific format, this code will return a JSON file that
will be used as input in other programs.


Inputs:
An example of STA_CSV_FILE:
LOC;SENSOR;KNETWK;KSTNM;STLA;STLO;STEL;FDAY;EDAY
RIO;1456;ON;OBS55;----;----;----;2060-01-27;2030-05-14
SAOPAULO;1456;ON;OBS97;----;----;----;2089-12-28;1920-01-15;

Header explanation:
		LOC: Location of the station (str)
		SENSOR: Serial number of the sensor (int)
		KNETWK: Network name (str)
		KSTNM: Network name (str)
		STLA: Latitude of the station (float)
		STLO: Longitude of the station (float)
		STEL: Elevation/Depth of the station (float)
		FDAY: Deployment day - First day (year-month-day)
		EDAY: Recovery day - End day (year-month-day)


Outputs:
JSON file with same structure of the input file


Examples of Usage (in command line):
   >> python get_STATTION_INFORMATION.py

--------------------------------------------------------------------------------
"""

from get_information_py import get_station_information
