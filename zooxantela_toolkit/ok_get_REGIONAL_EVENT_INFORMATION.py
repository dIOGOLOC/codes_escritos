#!/usr/bin/python -u
"""
--------------------------------------------------------------------------------
                     Collecting information of regional events
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2021


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
Given a starttime and endtime, this code will return a JSON file with a list of
events downloaded from Data Centers using OBSPY

More information in:
https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.get_events.html
and
https://docs.obspy.org/tutorial/code_snippets/retrieving_data_from_datacenters.html

Keep in mind that data centers and web services are constantly changing so this recommendation
might not be valid anymore at the time you read this.

Inputs:
INITIAL_DATE_EVENT: Initial date for looking for events
FINAL_DATE_EVENT: Final date for looking for events
EV_MAGNITUDE_MB: Event magnitude threshold


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
   >> python get_REGIONAL_EVENT_INFORMATION.py

"""

from get_information_py import get_events_information
