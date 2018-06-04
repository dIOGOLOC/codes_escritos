#!/usr/bin/python -u
"""
[Advice: run this script using python with unbuffered output:
`python -u crosscorrelation.py`]
This script reads seismic waveform data from a set of stations, and
calculates the cross-correlations between all pairs of stations. The
data (in miniseed format) must be located in folder *MSEED_DIR*. The
stations information (coordinates, instrument response) can be read
from dataless seed files (if *USE_DATALESSPAZ* = True) located in
folder *DATALESS_DIR*, and/or stationXML files (if *USE_STATIONXML* =
True) located in folder *STATIONXML_DIR*. Note that two different
stations MUST HAVE DIFFERENT NAMES, even if they do not belong to
the same network. Also, one given station cannot have several
sets of coordinates: if so, it will be skipped.
In the current version of the program, miniseed files MUST be
organized inside their directory as:
<year>-<month>/<network>.<station>.<channel>.mseed, e.g.:
1988-10/BL.JFOB.BHZ.mseed
So, there is one sub-directory per month, and inside it, one miniseed
file  per month and per station.
The implemented algorithm follows the lines of Bensen et al.,
"Processing seismic ambient noise data to obtain reliable broad-band
surface wave dispersion measurements", Geophys. J. Int. (2007).
The procedure consists in stacking daily cross-correlations between
pairs of stations, from *FIRSTDAY* to *LASTDAY* and, in each given day,
rejecting stations whose data fill is < *MINFILL*. Define a subset of
stations to cross-correlate in *CROSSCORR_STATIONS_SUBSET* (or let it
empty to cross-correlate all stations). Define a list of locations to
skip in *CROSSCORR_SKIPLOCS*, if any. The cross-correlations are
calculated between -/+ *CROSSCORR_TMAX* seconds.
Several pre-processing steps are applied to the daily seismic waveform
data, before the daily cross-correlation is calculated and stacked:
(1) removal of the instrument response, the mean and the trend;
(2) band-pass filter between *PERIODMIN* and *PERIODMAX* sec
(3) down-sampling to sampling step = *PERIOD_RESAMPLE* sec
(4) time-normalization:
    - if *ONEBIT_NORM* = False, normalization of the signal by its
      (smoothed) absolute amplitude in the earthquake period band,
      defined as *PERIODMIN_EARTHQUAKE* - *PERIODMIN_EARTHQUAKE* sec.
      The smoothing window is *PERIODMAX_EARTHQUAKE* / 2;
    - if *ONEBIT_NORM* = False, one-bit normalization, wherein
      only the sign of the signal is kept (+1 or -1);
(5) spectral whitening of the Fourier amplitude spectrum: the Fourier
    amplitude spectrum of the signal is divided by a smoothed version
    of itself. The smoonthing window is *WINDOW_FREQ*.
Note that all the parameters mentioned above are defined in the
configuration file.
When all the cross-correlations are calculated, the script exports
several files in dir *CROSSCORR_DIR*, whose name (without extension)
is:
xcorr[_<stations of subset>]_<first year>-<last year>[_1bitnorm] ...
      _[datalesspaz][+][xmlresponse][_<suffix>]
where <suffix> is provided by the user. For example:
"xcorr_1996-2012_xmlresponse"
The files, depending on their extension, contain the following data:
- .pickle       = set of all cross-correlations (instance of
                  pscrosscorr.CrossCorrelationCollection) exported in binary
                  format with module pickle;
- .txt          = all cross-correlations exported in ascii format
                  (one column per pair);
- .stats.txt    = general information on cross-correlations in ascii
                  format: stations coordinates, number of days, inter-
                  station distance etc.
- .stations.txt = general information on the stations: coordinates,
                  nb of cross-correlations in which it appears, total
                  nb of days it has been cross-correlated etc.
- .png          = figure showing all the cross-correlations (normalized to
                  unity), stacked as a function of inter-station distance.
"""

from time_scripts import time_P_Pds
from variable_scripts import mgconfig
import os
import sys
import warnings

# ====================================================
# parsing configuration file to import some parameters
# ====================================================

from variable_scripts.mgconfig import (RF_DIR,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,DIST_T_DIR)

print("- dir of receiver function data: " + RF_DIR)
print("- dir of dist time: " + DIST_T_DIR)
print("- velocity model file: " + MODEL_FILE_NPZ)
print("- dir of the program: " + PROG_MIGRATION_DIR)
