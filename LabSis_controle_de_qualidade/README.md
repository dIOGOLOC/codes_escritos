<p align="center">
  <img src="labsis_logo.png">
</p>

This project is dedicated to provide a Python framework for analysing the quality of
seismological data archived in  based on [ObsPy](https://github.com/obspy/obspy/wiki).

Version
---------
v0.2

Requirements
------------
The code is developped on Ubuntu with Python 3.7. I think that is not a problem use in other systems.

In addition to [Python 3.7](https://docs.python.org/3/), you need
to install the following packages: 

- [numpy](http://www.numpy.org/)
- [matplotlib](http://matplotlib.org/)
- [ObsPy](https://github.com/obspy/obspy/wiki)
- [json](https://docs.python.org/3/library/json.html)
- [os](https://docs.python.org/3/library/os.html)

> I suggest to use the [Anaconda Cloud](https://anaconda.org/) to install your packages.


***Shortcuts to install the required packages on UBUNTU:***

> Download the '.sh' file in [Anaconda Cloud](https://anaconda.org/) and enter the following command to install Anaconda3:

```shell
$ bash ~/Downloads/Anaconda3-2019.07-Linux-x86_64.sh (File location)
```

> After installing Anaconda, enter the following command to install Obspy via Anaconda:

```shell
$ conda config --add channels conda-forge
$ conda install obspy
```

> The another required packages are pre-installed with [Anaconda Cloud](https://anaconda.org/).


Brief explanation about the main code:
---------------------------------------

- You should start reading the configuration file (config_file.cnf), which contains global parameters and detailed instructions.

- You should create your own 'STA_LAT_LON.txt' file, which is the main archive to create the XML file (see more about XML file [here](https://docs.obspy.org/tutorial/code_snippets/stationxml_file_from_scratch.html)).

- For creating your own 'STA_LAT_LON.txt' file with your stations, you should pay attention in the following parameters:

	- NAME (name of the station)
	- LAT (latitude of the station)
	- LON (longitude of the station)
	- ELEV (elevation of the station)
	- SENSOR_KEYS (NRL parameters of the station sensor)
	- DATALOGGER_KEYS (NRL parameters of the station data logger)
	- ACCER_KEYS (NRL parameters, if the station has another sensor, like a accelerograph) .

> The SENSOR_KEYS, DATALOGGER_KEYS and ACCER_KEYS are filled according to Nominal Response Library keys (see [here](http://docs.obspy.org/packages/obspy.clients.nrl.html) how to fill these parameters). 

```
Examples:
	sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])
```
> If the same station has two sensors, you must to fill SENSOR_KEYS and ACCER_KEYS, else let ACCER_KEYS empty.

> The program will process stations that were discriminated in the 'STA_LAT_LON.txt', so you can comment (use '%') the line with useless stations.

- After the program creates the XML file, we estimate the completeness of the dataset and the Probabilistic Power Spectral Densities of the data via Obspy (more information [here](https://docs.obspy.org/tutorial/code_snippets/probabilistic_power_spectral_density.html)).  

*DATASET VIA HD*

- First of all, you must to get information of your the stations:

```shell
$ python get_STATION_INFORMATION.py
```

- You need to create your XML file:

```shell
$ python create_XML_network.py
```

- To check and pre-process your dataset:

```shell
$ python get_plot_DATA_AVAILABILITY.py
```

- Estimating the probabilistic power spectral densities of your data.

```shell
$ python estimate_plot_PPSD_TOTAL.py
```

- Plotting the probabilistic power spectral densities of your data.**

```shell
$ python estimate_plot_PPSD_WINDOWED.py
```

---------------------------------------
---------------------------------------

*DATASET VIA CLIENT:*

- First of all, you must to get information of your the stations:

```shell
$ python get_STATION_INFORMATION.py
```

- You need to create your XML file:

```shell
$ python create_XML_network.py
```

- Estimating the probabilistic power spectral densities of your data.

```shell
$ python estimate_plot_PPSD_TOTAL_via_client.py
```

- Plotting the probabilistic power spectral densities of your data.

```shell
$ python estimate_plot_PPSD_WINDOWED.py
```

**If you want to check your dataset completeness:**

PLUS: 

```shell
$ python get_plot_DATA_AVAILABILITY_via_client.py
```


---------------------------------------
---------------------------------------

Stand-alone Option
-------------

*If you want a fast solution, you can use stand-alone  codes in the folder Standalone_scripts_py.*


How to update
-------------
The code is still experimental so you should regularly check for (and pull) updates.

ToDo list
-------------
- Daily PQLX frequency plot. 🔨🔨🔨

References
----------

- McNamara, D. E. and Buland, R. P. (2004).
Ambient Noise Levels in the Continental United States.
*Bulletin of the Seismological Society of America*, **94 (4)**, 1517-1527.

- Peterson, J. (1993).
Observations and Modeling of Seismic Background Noise.
*U.S. Geological Survey open-file report*, **93-322**, Albuquerque, N.M.


- M. Beyreuther, R. Barsch, L. Krischer, T. Megies, Y. Behr and J. Wassermann (2010).
ObsPy: A Python Toolbox for Seismology.
*SRL*, **81(3)**, 530-533. DOI: 10.1785/gssrl.81.3.530


- L. Krischer, T. Megies, R. Barsch, M. Beyreuther, T. Lecocq, C. Caudron, J. Wassermann (2015).
ObsPy: a bridge for seismology into the scientific Python ecosystem.
*Computational Science & Discovery*, **8(1)**, 014003. DOI: 10.1088/1749-4699/8/1/014003


Inspiration
----------
The code is for [Labsis](http://www.labsis.ufrn.br/). 
