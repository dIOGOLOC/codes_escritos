<p align="center">
  <img src="labsis_logo.png">
</p>

This project is dedicated to provide a Python framework for analysing the quality of
seismological data based on [ObsPy](https://github.com/obspy/obspy/wiki).

Version
---------
v0.2

Requirements
------------
The code is developped on Ubuntu with Python Python 3.7.

In addition to [Python 3.7](https://docs.python.org/3/), you need
to install the following packages: 

- [numpy](http://www.numpy.org/)
- [matplotlib](http://matplotlib.org/)
- [ObsPy](https://github.com/obspy/obspy/wiki)
- [json](https://docs.python.org/3/library/json.html)
- [os](https://docs.python.org/3/library/os.html)

*I suggest to use the [Anaconda Cloud](https://anaconda.org/) to install your packages.


Brief explanation about the main code:
---------------------------------------

*DATASET VIA HD*

**First of all, you must to get information of your the stations:**

1) *python get_STATION_INFORMATION.py*

**You need to create your XML file:**

2) *python create_XML_network.py*

**To check and pre-process your dataset:**

3) *python get_plot_DATA_AVAILABILITY.py*

**Estimating the probabilistic power spectral densities of your data.**

4) *python estimate_plot_PPSD_TOTAL.py*

**Plotting the probabilistic power spectral densities of your data.**

5) *python estimate_plot_PPSD_WINDOWED.py*

---------------------------------------
---------------------------------------

*DATASET VIA CLIENT:*

**First of all, you must to get information of your the stations:**

1) *python get_STATION_INFORMATION.py*

**You need to create your XML file:**

2) *python create_XML_network.py*

**Estimating the probabilistic power spectral densities of your data.**

3) *python estimate_plot_PPSD_TOTAL_via_client.py*

**Plotting the probabilistic power spectral densities of your data.**

4) *python estimate_plot_PPSD_WINDOWED.py*

**If you want to check your dataset:**

PLUS: *python get_plot_DATA_AVAILABILITY_via_client.py*

**VIA CLIENT**



How to update
-------------
The code is still experimental so you should regularly check for (and pull) 
updates.

ToDo list
-------------
- ...

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
The code is for ([Labsis](http://www.labsis.ufrn.br/)). 
