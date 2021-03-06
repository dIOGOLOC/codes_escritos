<p align="center">
  <img src="caninana_toolkit_logo.png">
</p>

This project is dedicated to provide a Python framework for analysing the quality of
seismological data based on [ObsPy](https://github.com/obspy/obspy/wiki). Some functions
has the [multiprocessing](https://docs.python.org/3/library/multiprocessing.html) implemented.

Version
---------
v0.1

Requirements
------------
The code is developped on Ubuntu with Python Python 3.7.

In addition to [Python 3.7](https://www.python.org/downloads/release/python-370/), you need
to install the following packages: 

- [numpy](http://www.numpy.org/)
- [matplotlib](http://matplotlib.org/)
- [ObsPy](https://github.com/obspy/obspy/wiki)
- [json](https://docs.python.org/3/library/json.html)
- [os](https://docs.python.org/3/library/os.html)
- [multiprocessing](https://docs.python.org/3/library/multiprocessing.html)

*I suggest to use the [Anaconda Cloud](https://anaconda.org/) to install your packages.


Brief explanation about the main code:
---------------------------------------

**First of all, you must to get information of your the stations:**

1) *python get_STATION_EVENT_INFORMATION.py*

**To check your dataset:**

2) *python get_plot_DATA_AVAILABILITY.py*

**After that, you need to create your XML file:**

3) *python create_XML_network.py*

**Then, you can plot events data to compare the response of your stations:**

4) *python cut_plot_EVENT_DATA.py*

**Finally, you estimating and plotting the probabilistic power spectral densities of your data in the following scenarios.**

If you want to estimate and plot a small percentage of your dataset to check a statistical noise level (suggestion: 20%):

5) *python estimate_plot_PPSD_PERCENTAGE.py*

If you want to estimate, save and plot the whole dataset to check the noise level:

6) *python estimate_plot_PPSD_TOTAL.py*

If you want to estimate, save and plot a specific time window of your dataset to check the noise level:

7) *python estimate_plot_PPSD_WINDOWED.py*


How to update
-------------
The code is still experimental so you should regularly check for (and pull) 
updates.

ToDo list
-------------
- Create a code to check the orientation of the sensor.

References
----------

- M. Beyreuther, R. Barsch, L. Krischer, T. Megies, Y. Behr and J. Wassermann (2010).
ObsPy: A Python Toolbox for Seismology.
*SRL*, **81(3)**, 530-533. DOI: 10.1785/gssrl.81.3.530


- L. Krischer, T. Megies, R. Barsch, M. Beyreuther, T. Lecocq, C. Caudron, J. Wassermann (2015).
ObsPy: a bridge for seismology into the scientific Python ecosystem.
*Computational Science & Discovery*, **8(1)**, 014003. DOI: 10.1088/1749-4699/8/1/014003

- McNamara, D. E. and Buland, R. P. (2004).
Ambient Noise Levels in the Continental United States.
*Bulletin of the Seismological Society of America*, **94 (4)**, 1517-1527.

- Peterson, J. (1993).
Observations and Modeling of Seismic Background Noise.
*U.S. Geological Survey open-file report*, **93-322**, Albuquerque, N.M.


Inspiration
----------
The code logo is based in Spilotes pullatus([Caninana](https://en.wikipedia.org/wiki/Spilotes_pullatus)), a snake from the northeast part of Brazil fauna. 
