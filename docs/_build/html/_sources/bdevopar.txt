.. SpeX Prism Library Analysis Toolkit documentation master file, created by
   sphinx-quickstart on Sat Jul 11 20:07:28 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Brown Dwarf Evolutionary Model Package
===============================================================

.. toctree
   :maxdepth: 3

The SPLAT brown dwarf evolutionary model package serves routines that allow for the conversion between physical 
parameters (e.g., mass, age, metallicity) and observable parameters (temperature, luminosity, radius, surface gravity). Parameters are determined through logarithmic linear interpolation.

The evolutionary models currently available through this package are:

	- **burrows**: `Burrows et al. (2001) <http://adsabs.harvard.edu/abs/2001RvMP...73..719B>`_ for 1 Myr < age < 10 Gyr, 0.005 Msol < mass < 0.2 Msol, and solar metallicity
	- **baraffe**: `Baraffe et al. (2003) <http://adsabs.harvard.edu/abs/2003A&A...402..701B>`_ for 1 Myr < age < 10 Gyr, 0.005 Msol < mass < 0.1 Msol, and solar metallicity (COND dust prescription)
	- **saumon**: `Saumon et al. (2003) <http://adsabs.harvard.edu/abs/2003A&A...402..701B>`_ for 3 Myr < age < 10 Gyr, 0.002 Msol < mass < 0.085 Msol (although this varies, as the maximum temperature is 2500 K); in addition, there are options for
		- ``metallicity`` = ``solar``, ``+0.3``, ``-0.3``
		- ``cloud`` =  ``cloud-free``, ``hybrid``, ``f2``



.. automodule:: bdevopar
 :members:


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

