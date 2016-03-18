.. SpeX Prism Library Analysis Toolkit documentation master file, created by
   sphinx-quickstart on Sat Jul 11 20:07:28 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



Main SPLAT module
===============================================================

.. toctree
   :maxdepth: 3

The primary SPLAT module contains the definition of the core Spectrum object and associated function, 
standard classification and spectral analysis routines, and helper functions.


The SPLAT object
-----------------

Spectal data are manipulated through a SPLAT class object, which contains the relevant data (wavelength,
flux, uncertainty) and additional information on the source and/or observation.  A Spectrum object is 
the output of the various database access routines:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]	# note that getSpectrum returns a list
>>> sp = splat.getStandard('M0')[0]
>>> sp = splat.getModel(teff = 700, logg=4.5)
>>> splist = splat.getSpectrum(spt=['M7','L5'],jmag=[14.,99.])

You can also read in your own spectrum by passing a filename

>>> sp = splat.Spectrum(filename='PATH_TO/myspectrum.fits')

Note that this file must conform to the standard of the SPL data: the first column is
wavelength in microns, second column flux in f_lambda units, third column (optional) is 
flux uncertainty. The file can be a fits or ascii file.

There are many aspects of the Spectrum object that can be referenced or set, all of which are 
described in the API entry. Some primary examples:

``sp.wave``, ``sp.flux``, ``sp.noise``
	Wavelengths, flux density (per wavelength) values, and flux uncertainty of the spectrum (uncertainty = NaN for models)
``sp.wunit``, ``sp.funit``
	Astropy units for wavelength and flux, by default microns and erg/cm\ :sup: `2`\/s/micron
``sp.nu``, ``sp.fnu``, ``sp.fnu_unit``
	Frequencies, flux density (per frequency), and fnu units, by default Jy
``sp.snr``
	Median estimate of spectrum signal-to-noise
``sp.header``
	Fits header (dictionary) from original file, if present
``sp.teff``, ``sp.logg``, ``sp.z``, ``sp.fsed``, ``sp.cld``, ``sp.kzz``, 
	If Spectrum is a model, these values give the effective temperature (in K), log surface gravity (cm/s\ :sup:`2`),
	log metallicity (relative to Sun), sedimentation efficient f\ :sub: sed\ (for Saumon & Marley 2012 models),
	cloud coverage fraction (for Morley models) and non-equilibrium chemistry diffusion constant

The class ``info()`` command produces a summary of the Spectrum object's primary information:

>>>sp.info()


Built-in Commands
^^^^^^^^^^^^^^^^^

* To scale a spectrum



* To flux calibrate the spectrum, use the object's built in ``fluxCalibrate()`` method:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp.fluxCalibrate('2MASS J',14.0)

* To display the spectrum, use the Spectrum object's plot function 

>>> sp.plot()
 
You can save this display by adding a filename:

>>> sp.plot(sp,file='spectrum.png')


* Spectral math






Spectral Analysis
-----------------

SPLAT has several routines to do basic spectral analysis and combining of spectra.


Spectrophotometry
-----------------

SPLAT allows spectrophotometry of spectra using common filters in the red optical and near-infared. 



Classification
-----------------

SPLAT contains several different methods for classifying a spectrum, using both indices and 


Potentially Useful Program Constants
------------------------------------

``splat.DB_SOURCES``
	An Astropy Table object containing the Source Database
	
``splat.DB_SPECTRA``
	An Astropy Table object containing the Spectrum Database
	
``splat.SPEX_STDS``
	A dictionary containing Spectrum objects of the M0-T9 dwarf standards; this dictionary is 
	populated through calls to ``splat.getStandard``. A standard Spectrum object can be accessed
	by using the spectral type as the referring key:

>>> sp = splat.getStandard('M0')[0]		# both are Spectrum objects of Gliese 270
>>> sp = splat.SPEX_STDS['M0.0']		# note the mandatory decimal

	Available standards can be accessed through the command:

>>> splat.SPEX_STDS.keys()

	
``splat.SPEX_SD_STDS``
	Same as ``splat.SPEX_STDS`` for subdwarf standards

``splat.SPEX_ESD_STDS``
	Same as ``splat.SPEX_STDS`` for extreme subdwarf standards

``splat.FILTERS``
	A dictionary containing information on all of the filters used in SPLAT photometry. The command:

>>> splat.FILTERS.keys()

	will produce an (unordered) list of filters
	


Additional Programs
----------------------





* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

