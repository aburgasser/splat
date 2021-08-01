.. SpeX Prism Library Analysis Toolkit documentation master file, created by
   sphinx-quickstart on Sat Jul 11 20:07:28 2015.

.. _`SpeX Prism Library`: http://www.browndwarfs.org/spexprism
.. _`SPL`: http://www.browndwarfs.org/spexprism
.. _`pip`: https://pip.pypa.io/en/stable/


SPLAT: The SpeX Prism Library Analysis Toolkit
===============================================


SPLAT is a python-based spectral access and analysis package designed to interface  
with the `SpeX Prism Library`_ (`SPL`_), an online repository of over
3000 low-resolution, near-infrared spectra of low-temperature stars and brown dwarfs.
It is built on common python packages such as `astropy <http://www.astropy.org/>`_, `emcee <http://dan.iel.fm/emcee/current/>`_, `matplotlib <http://matplotlib.org/>`_, `numpy <http://www.numpy.org/>`_, `pandas <http://pandas.pydata.org/>`_, `scipy <https://www.scipy.org/>`_, and others.  


SPLAT tools allow you to:
    * search the SPL for spectral data and source information;
    * access the publically-available spectra contained in it;
    * compare your own near-infrared spectrum to these data;
    * make use of published empirical trends in absolute magnitudes, luminosity, effective temperatures, and others;
    * perform basic spectral analyses such as spectral classification, gravity classification, index measurement, spectrophotometry, reddening, blended light analysis, and basic math operations;
    * access atmosphere models and perform fits to spectral data;
    * transform observables to physical parameters using evolutionary models; 
    * simulate very low mass star and brown dwarf populations; and
    * plot, tabulate and publish your results.  

**Note that many of these features continue to be in development.**

Installation and Dependencies
-----------------------------

SPLAT should be cloned from the github site `https://github.com/aburgasser/splat <https://github.com/aburgasser/splat>`_.
which is updated on a regular basis. 
**Warning** at this time please do not try to install splat using `pip`_, as this is an outdated version that is no longer supported.

Once you've downloaded the code and data, you will need to add the SPLAT top-level directory to the environment variables ``SPLAT_PATH`` and ``PYTHONPATH`` (and optionally to your system ``PATH``).  More detailed instructions are on the `installation <installation.html>`_ page. 

SPLAT has core dependencies on the following packages:
    * `astropy <http://www.astropy.org/>`_
    * `astroquery <https://astroquery.readthedocs.io/en/latest/>`_
    * `bokeh <http://bokeh.pydata.org/en/latest/>`_ (for SPLAT web interface only)
    * `corner <http://corner.readthedocs.io/en/latest/>`_
    * `emcee <http://dan.iel.fm/emcee/current/>`_ (for emcee model fitting only)
    * `flask <http://flask.pocoo.org/>`_ (for SPLAT web interface only)
    * `matplotlib <http://matplotlib.org/>`_
    * `numpy <http://www.numpy.org/>`_
    * `pandas <http://pandas.pydata.org/>`_
    * `requests <http://docs.python-requests.org/en/master/>`_
    * `scipy <https://www.scipy.org/>`_

SPLAT has not yet reached v1.0, so bugs are common. Please help us squish them by 
sending bug reports to aburgasser@ucsd.edu or start an issue on the github site.

Using SPLAT
-----------

.. _`Spectrum`: splat.html?highlight=spectrum#the-splat-spectrum-object
.. _`getSpectrum()`: api.html#splat.getSpectrum
.. _`fluxCalibrate()`: api.html#splat.Spectrum.fluxCalibrate
.. _`plot()`: api.html#splat.Spectrum.plot
.. _`plotSpectrum()`: api.html#splat.plot.plotSpectrum
.. _`measureIndex()`: api.html#splat.measureIndex
.. _`measureIndexSet()`: api.html#splat.measureIndexSet
.. _`classifyGravity()`: api.html#splat.classifyGravity
.. _`classifyByXXX`: api.html#spectral-classification
.. _`compareSpectra()`: api.html#splat.compareSpectra
.. _`modelFitMCMC()`: api.html#splat.model.modelFitMCMC


SPLAT is organized into a series of modules based on core functionalities:
  * `splat.core <splat.html>`_: core functionalities, including index measurement, database access and classification
  * `splat.citations <splat_citations.html>`_: biblographic/bibtex routines
  * `splat.database <splat_database.html>`_: access the spectral and source databases, as well as online resources through astroquery
  * `splat.empirical <splat_empirical.html>`_: empirical conversion relations
  * `splat.evolve <splat_evolve.html>`_: access to evolutionary models
  * `splat.model <splat_model.html>`_: access to spectral models and model-fitting routines
  * `splat.photometry <splat_photometry.html>`_: spectrophotometry routines and filter access
  * `splat.plot <splat_plot.html>`_: plotting and visualization routines
  * `splat.simulate <splat_simulate.html>`_: population simulation routines
  * `splat.utilities <splat_utilities.html>`_: additional routines for general analysis
  * `splat.web <splat_web.html>`_: SPLAT's web interface

SPLAT has been tested on both Python 2.7 and 3.5, and is best used in the 
**ipython** or **ipython notebook**; all of the necessary data is
included in the github package, so you don't need to be online to run most programs.

Here are some examples:

* The best way to read in a spectrum is to use `getSpectrum()`_, which takes a number of search keywords and
returns a list of `Spectrum`_ objects:

>>> import splat
>>> splist = splat.getSpectrum(shortname='0415-0935')  
>>> splist = splat.getSpectrum(young=True)  
>>> splist = splat.getSpectrum(spt=['M7','L5'],jmag=[14.,99.])

In each case, splist is a list of `Spectrum`_ objects, each a container of various aspects of each spectrum and its source properties. For example, selecting the first spectrum,

>>> sp = splist[0]

``sp.wave`` gives the wavelengths of this spectrum, ``sp.flux`` the flux values, and ``sp.noise`` the 
flux uncertainty. A summary of the `Spectrum`_ object can be accessed using ``sp.info()``.

You can also read in your own spectrum by passing a filename

>>> sp = splat.Spectrum(filename='PATH_TO/myspectrum.fits')

Both fits and ascii (tab or csv) data formats are supported, but files 
should conform to the following data format standard: 
    * column 1: wavelength, assumed in microns
    * column 2: flux in f_lambda units
    * column 3: (optional) flux uncertainty in f_lambda units.

* To flux calibrate the spectrum, use the `Spectrum`_ object's built in `fluxCalibrate()`_ method:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp.fluxCalibrate('2MASS J',14.0)

* To display the spectrum, use the Spectrum object's `plot()`_ function 

>>> sp.plot()

or the splat.plot routine `plotSpectrum()`_ :

>>> import splat.plot as splot
>>> splot.plotSpectrum(sp)

You can save your spectrum by adding a filename:

>>> splot.plotSpectrum(sp,file='spectrum.pdf')

You can also compare multiple spectra:

>>> sp1 = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp2 = splat.getSpectrum(shortname='1217-0311')[0]
>>> splot.plotSpectrum(sp1,sp2,colors=['k','r'])

`plotSpectrum()`_ and related routines have many extras to label features, plot uncertainties, 
indicate telluric absorption regions, make multi-panel and multi-page plots
of lists of spectra, plot batches of spectra, etc. Be sure to look through the `splat.plot`_ 
subpackage for more details.


SPLAT's primary purpose is to allow the analysis of ultracool dwarf spectra.

* To measure spectral indices, use `measureIndex()`_ or `measureIndexSet()`_:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> value, error = splat.measureIndex(sp,[1.14,1.165],[1.21,1.235],method='integrate')
>>> indices = splat.measureIndexSet(sp,set='testi')

The last line returns a dictionary, whose value,error pair can be accessed by the name 
of the index:

>>> print(indices['sH2O-J'])		# returns value, error

* You can also determine the gravity classification of a source following `Allers & Liu (2013) <http://adsabs.harvard.edu/abs/2013ApJ...772...79A>`_ using `classifyGravity()`_:

>>> sp = splat.getSpectrum(young=True, lucky=True)[0]
>>> print(splat.classifyGravity(sp))   # returned 'VL-G'

* To classify a spectrum, use the various `classifyByXXX`_ methods:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> spt,unc = splat.classifyByIndex(sp,set='burgasser')
>>> spt,unc = splat.classifyByStandard(sp,spt=['T5','T9'])
>>> result = splat.classifyByTemplate(sp,spt=['T6','T9'],nbest=5)

The last line returns a dictionary containing the best 5 template matches.

* To compare a spectrum to another spectrum or a model, use `compareSpectra()`_ :

>>> import splat.model as spmod
>>> mdl = spmod.loadModel(teff=720,logg=4.8,set='btsettl')      # loads a BTSettl08 model 
>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> chi,scale = splat.compareSpectra(sp,mdl)
>>> mdl.scale(scale)
>>> splat.plotSpectrum(sp,mdl,colors=['k','r'],legend=[sp.name,mdl.name])

You can shortcut the last three lines using the ``plot`` keyword:

>>> chi,scale = splat.compareSpectra(sp,mdl,plot=True)


* There are also codes **still in development** to fit models directly to spectra: `modelFitGrid()`_, `modelFitMCMC()`_, and `modelFitEMCEE()`_:

>>> import splat.model as spmod
>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp.fluxCalibrate('2MASS J',14.49,absolute=True)
>>> nbest = 5
>>> result1 = splat.modelFitGrid(sp,set='btsettl')
>>> result2 = splat.modelFitMCMC(sp,set='btsettl',initial_guess=[800,5.0,0.],nsamples=300,step_sizes=[50.,0.5,0.])
>>> result3 = splat.modelFitEMCEE(sp,set='btsettl',initial_guess=[800,5.0,0.],nwalkers=12,nsamples=500)

The outputs of all of these fitting functions is a dictionary or list of dictionaries containing the parameters of the best-fitting models; there are also several diagnostic plots produced depending on the routine. View the model fitting page for more details.

All of these routines have many options worth exploring, and which are (increasingly) documented on this website. If there are capabilities
you need, please suggest them to aburgasser@ucsd.edu, or note it in the "Issues" link on our `github site <https://github.com/aburgasser/splat>`_.

Acknowledgements
----------------

SPLAT is an collaborative project of research students in the `UCSD Cool Star Lab <http://www.coolstarlab.org>`_, aimed at developing research through the building of spectral analysis tools.  Contributors to SPLAT have included Christian Aganze, Jessica Birky, Daniella Bardalez Gagliuffi, Adam Burgasser (PI), Caleb Choban, Andrew Davis, Ivanna Escala, Joshua Hazlett, Carolina Herrara Hernandez, Elizabeth Moreno Hilario, Aishwarya Iyer, Yuhui Jin, Mike Lopez, Dorsa Majidi, Diego Octavio Talavera Maya, Alex Mendez, Gretel Mercado, Niana Mohammed, Johnny Parra, Maitrayee Sahi, Adrian Suarez, Melisa Tallis, Tomoki Tamiya, Chris Theissen, and Russell van Linge.

This project has been supported by the National Aeronautics and Space Administration under Grant No. NNX15AI75G.

*Contents*

.. toctree::
   :maxdepth: 3

 
