.. SpeX Prism Library Analysis Toolkit documentation master file, created by
   sphinx-quickstart on Sat Jul 11 20:07:28 2015.


SPLAT: The SpeX Prism Library Analysis Toolkit
===============================================


SPLAT is a python-based spectral access and analysis package designed to interface  
with the SpeX Prism Library (SPL: http://www.browndwarfs.org/spexprism), 
an online repository of over
1500 low-resolution, near-infrared spectra of low-temperature stars and brown dwarfs.
It is built on common python packages such as astropy, matplotlib, numpy, pandas and scipy.  

SPLAT tools allow you to:
    * search the SPL for data and source information;
    * access the publically-available (published) spectra contained in it;
    * compare your near-infrared spectrum to these data;
    * make use of published empirical trends in absolute magnitudes and effective temperatures;
    * perform basic spectral analyses such as spectral classification, gravity classification, index measurement, spectrophotometry, reddening, robust comparison statistics, basic math operations;
    * perform advanced analyses such as MCMC spectral model fitting;
    * transform observables using empirical trends;
    * transform observable to physical parameters using evolutionary models; and
    * plot/tabulate/publish your results.  

**Note that many of these features are currently under development.**

Installation and Dependencies
-----------------------------

SPLAT is best forked from the github site http://github.org/aburgasser/splat, 
which is updated on a regular basis.
SPLAT has not yet reached v1.0, so bugs are common. Please help us squish them by 
sending bug reports to aburgasser@ucsd.edu or start an issue on the github site.

You may also obtain splat using pip_:

.._pip: https://pip.pypa.io/en/stable/

>>> pip install splat

Instructions on setting up and using SPLAT are maintained at http://www.browndwarfs.org/splat.

Copy the file ``.splat_access`` into your home directory - this is your access key
if you have priveleged access to unpublished data in the SPL.

Using SPLAT
-----------

SPLAT is best used in the **ipython** or **ipython notebook**; all of the necessary data is
included in the github/pip install, so you don't need to be online to run most programs.

Here are some examples:

* The best way to read in a spectrum is to use getSpectrum_:

.. _getSpectrum: api.html#splat.getSpectrum

>>> import splat
>>> splist = splat.getSpectrum(shortname='0415-0935')
>>> splist = splat.getSpectrum(young=True)
>>> splist = splat.getSpectrum(spt=['M7','L5'],jmag=[14.,99.])

In each case, splist is a list of Spectrum_ objects, which is the container of various 
aspects of the spectrum and it source properties. For example, selecting the first spectrum,

.. _Spectrum: splat.html?highlight=spectrum#the-splat-spectrum-object

>>> sp = splist[0]

``sp.wave`` gives the wavelengths of this spectrum, ``sp.flux`` the flux values, and ``sp.noise`` the 
flux uncertainty. There are several other elements to the Spectrum_ object that can be accessed using ``sp.info()``.

You can also read in your own spectrum by passing a filename

>>> sp = splat.Spectrum(filename='PATH_TO/myspectrum.fits')

Note that this file must conform to the following standard: the first column is
wavelength in microns, second column flux in f_lambda units, third column (optional) is 
flux uncertainty in f_lambda units.

* To flux calibrate the spectrum, use the object's built in fluxCalibrate_ method:

.. _fluxCalibrate: api.html#splat.Spectrum.fluxCalibrate

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp.fluxCalibrate('2MASS J',14.0)

* To display the spectrum, use the Spectrum object's plot_ function or plotSpectrum_ :

.. _plot: api.html#splat.Spectrum.plot
.. _plotSpectrum: api.html#splat_plot.plotSpectrum

>>> sp.plot()
>>> splat.plotSpectrum(sp)

which will pop up a window displaying flux vs. wavelength. 
You can save this display by adding a filename:

>>> splat.plotSpectrum(sp,file='spectrum.png')

You can also compare multiple spectra:

>>> sp1 = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp2 = splat.getSpectrum(shortname='1217-0311')[0]
>>> splat.plotSpectrum(sp1,sp2,colors=['black','red'])

You can add several extras to this to label features, plot uncertainties, 
indicate telluric absorption regions, make multi-panel and multi-page plots
of lists of spectra, etc. Be sure to look through the plotting 
subpackage for more details.


SPLAT can analyze and compare an arbitrary number of spectra.

* To measure spectral indices, use measureIndex_ or measureIndexSet_:

.. _measureIndex: api.html#splat.measureIndex
.. _measureIndexSet: api.html#splat.measureIndexSet

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> value, error = splat.measureIndex(sp,[1.14,1.165],[1.21,1.235],method='integrate')
>>> indices = splat.measureIndexSet(sp,set='testi')

The last line returns a dictionary, whose value,error pair can be accessed by the name 
of the index:

>>> print indices['sH2O-J']		# returns value, error

* You can also determine the gravity classification of a source via `Allers & Liu (2013) <http://adsabs.harvard.edu/abs/2013ApJ...772...79A>`_ using classifyGravity_:

.. _classifyGravity: api.html#splat.classifyGravity

>>> sp = splat.getSpectrum(young=True, lucky=True)[0]
>>> print splat.classifyGravity(sp)   # returned 'VL-G'


* To classify a spectrum, use the various classifyByXXX_ methods:

.. _classifyByXXX: api.html#spectral-classification

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> spt,unc = splat.classifyByIndex(sp,set='burgasser')
>>> spt,unc = splat.classifyByStandard(sp,spt=['T5','T9'])
>>> result = splat.classifyByTemplate(sp,spt=['T6','T9'],nbest=5)

The last line returns a dictionary containing the best 5 template matches to the Spectrum ``sp``.


* To compare a spectrum to another spectrum or a model, use compareSpectra_ :

.. _compareSpectra: api.html#splat.compareSpectra

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> mdl = splat.loadModel(teff=700,logg=5.0)			# loads a BTSettl08 model by default
>>> chi,scale = splat.compareSpectra(sp,mdl)
>>> mdl.scale(scale)
>>> splat.plotSpectrum(sp,mdl,colors=['black','red'],legend=[sp.name,mdl.name])


# There is also a basic Markov Chain Monte Carlo code to compare models to spectra called modelFitMCMC_:

.. _modelFitMCMC: api.html#splat_model.modelFitMCMC

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp.fluxCalibrate('2MASS J',14.49,absolute=True)
>>> table = splat.modelFitMCMC(sp,initial_guess=[800,5.0,0.],nsamples=300,step_sizes=[50.,0.5,0.])


All of these routines have many options worth exploring, and which are (increasingly) documented on this website. If there are capabilities
you need, please suggest them to aburgasser@ucsd.edu, or note it in the "Issues" link on our `github site <https://github.com/aburgasser/splat>`_.

Acknowledgements
----------------

SPLAT is an experimental, collaborative project of research students in `Adam Burgasser's
UCSD Cool Star Lab <http://www.coolstarlab.org>`_, aimed at teaching students how to do research by building their own analysis tools.  Contributors to SPLAT have included Christian Aganze, Jessica Birky, Daniella Bardalez Gagliuffi, Adam Burgasser (PI), Caleb Choban, Andrew Davis, Ivanna Escala, Aishwarya Iyer, Yuhui Jin, Mike Lopez, Alex Mendez, Gretel Mercado, Elizabeth Moreno, Johnny Parra, Maitrayee Sahi, Adrian Suarez, Melisa Tallis, Tomoki Tamiya, Chris Theissen and Russell van Linge.

This project is supported by the National Aeronautics and Space Administration under Grant No. NNX15AI75G.

*Contents*

.. toctree::
   :maxdepth: 3

   installation
   quickstart
   splat
   splat_db
   splat_plot
   splat_model
   splat_evolve
   api
   :ref:`genindex`
   
*Search*

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

