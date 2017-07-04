.. _`compareSpectra()` : api.html#splat.core.compareSpectra
.. _`classifyByIndex()` : api.html#splat.core.classifyByIndex
.. _`classifyByStandard()` : api.html#splat.core.classifyByStandard
.. _`classifyByTemplate()` : api.html#splat.core.classifyByTemplate
.. _`getSpectrum()` : api.html#splat.core.getSpectrum
.. _`measureIndex()` : api.html#splat.core.measureIndex
.. _`measureIndexSet()` : api.html#splat.core.measureIndexSet
.. _`plotSpectrum()` : api.html#splat.plot.plotSpectrum
.. _`fluxCalibrate()` : api.html#splat.core.Spectrum.fluxCalibrate
.. _`plot()` : api.html#splat.core.Spectrum.plot
.. _`Spectrum` : api.html#splat.core.Spectrum
.. _`modelFitGrid()` : api.html#splat.model.modelFitGrid
.. _`modelFitMCMC()` : api.html#splat.model.modelFitMCMC
.. _`modelFitEMCEE()` : api.html#splat.model.modelFitEMCEE

Quickstart
===============================================

SPLAT is best used in the ipython or jupyter notebook; all of the necessary data is
included in the github install, so you shouldn't need to be online to run anything
unless you are using proprietary data (these are not distributed with the package)
or the web interface.

All of the routines in SPLAT are organized into packages:

  * **splat.core**: core functionalities, including index measurement, database access and classification
  * **splat.citations**: biblographic/bibtex routines
  * **splat.database**: access the spectral and source databases, as well as online resources through astroquery
  * **splat.empirical**: empirical conversion relations
  * **splat.evolve**: access to evolutionary models and population synthesis routines
  * **splat.model**: access to spectral models and model-fitting routines
  * **splat.photometry**: spectrophotometry routines
  * **splat.plot**: plotting and visualization routines
  * **splat.utilities**: additional routines for general analysis
  * **splat.web**: SPLAT's web interface

To access routines, you should first import the appropriate package. 

Most routines also have help docstrings that can be accessed by appending a question mark to the end of the function call:

>>> splat.getSpectrum?


Examples
^^^^^^^^

* The best way to read in a spectrum from the SPLAT library is to use `getSpectrum()`_

>>> import splat
>>> splist = splat.getSpectrum(shortname='0415-0935')
>>> splist = splat.getSpectrum(young=True)
>>> splist = splat.getSpectrum(spt=['M7','L5'],jmag=[14.,99.])
>>> splist = splat.getSpectrum(lucky=True)

The last case returns a randome spectrum. 
In each case, ``splist`` is a list of Spectrum objects, which is the container of various 
aspects of the spectrum and its source properties. For example, selecting the first spectrum,


>>> sp = splist[0]

``sp.wave`` contains the wavelengths of this spectrum, ``sp.flux`` the flux values, and ``sp.noise`` the 
flux uncertainties.

You can also read in your own spectrum by passing a filename directly to the `Spectrum`_ class object:

>>> sp = splat.Spectrum(filename='PATH_TO/myspectrum.fits')

Note that this file (fits or ascii) must conform to the standard of the SPL data: the first column is
wavelength in microns, second column flux in f_lambda units, third column (optional) is 
flux uncertainty. You can preset the units for wavelength and flux density with optional keywords and using the 
`astropy <http://www.astropy.org/>`_ units package:

>>> import astropy.units as u
>>> sp = splat.Spectrum(filename='PATH_TO/myspectrum.fits',wunit=u.Angstrom,funit=u.erg/u.cm/u.cm/u.sec/u.Angstrom)

* To flux calibrate the spectrum, use the object's built in `fluxCalibrate()`_ method:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp.fluxCalibrate('2MASS J',14.0)

* To display the spectrum, use the Spectrum object's `plot()`_ function or `plotSpectrum()`_

>>> sp.plot()
>>> from splat.plot import plotSpectrum
>>> plotSpectrum(sp)

which will pop up a window displaying flux vs. wavelength. 
You can save this display by adding a filename:

>>> sp.plot(sp,file='spectrum.png')

You can also compare multiple spectra:

>>> sp1 = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp2 = splat.getSpectrum(shortname='1217-0311')[0]
>>> from splat.plot import plotSpectrum
>>> plotSpectrum(sp1,sp2,colors=['black','red'])

You can add several extras to this to label features, plot uncertainties, 
indicate telluric absorption regions, make multi-panel and multi-page plots
of lists of spectra, etc. Be sure to look through the `SPLAT plotting 
subpackage <splat_plot.html>`_ for more details.


SPLAT can analyze and compare an arbitrary number of spectra.

* To measure spectral indices, use `measureIndex()`_ or `measureIndexSet()`_:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> value, error = splat.measureIndex(sp,[1.14,1.165],[1.21,1.235],method='integrate')
>>> indices = splat.measureIndexSet(sp,set='testi')

The last line returns a dictionary, whose value,error pair can be accessed by the name 
of the index:

>>> print(indices['sH2O-J'])		# returns value, error

* You can also determine the gravity classification of a source via `Allers & Liu (2013) <http://adsabs.harvard.edu/abs/2013ApJ...772...79A>`_:

>>> sp = splat.getSpectrum(young=True, lucky=True)[0]
>>> splat.classifyGravity(sp,verbose=True)

This returns (depending on the source returned):

>>> Gravity Classification:
>>>   SpT = L1.0
>>>   VO-z: 1.193+/-0.018 => 1.0
>>>   FeH-z: 1.096+/-0.026 => 2.0
>>>   H-cont: 0.973+/-0.010 => 2.0
>>>   KI-J: 1.044+/-0.008 => 2.0
>>>   Gravity Class = VL-G


* To classify a spectrum, use the `classifyByStandard()`_, `classifyByIndex()`_, or `classifyByTemplate()`_ methods:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> spt,unc = splat.classifyByStandard(sp,spt=['T5','T9'])
>>> spt,unc = splat.classifyByIndex(sp,set='burgasser')
>>> bestMatches = splat.classifyByTemplate(sp,spt=['T6','T9'],nbest=5)

The last line returns a dictionary containing the best 5 template matches to the Spectrum sp. 
Note that comparing to the large template library of SPLAT can take a long time to run!


* To compare a spectrum to another spectrum or a model, use `compareSpectra()`_:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> from splat.model import loadModel			# loads in model reading package
>>> mdl = splat.loadModel(teff=700,logg=5.0)	# BTSettl08 model by default
>>> chi,scale = splat.compareSpectra(sp,mdl,plot=True,file='comparison.pdf')

Notice in the last line, the keyword ``plot`` allows you to visualize the comparison, which can be saved as a file

The available spectral models are 

	* *BTSettl2008* from `Allard et al. (2012) <http://adsabs.harvard.edu/abs/2012RSPTA.370.2765A>`_  (default)
	* *burrows06* from `Burrows et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...640.1063B>`_ 
	* *saumon12* from `Saumon et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...750...74S>`_ 
	* *drift* from `Witte et al. (2011) <http://adsabs.harvard.edu/abs/2011A%26A...529A..44W>`_ 
	* *madhusudhan* from `Madhusudhan et al. (2011) <http://adsabs.harvard.edu/abs/2011ApJ...737...34M>`_ 
	* *morley12* from `Morley et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...756..172M>`_ 
	* *morley14* from `Morley et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...787...78M>`_ 

* You can fit models to the spectral data using one of three built-in model-fitting routines: `modelFitGrid()`_, `modelFitMCMC()`_, and `modelFitEMCEE()`_:

>>> import splat.model as spmodel
>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> results = spmod.modelFitGrid(sp,teff_range=[1200,2500],model='Saumon',file='fit1507')
>>> results = spmod.modelFitMCMC(sp,initial_guess=[900,5.0,0.0],nsamples=1000)
>>> results = modelFitEMCEE(sp,t0=900,g0=5.0,nwalkers=50,nsamples=500,output='basefilename')


The output of each of these is a dictionary containing the best fit model parameters, average fit model parameters, and visualization plots. Refer to the API for these for more details.



All of these routines have many options worth exploring, and which are (partially) documented 
in the following pages. If there are other capabilities
you need, please suggest them, or note it in the "Issues" link on our github site



*Search*





* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

