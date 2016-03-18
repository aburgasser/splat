Quickstart
===============================================

SPLAT is best used in the ipython or ipython notebook; all of the necessary data is
included in the github install, so you shouldn't need to be online to run anything
unless you are using proprietary data (these are not distributed with the package).

Here are some examples:

* The best way to read in a spectrum is to use getSpectrum:

>>> import splat
>>> splist = splat.getSpectrum(shortname='0415-0935')
>>> splist = splat.getSpectrum(young=True)
>>> splist = splat.getSpectrum(spt=['M7','L5'],jmag=[14.,99.])

In each case, splist is a list of Spectrum objects, which is the container of various 
aspects of the spectrum and it source properties. For example, selecting the first spectrum,


>>> sp = splist[0]

``sp.wave`` gives the wavelengths of this spectrum, ``sp.flux`` the flux values, and ``sp.noise`` the 
flux uncertainty.

You can also read in your own spectrum by passing a filename

>>> sp = splat.Spectrum(filename='PATH_TO/myspectrum.fits')

Note that this file must conform to the standard of the SPL data: the first column is
wavelength in microns, second column flux in f_lambda units, third column (optional) is 
flux uncertainty.

* To flux calibrate the spectrum, use the object's built in ``fluxCalibrate()`` method:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp.fluxCalibrate('2MASS J',14.0)

* To display the spectrum, use the Spectrum object's plot function or plotSpectrum

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
of lists of spectra, etc. Be sure to look through the SPLAT plotting 
subpackage for more details.


SPLAT can analyze and compare an arbitrary number of spectra.

* To measure spectral indices, use measureIndex or measureIndexSet:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> value, error = splat.measureIndex(sp,[1.14,1.165],[1.21,1.235],method='integrate')
>>> indices = splat.measureIndexSet(sp,set='testi')

The last line returns a dictionary, whose value,error pair can be accessed by the name 
of the index:

>>> print indices['sH2O-J']		# returns value, error

* You can also determine the gravity classification of a source via Allers & Liu (2013):

>>> sp = splat.getSpectrum(young=True, lucky=True)[0]
>>> splat.classifyGravity(sp,verbose=True)

This returns:

>>> Gravity Classification:
>>>   SpT = L1.0
>>>   VO-z: 1.193+/-0.018 => 1.0
>>>   FeH-z: 1.096+/-0.026 => 2.0
>>>   H-cont: 0.973+/-0.010 => 2.0
>>>   KI-J: 1.044+/-0.008 => 2.0
>>>   Gravity Class = VL-G


* To classify a spectrum, use the classifyByXXX methods:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> spt,unc = splat.classifyByIndex(sp,set='burgasser')
>>> spt,unc = splat.classifyByStandard(sp,spt=['T5','T9'])
>>> bestMatches = splat.classifyByTemplate(sp,spt=['T6','T9'],nbest=5)

The last line returns a dictionary containing the best 5 template matches to the Spectrum sp. 
Note that this can take a long time to run!


* To compare a spectrum to another spectrum or a model, use compareSpectra:

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> mdl = splat.loadModel(teff=700,logg=5.0)			# BTSettl08 model by default
>>> chi,scale = splat.compareSpectra(sp,mdl)
>>> mdl.scale(scale)
>>> splat.plotSpectrum(sp,mdl,colors=['black','red'],legend=[sp.name,mdl.name])

The available spectral models are 

  * 'BTSettl08' (Allard et al. 2008)
  * 'drift' (Witte et al. 2008)
  * 'burrows06' (Burrows et al. 2006)
  * 'saumon12' (Saumon & Marley 2012)
  * 'morley12' (Morley et al. 2012)
  * 'morley14; (Morley et al. 2014)

* There is also a basic Markov Chain Monte Carlo code to compare models to spectra (Note: still in development)

>>> sp = splat.getSpectrum(shortname='0415-0935')[0]
>>> sp.fluxCalibrate('2MASS J',14.49,absolute=True)
>>> table = splat.modelFitMCMC(sp,initial_guess=[800,5.0,0.],nsamples=300,step_sizes=[50.,0.5,0.])


All of these routines have many options worth exploring, and which are (partially) documented 
in the following pages. If there are other capabilities
you need, please suggest them, or note it in the "Issues" link on our github site



*Search*





* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

