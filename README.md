# SPLAT: The SpeX Prism Library Analysis Toolkit

Access SPLAT's full documentation at https://splat.physics.ucsd.edu/splat.

## Preamble

SPLAT is a python-based spectral access and analysis package designed to interface  
with the **SpeX Prism Library (SPL)**, an online repository of over
3,000 low-resolution, near-infrared spectra, primarily 
of low-temperature stars and brown dwarfs.
It is built on common python packages such as 
[astropy](http://www.astropy.org), 
[astroquery](https://astroquery.readthedocs.io/en/latest), 
[emcee](http://dan.iel.fm/emcee/current), 
[matplotlib](http://matplotlib.org), 
[numpy](http://www.numpy.org), 
[pandas](http://pandas.pydata.org), 
[scipy](https://www.scipy.org), and others.  


SPLAT tools allow you to:
* Search the SpeX Prism Library for spectral data and source information;
* Access and analyze publically-available spectra contained in it;
* Analyze your own spectral data from SpeX and other instruments;
* Perform basic spectral analyses such as type classification, gravity classification, index measurement, spectrophotometry, reddening, blended light analysis, and basic math operations;
* Access atmosphere models and perform fits to spectral data;
* Transform observables to physical parameters using evolutionary models; 
* Use published empirical trends between spectral type, absolute magnitudes, colors, luminosities, effective temperatures, and others;
* Access online data repositories through wrappers to [astroquery] (https://astroquery.readthedocs.io/en/latest)
* Simulate very low mass star and brown dwarf populations by combining spatial, evolutionary, and observational properties; and
* Plot, tabulate, and publish your results.  

> Note: Many features in SPLAT continue to be in development.
    Help us improve the code by reporting bugs (and solutions!) to our github site,
    https://github.com/aburgasser/splat.


## Installation and Dependencies

The optimal installation method for SPLAT is cloning from the github site https://github.com/aburgasser/splat, which is updated on a (semi-)regular basis. The installations uses a pyproject.toml to install dependencies, so recommended set up is as follows:


    conda create -n splat
    conda activate splat
    git clone https://github.com/aburgasser/splat.git
    cd splat
    python -m pip install .


> Warning: At this time please do not install splat using `pip`, as this is an outdated version of SPLAT that is no longer supported.

Once you've downloaded the code and data, you should add the SPLAT top-level directory to the environment variable `PYTHONPATH`.  See the following links on how to update environment variables:

* bash or zsh (Unix/Mac): https://apple.stackexchange.com/questions/356441/how-to-add-permanent-environment-variable-in-zsh
* windows: https://www.computerhope.com/issues/ch000549.htm 



SPLAT has core dependencies on the following packages:
* [astropy](http://www.astropy.org)
* [astroquery](https://astroquery.readthedocs.io/en/latest)
* [matplotlib](http://matplotlib.org)
* [numpy](http://www.numpy.org)
* [pandas](http://pandas.pydata.org)
* [requests](http://docs.python-requests.org/en/master)
* [scipy](https://www.scipy.org)
* [corner](http://corner.readthedocs.io/en/latest)  (for model fitting only)
* [emcee](http://dan.iel.fm/emcee/current) (for model fitting only)
* [bokeh](http://bokeh.pydata.org/en/latest) (for SPLAT web interface only)
* [flask](http://flask.pocoo.org) (for SPLAT web interface only)

## Using SPLAT


SPLAT is organized into a series of modules based on core functionalities:
* `splat.core`: core functionalities, including index measurement, database access and classification
* `splat.citations`: biblographic/bibtex routines
* `splat.database`: access the spectral and source databases, as well as online resources through astroquery
* `splat.empirical`: empirical conversion relations
* `splat.evolve`: access to evolutionary models
* `splat.model`: access to spectral models and model-fitting routines
* `splat.photometry`: spectrophotometry routines and filter access
* `splat.plot`: plotting and visualization routines
* `splat.simulate`: population simulation routines
* `splat.utilities`: additional routines for general analysis
* `splat.web`: SPLAT's web interface (under development)

SPLAT has been tested on both Python 2.7 and 3.0-3.11, and is best used in 
`ipython` or `jupyter notebook`.
All of the necessary data is included in the github package, so you don't need to be online to run most programs.


### Reading in Spectra

The best way to read in a spectrum is to use `getSpectrum()`, which takes a number of search keywords and returns a list of `Spectrum` objects:

    import splat
    splist = splat.getSpectrum(shortname='0415-0935')  
> Retrieving 1 file

    splist = splat.getSpectrum(name='TWA30A')  
> Retrieving 3 files

    splist = splat.getSpectrum(opt_spt=['L2','L5'],jmag=[12,13])
> Retrieving 5 files

In each case, `splist` is a list of `Spectrum` objects, each a container of various aspects of each spectrum and its source properties. For example, selecting the first spectrum,

    sp = splist[0]
    sp

> SPEX-PRISM spectrum of 2MASSW J0036159+182110

The main elements of the `Spectrum` obejct are:
* `sp.wave`: wavelength array in default units of micron
* `sp.flux`: flux array in default units of erg/cm^2/s/micron
* `sp.noise`: flux uncertainty array in default units of erg/cm^2/s/micron

A summary of the `Spectrum` object can be accessed using `sp.info()`.

    sp.info()

> SPEX-PRISM spectrum of 2MASSW J0036159+182110 
>
>  Airmass = nan
> 
> Source designation = J00361617+1821104
>
> Median S/N = 274
>
> SpeX Classification = L2.0
>
> Spectrum key = 10249, Source key = 10068
>
> If you use these data, please cite:
>
>> Burgasser, A. J. et al. (2008, Astrophysical Journal, 681, 579-593)
>>
>> bibcode: 2008ApJ...681..579B
>
> History:
>
>> SPEX-PRISM spectrum successfully loaded

You can also read in your own spectrum by passing a filename

    sp = splat.Spectrum(file='PATH_TO/myspectrum.fits')

or a URL

    sp = splat.Spectrum(file='http://splat.physics.ucsd.edu/splat/spexprism/spectra/spex-prism_SO0253+1625_20040908_BUR08B.txt')


Both fits and ascii (tab or csv) data formats are supported, but files 
should ideally conform to the following data format standard: 
* column 1: wavelength, assumed in microns
* column 2: flux in flambda units
* column 3: (optional) flux uncertainty in flambda units.

There are a few built-in readers for specific data formats.

To flux calibrate a spectrum, use the `Spectrum` object's built in `fluxCalibrate()` method:

    sp = splat.getSpectrum(shortname='0415-0935')[0]
    sp.fluxCalibrate('2MASS J',14.0)


### Visualizing Spectra

To display the spectrum, use the Spectrum object's `plot()` function 

    sp.plot()

or the splat.plot routine `plotSpectrum()` :

    import splat.plot as splot
    splot.plotSpectrum(sp)

You can save your spectrum by adding a filename:

    splot.plotSpectrum(sp,file='spectrum.pdf')

You can also compare multiple spectra:

    sp1 = splat.getSpectrum(shortname='0415-0935')[0]
    sp2 = splat.getSpectrum(shortname='1217-0311')[0]
    splot.plotSpectrum(sp1,sp2,colors=['k','r'])

`plotSpectrum()` and related routines have many extras to label features, plot uncertainties, 
indicate telluric absorption regions, make multi-panel and multi-page plots
of lists of spectra, plot batches of spectra, etc. Be sure to look through the `splat.plot`
subpackage for more details.

### Analysis functions


SPLAT's primary purpose is to allow the analysis of ultracool dwarf spectra.

To measure spectral indices, use `measureIndex()` or `measureIndexSet()`:

    sp = splat.getSpectrum(shortname='0415-0935')[0]
    value, error = splat.measureIndex(sp,[1.14,1.165],[1.21,1.235],method='integrate')
    indices = splat.measureIndexSet(sp,set='testi')

The last line returns a dictionary, whose value,error pair can be accessed by the name 
of the index:

    print(indices['sH2O-J'])		# returns value, error

You can also determine the gravity classification of a source following [Allers & Liu (2013)] (http://adsabs.harvard.edu/abs/2013ApJ...772...79A) using `classifyGravity()`:

    sp = splat.getSpectrum(young=True, lucky=True)[0]
    print(splat.classifyGravity(sp))   # returned 'VL-G'

To classify a spectrum, use the various `classifyByXXX` methods:

    sp = splat.getSpectrum(shortname='0415-0935')[0]
    spt,unc = splat.classifyByIndex(sp,set='burgasser')
    spt,unc = splat.classifyByStandard(sp,spt=['T5','T9'])
    result = splat.classifyByTemplate(sp,spt=['T6','T9'],nbest=5)

The last line returns a dictionary containing the best 5 template matches.

To compare a spectrum to another spectrum or a model, use `compareSpectra()`:

    import splat.model as spmod
    mdl = spmod.loadModel(teff=720,logg=4.8,set='btsettl')      # loads a BTSettl08 model 
    sp = splat.getSpectrum(shortname='0415-0935')[0]
    chi,scale = splat.compareSpectra(sp,mdl)
    mdl.scale(scale)
    splat.plotSpectrum(sp,mdl,colors=['k','r'],legend=[sp.name,mdl.name])

You can shortcut the last three lines using the `plot` keyword:

    chi,scale = splat.compareSpectra(sp,mdl,plot=True)


There are also codes *still in development* to fit models directly to spectra: `modelFitGrid()`, `modelFitMCMC()`, and `modelFitEMCEE()`:

    import splat.model as spmod
    sp = splat.getSpectrum(shortname='0415-0935')[0]
    sp.fluxCalibrate('2MASS J',14.49,absolute=True)
    nbest = 5
    result1 = splat.modelFitGrid(sp,set='btsettl')
    result2 = splat.modelFitMCMC(sp,set='btsettl',initial_guess=[800,5.0,0.],nsamples=300,step_sizes=[50.,0.5,0.])
    result3 = splat.modelFitEMCEE(sp,set='btsettl',initial_guess=[800,5.0,0.],nwalkers=12,nsamples=500)

The outputs of all of these fitting functions is a dictionary or list of dictionaries containing the parameters of the best-fitting models; there are also several diagnostic plots produced depending on the routine. View the model fitting page for more details.

All of these routines have many options worth exploring, and which are (increasingly) documented at https://splat.physics.ucsd.edu/splat. If there are capabilities
you need, please suggest them to aburgasser@ucsd.edu, or note it in the "Issues" link on our `github site https://github.com/aburgasser/splat.

## Citing SPLAT and its data


If you use SPLAT tools for your research, please cite Burgasser et al. (2017, ASInC 14, 7), bibcode 2017ASInC..14....7B [NASA ADS] (https://ui.adsabs.harvard.edu/abs/2017ASInC..14....7B/abstract). 

In addition, if you use data contained in SPLAT or the SpeX Prism Library, please be sure to cite the original spectral data source, which can be accessed from the Spectrum object:

    sp = splat.getSpectrum(lucky=True)
    sp.citation().data_reference
> '2016ApJ...817..112S'

    import splat.citations as spcite
    spcite.shortRef(sp.data_reference)
> Schneider, A. C. et al. (2016, Astrophysical Journal, 817, 112)

## Acknowledgements


SPLAT is an collaborative project of research students in the [UCSD Cool Star Lab] (http://www.coolstarlab.org), aimed at developing research through the building of spectral analysis tools.  Contributors to SPLAT have included Christian Aganze, Jessica Birky, Daniella Bardalez Gagliuffi, Adam Burgasser (PI), Caleb Choban, Andrew Davis, Ivanna Escala, Joshua Hazlett, Carolina Herrara Hernandez, Elizabeth Moreno Hilario, Aishwarya Iyer, Yuhui Jin, Mike Lopez, Dorsa Majidi, Diego Octavio Talavera Maya, Alex Mendez, Gretel Mercado, Niana Mohammed, Johnny Parra, Maitrayee Sahi, Adrian Suarez, Melisa Tallis, Tomoki Tamiya, Chris Theissen, and Russell van Linge.

This project has been supported by the National Aeronautics and Space Administration under Grant No. NNX15AI75G.


