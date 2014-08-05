# SPLAT: The SpeX Prism Library Analysis Toolkit

## Preamble

SPLAT is a python package built upon numpy, scipy, astropy and matplotlib, as well as 
some other common packages.  SPLAT is
designed to interface specifically with the SpeX Prism Library (SPL: http://www.browndwarfs.org/spexprism), 
an online repository of over
1500 low-resolution, near-infrared spectra of low-temperature stars and brown dwarfs.
SPLAT tools allow you to search the library; read in spectra from it; perform basic spectral 
analyses such as classification, index measurement and spectrophotometry; perform
advanced analyses such as spectral model fitting and spectral binary analysis; and 
plot/tabulate your results.  

## Installation and Dependencies

SPLAT is best forked from this github site, which is updated on a semi-regular basis.
SPLAT has not yet reached v1.0, so bugs are rampant. Please help us knock them down by 
sending bug reports to aburgasser@ucsd.edu 

General instructions on setting up to run SPLAT are maintained at http://bit.ly/1AQuy9G

## Using SPLAT

The best place to start is the code documentation, housed at http://bit.ly/1zPZgi2

Here are some examples:

* The best way to read in a spectrum is to use getSpectrum:

```
import splat
splist = splat.getSpectrum(shortname='0415-0935')
splist = splat.getSpectrum(young=True)
splist = splat.getSpectrum(spt=['M7','L5'],jmag=[14.,99.])
```

In each case, splist is a list of Spectrum objects, which is a container of various 
aspects of the spectrum and it source. For example, selecting the first spectrum,

```
sp = splist[0]
```

sp.wave gives the wavelengths of this spectrum, sp.flux the flux values, and sp.noise the 
flux uncertainty.

You can also read in your own spectrum using the loadSpectrum function

```
sp = splat.loadSpectrum(filename='myspectrum.fits',local=True)
```

Note that this file must conform to the standard of the SPL data: the first column is
wavelength in microns, second column flux in f_lambda units, third column (optional) is 
flux uncertainty.

* To display the spectrum, use plotSpectrum

```
splat.plotSpectrum(sp)
```

which will pop up a window displaying flux and noise vs. wavelength. You can save this 
display by adding a filename:

```
splat.plotSpectrum(sp,file='spectrum.png')
```

You can also compare multiple spectra:

```
sp1 = splat.getSpectrum(shortname='0415-0935')[0]
sp2 = splat.getSpectrum(shortname='1217-0311')[1]
splat.plotSpectrum(sp1,sp2,colors=['k','r'])
```

SPLAT can compare an arbitrary number of spectra.

* To measure the indices of a spectrum, use measureIndex or measureIndexSet:

```
sp = splat.getSpectrum(shortname='0415-0935')[0]
value, error = splat.measureIndex(sp,[1.14,1.165],[1.21,1.235],method='integrate')
indices = splat.measureIndexSet(sp,set='burgasser')
```

Note that the latter is a dictionary, whose value,error pair can be accessed by the name 
of the index:

```
print indices['H2O-J']		# returns value, error
```

* To classify a spectrum, use the classifyByXXX methods:

```
sp = splat.getSpectrum(shortname='0415-0935')[0]
spt,unc = splat.classifyByIndex(sp,set='burgasser')
spt,unc = splat.classifyByStandard(sp)
spt,unc = splat.classifyByTemplate(sp)
```

* To compare a spectrum to another spectrum or a model, use compareSpectra:

```
sp = splat.getSpectrum(shortname='0415-0935')[0]
mdl = splat.loadModel(teff=700,logg=5.0)			# currently BTSettl08 only
chi,scale = splat.compareSpectra(sp,mdl)
mdl.scale(scale)
splat.plotSpectrum(sp,mdl,colors=['k','r'])
```

This can be placed in a for loop or MCMC chain to do best-fit parameter determination.

All of these routines have many options worth exploring, and if there are capabilities
you need, please suggest them or contribute code.

## Authors

SPLAT is an experimental, collaborative project of research students in Adam Burgasser's
UCSD Cool Star Lab, aimed at teaching students how to do research by building 
their own analysis tools.  Contributors to SPLAT include Christian Aganze, Daniella Bardalez Gagliuffi,
Adam Burgasser (PI), Caleb Choban, Ivanna Escala, Aishwarya Iyer, Yuhui Jin, Mike Lopez,
Alex Mendez, Johnny Parra, Julian Pilate-Hutcherson, Maitrayee Sahi and Melisa Tallis.



 







.