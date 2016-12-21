# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
.. note::
         These are the spectrophotometry functions for SPLAT 
"""

# imports - internal
import os

# imports - external
import numpy
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units
from scipy.integrate import trapz        # for numerical integration
from scipy.interpolate import interp1d

# splat functions and constants
from .initialize import *
from .utilities import *



#####################################################
###############   SPECTROPHOTOMETRY   ###############
#####################################################


def filterMag(sp,filt,*args,**kwargs):
    '''
    :Purpose: Determine the photometric magnitude of a source based on its
                spectrum. Spectral fluxes are convolved with the filter profile specified by
                the ``filter`` input.  By default this filter is also
                convolved with a model of Vega to extract Vega magnitudes,
                but the user can also specify AB magnitudes, photon flux or
                energy flux.

    :param sp: Spectrum class object, which should contain wave, flux and
                 noise array elements.
    :type sp: required
    :param filter: String giving name of filter, which can either be one of the predefined filters listed in splat.FILTERS.keys() or a custom filter name
    :type filter: required
    
    :param custom: A 2 x N vector array specifying the wavelengths and transmissions for a custom filter
    :type custom: optional, default = None
    :param notch: A 2 element array that specifies the lower and upper wavelengths for a notch filter (100% transmission within, 0% transmission without)
    :type notch: optional, default = None
    :param vega: compute Vega magnitudes
    :type vega: optional, default = True
    :param ab: compute AB magnitudes
    :type ab: optional, default = False
    :param energy: compute energy flux
    :type energy: optional, default = False
    :param photon: compute photon flux
    :type photon: optional, default = False
    :param filterFolder: folder containing the filter transmission files
    :type filterFolder: optional, default = splat.FILTER_FOLDER
    :param vegaFile: name of file containing Vega flux file, must be within ``filterFolder``
    :type vegaFile: optional, default = vega_kurucz.txt
    :param nsamples: number of samples to use in Monte Carlo error estimation
    :type nsamples: optional, default = 100
    :param info: List the predefined filter names available
    :type info: optional, default = False

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
    >>> sp.fluxCalibrate('2MASS J',14.5)
    >>> splat.filterMag(sp,'MKO J')
        (14.345894376898123, 0.027596454828421831)
    '''
# keyword parameters
    filterFolder = kwargs.get('filterFolder',SPLAT_PATH+FILTER_FOLDER)
    if not os.path.exists(filterFolder):
        filterFolder = SPLAT_URL+FILTER_FOLDER
    vegaFile = kwargs.get('vegaFile','vega_kurucz.txt')
    info = kwargs.get('info',False)
    custom = kwargs.get('custom',False)
    notch = kwargs.get('notch',False)
    vega = kwargs.get('vega',True)
    ab = kwargs.get('ab',False)
    photons = kwargs.get('photons',False)
    photons = kwargs.get('photon',photons)
    energy = kwargs.get('energy',False)
    if (photons or energy or ab):
        vega = False
    nsamples = kwargs.get('nsamples',100)


# check that requested filter is in list
    filter0 = filt
    filt = filt.replace(' ','_')
    filt.upper()
    if (filt not in FILTERS.keys() and isinstance(custom,bool) and isinstance(notch,bool)):
        print('\nFilter '+filt+' not currently available for SPLAT; contact '+SPLAT_EMAIL+'\n')
        info = True

# print out what's available
    if (info):
        filterInfo()
        return numpy.nan, numpy.nan

# Read in filter
    if isinstance(custom,bool) and isinstance(notch,bool):
        fwave,ftrans = numpy.genfromtxt(filterFolder+FILTERS[filt]['file'], comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
# notch filter
    elif isinstance(custom,bool) and isinstance(notch,list):
        d = (notch[1]-notch[0])/1000
        fwave = numpy.arange(notch[0]-5.*d,notch[1]+5.*d,d)
        ftrans = numpy.zeros(len(fwave))
        ftrans[numpy.where(numpy.logical_and(fwave >= notch[0],fwave <= notch[1]))] = 1.
# custom filter
    else:
        fwave,ftrans = custom[0],custom[1]
    fwave = fwave[~numpy.isnan(ftrans)]*u.micron   # temporary fix
    ftrans = ftrans[~numpy.isnan(ftrans)]

# check that spectrum and filter cover the same wavelength ranges
    if numpy.nanmax(fwave) < numpy.nanmin(sp.wave) or numpy.nanmin(fwave) > numpy.nanmax(sp.wave):
        print('\nWarning: no overlap between spectrum for {} and filter {}'.format(sp.name,filter0))
        return numpy.nan, numpy.nan

    if numpy.nanmin(fwave) < numpy.nanmin(sp.wave) or numpy.nanmax(fwave) > numpy.nanmax(sp.wave):
        print('\nWarning: spectrum for {} does not span full filter profile for {}'.format(sp.name,filter0))

# interpolate spectrum onto filter wavelength function
    wgood = numpy.where(~numpy.isnan(sp.noise))
    if len(sp.wave[wgood]) > 0:
        d = interp1d(sp.wave[wgood].value,sp.flux[wgood].value,bounds_error=False,fill_value=0.)
        n = interp1d(sp.wave[wgood].value,sp.noise[wgood].value,bounds_error=False,fill_value=0)
# catch for models
    else:
        print('\nWarning: data values in range of filter {} have no uncertainties'.format(filt))
        d = interp1d(sp.wave.value,sp.flux.value,bounds_error=False,fill_value=0.)
        n = interp1d(sp.wave.value,sp.flux.value*1.e-9,bounds_error=False,fill_value=0.)

    result = []
    if (vega):
# Read in Vega spectrum
        vwave,vflux = numpy.genfromtxt(filterFolder+vegaFile, comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
        vwave = vwave[~numpy.isnan(vflux)]*u.micron
        vflux = vflux[~numpy.isnan(vflux)]*(u.erg/(u.cm**2 * u.s * u.micron))
        vflux.to(sp.funit,equivalencies=u.spectral_density(vwave))
# interpolate Vega onto filter wavelength function
        v = interp1d(vwave.value,vflux.value,bounds_error=False,fill_value=0.)
        val = -2.5*numpy.log10(trapz(ftrans*d(fwave.value),fwave.value)/trapz(ftrans*v(fwave.value),fwave.value))
        for i in numpy.arange(nsamples):
#            result.append(-2.5*numpy.log10(trapz(ftrans*numpy.random.normal(d(fwave),n(fwave))*sp.funit,fwave)/trapz(ftrans*v(fwave)*sp.funit,fwave)))
            result.append(-2.5*numpy.log10(trapz(ftrans*(d(fwave.value)+numpy.random.normal(0,1.)*n(fwave.value)),fwave.value)/trapz(ftrans*v(fwave.value),fwave.value)))
        outunit = 1.

    elif (ab):
        nu = sp.wave.to('Hz',equivalencies=u.spectral())
        fnu = sp.flux.to('Jy',equivalencies=u.spectral_density(sp.wave))
        noisenu = sp.noise.to('Jy',equivalencies=u.spectral_density(sp.wave))
        filtnu = fwave.to('Hz',equivalencies=u.spectral())
        fconst = 3631*u.jansky
        d = interp1d(nu.value,fnu.value,bounds_error=False,fill_value=0.)
        n = interp1d(nu.value,noisenu.value,bounds_error=False,fill_value=0.)
        b = trapz((ftrans/filtnu.value)*fconst.value,filtnu.value)
        val = -2.5*numpy.log10(trapz(ftrans*d(filtnu.value)/filtnu.value,filtnu.value)/b)
        for i in numpy.arange(nsamples):
            a = trapz(ftrans*(d(filtnu.value)+numpy.random.normal(0,1)*n(filtnu.value))/filtnu.value,filtnu.value)
            result.append(-2.5*numpy.log10(a/b))
        outunit = 1.
    elif (energy):
        val = trapz(ftrans*d(fwave.value),fwave.value)
        for i in numpy.arange(nsamples):
            result.append(trapz(ftrans*(d(fwave.value)+numpy.random.normal(0,1.)*n(fwave.value)),fwave.value))
        outunit = u.erg/u.s/u.cm**2
    elif (photons):
        convert = const.h.to('erg s')*const.c.to('micron/s')
        val = trapz(ftrans*fwave.value*convert.value*d(fwave.value),fwave.value)
        for i in numpy.arange(nsamples):
            result.append(trapz(ftrans*fwave.value*convert.value*(d(fwave.value)+numpy.random.normal(0,1.)*n(fwave.value)),fwave.value))
        outunit = 1./u.s/u.cm**2
    else:
        raise NameError('\nfilterMag not given a correct physical quantity (vega, ab, energy, photons) to compute photometry\n\n')


#    val = numpy.nanmean(result)*outunit
    err = numpy.nanstd(result)*outunit
    if len(sp.wave[wgood]) == 0:
        err = 0.
    return val,err


def filterInfo(**kwargs):
    '''
    :Purpose: Prints out the current list of filters in the SPLAT reference library.
    '''

    fname = kwargs.get('filter',sorted(list(FILTERS.keys())))
    for k in fname:
        if k in list(FILTERS.keys()):
        	print('  '+k.replace('_',' ')+': '+FILTERS[k]['description'])
        else:
        	print('  Filter {} not in SPLAT filter list'.format(k))
    return


def filterProperties(filt,**kwargs):
    '''
    :Purpose: Returns a dictionary containing key parameters for a particular filter.

    :param filter: name of filter, must be one of the specifed filters given by splat.FILTERS.keys()
    :type filter: required
    :param verbose: print out information about filter to screen
    :type verbose: optional, default = True

    :Example:
    >>> import splat
    >>> data = splat.filterProperties('2MASS J')
    Filter 2MASS J: 2MASS J-band
    Zeropoint = 1594.0 Jy
    Pivot point: = 1.252 micron
    FWHM = 0.323 micron
    Wavelength range = 1.066 to 1.442 micron
    >>> data = splat.filterProperties('2MASS X')
    Filter 2MASS X not among the available filters:
      2MASS H: 2MASS H-band
      2MASS J: 2MASS J-band
      2MASS KS: 2MASS Ks-band
      BESSEL I: Bessel I-band
      FOURSTAR H: FOURSTAR H-band
      FOURSTAR H LONG: FOURSTAR H long
      FOURSTAR H SHORT: FOURSTAR H short
      ...
    '''
    filt = filt.replace(' ','_')
    filterFolder = kwargs.get('filterFolder',SPLAT_PATH+FILTER_FOLDER)
    if not os.path.exists(filterFolder):
        filterFolder = SPLAT_URL+FILTER_FOLDER

    if (filt not in FILTERS.keys()):
        print('Filter '+filt+' not among the available filters; run filterInfo() to get a list')
#        filterInfo()
        return None
    else:
        report = {}
        report['name'] = filt
        report['description'] = FILTERS[filt]['description']
        report['zeropoint'] = FILTERS[filt]['zeropoint']
        fwave,ftrans = numpy.genfromtxt(filterFolder+FILTERS[filt]['file'], comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
        fw = fwave[numpy.where(ftrans > 0.01*numpy.nanmax(ftrans))]
        ft = ftrans[numpy.where(ftrans > 0.01*numpy.nanmax(ftrans))]
        fw05 = fwave[numpy.where(ftrans > 0.5*numpy.nanmax(ftrans))]
#        print(trapz(ft,fw))
#        print(trapz(fw*ft,fw))
        report['lambda_mean'] = trapz(ft*fw,fw)/trapz(ft,fw)
        report['lambda_pivot'] = numpy.sqrt(trapz(fw*ft,fw)/trapz(ft/fw,fw))
        report['lambda_central'] = 0.5*(numpy.max(fw)+numpy.min(fw))
        report['lambda_fwhm'] = numpy.max(fw05)-numpy.min(fw05)
        report['lambda_min'] = numpy.min(fw)
        report['lambda_max'] = numpy.max(fw)
# report values out
        if kwargs.get('verbose',True):
            print('\nFilter '+filt+': '+report['description'])
            print('Zeropoint = {} Jy'.format(report['zeropoint']))
            print('Pivot point: = {:.3f} micron'.format(report['lambda_pivot']))
            print('FWHM = {:.3f} micron'.format(report['lambda_fwhm']))
            print('Wavelength range = {:.3f} to {:.3f} micron\n'.format(report['lambda_min'],report['lambda_max']))
        return report





