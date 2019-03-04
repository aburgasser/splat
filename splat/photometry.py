# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
.. note::
         These are the spectrophotometry functions for SPLAT 
"""

# imports - internal
import copy
import os

# imports - external
import numpy
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from scipy.integrate import trapz        # for numerical integration
from scipy.interpolate import interp1d

# splat functions and constants
from .initialize import *
from .utilities import *



#####################################################
###############   SPECTROPHOTOMETRY   ###############
#####################################################

# this function has been obseleted
def checkFilter(filt,verbose=True):
    output = False
    f = copy.deepcopy(filt)
    f = f.replace(' ','_').upper()
    for k in list(FILTERS.keys()):
        if f==k.upper() or f.lower() in FILTERS[k]['altnames']:
            output = k
    if output == False and verbose == True: 
        print('\nFilter '+filt+' not currently available for SPLAT; contact '+EMAIL+'\n')
        filterInfo()
    return output


def filterProfile(filt,**kwargs):
    '''
    :Purpose: Retrieve the filter profile for a SPLAT filter. Returns two arrays: the filter wavelength and filter transmission curve.

    :param filter: String giving the name of one of the predefined filters listed in splat.FILTERS.keys() (required)
    
    :param filterFolder: folder containing the filter transmission files (optional, default = splat.FILTER_FOLDER)

    :Example:
    >>> import splat
    >>> import splat.photometry as spphot
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
    >>> sp.fluxCalibrate('2MASS J',14.5)
    >>> spphot.filterMag(sp,'MKO J')
        (14.345894376898123, 0.027596454828421831)
    '''
# keyword parameters
    filterFolder = kwargs.get('filterFolder',SPLAT_PATH+FILTER_FOLDER)
    if not os.path.exists(filterFolder):
        filterFolder = SPLAT_URL+FILTER_FOLDER

# check that requested filter is in list
    f0 = checkFilterName(filt, verbose=True)
    if f0 == False: raise ValueError
    filt = f0

# read in filter
    fwave,ftrans = numpy.genfromtxt(os.path.normpath(filterFolder+FILTERS[filt]['file']), comments='#', unpack=True, missing_values = ('NaN','nan'), filling_values = (numpy.nan))
#    print(type(fwave),type(ftrans),isinstance(fwave,numpy.ndarray),isinstance(ftrans,numpy.ndarray),not isinstance(fwave,numpy.ndarray) or not isinstance(ftrans,numpy.ndarray))
    if not isinstance(fwave,numpy.ndarray) or not isinstance(ftrans,numpy.ndarray):
        raise ValueError('\nProblem reading in {}'.format(filterFolder+FILTERS[filt]['file']))
    fwave = fwave[~numpy.isnan(ftrans)]*u.micron
    ftrans = ftrans[~numpy.isnan(ftrans)]
    return fwave,ftrans



def filterMag(sp,filt,*args,**kwargs):
    '''
    :Purpose: 

    Determine the photometric magnitude of a source based on its
    spectrum. Spectral fluxes are convolved with the filter profile specified by
    the ``filter`` input.  By default this filter is also
    convolved with a model of Vega to extract Vega magnitudes,
    but the user can also specify AB magnitudes, photon flux or energy flux.

    :Required Parameters:

        **sp**: Spectrum class object, which should contain wave, flux and noise array elements.
        **filter**: String giving name of filter, which can either be one of the predefined filters listed in splat.FILTERS.keys() or a custom filter name

    :Optional Parameters:
    
        **custom** = None: A 2 x N vector array specifying the wavelengths and transmissions for a custom filter
        **notch** = None: A 2 element array that specifies the lower and upper wavelengths for a notch filter (100% transmission within, 0% transmission without)
        **vega** = True: compute Vega magnitudes (may be set by filter)
        **ab** = False: compute AB magnitudes (may be set by filter)
        **energy** = False: compute energy flux
        **photon** = False: compute photon flux
        **filterFolder** = splat.FILTER_FOLDER: folder containing the filter transmission files
        **vegaFile** = 'vega_kurucz.txt': name of file containing Vega flux file, must be within ``filterFolder``
        **nsamples** = 100: number of samples to use in Monte Carlo error estimation
        **info** = False: List the predefined filter names available
        **verbose** = True: List the predefined filter names available

    :Example:
    >>> import splat
    >>> import splat.photometry as spphot
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
    >>> sp.fluxCalibrate('2MASS J',14.5)
    >>> spphot.filterMag(sp,'MKO J')
        (14.345894376898123, 0.027596454828421831)
    '''
# keyword parameters
    filterFolder = kwargs.get('filterFolder',SPLAT_PATH+FILTER_FOLDER)
    if not os.path.exists(filterFolder):
        filterFolder = SPLAT_URL+FILTER_FOLDER
    vegaFile = kwargs.get('vegaFile',VEGAFILE)
    info = kwargs.get('info',False)
    custom = kwargs.get('custom',False)
    notch = kwargs.get('notch',False)
    vega = kwargs.get('vega',True)
    ab = kwargs.get('ab',not vega)
    rsr = kwargs.get('rsr',False)
    nsamples = kwargs.get('nsamples',100)
    verbose = kwargs.get('verbose',False)


# check that requested filter is in list
    if isinstance(custom,bool) and isinstance(notch,bool):
        f0 = checkFilterName(filt,verbose=True)
        if f0 == False: 
            return numpy.nan, numpy.nan
        filt = f0
        
# reset filter calculation methods based on filter design
        if 'ab' in FILTERS[filt]['method']: 
            ab = kwargs.get('ab',True)
            vega = not ab
        if 'vega' in FILTERS[filt]['method']: 
            vega = kwargs.get('vega',True)
            ab = not vega
        rsr = FILTERS[filt]['rsr']

# other possibilities
    photons = kwargs.get('photons',False)
    photons = kwargs.get('photon',photons)
    energy = kwargs.get('energy',False)
    energy = kwargs.get('flux',energy)
    if (photons or energy):
        vega = False
        ab = False
    if photons: energy = False
    if energy: photons = False


# Read in filter
    if isinstance(custom,bool) and isinstance(notch,bool):
        fwave,ftrans = filterProfile(filt,**kwargs)
# notch filter
    elif isinstance(custom,bool) and isinstance(notch,list):
        dn = (notch[1]-notch[0])/1000
        fwave = numpy.arange(notch[0]-5.*dn,notch[1]+5.*dn,dn)
        ftrans = numpy.zeros(len(fwave))
        ftrans[numpy.where(numpy.logical_and(fwave >= notch[0],fwave <= notch[1]))] = 1.
# custom filter
    else:
        fwave,ftrans = custom[0],custom[1]

# units
    if isinstance(fwave,u.quantity.Quantity) == True:
        fwave = fwave.to(u.micron)
    else:
        fwave = fwave*u.micron

# check that spectrum and filter cover the same wavelength ranges
    if numpy.nanmax(fwave) < numpy.nanmin(sp.wave) or numpy.nanmin(fwave) > numpy.nanmax(sp.wave):
        if verbose==True: print('\nWarning: no overlap between spectrum for {} and filter {}'.format(sp.name,filt))
        return numpy.nan, numpy.nan

    if numpy.nanmin(fwave) < numpy.nanmin(sp.wave) or numpy.nanmax(fwave) > numpy.nanmax(sp.wave):
        if verbose==True: print('\nWarning: spectrum for {} does not span full filter profile for {}'.format(sp.name,filt))

# interpolate spectrum onto filter wavelength function
    wgood = numpy.where(~numpy.isnan(sp.noise))
    if len(sp.wave[wgood]) > 0:
        d = interp1d(sp.wave[wgood].value,sp.flux[wgood].value,bounds_error=False,fill_value=0.)
        n = interp1d(sp.wave[wgood].value,sp.noise[wgood].value,bounds_error=False,fill_value=0)
# catch for models
    else:
        if verbose==True: print('\nWarning: data values in range of filter {} have no uncertainties'.format(filt))
        d = interp1d(sp.wave.value,sp.flux.value,bounds_error=False,fill_value=0.)
        n = interp1d(sp.wave.value,sp.flux.value*1.e-9,bounds_error=False,fill_value=0.)

    result = []
    if (vega):
# Read in Vega spectrum
        vwave,vflux = numpy.genfromtxt(os.path.normpath(filterFolder+vegaFile), comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
        vwave = vwave[~numpy.isnan(vflux)]*u.micron
        vflux = vflux[~numpy.isnan(vflux)]*(u.erg/(u.cm**2 * u.s * u.micron))
        vflux.to(sp.funit,equivalencies=u.spectral_density(vwave))
# interpolate Vega onto filter wavelength function
        v = interp1d(vwave.value,vflux.value,bounds_error=False,fill_value=0.)
        if rsr:
            val = -2.5*numpy.log10(trapz(ftrans*fwave.value*d(fwave.value),fwave.value)/trapz(ftrans*fwave.value*v(fwave.value),fwave.value))
        else:
            val = -2.5*numpy.log10(trapz(ftrans*d(fwave.value),fwave.value)/trapz(ftrans*v(fwave.value),fwave.value))
        for i in numpy.arange(nsamples):
#            result.append(-2.5*numpy.log10(trapz(ftrans*numpy.random.normal(d(fwave),n(fwave))*sp.funit,fwave)/trapz(ftrans*v(fwave)*sp.funit,fwave)))
            if rsr:
                result.append(-2.5*numpy.log10(trapz(ftrans*fwave.value*(d(fwave.value)+numpy.random.normal(0,1.)*n(fwave.value)),fwave.value)/trapz(ftrans*fwave.value*v(fwave.value),fwave.value)))
            else:
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
        outunit = u.erg/u.s/u.cm**2
        if rsr:
            a = trapz(ftrans*fwave.value*d(fwave.value),fwave.value)*sp.wave.unit*sp.flux.unit
            b = trapz(ftrans*fwave.value,fwave.value)*sp.wave.unit
            c = trapz(ftrans*fwave.value*fwave.value,fwave.value)*sp.wave.unit*sp.wave.unit
            val = (a/b * c/b).to(outunit).value
        else:
            a = trapz(ftrans*d(fwave.value),fwave.value)*sp.wave.unit*sp.flux.unit
            b = trapz(ftrans,fwave.value)*sp.wave.unit
            c = trapz(ftrans*fwave.value,fwave.value)*sp.wave.unit*sp.wave.unit
            val = (a/b * c/b).to(outunit).value
        for i in numpy.arange(nsamples):
            if rsr:
                result.append((trapz(ftrans*fwave.value*(d(fwave.value)+numpy.random.normal(0,1.)*n(fwave.value)),fwave.value)*sp.wave.unit*sp.flux.unit).to(outunit).value)
            else:
                result.append((trapz(ftrans*(d(fwave.value)+numpy.random.normal(0,1.)*n(fwave.value)),fwave.value)*sp.wave.unit*sp.flux.unit).to(outunit).value)

    elif (photons):
        outunit = 1./u.s/u.cm**2
        convert = const.h.to('erg s')*const.c.to('micron/s')
        val = (trapz(ftrans*fwave.value*convert.value*d(fwave.value),fwave.value)*sp.wave.unit*sp.flux.unit*convert.unit).to(outunit).value
        for i in numpy.arange(nsamples):
            result.append((trapz(ftrans*fwave.value*convert.value*(d(fwave.value)+numpy.random.normal(0,1.)*n(fwave.value)),fwave.value)*sp.wave.unit*sp.flux.unit*convert.unit).to(outunit).value)
    else:
        raise NameError('\nfilterMag not given a correct physical quantity (vega, ab, energy, photons) to compute photometry\n\n')


#    val = numpy.nanmean(result)*outunit
    err = numpy.nanstd(result)
    if len(sp.wave[wgood]) == 0:
        err = 0.
    return val*outunit,err*outunit


def filterInfo(*args,**kwargs):
    '''
    :Purpose: Prints out the current list of filters in the SPLAT reference library.
    '''

    verbose = kwargs.get('verbose',True)

    if len(args) > 0: 
        fname = list(args)
    elif kwargs.get('filter',False) != False: 
        fname = kwargs['filter']
    else: 
        fname = sorted(list(FILTERS.keys()))
    if isinstance(fname,list) == False: 
        fname = [fname]

    output = {}
    for k in fname:
        f = checkFilterName(k)
        if f != False:
            output[f] = {}
            output[f]['description'] = FILTERS[f]['description']
            output[f]['zeropoint'] = FILTERS[f]['zeropoint']
            fwave,ftrans = filterProfile(f,**kwargs)
            try:
                fwave = fwave.to(u.micron)
            except:
                fwave = fwave*u.micron
            fw = fwave[numpy.where(ftrans > 0.01*numpy.nanmax(ftrans))]
            ft = ftrans[numpy.where(ftrans > 0.01*numpy.nanmax(ftrans))]
            fw05 = fwave[numpy.where(ftrans > 0.5*numpy.nanmax(ftrans))]
            output[f]['lambda_mean'] = trapz(ft*fw,fw)/trapz(ft,fw)
            output[f]['lambda_pivot'] = numpy.sqrt(trapz(fw*ft,fw)/trapz(ft/fw,fw))
            output[f]['lambda_central'] = 0.5*(numpy.max(fw)+numpy.min(fw))
            output[f]['lambda_fwhm'] = numpy.max(fw05)-numpy.min(fw05)
            output[f]['lambda_min'] = numpy.min(fw)
            output[f]['lambda_max'] = numpy.max(fw)
            if verbose ==True: 
                print(f.replace('_',' ')+': '+output[f]['zeropoint'])
                print('Zeropoint = {} Jy'.format(output[f]['zeropoint']))
                print('Central wavelength: = {:.3f}'.format(output[f]['lambda_central']))
                print('Mean wavelength: = {:.3f}'.format(output[f]['lambda_mean']))
                print('Pivot point: = {:.3f}'.format(output[f]['lambda_pivot']))
                print('FWHM = {:.3f}'.format(output[f]['lambda_fwhm']))
                print('Wavelength range = {:.3f} to {:.3f}\n'.format(output[f]['lambda_min'],output[f]['lambda_max']))             
        else:
        	if verbose ==True: print('  Filter {} not in SPLAT filter list'.format(k))
    kys = list(output.keys())
    if len(kys) == 1: return output[kys[0]]
    else: return output


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
    filterFolder = kwargs.get('filterFolder',SPLAT_PATH+FILTER_FOLDER)
    if not os.path.exists(filterFolder):
        filterFolder = SPLAT_URL+FILTER_FOLDER

# check that requested filter is in list
    filt = checkFilterName(filt)
    if filt == False: return None

    report = {}
    report['name'] = filt
    report['description'] = FILTERS[filt]['description']
    report['zeropoint'] = FILTERS[filt]['zeropoint']
    report['method'] = FILTERS[filt]['method']
    report['rsr'] = FILTERS[filt]['rsr']
    fwave,ftrans = filterProfile(filt,**kwargs)
    try:
        fwave = fwave.to(u.micron)
    except:
        fwave = fwave*u.micron
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
    report['wave'] = fwave
    report['transmission'] = ftrans
# report values out
    if kwargs.get('verbose',False):
        print('\nFilter '+filt+': '+report['description'])
        print('Zeropoint = {} Jy'.format(report['zeropoint']))
        print('Pivot point: = {:.3f}'.format(report['lambda_pivot']))
        print('FWHM = {:.3f}'.format(report['lambda_fwhm']))
        print('Wavelength range = {:.3f} to {:.3f}\n'.format(report['lambda_min'],report['lambda_max']))
    return report


def magToFlux(mag,filt,**kwargs):
    '''
    :Purpose: Converts a magnitude into an energy, and vice versa.

    :param mag: magnitude on whatever system is defined for the filter or provided (required)
    :param filter: name of filter, must be one of the specifed filters given by splat.FILTERS.keys() (required)
    :param reverse: convert energy into magnitude instead (optional, default = False)
    :param ab: magnitude is on the AB system (optional, default = filter preference)
    :param vega: magnitude is on the Vega system (optional, default = filter preference)
    :param rsr: magnitude is on the Vega system (optional, default = filter preference)
    :param units: units for energy as an astropy.units variable; if this conversion does not work, the conversion is ignored (optional, default = erg/cm2/s)
    :param verbose: print out information about filter to screen (optional, default = False)

    WARNING: THIS CODE IS ONLY PARTIALLY COMPLETE
    '''

# keyword parameters
    filterFolder = kwargs.get('filterFolder',SPLAT_PATH+FILTER_FOLDER)
    if not os.path.exists(filterFolder):
        filterFolder = SPLAT_URL+FILTER_FOLDER
    vegaFile = kwargs.get('vegaFile','vega_kurucz.txt')
    vega = kwargs.get('vega',True)
    ab = kwargs.get('ab',not vega)
    rsr = kwargs.get('rsr',False)
    nsamples = kwargs.get('nsamples',100)
    custom = kwargs.get('custom',False)
    notch = kwargs.get('notch',False)
    base_unit = u.erg/(u.cm**2 * u.s)
    return_unit = kwargs.get('unit',base_unit)
    e_mag = kwargs.get('uncertainty',0.)
    e_mag = kwargs.get('unc',e_mag)
    e_mag = kwargs.get('e_mag',e_mag)
    if not isinstance(mag,u.quantity.Quantity): mag=mag*u.s/u.s
    if not isinstance(e_mag,u.quantity.Quantity): e_mag=e_mag*mag.unit

# check that requested filter is in list
    filt = checkFilterName(filt)
    if filt == False: return numpy.nan, numpy.nan

# reset filter calculation methods based on filter design
    if 'ab' in FILTERS[filt]['method']: 
        ab = kwargs.get('ab',True)
        vega = not ab
    if 'vega' in FILTERS[filt]['method']: 
        vega = kwargs.get('vega',True)
        ab = not vega
    if 'rsr' in FILTERS[filt]['method']: 
        rsr = kwargs.get('rsr',True)


# Read in filter
    if isinstance(custom,bool) and isinstance(notch,bool):
        fwave,ftrans = filterProfile(filt,**kwargs)
# notch filter
    elif isinstance(custom,bool) and isinstance(notch,list):
        dn = (notch[1]-notch[0])/1000
        fwave = numpy.arange(notch[0]-5.*dn,notch[1]+5.*dn,dn)*u.micron
        ftrans = numpy.zeros(len(fwave))
        ftrans[numpy.where(numpy.logical_and(fwave >= notch[0],fwave <= notch[1]))] = 1.
# custom filter
    else:
        fwave,ftrans = custom[0],custom[1]
    if isinstance(fwave,u.quantity.Quantity) == False: fwave=fwave*u.micron
    if isinstance(ftrans,u.quantity.Quantity) == True: ftrans=ftrans.value
    fwave = fwave[~numpy.isnan(ftrans)]
    ftrans = ftrans[~numpy.isnan(ftrans)]

    result = []
    err = 0.
# magnitude -> energy
    if kwargs.get('reverse',False) == False:
        
        if vega == True:
    # Read in Vega spectrum
            vwave,vflux = numpy.genfromtxt(os.path.normpath(filterFolder+vegaFile), comments='#', unpack=True, \
                missing_values = ('NaN','nan'), filling_values = (numpy.nan))
            vwave = vwave[~numpy.isnan(vflux)]*u.micron
            vflux = vflux[~numpy.isnan(vflux)]*(u.erg/(u.cm**2 * u.s * u.micron))
    # interpolate Vega onto filter wavelength function
            v = interp1d(vwave.value,vflux.value,bounds_error=False,fill_value=0.)
            if rsr: fact = trapz(ftrans*fwave.value*v(fwave.value),fwave.value)
            else: fact = trapz(ftrans*v(fwave.value),fwave.value)
            val = 10.**(-0.4*mag.value)*fact*u.erg/(u.cm**2 * u.s)
    # calculate uncertainty        
            if e_mag.value > 0.:
                for i in numpy.arange(nsamples): result.append(10.**(-0.4*(mag.value+numpy.random.normal(0,1.)*e_mag.value))*fact)
                err = (numpy.nanstd(result))*u.erg/(u.cm**2 * u.s)
            else: err = 0.*u.erg/(u.cm**2 * u.s)
        elif ab == True:
            fconst = 3631*u.jansky
            ftrans = (ftrans*fconst).to(u.erg/(u.cm**2 * u.s * u.micron),equivalencies=u.spectral_density(fwave))
            if rsr: fact = trapz(ftrans.value*fwave.value,fwave.value)
            else: fact = trapz(ftrans.value,fwave.value)
            val = (10.**(-0.4*mag.value)*fact)*u.erg/(u.cm**2 * u.s)
    # calculate uncertainty        
            if e_mag.value > 0.:
                for i in numpy.arange(nsamples): result.append(10.**(-0.4*(mag.value+numpy.random.normal(0,1.)*e_mag.value))*fact)
                err = (numpy.nanstd(result))*u.erg/(u.cm**2 * u.s)
            else: err = 0.*u.erg/(u.cm**2 * u.s)
        else:
            raise ValueError('\nmagToFlux needs vega or ab method specified')

# convert to desired energy units
#        try:
        val.to(return_unit)
        err.to(return_unit)
#        except:
#            print('\nWarning: unit {} is not an energy flux unit'.format(return_unit))
        try:
            val.to(base_unit)
            err.to(base_unit)
        except:
            print('\nWarning: cannot convert result to an energy flux unit'.format(base_unit))
            return numpy.nan, numpy.nan
        return val, err

# energy -> magnitude
# THIS NEEDS TO BE COMPLETED
    else:
        print('passed')
        pass
# check that input is an energy flux
#        try:
#            mag.to(base_unit)        
#            e_mag.to(base_unit)        
#        except:
#            raise ValueError('\nInput quantity unit {} is not a flux unit'.format(mag.unit))


def visualizeFilter(filters,verbose=True,xra=[],yra=[0,1.2],**kwargs):
    '''
    :Purpose: Plots a filter profile or set of filter profiles, optionally on top of a spectrum

    WARNING: THIS CODE IS CURRENTLY UNDER DEVELOPMENT, BUGS MAY BE COMMON

    '''
    filt = copy.deepcopy(filters)
    wunit = kwargs.get('wunit',DEFAULT_WAVE_UNIT)

# single filter name  
    if isinstance(filt,str):
        filt = [filt]

    if isinstance(filt,list):

# list of filter names
        if isinstance(filt[0],str):
            for f in filt:
                fc = checkFilterName(f)
                filt.remove(f)
                if fc == False: 
                    if verbose==True: print('Removed filter {}: not included in SPLAT'.format(f))
                else:
                    filt.insert(len(filt),fc)
            if len(filt) == 0:
                raise ValueError('Did not recognize any of the input filters {}'.format(filters))

# prep parameters
            fwave,ftrans = filterProfile(f,**kwargs)
            if isUnit(fwave): wunit = kwargs.get('wunit',fwave.unit)

            xl = kwargs.get('xlabel','Wavelength ({})'.format(wunit))
            yl = kwargs.get('ylabel','Transmission Curve')
            legend = []
            fig = plt.figure(figsize=kwargs.get('figsize',[5,4]))
            for i,f in enumerate(filt):
                fwave,ftrans = filterProfile(f,**kwargs)
                if isUnit(fwave): fwave.to(wunit)
                else: fwave = fwave*wunit
                if kwargs.get('normalize',False): ftrans = ftrans/numpy.nanmax(ftrans)
                plt.plot(fwave,ftrans)
                if len(xra) == 0: xra = [numpy.nanmin(fwave.value),numpy.nanmax(fwave.value)]
                xra = [numpy.nanmin([xra[0],numpy.nanmin(fwave.value)]),numpy.nanmax([xra[1],numpy.nanmax(fwave.value)])]
                yra = [yra[0],numpy.nanmax([yra[1],numpy.nanmax(ftrans)])]
                legend.append(FILTERS[f]['description'])
                if FILTERS[f]['rsr'] == True: yl = kwargs.get('ylabel','Transmission Curve')

# list of notch ranges
        if isinstance(filt[0],int) or isinstance(filt[0],float):
            filt = [filt]

# list of notch ranges
        if isinstance(filt[0],list):
            xl = kwargs.get('xlabel','Wavelength ({})'.format(wunit))
            yl = kwargs.get('ylabel','Transmission Curve')
            legend = []
            fig = plt.figure(figsize=kwargs.get('figsize',[5,4]))
            for i,f in enumerate(filt):
                fwave,ftrans = numpy.linspace(f[0],f[1],1000)*wunit,numpy.ones(1000)
                plt.plot(fwave,ftrans)
                if len(xra) == 0: xra = [numpy.nanmin(fwave.value),numpy.nanmax(fwave.value)]
                xra = [numpy.nanmin([xra[0],numpy.nanmin(fwave.value)]),numpy.nanmax([xra[1],numpy.nanmax(fwave.value)])]
                yra = [yra[0],numpy.nanmax([yra[1],numpy.nanmax(ftrans)])]
                legend.append('Filter {}'.format(i+1))

    else:
        raise ValueError('Could not parse input {}'.format(filt))
        

# add a comparison spectrum
    sp = kwargs.get('spectrum',None)
    sp = kwargs.get('comparison',sp)

    if isinstance(sp,splat.core.Spectrum) == True:
        print(xra)
        sp.normalize(xra)
        sp.scale(numpy.nanmax(ftrans)*kwargs.get('comparison_scale',0.8))
        plt.plot(sp.wave,sp.flux,color=kwargs.get('comparison_color','k'),alpha=kwargs.get('comparison_alpha',0.5))
        legend.append(sp.name)
        yra = [yra[0],yra[1]*1.1]

# finish up
    plt.xlim(xra)
    plt.ylim(yra)
    plt.xlabel(xl)
    plt.ylabel(yl)
    plt.legend(legend)

# save if desired
    file = kwargs.get('file','')
    file = kwargs.get('filename',file)
    file = kwargs.get('output',file)
    if file != '': plt.savefig(file)
    
    return fig


#########################################
########   SED FITTING TOOLS    #########
### WARNING: THESE ARE EXPERIMENTAL!! ###
#########################################

# plan:

def modelMagnitudes(verbose=True):
    '''
    this will be a code that calculates a set of magnitudes for a model set's SED models
    saves to file that could be uploaded
    pre-save some model magnitudes
    '''
    pass


def interpolateMagnitudes(verbose=True):
    '''
    produces an interpolated value for a grid set of model magnitudes
    '''
    pass

def compareMagnitudes(mags1,mags2,unc=None,unc2=None,ignore=[],verbose=True):
    '''
    this code compares a set of magnitudes using one of several statistics
    '''
    chi = 0.
    dm,em = [],[]
    for f in list(mags1.keys()):
        if f in list(mags2.keys()) and f in list(unc.keys()) and f not in ignore: 
            dm.append(mags1[f]-mags2[f])
            em.append(unc[f])
# find best scale factor
    dm = numpy.array(dm)
    em = numpy.array(em)
    offset = numpy.sum(dm/em**2)/numpy.sum (1./em**2)
    dmo = numpy.array([m-offset for m in dm])
    return numpy.sum((dmo/em)**2), offset

def SEDFitGrid(verbose=True):
    '''
    this code will compare a set of magnitudes to a grid of model magnitudes and choose the
    closest match based on various statistics
    '''
    pass

def SEDFitMCMC(verbose=True):
    '''
    this code will conduct a comparison of a set of magnitudes to model magnitudes using an
    MCMC wrapper, and choose best/average/distribution of parameters
    '''
    pass

def SEDFitAmoeba(verbose=True):
    '''
    this code will conduct a comparison of a set of magnitudes to model magnitudes using an
    Amoeba (Nelder-Mead) wrapper, and choose the closest match
    '''
    pass

def SEDVisualize(*args,e_mag=None,spectrum=None,file='',verbose=True):
    '''
    Visualizes magnitudes on SED scale (flux = lam x F_lam), with option of also comparing to SED spectrum
    '''




#####################################################
###############    MAGNITUDE CLASS    ###############
#####################################################

class Magnitude(object):
    '''
    :Description: 

        This is a class data structure for a magnitude value

    '''

    def __init__(self, magnitude, filt, uncertainty=0., magtype='apparent', verbose=False,**kwargs):
        self.magnitude = magnitude
        self.uncertainty = uncertainty
        self.type = magtype

# check filter and rename if necessary
        self.knownfilter = True
        fflag = checkFilterName(filt,verbose=verbose)
        if fflag == False:
            if verbose== True: print('filter {} is not a standard filter; some functions may not work'.format(filt))
            self.knownfilter = False
        else: filt = fflag
        self.filter = filt

# some things that are based on presets
        if self.knownfilter == True:
            self.wave,self.transmission = filterProfile(self.filter)
            info = filterProperties(self.filter)
            for k in info.keys(): setattr(self,k,info[k])

    def __copy__(self):
        '''
        :Purpose: Make a copy of a Magnitude object
        '''
        s = type(self)(self.magnitude,self.filter,uncertainty=self.uncertainty)
        s.__dict__.update(self.__dict__)
        return s

# backup version
    def copy(self):
        '''
        :Purpose: Make a copy of a Magnitude object
        '''
        s = type(self)(self.magnitude,self.filter,uncertainty=self.uncertainty)
        s.__dict__.update(self.__dict__)
        return s

    def __repr__(self):
        '''
        :Purpose: A simple representation of the Spectrum object
        '''
        if self.uncertainty != 0. and numpy.isfinite(self.uncertainty):
            return '{} magnitude of {}+/-{}'.format(self.filter,self.magnitude,self.uncertainty)
        else:
            return '{} magnitude of {}'.format(self.filter,self.magnitude)

    def __add__(self,other,samp=1000):
        '''
        :Purpose: 

            A representation of addition for Magnitude classes that takes into account uncertainties

        :Output: 

            A new Magnitude object equal to the sum of values
        '''

# make a copy and fill in combined magnitude
        out = copy.deepcopy(self)
        out.magnitude = self.magnitude+other.magnitude
        out.uncertainty = self.uncertainty+other.uncertainty

# combine noises
        if self.uncertainty != 0 and other.uncertainty != 0:
            m1 = numpy.random.normal(self.magnitude,self.uncertainty,samp)
            m2 = numpy.random.normal(other.magnitude,other.uncertainty,samp)
            val = m1+m2
            out.uncertainty = numpy.nanstd(val)

# check filter agreement
        if self.filter != other.filter:
            out.filter = '{}+{}'.format(self.filter,other.filter)

        return out

    def __sub__(self,other,samp=1000):
        '''
        :Purpose: 

            A representation of subtraction for Magnitude classes that takes into account uncertainties

        :Output: 

            A new Magnitude object equal to the diffence of values
        '''

# make a copy and fill in combined magnitude
        out = copy.deepcopy(self)
        out.magnitude = self.magnitude-other.magnitude
        out.uncertainty = self.uncertainty+other.uncertainty

# combine noises
        if self.uncertainty != 0 and other.uncertainty != 0:
            m1 = numpy.random.normal(self.magnitude,self.uncertainty,samp)
            m2 = numpy.random.normal(other.magnitude,other.uncertainty,samp)
            val = m1-m2
            out.uncertainty = numpy.nanstd(val)

# check filter agreement
        if self.filter != other.filter:
            out.filter = '{}-{}'.format(self.filter,other.filter)

        return out

    def flux(self,type='fnu',samp=1000):
        '''
        :Purpose: 

            Report the equivalent flux density of a magnitude

        :Output: 

            astropy quantity in flux density units (default = erg/cm2/s/micron)

        '''
        pass


    def addFlux(self,other,samp=1000):
        '''
        :Purpose: 

            A representation of addition for magnitudes (addition of fluxes)

        :Output: 

            A new magnitude object equal to the equivalent sum of fluxes

        '''
# check filter agreement
        if self.filter != other.filter:
            raise ValueError('magnitudes filters {} and {} are not the same'.format(self.filter,other.filter))

# make a copy and fill in combined magnitude
        out = copy.deepcopy(self)
        out.magnitude = self.magnitude-2.5*numpy.log10(1.+10.**(-0.4*(other.magnitude-self.magnitude)))
        out.uncertainty = self.uncertainty+other.uncertainty

# combine noises
        if self.uncertainty != 0 and other.uncertainty != 0:
            m1 = numpy.random.normal(self.magnitude,self.uncertainty,samp)
            m2 = numpy.random.normal(other.magnitude,other.uncertainty,samp)
            val = m1-2.5*numpy.log10(1.+10.**(-0.4*(m2-m1)))
            out.uncertainty = numpy.nanstd(val)

        return out


