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
from scipy.integrate import trapezoid as trapz        # for numerical integration
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
        if f==k.upper() or f.lower() in FILTERS[k]['altname']:
            output = k
    if output == False and verbose == True: 
        print('\nFilter {} not currently available for SPLAT\n\nTry {}'.format(filt,list(FILTERS.keys())))
        filterInfo()
    return output


def filterProfile(filt,filterfolder=FILTER_FOLDER,plot=False,file='',verbose=False,**kwargs):
    '''
    Purpose
    -------

    Retrieve the filter profile for a SPLAT filter. Returns two arrays: the filter wavelength and filter transmission curve.

    Parameters
    ----------

    filter : str
        String giving the name of one of the predefined filters listed in splat.FILTERS.keys()
    
    filterfolder : str, default = splat.FILTER_FOLDER
        Optional filepath to folder containing the filter transmission files

    plot : bool, default = False
        Set to True to plot the filter profile (using `visualizeFilter()_`)

    file : str, default = ''
        Full filepath for plot output; if empty string, no plot created

    verbose : bool, default = False
        Set to True to receive verbose feedback

    additional plotting options contained in matplotlib.pyplot can be passed as well, 
    including color, label, xlim, ylim, xlabel, ylabel, figsize, normalize, and fill

    Outputs
    -------

    wave,trans : numpy.array
        float arrays containing the wavelength and transmission profile

    Example
    -------
    
    >>> import splat
    >>> import splat.photometry as spphot
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
    >>> sp.fluxCalibrate('2MASS J',14.5)
    >>> spphot.filterMag(sp,'MKO J')
        (14.345894376898123, 0.027596454828421831)
 
    Dependencies
    ------------
    matplotlib.pyplot
    numpy
    os
    `checkFilterName()`_
    `isUnit()`_
    `visualizeFilter()`_

    .. _`checkFilterName()` : api.html#splat.photometry.checkFilterName
    .. _`isUnit()` : api.html#splat.utilities.isUnit
    .. _`visualizeFilter()` : api.html#splat.photometry.visualizeFilter
    '''
# keyword parameters
    if not os.path.exists(filterfolder): filterfolder = FILTER_FOLDER

# check that requested filter is in list
    f0 = checkFilterName(filt, verbose=True)
    if f0 == False: raise ValueError('Cannot find filter {} among default filters'.format(filt))
    filt = f0
    if verbose==True: print('Reading in filter profile for filter {}'.format(filt))

# read in filter
    fwave,ftrans = numpy.genfromtxt(os.path.normpath(filterfolder+FILTERS[filt]['file']), comments='#', unpack=True, missing_values = ('NaN','nan'), filling_values = (numpy.nan))
    if not isinstance(fwave,numpy.ndarray) or not isinstance(ftrans,numpy.ndarray):
        raise ValueError('\nProblem reading in {}'.format(filterfolder+FILTERS[filt]['file']))
    fwave = fwave[~numpy.isnan(ftrans)]*u.micron
    ftrans = ftrans[~numpy.isnan(ftrans)]

# optional plot
    if plot==True:
        visualizeFilter(filt,file=file,**kwargs)
        # if isUnit(fwave): x = fwave.value
        # else: x = fwave
        # if isUnit(ftrans): y = ftrans.value
        # else: y = ftrans
        # fig = plt.figure(figsize=kwargs.get('figsize',[5,4]))
        # if kwargs.get('normalize',False)==True: y = y/numpy.nanmax(y)
        # plt.plot(x,y,label=kwargs.get('label',filt.replace('_',' ')),color=kwargs.get('color','k'))
        # if kwargs.get('fill',False)==True: plt.fill_between(x,y,numpy.zeros(len(y)),color=kwargs.get('color','k'),alpha=0.1)
        # plt.xlim([numpy.nanmin(fwave.value),numpy.nanmax(fwave.value)])
        # xlim = kwargs.get('xlim',[])
        # if len(xlim)==2: plt.xlim(xlim)
        # ylim = kwargs.get('ylim',[])
        # if len(ylim)==2: plt.ylim(ylim)
        # plt.xlabel(kwargs.get('xlabel','Wavelength'),fontsize=14)
        # plt.ylabel(kwargs.get('ylabel','Transmission'),fontsize=14)
        # plt.xticks(fontsize=12)
        # plt.yticks(fontsize=12)
        # plt.legend(fontsize=14)
        # plt.tight_layout()
        # if file!='': 
        #     plt.savefig(file)
        #     if verbose==True: print('Plot file saved to {}'.format(file))
        # plt.show()

# output
    return fwave,ftrans



def filterMag(sp,filt,computevalue='vega',nsamples=100,custom=False,notch=False,rsr=False,vegafile=VEGAFILE,filterfolder=FILTER_FOLDER,info=False,verbose=False,**kwargs):
    '''
    Purpose 
    -------

    Determine the photometric magnitude of a source based on its
    spectrum. Spectral fluxes are convolved with the filter profile specified by
    the ``filter`` input.  By default this filter is also
    convolved with a model of Vega to extract Vega magnitudes,
    but the user can also specify AB magnitudes, photon flux or energy flux.

    Parameters
    ----------

    sp : Spectrum class object
        Spectrum to compute spectrophotometry, which should contain wave, flux and noise array elements.

    filter : str
        String giving name of filter, which can either be one of the predefined filters listed in splat.FILTERS.keys() or a custom filter name
    
    computevalue : str, default = 'vega'
        Type of value to compute; must be one of vega, ab, energy, photons; can also be set with vega, ab, energy, and photons keywords

    custom : array, default = None
        A 2 x N vector array specifying the wavelengths and transmissions for a custom filter

    notch : array, default = None: 
        A 2 element array that specifies the lower and upper wavelengths for a notch filter (100% transmission within, 0% transmission without)
    
    filterfolder : str, default = splat.FILTER_FOLDER
        Optional filepath to folder containing the filter transmission files

    vegafile :  str, default = 'vega_kurucz.txt'
        Filename containing Vega flux file, either full path or assumed to be within ``filterfolder``
    
    nsamples : int, default = 100
        Number of samples to use in Monte Carlo error estimation
    
    info : bool default = False: 
        Set to True to list the predefined filter names available
    
    verbose : bool, default = False
        Set to True to receive verbose feedback

    Outputs
    -------

    mag, err : float
        Float value of magnitude/energy/photons and its uncertainty

    Example
    -------
    
    >>> import splat
    >>> import splat.photometry as spphot
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
    >>> sp.fluxCalibrate('MKO J',14.5)
    >>> spphot.filterMag(sp,'MKO J')
        (14.499999979756499, 0.029712572838966562)
 
    Dependencies
    ------------
    astropy.units
    matplotlib.pyplot
    numpy
    os
    scipy.interpolate.interp1d
    scipy.integrate.trapz
    `checkFilterName()`_
    `filterProfile()`_
    `isUnit()`_

    .. _`checkFilterName()` : api.html#splat.photometry.checkFilterName
    .. _`filterProfile()` : api.html#splat.photometry.filterProfile
    .. _`isUnit()` : api.html#splat.utilities.isUnit
    .. _`maskFlux()` : api.html#splat.core.maskFlux

    :Example:
    '''
# keyword parameters
    evals = ['vega','ab','photons','energy']
    for x in evals:
        if kwargs.get(x,False) == True: computevalue = x
    if computevalue not in evals:
        if verbose==True: ('WARNING: do not recongize output type {}; using vega'.format(computevalue))
        computevalue = 'vega'

# check that requested filter is in list
    if isinstance(custom,bool) and isinstance(notch,bool):
        f0 = checkFilterName(filt,verbose=True)
        if f0 == False: 
            if verbose==True: print('\nWarning: SPLAT does not have filter {}; print splat.FILTERS.keys() to get list of current filters '.format(filt))
            return numpy.nan, numpy.nan
        fwave,ftrans = filterProfile(f0,**kwargs)
# notch filter
    elif isinstance(custom,bool) and isinstance(notch,list):
        dn = (notch[-1]-notch[0])/1000
        fwave = numpy.arange(notch[0]-5.*dn,notch[-1]+5.*dn,dn)
        ftrans = numpy.zeros(len(fwave))
        ftrans[numpy.where(numpy.logical_and(fwave >= notch[0],fwave <= notch[-1]))] = 1.
# custom filter
    else:
        fwave,ftrans = custom[0],custom[1]
# units
    if isinstance(fwave,u.quantity.Quantity) == True: fwave = fwave.to(sp.wave.unit)
    else: fwave = (fwave*u.micron).to(sp.wave.unit)


# check that spectrum and filter cover the same wavelength ranges
    if numpy.nanmax(fwave.value) < numpy.nanmin(sp.wave.value) or numpy.nanmin(fwave.value) > numpy.nanmax(sp.wave.value):
        if verbose==True: print('\nWarning: no overlap between spectrum for {} and filter {}'.format(sp.name,filt))
        return numpy.nan, numpy.nan

    if numpy.nanmin(fwave.value) < numpy.nanmin(sp.wave.value) or numpy.nanmax(fwave.value) > numpy.nanmax(sp.wave.value):
        if verbose==True: print('\nWarning: spectrum for {} does not span full filter profile for {}'.format(sp.name,filt))

# interpolate spectrum onto filter wavelength function
    wgood = numpy.where(numpy.logical_and(~numpy.isnan(sp.noise.value),~numpy.isnan(sp.flux.value)))
    if len(sp.wave[wgood]) > 0:
        d = interp1d(sp.wave[wgood].value,sp.flux[wgood].value,bounds_error=False,fill_value=0.)
        n = interp1d(sp.wave[wgood].value,sp.noise[wgood].value,bounds_error=False,fill_value=0)
# catch for models
    else:
        if verbose==True: print('\nWarning: data values in range of filter {} have no uncertainties'.format(filt))
        d = interp1d(sp.wave.value,sp.flux.value,bounds_error=False,fill_value=0.)
        n = interp1d(sp.wave.value,sp.flux.value*1.e-9,bounds_error=False,fill_value=0.)

    result = []
    if computevalue== 'vega':
# make sure vega file is there
        if os.path.exists(vegafile)==False:
            tmp = os.path.exists(filterfolder,vegafile)
            if os.path.exists(tmp)==False:
                raise ValueError('Cannot find vega file {} locally or in {}'.format(vegafile,FILTER_FOLDER))
            vegafile=tmp

# Read in Vega spectrum
        vwave,vflux = numpy.genfromtxt(os.path.normpath(vegafile), comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
        vwave = vwave[~numpy.isnan(vflux)]*u.micron
        vwave.to(sp.wave.unit)
        vflux = vflux[~numpy.isnan(vflux)]*(u.erg/(u.cm**2 * u.s * u.micron))
        vflux.to(sp.flux.unit,equivalencies=u.spectral_density(vwave))
# interpolate Vega onto filter wavelength function
        v = interp1d(vwave.value,vflux.value,bounds_error=False,fill_value=0.)
        if rsr == True:
            val = -2.5*numpy.log10(trapz(ftrans*fwave.value*d(fwave.value),fwave.value)/trapz(ftrans*fwave.value*v(fwave.value),fwave.value))
        else:
            val = -2.5*numpy.log10(trapz(ftrans*d(fwave.value),fwave.value)/trapz(ftrans*v(fwave.value),fwave.value))
        for i in numpy.arange(nsamples):
#            result.append(-2.5*numpy.log10(trapz(ftrans*numpy.random.normal(d(fwave),n(fwave))*sp.flux_unit,fwave)/trapz(ftrans*v(fwave)*sp.flux_unit,fwave)))
            if rsr == True:
                result.append(-2.5*numpy.log10(trapz(ftrans*fwave.value*numpy.random.normal(d(fwave.value),n(fwave.value)),fwave.value)/trapz(ftrans*fwave.value*v(fwave.value),fwave.value)))
            else:
                result.append(-2.5*numpy.log10(trapz(ftrans*numpy.random.normal(d(fwave.value),n(fwave.value)),fwave.value)/trapz(ftrans*v(fwave.value),fwave.value)))
        outunit = 1.

    elif computevalue== 'ab':
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
            a = trapz(ftrans*numpy.random.normal(d(filtnu.value),n(filtnu.value))/filtnu.value,filtnu.value)
            result.append(-2.5*numpy.log10(a/b))
        outunit = 1.

    elif computevalue== 'energy':
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
                result.append((trapz(ftrans*fwave.value*numpy.random.normal(d(fwave.value),n(fwave.value)),fwave.value)*sp.wave.unit*sp.flux.unit).to(outunit).value)
            else:
                result.append((trapz(ftrans*numpy.random.normal(d(fwave.value),n(fwave.value)),fwave.value)*sp.wave.unit*sp.flux.unit).to(outunit).value)

    elif computevalue== 'photons':
        outunit = 1./u.s/u.cm**2
        nu = sp.wave.to('Hz',equivalencies=u.spectral())
        fnu = sp.flux.to('Jy',equivalencies=u.spectral_density(sp.wave))
        noisenu = sp.noise.to('Jy',equivalencies=u.spectral_density(sp.wave))
        filtnu = fwave.to('Hz',equivalencies=u.spectral())
        d = interp1d(nu.value,fnu.value,bounds_error=False,fill_value=0.)
        n = interp1d(nu.value,noisenu.value,bounds_error=False,fill_value=0.)
#        convert = 1./((const.h.to('erg s')*numpy.nanmedian(nu.value)*u.Hz).to(u.erg))
        val = -1.*(trapz(ftrans*(d(filtnu.value)/(const.h.to('erg s').value*filtnu.value)),filtnu.value)*(u.Jy*u.Hz/u.erg)).to(outunit).value
        for i in numpy.arange(nsamples):
            result.append((trapz(ftrans*(numpy.random.normal(d(filtnu.value),n(filtnu.value))/(const.h.to('erg s').value*filtnu.value)),filtnu.value)*(u.Jy*u.Hz/u.erg)).to(outunit).value)
#            result.append((trapz(ftrans*fwave.value*convert.value*numpy.random.normal(d(fwave.value),n(fwave.value)),fwave.value)*sp.wave.unit*sp.flux.unit*convert.unit).to(outunit).value)
    else:
        raise NameError('\nfilterMag not given a correct physical quantity (vega, ab, energy, photons) to compute photometry\n\n')

#    val = numpy.nanmean(result)*outunit
    err = numpy.nanstd(result)
    if len(sp.wave[wgood]) == 0: err = 0.
    return val*outunit,err*outunit


def vegaToAB(filt,vegafile=VEGAFILE,filterfolder=FILTER_FOLDER,custom=False,notch=False,rsr=False,**kwargs):

# check that requested filter is in list
    if isinstance(custom,bool) and isinstance(notch,bool):
        f0 = checkFilterName(filt,verbose=True)
        if f0 == False: 
            return numpy.nan, numpy.nan
        filt = f0
        rsr = FILTERS[filt]['rsr']

# make sure vega file is there
    if os.path.exists(vegafile)==False:
        tmp = os.path.exists(filterfolder,vegafile)
        if os.path.exists(tmp)==False:
            raise ValueError('Cannot find vega file {} locally or in {}'.format(vegafile,FILTER_FOLDER))
        vegafile=tmp

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

# Read in Vega spectrum
    vwave,vflux = numpy.genfromtxt(os.path.normpath(vegafile), comments='#', unpack=True, \
        missing_values = ('NaN','nan'), filling_values = (numpy.nan))
    vwave = vwave[~numpy.isnan(vflux)]*u.micron
    vflux = vflux[~numpy.isnan(vflux)]*(u.erg/(u.cm**2 * u.s * u.micron))

# trim spectrum
    vflux = vflux[vwave>=numpy.nanmin(fwave)]
    vwave = vwave[vwave>=numpy.nanmin(fwave)]
    vflux = vflux[vwave<=numpy.nanmax(fwave)]
    vwave = vwave[vwave<=numpy.nanmax(fwave)]

# convert to fnu
    nu = vwave.to('Hz',equivalencies=u.spectral())
    fnu = vflux.to('Jy',equivalencies=u.spectral_density(vwave))
    filtnu = fwave.to('Hz',equivalencies=u.spectral())
    fconst = 3631*u.jansky
    d = interp1d(nu.value,fnu.value,bounds_error=False,fill_value=0.)
    b = trapz((ftrans/filtnu.value)*fconst.value,filtnu.value)
    return -2.5*numpy.log10(trapz(ftrans*d(filtnu.value)/filtnu.value,filtnu.value)/b)



def filterInfo(fname='',verbose=True,**kwargs):
    '''
    :Purpose: Prints out the current list of filters in the SPLAT reference library.
    '''

    if fname=='': fname = list(FILTERS.keys())
    if isinstance(fname,list) == False: fname = [fname]

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
            output[f]['vega_ab'] = vegaToAB(f)
            if verbose ==True: 
                print('\n{}: {}'.format(f.replace('_',' '),output[f]['description']))
                print('\tWavelength range = {:.3f} to {:.3f}\n'.format(output[f]['lambda_min'],output[f]['lambda_max']))             
                print('\tZeropoint = {:.2f} Jy'.format(output[f]['zeropoint']))
                print('\tCentral wavelength: = {:.3f}'.format(output[f]['lambda_central']))
                print('\tMean wavelength: = {:.3f}'.format(output[f]['lambda_mean']))
                print('\tPivot point: = {:.3f}'.format(output[f]['lambda_pivot']))
                print('\tFWHM = {:.3f}'.format(output[f]['lambda_fwhm']))
                print('\tVega to AB = {:.3f}'.format(output[f]['vega_ab']))
        else:
        	if verbose ==True: print('\tFilter {} not in SPLAT filter list'.format(k))
    return output
    # kys = list(output.keys())
    # if len(kys) == 1: return output[kys[0]]
    # else: return output

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
    filterfolder = kwargs.get('filterfolder',FILTER_FOLDER)
    if not os.path.exists(filterfolder):
        raise ValueError('Could not find filter folder {}'.format(filterfolder))
    vegafile = kwargs.get('vegafile',VEGAFILE)
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
    # make sure vega file is there
            if os.path.exists(vegafile)==False:
                tmp = os.path.exists(filterfolder,vegafile)
                if os.path.exists(tmp)==False:
                    raise ValueError('Cannot find vega file {} locally or in {}'.format(vegafile,FILTER_FOLDER))
                vegafile=tmp

    # Read in Vega spectrum
            vwave,vflux = numpy.genfromtxt(os.path.normpath(vegafile), comments='#', unpack=True, \
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




def visualizeFilter(fprof,spectrum=None,file='',color='k',linestyle='-',fill=False,legend=False,nsample=100,verbose=True,**kwargs):
    '''
    Purpose
    -------

    Plots a filter or set of filter transmission curves, and optionally a comparison spectrum.

    Parameters
    ----------

    fprof : str or list/array
        Either of string or list of strings of predefined filter names listed in splat.FILTERS.keys()
        or a list of numerical arrays: lists of pairs of numbers => notch filters; list of pairs of number arrays => wave,transmission profiles
    
    spectrum : Spectrum class object, default = None
        Assign to a spectrum class object you want to overplot; spectrum should be normalized in advance

    file : str, default = ''
        Full filepath for plot output; if empty string, no plot created

    color : str or list, default = 'k''
        single or list of color labels

    linestyle : str or list, default = '-''
        single or list of linestyles

    fill : bool, default = False
        Set to True to fill in the filter profile region

    legend : bool, default = False
        Set to True to include legend

    nsample : int, default = 100
        Number of samples for a notch filter

    verbose : bool, default = False
        Set to True to receive verbose feedback

    additional plotting options contained in matplotlib.pyplot can be passed as well, 
    including figsize, xlim, ylim, xlabel, ylabel, wave_unit, normalize, spectrum_scale, spectrum_color, spectrum_alpha

    Outputs
    -------

    fig : matplotlib figure object
        figure object for plot

    Example
    -------
    
    >>> import splat
    >>> import splat.photometry as spphot
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
    >>> sp.normalize()
    >>> spphot.visualizeFilter(['MKO Y','MKO J','MKO H','MKO K'],comparison=sp,color=['c','b','g','m'],fill=True) 
        [plot of filter tranmission curves and spectrum object] 

    Dependencies
    ------------
    copy
    matplotlib.pyplot
    numpy
    `checkFilterName()`_
    `isUnit()`_
    `filterProfile()`_

    .. _`checkFilterName()` : api.html#splat.photometry.checkFilterName
    .. _`filterProfile()` : api.html#splat.photometry.filterProfile
    .. _`isUnit()` : api.html#splat.utilities.isUnit
    '''
    filt = copy.deepcopy(fprof)
    wave_unit = kwargs.get('wave_unit',DEFAULT_WAVE_UNIT)

# filter names  
    if isinstance(filt,str): filt = [filt]
    if isinstance(filt,list):

# prep parameters
        fig = plt.figure(figsize=kwargs.get('figsize',[5,4]))
        if isinstance(color,str): color = [color]
        while len(color) < len(filt): color.append(color[-1])
        if isinstance(linestyle,str): linestyle = [linestyle]
        while len(linestyle) < len(filt): linestyle.append(linestyle[-1])

# filter names  
        if isinstance(filt[0],str):
            for i,f in enumerate(filt):
                fc = checkFilterName(f)
                if fc==False:
                    if verbose==True: print('Filter {} not included in SPLAT; skipping'.format(f))
                else:
                    fwave,ftrans = filterProfile(f,**kwargs)
                    if isUnit(fwave): fwave = fwave.to(wave_unit).value
    #                else: fwave = fwave*wave_unit
                    if isUnit(ftrans): ftrans = ftrans.value
                    if kwargs.get('normalize',False): ftrans = ftrans/numpy.nanmax(ftrans)
                    plt.plot(fwave,ftrans,label=f.replace('_',' '),color=color[i],ls=linestyle[i])
                    # legend.append(FILTERS[f]['description'])
                    # if FILTERS[f]['rsr'] == True: yl = kwargs.get('ylabel','Transmission Curve')
                    if fill==True: plt.fill_between(fwave,ftrans,numpy.zeros(len(ftrans)),color=color[i],alpha=0.1)

# number array
        if isinstance(filt[0],int) or isinstance(filt[0],float): filt = [filt]
        if isinstance(filt[0],list):
            for i,f in enumerate(filt):
# notch filters
                if isinstance(f[0],int) or isinstance(f[0],float):
                    fwave,ftrans = numpy.linspace(f[0],f[1],nsample)*wave_unit,numpy.ones(nsample)
                else:
                    fwave,ftrans = f[0],f[1]
                    if isUnit(fwave): fwave = fwave.to(wave_unit).value
    #                else: fwave = fwave*wave_unit
                    if isUnit(ftrans): ftrans = ftrans.value
                    if kwargs.get('normalize',False): ftrans = ftrans/numpy.nanmax(ftrans)
                plt.plot(fwave,ftrans,label='Filter {}'.format(i+1),color=colors[i],ls=linestyle[i])
                if fill==True: plt.fill_between(fwave,ftrans,numpy.zeros(len(ftrans)),color=color[i],alpha=0.1)
    else:
        raise ValueError('Could not parse input {}'.format(filt))

# add a comparison spectrum
    #spectrum = kwargs.get('comparison',spectrum)
    if isinstance(spectrum,splat.Spectrum) == True:
        # print(xra)
#        spectrum.normalize(xra)
        spectrum.scale(numpy.nanmax(ftrans)*kwargs.get('spectrum_scale',0.8))
        plt.plot(spectrum.wave,spectrum.flux,color=kwargs.get('spectrum_color','k'),alpha=kwargs.get('spectrumalpha',1.0),label=spectrum.name)
        plt.plot(spectrum.wave,spectrum.noise,color='k',ls='--')

# finish up
    if len(kwargs.get('xlim',[]))>0: plt.xlim(kwargs.get('xlim'))
    if len(kwargs.get('ylim',[]))>0: plt.ylim(kwargs.get('ylim'))
    plt.xlabel(kwargs.get('xlabel','Wavelength ({})'.format(wave_unit)),fontsize=14)
    plt.ylabel(kwargs.get('ylabel','Transmission'),fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    if legend==True: plt.legend(fontsize=14)
    plt.tight_layout()

# save if desired
    for x in ['file','filename','output']: file = kwargs.get(x,file)
    if file != '': plt.savefig(file)
    plt.show()
    
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

def SEDVisualize(verbose=True):
    '''
    Visualizes magnitudes on SED scale (flux = lam x F_lam), with option of also comparing to SED spectrum
    '''
    pass



#####################################################
###############    MAGNITUDE CLASS    ###############
##############    UNDER DEVELOPMENT    ##############
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


