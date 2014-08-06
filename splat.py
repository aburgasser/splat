# WORKING COPY OF SPLAT CODE LIBRARY
# based on routines developed by:
#    Christian Aganze
#    Daniella Bardalez Gagliuffi
#    Adam Burgasser
#    Caleb Choban
#    Ivanna Escala
#    Aishwarya Iyer
#    Yuhui Jin
#    Mike Lopez
#    Alex Mendez
#    Johnny Parra
#    Julian Pilate-Hutcherson
#    Maitrayee Sahi
#    Melisa Tallis

#
# CURRENT STATUS (5/6/2014)
# can now load up spectra and models from online sources
# source spectra can be selected by name, designation, young/subdwarf/red/blue/binary/spbin
# returned as array of spectra
#

# imports
import astropy
import base64
import copy
import os
import matplotlib.pyplot as plt
import numpy
import re
import scipy
import string
import sys
import urllib2
import warnings
from scipy import stats, signal
from scipy.integrate import trapz        # for numerical integration
from scipy.interpolate import interp1d, griddata
from astropy.io import ascii, fits            # for reading in spreadsheet
from astropy.table import Table, join            # for reading in table files
from astropy.coordinates import SkyCoord      # coordinate conversion
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units

# suppress warnings - probably not an entirely safe approach!
numpy.seterr(all='ignore')
warnings.simplefilter("ignore")
#from splat._version import __version__

############ CONSTANTS - THESE SHOULD STAY FIXED ############
SPLAT_URL = 'http://pono.ucsd.edu/~adam/splat/'
HOME_DIR = os.environ['HOME']
access_file = '.splat_access'
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
parameter_names = ['teff','logg','z','fsed','kzz']
spex_pixel_scale = 0.15            # spatial scale in arcseconds per pixel
spex_wave_range = [0.65,2.45]    # default wavelength range
max_snr = 1000.0                # maximum S/N ratio permitted
#############################################################

        
# helper functions from Alex
def lazyprop(fn):
     attr_name = '_lazy_' + fn.__name__
     @property
     def _lazyprop(self):
          if not hasattr(self, attr_name):
                setattr(self, attr_name, fn(self))
          return getattr(self, attr_name)
     return _lazyprop

def Show(fn):
     def _show(self, *args, **kwargs):
          noplot = kwargs.pop('noplot', False)
          quiet = kwargs.pop('quiet', False)
          tmp = fn(self, *args, **kwargs)
#         if not quiet:
#                self.info()
#          if not noplot:
#                self.plot(**kwargs)
          return tmp
     return _show

def Copy(fn):
     def _copy(self, *args, **kwargs):
          out = copy.copy(self)
          return fn(out, *args, **kwargs)
     return _copy


# define the Spectrum class which contains the relevant information
class Spectrum(object):
    @Show
    def __init__(self, filename, **kwargs):
        '''Load the file'''
        self.model = False
# various options for setting the filename
        self.filename = filename
        if kwargs.get('filename','') != '':
            self.filename = kwargs.get('filename','')
        if kwargs.get('file','') != '':
            self.filename = kwargs.get('file','')
        kwargs['filename'] = self.filename
        self.wlabel = kwargs.get('wlabel','Wavelength')
        self.wunit = kwargs.get('wunit',u.micron)
        self.flabel = kwargs.get('flabel','F_lambda')
        self.fscale = kwargs.get('fscale','')
        self.funit = kwargs.get('funit',u.erg/(u.cm**2 * u.s * u.micron))
        self.resolution = kwargs.get('resolution',150)    # default placeholder
        self.slitpixelwidth = kwargs.get('slitwidth',3.33)        # default placeholder
        self.slitwidth = self.slitpixelwidth*spex_pixel_scale
        self.simplefilename = os.path.basename(self.filename)
# wave and flux given
        if len(kwargs.get('wave','')) > 0 and len(kwargs.get('flux','')) > 0:
            self.wave = kwargs['wave']
            self.flux = kwargs['flux']
            if len(kwargs.get('noise','')) > 0:
                self.noise = kwargs['noise']
            else:
                self.noise = numpy.array([numpy.nan for i in self.wave])
        else:
# filename given
            try:
                self.wave, self.flux, self.noise = readSpectrum(**kwargs)
            except:
                raise NameError('\nCould not load up spectral file {:s}'.format(kwargs.get('filename','')))
        self.wave = numpy.array(self.wave)
        self.flux = numpy.array(self.flux)
        self.noise = numpy.array(self.noise)
# enforce positivity
        if (numpy.nanmin(self.flux) < 0):
            self.flux[numpy.where(self.flux < 0)] = 0.
# check on noise being too low
        if (numpy.nanmax(self.flux/self.noise) > max_snr):
            self.noise[numpy.where(self.flux/self.noise > max_snr)]=numpy.median(self.noise)

# some conversions
        self.nu = const.c.to('micron/s').value/self.wave
#        w = numpy.where(numpy.isnan(self.flux))
#        self.flux[w] = 0.
        self.flux[numpy.isnan(self.flux)] = 0.
# calculate variance
        self.variance = numpy.array([n**2 for n in self.noise])
        self.dof = numpy.round(len(self.wave)/self.slitpixelwidth)

# preserve original values
        self.wave_original = copy.deepcopy(self.wave)
        self.flux_original = copy.deepcopy(self.flux)
        self.noise_original = copy.deepcopy(self.noise)
        self.variance_original = copy.deepcopy(self.variance)
        self.resolution = copy.deepcopy(self.resolution)
        self.slitpixelwidth = copy.deepcopy(self.slitpixelwidth)
        
# information on source spectrum
        if not (kwargs.get('model',False)):
            x,y = filenameToNameDate(self.filename)
            self.name = kwargs.get('name',x)
            self.date = kwargs.get('date',y)
            try:
                self.caldate = dateToCaldate(self.date)
            except:
                self.caldate = ''
        else:
# information on model
            self.model = True
            self.teff = kwargs.get('teff',numpy.nan)
            self.logg = kwargs.get('logg',numpy.nan)
            self.z = kwargs.get('z',numpy.nan)
            self.cloud = kwargs.get('cloud',numpy.nan)
            self.modelset = kwargs.get('set','')
            self.name = self.modelset+' Teff='+str(self.teff)+' logg='+str(self.logg)
            self.fscale = 'Surface'
        self.history = ['Loaded']

                
    def __repr__(self):
        '''A simple representation of an object is to just give it a name'''
        return 'Spectra Object for {}'.format(self.name)

    def __add__(self,other):
        '''Adding two spectra '''
        sp = copy.deepcopy(self)
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.add(self.flux,f(self.wave))
        sp.variance = sp.variance+n(self.wave)
        sp.noise = numpy.array([v**0.5 for v in sp.variance])
        sp.flux_original=sp.flux
        sp.noise_original=sp.noise
        sp.variance_original=sp.variance
        return sp

    def __sub__(self,other):
        '''Subtracting two spectra '''
        sp = copy.deepcopy(self)
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.subtract(self.flux,f(self.wave))
        sp.variance = sp.variance+n(self.wave)
        sp.noise = numpy.array([v**0.5 for v in sp.variance])
        sp.flux_original=sp.flux
        sp.noise_original=sp.noise
        sp.variance_original=sp.variance
        return sp

    def __mul__(self,other):
        '''Multiplying two spectra'''
        sp = copy.deepcopy(self)
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.multiply(self.flux,f(self.wave))
        sp.variance = numpy.multiply(numpy.power(sp.flux,2),(\
            numpy.divide(self.variance,numpy.power(sp.flux,2))+\
            numpy.divide(n,numpy.power(f(self.wave),2))))
        sp.noise = numpy.array([v**0.5 for v in sp.variance])
        sp.flux_original=sp.flux
        sp.noise_original=sp.noise
        sp.variance_original=sp.variance
        return sp

    def __div__(self,other):
        '''Dividing two spectra'''
        sp = copy.deepcopy(self)
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.divide(self.flux,f(self.wave))
        sp.variance = numpy.multiply(numpy.power(sp.flux,2),(\
            numpy.divide(self.variance,numpy.power(sp.flux,2))+\
            numpy.divide(n,numpy.power(f(self.wave),2))))
        sp.noise = numpy.array([v**0.5 for v in sp.variance])
        sp.flux_original=sp.flux
        sp.noise_original=sp.noise
        sp.variance_original=sp.variance
        return sp

    def copy(self):
          '''Make a copy of the current spectrum'''
          other = copy.deepcopy(self)
          return other

    def info(self):
          '''Report some information about this spectrum'''
          if (self.model):
              print '''{0} model with Teff = {1} and log g = {2}'''.format(self.modelset, self.teff, self.logg)
          else:
              print '''Spectrum of {0} taken on {1}'''.format(self.name, self.date)
          return

    def flamToFnu(self):
         '''Convert flux density from F_lam to F_nu, the later in Jy'''
         self.funit = u.Jy
         self.flabel = 'F_nu'
         self.flux.to(self.funit,equivalencies=u.spectral_density(self.wave))
         self.noise.to(self.funit,equivalencies=u.spectral_density(self.wave))
         return

    def fluxCalibrate(self,filter,mag,**kwargs):
        '''Calibrate spectrum to input magnitude'''
        absolute = kwargs.get('absolute',False)
        apparent = kwargs.get('apparent',False)
        self.normalize()
        apmag = filterMag(self,filter,**kwargs)
        if (~numpy.isnan(apmag)):
            self.scale(10.**(0.4*(apmag-mag)))
            if (absolute):
                self.fscale = 'Absolute'
            if (apparent):
                self.fscale = 'Apparent'
        return

    def fluxMax(self):
        return numpy.nanmax(self.flux[numpy.where(\
            numpy.logical_and(self.wave > 0.8,self.wave < 2.3))])

    def fnuToFlam(self):
         '''Convert flux density from F_nu to F_lam, the later in erg/s/cm2/Hz'''
         self.funit = u.erg/(u.cm**2 * u.s * u.micron)
         self.flabel = 'F_lam'
         self.flux.to(self.funit,equivalencies=u.spectral_density(self.wave))
         self.noise.to(self.funit,equivalencies=u.spectral_density(self.wave))
         self.variance = [n**2 for n in self.noise]
         return

    def normalize(self):
        '''Normalize spectrum'''
        self.scale(1./self.fluxMax())
        self.fscale = 'Normalized'
        return
    def reset(self):
        '''Reset to original spectrum'''
        self.wave = copy.deepcopy(self.wave_original)
        self.flux = copy.deepcopy(self.flux_original)
        self.noise = copy.deepcopy(self.noise_original)
        self.variance = copy.deepcopy(self.variance_original)
        self.resolution = copy.deepcopy(self.resolution_original)
        self.slitpixelwidth = copy.deepcopy(self.slitpixelwidth_original)
        self.slitwidth = self.slitpixelwidth*spex_pixel_scale
        self.fscale = ''
        return

    def scale(self,factor):
         '''Scale spectrum and noise by a constant factor'''
         self.flux = self.flux*factor
         self.noise = numpy.array([n*factor for n in self.noise])
         self.variance = numpy.array([n**2 for n in self.noise])
         self.fscale = 'Scaled'
         return
        
    def smooth(self,**kwargs):
        '''Smooth spectrum to a constant slit width (smooth by pixels)'''
        method = kwargs.get('method','hanning')
        kwargs['method'] = method
        swargs = copy.deepcopy(kwargs)
        if (kwargs.get('slitPixelWidth','') != ''):
            del swargs['slitPixelWidth']
            self.smoothToSlitPixelWidth(kwargs['slitPixelWidth'],**swargs)
        elif (kwargs.get('resolution','') != ''):
            del swargs['resolution']
            self.smoothToResolution(kwargs['resolution'],**swargs)
        elif (kwargs.get('slitWidth','') != ''):
            del swargs['slitWidth']
            self.smoothToSlitWidth(kwargs['slitWidth'],**swargs)
        return

    def smoothToResolution(self,resolution,**kwargs):
        '''Smooth spectrum to a constant resolution'''
        overscale = kwargs.get('overscale',10.)
        method = kwargs.get('method','hanning')
        kwargs['method'] = method

# do nothing if requested resolution is higher than current resolution
        if (resolution < self.resolution):
# sample onto a constant resolution grid at 5x current resolution
             r = resolution*overscale
             waveRng = self.waveRange()
             npix = numpy.floor(numpy.log(waveRng[1]/waveRng[0])/numpy.log(1.+1./r))
             wave_sample = [waveRng[0]*(1.+1./r)**i for i in numpy.arange(npix)]
             f = interp1d(self.wave,self.flux,bounds_error=False,fill_value=0.)
             v = interp1d(self.wave,self.variance,bounds_error=False,fill_value=numpy.nan)
             flx_sample = f(wave_sample)
             var_sample = v(wave_sample)
# now convolve a function to smooth resampled spectrum
             window = signal.get_window(method,numpy.round(overscale))
             neff = numpy.sum(window)/numpy.nanmax(window)        # effective number of pixels
             flx_smooth = signal.convolve(flx_sample, window/numpy.sum(window), mode='same')
             var_smooth = signal.convolve(var_sample, window/numpy.sum(window), mode='same')/neff
# resample back to original wavelength grid
             f = interp1d(wave_sample,flx_smooth,bounds_error=False,fill_value=0.)
             v = interp1d(wave_sample,var_smooth,bounds_error=False,fill_value=0.)
             self.flux = f(self.wave)
             self.variance = v(self.wave)
             self.noise = numpy.array([n**0.5 for n in self.variance])
             self.slitpixelwidth = self.slitpixelwidth*self.resolution/resolution
             self.resolution = resolution
             self.slitwidth = self.slitpixelwidth*spex_pixel_scale
             self.history = ['Smoothed to constant resolution {}'.format(self.resolution)]        
        return

    def smoothToSlitPixelWidth(self,width,**kwargs):
        '''Smooth spectrum to a constant slit width (smooth by pixels)'''
        method = kwargs.get('method','hanning')
        kwargs['method'] = method
# do nothing if requested resolution is higher than current resolution
        if (width > self.slitpixelwidth):
# convolve a function to smooth spectrum
            window = signal.get_window(method,numpy.round(width))
            neff = numpy.sum(window)/numpy.nanmax(window)        # effective number of pixels
            self.flux = signal.convolve(self.flux, window/numpy.sum(window), mode='same')
            self.variance = signal.convolve(self.variance, window/numpy.sum(window), mode='same')/neff
            self.noise = numpy.array([n**0.5 for n in self.variance])
            self.resolution = self.resolution*self.slitpixelwidth/width
            self.slitpixelwidth = width
            self.slitwidth = self.slitpixelwidth*spex_pixel_scale
            self.history = ['Smoothed to slit width of {}'.format(self.slitwidth)]        
        return

    def smoothToSlitWidth(self,width,**kwargs):
        method = kwargs.get('method','hanning')
        kwargs['method'] = method
        '''Smooth spectrum to a constant slit width (smooth by pixels)'''
        pwidth = width/spex_pixel_scale
        self.smoothToSlitPixelWidth(pwidth,**kwargs)
        return

    def snr(self):
         '''Compute a representative S/N value'''
         pass
         return

    def surface(self,radius):
         '''Convert to surface fluxes given a radius, assuming at absolute fluxes'''
         pass
         return

    def waveRange(self):
        ii = numpy.where(self.flux > 0)
        return [numpy.nanmin(self.wave[ii]), numpy.nanmax(self.wave[ii])]
     
        
                             

# FUNCTIONS FOR SPLAT
def caldateToDate(d):
    '''Convert from numeric date to calendar date'''
    return d[:4]+str((months.index(d[5:8])+1)/100.)[2:4]+d[-2:]


def checkFile(filename,**kwargs):
    '''Check if a particular file is present in the online database'''
    url = kwargs.get('url',SPLAT_URL)+'/Spectra/'
    flag = checkOnline()
    if (flag):
        try:
            open(os.path.basename(filename), 'wb').write(urllib2.urlopen(url+filename).read())
        except urllib2.URLError:
            flag = False
    return flag


def checkAccess():
    try:
        bcode = urllib2.urlopen(SPLAT_URL+access_file).read()
        lcode = base64.b64encode(open(HOME_DIR+'/'+access_file,'r').read())
        if (lcode == bcode):
            return True
        else:
            return False
    except:
        return False
        

def checkOnline():
    '''Check if you are online'''
    try:
        urllib2.urlopen(SPLAT_URL)
        return True
    except urllib2.URLError:
        return False



def classifyByIndex(sp, *args, **kwargs):
    '''Classify a spectrum based on its spectral indices'''
    
    str_flag = kwargs.get('string', True)
    rnd_flag = kwargs.get('round', False)
    rem_flag = kwargs.get('remeasure', True)
    nsamples = kwargs.get('nsamples', 100)
    nloop = kwargs.get('nloop', 5)
    set = kwargs.get('set','burgasser')
    allowed_sets = ['burgasser','reid','testi','allers']
    if (set.lower() not in allowed_sets):
        print '\nWarning: index classification method {} not present; returning nan\n\n'.format(set)
        return numpy.nan, numpy.nan

# measure indices if necessary
    if (len(args) != 0):
        indices = args[0]

# Burgasser (2007, ApJ, 659, 655) calibration
    if (set.lower() == 'burgasser'):
        if (rem_flag or len(args) == 0):
            indices = measureIndexSet(sp, **kwargs)
        sptoffset = 20.
        sptfact = 1.
        coeffs = { \
            'H2O-J': {'fitunc': 0.8, 'range': [20,39], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [1.038e2, -2.156e2,  1.312e2, -3.919e1, 1.949e1]}, \
            'H2O-H': {'fitunc': 1.0, 'range': [20,39], 'spt': 0., 'sptunc': 99., 'mask': 1.,  \
            'coeff': [9.087e-1, -3.221e1, 2.527e1, -1.978e1, 2.098e1]}, \
            'CH4-J': {'fitunc': 0.7, 'range': [30,39], 'spt': 0., 'sptunc': 99., 'mask': 1.,  \
            'coeff': [1.491e2, -3.381e2, 2.424e2, -8.450e1, 2.708e1]}, \
            'CH4-H': {'fitunc': 0.3, 'range': [31,39], 'spt': 0., 'sptunc': 99., 'mask': 1.,  \
            'coeff': [2.084e1, -5.068e1, 4.361e1, -2.291e1, 2.013e1]}, \
            'CH4-K': {'fitunc': 1.1, 'range': [20,37], 'spt': 0., 'sptunc': 99., 'mask': 1.,  \
            'coeff': [-1.259e1, -4.734e0, 2.534e1, -2.246e1, 1.885e1]}}

# Reid et al. (2001, AJ, 121, 1710)
    elif (set.lower() == 'reid'):
        if (rem_flag or len(args) == 0):
            indices = measureIndexSet(sp, **kwargs)
        sptoffset = 20.
        sptfact = 1.
        coeffs = { \
            'H2O-A': {'fitunc': 1.18, 'range': [18,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [-32.1, 23.4]}, \
            'H2O-B': {'fitunc': 1.02, 'range': [18,28], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [-24.9, 20.7]}}

# Testi et al. (2001, ApJ, 522, L147)
    elif (set.lower() == 'testi'):
        if (rem_flag or len(args) == 0):
            indices = measureIndexSet(sp, **kwargs)
        sptoffset = 20.
        sptfact = 10.
        coeffs = { \
            'sHJ': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [-1.87, 1.67]}, \
            'sKJ': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [-1.20, 2.01]}, \
            'sH2O_J': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [1.54, 0.98]}, \
            'sH2O_H1': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [1.27, 0.76]}, \
            'sH2O_H2': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [2.11, 0.29]}, \
            'sH2O_K': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [2.36, 0.60]}}

# Allers et al. (2013, ApJ, 657, 511)
    elif (set.lower() == 'allers'):
        if (rem_flag or len(args) == 0):
            kwargs['set'] = 'mclean'
            i1 = measureIndexSet(sp, **kwargs)
            kwargs['set'] = 'slesnick'
            i2 = measureIndexSet(sp, **kwargs)
            kwargs['set'] = 'allers'
            i3 = measureIndexSet(sp, **kwargs)
            indices = dict(i1.items() + i2.items() + i3.items())
        sptoffset = 10.
        sptfact = 1.
        coeffs = { \
            'H2O': {'fitunc': 0.390, 'range': [15,25], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [24.0476, -104.424, 169.388,-83.5437]}, \
            'H2O-1': {'fitunc': 1.097, 'range': [14,25], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [28.5982, -80.7404, 39.3513, 12.1927]}, \
            'H2OD': {'fitunc': 0.757, 'range': [20,28], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [-97.230, 229.884, -202.245, 79.4477]}, \
            'H2O-2': {'fitunc': 0.501, 'range': [14,22], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [37.5013, -97.8144, 55.4580, 10.8822]}}

    else:
        sys.stderr.write('\nWarning: '+set.lower()+' SpT-index relation not in measureSpT code\n\n')
        return numpy.nan, numpy.nan


    for index in coeffs.keys():
        vals = numpy.polyval(coeffs[index]['coeff'],numpy.random.normal(indices[index][0],indices[index][1],nsamples))*sptfact
        coeffs[index]['spt'] = numpy.nanmean(vals)+sptoffset
        coeffs[index]['sptunc'] = (numpy.nanstd(vals)**2+coeffs[index]['fitunc']**2)**0.5
#        print index, indices[index], numpy.nanmean(vals), numpy.nanstd(vals), coeffs[index]['spt']
    
#    print indices[index][0], numpy.polyval(coeffs[index]['coeff'],indices[index][0]), coeffs[index]
#    mask = numpy.ones(len(coeffs.keys()))
#    result = numpy.zeros(2)
    for i in numpy.arange(nloop):
        wts = [coeffs[index]['mask']/coeffs[index]['sptunc']**2 for index in coeffs.keys()]
        if (numpy.nansum(wts) == 0.):
            sys.stderr.write('\nIndices do not fit within allowed ranges\n\n')
            return numpy.nan, numpy.nan            
        vals = [coeffs[index]['mask']*coeffs[index]['spt']/coeffs[index]['sptunc']**2 \
            for index in coeffs.keys()]
        sptn = numpy.nansum(vals)/numpy.nansum(wts)
        sptn_e = 1./numpy.nansum(wts)**0.5
        for index in coeffs.keys():
            coeffs[index]['mask'] = numpy.where( \
                coeffs[index]['range'][0] <= sptn <= coeffs[index]['range'][1],1,0)

# round off to nearest 0.5 subtypes if desired
    if (rnd_flag):
        sptn = 0.5*numpy.around(sptn*2.)

# change to string if desired
    if (str_flag):
        spt = typeToNum(sptn,uncertainty=sptn_e)
    else:
        spt = sptn

    return spt, sptn_e


def classifyByStandard(sp, *args, **kwargs):
    '''Classify a spectrum by comparing to spectral standards'''

    str_flag = kwargs.get('string', True)
    method = kwargs.get('method','')
    best_flag = kwargs.get('best',False)
    plot_flag = kwargs.get('plot',False)
    unc_sys = 0.5

    stdsptnum = numpy.arange(30)+10.
#    stdsptstr = ['M'+str(n) for n in numpy.arange(10)]
#    stdsptstr.append(['L'+str(n) for n in numpy.arange(10)])
#    stdsptstr.append(['T'+str(n) for n in numpy.arange(10)])
    stdfiles = [ \
        'spex_prism_Gl270_091203.fits',\
        'spex_prism_Gl424_091229.fits',\
        'spex_prism_Gl91_081012.fits',\
        'spex_prism_Gl752A_070704.fits',\
        'spex_prism_Gl213_071110dl.fits',\
        'spex_prism_Gl51_070728.fits',\
        'spex_prism_LHS1375_081012.fits',\
        'spex_prism_vB8_070704.fits',\
        'spex_prism_vB10_070704.fits',\
        'spex_prism_LHS2924_070704.fits',\
        'spex_prism_0345+2540_kc.txt',\
        'spex_prism_2130-0845_080713.txt',\
        'spex_prism_kelu-1_060411.txt',\
        'spex_prism_1506+1321_060410.txt',\
        'spex_prism_2158-1550_060901.fits',\
        'spex_prism_sdss0835+19_chiu06.txt',\
        'spex_prism_1010-0406_kc.txt',\
        'spex_prism_0103+1935_060902.fits',\
        'spex_prism_1632+1904_030905.txt',\
        'spex_prism_denis0255-4700_040908.txt',\
        'spex_prism_1207+0244_061221dl.txt',\
        'spex_prism_0837-0000_061221.fits',\
        'spex_prism_sdss1254-0122_030522.txt',\
        'spex_prism_1209-1004_030523.fits',\
        'spex_prism_2254+3123_030918.txt',\
        'spex_prism_1503+2525_030522.txt',\
        'spex_prism_sdss1624+0029_040312.txt',\
        'spex_prism_0727+1710_040310.txt',\
        'spex_prism_0415-0935_030917.txt',\
        'spex_prism_0722-0540_110404.fits']

    if (method == 'kirkpatrick'):
        comprng = [0.9,1.4]         # as prescribed in Kirkpatrick et al. 2010, ApJS, 
    else:
        comprng = [0.7,2.45]        # by default, compare whole spectrum

# compute fitting statistics
    stat = []
    for file in stdfiles:
        spstd = loadSpectrum(file=file)
        chisq,scale = compareSpectra(sp,spstd,fit_ranges=[comprng],chisqr=True,novar2=True)
        stat.append(chisq)

# list of sorted standard files and spectral types
    sorted_stdfiles = [x for (y,x) in sorted(zip(stat,stdfiles))]
    sorted_stdsptnum = [x for (y,x) in sorted(zip(stat,stdsptnum))]

# select either best match or an ftest-weighted average
    if best_flag:
        sptn = sorted_stdsptnum[0]
        sptn_e = unc_sys
    else:
        mean,var = weightedMeanVar(stdsptnum,stat,method='ftest',dof=sp.dof)
        if (var**0.5 < 1.):
            sptn = numpy.round(mean*2)*0.5
        else:
            sptn = numpy.round(mean)
        sptn_e = (unc_sys**2+var)**0.5

# string or not?
    if (str_flag):
        spt = typeToNum(sptn,uncertainty=sptn_e)
    else:
        spt = sptn

# plot spectrum compared to best spectrum
    if (plot_flag):
        spc = sp.copy()     # copy to avoid change original value
        spc.normalize()
        spstd = loadSpectrum(file=sorted_stdfiles[0])
        chisq,scale = compareSpectra(spc,spstd,fit_ranges=[comprng],chisqr=True)
        spstd.scale(scale)
        plotSpectrum(spc,spstd,colors=['k','r'],\
            title=sp.name+' vs '+typeToNum(sorted_stdsptnum[0])+' Standard',**kwargs)

    return spt, sptn_e
    

def classifyByTemplate(sp, *args, **kwargs):
    '''Classify a spectrum by comparing to spectral templates'''
    return numpy.nan, numpy.nan
    
    
def compareSpectra(sp1, sp2, *args, **kwargs):
    '''Compare two spectra against each other'''
    weights = kwargs.get('weights',numpy.zeros(len(sp1.wave)))
    mask = kwargs.get('mask',numpy.zeros(len(sp1.wave)))    # mask = 1 -> ignore
    fit_ranges = kwargs.get('fit_ranges',[spex_wave_range])
    mask_ranges = kwargs.get('mask_ranges',[])
    mask_telluric = kwargs.get('mask_telluric',False)
    mask_standard = kwargs.get('mask_standard',False)
    var_flag = kwargs.get('novar2',False)
    chisqr = kwargs.get('chisqr',True)
    stddev = kwargs.get('stddev',False)
    absdev = kwargs.get('absdev',False)
    minreturn = 1.e-9

    if (absdev or stddev):
        chisqr = False
    
    if (mask_standard == True):
        mask_telluric == True
        
# create interpolation function for second spectrum
    f = interp1d(sp2.wave,sp2.flux,bounds_error=False,fill_value=0.)
    if var_flag:
        v = interp1d(sp2.wave,sp2.variance*numpy.nan,bounds_error=False,fill_value=numpy.nan)
    else:        
        v = interp1d(sp2.wave,sp2.variance,bounds_error=False,fill_value=numpy.nan)
    
# total variance - funny form to cover for nans
    vtot = numpy.nanmax([sp1.variance,sp1.variance+v(sp1.wave)],axis=0)
 #   vtot = sp1.variance
    
# Mask certain wavelengths
# telluric absorption
    if (mask_telluric):
        mask_ranges.append([0.,0.65])        # meant to clear out short wavelengths
        mask_ranges.append([1.35,1.42])
        mask_ranges.append([1.8,1.92])
        mask_ranges.append([2.45,99.])        # meant to clear out long wavelengths

    if (mask_standard):
        mask_ranges.append([0.,0.8])        # standard short cut
        mask_ranges.append([2.35,99.])        # standard long cut

    for ranges in mask_ranges:
        mask[numpy.where(((sp1.wave >= ranges[0]) & (sp1.wave <= ranges[1])))] = 1

# set the weights    
    for ranges in fit_ranges:
        weights[numpy.where(((sp1.wave >= ranges[0]) & (sp1.wave <= ranges[1])))] = 1
     
# mask flux < 0
    mask[numpy.where(numpy.logical_or(sp1.flux < 0,f(sp1.wave) < 0))] = 1

    weights = weights*(1.-mask)

    
# comparison statistics

# chi^2
    if (chisqr):
# compute scale factor    
        scale = numpy.nansum(weights*sp1.flux*f(sp1.wave)/vtot)/ \
            numpy.nansum(weights*f(sp1.wave)*f(sp1.wave)/vtot)
# correct variance
        vtot = numpy.nanmax([sp1.variance,sp1.variance+v(sp1.wave)*scale**2],axis=0)
        stat = numpy.nansum(weights*(sp1.flux-f(sp1.wave)*scale)**2/vtot)

# standard deviation
    elif (stddev):
# compute scale factor    
        scale = numpy.nansum(weights*sp1.flux*f(sp1.wave))/ \
            numpy.nansum(weights*f(sp1.wave)*f(sp1.wave))
# correct variance
        vtot = numpy.nanmax([sp1.variance,sp1.variance+v(sp1.wave)*scale**2],axis=0)
        stat = numpy.nansum(weights*(sp1.flux-f(sp1.wave)*scale)**2)

# absolute deviation
    elif (absdev):
# compute scale factor    
        scale = numpy.nansum(weights*sp1.flux)/ \
            numpy.nansum(weights*f(sp1.wave))
# correct variance
        vtot = numpy.nanmax([sp1.variance,sp1.variance+v(sp1.wave)*scale**2],axis=0)
        stat = numpy.nansum(weights*abs(sp1.flux-f(sp1.wave)*scale))

# error
    else:
        print 'Error: statistic for compareSpectra not given'
        return numpy.nan, numpy.nan

    return numpy.nanmax([stat,minreturn]), scale
        

def coordinateToDesignation(c):
    '''Convert RA, Dec into designation string'''
# input is ICRS
    if isinstance(c,SkyCoord):
        cc = c
    else:
        cc = properCoordinates(c)
# input is [RA,Dec] pair in degrees
    return string.replace('J{0}{1}'.format(cc.ra.to_string(unit=u.hour, sep='', precision=2, pad=True), \
        cc.dec.to_string(unit=u.degree, sep='', precision=1, alwayssign=True, pad=True)),'.','')


def dateToCaldate(d):
    '''Convert from numeric date to calendar date'''
    return d[:4]+' '+months[int(d[5:6])-1]+' '+d[-2:]


def designationToCoordinate(value, **kwargs):
    '''Convert a designation into a RA, Dec tuple or ICRS'''
    icrsflag = kwargs.get('icrs',True)

    a = re.sub('[j.:hms]','',value.lower())
    fact = 1.
    spl = a.split('+')
    if len(spl) == 1:
        spl = a.split('-')
        fact = -1.
    ra = 15.*float(spl[0][0:2])
    if (len(spl[0]) > 2):
        ra+=15.*float(spl[0][2:4])/60.
    if (len(spl[0]) > 4):
        ra+=15.*float(spl[0][4:6])/3600.
    if (len(spl[0]) > 6):
        ra+=15.*float(spl[0][6:8])/360000.
    dec = float(spl[1][0:2])
    if (len(spl[0]) > 2):
        dec+=float(spl[1][2:4])/60.
    if (len(spl[0]) > 4):
        dec+=float(spl[1][4:6])/3600.
    if (len(spl[1]) > 6): 
        dec+=float(spl[1][6:8])/360000.
    dec = dec*fact
    if icrsflag:
        return SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    else:
        return [ra,dec]


def designationToShortName(value):
    '''Produce a shortened version of designation'''
    if isinstance(value,str):
        a = re.sub('[j.:hms]','',value.lower())
        mrk = '+'
        spl = a.split(mrk)
        if len(spl) == 1:
            mrk = '-'
            spl = a.split(mrk)
        if len(spl) == 2:
            return 'J'+spl[0][0:4]+mrk+spl[1][0:4]
        else:
            return value
    else:
        raise ValueError('\nMust provide a string value for designation\n\n')

        

def fetchDatabase(*args, **kwargs):    
    '''Get the SpeX Database from either online repository or local drive'''
    dataFile = kwargs.get('dataFile','db_spexprism.txt')
    folder = kwargs.get('folder','~/')
    url = kwargs.get('url',SPLAT_URL+'/Databases/')
    local = kwargs.get('local',False)

# check if online
    local = local or (not checkOnline())
        
# first try online
    if not local:
        try:
            open(os.path.basename(dataFile), 'wb').write(urllib2.urlopen(url+dataFile).read())
            data = ascii.read(dataFile, delimiter='\t',fill_values='-99.')
            os.remove(os.path.basename(dataFile))
#            return data
        except urllib2.URLError:
            sys.stderr.write('\nReading local '+dataFile+'\n\n')
            local = True
    
# now try local drive - NOT WORKING
    file = dataFile
    if local:
        if (os.path.exists(file) == False):
            file = folder+os.path.basename(file)
        if (os.path.exists(file) == False):
            raise NameError('\nCould not find '+dataFile+' or '+file+' locally\n\n')
        else:
            data = ascii.read(file, delimiter='    ')

# clean up blanks and convert numerical values to numbers
    data['ra'][numpy.where(data['ra'] == '')] = '0.'
#    data['ran'] = [float(x) for x in data['ra']]
    data['dec'][numpy.where(data['dec'] == '')] = '0.'
#    data['decn'] = [float(x) for x in data['dec']]
    data['jmag'][numpy.where(data['jmag'] == '')] = '99.'
#    data['jmagn'] = [float(x) for x in data['jmag']]
    data['hmag'][numpy.where(data['hmag'] == '')] = '99.'
#    data['hmagn'] = [float(x) for x in data['hmag']]
    data['kmag'][numpy.where(data['kmag'] == '')] = '99.'
#    data['kmagn'] = [float(x) for x in data['kmag']]
    data['jmag_error'][numpy.where(data['jmag_error'] == '')] = '99.'
#    data['jmag_errorn'] = [float(x) for x in data['jmag_error']]
    data['hmag_error'][numpy.where(data['hmag_error'] == '')] = '99.'
#    data['hmag_errorn'] = [float(x) for x in data['hmag_error']]
    data['kmag_error'][numpy.where(data['kmag_error'] == '')] = '99.'
#    data['kmag_errorn'] = [float(x) for x in data['kmag_error']]
    data['resolution'][numpy.where(data['resolution'] == '')] = '120'
#    data['resolutionn'] = [float(x) for x in data['resolution']]
    data['airmass'][numpy.where(data['airmass'] == '')] = '1.'
#    data['airmassn'] = [float(x) for x in data['airmass']]
    data['median_snr'][numpy.where(data['median_snr'] == '')] = '0'
#    data['median_snrn'] = [float(x) for x in data['median_snr']]

# convert coordinates to SkyCoord format
#    data['skycoords'] = data['ra']
    s = []
    for i in numpy.arange(len(data['ra'])):
        try:        # to deal with a blank string
            s.append(SkyCoord(ra=float(data['ra'][i])*u.degree,dec=float(data['dec'][i])*u.degree,frame='icrs'))
        except:
            s.append(SkyCoord(ra=0.*u.degree,dec=0.*u.degree,frame='icrs'))
    data['skycoords'] = s

# add in RA/Dec (TEMPORARY)
#    ra = []
#    dec = []
#    for x in data['designation']:
#        c = designationToCoordinate(x,ICRS=False)
#        ra.append(c[0])
#        dec.append(c[1])
#    data['ra'] = ra
#    data['dec'] = dec

    data['young'] = ['young' in x for x in data['library']]
    data['subdwarf'] = ['subdwarf' in x for x in data['library']]
    data['binary'] = ['binary' in x for x in data['library']]
    data['spbinary'] = ['sbinary' in x for x in data['library']]
    data['blue'] = ['blue' in x for x in data['library']]
    data['red'] = ['red' in x for x in data['library']]
    data['giant'] = ['giant' in x for x in data['library']]
    data['wd'] = ['wd' in x for x in data['library']]
    data['standard'] = ['std' in x for x in data['library']]

# add in shortnames
    data['shortname'] = [designationToShortName(x) for x in data['designation']]

    return data


def filenameToNameDate(filename):
    '''Extract from a SPLAT filename the source name and observation date'''
    ind = filename.rfind('.')
    base = filename[:ind]
    spl = base.split('_')
    if (len(spl) < 2):
        return '', ''
    else:
        name = spl[-2]
        d = spl[-1]
        try:
            float(d)
            date = '20'+d
        except ValueError:
#            print filename+' does not contain a date'
            date = ''
        
        return name, date



def filterMag(sp,filter,*args,**kwargs):
    '''Compute the magnitudes of a spectrum given a filter name'''
# keyword parameters
    filterFolder = kwargs.get('filterFolder',SPLAT_URL+'Filters/')
    vegaFile = kwargs.get('vegaFile','vega_kurucz.txt')
    info = kwargs.get('info',False)
    units = kwargs.get('units',sp.funit)
    custom = kwargs.get('custom',False)
    vega = kwargs.get('vega',True)
    ab = kwargs.get('ab',False)
    photons = kwargs.get('photons',False)
    energy = kwargs.get('energy',False)

# filter file assignments
    filters = { \
        '2MASS J': {'file': 'j_2mass.txt', 'description': '2MASS J-band'}, \
        '2MASS H': {'file': 'h_2mass.txt', 'description': '2MASS H-band'}, \
        '2MASS Ks': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band'}, \
        'MKO J': {'file': 'j_atm_mko.txt', 'description': 'MKO J-band + atmosphere'}, \
        'MKO H': {'file': 'h_atm_mko.txt', 'description': 'MKO H-band + atmosphere'}, \
        'MKO K': {'file': 'k_atm_mko.txt', 'description': 'MKO K-band + atmosphere'}, \
        'MKO Ks': {'file': 'k_atm_mko.txt', 'description': 'MKO Ks-band'} \
        }

# check that requested filter is in list
    if (filter not in filters.keys()):
        print 'Filter '+filter+' not included in filterMag'
        info = True
        
# print out what's available
    if (info):
        print 'Filter names:'
        for x in filters.keys():
            print x+': '+filters[x]['description']    
        return numpy.nan

# convert units to erg/cm2/s/um if needed - TO BE DONE

# Read in filter
    if (custom == False):
        fwave,ftrans = numpy.genfromtxt(filterFolder+filters[filter]['file'], comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
    else:
        fwave,ftrans = custom[0],custom[1]
    fnu = const.c.to('micron/s').value/fwave
    fwave = fwave[~numpy.isnan(ftrans)]            # temporary fix
    ftrans = ftrans[~numpy.isnan(ftrans)]
                
# interpolate spectrum and vega spectrum onto filter wavelength function
    d = interp1d(sp.wave,sp.flux,bounds_error=False,fill_value=0.)
    result = numpy.nan
    
# compute relevant quantity
    if (vega):
# Read in Vega spectrum
        vwave,vtrans = numpy.genfromtxt(filterFolder+vegaFile, comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
        vwave = vwave[~numpy.isnan(vtrans)]            # temporary fix
        vtrans = vtrans[~numpy.isnan(vtrans)]
        v = interp1d(vwave,vtrans,bounds_error=False,fill_value=0.)
        result = -2.5*numpy.log10(trapz(ftrans*d(fwave),fwave)/trapz(ftrans*v(fwave),fwave))
    if (energy):
        result = trapz(ftrans*d(fwave),fwave)
    if (photons):
        convert = const.h.to('erg s')*const.c.to('micron/s')
        result = trapz(ftrans*d(fwave)*fwave,fwave) / convert.value
    if (ab):
# NOTE: THIS IS CURRENTLY NOT WORKING
        sp.flamToFnu()
        dd = interp1d(sp.nu,sp.flux,bounds_error=False,fill_value=0.)
        a = trapz(ftrans*dd(fnu),numpy.log10(fnu))
        b = trapz(ftrans,numpy.log10(fnu))
        result = -2.5*numpy.log10(a/b)-48.6
        sp.fnuToFlam()
    
    return result



def getSourceKey(*args, **kwargs):
    '''Search the SpeX database to extract the key reference for that Source'''

    pass
    return


def getSpectrum(*args, **kwargs):
    '''Get specific spectra from online library'''

    result = []
    kwargs['output'] = 'all'
    kwargs['access'] = checkAccess()
    list = kwargs.get('list',False)
    search = searchLibrary(*args, **kwargs)

# only grab published data if no proprietary access
    if checkAccess() == False:
        search = search[:][numpy.where(search['public'] == 'Y')]
    
    files = search['data_file']
        
# return just the filenames
    if (list):
        return files
    
    if len(files) > 0:
        if (len(files) == 1):
            print '\nRetrieving {} file\n'.format(len(files))
        else:
            print '\nRetrieving {} files\n'.format(len(files))
        for x in files:
            result.append(loadSpectrum(x))
    else:
        if checkAccess() == False:
            sys.stderr.write('\nNo published files match search criteria\n\n')
        else:
            sys.stderr.write('\nNo files match search criteria\n\n')
    return result
        


# simple number checker
def isNumber(s):
    '''check if something is a number'''
    try:
        float(s)
        return True
    except ValueError:
        return False


def loadInterpolatedModel(*args,**kwargs):
# interpolates within a 2x2 grid of teff and logg ONLY
    set = kwargs.get('set','btsettl08')
    kwargs['set'] = set
    teff = float(kwargs.get('teff',1000.0))
    kwargs['teff'] = teff
    logg = float(kwargs.get('logg',5.0))
    kwargs['logg'] = logg
    url = kwargs.get('folder',SPLAT_URL+'/Models/')
    kwargs['url'] = url
#    local = kwargs.get('local',False)
    kwargs['model'] = True

# first get model parameters
    param = loadModelParameters(set)
#    print param
    
# check if input parameters are within range
    for t in param.colnames:
        flg = kwargs.get(t,False)
        if flg != False:
            if (flg < param[t][0] or flg > param[t][1]):
                raise ValueError('\n\nValue {}={} outside model range'.format(t,str(flg)))

# try simply loading model (on grid)
#    try:
#        return loadModel(**kwargs)
#    except:
#        sys.stderr.write('\nCalculating interpolated model\n')

# identify grid points around input parameters
# NOTE: THIS IS JUST FOR TEFF AND LOGG GRIDDING FOR NOW

    t = teff - teff%param['teff'][2]
    trng = [max(param['teff'][0],t),min(t+param['teff'][2],param['teff'][1])]
    g = logg - logg%param['logg'][2]
    grng = [max(param['logg'][0],g),min(g+param['logg'][2],param['logg'][1])]
    x,y = numpy.meshgrid(trng,grng)

    mkwargs = kwargs.copy()
#    try:
#    print trng, grng
    mkwargs['teff'] = trng[0]
    mkwargs['logg'] = grng[0]
    md11 = loadModel(**mkwargs)
    mkwargs['teff'] = trng[1]
    md21 = loadModel(**mkwargs)
    mkwargs['teff'] = trng[0]
    mkwargs['logg'] = grng[1]
    md12 = loadModel(**mkwargs)
    mkwargs['teff'] = trng[1]
    md22 = loadModel(**mkwargs)
    
#    except:
#        raise NameError('\nProblem loading grid of models in loadInterpolatedModel\n\n')

    mflx = numpy.zeros(len(md11.wave))
    val = numpy.array([[None]*2]*2)
    for i,w in enumerate(md11.wave):
        val = numpy.array([ \
            [numpy.log10(md11.flux[i]),numpy.log10(md21.flux[i])], \
            [numpy.log10(md12.flux[i]),numpy.log10(md22.flux[i])]])
        mflx[i] = 10.**(griddata((x.flatten(),y.flatten()),val.flatten(),(teff,logg),'linear'))
    
    return Spectrum(wave=md11.wave,flux=mflx,**kwargs)


def loadModel(*args, **kwargs):
    '''load up a model spectrum based on parameters'''
# keyword parameters
    set = kwargs.get('set','btsettl08')
    teff = float(kwargs.get('teff',1000.0))
    logg = float(kwargs.get('logg',5.0))
    z = kwargs.get('z',0.0)
    kzz = kwargs.get('kzz',0)
    fsed = kwargs.get('fsed',0)
    folder = kwargs.get('folder','')
    url = kwargs.get('folder',SPLAT_URL+'/Models/')
    local = kwargs.get('local',False)
    kwargs['model'] = True
    fileFlag = False
    
# check if online
    local = local or (not checkOnline())

# a filename has been passed - simply read this file
    if (len(args) > 0):
        kwargs['filename'] = args[0]
        fileFlag = True

# determine model set
    else:

# get model parameters
        param = loadModelParameters(set)

        if (set == 'btsettl08'):
            url = url+'/'+set+'/'
#            print teff, logg
            kwargs['filename'] = set+'_{:.0f}_{:.1f}_-0.0_nc_nc_eq_0.5.txt'.format(teff,logg)
#            kwargs['filename'] = 'lte'+'{:5.3f}'.format(teff/100000.)[2:]+'-'+str(logg)[0:3]+'-0.0.BT-Settl.7_r120.txt'        
        else: 
            raise NameError('\nCurrently only have BTSettl models\n\n')

        if (teff < param['teff'][0] or teff > param['teff'][1]):
            raise NameError('\nInput Teff out of range\n\n')

        if (logg < param['logg'][0] or logg > param['logg'][1]):
            raise NameError('\nInput logg out of range\n\n')

    
# if model is on grid read in that model
    if (teff in numpy.arange(param['teff'][0],param['teff'][1]+param['teff'][2],param['teff'][2]) and \
        logg in numpy.arange(param['logg'][0],param['logg'][1]+param['logg'][2],param['logg'][2])):
        fileFlag = True

# simple read in a file
    if fileFlag:
    
# first try online
        if not local:
            try:
                open(os.path.basename(kwargs['filename']), 'wb').write(urllib2.urlopen(url+kwargs['filename']).read()) 
                kwargs['filename'] = os.path.basename(kwargs['filename'])
                sp = Spectrum(**kwargs)
#                os.close(kwargs['filename'])
                os.remove(kwargs['filename'])
                return sp
#fh, filename = tempfile.mkstemp()    do i need this?
#os.close(fh)
#os.remove(filename)
            except urllib2.URLError:
                sys.stderr.write('\nCould not find model file '+kwargs['filename']+' on SPLAT website\n\n')
                os.remove(os.path.basename(kwargs['filename']))
                local = True


    
    # now try local drive
        if (os.path.exists(kwargs['filename']) == False):
            kwargs['filename'] = folder+os.path.basename(kwargs['filename'])
            if (os.path.exists(kwargs['filename']) == False):
                raise NameError('\nCould not find '+kwargs['filename']+' locally\n\n')
            else:
                return Spectrum(**kwargs)
        else:
            return Spectrum(**kwargs)

# or do an interpolated Model
    else:
        return loadInterpolatedModel(**kwargs)
        


def loadModelParameters(*args, **kwargs):
    '''Load up Model Parameters based on help file'''
# keyword parameters
    set = args[0]
    pfile = kwargs.get('parameterFile','parameters.txt')
    
# check if online
    if not checkOnline():
        raise urllib2.URLError('\n\nCannot access online models\n\n')

    url = kwargs.get('folder',SPLAT_URL+'/Models/'+set+'/')
    try:
        open(os.path.basename(pfile), 'wb').write(urllib2.urlopen(url+pfile).read())
        parameters = ascii.read(pfile)
        os.remove(os.path.basename(pfile))
    except urllib2.URLError:
        raise urllib2.URLError('\n\nCannot access models from SPLAT website\n\n')
    
    return parameters



def loadSpectrum(*args, **kwargs):
    '''load up a SpeX spectrum based name, shortname and/or date'''
# keyword parameters
#    name = kwargs.get('name','')
#    shname = kwargs.get('shname','')
#    date = kwargs.get('date','')
    local = kwargs.get('local',False)
    folder = kwargs.get('folder','')
    url = kwargs.get('folder',SPLAT_URL+'/Spectra/')
    kwargs['model'] = False
    file = kwargs.get('filename','')
    file = kwargs.get('file',file)
    if (len(args) > 0):
        file = args[0]
    kwargs['filename'] = file

# check if online
    local = local or (not checkOnline())
    
# CODE NEEDED HERE TO SET UP FILE NAME; FOR NOW JUST ERROR
    
# a filename has been passed
    if (kwargs['filename'] == ''):
        raise NameError('\nNeed to pass in filename for spectral data')

# first try online
    if not local:
        try:
            open(os.path.basename(kwargs['filename']), 'wb').write(urllib2.urlopen(url+kwargs['filename']).read())
            kwargs['filename'] = os.path.basename(kwargs['filename'])
            sp = Spectrum(**kwargs)
            os.remove(os.path.basename(kwargs['filename']))
            return sp
        except urllib2.URLError:
            sys.stderr.write('\nCould not find data file '+kwargs['filename']+' at '+url+'\n\n')
            local = True
    
# now try local drive
    if (os.path.exists(kwargs['filename']) == False):
        kwargs['filename'] = folder+os.path.basename(kwargs['filename'])
        if (os.path.exists(kwargs['filename']) == False):
            raise NameError('\nCould not find '+kwargs['filename']+' locally\n\n')
        else:
            return Spectrum(**kwargs)
    else:
        return Spectrum(**kwargs)


# code to measure a defined index from a spectrum using Monte Carlo noise estimate
# measure method can be mean, median, integrate
# index method can be ratio = 1/2, valley = 1-2/3, OTHERS
# output is index value and uncertainty
def measureIndex(sp,*args,**kwargs):
    '''Measure an index on a spectrum based on defined methodology'''

# keyword parameters
    method = kwargs.get('method','ratio')
    sample = kwargs.get('sample','integrate')
    nsamples = kwargs.get('nsamples',100)
    noiseFlag = kwargs.get('nonoise',False)
            
# create interpolation functions
    w = numpy.where(numpy.isnan(sp.flux) == False)
    f = interp1d(sp.wave[w],sp.flux[w],bounds_error=False,fill_value=0.)
    w = numpy.where(numpy.isnan(sp.noise) == False)
    if (numpy.size(w) != 0):
        s = interp1d(sp.wave[w],sp.noise[w],bounds_error=False,fill_value=numpy.nan)
        noiseFlag = False
    else:
        s = interp1d(sp.wave,sp.noise,bounds_error=False,fill_value=numpy.nan)
        noiseFlag = True
                
# error checking on number of arguments provided
    if (len(args) < 2):
        print 'measureIndex needs at least two samples to function'
        return numpy.nan, numpy.nan
    elif (len(args) < 3 and (method == 'line' or method == 'allers')):
        print method+' requires at least 3 sample regions'
        return numpy.nan, numpy.nan

# define the sample vectors
    values = numpy.zeros((len(args),nsamples))

# loop over all sampling regions
    for i,waveRng in enumerate(args):
        xNum = (numpy.arange(0,nsamples+1.0)/nsamples)* \
            (numpy.nanmax(waveRng)-numpy.nanmin(waveRng))+numpy.nanmin(waveRng)
        yNum = f(xNum)
        yNum_e = s(xNum)

# now do MonteCarlo measurement of value and uncertainty
        for j in numpy.arange(0,nsamples):

# doing this for noise = nan
            if (numpy.isnan(yNum_e[0]) == False):
                yVar = numpy.random.normal(yNum,yNum_e)
            else:
                yVar = yNum

# choose function for measuring indices
            if (sample == 'integrate'):
                values[i,j] = trapz(yVar,xNum)
            elif (sample == 'average'):
                values[i,j] = numpy.nanmean(yVar)
            elif (sample == 'median'):
                values[i,j] = numpy.median(yVar)
            elif (sample == 'maximum'):
                values[i,j] = numpy.nanmax(yVar)
            elif (sample == 'minimum'):
                values[i,j] = numpy.nanmin(yVar)
            else:
                values[i,j] = numpy.nanmean(yVar)

# compute index based on defined method
# default is a simple ratio
    if (method == 'ratio'):
        vals = values[0,:]/values[1,:]
    elif (method == 'line'):
        vals = (values[0,:]+values[1,:])/values[2,:]
    elif (method == 'change'):
        vals = 2.*(values[0,:]-values[1,:])/(values[0,:]+values[1,:])
    elif (method == 'allers'):
        vals = (((numpy.mean(args[0])-numpy.mean(args[1]))/(numpy.mean(args[2])-numpy.mean(args[1])))*values[2,:] \
            + ((numpy.mean(args[2])-numpy.mean(args[0]))/(numpy.mean(args[2])-numpy.mean(args[1])))*values[1,:]) \
            /values[0,:]
    else:
        vals = values[0,:]/values[1,:]
            
# output mean, standard deviation
    if (noiseFlag):
        return numpy.mean(vals), numpy.nan
    else:
        return numpy.mean(vals), numpy.std(vals)


# wrapper function for measuring specific sets of indices

def measureIndexSet(sp,**kwargs):

# keyword parameters
    set = kwargs.get('set','burgasser')

# determine combine method
    if (set.lower() == 'burgasser'):
        reference = 'Indices from Burgasser et al. (2006)'
        names = ['H2O-J','CH4-J','H2O-H','CH4-H','H2O-K','CH4-K','K-J']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.14,1.165],[1.26,1.285],method='ratio',sample='integrate',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[1.315,1.335],[1.26,1.285],method='ratio',sample='integrate',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[1.48,1.52],[1.56,1.60],method='ratio',sample='integrate',**kwargs)
        inds[3],errs[3] = measureIndex(sp,[1.635,1.675],[1.56,1.60],method='ratio',sample='integrate',**kwargs)
        inds[4],errs[4] = measureIndex(sp,[1.975,1.995],[2.08,2.12],method='ratio',sample='integrate',**kwargs)
        inds[5],errs[5] = measureIndex(sp,[2.215,2.255],[2.08,2.12],method='ratio',sample='integrate',**kwargs)
        inds[6],errs[6] = measureIndex(sp,[2.06,2.10],[1.25,1.29],method='ratio',sample='integrate',**kwargs)
    elif (set.lower() == 'tokunaga'):
        reference = 'Indices from Tokunaga & Kobayashi (1999)'
        names = ['K1','K2']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[2.1,2.18],[1.96,2.04],method='change',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[2.2,2.28],[2.1,2.18],method='change',sample='average',**kwargs)
    elif (set.lower() == 'reid'):
        reference = 'Indices from Reid et al. (2001)'
        names = ['H2O-A','H2O-B']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.33,1.35],[1.28,1.30],method='ratio',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[1.47,1.49],[1.59,1.61],method='ratio',sample='average',**kwargs)
    elif (set.lower() == 'geballe'):
        reference = 'Indices from Geballe et al. (2002)'
        names = ['H2O-1.2','H2O-1.5','CH4-2.2']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.26,1.29],[1.13,1.16],method='ratio',sample='integrate',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[1.57,1.59],[1.46,1.48],method='ratio',sample='integrate',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[2.08,2.12],[2.215,2.255],method='ratio',sample='integrate',**kwargs)
    elif (set.lower() == 'allers'):
        reference = 'Indices from Allers et al. (2007), Allers & Liu (2013)'
        names = ['H2O','FeH-z','VO-z','FeH-J','KI-J','H-cont']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.55,1.56],[1.492,1.502],method='ratio',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[0.99135,1.00465],[0.97335,0.98665],[1.01535,1.02865],method='allers',sample='average',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[1.05095,1.06505],[1.02795,1.04205],[1.07995,1.09405],method='allers',sample='average',**kwargs)
        inds[3],errs[3] = measureIndex(sp,[1.19880,1.20120],[1.19080,1.19320],[1.20680,1.20920],method='allers',sample='average',**kwargs)
        inds[4],errs[4] = measureIndex(sp,[1.23570,1.25230],[1.21170,1.22830],[1.26170,1.27830],method='allers',sample='average',**kwargs)
        inds[5],errs[5] = measureIndex(sp,[1.54960,1.57040],[1.45960,1.48040],[1.65960,1.68040],method='allers',sample='average',**kwargs)
    elif (set.lower() == 'slesnick'):
        reference = 'Indices from Slesnick et al. (2004)'
        names = ['H2O-1','H2O-2','FeH']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.335,1.345],[1.295,1.304],method='ratio',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[2.035,2.045],[2.145,2.155],method='ratio',sample='average',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[1.1935,1.2065],[1.2235,1.2365],method='ratio',sample='average',**kwargs)
    elif (set.lower() == 'mclean'):
        reference = 'Indices from McLean et al. (2003)'
        names = ['H2OD']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.951,1.977],[2.062,2.088],method='ratio',sample='average',**kwargs)
    else:
        print '{} is not one of the sets used for measureIndexSet'.format(set)
        return numpy.nan

# output dictionary of indices
    result = {names[i]: (inds[i],errs[i]) for i in numpy.arange(len(names))}
#    result['reference'] = reference
#    return inds,errs,names

    return result



# To do:
#     masking telluric regions
#    labeling features
#     add legend
def plotSpectrum(*args, **kwargs):

# error check - make sure you're plotting something
    if (len(args) < 1):
        print 'plotSpectrum needs at least on Spectrum object to plot'
        return

# this loop solves issues if a list of spectrum objects is provided accidentally
    sp = []
    for i,a in enumerate(args):
        if (type(a) == list):
            sp.append(a[0])
        else:
            sp.append(a)

# keyword parameters
    title = kwargs.get('title','')
    zeropoint = kwargs.get('zeropoint',[0. for x in range(len(sp))])
    xlabel = kwargs.get('xlabel','{} ({})'.format(sp[0].wlabel,sp[0].wunit))
    ylabel = kwargs.get('ylabel','{} {} ({})'.format(sp[0].fscale,sp[0].flabel,sp[0].funit))
    xrange = kwargs.get('xrange',sp[0].waveRange())
    bound = xrange
    ymax = [s.fluxMax() for s in sp]
    yrange = kwargs.get('yrange',[0,numpy.nanmax(ymax)+numpy.nanmax(zeropoint)])
    bound.extend(yrange)
    grid = kwargs.get('grid',False)
    colors = kwargs.get('colors',['k' for x in range(len(sp))])
    if (len(colors) < len(args)):
        colors.extend(['k' for x in range(len(sp)-len(colors))])
    colorsUnc = kwargs.get('colors',['k' for x in range(len(sp))])
    if (len(colorsUnc) < len(sp)):
        colorsUnc.extend(['k' for x in range(len(sp)-len(colorsUnc))])
    linestyle = kwargs.get('linestyle',['steps' for x in range(len(sp))])
    if (len(linestyle) < len(sp)):
        linestyle.extend(['steps' for x in range(len(sp)-len(linestyle))])
    filename = kwargs.get('filename','')
    format = kwargs.get('format',filename.split('.')[-1])
    showNoise = kwargs.get('showNoise',[False for x in range(len(sp))])
    if not isinstance(showNoise, list):
        showNoise = [showNoise]
    if (len(showNoise) < len(sp)):
        showNoise.extend([True for x in range(len(sp)-len(showNoise))])
    showZero = kwargs.get('showZero',[False for x in numpy.arange(len(sp))])
    if not isinstance(showZero, list):
        showZero = [showZero]
    if (len(showZero) < len(sp)):
        showZero.extend([True for x in range(len(sp)-len(showZero))])
#    mask = kwargs.get('mask',False)                # not yet implemented
#    labels = kwargs.get('labels','')            # not yet implemented
#    features = kwargs.get('features','')        # not yet implemented


#    plt.clf()
# loop through sources
    plt.subplots(1)
    for ii,a in enumerate(sp):
        flx = [i+zeropoint[ii] for i in a.flux]
        plt.plot(a.wave,flx,color=colors[ii],linestyle=linestyle[ii])
# show noise
        if (showNoise[ii]):
            ns = [i+zeropoint[ii] for i in a.noise]
            plt.plot(a.wave,ns,color=colorsUnc[ii],linestyle=linestyle[ii],alpha=0.3)
# zeropoint
        if (showZero[ii]):
            ze = numpy.ones(len(a.flux))*zeropoint[ii]
            plt.plot(a.wave,ze,color=colors[ii],linestyle=':',alpha=0.3)

# grid
    if (grid):
        plt.grid()
# labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axis(bound)
    plt.title(title)
    
# save to filen or display
    if (len(filename) > 0): 
        plt.savefig(filename, format=format)
    
    else:
        plt.show()
        plt.ion()        # make window interactive by default
    
    return


def properCoordinates(c):
    '''converts various coordinate forms to the proper SkyCoord format'''
    if isinstance(c,SkyCoord):
        return c
    elif isinstance(c,list):
        return SkyCoord(c[0]*u.deg,c[1]*u.deg,frame='icrs')
# input is sexigessimal string
    elif isinstance(c,str):
        return SkyCoord(c,'icrs', unit=(u.hourangle, u.deg))
    else:
        raise ValueError('\nCould not parse input format\n\n')


def readSpectrum(**kwargs):

# TO BE DONE:
# FIX IF THERE IS NO NOISE CHANNEL
# PRODUCE AND RETURN HEADER => CHANGE OUTPUT TO DICTIONARY?
    
# keyword parameters
    folder = kwargs.get('folder','./')
    catchSN = kwargs.get('catchSN',True)
#    model = kwargs.get('model',False)
#    uncertainty = kwargs.get('uncertainty',not model)
    file = kwargs.get('filename','')
    if (os.path.exists(file) == False):
        file = folder+os.path.basename(file)
    if (os.path.exists(file) == False):
        raise NameError('\nCould not find ' + file+'\n\n')
    
# determine which type of file
    ftype = file.split('.')[-1]

# fits file    
    if (ftype == 'fit' or ftype == 'fits'):
        data = fits.open(file)
        if 'NAXIS3' in data[0].header.keys():
            d = data[0].data[0,:,:]
        else:
            d = data[0].data
        data.close()

# ascii file    
    else:
        try:
            d = numpy.genfromtxt(file, comments='#', unpack=False, \
                missing_values = ('NaN','nan'), filling_values = (numpy.nan)).transpose()
        except ValueError:
            d = numpy.genfromtxt(file, comments=';', unpack=False, \
                 missing_values = ('NaN','nan'), filling_values = (numpy.nan)).transpose()

# assign arrays to wave, flux, noise
    wave = d[0,:]
    flux = d[1,:]
    if len(d[:,0]) > 2:
        noise = d[2,:]
    else:
        noise = numpy.zeros(len(flux))
        noise[:] = numpy.nan        

# fix places where noise is claimed to be 0
    w = numpy.where(noise == 0.)
    noise[w] = numpy.nan

# fix nans in flux
    w = numpy.where(numpy.isnan(flux) == True)
    flux[w] = 0.

# fix to catch badly formatted files where noise column is S/N                 
#    print flux, numpy.median(flux)
    if (catchSN):
          w = numpy.where(flux > stats.nanmedian(flux))
          if (stats.nanmedian(flux[w]/noise[w]) < 1.):
              noise = flux/noise
              w = numpy.where(numpy.isnan(noise))
              noise[w] = stats.nanmedian(noise)

    return wave, flux, noise


# THIS HAS BEEN SUPERCEDED BY NEW FXN; KEEPING AROUND JUST IN CASE
def searchLibrary_OLD(*args, **kwargs):
    '''Search the SpeX database to extract the key reference for that Spectrum
        Note that this is currently only and AND search - need to figure out
        how to a full SQL style search'''

# get database
    data = fetchDatabase(**kwargs)
    sql = Table()
    ref = kwargs.get('output','data_file')

# search parameters
    if kwargs.get('name',False) != False:
        nm = kwargs['name']
        if isinstance(nm,str):
            nm = [nm]
        sql['name'] = nm
    if kwargs.get('designation',False) != False:
        desig = kwargs['designation']
        if isinstance(desig,str):
            desig = [desig]
        sql['designation'] = desig
    if kwargs.get('shortname',False) != False:
        sname = kwargs['shortname']
        if isinstance(sname,str):
            sname = [sname]
        for i,sn in enumerate(sname):
            if sn[0].lower() != 'j':
                sname[i] = 'J'+sname[i]
        sql['shortname'] = sname
    if kwargs.get('date',False) != False:
        sql['observation_date'] = [kwargs['date']]
    if kwargs.get('young',False) != False:
        sql['young'] = [True]
    if kwargs.get('subdwarf',False) != False:
        sql['subdwarf'] = [True]
    if kwargs.get('binary',False) != False:
        sql['binary'] = [True]
    if kwargs.get('spbin',False) != False:
        sql['spbin'] = [True]
    if kwargs.get('red',False) != False:
        sql['red'] = [True]
    if kwargs.get('blue',False) != False:
        sql['blue'] = [True]

# NEED TO ADD IN SEARCH IN AREA, MAGNITUDE RANGE, OBS DATE RANGE, REFERENCES

    if len(sql) > 0:
        result = join(data,sql)
    else:
        result = data
        
    return result[ref]



def searchLibrary(*args, **kwargs):
    '''Search the SpeX database to extract the key reference for that Spectrum
        Note that this is currently only and AND search - need to figure out
        how to a full SQL style search'''
# program parameters
    ref = kwargs.get('output','all')
    radius = kwargs.get('radius',10.)      # search radius in arcseconds
    classes = ['young','subdwarf','binary','spbinary','red','blue','giant','wd','standard']

# get database
    data = fetchDatabase(**kwargs)
    if (ref not in data.colnames and ref != 'all'):
        print '\nWarning: searchLibrary cannot return unknown column {}; returning filename instead\n\n'.format(ref)
        ref = 'all'
    data['select'] = numpy.zeros(len(data['ra']))
    count = 0.
    
# search by name
    if kwargs.get('name',False) != False:
        nm = kwargs['name']
        if isinstance(nm,str):
            nm = [nm]
        for n in nm:
            data['select'][numpy.where(data['name'] == n)] += 1
        count+=1.
# search by shortname
    if kwargs.get('shortname',False) != False:
        sname = kwargs['shortname']
        if isinstance(sname,str):
            sname = [sname]
        for sn in sname:
            if sn[0].lower() != 'j':
                sn = 'J'+sn
            data['select'][numpy.where(data['shortname'] == sn)] += 1
        count+=1.
# search by reference list
    if kwargs.get('reference',False) != False:
        refer = kwargs['reference']
        if isinstance(ref,str):
            refer = [refer]
        for r in refer:
            data['select'][numpy.where(data['data_reference'] == r)] += 1
        count+=1.
# search by designation
    if kwargs.get('designation',False) != False:
        desig = kwargs['designation']
        if isinstance(desig,str):
            desig = [desig]
        for d in desig:
            data['select'][numpy.where(data['designation'] == d)] += 1
        count+=1.
# search by observation date range
    if kwargs.get('date',False) != False:
        date = kwargs['date']
        if isinstance(date,str) or isinstance(date,long) or isinstance(date,float) or isinstance(date,int):
            date = [float(date),float(date)]
        elif isinstance(date,list):
            date = [float(date[0]),float(date[-1])]
        else:
            raise ValueError('\nCould not parse date input {}\n\n'.format(date))
        data['daten'] = [float(x) for x in data['observation_date']]
        data['select'][numpy.where(numpy.logical_and(data['daten'] >= date[0],data['daten'] <= date[1]))] += 1
        count+=1.
# search by coordinate - NOTE: THIS IS VERY SLOW RIGHT NOW
    if kwargs.get('coordinate',False) != False:
        coord = kwargs['coordinate']
        if isinstance(coord,SkyCoord):
            cc = coord
        else:
            cc = properCoordinates(coord)
        data['separation'] = [cc.separation(data['skycoords'][i]).arcsecond for i in numpy.arange(len(data['skycoords']))]
        data['select'][numpy.where(data['separation'] <= radius)] += 1
        count+=1.
# search by spectral type
    if kwargs.get('spt',False) != False:
        spt = kwargs['spt']
        if not isinstance(spt,list):        # one value = only this type
            spt = [spt,spt]
        if isinstance(spt[0],str):          # convert to numerical spt
            spt = [typeToNum(spt[0]),typeToNum(spt[1])]
        data['sptn'] = [typeToNum(x) for x in data['spex_type']]
        data['select'][numpy.where(numpy.logical_and(data['sptn'] >= spt[0],data['sptn'] <= spt[1]))] += 1
        count+=1.
    if kwargs.get('spex_spt',False) != False:
        spt = kwargs['spex_spt']
        if not isinstance(spt,list):        # one value = only this type
            spt = [spt,spt]
        if isinstance(spt[0],str):          # convert to numerical spt
            spt = [typeToNum(spt[0]),typeToNum(spt[1])]
        data['sptn'] = [typeToNum(x) for x in data['spex_type']]
        data['select'][numpy.where(numpy.logical_and(data['sptn'] >= spt[0],data['sptn'] <= spt[1]))] += 1
        count+=1.
    if kwargs.get('nir_spt',False) != False:
        spt = kwargs['nir_spt']
        if not isinstance(spt,list):        # one value = only this type
            spt = [spt,spt]
        if isinstance(spt[0],str):          # convert to numerical spt
            spt = [typeToNum(spt[0]),typeToNum(spt[1])]
        data['sptn'] = [typeToNum(x) for x in data['nir_type']]
        data['select'][numpy.where(numpy.logical_and(data['sptn'] >= spt[0],data['sptn'] <= spt[1]))] += 1
        count+=1.
    if kwargs.get('opt_spt',False) != False:
        spt = kwargs['opt_spt']
        if not isinstance(spt,list):        # one value = only this type
            spt = [spt,spt]
        if isinstance(spt[0],str):          # convert to numerical spt
            spt = [typeToNum(spt[0]),typeToNum(spt[1])]
        data['sptn'] = [typeToNum(x) for x in data['opt_type']]
        data['select'][numpy.where(numpy.logical_and(data['sptn'] >= spt[0],data['sptn'] <= spt[1]))] += 1
        count+=1.
# search by magnitude range
    if kwargs.get('jmag',False) != False:
        mag = kwargs['jmag']
        if not isinstance(mag,list):        # one value = faint limit
            mag = [0,mag]
        data['magn'] = [float(x) for x in data['jmag']]
        data['select'][numpy.where(numpy.logical_and(data['magn'] >= mag[0],data['magn'] <= mag[1]))] += 1
        count+=1.
    if kwargs.get('hmag',False) != False:
        mag = kwargs['hmag']
        if not isinstance(mag,list):        # one value = faint limit
            mag = [0,mag]
        data['magn'] = [float(x) for x in data['hmag']]
        data['select'][numpy.where(numpy.logical_and(data['magn'] >= mag[0],data['magn'] <= mag[1]))] += 1
        count+=1.
    if kwargs.get('kmag',False) != False:
        mag = kwargs['kmag']
        if not isinstance(mag,list):        # one value = faint limit
            mag = [0,mag]
        data['magn'] = [float(x) for x in data['kmag']]
        data['select'][numpy.where(numpy.logical_and(data['magn'] >= mag[0],data['magn'] <= mag[1]))] += 1
        count+=1.
# search by S/N range
    if kwargs.get('snr',False) != False:
        snr = kwargs['snr']
        if not isinstance(snr,list):        # one value = minimum S/N
            snr = [float(snr),1.e9]
        data['snrn'] = [float(x) for x in data['median_snr']]
        data['select'][numpy.where(numpy.logical_and(data['snrn'] >= snr[0],data['snrn'] <= snr[1]))] += 1
        count+=1.
# search by class
    for c in classes:
        if kwargs.get(c,'n') != 'n':
            test = kwargs.get(c)
            if isinstance(test,bool):
                data['select'][numpy.where(data[c] == test)]+=1
                count+=1.

# limit access to public data for most users
    if checkAccess() == True:
        data['select'][numpy.where(data['public'] != 'Y')] = 0

# logic of search
    logic = 'and'         # default combination
    logic = kwargs.get('combine',logic).lower()
    logic = kwargs.get('logic',logic).lower()
    if (logic == 'and'):
        data['select'] = numpy.floor(data['select']/count)
    elif (logic == 'or'):
        data['select'] = numpy.ceil(data['select']/count)
    else:
        raise NameError('\nDo not recognize logical operation {}; use logic = and/or\n\n'.format(logic))

# return sorted by ra by default
    if kwargs.get('sort',True) != False:
        data.sort('ra')
    if (ref == 'all'):
        return data[:][numpy.where(data['select']==1)]
    else:
        return data[ref][numpy.where(data['select']==1)]


def test():
    test_src = 'J1503+2525'

    sys.stderr.write('\n\n>>>>>>>>>>>> TESTING SPLAT CODE <<<<<<<<<<<<\n')
# check you are online
    if checkOnline():
        sys.stderr.write('\n...you are online and can see SPLAT website\n')
    else:
        sys.stderr.write('\n...you are NOT online or cannot see SPLAT website\n')

# check your access
    if checkAccess():
        sys.stderr.write('\n...you currently HAVE access to unpublished spectra\n')
    else:
        sys.stderr.write('\n...you currently DO NOT HAVE access to unpublished spectra\n')
        
# check you can search for and load a spectrum
    sp = getSpectrum(shortname=test_src)[0]
    sp.info()
    sys.stderr.write('\n...getSpectrum and loadSpectrum successful\n')

# check searchLibrary
    list = searchLibrary(young=True,output='data_file')
    sys.stderr.write('\n{} young spectra in the SPL  ...searchLibrary successful\n'.format(len(list)))

# check index measurement
    ind = measureIndexSet(sp,set='burgasser')
    sys.stderr.write('\nSpectral indices for '+test_src+':\n')
    for k in ind.keys():
        print '\t{:s}: {:4.3f}+/-{:4.3f}'.format(k, ind[k][0], ind[k][1])
    sys.stderr.write('...index measurement successful\n')

# check classification
    spt, spt_e = classifyByIndex(sp,set='burgasser')
    sys.stderr.write('\n...index classification of '+test_src+' = {:s}+/-{:2.1f}; successful\n'.format(spt,spt_e))

    spt, spt_e = classifyByStandard(sp,method='kirkpatrick')
    sys.stderr.write('\n...standard classification of '+test_src+' = {:s}+/-{:2.1f}; successful\n'.format(spt,spt_e))

# check SpT -> Teff
    teff, teff_e = typeToTeff(spt,unc=spt_e)
    sys.stderr.write('\n...Teff of '+test_src+' = {:.1f}+/-{:.1f} K; successful\n'.format(teff,teff_e))

# check flux calibration
    sp.fluxCalibrate('2MASS J',15.0)
    mag = filterMag(sp,'MKO J')
    sys.stderr.write('\n...apparent magnitude MKO J = {:3.2f} from 2MASS J = 15.0; filter calibration successful\n'.format(mag))

# check models
    mdl = loadModel(teff=teff,logg=5.0,set='btsettl08')
    sys.stderr.write('\n...interpolated model generation successful\n')

# check normalization
    mdl.normalize()
    sp.normalize()
    sys.stderr.write('\n...normalization successful\n')

# check compareSpectrum
    chi, scale = compareSpectra(sp,mdl,mask_standard=True)
    sys.stderr.write('\nFit to model: chi^2 = {}, scale = {}'.format(chi,scale))
    sys.stderr.write('\n...compareSpectra successful\n'.format(chi,scale))

# check plotSpectrum
    mdl.scale(scale)
    plotSpectrum(sp,mdl,colors=['k','r'],title='If this appears everything is OK: close window')
    sys.stderr.write('\n...plotSpectrum successful\n')
    sys.stderr.write('\n>>>>>>>>>>>> SPLAT TEST SUCCESSFUL; HAVE FUN! <<<<<<<<<<<<\n\n')



def typeToNum(input, **kwargs):
    '''convert between string and numeric spectral types'''
# keywords     
    error = kwargs.get('error','')
    unc = kwargs.get('uncertainty',0.)
    subclass = kwargs.get('subclass','')
    lumclass = kwargs.get('lumclass','')
    ageclass = kwargs.get('ageclass','')
    colorclass = kwargs.get('colorclass','')
    peculiar = kwargs.get('peculiar',False)
    spletter = 'KMLTY'

# number -> spectral type
    if (isNumber(input)):
        spind = int(abs(input/10))
        spdec = numpy.around(input,1)-spind*10.
        pstr = ''
        if (unc > 1.):
            error = ':'
        if (unc > 2.):
            error = '::'
        if (peculiar):
            pstr = 'p'
        if (0 <= spind < len(spletter)):
            output = colorclass+subclass+spletter[spind]+'{:3.1f}'.format(spdec)+ageclass+lumclass+pstr+error
        else:
            print 'Spectral type number must be between 0 ({}0) and {} ({}9)'.format(spletter[0],len(spletter)*10.-1.,spletter[-1])
            output = numpy.nan

# spectral type -> number
    else:
        input = string.split(input,sep='+')[0]    # remove +/- sides
        input = string.split(input,sep='-')[0]    # remove +/- sides
        sptype = re.findall('[{}]'.format(spletter),input)
        if (len(sptype) == 1):
            output = spletter.find(sptype[0])*10.
            spind = input.find(sptype[0])+1
            if (spind < len(input)):
                if (input.find('.') < 0):
                    if (isNumber(input[spind])):
                        output = output+float(input[spind])
                else:
                    output = output+float(input[spind:spind+3])
                    spind = spind+3
            ytype = re.findall('[abcd]',input.split('p')[-1])
            if (len(ytype) == 1):
                ageclass = ytype[0]
            if (input.find('p') != -1):
                 peculiar = True
            if (input.find('sd') != -1):
                 subclass = 'sd'
            if (input.find('esd') != -1):
                 subclass = 'esd'
            if (input.find('usd') != -1):
                 subclass = 'usd'
            if (input.count('I') > 0):
                 lumclass = ''.join(re.findall('I',input))
            if (input.count(':') > 0):
                 error = ''.join(re.findall(':',input))
            if (input[0] == 'b' or input[0] == 'r'):
                 colorclass = input[0]
        if (len(sptype) != 1):
#            print 'Only spectral classes {} are handled with this routine'.format(spletter)
            output = numpy.nan
    return output


def typeToTeff(input, **kwargs):
    '''return Teff for a given SpT'''
# keywords     
    nsamples = kwargs.get('nsamples',100)
    unc = kwargs.get('uncertainty',0.)
    unc = kwargs.get('unc',unc)
    unc = kwargs.get('spt_e',unc)
    ref = kwargs.get('ref','looper2008')
    ref = kwargs.get('set',ref)
    ref = kwargs.get('method',ref)

# convert spectral type string to number
    if (type(input) == str):
        spt = typeToNum(input,uncertainty=unc)
    else:
        spt = copy.deepcopy(input)

# choose among possible options
# Stephens et al. (2009, ApJ, 702, 1545); using OPT/IR relation for M6-T8
# plus alternate coefficients for L3-T8
    if (ref.lower() == 'stephens2009'):
        reference = 'Teff/SpT relation from Stephens et al. (2009)'
        sptoffset = 10.
        coeff = [-0.0025492,0.17667,-4.4727,54.67,-467.26,4400.]
        range = [16.,38.]
        fitunc = 100.
        coeff_alt = [-0.011997,1.2315,-50.472,1031.9,-10560.,44898.]
        range_alt = [23.,38.]

# Looper et al. (2008, ApJ, 685, 1183)
    elif (ref.lower() == 'looper2008'):
        reference = 'Teff/SpT relation from Looper et al. (2008)'
        sptoffset = 20.
        coeff = [9.084e-4,-4.255e-2,6.414e-1,-3.101,1.950,-108.094,2319.92]
        range = [20.,38.]
        fitunc = 87.
# Looper is default
    else:
        sys.stderr.write('\nInvalid Teff/SpT relation given ({})\n'.format(ref))
        return numpy.nan, numpy.nan

    if (range[0] <= spt <= range[1]):
        vals = numpy.polyval(coeff,numpy.random.normal(spt-sptoffset,unc,nsamples))
        if (ref.lower() == 'stephens' and (range_alt[0] <= spt <= range_alt[1])):
            vals = numpy.polyval(coeff_alt,numpy.random.normal(spt-sptoffset,unc,nsamples))
        teff = numpy.nanmean(vals)
        teff_e = (numpy.nanstd(vals)**2+fitunc**2)**0.5
        return teff, teff_e
    else:
        sys.stderr.write('\nSpectral Type is out of range for {:s} Teff/SpT relation\n'.format(reference))
        return numpy.nan, numpy.nan


def weightedMeanVar(vals, winput, *args, **kwargs):
    '''Compute weighted mean of an array of values through various methods'''
    
    method = kwargs.get('method','')
    minwt = kwargs.get('weight_minimum',0.)
    dof = kwargs.get('dof',len(vals)-1)
    if (numpy.nansum(winput) <= 0.):
        weights = numpy.ones(len(vals))
    
# uncertainty weighting: input is unceratinties   
    if (method == 'uncertainty'):
        weights = [w**(-2) for w in winput]
# ftest weighting: input is chisq values, extra dof value is required   
    elif (method == 'ftest'):
# fix issue of chi^2 = 0
        minchi = numpy.nanmin(winput)
        weights = numpy.array([stats.f.pdf(w/minchi,dof,dof) for w in winput])
# just use the input as the weights
    else:
        weights = [w for w in winput]
    
    weights = weights/numpy.nanmax(weights)
    weights[numpy.where(weights < minwt)] = 0.
    mn = numpy.nansum(vals*weights)/numpy.nansum(weights)
    var = numpy.nansum(weights*(vals-mn)**2)/numpy.nansum(weights)
    
    return mn,var        

