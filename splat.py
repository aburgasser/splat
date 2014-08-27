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

# local application/library specific imports
import bdevopar # still in development

# suppress warnings - probably not an entirely safe approach!
numpy.seterr(all='ignore')
warnings.simplefilter("ignore")
#from splat._version import __version__

# CONSTANTS 
SPLAT_URL = 'http://pono.ucsd.edu/~adam/splat/'
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
parameter_names = ['teff','logg','z','fsed','kzz']
spex_pixel_scale = 0.15            # spatial scale in arcseconds per pixel
spex_wave_range = [0.65,2.45]    # default wavelength range
max_snr = 1000.0                # maximum S/N ratio permitted

spex_stdfiles = { \
    'M0.0': 'spex_prism_Gl270_091203.fits',\
    'M1.0': 'spex_prism_Gl424_091229.fits',\
    'M2.0': 'spex_prism_Gl91_081012.fits',\
    'M3.0': 'spex_prism_Gl752A_070704.fits',\
    'M4.0': 'spex_prism_Gl213_071110dl.fits',\
    'M5.0': 'spex_prism_Gl51_070728.fits',\
    'M6.0': 'spex_prism_LHS1375_081012.fits',\
    'M7.0': 'spex_prism_vB8_070704.fits',\
    'M8.0': 'spex_prism_vB10_070704.fits',\
    'M9.0': 'spex_prism_LHS2924_070704.fits',\
    'L0.0': 'spex_prism_0345+2540_kc.txt',\
    'L1.0': 'spex_prism_2130-0845_080713.txt',\
    'L2.0': 'spex_prism_kelu-1_060411.txt',\
    'L3.0': 'spex_prism_1506+1321_060410.txt',\
    'L4.0': 'spex_prism_2158-1550_060901.fits',\
    'L5.0': 'spex_prism_sdss0835+19_chiu06.txt',\
    'L6.0': 'spex_prism_1010-0406_kc.txt',\
    'L7.0': 'spex_prism_0103+1935_060902.fits',\
    'L8.0': 'spex_prism_1632+1904_030905.txt',\
    'L9.0': 'spex_prism_denis0255-4700_040908.txt',\
    'T0.0': 'spex_prism_1207+0244_061221dl.txt',\
    'T1.0': 'spex_prism_0837-0000_061221.fits',\
    'T2.0': 'spex_prism_sdss1254-0122_030522.txt',\
    'T3.0': 'spex_prism_1209-1004_030523.fits',\
    'T4.0': 'spex_prism_2254+3123_030918.txt',\
    'T5.0': 'spex_prism_1503+2525_030522.txt',\
    'T6.0': 'spex_prism_sdss1624+0029_040312.txt',\
    'T7.0': 'spex_prism_0727+1710_040310.txt',\
    'T8.0': 'spex_prism_0415-0935_030917.txt',\
    'T9.0': 'spex_prism_0722-0540_110404.fits'}

        
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
        self.header = kwargs.get('header',Table())
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
#                self.wave, self.flux, self.noise = readSpectrum(**kwargs)
                rs = readSpectrum(**kwargs)
                self.wave = rs['wave']
                self.flux = rs['flux']
                self.noise = rs['noise']
                self.header = rs['header']
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
# load in header from splat database
#        self.header = kwargs.get('header',searchLibrary(file=self.filename))

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
    '''
    :Purpose: ``Convert from numeric date to calendar date, and vice-versa.``
    :param d: ``A numeric date of the format '20050412', or a date in the 
                calendar format '2005 Jun 12'``
    :Example:
       >>> import splat
       >>> caldate = splat.dateToCaldate('20050612')
       >>> print caldate
       2005 Jun 12
       >>> date = splat.caldateToDate('2005 June 12')
       >>> print date
       20050612
    '''
    return d[:4]+str((months.index(d[5:8])+1)/100.)[2:4]+d[-2:]


def checkFile(filename,**kwargs):
    '''
    :Purpose: ``Checks if a spectrum file exists in the SPLAT's library.``
    :param filename: ``A string containing the spectrum's filename.``
    :Example: 
       >>> import splat
       >>> spectrum1 = 'spex_prism_1315+2334_110404.fits'
       >>> print splat.checkFile(spectrum1)
       True
       >>> spectrum2 = 'fake_name.fits'
       >>> print splat.checkFile(spectrum2)
       False
    '''
    url = kwargs.get('url',SPLAT_URL)+'/Spectra/'
    flag = checkOnline()
    if (flag):
        try:
            open(os.path.basename(filename), 'wb').write(urllib2.urlopen(url+filename).read())
        except urllib2.URLError:
            flag = False
    return flag


def checkAccess(**kwargs):
    access_file = '.splat_access'
    result = False

    try:
        home = os.environ.get('HOME')
        if home == None:
            home = './'
        bcode = urllib2.urlopen(SPLAT_URL+access_file).read()
        lcode = base64.b64encode(open(home+'/'+access_file,'r').read())
        if (bcode in lcode):        # changed to partial because of EOL variations
            result = True
    except:
        result = False

    if (kwargs.get('report','') != ''):
        if result == True:
            print 'You have full access to all SPLAT data'
        else:
            print 'You have access only to published data'
    return result
    

def checkOnline():
    '''
    :Purpose: ``Checks if SPLAT's URL is accessible from your machine--
                that is, checks if you and the host are online.``
    :Example:
       >>> import splat
       >>> splat.checkOnline()
       True  # SPLAT's URL was detected.
       >>> splat.checkOnline()
       False # SPLAT's URL was not detected.
    '''
    try:
        urllib2.urlopen(SPLAT_URL)
        return True
    except urllib2.URLError:
        return False



def classifyByIndex(sp, *args, **kwargs):
    '''
    :Purpose: ``Determine the spectral type and uncertainty for a spectrum 
                based on indices. Makes use of published index-SpT relations
                from Reid et al. (2001); Testi et al. (2001); Allers et al. 
                (2007); and Burgasser (2007). Returns 2-element tuple 
                containing spectral type (numeric or string) and 
                uncertainty.``
    :param sp: ``Spectrum class object, which should contain wave, flux and 
                 noise array elements.``
    :param \**kwargs (optional): - ``set = 'burgasser': named set of indices to measure and compute spectral type:``
                          * ``'allers': H2O from Allers et al.``
                          * ``'burgasser': H2O-J, CH4-J, H2O-H, CH4-H, 
                            CH4-K from Burgasser (2007)``
                          * ``'reid':H2O-A and H2O-B from Reid et al.(2001)``
                          * ``'testi': sHJ, sKJ, sH2O_J, sH2O_H1, sH2O_H2, 
                            sH2O_K from Testi et al. (2001)``
                      - ``string = False: return spectral type as a string 
                        (uses typeToNum)``
                      - ``round = False: rounds off to nearest 0.5 subtypes``
                      - ``remeasure = True: force remeasurement of indices``
                      - ``nsamples = 100: number of Monte Carlo samples for 
                        error computation``
                      - ``nloop = 5: number of testing loops to see if 
                        spectral type is within a certain range``

    :Example:
       >>> import splat
       >>> spc = splat.loadSpectrum('spex_prism_gl570d_030522.txt')
       >>> print splat.classifyByIndex(spc, string=True, set='burgasser', round=True)
       ('T7.5', 0.25285169510990341)
    '''
    
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

    verbose = kwargs.get('verbose',True)
    method = kwargs.get('method','')
    best_flag = kwargs.get('best',False)
    sptrange = kwargs.get('sptrange',[10,39])
    if (isinstance(sptrange[0],str) != False):
        sptrange = [typeToNum(sptrange[0]),typeToNum(sptrange[1])]
    unc_sys = 0.5

# if you just want to compare to one standard
    cspt = kwargs.get('compareto',False)
    if (cspt != False):
        if (isinstance(cspt,str) == False):
            cspt = typeToNum(cspt)
# round off
        cspt = typeToNum(numpy.round(typeToNum(cspt)))
        mkwargs = copy.deepcopy(kwargs)
        mkwargs['compareto'] = False
        mkwargs['sptrange'] =[cspt,cspt]
        return classifyByStandard(sp,**mkwargs)
            
#    stdsptnum = numpy.arange(30)+10.
#    stdsptstr = ['M'+str(n) for n in numpy.arange(10)]
#    stdsptstr.append(['L'+str(n) for n in numpy.arange(10)])
#    stdsptstr.append(['T'+str(n) for n in numpy.arange(10)])

    if (method == 'kirkpatrick'):
        comprng = [0.9,1.4]         # as prescribed in Kirkpatrick et al. 2010, ApJS, 
    else:
        comprng = [0.7,2.45]        # by default, compare whole spectrum

# compute fitting statistics
    stat = []
    sspt = []
    sfile = []
    for t in numpy.arange(sptrange[0],sptrange[1]+1):
        spstd = loadSpectrum(file=spex_stdfiles[typeToNum(t)])
        chisq,scale = compareSpectra(sp,spstd,fit_ranges=[comprng],chisqr=True,novar2=True)
        stat.append(chisq)
        sspt.append(t)
        sfile.append(spex_stdfiles[typeToNum(t)])
        if (verbose):
            print spex_stdfiles[typeToNum(t)], typeToNum(t), chisq, scale

# list of sorted standard files and spectral types
    sorted_stdfiles = [x for (y,x) in sorted(zip(stat,sfile))]
    sorted_stdsptnum = [x for (y,x) in sorted(zip(stat,sspt))]

# select either best match or an ftest-weighted average
    if (best_flag or len(stat) == 1):
        sptn = sorted_stdsptnum[0]
        sptn_e = unc_sys
    else:
        mean,var = weightedMeanVar(sspt,stat,method='ftest',dof=sp.dof)
        if (var**0.5 < 1.):
            sptn = numpy.round(mean*2)*0.5
        else:
            sptn = numpy.round(mean)
        sptn_e = (unc_sys**2+var)**0.5

# string or not?
    if (kwargs.get('string', True) == True):
        spt = typeToNum(sptn,uncertainty=sptn_e)
    else:
        spt = sptn

# plot spectrum compared to best spectrum
    if (kwargs.get('plot',False) != False):
        spc = sp.copy()     # copy to avoid change original value
        spc.normalize()
        spstd = loadSpectrum(file=sorted_stdfiles[0])
        chisq,scale = compareSpectra(spc,spstd,fit_ranges=[comprng],chisqr=True)
        spstd.scale(scale)
        plotSpectrum(spc,spstd,colors=['k','r'],\
            title=sp.name+' vs '+typeToNum(sorted_stdsptnum[0])+' Standard',**kwargs)

    return spt, sptn_e
    

def classifyByTemplate(sp, *args, **kwargs):
    '''Classify a spectrum by comparing to other spectra in the library'''
    spt_return = kwargs.get('spt_return','spex')
    verbose = kwargs.get('verbose',True)
    set = kwargs.get('set','')
    unc_sys = 0.5
    if (kwargs.get('method','') == 'kirkpatrick'):
        comprng = [0.9,1.4]         # as prescribed in Kirkpatrick et al. 2010, ApJS, 
    else:
        comprng = [0.7,2.45]        # by default, compare whole spectrum

#  canned searches
#  constrain spectral types
    spt = [10.,39.9]
    if ('m dwarf' in set):
        spt = [10,19.9]
    if ('l dwarf' in set):
        spt = [20,29.9]
    if ('t dwarf' in set):
        spt = [30,39.9]
    if ('vlm' in set):
        spt = [numpy.max([17,spt[0]]),numpy.min([39.9,spt[-1]])]
    spt = kwargs.get('spt',spt)
    
#  constrain S/N
    snr = 0.
    if ('high sn' in set):
        snr = 100.
    snr = kwargs.get('snr',snr)
    
#  don't compare to same spectrum
    try:
        excludefile = [sp.filename]
    except:
        excludefile = ['']
    excludefile = kwargs.get('excludefile',excludefile)
        
    if ('optical' in set):
        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,opt_type=spt,giant=False,logic='and')
    elif ('standard' in set):
        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,standard=True,logic='and')
    elif ('companion' in set):
        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,companion=True,giant=False,logic='and')
    elif ('young' in set):
        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,young=True,logic='and')
    elif ('subdwarf' in set):
        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,subdwarf=True,logic='and')
    elif ('single' in set):
        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,spbinary=False,binary=False,giant=False,logic='and')
    elif ('spectral binaries' in set):
        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,spbinary=True,logic='and')
    elif (set != ''):
        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt)
    else:
        lib = searchLibrary(output='all',excludefile=excludefile,**kwargs)
        
# first search for the spectra desired - parameters are set by user
    files = lib['data_file']

# which spectral type to return
    if ('spex' in spt_return):
        sptref = 'spex_type'
    elif ('opt' in spt_return):
        sptref = 'opt_type'
    elif ('nir' in spt_return):
        sptref = 'nir_type'
    else:
        sptref = 'lit_type'

    lib = lib[:][numpy.where(lib[sptref] != '')]
    files = lib['data_file']
    sspt = [typeToNum(s) for s in lib[sptref]]

    if len(files) == 0:
        print '\nNo templates available for comparison\n\n'
        return numpy.nan, numpy.nan
    if (verbose):
        print '\nComparing to {} templates\n\n'.format(len(files))
    if len(files) > 100:
        if (verbose):
            print 'This may take some time!\n\n'.format(len(files))

# do comparison
    stat = []
    for i,f in enumerate(files):
        s = loadSpectrum(file=f)
        chisq,scale = compareSpectra(sp,s,fit_ranges=[comprng],chisqr=True,novar2=True)
        stat.append(chisq)
        if (verbose):
            print f, typeToNum(sspt[i]), chisq, scale
        
# list of sorted standard files and spectral types
    sorted_files = [x for (y,x) in sorted(zip(stat,files))]
    sorted_spt = [x for (y,x) in sorted(zip(stat,sspt))]

# select either best match or an ftest-weighted average
    if (kwargs.get('best',False) or len(stat) == 1):
        sptn = sorted_spt[0]
        sptn_e = unc_sys
    else:
        mean,var = weightedMeanVar(sspt,stat,method='ftest',dof=sp.dof)
        if (var**0.5 < 1.):
            sptn = numpy.round(mean*2.)*0.5
        else:
            sptn = numpy.round(mean)
        sptn_e = (unc_sys**2+var)**0.5

# plot spectrum compared to best spectrum
    if (kwargs.get('plot',False) != False):
        s = loadSpectrum(file=sorted_files[0])
        chisq,scale = compareSpectra(sp,s,fit_ranges=[comprng],chisqr=True)
        s.scale(scale)
        plotSpectrum(sp,s,colors=['k','r'],title=sp.name+' vs '+s.name,**kwargs)

# string or not?
    if (kwargs.get('string', True) == True):
        spt = typeToNum(sptn,uncertainty=sptn_e)
    else:
        spt = sptn

    return spt, sptn_e
 

def classifyGravity(sp, *args, **kwargs):
    '''Determine the gravity classification of a brown dwarf
    using the method of Allers & Liu 2013'''
    
    verbose = kwargs.get('verbose',False)
    
# Chart for determining gravity scores based on gravity sensitive 
# indices as described in the Allers and Liu paper.
# The key to the overall indices dictionary is each index name.
# The key to each index dictionary are the spectral types, which 
# contain the limits for each gravity score.
# To access a value do the following: print grav['FeH-z']['M5'][0] 
# which should return 'nan'
    grav = {\
        'FeH-z':{'M5.0':[numpy.nan,numpy.nan],'M6.0':[1.068,1.039],'M7.0':[1.103,1.056],'M8.0':[1.146,1.074],'M9.0': [1.167,1.086],'L0.0': [1.204,1.106],'L1.0':[1.252,1.121],'L2.0':[1.298,1.142],'L3.0': [1.357,1.163],'L4.0': [1.370,1.164],'L5.0': [1.258,1.138],'L6.0': [numpy.nan,numpy.nan],'L7.0': [numpy.nan,numpy.nan]},\
        'VO-z': {'M5.0':[numpy.nan,numpy.nan],'M6.0':[numpy.nan,numpy.nan],'M7.0': [numpy.nan,numpy.nan],'M8.0': [numpy.nan,numpy.nan],'M9.0': [numpy.nan,numpy.nan],'L0.0': [1.122,1.256],'L1.0': [1.112,1.251],'L2.0': [1.110,1.232],'L3.0': [1.097,1.187],'L4.0': [1.073,1.118],'L5.0': [numpy.nan,numpy.nan],'L6.0': [numpy.nan,numpy.nan],'L7.0': [numpy.nan,numpy.nan]},\
        'KI-J': {'M5.0': [numpy.nan,numpy.nan], 'M6.0': [1.042,1.028], 'M7.0': [1.059,1.036],'M8.0': [1.077,1.046],'M9.0': [1.085,1.053],'L0.0': [1.098,1.061],'L1.0': [1.114,1.067],'L2.0': [1.133,1.073],'L3.0': [1.135,1.075],'L4.0': [1.126,1.072],'L5.0': [1.094,1.061],'L6.0': [numpy.nan,numpy.nan],'L7.0': [numpy.nan,numpy.nan]},\
        'H-cont': {'M5.0': [numpy.nan,numpy.nan], 'M6.0': [.988,.994], 'M7.0': [.981,.990],'M8.0': [.963,.984],'M9.0': [.949,.979],'L0.0': [.935,.972],'L1.0': [.914,.968],'L2.0': [.906,.964],'L3.0': [.898,.960],'L4.0': [.885,.954],'L5.0': [.869,.949],'L6.0': [.874,.950],'L7.0': [numpy.nan,numpy.nan]}}

# Calculate Allers indices and their uncertainties 
    ind = kwargs.get('indices',False)
    if ind == False:
        ind = measureIndexSet(sp,set='allers')

# Determine the object's NIR spectral type and its uncertainty 
    sptn = kwargs.get('spt',False)
    if sptn == False:
        sptn, spt_e = classifyByIndex(sp,string=False,set='allers')
        if numpy.isnan(sptn):
            print 'Spectral type could not be determined from indices'
            return numpy.nan
    if isinstance(sptn,str):
        sptn = typeToNum(sptn)
    Spt = typeToNum(numpy.floor(sptn))

#Check whether the NIR SpT is within gravity sensitive range values
    if ((sptn < 16.0) or (sptn > 27.0)):
        print 'Spectral type '+typeToNum(sptn)+' outside range for gravity classification'
        return numpy.nan

#Creates an empty array with dimensions 4x1 to fill in later with 5 gravscore values
    gravscore = {}
    medgrav = []

# Use the spt to pick the column that contains the 
# values we want to compare our indices with. 
    for k in grav.keys():
        val = 0.0
        if k == 'VO-z' or k=='H-cont':
            if numpy.isnan(grav[k][Spt][0]):
                val = numpy.nan
            if ind[k][0] >= grav[k][Spt][0]:
                val = 1.0
            if ind[k][0] >= grav[k][Spt][1]:
                val = 2.0
            if verbose:
                print k,ind[k][0], ind[k][1], val 
        if k == 'FeH-z' or k=='KI-J':
            if numpy.isnan(grav[k][Spt][0]):
                val = numpy.nan
            if ind[k][0] <= grav[k][Spt][0]:
                val = 1.0
            if ind[k][0] <= grav[k][Spt][1]:
                val = 2.0
            if verbose:
                print k,ind[k][0], ind[k][1], val 
        gravscore[k] = val
        medgrav.append(val)

# determine median score, or mean if even					
    if (len(numpy.where(numpy.isnan(medgrav) == False))%2 == 0):
        gravscore['score'] = scipy.stats.nanmean(medgrav)        
    else:        
        gravscore['score'] = scipy.stats.nanmedian(medgrav)        

    if gravscore['score'] <= 0.5:
       gravscore['gravity_class'] = 'FLD-G'
    elif gravscore['score'] > 0.5 and gravscore['score'] < 1.5:
       gravscore['gravity_class'] = 'INT-G'
    elif gravscore['score'] >= 1.5:
       gravscore['gravity_class'] = 'VL-G'

# plot spectrum against standard
    if (kwargs.get('plot',False) != False):
        spt,unc = classifyByStandard(sp,compareto=Spt,method='kirkpatrick',**kwargs)
        
# return gravity class or entire dictionary
    if (kwargs.get('allscores',False) == False):
        return gravscore['gravity_class']
    else:
        return gravscore

    
    
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
    data['companion'] = ['companion' in x for x in data['library']]

# add in shortnames
    data['shortname'] = [designationToShortName(x) for x in data['designation']]

# create literature spt 
    data['lit_type'] = data['opt_type']
    w = numpy.where(numpy.logical_and(data['lit_type'] == '',data['nir_type'] != ''))
    data['lit_type'][w] = data['nir_type'][w]
    sptn = [typeToNum(x) for x in data['lit_type']]
    w = numpy.where(numpy.logical_and(sptn > 29.,data['nir_type'] != ''))
    data['lit_type'][w] = data['nir_type'][w]
#    w = numpy.where(numpy.logical_and(data['lit_type'] == '',typeToNum(data['spex_type']) > 17.))
#    data['lit_type'][w] = data['spex_type'][w]

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



def getSpectrum(*args, **kwargs):
    '''Get specific spectra from online library'''

    result = []
    kwargs['output'] = 'all'
    search = searchLibrary(*args, **kwargs)

# limit access to most users
    if checkAccess() == False:
        search = search[:][numpy.where(search['public'] == 'Y')]
    
    files = search['data_file']
        
# return just the filenames
    if (kwargs.get('list',False) != False):
        return files
    
    if len(files) > 0:
        if (len(files) == 1):
            print '\nRetrieving {} file\n'.format(len(files))
        else:
            print '\nRetrieving {} files\n'.format(len(files))
        for i,x in enumerate(files):
            result.append(loadSpectrum(x,header=search[i:i+1]))
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
#    header = kwargs.get('header',False)
#    kwargs['header'] = header
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
#            if kwargs['header'] == False:
#                kwargs['header'] = searchLibrary(file=kwargs['filename'])
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
        header = data[0].header
        data.close()

# ascii file    
    else:
        try:
            d = numpy.genfromtxt(file, comments='#', unpack=False, \
                missing_values = ('NaN','nan'), filling_values = (numpy.nan)).transpose()
        except ValueError:
            d = numpy.genfromtxt(file, comments=';', unpack=False, \
                 missing_values = ('NaN','nan'), filling_values = (numpy.nan)).transpose()
        header = fits.Header()
        
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

    return {'wave':wave,'flux':flux,'noise':noise,'header':header}



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
    
# search by filename
    file = kwargs.get('file','')
    file = kwargs.get('filename',file)
    if (file != ''):
        if isinstance(file,str):
            file = [file]
        for f in file:
            data['select'][numpy.where(data['data_file'] == f)] += 1
        count+=1.
# exclude by filename
    if kwargs.get('excludefile',False) != False:
        file = kwargs['excludefile']
        if isinstance(file,str):
            file = [file]
        for f in file:
            data['select'][numpy.where(data['data_file'] != f)] += 1
        count+=1.
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
# exclude by shortname
    if kwargs.get('excludesource',False) != False:
        sname = kwargs['excludesource']
        if isinstance(sname,str):
            sname = [sname]
        for sn in sname:
            if sn[0].lower() != 'j':
                sn = 'J'+sn
            data['select'][numpy.where(data['shortname'] != sn)] += 1
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
    sref = ''
    if (kwargs.get('spt',False) != False):
        sref = 'spex_type'
        spt = kwargs['spt']
    if (kwargs.get('spex_spt',False) != False):
        sref = 'spex_type'
        spt = kwargs['spex_spt']
    if (kwargs.get('spex_type',False) != False):
        sref = 'spex_type'
        spt = kwargs['spex_type']
    if (kwargs.get('opt_spt',False) != False):
        sref = 'opt_type'
        spt = kwargs['opt_spt']
    if (kwargs.get('opt_type',False) != False):
        sref = 'opt_type'
        spt = kwargs['opt_type']
    if (kwargs.get('nir_spt',False) != False):
        sref = 'nir_type'
        spt = kwargs['nir_spt']
    if (kwargs.get('nir_type',False) != False):
        sref = 'nir_type'
        spt = kwargs['nir_type']
    if sref != '':
        if not isinstance(spt,list):        # one value = only this type
            spt = [spt,spt]
        if isinstance(spt[0],str):          # convert to numerical spt
            spt = [typeToNum(spt[0]),typeToNum(spt[1])]
        data['sptn'] = [typeToNum(x) for x in data[sref]]
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
    if checkAccess() == False:
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
    test_src = 'J1507-1627'

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

    grav = classifyGravity(sp)
    sys.stderr.write('\n...gravity class of '+test_src+' = {:s}; successful\n'.format(grav))

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


def typeToMag(spt, filt, **kwargs):
    """
    Takes a spectral type, its uncertainty, and a filter, returns absolute magnitude
    Filter options = ['MKO K', 'MKO H', 'MKO J', 'MKO Y', 'MKO LP', '2MASS J', '2MASS K', '2MASS H']
    References = ['faherty', 'dupuy', 'burgasser']
    """

#Keywords
    nsamples = kwargs.get('nsamples', 100)
    ref = kwargs.get('ref', 'dupuy')
    unc = kwargs.get('unc', 0.)

#Convert spectral type string to number
  
    if (type(spt) == str):
        spt = typeToNum(spt, uncertainty=unc)
    else:
        spt = copy.deepcopy(spt)

#Faherty
    if (ref.lower() == 'faherty'):
        sptoffset = 10.
        reference = 'Abs Mag/SpT relation from Faherty et al. (2012)'
        coeffs = { \
            'MKO J': {'fitunc' : 0.30, 'range' : [20., 38.],  \
                'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, \
            'MKO H': {'fitunc' : 0.27, 'range' : [20., 38.], \
                'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, \
            'MKO K': {'fitunc' : 0.28, 'range' : [20., 38.], \
                'coeff' : [.0000816516, -.00469032, .0940816, -.485519, 9.76100]}}

# Burgasser
    elif (ref.lower() == 'burgasser'):
        sptoffset = 20.
        reference = 'Abs Mag/SpT relation from Burgasser (2007)'
        coeffs = { \
            'MKO K': {'fitunc' : 0.26, 'range' : [20., 38.], \
                'coeff': [.0000001051, -.000006985, .0001807, -.002271, .01414, -.04024, .05129, .2322, 10.45]}}
 
# Dupuy & Liu, default reference  
    elif (ref.lower() == 'dupuy'):
        reference = 'Abs Mag/SpT relation from Dupuy & Liu (2012)'
        sptoffset = 10.
        coeffs = { \
            'MKO J': {'fitunc' : 0.39, 'range' : [16., 39.], \
                'coeff' : [-.00000194920, .000227641, -.0103332, .232771, -2.74405, 16.3986, -28.3129]}, \
            'MKO Y': {'fitunc': 0.40, 'range' : [16., 39.], \
                'coeff': [-.00000252638, .000285027, -.0126151, .279438, -3.26895, 19.5444, -35.1560]}, \
            'MKO H': {'fitunc': 0.38, 'range' : [16., 39.], \
                'coeff': [-.00000224083, .000251601, -.0110960, .245209, -2.85705, 16.9138, -29.7306]}, \
            'MKO K': {'fitunc': 0.40, 'range' : [16., 39.], \
                'coeff': [-.00000104935, .000125731, -.00584342, .135177, -1.63930, 10.1248, -15.2200]}, \
            'MKO LP': {'fitunc': 0.28, 'range': [16., 39.], \
                'coeff': [0.00000, 0.00000, .0000546366, -.00293191, .0530581,  -.196584, 8.89928]}, \
            '2MASS J': {'fitunc': 0.40, 'range': [16., 39.], \
                'coeff': [-.000000784614, .000100820, -.00482973, .111715, -1.33053, 8.16362, -9.67994]}, \
            '2MASS H': {'fitunc': 0.40, 'range': [16., 38.5], \
                'coeff': [-.00000111499, .000129363, -.00580847, .129202, -1.50370, 9.00279, -11.7526]}, \
            '2MASS K': {'fitunc': 0.43, 'range':[16., 38.5], \
                'coeff': [0.00000, 0.00000, .000106693, -.00642118, .134163, -.867471, 11.0114]}}

    else:
        sys.stderr.write('\nInvalid Abs Mag/SpT relation given: %s\n' % ref)
        return numpy.nan, numpy.nan

    if (filt.upper() in coeffs.keys()) == 1:
        for filter in coeffs.keys():
            if filt.upper() == filter:
                coeff = coeffs[filter]['coeff']
                fitunc = coeffs[filter]['fitunc']
                range = coeffs[filter]['range']
    else:
        sys.stderr.write('\n Invalid filter given for %s\n' % reference)
        return numpy.nan, numpy.nan

# compute magnitude if its in the right spectral type range
    if (range[0] <= spt <= range[1]):
        if (unc > 0.):
            vals = numpy.polyval(coeff, numpy.random.normal(spt - sptoffset, unc, nsamples))
            abs_mag = numpy.nanmean(vals)
            abs_mag_error = (numpy.nanstd(vals)**2+fitunc**2)**0.5
            return abs_mag, abs_mag_error
        else:
            abs_mag = numpy.polyval(coeff, spt-sptoffset)
            return abs_mag, fitunc
    else:
        sys.stderr.write('\nSpectral Type is out of range for %s Abs Mag/SpT relation\n' % reference)
        return numpy.nan, numpy.nan


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
    unc = kwargs.get('uncertainty',0.001)
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
        if (ref.lower() == 'stephens'):
            if (range_alt[0] <= spt <= range_alt[1]):
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

