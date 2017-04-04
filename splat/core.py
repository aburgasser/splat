# -*- coding: utf-8 -*-
from __future__ import print_function

# WORKING COPY OF SPLAT CODE LIBRARY
# based on routines developed by:
#    Christian Aganze
#    Daniella Bardalez Gagliuffi
#    Jessica Birky
#    Adam Burgasser
#    Caleb Choban
#    Andrew Davis
#    Ivanna Escala
#    Aishwarya Iyer
#    Yuhui Jin
#    Michael Lopez
#    Alex Mendez
#    Gretel Mercado
#    Elizabeth Moreno
#    Jonathan Parra
#    Maitrayee Sahi
#    Adrian Suarez
#    Melisa Tallis
#    Chris Theissen
#    Tomoki Tamiya
#    Russell van Linge

# imports - internal
import copy
import os
import sys
import warnings

# imports - external
import matplotlib.pyplot as plt
import numpy
import requests
from astropy.io import ascii, fits            # for reading in spreadsheet
from astropy.table import Table, join            # for reading in table files
from astropy.coordinates import SkyCoord      # coordinate conversion
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units
from scipy import stats, signal
from scipy.integrate import trapz        # for numerical integration
from scipy.interpolate import interp1d

if sys.version_info.major != 2 and sys.version_info.major != 3:
    raise NameError('\nSPLAT only works on Python 2.7 and 3.X\n')
if sys.version_info.major == 2:     # switch for those using python 3
    import string

# splat functions and constants
from splat.initialize import *
from splat.utilities import *
import splat.citations as spbib
from splat.photometry import filterMag
#from splat.database import searchLibrary, keySpectrum

__version__ = VERSION

# holding arrays
SPECTRA_READIN = {}
STDS_DWARF_SPEX = {}
STDS_SD_SPEX = {}
STDS_ESD_SPEX = {}

# databases - using the .txt files for now, will need to change to SQL at a future date
DB_SOURCES = ascii.read(SPLAT_PATH+DB_FOLDER+DB_SOURCES_FILE)
DB_SPECTRA = ascii.read(SPLAT_PATH+DB_FOLDER+DB_SPECTRA_FILE)


# suppress warnings - probably not an entirely safe approach!
numpy.seterr(all='ignore')
warnings.simplefilter("ignore")

# temporary constants - will be removed
spex_pixel_scale = 0.15            # spatial scale in arcseconds per pixel
uspex_pixel_scale = 0.10            # spatial scale in arcseconds per pixel
spex_wave_range = [0.65,2.45]*u.micron    # default wavelength range
max_snr = 1000.0                # maximum S/N ratio permitted


#####################################################
###############   Spectrum class   ##################
#####################################################

class Spectrum(object):
    '''
    :Description: Primary class for containing spectral and source data for SpeX Prism Library.

    Optional Inputs:

    :param ismodel: Set to True to specify spectrum as a model (default = False)
    :param wlabel: label of wavelength (default = 'Wavelength')
    :param wunit: unit in which wavelength is measured (default = u.micron)
    :param wunit_label: label of the unit of wavelength (default = 'micron')
    :param flabel: label of flux density (default = 'F\_lambda')
    :param fscale: string describing how flux density is scaled (default = '')
    :param funit: unit in which flux density is measured (default = u.erg/(u.cm**2 * u.s * u.micron)
    :param funit_label: label of the unit of flux density (default = 'erg cm\^-2 s\^-1 micron\^-1')
    :param resolution: Resolution of spectrum (default = 150)
    :param slitpixelwidth: Width of the slit measured in subpixel values (default = 3.33)
    :param slitwidth: Actual width of the slit, measured in arcseconds. Default value is the ``slitpixelwidth`` multiplied an assumed (for SpeX) spectrograph pixel scale of 0.15 arcseconds 
    :param header: header info of the spectrum (default = Table())
    :param filename: a string containing the spectrum's filename (default = '')
    :param file: same as filename (default = '')
    :param idkey: spectrum key of the desired spectrum (default = False)

    :Example:
       >>> import splat
       >>> sp = splat.Spectrum(filename='myspectrum.fits')      # read in a file
       >>> sp = splat.Spectrum('myspectrum.fits')               # same
       >>> sp = splat.Spectrum(10002)                           # read in spectrum with idkey = 10002
       >>> sp = splat.Spectrum(wave=wavearray,flux=fluxarray)   # create objects with wavelength & flux arrays
    '''

    def __init__(self, *args, **kwargs):
# some presets
        sdb = False
        self.ismodel = kwargs.get('ismodel',False)
        self.wlabel = kwargs.get('wlabel',r'Wavelength')
        self.wunit = kwargs.get('wunit',u.micron)
        self.wunit_label = kwargs.get('wunit_label',r'$\mu$m')
        self.flabel = kwargs.get('flabel',r'F$_{\lambda}$')
        self.fscale = kwargs.get('fscale','')
        self.funit = kwargs.get('funit',u.erg/(u.cm**2 * u.s * u.micron))
        self.funit_label = kwargs.get('funit_label',r'erg~cm$^{-2}$~s$^{-1}$~$\mu$m$^{-1}$')
#        self.header = kwargs.get('header',fits.PrimaryHDU())
        self.header = kwargs.get('header',{})
        self.filename = kwargs.get('file','')
        self.filename = kwargs.get('filename',self.filename)
        self.idkey = kwargs.get('idkey',False)
        self.instrument = kwargs.get('instrument','SPEX_PRISM')
        self.history = []

# breakouts for specific instruments
        if kwargs.get('APOGEE') == True and kwargs.get('file',False) != False:
            rs = _readAPOGEE(kwargs['file'],**kwargs)
            self.instrument = 'APOGEE'
            self.wave = rs['wave']
            self.flux = rs['flux']
            self.noise = rs['noise']
            self.model = rs['model']
            self.header = rs['header']        
            self.history.append('Spectrum successfully loaded')
    # create a copy to store as the original
            self.original = copy.deepcopy(self)
            return


# option 1: a filename is given
        if (len(args) > 0):
            if isinstance(args[0],str):
                self.filename = args[0]
        if kwargs.get('file',self.filename) != '':
            self.filename = kwargs.get('file',self.filename)
        if kwargs.get('filename',self.filename) != '':
            self.filename = kwargs.get('filename',self.filename)

# option 2: a spectrum ID is given
        if (len(args) > 0):
            if isinstance(args[0],int):
                self.idkey = args[0]


        if self.idkey != False:
#            self.idkey = kwargs.get('idkey',self.idkey)
            try:
                sdb = keySpectrum(self.idkey)
                if sdb != False:
                    self.filename = sdb['DATA_FILE'][0]
            except:
                print('Warning: problem reading in spectral database, a known problem for Python 3.X')
        elif self.ismodel == False and self.filename != '':
            kwargs['filename']=self.filename
            kwargs['silent']=True
            try:
                t = searchLibrary(**kwargs)
                if len(t) > 0:
                    sdb = t
            except:
                print('Warning: problem reading in source or spectral database, a known problem for Python 3.X')            
        else:
            sdb = False

# set up folder - by default this is local data directory
        kwargs['folder'] = kwargs.get('folder',SPLAT_PATH+DATA_FOLDER)
        self.simplefilename = os.path.basename(self.filename)
        self.file = self.filename
        self.name = kwargs.get('name',self.simplefilename)
        kwargs['filename'] = self.filename

# option 3: wave and flux are given
        if len(kwargs.get('wave','')) > 0 and len(kwargs.get('flux','')) > 0:
            self.wave = kwargs['wave']
            self.flux = kwargs['flux']
            if len(kwargs.get('noise','')) > 0:
                self.noise = kwargs['noise']
            else:
                self.noise = numpy.array([numpy.nan for i in self.wave])

# read in data from file
        elif self.filename != '':
            try:
                rs = readSpectrum(self.filename,**kwargs)
                self.wave = rs['wave']
                self.flux = rs['flux']
                self.noise = rs['noise']
                self.header = rs['header']
            except:
                raise NameError('\nCould not load spectral file {:s}, recheck the filename and path'.format(kwargs.get('filename','')))

# empty spectrum vessel (used for copying)
        else:
            self.wave = []
            self.flux = []
            self.noise = []
            print ('Warning: Creating an empty Spectrum object')

# process spectral data
        if len(self.wave) > 0:
# convert to numpy arrays
            self.wave = numpy.array(self.wave)
            self.flux = numpy.array(self.flux)
            self.noise = numpy.array(self.noise)
# enforce positivity and non-nan
            if (numpy.nanmin(self.flux) < 0):
                self.flux[numpy.where(self.flux < 0)] = 0.
            self.flux[numpy.isnan(self.flux)] = 0.
# check on noise being too low
            if (numpy.nanmax(self.flux/self.noise) > max_snr):
                self.noise[numpy.where(self.flux/self.noise > max_snr)]=numpy.median(self.noise)
# convert to astropy quantities with units
# assuming input is flam in erg/s/cm2/micron
            if ~isinstance(self.wave,u.quantity.Quantity):
                self.wave = numpy.array(self.wave)*self.wunit
            if ~isinstance(self.flux,u.quantity.Quantity):
                self.flux = numpy.array(self.flux)*self.funit
            if ~isinstance(self.wave,u.quantity.Quantity):
                self.noise = numpy.array(self.noise)*self.funit
# some conversions
            self.flam = self.flux
            self.nu = self.wave.to('Hz',equivalencies=u.spectral())
            self.fnu = self.flux.to('Jy',equivalencies=u.spectral_density(self.wave))
            self.noisenu = self.noise.to('Jy',equivalencies=u.spectral_density(self.wave))
            self.fnu_unit = u.Jansky
# calculate variance
            self.variance = self.noise**2
# signal to noise
            self.snr = self.computeSN()
#            self.wave_original = copy.deepcopy(self.wave)
#            self.flux_original = copy.deepcopy(self.flux)
#            self.noise_original = copy.deepcopy(self.noise)
#            self.variance_original = copy.deepcopy(self.variance)
#            self.resolution_original = copy.deepcopy(self.resolution)
#            self.slitpixelwidth_original = copy.deepcopy(self.slitpixelwidth)


# populate information on source and spectrum from database
        if sdb != False:
            for k in sdb.keys():
                setattr(self,k.lower(),str(sdb[k][0]))
            self.shortname = designationToShortName(self.designation)
            self.date = str(self.observation_date)
# convert some data into numbers
            kconv = ['ra','dec','julian_date','median_snr','resolution','airmass',\
            'jmag','jmag_error','hmag','hmag_error','kmag','kmag_error','source_key']
            for k in kconv:
                try:
                    setattr(self,k,float(getattr(self,k)))
                except:
                    setattr(self,k,numpy.nan)
# this is to make sure the database resolution is the default value
            if kwargs.get('resolution',False) == False:
                kwargs['resolution'] = self.resolution
            if kwargs.get('instrument',False) == False:
                kwargs['resolution'] = self.resolution

# instrument specific information
        hkys = list(self.header.keys())

# automated stuff for spex data
        if 'INSTRUME' in hkys:
            if 'spex' in self.header['INSTRUME'].lower() and 'GRAT' in hkys:
                if 'lowres15' in self.header['GRAT'].lower() or 'prism' in self.header['GRAT'].lower(): self.instrument = 'SPEX_PRISM'
                if 'shortxd' in self.header['GRAT'].lower() or 'sxd' in self.header['GRAT'].lower(): self.instrument = 'SPEX_SXD'
        if 'INSTR' in hkys:
            if 'spex' in self.header['INSTR'].lower() and 'GRAT' in hkys:
                if 'lowres15' in self.header['GRAT'].lower() or 'prism' in self.header['GRAT'].lower(): self.instrument = 'SPEX_PRISM'
                if 'shortxd' in self.header['GRAT'].lower() or 'sxd' in self.header['GRAT'].lower(): self.instrument = 'SPEX_SXD'
        if 'spex' in self.instrument.lower():
            if 'observation_date' in list(self.__dict__.keys()): dt = self.observation_date
            elif 'OBS-DATE' in hkys: dt = self.header['OBS-DATE'].replace('-','')
            elif 'OBS_DATE' in hkys: dt = self.header['OBS_DATE'].replace('-','')
            elif 'DATE_OBS' in hkys: dt = self.header['DATE_OBS'].replace('-','')
            elif 'DATE-OBS' in hkys: dt = self.header['DATE-OBS'].replace('-','')
            else: dt = '20000101'
            if int(dt) > 20140800:
                self.instrument = self.instrument.replace('SPEX','USPEX')

# populate defaults
        self.instrument = self.instrument.upper().replace(' ','_')
        if self.instrument in list(INSTRUMENTS.keys()):
            for k in list(INSTRUMENTS[self.instrument].keys()): setattr(self,k,INSTRUMENTS[self.instrument][k])
            self.slitpixelwidth = (self.slitwidth/self.pixelscale).value
            if isNumber(kwargs.get('resolution',False)) != False:
                if kwargs['resolution'] > 0.:
                    self.slitwidth *= self.resolution/kwargs['resolution']
                    self.slitpixelwidth *= self.resolution/kwargs['resolution']
                    self.resolution = kwargs['resolution']
            if isNumber(kwargs.get('slitwidth',False)) != False:
                sl = kwargs['slitwidth']
                if isinstance(sl,u.quantity.Quantity): sl = sl.value
                if s1 > 0.:
                    self.resolution *= (self.slitwidth.value)/sl
                    self.slitpixelwidth *= s1/(self.slitwidth.value)
                    self.slitwidth = kwargs['slitwidth']
                    if not isinstance(self.slitwidth,u.quantity.Quantity): sef.slitwidth *= u.arcsec
            if isNumber(kwargs.get('slitpixelwidth',False)) != False:
                if kwargs['slitpixelwidth'] > 0.:
                    self.resolution *= self.slitpixelwidth/kwargs['slitpixelwidth']
                    self.slitwidth *= kwargs['slitpixelwidth']/self.slitpixelwidth
                    self.slitpixelwidth = kwargs['slitpixelwidth']
        else:
            kys = list(splat.INSTRUMENTS.keys())
            for k in list(INSTRUMENTS[kys[0]].keys()): setattr(self,k,'UNKNOWN')
            self.instrument_name = self.instrument
            self.slitpixelwidth = 3.
            self.resolution = (self.wave[int(0.5*len(self.wave))]/(self.wave[int(0.5*len(self.wave))+2]-self.wave[int(0.5*len(self.wave))-1])).value
            self.slitwidth = 1.0*u.arcsec # not sure I should set this
            print('Warning: {} is not one of the instruments defined for SPLAT'.format(self.instrument))
            print('Setting slit width to {}, slit pixel width to {} and resolution to {}'.format(self.slitwidth,self.slitpixelwidth,self.resolution))

        self.dof = numpy.round(len(self.wave)/self.slitpixelwidth)


# published? assume not
        if not hasattr(self,'published'):
            self.published = kwargs.get('published','N')
        

# information on model
        if self.ismodel == True:
            self.teff = kwargs.get('teff',numpy.nan)
            self.logg = kwargs.get('logg',numpy.nan)
            self.z = kwargs.get('z',numpy.nan)
            self.fsed = kwargs.get('fsed',numpy.nan)
            self.cld = kwargs.get('cld',numpy.nan)
            self.kzz = kwargs.get('kzz',numpy.nan)
            self.slit = kwargs.get('slit',numpy.nan)
            self.modelset = kwargs.get('model','')
# temporary fix of incorrect units in SED spectra            
            if kwargs.get('sed',False):
                self.scale(1.e4)
            try:
                self.name = DEFINED_MODEL_NAMES[self.modelset]
            except:
                self.name = self.modelset
            self.shortname = self.name
            self.name = self.name+' Teff='+str(self.teff)+' logg='+str(self.logg)+' [M/H]='+str(self.z)
            self.fscale = 'Surface'

# populate header            
        else:
            kconv = {'designation': 'DESIG','name': 'NAME','shortname': 'SNAME','ra': 'RA_DEC','dec': 'DEC_DEC','slitwidth': 'SLTW_ARC','source_key': 'SRC_KEY','data_key': 'DATA_KEY','observer': 'OBSERVER', 'data_reference': 'BIB_DATA','program_pi': 'PI','program_number': 'PROGRAM','airmass': 'AIRMASS','reduction_spextool_version': 'VERSION','reduction_person': 'RED_PERS','reduction_date': 'RED_DATE','observation_date': 'OBSDATE','julian_date': 'JDATE','median_snr': 'SNR','resolution': 'RES', 'instrument': 'INSTRUME','wunit': 'XUNITS','funit': 'YUNITS', 'wlabel': 'XTITLE', 'flabel': 'YTITLE', 'opt_type': 'SPT_OPT', 'lit_type': 'SPT_LIT', 'nir_type': 'SPT_NIR', 'spex_type': 'SPT_SPEX', 'gravity_class_nir': 'GRAV_NIR', 'gravity_class_opt': 'GRAV_OPT', 'metallicity_class': 'ZCLASS', 'luminosity_class': 'LUMCLASS', 'color_extremity': 'COLOREX','mu': 'MU','mu_e': 'E_MU', 'mu_ra': 'MU_RA', 'mu_dec': 'MU_DEC', 'parallax': 'PARALLAX', 'parallax_e': 'E_PARALL', 'vtan': 'VTAN','vtan_e': 'E_VTAN','rv': 'RV','rv_e': 'E_RV', 'vsini': 'VSINI', 'vsini_e': 'E_VSINI','distance': 'DISTANCE', 'distance_e': 'E_DISTAN','j_2mass': 'J_2MASS', 'h_2mass': 'H_2MASS', 'ks_2mass': 'K_2MASS', 'j_2mass_e': 'E_J_2MAS', 'h_2mass_e': 'E_H_2MAS', 'ks_2mass_e': 'E_K_2MAS', 'object_type': 'OBJ_TYPE', 'binary': 'BINARY','sbinary': 'SPBINARY', 'companion_name': 'COMPNAME', 'cluster': 'CLUSTER' }
            for k in list(kconv.keys()):
                if kconv[k].upper() not in list(self.header.keys()):
                    try:
                        self.header[kconv[k]] = getattr(self,k)
                    except:
                        self.header[kconv[k]] = ''
            if 'DATE_OBS' not in list(self.header.keys()) and 'observation_date' in list(self.__dict__.keys()):
                self.header['DATE_OBS'] = '{}-{}-{}'.format(self.observation_date[:4],self.observation_date[4:6],self.observation_date[6:])
            if 'TIME_OBS' not in list(self.header.keys()) and 'observation_time' in list(self.__dict__.keys()):
                self.header['TIME_OBS'] = self.observation_time.replace(' ',':')

        self.history.append('Spectrum successfully loaded')
# create a copy to store as the original
        self.original = copy.deepcopy(self)
        return

    def __copy__(self):
        '''
        :Purpose: Make a copy of a Spectrum object
        '''
        s = type(self)()
        s.__dict__.update(self.__dict__)
        return s

# backup version
    def copy(self):
        '''
        :Purpose: Make a copy of a Spectrum object
        '''
        s = type(self)()
        s.__dict__.update(self.__dict__)
        return s

    def __repr__(self):
        '''
        :Purpose: A simple representation of an object is to just give it a name
        '''
        return '{} spectrum of {}'.format(self.instrument,self.name)

    def __add__(self,other):
        '''
        :Purpose: A representation of addition for Spectrum objects which correctly interpolates as a function of wavelength and combines variances

        :Output: a new Spectrum object equal to the spectral sum of the inputs

        :Example:
           >>> import splat
           >>> sp1 = splat.getSpectrum(lucky=True)[0]
           >>> sp2 = splat.getSpectrum(lucky=True)[0]
           >>> sp3 = sp1 + sp2
           >>> sp3
            Spectrum of 2MASS J17373467+5953434 + WISE J174928.57-380401.6
        '''
        sp = copy.deepcopy(self)
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.add(self.flux,f(self.wave)*other.funit)
        sp.variance = sp.variance+n(self.wave)*(other.funit**2)
        sp.noise = sp.variance**0.5
        sp.snr = sp.computeSN()
# update information
        sp.name = self.name+' + '+other.name
        ref = ['date','observer','airmass','designation','source_key','data_key']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,'{} and {}'.format(getattr(self,r),getattr(other,r)))
        sp.history.append('Sum of {} and {}'.format(self.name,other.name))
# reset original
        sp.original = copy.deepcopy(sp)
        return sp

    def __sub__(self,other):
        '''
        :Purpose: A representation of subtraction for Spectrum objects which correctly interpolates as a function of wavelength and combines variances

        :Output: a new Spectrum object equal to the spectral difference of the inputs

        :Example:
           >>> import splat
           >>> sp1 = splat.getSpectrum(lucky=True)[0]
           >>> sp2 = splat.getSpectrum(lucky=True)[0]
           >>> sp3 = sp1 - sp2
           >>> sp3
            Spectrum of 2MASS J17373467+5953434 - WISE J174928.57-380401.6
        '''
        sp = copy.deepcopy(self)
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.subtract(self.flux,f(self.wave)*other.funit)
        sp.variance = sp.variance+n(self.wave)*(other.funit**2)
        sp.noise = sp.variance**0.5
        sp.snr = sp.computeSN()
# update information
        sp.name = self.name+' - '+other.name
        ref = ['date','observer','airmass','designation','source_key','data_key']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,'{} and {}'.format(getattr(self,r),getattr(other,r)))
        sp.history.append('Subtraction of {} by {}'.format(self.name,other.name))
# reset original
        sp.original = copy.deepcopy(sp)
        return sp

    def __mul__(self,other):
        '''
        :Purpose: A representation of multiplication for Spectrum objects which correctly interpolates as a function of wavelength and combines variances

        :Output: a new Spectrum object equal to the spectral product of the inputs

        :Example:
           >>> import splat
           >>> sp1 = splat.getSpectrum(lucky=True)[0]
           >>> sp2 = splat.getSpectrum(lucky=True)[0]
           >>> sp3 = sp1 * sp2
           >>> sp3
            Spectrum of 2MASS J17373467+5953434 x WISE J174928.57-380401.6
        '''
        sp = copy.deepcopy(self)
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.multiply(self.flux,f(self.wave)*other.funit)
        sp.variance = numpy.multiply(sp.flux**2,(\
            numpy.divide(self.variance.value,sp.flux.value**2)+numpy.divide(n(self.wave),f(self.wave)**2)))
        sp.noise = sp.variance**0.5
        sp.snr = sp.computeSN()
# reset originals
#        sp.flux_original=sp.flux
#        sp.noise_original=sp.noise
#        sp.variance_original=sp.variance
#        sp.funit = sp.flux.unit
# update information
        sp.name = self.name+' x '+other.name
        ref = ['date','observer','airmass','designation','source_key','data_key']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,'{} and {}'.format(getattr(self,r),getattr(other,r)))
        sp.history.append('Product of {} by {}'.format(self.name,other.name))
# reset original
        sp.original = copy.deepcopy(sp)
        return sp


    def __div__(self,other):
        '''
        :Purpose: A representation of division for Spectrum objects which correctly interpolates as a function of wavelength and combines variances

        :Output: a new Spectrum object equal to the spectral ratio of the inputs

        :Example:
           >>> import splat
           >>> sp1 = splat.getSpectrum(lucky=True)[0]
           >>> sp2 = splat.getSpectrum(lucky=True)[0]
           >>> sp3 = sp1/sp2
           >>> sp3
            Spectrum of 2MASS J17373467+5953434 + WISE J174928.57-380401.6
        '''
        sp = copy.deepcopy(self)
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.divide(self.flux,f(self.wave)*other.funit)
        sp.variance = numpy.multiply(sp.flux**2,(\
            numpy.divide(self.variance.value,sp.flux.value**2)+numpy.divide(n(self.wave),f(self.wave)**2)))
        sp.noise = sp.variance**0.5
# clean up infinities
        sp.flux = numpy.where(numpy.absolute(sp.flux) == numpy.inf, numpy.nan, sp.flux)*u.erg/u.erg
        sp.noise = numpy.where(numpy.absolute(sp.noise) == numpy.inf, numpy.nan, sp.noise)*u.erg/u.erg
        sp.variance = numpy.where(numpy.absolute(sp.variance) == numpy.inf, numpy.nan, sp.variance)*u.erg/u.erg
# update information
        sp.name = self.name+' / '+other.name
        ref = ['date','observer','airmass','designation','source_key','data_key']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,'{} and {}'.format(getattr(self,r),getattr(other,r)))
        sp.history.append('Division of {} by {}'.format(self.name,other.name))
# reset original
        sp.original = copy.deepcopy(sp)
        return sp

    def __truediv__(self,other):
        '''
        :Purpose: A representation of division for Spectrum objects which correctly interpolates as a function of wavelength and combines variances

        :Output: a new Spectrum object equal to the spectral ratio of the inputs

        :Example:
           >>> import splat
           >>> sp1 = splat.getSpectrum(lucky=True)[0]
           >>> sp2 = splat.getSpectrum(lucky=True)[0]
           >>> sp3 = sp1/sp2
           >>> sp3
            Spectrum of 2MASS J17373467+5953434 + WISE J174928.57-380401.6
        '''
        sp = copy.deepcopy(self)
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.divide(self.flux,f(self.wave)*other.funit)
        sp.variance = numpy.multiply(sp.flux**2,(\
            numpy.divide(self.variance.value,sp.flux.value**2)+numpy.divide(n(self.wave),f(self.wave)**2)))
        sp.noise = sp.variance**0.5
# clean up infinities
        sp.flux = numpy.where(numpy.absolute(sp.flux) == numpy.inf, numpy.nan, sp.flux)*u.erg/u.erg
        sp.noise = numpy.where(numpy.absolute(sp.noise) == numpy.inf, numpy.nan, sp.noise)*u.erg/u.erg
        sp.variance = numpy.where(numpy.absolute(sp.variance) == numpy.inf, numpy.nan, sp.variance)*u.erg/u.erg
# update information
        sp.name = self.name+' / '+other.name
        ref = ['date','observer','airmass','designation','source_key','data_key']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,'{} and {}'.format(getattr(self,r),getattr(other,r)))
        sp.history.append('Division of {} by {}'.format(self.name,other.name))
# reset original
        sp.original = copy.deepcopy(sp)
        return sp

    def computeSN(self):
        '''
        :Purpose: Compute a representative S/N value as the median value of S/N among the top 50% of flux values
        
        :Output: the S/N value

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.computeSN()
           115.96374031163553
        '''
        w = numpy.where(self.flux.value > numpy.median(self.flux.value))
        return numpy.nanmedian(self.flux.value[w]/self.noise.value[w])

    def info(self,**kwargs):
        '''
        :Purpose: Returns a summary of properties for the Spectrum object
        
        :Output: Text summary

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.info()
            Spectrum of NLTT 184
            Observed on 20071012
            at an airmass of 1.145
            Full source designation is J00054517+0723423
            Median S/N = 97.0
            SPLAT source key is 10012.0
            SPLAT spectrum key is 10857
            Data published in Kirkpatrick, J. D. et al. (2010, ApJS, 190, 100-146)
            History:
                Spectrum successfully loaded
        '''
        if (self.ismodel):
            f = '\n{} atmosphere model with the following parmeters:'.format(self.name)
            f+='\nTeff = {} {}'.format(self.teff,SPECTRAL_MODEL_PARAMETERS['teff']['unit'])
            f+='\nlogg = {} {}'.format(self.logg,SPECTRAL_MODEL_PARAMETERS['logg']['unit'])
            f+='\nz = {} {}'.format(self.z,SPECTRAL_MODEL_PARAMETERS['z']['unit'])
            f+='\nfsed = {}'.format(self.fsed)
            f+='\ncld = {}'.format(self.cld)
            f+='\nkzz = {}'.format(self.kzz)
#            f+='\nSmoothed to slit width {} {}'.format(self.slit,SPECTRAL_MODEL_PARAMETERS['slit']['unit'])
            f+='\n\nIf you use this model, please cite {}'.format(spbib.shortRef(SPECTRAL_MODELS[self.modelset]['bibcode']))
            f+='\nbibcode = {}\n'.format(SPECTRAL_MODELS[self.modelset]['bibcode'])
            if kwargs.get('printout',True):
                print(f)
                return
            else:
                return f
        else:
            f = '\n'
            if hasattr(self,'instrument_name'): f+='{} '.format(self.instrument_name)
            if hasattr(self,'name'): f+='spectrum of {}'.format(self.name)
            if hasattr(self,'observer') and hasattr(self,'date'): f+='\nObserved by {} on {}'.format(self.observer,properDate(self.date,output='YYYY MMM DD'))
            if hasattr(self,'airmass'): f+='\nAirmass = {}'.format(self.airmass)
            if hasattr(self,'designation'): f+='\nSource designation = {}'.format(self.designation)
            if hasattr(self,'median_snr'): f+='\nMedian S/N = {}'.format(self.median_snr)
            if hasattr(self,'spex_type'): f+='\nSpeX Classification = {}'.format(self.spex_type)
            if hasattr(self,'lit_type'): f+='\nLiterature Classification = {} from {}'.format(self.lit_type,spbib.shortRef(self.lit_type_ref))
            if hasattr(self,'source_key') and hasattr(self,'data_key'): f+='\nSpectrum key = {}, Source key = {}'.format(int(self.data_key),int(self.source_key))
            if self.published == 'Y':
                f+='\n\nIf you use these data, please cite:'
                f+='\n\t{}'.format(spbib.shortRef(self.data_reference))
                f+='\n\tbibcode: {}'.format(self.data_reference)
            else:
                f+='\n\nUNPUBLISHED DATA'
            f+='\n\nHistory:'
            for h in self.history:
                f+='\n\t{}'.format(h)
            if kwargs.get('printout',True):
                print(f)
                return
            else:
                return f

    def export(self,*args,**kwargs):
        '''
        :Purpose: Exports a Spectrum object to either a fits or ascii file, depending on file extension given.  If no filename is explicitly given, the Spectrum.filename attribute is used. If the filename does not include the full path, the file is saved in the current directory.  Spectrum.export and Spectrum.save_ function in the same manner.

        .. _Spectrum.save : api.html#splat.Spectrum.save

        :param filename: String specifying the filename to save
        :type filename: optional, default = Spectrum.simplefilename

        :Output: An ascii or fits file with the data

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.export('/Users/adam/myspectrum.txt')
           >>> from astropy.io import ascii
           >>> data = ascii.read('/Users/adam/myspectrum.txt',format='tab')
           >>> data
            <Table length=564>
              wavelength          flux          uncertainty   
               float64          float64           float64     
            -------------- ----------------- -----------------
            0.645418405533               0.0               nan
            0.647664904594 6.71920214475e-16 3.71175052033e-16
            0.649897933006 1.26009925777e-15 3.85722895842e-16
            0.652118623257 7.23781818374e-16 3.68178778862e-16
            0.654327988625 1.94569566622e-15 3.21007116982e-16
            ...
        '''

#
# NOTE: EXPORT CURRENTLY DOES NOT RETURN A PROPER HEADER
#
        filename = self.simplefilename
        if len(args) > 0:
            filename = args[0]
        filename = kwargs.get('filename',filename)
        filename = kwargs.get('file',filename)

# determine which type of file
        ftype = filename.split('.')[-1]
#        print(filename,ftype)
# fits file
        if (ftype == 'fit' or ftype == 'fits'):
            try:
                data = numpy.vstack((self.wave.value,self.flux.value,self.noise.value))
                hdu = fits.PrimaryHDU(data)
                for k in list(self.header.keys()):
                    if k != 'HISTORY' and k != 'COMMENT' and k.replace('#','') != '':
                        hdu.header[k] = self.header[k]
                hdu.writeto(filename,clobber=kwargs.get('clobber',True))
            except:
                raise NameError('Problem saving spectrum object to file {}'.format(filename))

# ascii file - by default space delimited (could make this more flexible)
        else:
            try:
                if kwargs.get('model',True) == True:
                    t = Table([self.wave.value,self.flux.value],names=['#wavelength ({})'.format(self.wave.unit),'surface flux ({})'.format(self.flux.unit)])
                else:
                    t = Table([self.wave.value,self.flux.value,self.noise.value],names=['#wavelength','flux','uncertainty'])
                if kwargs.get('header',True):
                    kys = self.header.keys()
                    while 'HISTORY' in kys:
                        kys.remove('HISTORY')
                    hd = ['{} = {}'.format(k,self.header[k]) for k in kys]
                    hd = list(set(hd))
                    hd.sort()
                    t.meta['comments']=hd
#                ascii.write(t,output=filename,format=kwargs.get('format','commented_header'))
                t.write(filename,format='ascii.tab')
            except:
                raise NameError('Problem saving spectrum object to file {}'.format(filename))
        self.history.append('Spectrum saved to {}'.format(filename))
        return


    def save(self,*args,**kwargs):
        '''
        :Purpose: Exports a Spectrum object to either a fits or ascii file, depending on file extension given.  If no filename is explicitly given, the Spectrum.filename attribute is used. If the filename does not include the full path, the file is saved in the current directory.  

        `Spectrum.export()`_ and `Spectrum.save` function in the same manner.

        .. _`Spectrum.export()` : api.html#splat.Spectrum.export
        '''
        self.export(*args,**kwargs)


    def flamToFnu(self):
        return self.toFnu()

    def toFnu(self):
        '''
        :Purpose: Converts flux density F\_nu in units of Jy. This routine changes the underlying Spectrum object. There is no change if the spectrum is already in F\_nu units.
        
        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.toFnu()
           >>> sp.flux.unit
            Unit("Jy")
        '''
        if self.fscale == 'Temperature':
            self.reset()
        self.funit = u.Jy
        self.flabel = 'F_nu'
        self.flux = self.flux.to(self.funit,equivalencies=u.spectral_density(self.wave))
        self.noise = self.noise.to(self.funit,equivalencies=u.spectral_density(self.wave))
        self.snr = self.computeSN()
        self.history.append('Converted to Fnu units of {}'.format(self.funit))
        return

    def fnuToFlam(self):
        return self.toFlam()

    def toFlam(self):
        '''
        :Purpose: Converts flux density to F\_lambda in units of erg/s/cm\^2/Hz. This routine changes the underlying Spectrum object. There is no change if the spectrum is already in F\_lambda units.
        
        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.toFnu()
           >>> sp.flux.unit
            Unit("Jy")
           >>> sp.toFlam()
           >>> sp.flux.unit
            Unit("erg / (cm2 micron s)")
        '''
        if self.fscale == 'Temperature':
            self.reset()
        self.funit = u.erg/(u.cm**2 * u.s * u.micron)
        self.flabel = 'F_lam'
        self.flux = self.flux.to(self.funit,equivalencies=u.spectral_density(self.wave))
        self.noise = self.noise.to(self.funit,equivalencies=u.spectral_density(self.wave))
        self.variance = self.noise**2
        self.snr = self.computeSN()
        self.history.append('Converted to Flam units of {}'.format(self.funit))
        return

    def toAngstrom(self):
        '''
        :Purpose: Converts wavelength to Angstrom. This routine changes the underlying Spectrum object.
        
        :Example:
           >>> import splat
           >>> import astropy.units as u
           >>> sp = splat.getSpectrum(lucky=True)
           >>> sp.wave.unit
            Unit("micron")
           >>> sp.toAngstrom()
           >>> sp.wave.unit
            Unit("Angstrom")
        '''
        self.wunit = u.Angstrom
        self.wave = self.wave.to(self.wunit)
        self.history.append('Converted wavelength to units of {}'.format(self.wunit))
        return

    def toMicron(self):
        '''
        :Purpose: Converts wavelength to microns. This routine changes the underlying Spectrum object.
        
        :Example:
           >>> import splat
           >>> import astropy.units as u
           >>> sp = splat.Spectrum(file='somespectrum',wunit=u.Angstrom)
           >>> sp.wave.unit
            Unit("Angstrom")
           >>> sp.toMicron()
           >>> sp.wave.unit
            Unit("micron")
        '''
        self.wunit = u.micron
        self.wave = self.wave.to(self.wunit)
        self.history.append('Converted wavelength to units of {}'.format(self.wunit))
        return

    def waveUnit(self,wunit):
        '''
        :Purpose: Converts wavelength to microns. This routine changes the underlying Spectrum object.
        
        :Example:
           >>> import splat
           >>> import astropy.units as u
           >>> sp = splat.Spectrum(file='somespectrum',wunit=u.Angstrom)
           >>> sp.wave.unit
            Unit("Angstrom")
           >>> sp.toMicron()
           >>> sp.wave.unit
            Unit("micron")
        '''
        try:
            self.wave = self.wave.to(wunit)
            self.wunit = wunit
            self.history.append('Converted wavelength to units of {}'.format(self.wunit))
        except:
            print('Warning! failed to convert to wavelength unit {}'.format(wunit))            
        return

    def fluxUnit(self,funit):
        '''
        :Purpose: Converts wavelength to microns. This routine changes the underlying Spectrum object.
        
        :Example:
           >>> import splat
           >>> import astropy.units as u
           >>> sp = splat.Spectrum(file='somespectrum',wunit=u.Angstrom)
           >>> sp.wave.unit
            Unit("Angstrom")
           >>> sp.toMicron()
           >>> sp.wave.unit
            Unit("micron")
        '''
        try:
            self.flux = self.flux.to(funit,equivalencies=u.spectral_density(self.wave))
            self.noise = self.noise.to(funit,equivalencies=u.spectral_density(self.wave))
            self.variance = self.noise**2
            self.snr = self.computeSN()
            self.funit = funit
            self.history.append('Converted to flux units of {}'.format(self.funit))
        except:
            print('Warning! failed to convert to flux unit {}'.format(funit))
        return

    def rvShift(self,rv):
        '''
        :Purpose: Shifts the wavelength scale by a radial velocity factor. This routine changes the underlying Spectrum object.
        
        :Example:
           >>> import splat
           >>> import astropy.units as u
           >>> sp = splat.Spectrum(file='somespectrum',wunit=u.Angstrom)
           >>> sp.wave.unit
            Unit("Angstrom")
           >>> sp.toMicron()
           >>> sp.wave.unit
            Unit("micron")
        '''
        if not isinstance(rv,u.quantity.Quantity):
            rv*=(u.km/u.s)
        rv.to(u.km/u.s)
        self.wave = self.wave*(1.+(rv/const.c).to(u.m/u.m))
        return

    def fluxCalibrate(self,filt,mag,**kwargs):
        '''
        :Purpose: Flux calibrates a spectrum given a filter and a magnitude. The filter must be one of those listed in `splat.FILTERS.keys()`. It is possible to specifically set the magnitude to be absolute (by default it is apparent).  This function changes the Spectrum object's flux, noise and variance arrays.
        
        Required Inputs:

        :param filter: string specifiying the name of the filter
        :param mag: number specifying the magnitude to scale to 

        Optional Inputs:

        :param absolute: set to True to specify that the given magnitude is an absolute magnitude, which sets the ``fscale`` keyword in the Spectrum object to 'Absolute' (default = False)
        :param apparent: set to True to specify that the given magnitude is an apparent magnitude, which sets the ``fscale`` flag in the Spectrum object to 'Apparent' (default = False)

        Outputs:

        None, Spectrum object is changed to a flux calibrated spectrum

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.fluxCalibrate('2MASS J',15.0)
           >>> splat.filterMag(sp,'2MASS J')
            (15.002545668628173, 0.017635234089677564)
        '''

        if self.fscale == 'Temperature':
            self.reset()
        if self.funit != u.erg/(u.cm**2 * u.s * u.micron):
            self.toFlam()
        absolute = kwargs.get('absolute',False)
        apparent = kwargs.get('apparent',not absolute)
        apmag,apmag_e = filterMag(self,filt,**kwargs)
# NOTE: NEED TO INCORPORATE UNCERTAINTY INTO SPECTRAL UNCERTAINTY
        if (~numpy.isnan(apmag)):
            self.scale(10.**(0.4*(apmag-mag)),silent=True)
            if (absolute):
                self.fscale = 'Absolute'
                self.history.append('Flux calibrated with {} filter to an absolute magnitude of {}'.format(filt,mag))
            if (apparent):
                self.fscale = 'Apparent'
                self.history.append('Flux calibrated with {} filter to an apparent magnitude of {}'.format(filt,mag))
        self.snr = self.computeSN()

        return

# determine maximum flux, by default in non telluric regions
    def fluxMax(self,**kwargs):
        '''
        :Purpose: Reports the maximum flux of a Spectrum object ignoring nan's.

        :param maskTelluric: masks telluric regions
        :type maskTelluric: optional, default = True

        :Output: maximum flux (with units)

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.normalize()
           >>> sp.fluxMax()
           <Quantity 1.0 erg / (cm2 micron s)>
        '''
        if kwargs.get('maskTelluric',True):            
            try:
                fl = self.flux[numpy.where(numpy.logical_or(\
                    numpy.logical_and(self.wave > 0.9*u.micron,self.wave < 1.35*u.micron),
                    numpy.logical_and(self.wave > 1.42*u.micron,self.wave < 1.8*u.micron),
                    numpy.logical_and(self.wave > 1.92*u.micron,self.wave < 2.3*u.micron)))]
                if isinstance(fl[0],u.quantity.Quantity):
                    fl = [f.value for f in fl]
                return numpy.nanmax(fl)*self.funit
            except:
                pass
        
        fl = self.flux[numpy.where(\
                numpy.logical_and(self.wave > numpy.nanmin(self.wave)+0.1*(numpy.nanmax(self.wave)-numpy.nanmin(self.wave)),self.wave < numpy.nanmax(self.wave)-0.1*(numpy.nanmax(self.wave)-numpy.nanmin(self.wave))))]
        if isinstance(fl[0],u.quantity.Quantity):
            fl = [f.value for f in fl]
        return numpy.nanmax(fl)*self.funit


    def normalize(self,*args,**kwargs):
        '''
        :Purpose: Normalize a spectrum to a maximum value of 1 (in its current units)

        :param waveRange: choose the wavelength range to normalize; can be a list specifying minimum and maximum or a single number to normalize around a particular point
        :type waveRange: optional, default = None

        :Output: maximum flux (with units)

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.normalize()
           >>> sp.fluxMax()
           <Quantity 1.0 erg / (cm2 micron s)>
           >>> sp.normalize(waverange=[2.25,2.3])
           >>> sp.fluxMax()
           <Quantity 1.591310977935791 erg / (cm2 micron s)>
        '''
        rng = kwargs.get('waverange',False)
        rng = kwargs.get('waveRange',rng)
        rng = kwargs.get('range',rng)
        if len(args) > 0:
            rng = args[0]
        if rng != False:
            if not isinstance(rng,list):
                rng = [rng]
            if len(rng) < 2:
                rng = [rng[0]-0.02,rng[0]+0.02]
            self.scale(1./numpy.nanmax(self.flux.value[numpy.where(numpy.logical_and(self.wave > rng[0]*u.micron,self.wave < rng[1]*u.micron))]))
        else:
            self.scale(1./self.fluxMax(**kwargs).value,silent=True)
        self.fscale = 'Normalized'
        if not kwargs.get('silent',False):
            self.history.append('Spectrum normalized')
        self.snr = self.computeSN()
        return

    def plot(self,**kwargs):
        '''
        :Purpose: calls the plotSpectrum_ function, by default showing the noise spectrum and zeropoints. See the plotSpectrum_ API listing for details.

        .. _plotSpectrum: api.html#splat_plot.plotSpectrum

        :Output: A plot of the Spectrum object

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.plot()
        '''
        kwargs['legend'] = kwargs.get('legend',self.name)
        kwargs['showNoise'] = kwargs.get('showNoise',True)
        kwargs['showZero'] = kwargs.get('showZero',True)
        from .plot import plotSpectrum
        f = plotSpectrum(self,**kwargs)
        return f


    def reset(self):
        '''
        :Purpose: Restores a Spectrum to its original read-in state, removing scaling and smoothing. This routine changes the Spectrum object directly and there is no output.
        
        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.fluxMax()
           <Quantity 4.561630292384622e-15 erg / (cm2 micron s)>
           >>> sp.normalize()
           >>> sp.fluxMax()
           <Quantity 0.9999999403953552 erg / (cm2 micron s)>
           >>> sp.reset()
           >>> sp.fluxMax()
           <Quantity 4.561630292384622e-15 erg / (cm2 micron s)>
        '''
        for k in self.original.__dict__.keys():
            if k != 'history':
                setattr(self,k,getattr(self.original,k))

#        self = self.original
#        self.wave = copy.deepcopy(self.original.wave)
#        self.flux = copy.deepcopy(self.original.flux)
#        self.noise = copy.deepcopy(self.original.noise)
#        self.variance = copy.deepcopy(self.original.variance)
#        self.resolution = copy.deepcopy(self.original.resolution)
#        self.slitpixelwidth = copy.deepcopy(self.original.slitpixelwidth)
#        self.slitwidth = self.slitpixelwidth*spex_pixel_scale
#        self.snr = self.computeSN()
#        self.fscale = copy.deepcopy(self.original.fscale)

        self.history.append('Returned to original state')
        self.original = copy.deepcopy(self)
        return



    def scale(self,factor,**kwargs):
        '''
        :Purpose: Scales a Spectrum object's flux and noise values by a constant factor. This routine changes the Spectrum object directly.

        :param factor: A floating point number used to scale the Spectrum object
        :type factor: required, default = None

        :Output: maximum flux (with units)

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.fluxMax()
           <Quantity 1.0577336634332284e-14 erg / (cm2 micron s)>
           >>> sp.computeSN()
           124.5198
           >>> sp.scale(1.e15)
           >>> sp.fluxMax()
           <Quantity 1.0577336549758911 erg / (cm2 micron s)>
           >>> sp.computeSN()
           124.51981
        '''
        self.flux = self.flux*factor
        self.noise = self.noise*factor
        self.variance = self.noise**2
        self.snr = self.computeSN()
        self.fscale = 'Scaled'
        if not kwargs.get('silent',False):
            self.history.append('Spectrum scaled by a factor of {}'.format(factor))
        return

    def showHistory(self):
        '''
        :Purpose: Report history of actions taken on a Spectrum object. This can also be retrieved by printing the attribute Spectrum.history

        :Output: List of actions taken on spectrum
        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.normalize()
           >>> sp.fluxCalibrate('2MASS J',15.0)
           >>> sp.showHistory()
            Spectrum successfully loaded
            Spectrum normalized
            Flux calibrated with 2MASS J filter to an apparent magnitude of 15.0
        '''
        for h in self.history:
            print(h)
        return

    def smooth(self,**kwargs):
        '''
        :Purpose: Smoothes a spectrum either by selecting a constant slit width (smooth in spectral dispersion space), pixel width (smooth in pixel space) or resolution (smooth in velocity space). One of these options must be selected for any smoothing to happen. Changes spectrum directly.

        :param method: the type of smoothing window to use. See http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details.
        :type method: optional, default = Hanning
        :param resolution: Constant resolution to smooth toe(see smoothToResolution_)
        :type resolution: optional, default = None
        :param slitPixelWidth: Number of pixels to smooth in pixel space (see smoothToSlitPixelWidth_)
        :type slitPixelWidth: optional, default = None
        :param slitWidth: Number of pixels to smooth in angular space (see smoothToPixelWidth_)
        :type slitWidth: optional, default = None

        .. _smoothToResolution : api.html#splat.Spectrum.smoothToResolution
        .. _smoothToPixelWidth : api.html#splat.Spectrum.smoothToPixelWidth
        .. _smoothToSlitPixelWidth : api.html#splat.Spectrum.smoothToSlitPixelWidth

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.smoothfluxMax()
           <Quantity 1.0577336634332284e-14 erg / (cm2 micron s)>
           >>> sp.computeSN()
           124.5198
           >>> sp.scale(1.e15)
           >>> sp.fluxMax()
           <Quantity 1.0577336549758911 erg / (cm2 micron s)>
           >>> sp.computeSN()
           124.51981
        '''
        method = kwargs.get('method','hanning')
        kwargs['method'] = method
        swargs = copy.deepcopy(kwargs)
        smv = kwargs.get('slitPixelWidth',None)
        smv = kwargs.get('slitpixelwidth',smv)
        smv = kwargs.get('pixelwidth',smv)
        smv = kwargs.get('pixel',smv)
        if smv != None:
#            del swargs['slitPixelWidth']
            self._smoothToSlitPixelWidth(smv,**swargs)
            return
        smv = kwargs.get('slitWidth',None)
        smv = kwargs.get('slitwidth',smv)
        smv = kwargs.get('slit',smv)
        if smv != None:
#            del swargs['slitWidth']
            self._smoothToSlitWidth(smv,**swargs)
            return
        smv = kwargs.get('resolution',None)
        smv = kwargs.get('res',smv)
        if smv != None:
#            del swargs['resolution']
            self._smoothToResolution(smv,**swargs)
            return
        return

    def _smoothToResolution(self,res,**kwargs):
        '''
        :Purpose: Smoothes a spectrum to a constant or resolution (smooth in velocity space). Changes spectrum directly.  Note that no smoothing is done if requested resolution is greater than the current resolution

        :param resolution: number giving the desired resolution
        :type resolution: required
        :param method: the type of smoothing window to use. See http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details.
        :type method: optional, default = Hanning
        :param overscale: used for computing number of samples in the window
        :type overscale: optional, default = 10.

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.resolution()
           120
           >>> sp.computeSN()
           21.550974
           >>> sp.smoothToResolution(50)
           >>> sp.resolution()
           50
           >>> sp.computeSN()
           49.459522314460855
        '''

        overscale = kwargs.get('overscale',10.)
        method = kwargs.get('method','hamming')
        kwargs['method'] = method

# do nothing if requested resolution is higher than current resolution
        if (res < self.resolution):
# sample onto a constant resolution grid at 5x current resolution
            r = res*overscale
            waveRng = self.waveRange()
            npix = numpy.floor(numpy.log(waveRng[1]/waveRng[0])/numpy.log(1.+1./r))
            wave_sample = [waveRng[0].value*(1.+1./r)**i for i in numpy.arange(npix)]
            f = interp1d(self.wave,self.flux,bounds_error=False,fill_value=0.)
            v = interp1d(self.wave,self.variance,bounds_error=False,fill_value=numpy.nan)
            flx_sample = f(wave_sample)*self.funit
            var_sample = v(wave_sample)*self.funit**2
# now convolve a function to smooth resampled spectrum
            window = signal.get_window(method,numpy.round(overscale))
            neff = numpy.sum(window)/numpy.nanmax(window)        # effective number of pixels
            flx_smooth = signal.convolve(flx_sample, window/numpy.sum(window), mode='same')
            var_smooth = signal.convolve(var_sample, window/numpy.sum(window), mode='same')/neff
# resample back to original wavelength grid
            f = interp1d(wave_sample,flx_smooth,bounds_error=False,fill_value=0.)
            v = interp1d(wave_sample,var_smooth,bounds_error=False,fill_value=0.)
            self.flux = f(self.wave.value)*self.funit
            self.variance = v(self.wave.value)*self.funit**2
            self.noise = [ns**0.5 for ns in self.variance.value]*self.funit
            self.snr = self.computeSN()
            self.slitpixelwidth *= self.resolution/res
            self.slitwidth *= self.resolution/res
            self.resolution = res
            self.history.append('Smoothed to a constant resolution of {}'.format(self.resolution))
        else:
            print('\nTarget resolution {} less than current resolution {}; no change made'.format(res,self.resolution))
        return

    def _smoothToSlitPixelWidth(self,width,**kwargs):
        '''
        :Purpose: Smoothes a spectrum to a constant slit pixel width (smooth in pixel space). Changes spectrum directly.  Note that no smoothing is done if requested width is greater than the current slit width.

        :param width: number giving the desired smoothing scale in pixels
        :type width: required
        :param method: the type of smoothing window to use. See http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details.
        :type method: optional, default = Hanning

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.slitpixelwidth
           3.33
           >>> sp.resolution
           120
           >>> sp.computeSN()
           105.41789
           >>> sp.smoothToSlitPixelWidth(10)
           >>> sp.slitpixelwidth
           10
           >>> sp.resolution
           39.96
           >>> sp.computeSN()
           235.77536310249229
        '''
        method = kwargs.get('method','hanning')
        kwargs['method'] = method
# do nothing if requested resolution is higher than current resolution
        if (width > self.slitpixelwidth):
# convolve a function to smooth spectrum
            window = signal.get_window(method,numpy.round(width))
            neff = numpy.sum(window)/numpy.nanmax(window)        # effective number of pixels
            self.flux = signal.convolve(self.flux.value, window/numpy.sum(window), mode='same')*self.funit
            self.variance = signal.convolve(self.variance.value, window/numpy.sum(window), mode='same')/neff*(self.funit**2)
            self.noise = [n**0.5 for n in self.variance.value]*self.funit
            self.snr = self.computeSN()
            self.resolution *= self.slitpixelwidth/width
            self.slitwidth *= width/self.slitpixelwidth
            self.slitpixelwidth = width
            self.history.append('Smoothed to slit pixel width of {}'.format(self.slitpixelwidth))
        else:
            print('\nTarget slit width {} is less than current slit width {}; no change made'.format(width,self.slitpixelwidth))
        return

    def _smoothToSlitWidth(self,width,**kwargs):
        '''
        :Purpose: Smoothes a spectrum to a constant slit angular width (smooth in dispersion space). Changes spectrum directly.  Note that no smoothing is done if requested width is greater than the current slit width.

        :param width: number giving the desired smoothing scale in arcseconds
        :type width: required
        :param method: the type of smoothing window to use. See http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details.
        :type method: optional, default = Hanning

        :Output: maximum flux (with units)

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.slitwidth
           0.4995
           >>> sp.resolution
           120
           >>> sp.computeSN()
           105.41789
           >>> sp.smoothToSlitWidth(2.0)
           >>> sp.slitwidth
           2.0
           >>> sp.resolution
           29.97
           >>> sp.computeSN()
           258.87135134070593
        '''
        method = kwargs.get('method','hanning')
        kwargs['method'] = method
        if not isinstance(width,u.quantity.Quantity):
            width*=u.arcsec
        pwidth = self.slitpixelwidth*(width/self.slitwidth).value
        self._smoothToSlitPixelWidth(pwidth,**kwargs)
        return

    def surface(self,radius):
        '''
        :Purpose: Convert to surface fluxes given a radius, assuming at absolute fluxes
        .. note:: UNTESTED, NEED TO ADD IN UNCERTAINTIES
        '''
        if self.fscale == 'Surface':
            return
        if self.fscale != 'Absolute':
            print('To convert to surface fluxes you must first scale spectrum to absolute (10 pc) flux units')
            return
        r = copy.deepcopy(radius)
        if ~isinstance(r,u.quantity.Quantity):
            r*=const.R_sun
        self.scale(((10.*u.pc/r).to(u.m/u.m).value)**2,silent=True)
        self.history.append('Converted to surface fluxes assuming a radius of {} solar radii'.format((r/const.R_sun).to(u.m/u.m)))
        self.fscale = 'Surface'
        return

    def brightnessTemperature(self):
        '''
        :Purpose: Convert to surface fluxes given a radius, assuming at absolute fluxes
        .. note: UNTESTED
        '''
        if self.fscale == 'Temperature':
            return
        if self.fscale != 'Surface':
            print('To convert to brightness temperature you must first scale spectrum to surface flux units')
            return
        fs = copy.deepcopy(self.flux).to(u.erg/u.s/u.cm**3)
        fse = copy.deepcopy(self.noise).to(u.erg/u.s/u.cm**3)
        w = copy.deepcopy(self.wave).to(u.cm)
        x = (2.*numpy.pi*const.h.to(u.erg*u.s)*(const.c.to(u.cm/u.s)**2)/(fs*(w**5))).to(u.m/u.m).value
        self.flux = (const.h*const.c/(const.k_B*w)).to(u.K)/numpy.log(1.+x)
        self.noise = self.flux*x/((1.+x)*numpy.log(1.+x))*((fs/fse).to(u.m/u.m))
        self.variance = self.noise**2
        self.snr = self.computeSN()
        self.history.append('Converted to brightness temperature')
        self.fscale = 'Temperature'
        self.funit = u.K
        return

    def temperature(self):
        self.brightnessTemperature()


    def trim(self,rng,**kwargs):
        '''
        :Purpose: Trims a spectrum to be within a certain wavelength range or set of ranges. Data outside of these ranges are excised from the wave, flux and noise arrays. The full spectrum can be restored with the reset() procedure.

        :param range: the range(s) over which the spectrum is retained - a series of nested 2-element arrays

        .. _smoothToResolution : api.html#splat.Spectrum.smoothToResolution
        .. _smoothToPixelWidth : api.html#splat.Spectrum.smoothToPixelWidth
        .. _smoothToSlitPixelWidth : api.html#splat.Spectrum.smoothToSlitPixelWidth

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.smoothfluxMax()
           <Quantity 1.0577336634332284e-14 erg / (cm2 micron s)>
           >>> sp.computeSN()
           124.5198
           >>> sp.scale(1.e15)
           >>> sp.fluxMax()
           <Quantity 1.0577336549758911 erg / (cm2 micron s)>
           >>> sp.computeSN()
           124.51981
        '''

        if isinstance(rng,float):
            rng = [rng-0.1,rng+0.1]
        if ~isinstance(rng[0],list):
            rng = [rng]
        mask = numpy.zeros(len(self.wave))
        for r in rng:
            if ~isinstance(r[0],u.quantity.Quantity):
                r*=u.micron
            mask[numpy.where(((self.wave.value >= r[0].value) & (self.wave.value <= r[1].value)))] = 1
        w = numpy.where(mask == 1)
        self.wave = self.wave[w]
        self.flux = self.flux[w]
        self.noise = self.noise[w]
        self.variance = self.variance[w]
        self.flam = self.flux
        self.nu = self.wave.to('Hz',equivalencies=u.spectral())
        self.fnu = self.flux.to('Jy',equivalencies=u.spectral_density(self.wave))
        self.noisenu = self.noise.to('Jy',equivalencies=u.spectral_density(self.wave))
        self.snr = self.computeSN()
        self.history.append('Spectrum trimmed to ranges {}'.format(rng))
        return


    def waveRange(self):
        '''
        :Purpose: Return the wavelength range of the current Spectrum object.

        :Output: 2-element array giving minimum and maximum of wavelength range

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.slitwidth
           [<Quantity 0.6447611451148987 micron>, <Quantity 2.5517737865448 micron>]
        '''

        ii = numpy.where(self.flux.value > 0)
        return [numpy.nanmin(self.wave[ii]), numpy.nanmax(self.wave[ii])]


# stitch spectrum
def stitch(s1,s2,rng = [0.8,0.9],**kwargs):
    '''
    :Purpose: Stitches together two spectra covering different wavelength scales.

    :Output: New spectrum object of stitched spectrum

    :Example:
       >>> import splat
       >>> spopt = splat.Spectrum(file='myopticalspectrum.fits')
       >>> spnir = splat.Spectrum(file='myopticalspectrum.fits')
       >>> sp = splat.stitch(spopt,spnir,rng=[0.8,0.9],trim=[0.35,2.4])
    '''

# generate copies of spectrum objects
    sp1 = copy.deepcopy(s1)
    sp2 = copy.deepcopy(s2)
    
# assert common units
    wunit = kwargs.get('wunit',SPECTRAL_MODEL_WAVE_UNIT)
    sp1.waveUnit(wunit)
    sp2.waveUnit(wunit)
    funit = kwargs.get('funit',SPECTRAL_MODEL_FLUX_UNIT)
    sp1.fluxUnit(funit)
    sp2.fluxUnit(funit)

# scale to overlap region - this assumes full overlap over the stitch region
# for now doing simple normalization - better to do uncertainty-weighted scaling
    if numpy.max(sp1.wave.value) >= rng[1] and numpy.min(sp2.wave.value) <= rng[0]:
        sp1.normalize(rng)  
        sp2.normalize(rng)
    else:
        raise ValueError('Stich region {} to {} does not overlap both spectra'.format(rng[0],rng[1]))

# interpolation of second spectrum
    f2r = interp1d(sp2.wave.value,sp2.flux.value)
    v2r = interp1d(sp2.wave.value,sp2.variance.value)

# start with first spectrum
    wave = sp1.wave.value[numpy.where(sp1.wave.value < rng[0])]
    flux = sp1.flux.value[numpy.where(sp1.wave.value < rng[0])]
    variance = sp1.variance.value[numpy.where(sp1.wave.value < rng[0])]

# mixed region
    w1 = sp1.wave.value[numpy.where(numpy.logical_and(sp1.wave.value >= rng[0],sp1.wave.value < rng[1]))]
    f1 = sp1.flux.value[numpy.where(numpy.logical_and(sp1.wave.value >= rng[0],sp1.wave.value < rng[1]))]
    v1 = sp1.variance.value[numpy.where(numpy.logical_and(sp1.wave.value >= rng[0],sp1.wave.value < rng[1]))]
    f12 = numpy.array([((f1[i]/v1[i])+(f2r(w1[i])/v2r(w1[i])))/(1./v1[i]+1./v2r(w1[i])) for i in range(len(w1))])
    v12 = numpy.array([(1./(1./v1[i]+1./v2r(w1[i]))) for i in range(len(w1))])
    wave = numpy.append(wave,w1)
    flux = numpy.append(flux,f12)
    variance = numpy.append(variance,v12)

# tail of second spectrum
    wave = numpy.append(wave,sp2.wave.value[numpy.where(sp2.wave.value >= rng[1])])
    flux = numpy.append(flux,sp2.flux.value[numpy.where(sp2.wave.value >= rng[1])])
    variance = numpy.append(variance,sp2.variance.value[numpy.where(sp2.wave.value >= rng[1])])

# put back units
    wave = wave*wunit
    flux = flux*funit
    variance = variance*(funit**2)

# create new spectrum object containing combined spectrum
    sp = copy.deepcopy(sp1)
    sp.wave = wave
    sp.flux = flux
    sp.variance = variance
    sp.noise = sp.variance**0.5
    sp.snr = sp.computeSN()
    sp.name = 'Stitched spectrum of {} and {}'.format(sopt.name,snir.name)

# trim if desired
    if kwargs.get('trim',False) != False:
        sp.trim(kwargs['trim'])
        
    return sp


#####################################################
#################   DATA ACCESS   ###################
#####################################################


def getSpectrum(*args, **kwargs):
    '''
    :Purpose: Gets a spectrum from the SPLAT library using various selection criteria. Calls searchLibrary_ to select spectra; if any found it routines an array of Spectrum objects, otherwise an empty array. 

    .. _searchLibrary : api.html#splat_db.searchLibrary

    :Output: An array of Spectrum objects that satisfy the search criteria

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
        Retrieving 1 file
    >>> sparr = splat.getSpectrum(spt='M7')
        Retrieving 120 files
    >>> sparr = splat.getSpectrum(spt='T5',young=True)
        No files match search criteria
    '''

    if kwargs.get('lucky',False) == True:
        kwargs['published'] = True
    result = []
    kwargs['output'] = 'all'
    search = searchLibrary(*args, **kwargs)

    if len(search) > 0:
        files = []
        for i,x in enumerate(search['DATA_KEY']):
            files.append(str(search['DATA_KEY'][i])+'_'+str(search['SOURCE_KEY'][i])+'.fits')

# return just the filenames
        if (kwargs.get('list',False) != False):
            return files

        if (len(files) == 1):
            print('\nRetrieving 1 file\n')
            result.append(Spectrum(files[0],header=search[0]))
        else:
            if (kwargs.get('lucky',False) == True):
                print('\nRetrieving 1 lucky file\n')
                ind = numpy.random.choice(numpy.arange(len(files)))
                result.append(Spectrum(files[ind],header=search[ind]))
            else:
                print('\nRetrieving {} files\n'.format(len(files)))
                for i,x in enumerate(files):
                    result.append(Spectrum(x,header=search[i:i+1]))

    else:
        if checkAccess() == False:
            sys.stderr.write('\nNo published files match search criteria\n\n')
        else:
            sys.stderr.write('\nNo files match search criteria\n\n')

    return result



def getStandard(spt, **kwargs):
    '''
    :Purpose: Gets one of the pre-defined spectral standards from the SPLAT library.

    :param spt: Spectral type of standard desired, either string ('M7') or numberic (17)
    :type spt: required
    :param sd: Set to True to get a subdwarf standard
    :type sd: optional, default = False
    :param esd: Set to True to get an extreme subdwarf standard
    :type esd: optional, default = False

    :Example:
    >>> import splat
    >>> sp = splat.getStandard('M7')[0]
        Spectrum of VB 8
    >>> sparr = splat.getStandard('T5',esd=True)
        Type esdT5.0 is not in esd standards: try one of the following:
        ['esdM5.0', 'esdM7.0', 'esdM8.5']
    '''


# make sure standards are read in
#    initiateStandards(**kwargs)

# set up subtype to use, convert to number then back to string
    if isinstance(spt,str):
        sptstr = copy.deepcopy(spt)
        spt = typeToNum(spt)
    else: sptstr = ''

# get standards
    if kwargs.get('esd',False) or 'esd' in sptstr:
        stds = STDS_ESD_SPEX
        kys = STDS_ESD_SPEX_KEYS
        subclass = 'esd'
    elif kwargs.get('sd',False) or 'sd' in sptstr:
        stds = STDS_SD_SPEX
        kys = STDS_SD_SPEX_KEYS
        subclass = 'sd'
    else:
        stds = STDS_DWARF_SPEX
        kys = STDS_DWARF_SPEX_KEYS
        subclass = ''

    spt = typeToNum(spt,subclass=subclass)

# not a valid subtype
    if spt not in list(kys.keys()):
        print('Type {} is not in {} standards: try one of the following:'.format(spt,subclass))
        print(sorted(list(kys.keys())))
        return Spectrum()

# not yet read in
    if spt not in list(stds.keys()):
        stds[spt] = Spectrum(kys[spt])
    
    return stds[spt]



def initiateStandards(**kwargs):
    '''
    :Purpose: Initiates the spectral standards in the SpeX library. By default this loads the dwarfs standards, but you can also specify loading of subdwarf and extreme subdwarf standards as well. Once loaded, these standards remain in memory.

    :param sd: Set equal to True to load subdwarf standards
    :type sd: optional, default = False
    :param esd: Set equal to True to load extreme subdwarf standards
    :type esd: optional, default = False

    :Example:
    >>> import splat
    >>> splat.initiateStandards()
    >>> splat.SPEX_STDS['M5.0']
    Spectrum of Gl51
    '''

# choose what kind of standards desired - d, sd, esd
# and read in standards into dictionary if they haven't been already
    if kwargs.get('all',False):
        swargs = copy.deepcopy(kwargs)
        del swargs['all']
        initiateStandards()
        initiateStandards(sd=True)
        initiateStandards(esd=True)
        return

    elif kwargs.get('sd',False):
        stds = STDS_SD_SPEX
        kys = STDS_SD_SPEX_KEYS
    elif kwargs.get('esd',False):
        stds = STDS_ESD_SPEX
        kys = STDS_ESD_SPEX_KEYS
    else:
        stds = STDS_DWARF_SPEX
        kys = STDS_DWARF_SPEX_KEYS
    for t in list(kys.keys()):
        if t not in list(stds.keys()):
            stds[t] = Spectrum(kys[t])
            stds[t].normalize()
            stds[t].name += ' ({})'.format(t)

    return



def keySource(keys, **kwargs):
    '''
    :Purpose: Takes a source key and returns a table with the source information
    :param keys: source key or a list of source keys
    :Example:
    >>> import splat
    >>> print spl.keySource(10001)
        SOURCE_KEY           NAME              DESIGNATION    ... NOTE SELECT
        ---------- ------------------------ ----------------- ... ---- ------
             10001 SDSS J000013.54+255418.6 J00001354+2554180 ...        True
    >>> print spl.keySource([10105, 10623])
        SOURCE_KEY          NAME             DESIGNATION    ... NOTE SELECT
        ---------- ---------------------- ----------------- ... ---- ------
             10105 2MASSI J0103320+193536 J01033203+1935361 ...        True
             10623 SDSS J09002368+2539343 J09002368+2539343 ...        True
    >>> print spl.keySource(1000001)
        No sources found with source key 1000001
        False
    '''

# vectorize
    if isinstance(keys,list) == False:
        keys = [keys]

#    sdb = ascii.read(SPLAT_PATH+DB_FOLDER+SOURCES_DB, delimiter='\t',fill_values='-99.',format='tab')
#    sdb = fetchDatabase(SPLAT_PATH+DB_FOLDER+SOURCES_DB)
    sdb = copy.deepcopy(DB_SOURCES)
    sdb['SELECT'] = [x in keys for x in sdb['SOURCE_KEY']]

    if sum(sdb['SELECT']) == 0.:
        print('No sources found with source key {}'.format(keys[0]))
        return False
    else:
        db = sdb[:][numpy.where(sdb['SELECT']==1)]
        return db


def keySpectrum(keys, **kwargs):
    '''
    :Purpose: Takes a spectrum key and returns a table with the spectrum and source information
    :param keys: spectrum key or a list of source keys
    :Example:
    >>> import splat
    >>> print spl.keySpectrum(10001)
        DATA_KEY SOURCE_KEY    DATA_FILE     ... COMPANION COMPANION_NAME NOTE_2
        -------- ---------- ---------------- ... --------- -------------- ------
           10001      10443 10001_10443.fits ...
    >>> print spl.keySpectrum([10123, 11298])
        DATA_KEY SOURCE_KEY    DATA_FILE     ... COMPANION COMPANION_NAME NOTE_2
        -------- ---------- ---------------- ... --------- -------------- ------
           11298      10118 11298_10118.fits ...
           10123      10145 10123_10145.fits ...
    >>> print spl.keySpectrum(1000001)
        No spectra found with spectrum key 1000001
        False
    '''

    verbose = kwargs.get('verbose',False)

# vectorize
    if isinstance(keys,list) == False:
        keys = [keys]

#    sdb = ascii.read(SPLAT_PATH+DB_FOLDER+SPECTRA_DB, delimiter='\t',fill_values='-99.',format='tab')
#    sdb = fetchDatabase(SPLAT_PATH+DB_FOLDER+SPECTRA_DB)
    sdb = copy.deepcopy(DB_SPECTRA)
    sdb['SELECT'] = [x in keys for x in sdb['DATA_KEY']]

    if sum(sdb['SELECT']) == 0.:
        if verbose: print('No spectra found with spectrum key {}'.format(keys[0]))
        return False
    else:
#        s2db = ascii.read(SPLAT_PATH+DB_FOLDER+SOURCES_DB, delimiter='\t',fill_values='-99.',format='tab')
#        s2db = fetchDatabase(SPLAT_PATH+DB_FOLDER+SOURCES_DB)
        s2db = copy.deepcopy(DB_SOURCES)
        db = join(sdb[:][numpy.where(sdb['SELECT']==1)],s2db,keys='SOURCE_KEY')
        return db



def searchLibrary(*args, **kwargs):
    '''
    :Purpose: Search the SpeX database to extract the key reference for that Spectrum

    :param optional name: search by source name (e.g., ``name = 'Gliese 570D'``)
    :param optional shortname: search be short name (e.g. ``shortname = 'J1457-2124'``)
    :param optional designation: search by full designation (e.g., ``designation = 'J11040127+1959217'``)
    :param optional coordinate: search around a coordinate by a radius specified by radius keyword (e.g., ``coordinate = [180.,+30.], radius = 10.``)
    :param radius: search radius in arcseconds for coordinate search
    :type radius: optional, default = 10
    :param optional spt: search by SpeX spectral type; single value is exact, two-element array gives range (e.g., ``spt = 'M7'`` or ``spt = [24,39]``)
    :param optional spex_spt: same as ``spt``
    :param optional opt_spt: same as ``spt`` for literature optical spectral types
    :param optional nir_spt: same as ``spt`` for literature NIR spectral types
    :param optional jmag, hmag, kmag: select based on faint limit or range of J, H or Ks magnitudes (e.g., ``jmag = [12,15]``)
    :param optional snr: search on minimum or range of S/N ratios (e.g., ``snr = 30.`` or ``snr = [50.,100.]``)
    :param optional subdwarf, young, binary, spbinary, red, blue, giant, wd, standard: classes to search on (e.g., ``young = True``)
    :param logic: search logic, can be ``and`` or ``or``
    :type logic: optional, default = 'and'
    :param combine: same as logic
    :type combine: optional, default = 'and'
    :param optional date: search by date (e.g., ``date = '20040322'``) or range of dates (e.g., ``date=[20040301,20040330]``)
    :param optional reference: search by list of references (bibcodes) (e.g., ``reference = '2011ApJS..197...19K'``)
    :param sort: by default returned table is sorted by designation; set this parameter to a column name to sort on a different parameter
    :type sort: optional, default = 'DESIGNATION'
    :param reverse: set to True to do a reverse sort
    :type reverse: optional, default = False
    :param list: if True, return just a list of the data files
    :type list: optional, default = False
    :param lucky: if True, return one randomly selected spectrum from the selected sample
    :type lucky: optional, default = False


    :param output: returns desired output of selected results
    :type output: optional, default = 'all'
    :param logic: search logic, can be and`` or ``or``
    :type logic: optional, default = 'and'
    :param combine: same as logic
    :type combine: optional, default = 'and'
    :Example:
    >>> import splat
    >>> print SearchLibrary(shortname = '2213-2136')
        DATA_KEY SOURCE_KEY    DATA_FILE     ... SHORTNAME  SELECT_2
        -------- ---------- ---------------- ... ---------- --------
           11590      11586 11590_11586.fits ... J2213-2136      1.0
           11127      11586 11127_11586.fits ... J2213-2136      1.0
           10697      11586 10697_11586.fits ... J2213-2136      1.0
           10489      11586 10489_11586.fits ... J2213-2136      1.0
    >>> print SearchLibrary(shortname = '2213-2136', output = 'OBSERVATION_DATE')
        OBSERVATION_DATE
        ----------------
                20110908
                20080829
                20060902
                20051017

    .. note:: Note that this is currently only and AND search - need to figure out how to a full SQL style search
    '''

# program parameters
    ref = kwargs.get('output','all')
    radius = kwargs.get('radius',10.)      # search radius in arcseconds
    classes = ['YOUNG','SUBDWARF','BINARY','SPBINARY','RED','BLUE','GIANT','WD','STANDARD','COMPANION']
    verbose = kwargs.get('verbose',False)

# logic of search
    logic = 'and'         # default combination
    logic = kwargs.get('combine',logic).lower()
    logic = kwargs.get('logic',logic).lower()
    if (logic != 'and' and logic != 'or'):
        raise ValueError('\nLogical operator '+logic+' not supported\n\n')

# read in source database and add in shortnames and skycoords
#    source_db = ascii.read(SPLAT_PATH+DB_FOLDER+SOURCES_DB, delimiter='\t', fill_values='-99.', format='tab')
#    source_db = fetchDatabase(SOURCES_DB)
    source_db = copy.deepcopy(DB_SOURCES)

# first search by source parameters
    source_db['SELECT'] = numpy.zeros(len(source_db['RA']))
    count = 0.

# search by source key
    idkey = kwargs.get('sourcekey',False)
    idkey = kwargs.get('source_key',idkey)
    idkey = kwargs.get('source',idkey)
    idkey = kwargs.get('idkey',idkey)
    idkey = kwargs.get('id',idkey)
    if idkey != False:
        if not isinstance(idkey,list):
            idkey = [idkey]
        if isinstance(idkey[0],str):
            idkey = [int(i) for i in idkey]
        for s in idkey:
            source_db['SELECT'][numpy.where(source_db['SOURCE_KEY'] == s)] += 1
        count+=1.
# search by name
    if kwargs.get('name',False) != False:
        nm = kwargs['name']
        if isinstance(nm,str):
            nm = [nm]
        for n in nm:
            source_db['SELECT'][numpy.where(source_db['NAME'] == n)] += 1
        count+=1.
# search by shortname
    if kwargs.get('shortname',False) != False:
        if 'SHORTNAME' not in source_db.keys():
            source_db['SHORTNAME'] = [designationToShortName(x) for x in source_db['DESIGNATION']]
        sname = kwargs['shortname']
        if isinstance(sname,str):
            sname = [sname]
        for sn in sname:
            if sn[0].lower() != 'j':
                sn = 'J'+sn
            source_db['SELECT'][numpy.where(source_db['SHORTNAME'] == sn)] += 1
        count+=1.
# exclude by shortname
    sname = kwargs.get('excludesource',False)
    sname = kwargs.get('excludeshortname',sname)
    if sname != False and len(sname) > 0:
        if isinstance(sname,str):
            sname = [sname]
        for sn in sname:
            if sn[0].lower() != 'j':
                sn = 'J'+sn
#            t = numpy.sum(source_db['SELECT'][numpy.where(source_db['SHORTNAME'] != sn)])
            source_db['SELECT'][numpy.where(source_db['SHORTNAME'] != sn)] += 1
#            if numpy.sum(source_db['SELECT'][numpy.where(source_db['SHORTNAME'] != sn)]) > t:
#                print('rejected '+sn)
        count+=1.
# search by reference list
    if kwargs.get('reference',False) != False:
        refer = kwargs['reference']
        if isinstance(refer,str):
            refer = [refer]
        for r in refer:
            source_db['SELECT'][numpy.where(source_db['DISCOVERY_REFERENCE'] == r)] += 1
        count+=1.
# search by designation
    if kwargs.get('designation',False) != False:
        desig = kwargs['designation']
        if isinstance(desig,str):
            desig = [desig]
        for d in desig:
            source_db['SELECT'][numpy.where(source_db['DESIGNATION'] == d)] += 1
        count+=1.
# search by coordinate - NOTE: THIS IS VERY SLOW RIGHT NOW
    if kwargs.get('coordinate',False) != False:
        print('\nWarning, search by coordinates may take a few minutes\n')
        coord = kwargs['coordinate']
        if isinstance(coord,SkyCoord):
            cc = coord
        else:
            cc = properCoordinates(coord)
# calculate skycoords
        if 'SKYCOORD' not in source_db.keys():
            s = []
            for i in numpy.arange(len(source_db['RA'])):
                try:        # to deal with a blank string
                    s.append(SkyCoord(ra=float(source_db['RA'][i])*u.degree,dec=float(source_db['DEC'][i])*u.degree,frame='icrs'))
                except:
                    s.append(SkyCoord(ra=numpy.nan*u.degree,dec=numpy.nan*u.degree,frame='icrs'))
#                if numpy.mod(i,len(source_db['RA'])/10.) < 1 and i != 0:
#                    print('\b{:.0f}%...'.format(100*i/len(source_db['RA'])))
            source_db['SKYCOORD'] = s
#        print('measuring separations')
#        source_db['SEPARATION'] = [cc.separation(source_db['SKYCOORDS'][i]).arcsecond for i in numpy.arange(len(source_db['SKYCOORDS']))]
        source_db['SEPARATION'] = [cc.separation(c).arcsecond for c in source_db['SKYCOORD']]
#        print('done')
        source_db['SELECT'][numpy.where(source_db['SEPARATION'] <= radius)] += 1
        count+=1.
#        print(count,numpy.max(source_db['SELECT']))

# search by spectral type
# THIS COULD USE SOME CLEAN UP
    

    spt_range = kwargs.get('spt_range',False)
    spt_range = kwargs.get('spt',spt_range)
    spt_type = kwargs.get('spt_type','LIT_TYPE')
    if kwargs.get('opt_spt',False) != False:
        spt_type = 'OPT_TYPE'
        spt_range = kwargs['opt_spt']
    if kwargs.get('optspt',False) != False:
        spt_type = 'OPT_TYPE'
        spt_range = kwargs['optspt']
    if kwargs.get('spex_spt',False) != False:
        spt_type = 'SPEX_TYPE'
        spt_range = kwargs['spex_spt']
    if kwargs.get('spexspt',False) != False:
        spt_type = 'SPEX_TYPE'
        spt_range = kwargs['spexspt']
    if kwargs.get('nir_spt',False) != False:
        spt_type = 'NIR_TYPE'
        spt_range = kwargs['nir_spt']
    if kwargs.get('nirspt',False) != False:
        spt_type = 'NIR_TYPE'
        spt_range = kwargs['nirspt']
    if kwargs.get('lit_spt',False) != False:
        spt_type = 'LIT_TYPE'
        spt_range = kwargs['lit_spt']
    if kwargs.get('litspt',False) != False:
        spt_type = 'LIT_TYPE'
        spt_range = kwargs['litspt']

    if spt_type.lower() == 'lit_type' or spt_type.lower() == 'lit' or spt_type.lower() == 'literature' or spt_type.lower() == 'pub' or spt_type.lower() == 'published':
        spt_type = 'LIT_TYPE'
    elif spt_type.lower() == 'spex_type' or spt_type.lower() == 'spex':
        spt_type = 'SPEX_TYPE'
    elif spt_type.lower() == 'opt_type' or spt_type.lower() == 'optical_type' or spt_type.lower() == 'optical' or spt_type.lower() == 'opt':
        spt_type = 'OPT_TYPE'
    elif spt_type.lower() == 'nir_type' or spt_type.lower() == 'nir' or spt_type.lower() == 'infrared' or spt_type.lower() == 'near-infrared':
        spt_type = 'NIR_TYPE'
    else:
        spt_type = 'LIT_TYPE'
    if spt_range != False and spt_type != 'SPEX_TYPE':
        if not isinstance(spt_range,list):        # one value = only this type
            spt_range = [spt_range,spt_range]
        if isinstance(spt_range[0],str):          # convert to numerical spt
            spt_range = [typeToNum(spt_range[0]),typeToNum(spt_range[1])]
        source_db['SPTN'] = [typeToNum(x) for x in source_db[spt_type]]
        source_db['SELECT'][numpy.where(numpy.logical_and(source_db['SPTN'] >= spt_range[0],source_db['SPTN'] <= spt_range[1]))] += 1
        count+=1.

# search by magnitude range
    if kwargs.get('jmag',False) != False:
        mag = kwargs['jmag']
        if not isinstance(mag,list):        # one value = faint limit
            mag = [0,mag]
        source_db['JMAGN'] = [float('0'+x) for x in numpy.ma.filled(source_db['J_2MASS'],'')]
        source_db['SELECT'][numpy.where(numpy.logical_and(source_db['JMAGN'] >= mag[0],source_db['JMAGN'] <= mag[1]))] += 1
        count+=1.
    if kwargs.get('hmag',False) != False:
        mag = kwargs['hmag']
        if not isinstance(mag,list):        # one value = faint limit
            mag = [0,mag]
        source_db['HMAGN'] = [float('0'+x) for x in numpy.ma.filled(source_db['H_2MASS'],'')]
        source_db['SELECT'][numpy.where(numpy.logical_and(source_db['HMAGN'] >= mag[0],source_db['HMAGN'] <= mag[1]))] += 1
        count+=1.
    if kwargs.get('kmag',False) != False:
        mag = kwargs['kmag']
        if not isinstance(mag,list):        # one value = faint limit
            mag = [0,mag]
        source_db['KMAGN'] = [float('0'+x) for x in numpy.ma.filled(source_db['KS_2MASS'],'')]
        source_db['SELECT'][numpy.where(numpy.logical_and(source_db['KMAGN'] >= mag[0],source_db['KMAGN'] <= mag[1]))] += 1
        count+=1.

# low surface gravity
    if (kwargs.get('lowg','') != ''):
#        source_db['LOWG'] = [not numpy.ma.is_masked(i) for i in source_db['GRAVITY_CLASS_OPTICAL']] or [not numpy.ma.is_masked(i) for i in source_db['GRAVITY_CLASS_NIR']]
        source_db['LOWG'] = [source_db['GRAVITY_CLASS_OPTICAL'][i]=='alpha' or source_db['GRAVITY_CLASS_OPTICAL'][i]=='beta' or source_db['GRAVITY_CLASS_OPTICAL'][i]=='gamma' or source_db['GRAVITY_CLASS_OPTICAL'][i]=='delta' or source_db['GRAVITY_CLASS_NIR'][i]=='INT-G' or source_db['GRAVITY_CLASS_NIR'][i]=='VL-G' or source_db['GRAVITY_CLASS_NIR'][i]=='LOW-G' for i in range(len(source_db))]
        source_db['SELECT'][numpy.where(source_db['LOWG'] == kwargs.get('lowg'))] += 1
        count+=1.

# specific gravity class
    flag = kwargs.get('gravity_class','')
    flag = kwargs.get('gravity',flag)
    if (flag != ''):
        source_db['SELECT'][numpy.where(numpy.ma.filled(source_db['GRAVITY_CLASS_OPTICAL'],'') == flag)] += 1
        source_db['SELECT'][numpy.where(numpy.ma.filled(source_db['GRAVITY_CLASS_NIR'],'') == flag)] += 1
        count+=1.

# young => member of a young cluster
    if (kwargs.get('young','') != ''):
        source_db['INCLUSTER'] = [not numpy.ma.is_masked(i) for i in source_db['CLUSTER']]
        source_db['SELECT'][numpy.where(source_db['INCLUSTER'] == kwargs.get('young'))] += 1
        count+=1.

# young => member of a young cluster
    if (kwargs.get('cluster','') != '' and isinstance(kwargs.get('cluster'),bool)):
        source_db['INCLUSTER'] = [not numpy.ma.is_masked(i) for i in source_db['CLUSTER']]
        source_db['SELECT'][numpy.where(source_db['INCLUSTER'] == kwargs.get('cluster'))] += 1
        count+=1.

# specific cluster
    if (kwargs.get('cluster','') != '' and isinstance(kwargs.get('cluster'),str)):
        source_db['CLUSTER_FLAG'] = [i.lower() == kwargs.get('cluster').lower() for i in numpy.ma.filled(source_db['CLUSTER'],'')]
        source_db['SELECT'][numpy.where(source_db['CLUSTER_FLAG'] == True)] += 1
        count+=1.

# giant
    if (kwargs.get('giant','') != ''):
#        kwargs['vlm'] = False
        source_db['GIANT'] = [not numpy.ma.is_masked(i) for i in source_db['LUMINOSITY_CLASS']]
        source_db['SELECT'][numpy.where(source_db['GIANT'] == kwargs.get('giant'))] += 1
        count+=1.

# luminosity class - this is not quite right
    if (kwargs.get('giant_class','') != ''):
#        if 'GIANT' not in source_db.keys():
#            source_db['GIANT'] = [not numpy.ma.is_masked(i) for i in source_db['LUMINOSITY_CLASS']]
        source_db['GIANT_FLAG'] = [i.lower() == kwargs.get('giant_class').lower() for i in numpy.ma.filled(source_db['LUMINOSITY_CLASS'],'')]
        source_db['SELECT'][numpy.where(source_db['GIANT_FLAG'] == True)] += 1
        count+=1.

# subdwarf
    if (kwargs.get('subdwarf','') != ''):
        source_db['SUBDWARF'] = [not numpy.ma.is_masked(i) for i in source_db['METALLICITY_CLASS']]
        source_db['SELECT'][numpy.where(source_db['SUBDWARF'] == kwargs.get('subdwarf'))] += 1
        count+=1.

# metallicity class
    if (kwargs.get('subdwarf_class','') != ''):
        source_db['SD_FLAG'] = [i.lower() == kwargs.get('subdwarf_class').lower() for i in numpy.ma.filled(source_db['METALLICITY_CLASS'],'')]
        source_db['SELECT'][numpy.where(source_db['SD_FLAG'] == True)] += 1
        count+=1.

# red - THIS NEEDS TO BE CHANGED
    if (kwargs.get('red','') != ''):
        source_db['RED'] = ['red' in i for i in numpy.ma.filled(source_db['LIBRARY'],'')]
        source_db['SELECT'][numpy.where(source_db['RED'] == kwargs.get('red'))] += 1
        count+=1.

# blue - THIS NEEDS TO BE CHANGED
    if (kwargs.get('blue','') != ''):
        source_db['BLUE'] = ['blue' in i for i in numpy.ma.filled(source_db['LIBRARY'],'')]
        source_db['SELECT'][numpy.where(source_db['BLUE'] == kwargs.get('blue'))] += 1
        count+=1.

# binaries
    if (kwargs.get('binary','') != ''):
        source_db['BINARY_FLAG'] = [not numpy.ma.is_masked(i) for i in source_db['BINARY']]
        source_db['SELECT'][numpy.where(source_db['BINARY_FLAG'] == kwargs.get('binary'))] += 1
        count+=1.

# spectral binaries
    if (kwargs.get('sbinary','') != ''):
        source_db['SBINARY_FLAG'] = [not numpy.ma.is_masked(i) for i in source_db['SBINARY']]
        source_db['SELECT'][numpy.where(source_db['SBINARY_FLAG'] == kwargs.get('sbinary'))] += 1
        count+=1.

# companions
    if (kwargs.get('companion','') != ''):
        source_db['COMPANION_FLAG'] = [not numpy.ma.is_masked(i) for i in source_db['COMPANION_NAME']]
        source_db['SELECT'][numpy.where(source_db['COMPANION_FLAG'] == kwargs.get('companion'))] += 1
        count+=1.

# white dwarfs
    if (kwargs.get('wd','') != ''):
        kwargs['vlm'] = False
        source_db['WHITEDWARF'] = [i == 'WD' for i in numpy.ma.filled(source_db['OBJECT_TYPE'],'')]
        source_db['SELECT'][numpy.where(source_db['WHITEDWARF'] == kwargs.get('wd'))] += 1
        count+=1.

# galaxies
    if (kwargs.get('galaxy','') != ''):
        kwargs['vlm'] = False
        source_db['GALAXY'] = [i == 'GAL' for i in numpy.ma.filled(source_db['OBJECT_TYPE'],'')]
        source_db['SELECT'][numpy.where(source_db['GALAXY'] == kwargs.get('galaxy'))] += 1
        count+=1.

# carbon stars
    if (kwargs.get('carbon','') != ''):
        kwargs['vlm'] = False
        source_db['CARBON'] = [i == 'C' for i in numpy.ma.filled(source_db['OBJECT_TYPE'],'')]
        source_db['SELECT'][numpy.where(source_db['CARBON'] == kwargs.get('carbon'))] += 1
        count+=1.

# peculiars
    if (kwargs.get('peculiar','') != ''):
#        kwargs['vlm'] = False
        source_db['PECULIAR'] = ['p' in i for i in numpy.ma.filled(source_db['LIT_TYPE'],'')]
        source_db['SELECT'][numpy.where(source_db['PECULIAR'] == kwargs.get('peculiar'))] += 1
        count+=1.

# VLM dwarfs
    if (kwargs.get('vlm','') != ''):
        if (kwargs.get('vlm') == True):
            source_db['SELECT'][numpy.where(source_db['OBJECT_TYPE'] == 'VLM')] += 1
            count+=1.
        if (kwargs.get('vlm') == False):
            source_db['SELECT'][numpy.where(source_db['OBJECT_TYPE'] != 'VLM')] += 1
            count+=1.

# select source keys
    if (count > 0):
        if (logic == 'and'):
            source_db['SELECT'] = numpy.floor(source_db['SELECT']/count)
        elif (logic == 'or'):
            source_db['SELECT'] = numpy.ceil(source_db['SELECT']/count)

        source_keys = source_db['SOURCE_KEY'][numpy.where(source_db['SELECT']==1)]
# no selection made on sources - choose everything
    else:
        source_keys = source_db['SOURCE_KEY']

#    print(count,numpy.max(source_db['SELECT']),len(source_db[:][numpy.where(source_db['SELECT']==1)]),len(source_keys))


# read in spectral database
#    spectral_db = ascii.read(SPLAT_PATH+DB_FOLDER+SPECTRA_DB, delimiter='\t',fill_values='-99.',format='tab')
#    spectral_db = fetchDatabase(SPLAT_PATH+DB_FOLDER+SPECTRA_DB)
    spectral_db = copy.deepcopy(DB_SPECTRA)
# having to force dtype here so SELECT remains an integer
    spectral_db['SELECT'] = Table.Column(numpy.zeros(len(spectral_db['DATA_KEY'])),dtype=int)
    count = 0.

    spectral_db['SOURCE_SELECT'] = [x in source_keys for x in spectral_db['SOURCE_KEY']]
#    print(spectral_db['SOURCE_KEY'][numpy.where(spectral_db['SOURCE_SELECT']==True)])

# search by source key
    datakey = kwargs.get('datakey',False)
    datakey = kwargs.get('data_key',datakey)
    if datakey != False:
        if not isinstance(datakey,list):
            datakey = [datakey]
        if isinstance(datakey[0],str):
            datakey = [int(i) for i in datakey]
        for s in datakey:
            spectral_db['SELECT'][numpy.where(spectral_db['DATA_KEY'] == s)] += 1
        count+=1.

# search by filename
    file = kwargs.get('file','')
    file = kwargs.get('filename',file)
    if (file != ''):
        if isinstance(file,str):
            file = [file]
        for f in file:
            spectral_db['SELECT'][numpy.where(spectral_db['DATA_FILE'] == f)] += 1
        count+=1.

# exclude by data key
    if kwargs.get('excludekey',False) != False:
        exkey = kwargs['excludekey']
        if len(exkey) > 0:
            if isinstance(exkey,str):
                exkey = [exkey]
            for f in exkey:
                spectral_db['SELECT'][numpy.where(spectral_db['DATA_KEY'] != f)] += 1
            count+=1.

# exclude by filename
    if kwargs.get('excludefile',False) != False:
        file = kwargs['excludefile']
        if len(file) > 0:
            if isinstance(file,str):
                file = [file]
            for f in file:
                spectral_db['SELECT'][numpy.where(spectral_db['DATA_FILE'] != f)] += 1
            count+=1.

# search by observation date range
    if kwargs.get('date',False) != False:
        date = kwargs['date']
        if not isinstance(date,list):
            date = [date,date]
        try:
            date = [float(properDate(x,output='YYYYMMDD')) for x in date]
#        if isinstance(date,str) or isinstance(date,long) or isinstance(date,float) or isinstance(date,int):
#            date = [float(date),float(date)]
#        elif isinstance(date,list):
#            date = [float(date[0]),float(date[-1])]
#        else:
        except:
            raise ValueError('\nCould not parse date input {}\n\n'.format(date))
        spectral_db['DATEN'] = [float(x) for x in spectral_db['OBSERVATION_DATE']]
        spectral_db['SELECT'][numpy.where(numpy.logical_and(spectral_db['DATEN'] >= date[0],spectral_db['DATEN'] <= date[1]))] += 1
        count+=1.

# search by S/N range
    if kwargs.get('snr',False) != False:
        snr = kwargs['snr']
        if not isinstance(snr,list):        # one value = minimum S/N
            snr = [float(snr),1.e9]
        spectral_db['SNRN'] = [float('0'+x) for x in spectral_db['MEDIAN_SNR']]
        spectral_db['SELECT'][numpy.where(numpy.logical_and(spectral_db['SNRN'] >= snr[0],spectral_db['SNRN'] <= snr[1]))] += 1
        count+=1.

# search by reference list
    if kwargs.get('data_reference',False) != False:
        drefer = kwargs['data_reference']
        if isinstance(drefer,str):
            drefer = [drefer]
        for r in drefer:
            spectral_db['SELECT'][numpy.where(spectral_db['DATA_REFERENCE'] == r)] += 1
        count+=1.

# search by spex type
    if spt_range != False and spt_type == 'SPEX_TYPE':
        if not isinstance(spt_range,list):        # one value = only this type
            spt_range = [spt_range,spt_range]
        if isinstance(spt_range[0],str):          # convert to numerical spt
            spt_range = [typeToNum(spt_range[0]),typeToNum(spt_range[1])]
        spectral_db['SPTN'] = [typeToNum(x) for x in spectral_db['SPEX_TYPE']]
        spectral_db['SELECT'][numpy.where(numpy.logical_and(spectral_db['SPTN'] >= spt_range[0],spectral_db['SPTN'] <= spt_range[1]))] += 1
        count+=1.

# exclude by filename
    if kwargs.get('ok',False) != False or kwargs.get('quality','').upper() == 'OK':
        spectral_db['SELECT'][numpy.where(spectral_db['QUALITY_FLAG'] == 'OK')] += 1
        count+=1.

# combine selection logically
    if (count > 0):
        if (logic == 'and'):
            spectral_db['SELECT'] = numpy.floor(spectral_db['SELECT']/count)
        else:
            spectral_db['SELECT'] = numpy.ceil(spectral_db['SELECT']/count)

    else:
        spectral_db['SELECT'] = numpy.ones(len(spectral_db['DATA_KEY']))


# limit access to public data for most users
#    print(count,numpy.max(spectral_db['SOURCE_SELECT']),numpy.max(spectral_db['SELECT']))
#    print(len(spectral_db[:][numpy.where(spectral_db['SELECT']==1)]))
#    print(len(spectral_db[:][numpy.where(spectral_db['SOURCE_SELECT']==True)]))
#    print(len(spectral_db[:][numpy.where(numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==True))]))
    if (not checkAccess() or kwargs.get('published',False) or kwargs.get('public',False)):
        spectral_db['SELECT'][numpy.where(spectral_db['PUBLISHED'] != 'Y')] = 0.

#    print(spectral_db['SOURCE_KEY'][numpy.where(spectral_db['SELECT']==1)])
#    print(spectral_db['SOURCE_KEY'][numpy.where(spectral_db['SOURCE_SELECT']==True)])

# no matches
#    print(count,numpy.max(spectral_db['SOURCE_SELECT']),numpy.max(spectral_db['SELECT']))
#    print(len(spectral_db[:][numpy.where(spectral_db['SELECT']==1)]))
#    print(len(spectral_db[:][numpy.where(spectral_db['SOURCE_SELECT']==True)]))
#    print(len(spectral_db[:][numpy.where(numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==True))]))
    if len(spectral_db[:][numpy.where(numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==True))]) == 0:
        if verbose: print('No spectra in the SPL database match the selection criteria')
        return Table()
    else:

# merge databases
        db = join(spectral_db[:][numpy.where(numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==True))],source_db,keys='SOURCE_KEY')

# sort output - default is by designation
        db.sort('DESIGNATION')
        if kwargs.get('sort',False) != False:
            if kwargs['sort'].upper() in db.colnames:
                db.sort(kwargs['sort'].upper())
            if kwargs['sort'].upper() == 'SNR': db.sort('MEDIAN_SNR')
        if kwargs.get('reverse',False) != False: db.reverse()

# select what to return
        if (ref == 'all'):
            outdb = db
        else:
            outdb = db[ref]

        return outdb


def _readAPOGEE(file,**kwargs):
    '''
    This assumes you have an ascap data file with data dimensions of [data],[uncertainty],[]
    '''

# make sure file is there
    if not os.path.exists(file):
        raise NameError('\nCould not find APOGEE file {}'.format(file))

    hdulist = fits.open(file)
    header = hdulist[0].header
    wave = 10.**(numpy.linspace(hdulist[0].header['CRVAL1'],hdulist[0].header['CRVAL1']+hdulist[0].header['NAXIS1']*hdulist[0].header['CDELT1'],num=hdulist[0].header['NAXIS1'],endpoint=False)-4.)
    flx = hdulist[1].data
    unc = hdulist[2].data
    mdl = hdulist[3].data
    return Spectrum()





def readSpectrum(*args,**kwargs):
    '''
    .. DOCS: will come back to this one
    '''
# keyword parameters
    folder = kwargs.get('folder','')
    catchSN = kwargs.get('catchSN',True)
    local = kwargs.get('local',True)
    online = kwargs.get('online',not local and checkOnline())
    local = not online
    url = kwargs.get('url',SPLAT_URL+DATA_FOLDER)
    kwargs['model'] = False

# filename
    file = kwargs.get('file','')
    file = kwargs.get('filename',file)
    if (len(args) > 0):
        file = args[0]
    kwargs['filename'] = file
    kwargs['model'] = False

# a filename must be passed
    if (kwargs['filename'] == ''):
        raise NameError('\nNeed to pass in filename to read in spectral data (readSpectrum)\n\n')

# first pass: check if file is local
    if online == False:
        file = kwargs['filename']
        if not os.path.exists(file):
            file = folder+os.path.basename(kwargs['filename'])
            if not os.path.exists(file):
#                print('Cannot find '+kwargs['filename']+' locally, trying online\n\n')
                local = False

# second pass: download file if necessary
    online = not local
    if online == True:
        file = checkOnline(url+kwargs['filename'])
        if file=='':
            raise NameError('\nCannot find file '+kwargs['filename']+' on SPLAT website\n\n')
        else:
# read in online file
            file = kwargs['filename']
            try:
#                file = TMPFILENAME+'.'+ftype
                if os.path.exists(os.path.basename(file)):
                    os.remove(os.path.basename(file))
#                open(os.path.basename(file), 'wb').write(urllib2.urlopen(url+file).read())
                open(os.path.basename(file), 'wb').write(requests.get(url+file).content)
#                print(file)
#                kwargs['filename'] = os.path.basename(file)
#               sp = Spectrum(**kwargs)
#                os.remove(os.path.basename(tmp))
#                return sp
            except:
                raise NameError('\nProblem reading in '+file+' from SPLAT website\n\n')

# determine which type of file
    ftype = file.split('.')[-1]

# fits file
    if (ftype == 'fit' or ftype == 'fits'):
        df = fits.open(file)
        with fits.open(file) as data:
            data.verify('silentfix+warn')
            if 'NAXIS3' in list(data[0].header.keys()):
                d = numpy.copy(data[0].data[0,:,:])
            else:
                d =  numpy.copy(data[0].data)
            header = data[0].header

# ascii file
    else:
        try:
            d = numpy.genfromtxt(file, comments='#', unpack=False, \
                missing_values = ('NaN','nan'), filling_values = (numpy.nan)).transpose()
        except ValueError:
            d = numpy.genfromtxt(file, comments=';', unpack=False, \
                 missing_values = ('NaN','nan'), filling_values = (numpy.nan)).transpose()
        header = fits.Header()      # blank header

# delete file if this was an online read
    if online and not local and os.path.exists(os.path.basename(file)):
        os.remove(os.path.basename(file))

# assign arrays to wave, flux, noise
    if len(d[:,0]) > len(d[0,:]): d = d.transpose()  # array is oriented wrong

# SDSS format for wavelength scale - in header and log format
    if kwargs.get('sdss',False) == True or (kwargs.get('waveheader',False) == True and kwargs.get('wavelog',False) == True):
        flux = d[0,:]
        if 'CRVAL1' in list(data[0].header.keys()) and 'CDELT1' in list(data[0].header.keys()):
            wave = 10.**(numpy.linspace(float(data[0].header['CRVAL1']),float(data[0].header['CRVAL1'])+len(flux)*float(data[0].header['CDELT1']),num=len(flux)))
        else: raise ValueError('\nCannot find CRVAL1 and CDELT1 keywords in header of fits file {}'.format(file))
        if len(d[:,0]) > 1:
            noise = d[1,:]
        else:
            noise = numpy.zeros(len(flux))
            noise[:] = numpy.nan
#  wavelength scale in header and linear format
    elif (kwargs.get('waveheader',False) == True and kwargs.get('wavelinear',False) == True):
        flux = d[0,:]
        if 'CRVAL1' in list(data[0].header.keys()) and 'CDELT1' in list(data[0].header.keys()):
            wave = numpy.linspace(float(data[0].header['CRVAL1']),float(data[0].header['CRVAL1'])+len(flux)*float(data[0].header['CDELT1']),num=len(flux))
        else: raise ValueError('\nCannot find CRVAL1 and CDELT1 keywords in header of fits file {}'.format(file))
        if len(d[:,0]) > 1:
            noise = d[1,:]
        else:
            noise = numpy.zeros(len(flux))
            noise[:] = numpy.nan
# wavelength is explicitly in data array 
    else:
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
#    w = numpy.where(numpy.isnan(flux) == True)
#    flux[w] = 0.

# remove all parts of spectrum that are nans
    w = numpy.where(numpy.logical_and(numpy.isnan(wave) == False,numpy.isnan(flux) == False))
    wave = wave[w]
    flux = flux[w]
    noise = noise[w]


# fix to catch badly formatted files where noise column is S/N
#    print(flux, numpy.median(flux))
    if (catchSN):
          w = numpy.where(flux > numpy.nanmedian(flux))
          if (numpy.nanmedian(flux[w]/noise[w]) < 1.):
              noise = flux/noise
              w = numpy.where(numpy.isnan(noise))
              noise[w] = numpy.nanmedian(noise)

# clean up
#    if url != '' and not local:
#        os.remove(os.path.basename(TMPFILENAME))


    return {'wave':wave,'flux':flux,'noise':noise,'header':header}





#####################################################
###############   CLASSIFICATION   ##################
#####################################################


def classifyByIndex(sp, *args, **kwargs):
    '''
    :Purpose: 

    Determine the spectral type and uncertainty for a spectrum based on indices. Makes use of published index-SpT relations from `Reid et al. (2001) <http://adsabs.harvard.edu/abs/2001AJ....121.1710R>`_; `Testi et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...552L.147T>`_; `Allers et al. (2007) <http://adsabs.harvard.edu/abs/2007ApJ...657..511A>`_; and `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_. Returns 2-element tuple containing spectral type (numeric or string) and uncertainty.

    Required Inputs:

    :param sp: Spectrum class object, which should contain wave, flux and noise array elements.

    Optional Inputs:

    :param set: named set of indices to measure and compute spectral type

        - *'burgasser'*: H2O-J, CH4-J, H2O-H, CH4-H, CH4-K from `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_ (default)
        - *'allers'*: H2O from `Allers et al. (2007) <http://adsabs.harvard.edu/abs/2007ApJ...657..511A>`_
        - *'reid'*:H2O-A and H2O-B from `Reid et al. (2001) <http://adsabs.harvard.edu/abs/2001AJ....121.1710R>`_
        - *'testi'*: sHJ, sKJ, sH2O_J, sH2O_H1, sH2O_H2, sH2O_K from `Testi et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...552L.147T>`_

    :param string: return spectral type as a string using typeToNum (default = False)
    :param round: rounds off to nearest 0.5 subtypes (default = False)
    :param allmeasures: Set to True to return all of the index values and individual subtypes  (default = False)
    :param remeasure: force remeasurement of indices (default = True)
    :param nsamples: number of Monte Carlo samples for error computation (default = 100)
    :param nloop: number of testing loops to see if spectral type is within a certain range (default = 5)

    :Example:
    >>> import splat
    >>> spc = splat.getSpectrum(shortname='0559-1404')[0]
    >>> splat.classifyByIndex(spc, string=True, set='burgasser', round=True)
        ('T4.5', 0.2562934083414341)

    '''

    str_flag = kwargs.get('string', True)
    verbose = kwargs.get('verbose', False)
    rnd_flag = kwargs.get('round', False)
    rem_flag = kwargs.get('remeasure', True)
    nsamples = kwargs.get('nsamples', 100)
    nloop = kwargs.get('nloop', 5)
    set = kwargs.get('set','burgasser')
    set = kwargs.get('ref',set)
    kwargs['set'] = set
    allowed_sets = ['burgasser','reid','testi','allers']
    if (set.lower() not in allowed_sets):
        print('\nWarning: index classification method {} not present; returning nan\n\n'.format(set))
        return numpy.nan, numpy.nan

# measure indices if necessary
    if (len(args) != 0):
        indices = args[0]

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

# Burgasser (2007, ApJ, 659, 655) calibration
    elif (set.lower() == 'burgasser'):
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

# Allers et al. (2013, ApJ, 657, 511)
    elif (set.lower() == 'allers'):
        if (rem_flag or len(args) == 0):
            kwargs['set'] = 'mclean'
            i1 = measureIndexSet(sp, **kwargs)
            kwargs['set'] = 'slesnick'
            i2 = measureIndexSet(sp, **kwargs)
            kwargs['set'] = 'allers'
            i3 = measureIndexSet(sp, **kwargs)
            if sys.version_info.major == 2:
                indices = dict(i1.items() + i2.items() + i3.items())
            else:
                indices = dict(i1.items() | i2.items() | i3.items())
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

# Aganze et al. 2015 (in preparation)
#    elif (set.lower() == 'aganze'):
#        if (rem_flag or len(args) == 0):
#            kwargs['set'] = 'geballe'
#            i1 = measureIndexSet(sp, **kwargs)
#            kwargs['set'] = 'slesnick'
#            i2 = measureIndexSet(sp, **kwargs)
#            kwargs['set'] = 'allers'
#            i3 = measureIndexSet(sp, **kwargs)
#            kwargs['set'] = 'burgasser'
#            i4 = measureIndexSet(sp, **kwargs)
#            kwargs['set'] = 'reid'
#            i5 = measureIndexSet(sp, **kwargs)
#            kwargs['set'] = 'tokunaga'
#            i6 = measureIndexSet(sp, **kwargs)
#            if sys.version_info.major == 2:
#                indices = dict(i1.items() + i2.items() + i3.items()+ i4.items() + i5.items() + i6.items())
#            else:
#                indices = dict(i1.items() | i2.items() | i3.items()| i4.items() | i5.items() | i6.items())
#        sptoffset = 0.0
#        sptfact = 1.0
#        coeffs = { \
#            'H2O': {'fitunc': 0.863, 'range': [15,23], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#            'coeff': [ -361.25130485, 1663.93768276, -2870.50724103,  2221.99873698, -638.03203556]}, \
#            'H2O-J': {'fitunc': 0.902, 'range': [15,23], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#            'coeff': [ -146.21144969 ,  632.34633568,  -1008.79681307,   678.80156994 , -137.92921741]}, \
#            'H2O-K': {'fitunc':0.973, 'range': [15,23], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#            'coeff': [-21366.79781425,  38630.25299752,  -25984.2424891  ,  7651.46728497,  -805.79462608]}, \
#            'K1': {'fitunc':0.878, 'range': [15,23], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#            'coeff': [  10.29493194 ,  -62.71016723 ,  115.76162692,   -60.72606292 ,  15.1905955 ]}, \
#            'K2': {'fitunc':0.934, 'range': [15,23], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#            'coeff': [ -44.78083424 , 225.58312733 ,-428.98225919 ,379.28205312 , -114.74469746]}, \
#            'H2O-1': {'fitunc':1.035, 'range': [15,23], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#            'coeff': [ -2999.69506898 , 11118.42653046 , -15340.87706264  ,9307.5183138, -2068.63608393]}, \
#            'H2O-B': {'fitunc':1.096, 'range': [15,23], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#            'coeff': [ -458.07448646 , 1547.35113353 , -1936.51451632 , 1041.95275566  , -178.50240834]}, \
#            'H2O-H': {'fitunc':1.041, 'range': [15,23], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#            'coeff': [ -767.21126974 , 2786.26168556 , -3762.93498987,   2211.62680244,  -451.54693932]}, \
#            'CH4-2.2': {'fitunc': 0.932, 'range': [15,23], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#            'coeff': [-331.74150369, 133.08406514  , -0.84614999  , 19.78717161  , 17.18479766]}}

    else:
        sys.stderr.write('\nWarning: '+set.lower()+' SpT-index relation not in classifyByIndex code\n\n')
        return numpy.nan, numpy.nan

    for index in coeffs.keys():
        if indices[index][1] > 0.:
            vals = numpy.polyval(coeffs[index]['coeff'],numpy.random.normal(indices[index][0],indices[index][1],nsamples))*sptfact
            coeffs[index]['spt'] = numpy.polyval(coeffs[index]['coeff'],indices[index][0])+sptoffset
            coeffs[index]['sptunc'] = (numpy.nanstd(vals)**2+coeffs[index]['fitunc']**2)**0.5
        else:
            coeffs[index]['spt'] = numpy.nan
            coeffs[index]['sptunc'] = numpy.nan

        if (coeffs[index]['spt'] < coeffs[index]['range'][0] or coeffs[index]['spt'] > coeffs[index]['range'][1] or numpy.isnan(coeffs[index]['spt'])):
            coeffs[index]['mask'] = 0.
        else:
            coeffs[index]['mask'] = 1.

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
            if (sptn < coeffs[index]['range'][0] or sptn > coeffs[index]['range'][1]):
                coeffs[index]['mask'] = 0

# report individual subtypes
    if verbose:
        for i in coeffs.keys():
            flg = '*'
            if coeffs[i]['mask'] == 0:
                flg = ''
            print('{}{} = {:.3f}+/-{:.3f} = SpT = {}+/-{}'.format(flg,i,indices[i][0],indices[i][1],typeToNum(coeffs[i]['spt']),coeffs[i]['sptunc']))

# round off to nearest 0.5 subtypes if desired
    if (rnd_flag):
        sptn = 0.5*numpy.around(sptn*2.)

# change to string if desired
    if (str_flag):
        spt = typeToNum(sptn,uncertainty=sptn_e)
    else:
        spt = sptn

    if kwargs.get('allmeasures',False):
        output = {}
        for k in coeffs.keys():
            output[k] = {'spt': coeffs[k]['spt'], 'spt_e': coeffs[k]['sptunc'], 'index': indices[k][0], 'index_e': indices[k][0]}
        output['result'] = (spt,sptn_e)
        return output
    else:
        return spt, sptn_e



def classifyByStandard(sp, *args, **kwargs):
    '''
    :Purpose: Determine the spectral type and uncertainty for a
                spectrum by direct comparison to defined spectral standards.  
                Dwarf standards span M0-T9 and include the standards listed in
                `Burgasser et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...637.1067B>`_, `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_ and `Cushing et al. (2011) <http://adsabs.harvard.edu/abs/2011ApJ...743...50C>`_. 
                Comparison to subdwarf and extreme subdwarf standards may also be done.
                Returns the best
                match or an F-test weighted mean and uncertainty. There is an option
                to follow the procedure of `Kirkpatrick et al. (2010)
                <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_, fitting only in
                the 0.9-1.4 micron region.

    :Output: A tuple listing the best match standard and uncertainty based on F-test weighting and systematic uncertainty of 0.5 subtypes

    :param sp: Spectrum class object, which should contain wave, flux and
               noise array elements.
    :param sp: required
    :param sptrange: Set to the spectral type range over which comparisons should be made, can be a two-element array of strings or numbers
    :type sptrange: optional, default = ['M0','T9']
    :param statistic: string defining which statistic to use in comparison; available options are:

            - *'chisqr'*: compare by computing chi squared value (requires spectra with noise values)
            - *'stddev'*: compare by computing standard deviation
            - *'stddev_norm'*: compare by computing normalized standard deviation
            - *'absdev'*: compare by computing absolute deviation

    :type statistic: optional, default = 'chisqr'
    :param method: set to ``'kirkpatrick'`` to follow the `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_ method, fitting only to the 0.9-1.4 micron band
    :type method: optional, default = ''
    :param best: Set to True to return the best fit standard type
    :type best: optional, default = True
    :param average: Set to True to return an chi-square weighted type only
    :type average: optional, default = False
    :param compareto: Set to the single standard (string or number) you want to compare to 
    :type compareto: optional, default = None
    :param plot: Set to True to generate a plot comparing best fit template to source; can also set keywords associated with plotSpectrum_ routine 
    :type plot: optional, default = False
    :param string: return spectral type as a string
    :type string: optional, default = True
    :param verbose: Set to True to give extra feedback
    :type verbose: optional, default = False

    Users can also set keyword parameters defined in plotSpectrum_ and compareSpectra_ routine.

    .. _compareSpectra : api.html#splat.compareSpectra

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(lucky=True)[0]
    >>> result = splat.classifyByStandard(sp,verbose=True)
        Using dwarf standards
        Type M3.0: statistic = 5763368.10355, scale = 0.000144521824721
        Type M2.0: statistic = 5613862.67356, scale = 0.000406992798674
        Type T8.0: statistic = 18949835.2087, scale = 9.70960919364
        Type T9.0: statistic = 21591485.163, scale = 29.1529786804
        Type L8.0: statistic = 3115605.62687, scale = 1.36392504072
        Type L9.0: statistic = 2413450.79206, scale = 0.821131769522
        ...
        Best match to L1.0 spectral standard
        Best spectral type = L1.0+/-0.5
    >>> result
        ('L1.0', 0.5)
    >>> splat.classifyByStandard(sp,sd=True,average=True)
        ('sdL0.0:', 1.8630159149200021)
    '''

    verbose = kwargs.get('verbose',False)
    best_flag = kwargs.get('best',True)
    average_flag = kwargs.get('average',not best_flag)
    best_flag = not average_flag
    statistic = kwargs.get('statistic','chisqr')
    statistic = kwargs.get('stat',statistic)
    sptrange = kwargs.get('sptrange',[10,39])
    sptrange = kwargs.get('range',sptrange)
    sptrange = kwargs.get('spt',sptrange)
    if not isinstance(sptrange,list):
        sptrange = [sptrange,sptrange]
    if (isinstance(sptrange[0],str) != False):
        sptrange = [typeToNum(sptrange[0]),typeToNum(sptrange[1])]
    unc_sys = 0.5       # assumed systematic uncertainty


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

# assign subclasses
    if kwargs.get('sd',False):
        stds = STDS_SD_SPEX
        subclass = 'sd'
        initiateStandards(sd=True)
        if verbose:
            print('Using subdwarf standards')
    elif kwargs.get('esd',False):
        stds = STDS_ESD_SPEX
        subclass = 'esd'
        initiateStandards(esd=True)
        if verbose:
            print('Using extreme subdwarf standards')
    else:
        stds = STDS_DWARF_SPEX
        initiateStandards()
        subclass = ''
        if verbose:
            print('Using dwarf standards')

# select desired spectral range
    spt_allowed = numpy.array([typeToNum(s) for s in stds.keys()])
    spt_sample = spt_allowed[numpy.where(spt_allowed >= sptrange[0])]
    spt_sample = spt_sample[numpy.where(spt_sample <= sptrange[1])]

# determine comparison range based on method
    if (kwargs.get('method','') == 'kirkpatrick'):
        fit_ranges = [[0.9,1.4]]         # as prescribed in Kirkpatrick et al. 2010, ApJS,
    else:
        fit_ranges = [[0.7,2.45]]       # by default, compare whole spectrum
    fit_ranges = kwargs.get('fitrange',fit_ranges)
    fit_ranges = kwargs.get('fitrng',fit_ranges)
    fit_ranges = kwargs.get('comprange',fit_ranges)
    fit_ranges = kwargs.get('comprng',fit_ranges)
    if not isinstance(fit_ranges[0],list):
        fit_ranges = [fit_ranges]


# compute fitting statistics
    stat = []
    sspt = []

    for t in spt_sample:
        chisq,scale = compareSpectra(sp,stds[typeToNum(t,subclass=subclass)],fit_ranges=fit_ranges,statistic=statistic,novar2=True)
        stat.append(chisq)
        sspt.append(t)
        if (verbose):
            print('Type {}: statistic = {}, scale = {}'.format(typeToNum(t,subclass=subclass), chisq, scale))

# list of sorted standard files and spectral types
    sorted_stdsptnum = [x for (y,x) in sorted(zip(stat,sspt))]

# select either best match or an ftest-weighted average
# note that these are NUMBERS
    if (best_flag or len(stat) == 1):
        sptn = sorted_stdsptnum[0]
        sptn_e = unc_sys
    else:
        try:
            st = stat.value
        except:
            st = stat
        if numpy.isnan(numpy.median(sp.noise)):
            mean,var = weightedMeanVar(sspt,st)
        else:
            mean,var = weightedMeanVar(sspt,st,method='ftest',dof=sp.dof)
        if (var**0.5 < 1.):
            sptn = numpy.round(mean*2)*0.5
        else:
            sptn = numpy.round(mean)
        sptn_e = (unc_sys**2+var)**0.5

# string or not?
    if (kwargs.get('string', True) == True):
        output_spt = typeToNum(sptn,uncertainty=sptn_e,subclass=subclass)
    else:
        output_spt = sptn

    if verbose:
        print('\nBest match to {} spectral standard'.format(typeToNum(sorted_stdsptnum[0],subclass=subclass)))
        print('Best spectral type = {}+/-{}'.format(output_spt,sptn_e))

# plot spectrum compared to best spectrum
    if (kwargs.get('plot',False) != False):
#        spstd = Spectrum(file=sorted_stdfiles[0])
#        print(typeToNum(sorted_stdsptnum[0],subclass=subclass))
        spstd = getStandard(typeToNum(sorted_stdsptnum[0],subclass=subclass))
        chisq,scale = compareSpectra(sp,spstd,fit_ranges=fit_ranges,statistic=statistic)
        spstd.scale(scale)
        if kwargs.get('colors',False) == False:
            kwargs['colors'] = ['k','r','b']
        if kwargs.get('labels',False) == False:
            kwargs['labels'] = [sp.name,'{} Standard'.format(typeToNum(sorted_stdsptnum[0],subclass=subclass)),'Difference']
        from .plot import plotSpectrum
        if kwargs.get('difference',True):
            kwargs['labels'].append('Difference')
            pl = plotSpectrum(sp,spstd,sp-spstd,**kwargs)
        else:
            pl = plotSpectrum(sp,spstd,**kwargs)

    return output_spt, sptn_e




def classifyByTemplate(sp, *args, **kwargs):
    '''
    :Purpose: Determine the spectral type and uncertainty for a
                spectrum by direct comparison to a large set of spectra in
                the library. Returns a dictionary with the best spectral type (F-test weighted mean and
                uncertainty), and arrays for the N best-matching Spectrum objects, scale factors, spectral types and comparison statistics. 
                There is an option to follow the procedure of
                `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_,
                fitting only in the 0.9-1.4 micron region.
                It is strongly encouraged that users winnow down the templates used in the comparison
                by selecting templates using the searchLibrary_ options or optionally the ``set`` parameter. 

    :Output: A dictionary containing the following keys:

                    - **result**: a tuple containing the spectral type and its uncertainty based on F-test statistic
                    - **statistic**: array of N best statistical comparison values
                    - **scale**: array of N best optimal scale factors
                    - **spectra**: array of N best Spectrum objects
                    - **spt**: array of N best spectral types

    :param sp: Spectrum class object, which should contain wave, flux and
               noise array elements.
    :param sp: required
    :param statistic: string defining which statistic to use in comparison; available options are:

            - *'chisqr'*: compare by computing chi squared value (requires spectra with noise values)
            - *'stddev'*: compare by computing standard deviation
            - *'stddev_norm'*: compare by computing normalized standard deviation
            - *'absdev'*: compare by computing absolute deviation

    :type statistic: optional, default = 'chisqr'
    :param select: string defining which spectral template set you want to compare to; several options which can be combined:

            - *m dwarf*: fit to M dwarfs only
            - *l dwarf*: fit to M dwarfs only
            - *t dwarf*: fit to M dwarfs only
            - *vlm*: fit to M7-T9 dwarfs
            - *optical*: only optical classifications
            - *high sn*: median S/N greater than 100
            - *young*: only young/low surface gravity dwarfs
            - *companion*: only companion dwarfs
            - *subdwarf*: only subdwarfs
            - *single*: only dwarfs not indicated a binaries
            - *spectral binaries*: only dwarfs indicated to be spectral binaries
            - *standard*: only spectral standards (Note: use classifyByStandard_ instead)

    :type select: optional, default = ''
    :param method: set to ``'kirkpatrick'`` to follow the `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_ method, fitting only to the 0.9-1.4 micron band
    :type method: optional, default = ''
    :param best: Set to True to return only the best fit template type
    :type best: optional, default = False
    :param nbest: Set to the number of best fitting spectra to return
    :type nbest: optional, default = 1
    :param maxtemplates: Set to the maximum number of templates that should be fit
    :type maxtemplates: optional, default = 100
    :param force: By default, classifyByTemplate won't proceed if you have more than 100 templates; set this parameter to True to ignore that constraint
    :type force: optional, default = False
    :param plot: Set to True to generate a plot comparing best fit template to source; can also set keywords associated with plotSpectrum_ routine 
    :type plot: optional, default = False
    :param string: return spectral type as a string
    :type string: optional, default = True
    :param verbose: give lots of feedback
    :type verbose: optional, default = False

    Users can also set keyword parameters defined in plotSpectrum_ and searchLibrary_ routines

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
    >>> result = splat.classifyByTemplate(sp,string=True,spt=[24,26],nbest=5)
        Too many templates (171) for classifyByTemplate; set force=True to override this
    >>> result = splat.classifyByTemplate(sp,string=True,spt=[24,26],snr=50,nbest=5)
        Comparing to 98 templates
        LHS 102B L5.0 10488.1100432 11.0947838116
        2MASSI J0013578-223520 L4.0 7037.37441677 136.830522173
        SDSS J001608.44-004302.3 L5.5 15468.6209466 274.797693706
        2MASSI J0028394+150141 L4.5 63696.1897668 187.266152375
        ...
        Best match = DENIS-P J153941.96-052042.4 with spectral type L4:
        Mean spectral type = L4.5+/-0.718078660103
    >>> result
        {'result': ('L4.5', 0.71807866010293797),
         'scale': [3.0379089778408642e-14,
          96.534933767992072,
          3.812718429200959,
          2.9878801833735986e-14,
          3.0353579048704484e-14],
         'spectra': [Spectrum of DENIS-P J153941.96-052042.4,
          Spectrum of 2MASSI J0443058-320209,
          Spectrum of SDSSp J053951.99-005902.0,
          Spectrum of 2MASSI J1104012+195921,
          Spectrum of 2MASS J17502484-0016151],
         'spt': [24.0, 25.0, 25.0, 24.0, 25.5],
         'statistic': [<Quantity 2108.997879536768>,
          <Quantity 2205.640664932956>,
          <Quantity 2279.316858783139>,
          <Quantity 2579.0089210846527>,
          <Quantity 2684.003187310027>]}

    .. _classifyByStandard : api.html#splat.classifyByStandard
    .. _searchLibrary : api.html#splat_db.searchLibrary
    .. _plotSpectrum : api.html#splat_plot.plotSpectrum

    '''

#
    spt_type = kwargs.get('spt_type','literature')
    spt = kwargs.get('spt',[10.,39.9])
    spt = kwargs.get('spt_range',spt)
    nbest = kwargs.get('nbest',1)
    verbose = kwargs.get('verbose',False)
    published = kwargs.get('published','')
    published = kwargs.get('public',published)
    statistic = kwargs.get('statistic','chisqr')
    statistic = kwargs.get('stat',statistic)
    force = kwargs.get('force',False)
    maxtemplates = kwargs.get('maxtemplates',100)
    select = kwargs.get('select','')
    select = kwargs.get('set',select)
#   placeholder for a systematic uncertainty term
    unc_sys = 0.
    if (kwargs.get('method','') == 'kirkpatrick'):
        fit_ranges = [[0.9,1.4]]         # as prescribed in Kirkpatrick et al. 2010, ApJS,
    else:
        fit_ranges = [[0.7,2.45]]       # by default, compare whole spectrum
    fit_ranges = kwargs.get('fitrange',fit_ranges)
    fit_ranges = kwargs.get('fitrng',fit_ranges)
    fit_ranges = kwargs.get('comprange',fit_ranges)
    fit_ranges = kwargs.get('comprng',fit_ranges)
    if not isinstance(fit_ranges[0],list):
        fit_ranges = [fit_ranges]

#  canned searches
#  constrain spectral types
    if ('lit' in spt_type.lower()):
        spt_type = 'LIT_TYPE'
    elif ('opt' in spt_type.lower() or 'optical' in select):
        spt_type = 'OPT_TYPE'
    elif ('nir' in spt_type.lower()):
        spt_type = 'NIR_TYPE'
    else:
        spt_type = 'LIT_TYPE'

    if ('m dwarf' in select.lower() or kwargs.get('mdwarf',False)):
        spt = [numpy.max([10,spt[0]]),numpy.min([19.9,spt[-1]])]
    if ('l dwarf' in select.lower() or kwargs.get('ldwarf',False)):
        spt = [numpy.max([20,spt[0]]),numpy.min([29.9,spt[-1]])]
    if ('t dwarf' in select.lower() or kwargs.get('tdwarf',False)):
        spt = [numpy.max([30,spt[0]]),numpy.min([39.9,spt[-1]])]
    if ('vlm' in select.lower() or kwargs.get('vlm',False)):
        spt = [numpy.max([17,spt_range[0]]),numpy.min([39.9,spt_range[-1]])]

#  constrain S/N
    snr = 0.
    if ('high sn' in select.lower()):
        snr = 100.
    snr = kwargs.get('snr',snr)

#  don't compare to same spectrum
    try:
        excludefile = [sp.filename]
    except:
        excludefile = []
    if kwargs.get('excludefile',False) != False:
        e = kwargs.get('excludefile')
        if isinstance(e,list):
            excludefile.extend(e)
        else:
            excludefile.append(e)
    try:
        excludekey = [sp.data_key]
    except:
        excludekey = []
    if kwargs.get('excludekey',False) != False:
        e = kwargs.get('excludekey')
        if isinstance(e,list):
            excludekey.extend(e)
        else:
            excludekey.append(e)
    try:
        excludeshortname = [sp.shortname]
    except:
        excludeshortname = []
    if kwargs.get('excludeshortname',False) != False:
        e = kwargs.get('excludeshortname')
        if isinstance(e,list):
            excludeshortname.extend(e)
        else:
            excludeshortname.append(e)
#    print(excludefile, excludekey, excludeshortname)

# other classes
    giant = ''
    if 'giant' in select.lower() or kwargs.get('giant',False):
        giant = True
    if 'not giant' in select.lower():
        giant = False
    companion = ''
    if 'companion' in select.lower() or kwargs.get('companion',False):
        companion = True
    if 'not companion' in select.lower():
        companion = False
    young = ''
    if 'young' in select.lower() or kwargs.get('young',False):
        young = True
    if 'not young' in select.lower():
        young = False
    binary = ''
    if 'binary' in select.lower() or kwargs.get('binary',False):
        binary = True
    if 'not binary' in select.lower():
        binary = False
    spbinary = ''
    if 'spectral binary' in select.lower() or kwargs.get('sbinary',False):
        spbinary = True
    if 'not spectral binary' in select.lower():
        spbinary = False

# REARRANGE THIS - SEND IN KWARGS WITH OUTPUT, LOGIC SET, AND THE REST ARE UP TO USER?

    lib = searchLibrary(excludefile=excludefile,excludekey=excludekey,excludeshortname=excludeshortname, \
        snr=snr,spt_type=spt_type,spt=spt,published=published, \
        giant=giant,companion=companion,young=young,binary=binary,spbinary=spbinary,output='all',logic='and')

# first search for the spectra desired - parameters are set by user
    if len(lib) == 0:
        print('\nNo templates available for comparison\n\n')
        return numpy.nan, numpy.nan

    if len(lib) > maxtemplates and force == False:
        print('\nToo many templates ({}) for classifyByTemplate; set force=True to override this\n\n'.format(len(lib)))
        return numpy.nan, numpy.nan

    files = lib['DATA_FILE']
    dkey = lib['DATA_KEY']
    sspt = [typeToNum(s) for s in lib[spt_type]]

    if (verbose):
        print('\nComparing to {} templates\n'.format(len(files)))
        if len(files) > 100:
            print('This may take some time!\n\n'.format(len(files)))

# do comparison
    stat = []
    scl = []
    for i,d in enumerate(dkey):

# INSERT TRY STATEMNT HERE?

        s = Spectrum(idkey=d)
        stt,scale = compareSpectra(sp,s,fit_ranges=fit_ranges,statistic=statistic,novar2=True,*kwargs)
        stat.append(stt)
        scl.append(scale)
        if (verbose):
            print(keySpectrum(d)['NAME'][0], typeToNum(sspt[i]), stt, scale)

# list of sorted standard files and spectral types
    sorted_dkey = [x for (y,x) in sorted(zip(stat,dkey))]
    sorted_spt = [x for (y,x) in sorted(zip(stat,sspt))]
    sorted_scale = [x for (y,x) in sorted(zip(stat,scl))]

# select either best match or an ftest-weighted average
    if (kwargs.get('best',False) or len(stat) == 1):
        sptn = sorted_spt[0]
        sptn_e = unc_sys
    else:
        mean,var = weightedMeanVar(sspt,stat,method='ftest',dof=sp.dof)
# allow 1/2 subtypes if uncertainty is less than 1.0
        if (var**0.5 < 1.):
            sptn = numpy.round(mean*2.)*0.5
        else:
            sptn = numpy.round(mean)
        sptn_e = (unc_sys**2+var)**0.5

# plot spectrum compared to best spectrum
    if (kwargs.get('plot',False) != False):
        s = Spectrum(idkey=sorted_dkey[0])
#        chisq,scale = compareSpectra(s,sp,fit_ranges=[comprng],stat='chisqr',novar2=True)
        s.scale(sorted_scale[0])
        kwargs['legend'] = [sp.name,s.name]
        kwargs['colors'] = ['k','r','b']
        from .plot import plotSpectrum
        plotSpectrum(sp,s,sp-s,**kwargs)

# string or not?
    if (kwargs.get('string', True) == True):
        output_spt = typeToNum(sptn,uncertainty=sptn_e)
    else:
        output_spt = sptn

    if verbose:
        s = Spectrum(idkey=sorted_dkey[0])
        print('\nBest match = {} with spectral type {}'.format(s.name,s.lit_type))
        print('Mean spectral type = {}+/-{}'.format(output_spt,sptn_e))

# return dictionary of results
    return {'result': (output_spt,sptn_e), \
        'statistic': sorted(stat)[0:nbest], 'spt': sorted_spt[0:nbest], \
        'scale': sorted_scale[0:nbest], \
        'spectra': [Spectrum(idkey=d) for d in sorted_dkey[0:nbest]]}



def classifyGravity(sp, *args, **kwargs):
    '''
    :Purpose: Determine the gravity classification of a brown dwarf using the method of `Allers & Liu (2013) <http://adsabs.harvard.edu/abs/2013ApJ...772...79A>`_. 

    :param sp: Spectrum class object, which should contain wave, flux and
               noise array elements. Must be between M6.0 and L7.0.
    :type sp: required
    :param spt: spectral type of ``sp``. Must be between M6.0 and L7.0
    :type spt: optional, default = False
    :param indices: specify indices set using ``measureIndexSet``.
    :type indices: optional, default = 'allers'
    :param plot: Set to True to plot sources against closest dwarf spectral standard
    :type plot: optional, default = False
    :param allscores: Set to True to return a dictionary containing the gravity scores from individual indices
    :type allscores: optional, default = False
    :param verbose: Give feedback while computing
    :type verbose: optional, default = False

    :Output: Either a string specifying the gravity classification or a dictionary specifying the gravity scores for each index

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
    >>> splat.classifyGravity(sp)
        FLD-G
    >>> result = splat.classifyGravity(sp, allscores = True, verbose=True)
        Gravity Classification:
            SpT = L4.0
            VO-z: 1.012+/-0.029 => 0.0
            FeH-z: 1.299+/-0.031 => 1.0
            H-cont: 0.859+/-0.032 => 0.0
            KI-J: 1.114+/-0.038 => 1.0
            Gravity Class = FLD-G
    >>> result
        {'FeH-z': 1.0,
         'H-cont': 0.0,
         'KI-J': 1.0,
         'VO-z': 0.0,
         'gravity_class': 'FLD-G',
         'score': 0.5,
         'spt': 'L4.0'}
    '''

    verbose = kwargs.get('verbose',False)

# Chart for determining gravity scores based on gravity sensitive
# indices as described in the Allers and Liu paper.
# The key to the overall indices dictionary is each index name.
# The key to each index dictionary are the spectral types, which
# contain the limits for each gravity score.
# To access a value do the following: print grav['FeH-z']['M5'][0]
# which should return 'nan'

# Note: alternate method is Canty et al. (2013, MNRAS, 435, 2650)
# H2(K) index: median[2.16-2.18]/median[2.23-2.25]

    grav = {\
        'FeH-z':{'M5.0':[numpy.nan,numpy.nan],'M6.0':[1.068,1.039],'M7.0':[1.103,1.056],'M8.0':[1.146,1.074],'M9.0': [1.167,1.086],'L0.0': [1.204,1.106],'L1.0':[1.252,1.121],'L2.0':[1.298,1.142],'L3.0': [1.357,1.163],'L4.0': [1.370,1.164],'L5.0': [1.258,1.138],'L6.0': [numpy.nan,numpy.nan],'L7.0': [numpy.nan,numpy.nan]},\
        'VO-z': {'M5.0':[numpy.nan,numpy.nan],'M6.0':[numpy.nan,numpy.nan],'M7.0': [numpy.nan,numpy.nan],'M8.0': [numpy.nan,numpy.nan],'M9.0': [numpy.nan,numpy.nan],'L0.0': [1.122,1.256],'L1.0': [1.112,1.251],'L2.0': [1.110,1.232],'L3.0': [1.097,1.187],'L4.0': [1.073,1.118],'L5.0': [numpy.nan,numpy.nan],'L6.0': [numpy.nan,numpy.nan],'L7.0': [numpy.nan,numpy.nan]},\
        'KI-J': {'M5.0': [numpy.nan,numpy.nan], 'M6.0': [1.042,1.028], 'M7.0': [1.059,1.036],'M8.0': [1.077,1.046],'M9.0': [1.085,1.053],'L0.0': [1.098,1.061],'L1.0': [1.114,1.067],'L2.0': [1.133,1.073],'L3.0': [1.135,1.075],'L4.0': [1.126,1.072],'L5.0': [1.094,1.061],'L6.0': [numpy.nan,numpy.nan],'L7.0': [numpy.nan,numpy.nan]},\
        'H-cont': {'M5.0': [numpy.nan,numpy.nan], 'M6.0': [.988,.994], 'M7.0': [.981,.990],'M8.0': [.963,.984],'M9.0': [.949,.979],'L0.0': [.935,.972],'L1.0': [.914,.968],'L2.0': [.906,.964],'L3.0': [.898,.960],'L4.0': [.885,.954],'L5.0': [.869,.949],'L6.0': [.874,.950],'L7.0': [0.888,0.952]}}

# Calculate Allers indices and their uncertainties
    ind = kwargs.get('indices',False)
    if ind == False:
        ind = measureIndexSet(sp,set='allers')

# Determine the object's NIR spectral type and its uncertainty
    sptn = kwargs.get('spt',False)
    if sptn == False:
        sptn, spt_e = classifyByIndex(sp,string=False,set='allers')
        if numpy.isnan(sptn):
            print('Spectral type could not be determined from indices')
            return ''
    if isinstance(sptn,str):
        sptn = typeToNum(sptn)
    Spt = typeToNum(numpy.round(sptn))

#Check whether the NIR SpT is within gravity sensitive range values
    if ((sptn < 16.0) or (sptn > 27.0)):
        print('Spectral type '+typeToNum(sptn)+' outside range for gravity classification')
        return ''

# print spt if verbose
    if verbose:
        print('\nGravity Classification:\n\tSpT = {}'.format(Spt))

#Creates an empty array with dimensions 4x1 to fill in later with 5 gravscore values
    gravscore = {'spt': Spt}
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
                print('\t{}: {:.3f}+/-{:.3f} => {}'.format(k,ind[k][0], ind[k][1], val))
        if k == 'FeH-z' or k=='KI-J':
            if numpy.isnan(grav[k][Spt][0]):
                val = numpy.nan
            if ind[k][0] <= grav[k][Spt][0]:
                val = 1.0
            if ind[k][0] <= grav[k][Spt][1]:
                val = 2.0
            if verbose:
                print('\t{}: {:.3f}+/-{:.3f} => {}'.format(k,ind[k][0], ind[k][1], val))
        gravscore[k] = val
        medgrav.append(val)

# determine median score, or mean if even
    if (len(numpy.where(numpy.isnan(medgrav) == False))%2 == 0):
        gravscore['score'] = numpy.nanmean(medgrav)
    else:
        gravscore['score'] = numpy.nanmedian(medgrav)

    if gravscore['score'] <= 0.5:
       gravscore['gravity_class'] = 'FLD-G'
    elif gravscore['score'] > 0.5 and gravscore['score'] < 1.5:
       gravscore['gravity_class'] = 'INT-G'
    elif gravscore['score'] >= 1.5:
       gravscore['gravity_class'] = 'VL-G'
    else:
       gravscore['gravity_class'] = 'UNKNOWN'

# print spt if verbose
    if verbose:
        print('\tGravity Class = {}\n'.format(gravscore['gravity_class']))


# plot spectrum against standard
    if (kwargs.get('plot',False) != False):
        spt,unc = classifyByStandard(sp,compareto=Spt,method='kirkpatrick',**kwargs)

# return gravity class or entire dictionary
    if (kwargs.get('allscores',False) == False):
        return gravscore['gravity_class']
    else:
        return gravscore



def compareSpectra(sp1, sp2, *args, **kwargs):
    '''
    :Purpose: Compare two spectra against each other using a pre-selected statistic. Returns the value of the desired statistic as well as the optimal scale factor. Minimum possible value for statistic is 1.e-9.

    :param sp1: First spectrum class object, which sets the wavelength scale
    :type sp1: required
    :param sp2: Second spectrum class object, interpolated onto the wavelength scale of sp1
    :type sp2: required
    :param statistic: string defining which statistic to use in comparison; available options are:

            - *'chisqr'*: compare by computing chi squared value (requires spectra with noise values)
            - *'stddev'*: compare by computing standard deviation
            - *'stddev_norm'*: compare by computing normalized standard deviation
            - *'absdev'*: compare by computing absolute deviation

    :type statistic: optional, default = 'chisqr'
    :param fit_ranges: 2-element array or nested array of 2-element arrays specifying the wavelength ranges to be used for the fit, assumed to be measured in microns. This is effectively the opposite of mask_ranges.
    :type fit_ranges: optional, default = [0.65,2.45]
    :param weights: Array specifying the weights for individual wavelengths; must be an array with length equal to the wavelength scale of ``sp1``; need not be normalized
    :type weights: optional, default = [1, ..., 1] for len(sp1.wave)
    :param mask_ranges: Multi-vector array setting wavelength boundaries for masking data, assumed to be in microns
    :type mask_ranges: optional, default = None
    :param mask: Array specifiying which wavelengths to mask; must be an array with length equal to the wavelength scale of ``sp1`` with only 0 (OK) or 1 (mask).
    :type mask: optional, default = [0, ..., 0] for len(sp1.wave)
    :param mask_telluric: Set to True to mask pre-defined telluric absorption regions
    :type mask_telluric: optional, default = False
    :param mask_standard: Like ``mask_telluric``, with a slightly tighter cut of 0.80-2.35 micron
    :type mask_standard: optional, default = False
    :param novar2: Set to True to compute statistic without considering variance of ``sp2``
    :type novar2: optional, default = True
    :param plot: Set to True to plot ``sp1`` with scaled ``sp2`` and difference spectrum overlaid
    :type plot: optional, default = False
    :param verbose: Set to True to report things as you're going along
    :type verbose: optional, default = False

    :Example:
    >>> import splat
    >>> import numpy
    >>> sp1 = splat.getSpectrum(shortname = '2346-3153')[0]
        Retrieving 1 file
    >>> sp2 = splat.getSpectrum(shortname = '1421+1827')[0]
        Retrieving 1 file
    >>> sp1.normalize()
    >>> sp2.normalize()    
    >>> splat.compareSpectra(sp1, sp2, statistic='chisqr')
        (<Quantity 19927.74527822856>, 0.94360732593223595)
    >>> splat.compareSpectra(sp1, sp2, statistic='stddev')
        (<Quantity 3.0237604611215705 erg2 / (cm4 micron2 s2)>, 0.98180983971456637)
    >>> splat.compareSpectra(sp1, sp2, statistic='absdev')
        (<Quantity 32.99816249949072 erg / (cm2 micron s)>, 0.98155779612333172)
    >>> splat.compareSpectra(sp1, sp2, statistic='chisqr', novar2=False)
        (<Quantity 17071.690727945213>, 0.94029474635786015)
    '''
    weights = kwargs.get('weights',numpy.zeros(len(sp1.wave))) # these will be set to 1 later
    fit_ranges = kwargs.get('fit_ranges',[spex_wave_range.value])
    fit_ranges = kwargs.get('fit_range',fit_ranges)
    fit_ranges = kwargs.get('fitrange',fit_ranges)
    fit_ranges = kwargs.get('fitrng',fit_ranges)
    fit_ranges = kwargs.get('comprange',fit_ranges)
    fit_ranges = kwargs.get('comprng',fit_ranges)
#    mask_ranges = kwargs.get('mask_ranges',[])
#    mask_standard = kwargs.get('mask_standard',False)
#    mask_telluric = kwargs.get('mask_telluric',mask_standard)
    var_flag = kwargs.get('novar2',True)
    statistic = kwargs.get('statistic','chisqr')
    minreturn = 1.e-60
# THERE IS A MAJOR FLAW HERE
    if ~isinstance(fit_ranges[0],u.quantity.Quantity):
        fit_ranges*=u.micron
#    if ~isinstance(fit_ranges[0],list) or ~isinstance(fit_ranges[0],'numpy.ndarray'):
#        print(type(fit_ranges[0]))
#        fit_ranges = [fit_ranges]
#    print(fit_ranges)


# create interpolation function for second spectrum
    f = interp1d(sp2.wave,sp2.flux,bounds_error=False,fill_value=0.)
    if var_flag:
        v = interp1d(sp2.wave,sp2.variance*numpy.nan,bounds_error=False,fill_value=numpy.nan)
    else:
        v = interp1d(sp2.wave,sp2.variance,bounds_error=False,fill_value=numpy.nan)

# total variance - funny form to cover for nans
    vtot = numpy.nanmax([sp1.variance.value,sp1.variance.value+v(sp1.wave.value)],axis=0)
 #   vtot = sp1.variance

# Mask certain wavelengths
    mask = _generateMask(sp1.wave,**kwargs)
# mask flux < 0
    mask[numpy.where(numpy.logical_or(sp1.flux < 0,f(sp1.wave) < 0))] = 1

# set the weights
    for ranges in fit_ranges:
#        print(ranges)
        weights[numpy.where(numpy.logical_and(sp1.wave >= ranges[0],sp1.wave <= ranges[1]))] = 1

# combine weights and mask together to one array
    weights = weights*(1.-mask)

# comparison statistics
# switch to standard deviation if no uncertainty
    if numpy.isnan(numpy.nanmax(vtot)):
        statistic = 'stddev'
        if kwargs.get('verbose',False):
            print('No uncertainties provided; using the {} statistic by default'.format(statistic))
    else:
        if kwargs.get('verbose',False):
            print('Comparing spectra using the {} statistic'.format(statistic))

# chi^2
    if (statistic == 'chisqr' or statistic == 'chisq' or statistic == 'chi'):
# compute scale factor
        scale = numpy.nansum(weights*sp1.flux.value*f(sp1.wave)/vtot)/ \
            numpy.nansum(weights*f(sp1.wave)*f(sp1.wave)/vtot)
# correct variance
        vtot = numpy.nanmax([sp1.variance.value,sp1.variance.value+v(sp1.wave)*scale**2],axis=0)
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave)*scale)**2/vtot)
        unit = sp1.funit/sp1.funit

# normalized standard deviation
    elif (statistic == 'stddev_norm' or statistic == 'stdev_norm'):
# compute scale factor
        scale = numpy.nansum(weights*sp1.flux.value)/ \
            numpy.nansum(weights*f(sp1.wave))
# correct variance
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave)*scale)**2)/ \
            numpy.median(sp1.flux.value)**2
        unit = sp1.funit/sp1.funit

# standard deviation
    elif (statistic == 'stddev' or statistic == 'stdev'):
# compute scale factor
        scale = numpy.nansum(weights*sp1.flux.value*f(sp1.wave))/ \
            numpy.nansum(weights*f(sp1.wave)*f(sp1.wave))
# correct variance
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave)*scale)**2)
        unit = sp1.funit**2

# absolute deviation
    elif (statistic == 'absdev'):
# compute scale factor
        scale = numpy.nansum(weights*sp1.flux.value)/ \
            numpy.nansum(weights*f(sp1.wave))
# correct variance
        stat = numpy.nansum(weights*abs(sp1.flux.value-f(sp1.wave)*scale))
        unit = sp1.funit

# error
    else:
        print('Error: statistic {} for compareSpectra not available'.format(statistic))
        return numpy.nan, numpy.nan

# plot spectrum compared to best spectrum
    if (kwargs.get('plot',False) != False):
        spcomp = sp2.copy()
        spcomp.scale(scale)
        kwargs['colors'] = kwargs.get('colors',['k','r','b'])
        kwargs['title'] = kwargs.get('title',sp1.name+' vs '+sp2.name)
        from .plot import plotSpectrum
        plotSpectrum(sp1,spcomp,sp1-spcomp,**kwargs)

    return numpy.nanmax([stat,minreturn])*unit, scale





def _generateMask(wave,**kwargs):
    '''
    :Purpose: Generates a mask array based on wavelength vector and optional inputs on what to mask.

    :Output: A mask array, where 0 = OK and 1 = ignore

    :Example:
    '''

# generate initial mask structure
    mask = numpy.zeros(len(wave))
    mask_ranges = kwargs.get('mask_ranges',[])

# mask telluric bands
    if kwargs.get('mask_telluric',False):
        mask_ranges.append([0.,0.65])        # meant to clear out short wavelengths
        mask_ranges.append([1.35,1.42])
        mask_ranges.append([1.8,1.92])
        mask_ranges.append([2.45,99.])        # meant to clear out long wavelengths

# a standardized masking
    if kwargs.get('mask_standard',False):
        mask_ranges.append([0.,0.8])        # standard short cut
        mask_ranges.append([1.35,1.42])
        mask_ranges.append([1.8,1.92])
        mask_ranges.append([2.35,99.])        # standard long cut

# make sure mask_ranges array has same units as wave array
    if ~isinstance(wave[0],u.quantity.Quantity):
        mask_ranges*=wave.unit

# generate mask
    for ranges in mask_ranges:
        mask[numpy.where(numpy.logical_and(wave >= ranges[0],wave <= ranges[1]))]= 1
    return mask



def measureEW(sp, *args, **kwargs):
    '''
    :Purpose: Measures equivalent widths (EWs) of specified lines
    :param sp: Spectrum class object, which should contain wave, flux and noise array elements
    :param args: wavelength arrays. Needs at least two arrays to measure line and continuum regions.
    :type nsamples: optional, default = 100
    :param nonoise:
    :type nonoise: optional, default = False
    :param line:
    :type nonoise: optional, default = ''

    .. not too sure about how this one works; will come back later.
    '''

# presets
    nsamples = kwargs.get('nsamples',100)
    noiseFlag = kwargs.get('nonoise',False)

# predefined lines
    specline = kwargs.get('line','').replace(' ','').lower()
    if 'nai' in specline:
        if '2.2' in specline:
            wave_line = [2.2020, 2.2120]
            wave_cont = [2.1965, 2.2125, 2.2175]
    elif 'cai' in specline:
        if '2.2' in specline:
            wave_line = [2.2580, 2.2690]
            wave_cont = [2.2510, 2.2580, 2.2705, 2.2760]
    else:
        if len(args) < 2:
            print('measureEW needs at least two wavelength arrays to measure line and continuum regions')
            return numpy.nan, numpy.nan
        else:
            wave_line = args[0]
            wave_cont = args[1]


# create interpolation routines
    w = numpy.where(numpy.isnan(sp.flux) == False)
    f = interp1d(sp.wave.value[w],sp.flux.value[w],bounds_error=False,fill_value=0.)
    w = numpy.where(numpy.isnan(sp.noise) == False)

# note that units are stripped out
    if (numpy.size(w) != 0):
        n = interp1d(sp.wave.value[w],sp.noise.value[w],bounds_error=False,fill_value=numpy.nan)
        noiseFlag = False or noiseFlag
    else:
        n = interp1d(sp.wave.value[:],sp.noise.value[:],bounds_error=False,fill_value=numpy.nan)
        noiseFlag = True or noiseFlag

    wLine = (numpy.arange(0,nsamples+1.0)/nsamples)* \
            (numpy.nanmax(wave_line)-numpy.nanmin(wave_line))+numpy.nanmin(wave_line)
    fLine = f(wLine)
    nLine = n(wLine)
    fCont = f(wave_cont)
    nCont = n(wave_cont)

# first compute value
    pCont = numpy.poly1d(numpy.polyfit(wave_cont,fCont,1))
    fContFit = pCont(wLine)
    ew = trapz((numpy.ones(len(fLine))-(fLine/fContFit)), wLine)*1e4
#monte carlo
    if noiseFlag == False:
        ews=[]
        for i in numpy.arange(nsamples):
#generate simulated fluxes
#            fContVar = fCont+numpy.random.normal(0.,1.)*nCont
#            fLineVar = fLine+numpy.random.normal(0.,1.)*nLine
            fContVar = numpy.random.normal(fCont,nCont)
            fLineVar = numpy.random.normal(fLine,nLine)

#linear fit to continuum
            pCont = numpy.poly1d(numpy.polyfit(wave_cont,fContVar,1))
            fContFit = pCont(wLine)
            ews.append(trapz((numpy.ones(len(fLineVar))-(fLineVar/fContFit)), wLine)*1e4)

# some error checking
#            plt.plot(wLine,fContFit,color='r')
#            plt.plot(wLine,fLine,color='k')
#            plt.show()

# following line is more correct but having problem with output
#       return numpy.nanmean(ew)*u.angstrom, numpy.nanstd(ew)*u.angstrom
        return ew*u.angstrom, numpy.nanstd(ew)*u.angstrom
    else:
        return ew*u.angstrom, numpy.nan


def measureEWSet(sp,*args,**kwargs):
    '''
    :Purpose: Measures equivalent widths (EWs) of lines from specified sets. Returns dictionary of indices.
    :param sp: Spectrum class object, which should contain wave, flux and noise array elements
    :param set: string defining which EW measurement set you want to use; options include:

            - *rojas*: EW measures from `Rojas-Ayala et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...748...93R>`_;
              uses Na I 2.206/2.209 Ca I 2.26 micron lines.

    :type set: optional, default = 'rojas'

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1555+0954')[0]
    >>> print splat.measureEWSet(sp, set = 'rojas')
        {'Na I 2.206/2.209': (1.7484002652013144, 0.23332441577025356), 'Ca I 2.26': (1.3742491939667159, 0.24867705962337672), 'names': ['Na I 2.206/2.209', 'Ca I 2.26'], 'reference': 'EW measures from Rojas-Ayala et al. (2012)'}
    '''
    set = kwargs.get('set','rojas')

# determine combine method
    if ('rojas' in set.lower()):
        reference = 'EW measures from Rojas-Ayala et al. (2012)'
        names = ['Na I 2.206/2.209','Ca I 2.26']
        ews = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        ews[0],errs[0] = measureEW(sp,[2.2020, 2.2120],[2.1965, 2.2125, 2.2175],**kwargs)
        ews[1],errs[1] = measureEW(sp,[2.2580, 2.2690],[2.2510, 2.2580, 2.2705, 2.2760],**kwargs)
    else:
        print('{} is not one of the sets used for measureIndexSet'.format(set))
        return numpy.nan

# output dictionary of indices
    result = {names[i]: (ews[i],errs[i]) for i in numpy.arange(len(names))}
    result['reference'] = reference
    result['names'] = names
#    result['reference'] = reference
#    return inds,errs,names

    return result


def measureIndex(sp,*args,**kwargs):
    '''
    :Purpose: Measure an index on a spectrum based on defined methodology
                measure method can be mean, median, integrate
                index method can be ratio = 1/2, valley = 1-2/3, OTHERS
                output is index value and uncertainty
    .. will also come back to this one
    '''

# keyword parameters
    method = kwargs.get('method','ratio')
    sample = kwargs.get('sample','integrate')
    nsamples = kwargs.get('nsamples',1000)
    noiseFlag = kwargs.get('nonoise',False)

# create interpolation functions
    w = numpy.where(numpy.isnan(sp.flux) == False)
    f = interp1d(sp.wave.value[w],sp.flux.value[w],bounds_error=False,fill_value=0.)
    w = numpy.where(numpy.isnan(sp.noise) == False)
# note that units are stripped out
    if (numpy.size(w) != 0):
        s = interp1d(sp.wave.value[w],sp.noise.value[w],bounds_error=False,fill_value=numpy.nan)
        noiseFlag = False
    else:
        s = interp1d(sp.wave.value[:],sp.noise.value[:],bounds_error=False,fill_value=numpy.nan)
        noiseFlag = True

# error checking on number of arguments provided
    if (len(args) < 2):
        print('measureIndex needs at least two samples to function')
        return numpy.nan, numpy.nan
    elif (len(args) < 3 and (method == 'line' or method == 'allers' or method == 'inverse_line')):
        print(method+' requires at least 3 sample regions')
        return numpy.nan, numpy.nan

# define the sample vectors
    value = numpy.zeros(len(args))
    value_sim = numpy.zeros((len(args),nsamples))

# loop over all sampling regions
    for i,waveRng in enumerate(args):
        xNum = (numpy.arange(0,nsamples+1.0)/nsamples)* \
            (numpy.nanmax(waveRng)-numpy.nanmin(waveRng))+numpy.nanmin(waveRng)
        yNum = f(xNum)
        yNum_e = s(xNum)

# first compute the actual value
        if (sample == 'integrate'):
            value[i] = trapz(yNum,xNum)
        elif (sample == 'average'):
            value[i] = numpy.nanmean(yNum)
        elif (sample == 'median'):
            value[i] = numpy.median(yNum)
        elif (sample == 'maximum'):
            value[i] = numpy.nanmax(yNum)
        elif (sample == 'minimum'):
            value[i] = numpy.nanmin(yNum)
        else:
            value[i] = numpy.nanmean(yNum)

# now do MonteCarlo measurement of value and uncertainty
        for j in numpy.arange(0,nsamples):

# sample variance
            if (numpy.isnan(yNum_e[0]) == False):
                yVar = yNum+numpy.random.normal(0.,1.)*yNum_e
# NOTE: I'M NOT COMFORTABLE WITH ABOVE LINE - SEEMS TO BE TOO COARSE OF UNCERTAINTY
# BUT FOLLOWING LINES GIVE UNCERTAINTIES THAT ARE WAY TOO SMALL
#                yVar = numpy.random.normal(yNum,yNum_e)
#                yVar = yNum+numpy.random.normal(0.,1.,len(yNum))*yNum_e
            else:
                yVar = yNum

# choose function for measuring indices
            if (sample == 'integrate'):
                value_sim[i,j] = trapz(yVar,xNum)
            elif (sample == 'average'):
                value_sim[i,j] = numpy.nanmean(yVar)
            elif (sample == 'median'):
                value_sim[i,j] = numpy.median(yVar)
            elif (sample == 'maximum'):
                value_sim[i,j] = numpy.nanmax(yVar)
            elif (sample == 'minimum'):
                value_sim[i,j] = numpy.nanmin(yVar)
            else:
                value_sim[i,j] = numpy.nanmean(yVar)

# compute index based on defined method
# default is a simple ratio
    if (method == 'ratio'):
        val = value[0]/value[1]
        vals = value_sim[0,:]/value_sim[1,:]
    elif (method == 'line'):
        val = 0.5*(value[0]+value[1])/value[2]
        vals = 0.5*(value_sim[0,:]+value_sim[1,:])/value_sim[2,:]
    elif (method == 'inverse_line'):
        val = 2.*value[0]/(value[1]+value[2])
        vals = 2.*value_sim[0,:]/(value_sim[1,:]+value_sim[2,:])
    elif (method == 'change'):
        val = 2.*(value[0]-value[1])/(value[0]+value[1])
        vals = 2.*(value_sim[0,:]-value_sim[1,:])/(value_sim[0,:]+value_sim[1,:])
    elif (method == 'allers'):
        val = (((numpy.mean(args[0])-numpy.mean(args[1]))/(numpy.mean(args[2])-numpy.mean(args[1])))*value[2] \
            + ((numpy.mean(args[2])-numpy.mean(args[0]))/(numpy.mean(args[2])-numpy.mean(args[1])))*value[1]) \
            /value[0]
        vals = (((numpy.mean(args[0])-numpy.mean(args[1]))/(numpy.mean(args[2])-numpy.mean(args[1])))*value_sim[2,:] \
            + ((numpy.mean(args[2])-numpy.mean(args[0]))/(numpy.mean(args[2])-numpy.mean(args[1])))*value_sim[1,:]) \
            /value_sim[0,:]
    else:
        val = value[0]/value[1]
        vals = value_sim[0,:]/value_sim[1,:]


# output mean, standard deviation
    if (noiseFlag):
        return val, numpy.nan
    else:
        return val, numpy.nanstd(vals)


# wrapper function for measuring specific sets of indices

def measureIndexSet(sp,**kwargs):
    '''
    :Purpose: Measures indices of ``sp`` from specified sets. Returns dictionary of indices.
    :param sp: Spectrum class object, which should contain wave, flux and noise array elements
    :param set: string defining which indices set you want to use; options include:

            - *bardalez*: H2O-J, CH4-J, H2O-H, CH4-H, H2O-K, CH4-K, K-J, H-dip, K-slope, J-slope, J-curve, H-bump, H2O-Y from `Bardalez Gagliuffi et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...794..143B>`_
            - *burgasser*: H2O-J, CH4-J, H2O-H, CH4-H, H2O-K, CH4-K, K-J from `Burgasser et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...637.1067B>`_
            - *tokunaga*: K1, K2 from `Tokunaga & Kobayashi (1999) <http://adsabs.harvard.edu/abs/1999AJ....117.1010T>`_
            - *reid*: H2O-A, H2O-B from `Reid et al. (2001) <http://adsabs.harvard.edu/abs/2001AJ....121.1710R>`_
            - *geballe*: H2O-1.2, H2O-1.5, CH4-2.2 from `Geballe et al. (2002) <http://adsabs.harvard.edu/abs/2002ApJ...564..466G>`_
            - *allers*: H2O, FeH-z, VO-z, FeH-J, KI-J, H-cont from `Allers et al. (2007) <http://adsabs.harvard.edu/abs/2007ApJ...657..511A>`_, `Allers & Liu (2013) <http://adsabs.harvard.edu/abs/2013ApJ...772...79A>`_
            - *testi*: sHJ, sKJ, sH2O-J, sH2O-H1, sH2O-H2, sH2O-K from `Testi et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...552L.147T>`_
            - *slesnick*: H2O-1, H2O-2, FeH from `Slesnick et al. (2004) <http://adsabs.harvard.edu/abs/2004ApJ...610.1045S>`_
            - *mclean*: H2OD from `McLean et al. (2003) <http://adsabs.harvard.edu/abs/2003ApJ...596..561M>`_
            - *rojas*: H2O-K2 from `Rojas-Ayala et al.(2012) <http://adsabs.harvard.edu/abs/2012ApJ...748...93R>`_

    :type set: optional, default = 'burgasser'

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1555+0954')[0]
    >>> print splat.measureIndexSet(sp, set = 'reid')
        {'H2O-B': (1.0531856077273236, 0.0045092074790538221), 'H2O-A': (0.89673318593633422, 0.0031278302105038594)}
    '''
# keyword parameters
    set = kwargs.get('set','burgasser')

    if ('allers' in set.lower()):
        reference = 'Indices from Allers et al. (2007), Allers & Liu (2013)'
        refcode = '2013ApJ...772...79A'
        names = ['H2O','FeH-z','VO-z','FeH-J','KI-J','H-cont']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.55,1.56],[1.492,1.502],method='ratio',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[0.99135,1.00465],[0.97335,0.98665],[1.01535,1.02865],method='allers',sample='average',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[1.05095,1.06505],[1.02795,1.04205],[1.07995,1.09405],method='allers',sample='average',**kwargs)
        inds[3],errs[3] = measureIndex(sp,[1.19880,1.20120],[1.19080,1.19320],[1.20680,1.20920],method='allers',sample='average',**kwargs)
        inds[4],errs[4] = measureIndex(sp,[1.23570,1.25230],[1.21170,1.22830],[1.26170,1.27830],method='allers',sample='average',**kwargs)
        inds[5],errs[5] = measureIndex(sp,[1.54960,1.57040],[1.45960,1.48040],[1.65960,1.68040],method='allers',sample='average',**kwargs)
    elif ('bardalez' in set.lower()):
        reference = 'Indices from Bardalez Gagliuffi et al. (2014)'
        refcode = '2014ApJ...794..143B'
        names = ['H2O-J','CH4-J','H2O-H','CH4-H','H2O-K','CH4-K','K-J','H-dip','K-slope','J-slope','J-curve','H-bump','H2O-Y']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.14,1.165],[1.26,1.285],method='ratio',sample='integrate',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[1.315,1.335],[1.26,1.285],method='ratio',sample='integrate',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[1.48,1.52],[1.56,1.60],method='ratio',sample='integrate',**kwargs)
        inds[3],errs[3] = measureIndex(sp,[1.635,1.675],[1.56,1.60],method='ratio',sample='integrate',**kwargs)
        inds[4],errs[4] = measureIndex(sp,[1.975,1.995],[2.08,2.12],method='ratio',sample='integrate',**kwargs)
        inds[5],errs[5] = measureIndex(sp,[2.215,2.255],[2.08,2.12],method='ratio',sample='integrate',**kwargs)
        inds[6],errs[6] = measureIndex(sp,[2.06,2.10],[1.25,1.29],method='ratio',sample='integrate',**kwargs)
        inds[7],errs[7] = measureIndex(sp,[1.61,1.64],[1.56,1.59],[1.66,1.69] ,method='inverse_line',sample='integrate',**kwargs)
        inds[8],errs[8] = measureIndex(sp,[2.06,2.10],[2.10,2.14],method='ratio',sample='integrate',**kwargs)
        inds[9],errs[9] = measureIndex(sp,[1.27,1.30],[1.30,1.33],method='ratio',sample='integrate',**kwargs)
        inds[10],errs[10] = measureIndex(sp,[1.04,1.07],[1.26,1.29],[1.14,1.17],method='line',sample='integrate',**kwargs)
        inds[11],errs[11] = measureIndex(sp,[1.54,1.57],[1.66,1.69],method='ratio',sample='integrate',**kwargs)
        inds[12],errs[12] = measureIndex(sp,[1.04,1.07],[1.14,1.17],method='ratio',sample='integrate',**kwargs)
    elif ('burgasser' in set.lower()):
        reference = 'Indices from Burgasser et al. (2006)'
        refcode = '2006ApJ...637.1067B'
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
    elif ('geballe' in set.lower()):
        reference = 'Indices from Geballe et al. (2002)'
        refcode = '2002ApJ...564..466G'
        names = ['H2O-1.2','H2O-1.5','CH4-2.2']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.26,1.29],[1.13,1.16],method='ratio',sample='integrate',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[1.57,1.59],[1.46,1.48],method='ratio',sample='integrate',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[2.08,2.12],[2.215,2.255],method='ratio',sample='integrate',**kwargs)
    elif ('mclean' in set.lower()):
        reference = 'Indices from McLean et al. (2003)'
        refcode = '2003ApJ...596..561M'
        names = ['H2OD']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.951,1.977],[2.062,2.088],method='ratio',sample='average',**kwargs)
    elif ('reid' in set.lower()):
        reference = 'Indices from Reid et al. (2001)'
        refcode = '2001AJ....121.1710R'
        names = ['H2O-A','H2O-B']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.33,1.35],[1.28,1.30],method='ratio',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[1.47,1.49],[1.59,1.61],method='ratio',sample='average',**kwargs)
    elif ('rojas' in set.lower()):
        reference = 'Indices from Rojas-Ayala et al.(2012)'
        refcode = '2012ApJ...748...93R'
        names = ['H2O-K2']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        num, er1= measureIndex(sp,[2.070,2.090],[2.235,2.255],method='ratio',sample='average',**kwargs)
        den, er2= measureIndex(sp,[2.235,2.255],[2.360,2.380],method='ratio',sample='average',**kwargs)
        inds[0]= num/den
        errs[0]= inds[0]*numpy.sqrt((er1/num)**2+(er2/den)**2)
    elif ('slesnick' in set.lower()):
        reference = 'Indices from Slesnick et al. (2004)'
        refcode = '2004ApJ...610.1045S'
        names = ['H2O-1','H2O-2','FeH']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.335,1.345],[1.295,1.305],method='ratio',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[2.035,2.045],[2.145,2.155],method='ratio',sample='average',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[1.1935,1.2065],[1.2235,1.2365],method='ratio',sample='average',**kwargs)
    elif ('testi' in set.lower()):
        reference = 'Indices from Testi et al. (2001)'
        refcode = '2001ApJ...552L.147T'
        names = ['sHJ','sKJ','sH2O_J','sH2O_H1','sH2O_H2','sH2O_K']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.265,1.305],[1.6,1.7],method='change',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[1.265,1.305],[2.12,2.16],method='change',sample='average',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[1.265,1.305],[1.09,1.13],method='change',sample='average',**kwargs)
        inds[3],errs[3] = measureIndex(sp,[1.60,1.70],[1.45,1.48],method='change',sample='average',**kwargs)
        inds[4],errs[4] = measureIndex(sp,[1.60,1.70],[1.77,1.81],method='change',sample='average',**kwargs)
        inds[5],errs[5] = measureIndex(sp,[2.12,2.16],[1.96,1.99],method='change',sample='average',**kwargs)
    elif ('tokunaga' in set.lower()):
        reference = 'Indices from Tokunaga & Kobayashi (1999)'
        refcode = '1999AJ....117.1010T'
        names = ['K1','K2']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[2.1,2.18],[1.96,2.04],method='change',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[2.2,2.28],[2.1,2.18],method='change',sample='average',**kwargs)
    else:
        print('{} is not one of the sets used for measureIndexSet'.format(set))
        return numpy.nan

# output dictionary of indices
    result = {names[i]: (inds[i],errs[i]) for i in numpy.arange(len(names))}
#    result['reference'] = reference
#    return inds,errs,names

    return result


def metallicity(sp,**kwargs):
    '''
    :Purpose: Metallicity measurement using Na I and Ca I lines and H2O-K2 index as described in `Rojas-Ayala et al.(2012) <http://adsabs.harvard.edu/abs/2012ApJ...748...93R>`_
    :param sp: Spectrum class object, which should contain wave, flux and noise array elements
    :param nsamples: number of Monte Carlo samples for error computation
    :type nsamples: optional, default = 100

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='0559-1404')[0]
    >>> print splat.metallicity(sp)
        (-0.50726104530066363, 0.24844773591243882)
    '''
    nsamples = kwargs.get('nsamples',1000)

    coeff_feh = [-1.039,0.092,0.119]
    coeff_feh_e = [0.17,0.023,0.033]
    feh_unc = 0.100
    coeff_mh = [-0.731,0.066,0.083]
    coeff_mh_e = [0.12,0.016,0.023]
    mh_unc = 0.100

    h2ok2,h2ok2_e = measureIndexSet(sp, set='rojas')['H2O-K2']
    ew = measureEWSet(sp,set='rojas')
    nai = kwargs.get('nai',False)
    nai_e = kwargs.get('nai_e',0.)
    if nai is False:
        nai, nai_e = ew['Na I 2.206/2.209']
    cai = kwargs.get('cai',False)
    cai_e = kwargs.get('cai_e',0.)
    if cai is False:
        cai, cai_e = ew['Ca I 2.26']

    mh = coeff_mh[0]+(nai/h2ok2)*coeff_mh[1]+(cai/h2ok2)*coeff_mh[2]

# simulate uncertainties
    mhsim = numpy.ones(nsamples)*coeff_mh[0]+\
        (numpy.random.normal(nai,nai_e,nsamples)/numpy.random.normal(h2ok2,h2ok2_e,nsamples))*coeff_mh[1]+\
        (numpy.random.normal(cai,cai_e,nsamples)/numpy.random.normal(h2ok2,h2ok2_e,nsamples))*coeff_mh[2]

    return mh, numpy.sqrt(numpy.nanstd(mhsim)**2+mh_unc**2)




# run test program if calling from command line
if __name__ == '__main__':
    pass
#    test_gooddata()
#    splat.test()
#    test_info()
