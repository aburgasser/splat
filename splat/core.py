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
#    Joshua Hazlett
#    Carolina Herrara Hernandez
#    Elizabeth Moreno Hilario
#    Aishwarya Iyer
#    Yuhui Jin
#    Michael Lopez
#    Ryan Low
#    Dorsa Majidi
#    Diego Octavio Talavera Maya
#    Alex Mendez
#    Gretel Mercado
#    Niana Mohammed
#    Jonathan Parra
#    Maitrayee Sahi
#    Adrian Suarez
#    Melisa Tallis
#    Chris Theissen
#    Tomoki Tamiya
#    Steven Truong
#    Russell van Linge

# imports - internal
import copy
import bz2
import gzip
import os
import shutil
import sys
import warnings

# imports - external
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy
import pandas
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

# holding arrays
SPECTRA_READIN = {}
STDS_DWARF_SPEX = {}
STDS_SD_SPEX = {}
STDS_DSD_SPEX = {}
STDS_ESD_SPEX = {}
STDS_VLG_SPEX = {}
STDS_INTG_SPEX = {}

# databases - using the .txt files for now, will need to change to SQL at a future date
# these are now done using pandas
#DB_SOURCES = ascii.read(SPLAT_PATH+DB_FOLDER+DB_SOURCES_FILE)
#DB_SPECTRA = ascii.read(SPLAT_PATH+DB_FOLDER+DB_SPECTRA_FILE)
DB_SOURCES = pandas.read_csv(os.path.normpath(SPLAT_PATH+DB_FOLDER+DB_SOURCES_FILE),delimiter='\t')
DB_SPECTRA = pandas.read_csv(os.path.normpath(SPLAT_PATH+DB_FOLDER+DB_SPECTRA_FILE),delimiter='\t')

# suppress warnings - probably not an entirely safe approach!
numpy.seterr(all='ignore')
warnings.simplefilter('ignore')

# temporary constants - will be removed
max_snr = 1.e6                # maximum S/N ratio permitted

#######################################################
#######################################################
##################   DATA LOADING  ####################
#######################################################
#######################################################

DB_ALL_SPECTRA = pandas.DataFrame()
for k in (DB_SPECTRA_DEFAULT_PARAMETERS.keys()): DB_ALL_SPECTRA[k] = []

def _processNewData(verbose=True,**kwargs):
    pass

def _initializeAllData(override=False,verbose=False,**kwargs):
    '''
    :Purpose:

        Initializes the spectral data available for analysis by adding to splat.SPECTRAL_DATA global variable

    :Required Inputs:

        None

    :Optional Inputs:

        * :param verbose = False: provide verbose feedback

    :Outputs:
        
        None

    '''
##### ISSUES #####
# need to figure out how to properly read in quantities

# default information for a new data set    
    DATA_FOLDERS = []
# structure:
#   instrument name
#       folder

# folders from which data are to be found
    mfolders = [SPLAT_PATH+DATA_FOLDER+'../']

# as specified in .splat_spectral_data
    if os.path.exists(EXTERNAL_DATA_FILE):
        with open(EXTERNAL_DATA_FILE, 'r') as frd: x = frd.read()
        mfolders.extend(x.split('\n'))
    if os.path.exists(HOME_FOLDER+'/'+EXTERNAL_DATA_FILE):
        with open(HOME_FOLDER+'/'+EXTERNAL_DATA_FILE, 'r') as frd: x = frd.read()
        mfolders.extend(x.split('\n'))
# specified in environmental variable SPLAT_DATA
    if os.environ.get('SPLAT_DATA') != None:
        mfolders.extend(str(os.environ['SPLAT_DATA']).split(':'))

# check the folders
    if '' in mfolders: mfolders.remove('')
    rm = []
    for m in mfolders:
        if os.path.exists(m) == False: rm.append(m)
    if len(rm) > 0:
        for m in rm: mfolders.remove(m)
    if len(mfolders) == 0:
        if verbose == True: print('\nNo folders containing spectral data were found to be present')
        return
    mfolders = list(set(mfolders))

# go through each folder and check for subfolders = instrument names
    for i,f in enumerate(mfolders):
        instruments = os.listdir(f)
# only folders
        rm = []
        for m in instruments: 
            if os.path.isdir(os.path.join(f,m))==False: rm.append(m)
        if len(rm) > 0: 
            for m in rm: instruments.remove(m)
        if len(instruments) > 0:
            for inst in instruments:
                DATA_FOLDERS.append(os.path.join(f,inst))
                instcheck = checkInstrument(inst)
# new instrument; set defaults, then drew from kwargs or instrument definition file                
                if instcheck == False:
                    ikey = (inst.upper()).replace(' ','-').replace('_','-')
                    idict = {}
                    for k in list(INSTRUMENT_DEFAULT_PARAMETERS.keys()): idict[k] = kwargs.get(k,INSTRUMENT_DEFAULT_PARAMETERS[k]['default'])
                    if INSTRUMENT_DEFINITION_FILE in os.listdir(DATA_FOLDERS[-1]):
                        idict_read = readDictFromFile(os.path.join(DATA_FOLDERS[-1],INSTRUMENT_DEFINITION_FILE))
                        for k in list(idict_read.keys()): idict[k] = idict_read[k]
                    idict['altname'].append(inst)
                    INSTRUMENTS[ikey] = idict
                    if verbose == True:
                        print('Instrument {} was not one of the SPLAT defaults; added as a new instrument with following parameters:'.format(inst))
                        for k in list(idict.keys()): print('\t{} = {}'.format(k,idict[k]))
                else: 
                    inst = instcheck
                    idict = INSTRUMENTS[inst]
# save instrument information if not present
                if INSTRUMENT_DEFINITION_FILE not in os.listdir(DATA_FOLDERS[-1]) or override == True:
                    writeDictToFile(idict,os.path.join(DATA_FOLDERS[-1],INSTRUMENT_DEFINITION_FILE))

# read in or generate spectral_data.txt file, and integrate into a "super" database
        if len(DATA_FOLDERS) == 0:
            if verbose == True: print('No data uploaded')
            return

# STOPPED HERE
#         for i,inst in enumerate(list(SPLAT_DATA_INVENTORY.keys()))


# # now go through all data folders and add data
#                 fnm = os.path.join(f,nm)
#                 instruments = os.listdir(fnm)
#                 name = checkSpectralModelName(nm)
# # new model name, add to global variable
# # using info.txt data if available                    
#                 if name == False:
#                     name = nm
#                     adddict = {'name': name}
#                     definfo = copy.deepcopy(default_info)
#                     if 'info.txt' in instruments:
#                         with open(os.path.join(fnm,'info.txt'), 'r') as frd: x = frd.read()
#                         lines = x.split('\n')
#                         if '' in lines: lines.remove('')
#                         lines = [x.split('\t') for x in lines]
#                         adddict = dict(lines)
#                         if 'altnames' in list(adddict.keys()): adddict['altnames'] = adddict['altnames'].split(',')
#                    for k in list(SPECTRAL_MODELS[list(SPECTRAL_MODELS.keys())[0]].keys()):
#                        if k not in list(adddict.keys()):
#                            if k in list(default_info.keys()): adddict[k] = definfo[k]
#                            else: adddict[k] = ''
#                    for k in list(default_info.keys()):
#                        if k not in list(minfo.keys()): minfo[k] = default_info[k]

# #  this sets the default values - it would be better to just grab one file and set the defaults that way                    
#                     if 'default' not in list(adddict.keys()): adddict['default'] = {}
#                     for k in list(SPECTRAL_MODEL_PARAMETERS.keys()):
#                         if k in list(adddict.keys()): adddict['default'][k] = adddict[k]
#                         if 'default_'+k in list(adddict.keys()): adddict['default'][k] = adddict['default_'+k]
# #                        if k in list(adddict['default'].keys()): print(k,adddict['default'][k])
# #                    print('\nWarning: did not find info.txt file in {}; using default values for model information'.format(minfo['folder']))
# #                    adddict['name'] = nm
#                     if 'name' not in list(adddict.keys()): adddict['name'] = name
#                     if 'instruments' not in list(adddict.keys()): adddict['instruments'] = {}
#                     if 'bibcode' not in list(adddict.keys()): adddict['bibcode'] = ''
#                     SPECTRAL_MODELS[name] = adddict
#                     if verbose==True: print('\nAdded a new model {} with parameters {}'.format(name,adddict))
#                     del adddict, definfo
# # go through instruments                
#                 rm = []
#                 for m in instruments:
#                     if os.path.isdir(os.path.join(fnm,m))==False: rm.append(m)
#                 if len(rm) > 0:
#                     for m in rm: instruments.remove(m)
#                 if len(instruments) > 0:
#                     for inst in instruments:
# # make sure there are files in this folder
#                         fnmi = os.path.join(fnm,inst)
#                         mfiles = os.listdir(fnmi)
#                         if len(mfiles) > 0:
#                             instrument = checkInstrument(inst)
# # unknown instrument; just add for now
#                             if instrument == False:
#                                 instrument = (inst.replace(' ','-').replace('_','-')).upper()
#                             if instrument not in list(SPECTRAL_MODELS[name]['instruments'].keys()):
#                                 SPECTRAL_MODELS[name]['instruments'][instrument] = fnmi
#                                 if verbose == True: print('\nAdding model {} and instrument {} from {}'.format(name,instrument,fnmi))
#                             else:
#                                 if verbose == True: print('\nModel {} and instrument {}: ignoring {} as these already exists in {}'.format(name,instrument,fnmi,SPECTRAL_MODELS[name]['instruments'][instrument]))

#     return
    if verbose == True:
        print('Imported spectral data folders:')
        for d in DATA_FOLDERS: print('\t{}'.format(d))
    pass

_initializeAllData()


#####################################################
###############   Spectrum class   ##################
#####################################################

class Spectrum(object):
    '''
    :Description: 
        Class for containing spectral and source data from SpeX Prism Library.
        This is a temporary structure until astropy.specutils is completed

    Optional Inputs:

    :param ismodel: Set to True to specify spectrum as a model (default = False)
    :param wlabel: label of wavelength (default = 'Wavelength')
    :param wunit: unit in which wavelength is measured (default = u.micron)
    :param wunit_label: label of the unit of wavelength (default = 'micron')
    :param flabel: label of flux density (default = 'F\_lambda')
    :param fscale: string describing how flux density is scaled (default = '')
    :param funit: unit in which flux density is measured (default = u.erg/(u.cm**2 * u.s * u.micron)
    :param funit_label: label of the unit of flux density (default = 'erg cm\^-2 s\^-1 micron\^-1')
    :param resolution: Resolution of spectrum (default = median lam/lam step/2.)
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
       >>> sp = splat.Spectrum(10002)                           # read in spectrum with data_key = 10002
       >>> sp = splat.Spectrum(wave=wavearray,flux=fluxarray)   # create objects with wavelength & flux arrays
    '''

    def __init__(self, *args, **kwargs):
# some presets
        sdb = False
        self.ismodel = kwargs.get('ismodel',False)
        self.istransmission = kwargs.get('istransmission',False)
        self.wlabel = kwargs.get('wlabel',r'Wavelength')
        self.wunit = kwargs.get('wunit',DEFAULT_WAVE_UNIT)
        self.wunit_label = kwargs.get('wunit_label',self.wunit)
        self.flabel = kwargs.get('flabel',r'F$_{\lambda}$')
        self.fscale = kwargs.get('fscale','Arbitrary')
        if kwargs.get('surface',False) == True: self.fscale = 'Surface'
        if kwargs.get('apparent',False) == True: self.fscale = 'Apparent'
        if kwargs.get('absolute',False) == True: self.fscale = 'Absolute'
        self.funit = kwargs.get('funit',DEFAULT_FLUX_UNIT)
        self.funit_label = kwargs.get('funit_label',self.funit)
#        self.header = kwargs.get('header',fits.PrimaryHDU())
        self.header = kwargs.get('header',{})
        self.filename = kwargs.get('file','')
        self.filename = kwargs.get('filename',self.filename)
        self.name = kwargs.get('name','')
        self.idkey = kwargs.get('idkey',False)
# instrument
        self.instrument = kwargs.get('instrument','')
        inst = checkInstrument(self.instrument) 
        if inst != False: 
            self.instrument = inst
            for k in list(INSTRUMENTS[inst].keys()): setattr(self,k,kwargs.get(k,INSTRUMENTS[inst][k]))
        self.instrument_mode = kwargs.get('instrument_mode','')
#        self.runfast = kwargs.get('runfast',True)
        self.published = kwargs.get('published','N')
        self.bibcode = kwargs.get('bibcode','')
        self.history = []
        self.wave = []
        self.flux = []
        self.noise = []
        self.variance = []

# process arguments
# option 1: a filename is given
        if len(args) > 0:
            if isinstance(args[0],str):
                self.filename = args[0]

# option 2: a spectrum ID is given
            elif isinstance(args[0],int) or isinstance(args[0],float) or isinstance(args[0],numpy.int64) or isinstance(args[0],numpy.float64):
                self.idkey = int(args[0])
                try:
                    sdb = keySpectrum(self.idkey)
                    if isinstance(sdb,bool) == False:
                        self.filename = sdb['DATA_FILE'].iloc[0]
                except:
                    print('Warning: problem reading in spectral database')

# option 3: arrays are given - interpret as wave, flux, and optionally noise
        if len(args) > 1:
            if (isinstance(args[0],list) or isinstance(args[0],numpy.ndarray)) and (isinstance(args[1],list) or isinstance(args[1],numpy.ndarray)):
                kwargs['wave'] = kwargs.get('wave',args[0])
                kwargs['flux'] = kwargs.get('flux',args[1])
        if len(args) > 2:
            if isinstance(args[2],list) or isinstance(args[2],numpy.ndarray):
                kwargs['noise'] = kwargs.get('noise',args[2])

        if len(kwargs.get('wave','')) > 0 and len(kwargs.get('flux','')) > 0:
            self.wave = kwargs['wave']
            self.flux = kwargs['flux']
            if len(kwargs.get('noise','')) > 0:
                self.noise = kwargs['noise']
            else:
                self.noise = numpy.zeros(len(self.wave))
# some extras
            others = ['pixel','mask','flag','flags','model']
            for o in others:
                if len(kwargs.get(o,'')) > 0:
                    setattr(self,o,kwargs[o])

# read in file (NOTE: this overrules passing wave, flux, noise arrays)
        elif self.filename != '':

# set up parameters
            mkwargs = {}
            mkwargs['instrument']=self.instrument
            mkwargs['instrument_mode']=self.instrument_mode
            mkwargs['filename']=self.filename
            self.simplefilename = os.path.basename(self.filename)
#            self.file = self.filename
            self.name = kwargs.get('name',self.simplefilename)

# folder is by default the current directory
            mkwargs['folder'] = kwargs.get('folder','./')
            mkwargs['wunit'] = self.wunit
            mkwargs['funit'] = self.funit

# is this in the SPLAT database? if so use the default folder
# NOTE: NEED TO MAKE THIS INSTRUMENT FLEXIBLE
            if self.filename in list(DB_SPECTRA['DATA_FILE']): 
                mkwargs['folder'] = SPLAT_PATH+DATA_FOLDER                
#                mkwargs['data_key'] = int(DB_SPECTRA['DATA_KEY'][DB_SPECTRA['DATA_FILE']==self.filename])                
                mkwargs['silent']=True
                sdb = searchLibrary(**mkwargs)

# return prior spectrum - THIS IS NOT WORKING SO COMMENTED OUT
#            if self.filename in list(SPECTRA_READIN.keys()) and self.runfast == True:
#                self = SPECTRA_READIN[self.filename]
#                return

#            try:

        # breakouts for specific instruments
            if (kwargs.get('APOGEE') == True or kwargs.get('apogee') == True or kwargs.get('instrument','SPEX-PRISM').upper() == 'APOGEE') and self.filename != '':
                rs = _readAPOGEE(self.filename,**kwargs)
                self.instrument = 'APOGEE'
#                for k in list(rs.keys()): setattr(self,k.lower(),rs[k])
                self.history.append('Spectrum successfully loaded')
        # create a copy to store as the original
                self.original = copy.deepcopy(self)

            elif (kwargs.get('BOSS',False) == True or kwargs.get('boss',False) == True or kwargs.get('eboss',False) == True or kwargs.get('EBOSS',False) == True or kwargs.get('instrument','SPEX-PRISM').upper() == 'BOSS' or kwargs.get('instrument','SPEX-PRISM').upper() == 'EBOSS') and self.filename != '':
                rs = _readBOSS(self.filename)
#                for k in list(rs.keys()): setattr(self,k.lower(),rs[k])
                self.wunit = kwargs.get('wunit',u.Angstrom)
                self.funit = kwargs.get('funit',u.erg/(u.cm**2 * u.s * u.Angstrom))
                self.history.append('Spectrum successfully loaded')
        # create a copy to store as the original
                self.original = copy.deepcopy(self)


            else:
                rs = readSpectrum(self.filename,**mkwargs)
            for k in list(rs.keys()): setattr(self,k.lower(),rs[k])


# None of this worked; create an empty Spectrum object (used for copying)
        else:
#            print('Warning: Creating an empty Spectrum object')
            return

# process spectral data
        if len(self.wave) > 0:
# convert to numpy arrays
            self.wave = numpy.array(self.wave)
            self.flux = numpy.array(self.flux)
            self.noise = numpy.array(self.noise)
# enforce positivity and non-nan
# THESE HAVE BEEN REMOVED
#            if (numpy.nanmin(self.flux) < 0):
#                self.flux[numpy.where(self.flux < 0)] = 0.
#            self.flux[numpy.isnan(self.flux)] = 0.
            try:
                if numpy.nanmin(self.noise) < 0.:
                    self.noise[numpy.where(self.noise < 0.)] = 0.
            except:
                pass
# check on noise being too low - a problem with old SpeX files
# THIS HAS BEEN REMOVED
#            if (numpy.nanmax(self.flux/self.noise) > max_snr):
#                self.noise[numpy.where(self.flux/self.noise > max_snr)]=numpy.median(self.noise)
# convert to astropy quantities with units
            if not isUnit(self.wave):
                self.wave = numpy.array(self.wave)*self.wunit
            if not isUnit(self.flux):
                self.flux = numpy.array(self.flux)*self.funit
            if not isUnit(self.noise):
                self.noise = numpy.array(self.noise)*self.funit

# some conversions
            self.flam = self.flux
            try:
                self.nu = self.wave.to('Hz',equivalencies=u.spectral())
            except:
                pass
            try:
                self.fnu = self.flux.to('Jy',equivalencies=u.spectral_density(self.wave))
                self.fnu_unit = u.Jansky
            except:
                pass
            try:
                self.noisenu = self.noise.to('Jy',equivalencies=u.spectral_density(self.wave))
            except:
                pass
            self.temperature = numpy.zeros(len(self.flux))
# calculate variance & S/N
            self.variance = numpy.array([n**2 for n in self.noise.value])*self.noise.unit*self.noise.unit
            self.snr = self.computeSN()

# estimate resolution - be default central lam/lam spacing/3
            i = int(0.5*len(self.wave))
            self.resolution = kwargs.get('resolution',self.wave.value[i]/numpy.absolute(self.wave.value[i+1]-self.wave.value[i])/2.)

# populate information on source and spectrum from database
#        print(sdb)
        if isinstance(sdb,bool) == False :
            if isinstance(sdb,pandas.core.frame.DataFrame) and len(sdb) != 0:
                for k in list(sdb.columns):
                    setattr(self,k.lower(),str(sdb[k].iloc[0]))
#            elif isinstance(sdb,dict) == True: 
#                for k in list(sdb.keys()):
#                    setattr(self,k.lower(),str(sdb[k][0]))
            else:
                try:
                    for k in list(sdb.keys()):
                        setattr(self,k.lower(),str(sdb[k][0]))
                except:
                    pass
        try:
            self.shortname = designationToShortName(self.designation)
        except:
            pass
        try:
            self.date = str(self.observation_date)
        except:
            pass
# convert some data into numbers
        kconv = ['ra','dec','julian_date','median_snr','resolution','airmass',\
        'jmag','jmag_error','hmag','hmag_error','kmag','kmag_error','source_key']
        for k in kconv:
            try:
                setattr(self,k,float(getattr(self,k)))
            except:
                pass
# this is to make sure the database resolution is the default value
        try:
            if kwargs.get('resolution',False) == False:
                kwargs['resolution'] = self.resolution
            if kwargs.get('instrument',False) == False:
                kwargs['resolution'] = self.resolution
        except:
            pass

# instrument specific information
        hkys = list(self.header.keys())

# automated stuff for spex data
        if 'INSTRUME' in hkys:
            if 'spex' in self.header['INSTRUME'].lower() and 'GRAT' in hkys:
                if 'lowres15' in self.header['GRAT'].lower() or 'prism' in self.header['GRAT'].lower(): self.instrument = 'SPEX-PRISM'
                if 'shortxd' in self.header['GRAT'].lower() or 'sxd' in self.header['GRAT'].lower(): self.instrument = 'SPEX-SXD'
        if 'INSTR' in hkys:
            if 'spex' in self.header['INSTR'].lower() and 'GRAT' in hkys:
                if 'lowres15' in self.header['GRAT'].lower() or 'prism' in self.header['GRAT'].lower(): self.instrument = 'SPEX-PRISM'
                if 'shortxd' in self.header['GRAT'].lower() or 'sxd' in self.header['GRAT'].lower(): self.instrument = 'SPEX-SXD'
        if 'spex' in self.instrument.lower():
            if 'observation_date' in list(self.__dict__.keys()): dt = self.observation_date
            elif 'OBS-DATE' in hkys: dt = self.header['OBS-DATE'].replace('-','')
            elif 'OBS_DATE' in hkys: dt = self.header['OBS_DATE'].replace('-','')
            elif 'DATE_OBS' in hkys: dt = self.header['DATE_OBS'].replace('-','')
            elif 'DATE-OBS' in hkys: dt = self.header['DATE-OBS'].replace('-','')
            else: dt = '20000101'
#            if int(dt) > 20140800:
#                self.instrument = self.instrument.replace('SPEX','USPEX')

# populate defaults
        inst = checkInstrument(self.instrument)
        if inst == False: 
#            print(self.instrument)
            inst = 'UNKNOWN'
        self.instrument = inst
        for k in list(INSTRUMENTS[self.instrument].keys()): setattr(self,k,INSTRUMENTS[self.instrument][k])
        self.slitpixelwidth = (self.slitwidth/self.pixelscale).value
        if kwargs.get('resolution',False) != False:
            if kwargs['resolution'] > 0.:
                self.slitwidth = self.slitwidth*self.resolution/kwargs['resolution']
                self.slitpixelwidth = self.slitpixelwidth*self.resolution/kwargs['resolution']
                self.resolution = kwargs['resolution']
        if kwargs.get('slitwidth',False) != False:
            sl = kwargs['slitwidth']
            if isUnit(sl): sl = sl.value
            if s1 > 0.:
                self.resolution = self.resolution*(self.slitwidth.value)/sl
                self.slitpixelwidth = self.slitpixelwidth*s1/(self.slitwidth.value)
                self.slitwidth = kwargs['slitwidth']
                if not isUnit(self.slitwidth):
                    self.slitwidth = self.slitwidth*u.arcsec
        if kwargs.get('slitpixelwidth',False) != False:
            if kwargs['slitpixelwidth'] > 0.:
                self.resolution = self.resolution*self.slitpixelwidth/kwargs['slitpixelwidth']
                self.slitwidth = self.slitwidth*kwargs['slitpixelwidth']/self.slitpixelwidth
                self.slitpixelwidth = kwargs['slitpixelwidth']
#        else:
#            kys = list(splat.INSTRUMENTS.keys())
#            for k in list(INSTRUMENTS[kys[0]].keys()): setattr(self,k,'UNKNOWN')
#            self.instrument_name = self.instrument
#            self.slitpixelwidth = 3.
#            self.resolution = (self.wave[int(0.5*len(self.wave))]/(self.wave[int(0.5*len(self.wave))+2]-self.wave[int(0.5*len(self.wave))-1])).value
#            self.slitwidth = 1.0*u.arcsec # not sure I should set this
#            print('Warning: {} is not one of the instruments defined for SPLAT'.format(self.instrument))
#            print('Setting slit width to {}, slit pixel width to {} and resolution to {}'.format(self.slitwidth,self.slitpixelwidth,self.resolution))

        self.dof = numpy.round(len(self.wave)/self.slitpixelwidth)


# information on model
        if self.ismodel == True:
            self.teff = kwargs.get('teff',numpy.nan)
            self.logg = kwargs.get('logg',numpy.nan)
            self.z = kwargs.get('z',numpy.nan)
            self.fsed = kwargs.get('fsed',numpy.nan)
            self.cld = kwargs.get('cld',numpy.nan)
            self.kzz = kwargs.get('kzz',numpy.nan)
            self.slit = kwargs.get('slit',numpy.nan)
            self.model = kwargs.get('model','')
            mset = checkSpectralModelName(self.model)
#            print(self.model,mset)
            if mset != False:
                self.model = mset
                self.modelset = mset
                for k in list(SPECTRAL_MODELS[mset].keys()):
#                    print(k,SPECTRAL_MODELS[mset][k])
                    setattr(self,k.lower(),SPECTRAL_MODELS[mset][k])
#            self.name = self.model+' model'
#            self.shortname = self.name
            self.fscale = 'Surface'
            self.published = 'Y'
#            print(self.name)

# information on transmission/filter spectrum
        elif self.istransmission == True:
            self.shortname = self.name
            self.fscale = 'Normalized'
            self.published = 'Y'
# populate header            
        else:
            kconv = {'designation': 'DESIG','name': 'NAME','shortname': 'SNAME','ra': 'RA_DEC','dec': 'DEC_DEC','slitwidth': 'SLTW_ARC','source_key': 'SRC_KEY','data_key': 'DATA_KEY','observer': 'OBSERVER', 'bibcode': 'BIB_DATA','program_pi': 'PI','program_number': 'PROGRAM','airmass': 'AIRMASS','reduction_spextool_version': 'VERSION','reduction_person': 'RED_PERS','reduction_date': 'RED_DATE','observation_date': 'OBSDATE','julian_date': 'JDATE','median_snr': 'SNR','resolution': 'RES', 'instrument': 'INSTRUME','wunit': 'XUNITS','funit': 'YUNITS', 'wlabel': 'XTITLE', 'flabel': 'YTITLE', 'opt_type': 'SPT_OPT', 'lit_type': 'SPT_LIT', 'nir_type': 'SPT_NIR', 'spex_type': 'SPT_SPEX', 'gravity_class_nir': 'GRAV_NIR', 'gravity_class_opt': 'GRAV_OPT', 'metallicity_class': 'ZCLASS', 'luminosity_class': 'LUMCLASS', 'color_extremity': 'COLOREX','mu': 'MU','mu_e': 'E_MU', 'mu_ra': 'MU_RA', 'mu_dec': 'MU_DEC', 'parallax': 'PARALLAX', 'parallax_e': 'E_PARALL', 'vtan': 'VTAN','vtan_e': 'E_VTAN','rv': 'RV','rv_e': 'E_RV', 'vsini': 'VSINI', 'vsini_e': 'E_VSINI','distance': 'DISTANCE', 'distance_e': 'E_DISTAN','j_2mass': 'J_2MASS', 'h_2mass': 'H_2MASS', 'ks_2mass': 'K_2MASS', 'j_2mass_e': 'E_J_2MAS', 'h_2mass_e': 'E_H_2MAS', 'ks_2mass_e': 'E_K_2MAS', 'object_type': 'OBJ_TYPE', 'binary': 'BINARY','sbinary': 'SPBINARY', 'companion_name': 'COMPNAME', 'cluster': 'CLUSTER' }
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


        self.history.append('{} spectrum successfully loaded'.format(self.instrument))

# create a copy to store as the original
        self.original = copy.deepcopy(self)

# add to previous read spectra
        if self.filename != '' and self.ismodel == False:
            SPECTRA_READIN[self.filename] = self

        return


    def mapTo(self,other,overhang=0.1):
        '''
        Purpose: maps spectrum onto the wavelength scale of another spectrum

        '''
        self.toWaveUnit(other.wave.unit)
        funit = self.flux.unit
        trng = [numpy.nanmin(other.wave.value),numpy.nanmax(other.wave.value)]
        dt = numpy.abs(trng[1]-trng[0])
        trng = [trng[0]-overhang*dt,trng[1]+overhang*dt]*other.wave.unit
        self.trim(trng)
        self.flux = reMap(self.wave.value,self.flux.value,other.wave.value)*funit
        self.noise = reMap(self.wave.value,self.noise.value,other.wave.value)*funit
        self.wave = other.wave
        self.variance = [x**2 for x in self.noise.value]*funit**2
        self.history.append('Mapped onto wavelength grid of {}'.format(other))
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
        :Purpose: A simple representation of the Spectrum object
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
# convert wavelength and flux units
        other.toWaveUnit(self.wunit)
        other.toFluxUnit(self.funit)

# make a copy and fill in wavelength to be overlapping
        sp = copy.deepcopy(self)
        sp.wave = self.wave.value[numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))]
        sp.wave=sp.wave*self.wunit

# generate interpolated axes
        f1 = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        f2 = interp1d(other.wave.value,other.flux.value,bounds_error=False,fill_value=0.)
        n1 = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
        n2 = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)

# add
        sp.flux = (f1(sp.wave.value)+f2(sp.wave.value))*self.funit

# uncertainty
        sp.variance = (n1(sp.wave.value)+n2(sp.wave.value))*(self.funit**2)
        sp.noise = sp.variance**0.5
        sp.snr = sp.computeSN()

# update information
        sp.name = self.name+' + '+other.name
# remove these attributes
        ref = ['date','observer','airmass','designation','source_key','data_key']
        for r in ref:
            if r in sp.__dict__.keys():
                delattr(sp,r)
# combine these attributes
        ref = ['bibcode','date','observer','airmass','designation','source_key','data_key']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,[getattr(self,r),getattr(other,r)])
#        ref = ['source_key','data_key']
#        for r in ref:
#            setattr(sp,r,0)
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
# convert wavelength and flux units
        other.toWaveUnit(self.wave.unit)
        other.toFluxUnit(self.flux.unit)

# make a copy and fill in wavelength to be overlapping
        sp = copy.deepcopy(self)
        sp.wave = self.wave.value[numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))]
# this fudge is for astropy 1.*
        if not isUnit(sp.wave):
            sp.wave=sp.wave*self.wave.unit

# generate interpolated axes
        f1 = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        f2 = interp1d(other.wave.value,other.flux.value,bounds_error=False,fill_value=0.)
        n1 = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
        n2 = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)

# subtract
        sp.flux = (f1(sp.wave.value)-f2(sp.wave.value))*self.flux.unit

# uncertainty
        sp.variance = (n1(sp.wave.value)+n2(sp.wave.value))*(self.flux.unit**2)
        sp.noise = sp.variance**0.5
        sp.snr = sp.computeSN()

# update information
        sp.name = self.name+' - '+other.name
        ref = ['date','observer','airmass','designation']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,'{} and {}'.format(getattr(self,r),getattr(other,r)))
        ref = ['source_key','data_key']
        for r in ref:
            setattr(sp,r,0)
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
# convert wavelength units
        other.toWaveUnit(self.wave.unit)

# make a copy and fill in wavelength to be overlapping
        sp = copy.deepcopy(self)
        sp.wave = self.wave.value[numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))]
        sp.wave=sp.wave*self.wave.unit

# generate interpolated axes
        f1 = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        f2 = interp1d(other.wave.value,other.flux.value,bounds_error=False,fill_value=0.)
        n1 = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
        n2 = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)

# multiply
        sp.flux = numpy.multiply(numpy.array(f1(sp.wave.value)),numpy.array(f2(sp.wave.value)))*self.flux.unit*other.flux.unit

# uncertainty
        sp.variance = numpy.multiply(sp.flux**2,((numpy.divide(n1(sp.wave.value),f1(sp.wave.value))**2)+(numpy.divide(n2(sp.wave.value),f2(sp.wave.value))**2)))
        sp.variance=sp.variance*((self.flux.unit*other.flux.unit)**2)
        sp.noise = sp.variance**0.5
        sp.cleanNoise()
        sp.snr = sp.computeSN()

# update information
        sp.name = self.name+' x '+other.name
        sp.funit = self.flux.unit*other.flux.unit
        sp.funit_label = str(self.flux.unit*other.flux.unit)
        ref = ['date','observer','airmass','designation']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,'{} and {}'.format(getattr(self,r),getattr(other,r)))
        ref = ['source_key','data_key']
        for r in ref:
            setattr(sp,r,0)
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
# convert wavelength units
        other.toWaveUnit(self.wave.unit)

# make a copy and fill in wavelength to be overlapping
        sp = copy.deepcopy(self)
        sp.wave = self.wave.value[numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))]
        sp.wave=sp.wave*self.wave.unit

# generate interpolated axes
        f1 = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        f2 = interp1d(other.wave.value,other.flux.value,bounds_error=False,fill_value=0.)
        n1 = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
        n2 = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)

# divide
        sp.flux = numpy.divide(numpy.array(f1(sp.wave.value)),numpy.array(f2(sp.wave.value)))*(self.flux.unit/other.flux.unit)

# uncertainty
        sp.variance = numpy.multiply(sp.flux**2,((numpy.divide(n1(sp.wave.value),f1(sp.wave.value))**2)+(numpy.divide(n2(sp.wave.value),f2(sp.wave.value))**2)))
        sp.variance=sp.variance*((self.flux.unit/other.flux.unit)**2)
        sp.noise = sp.variance**0.5
        sp.snr = sp.computeSN()

# clean up infinities
        sp.flux = (numpy.where(numpy.absolute(sp.flux.value) == numpy.inf, numpy.nan, sp.flux.value))*self.funit/other.flux.unit
        sp.cleanNoise()

# update information
        sp.name = self.name+' / '+other.name
        sp.funit = self.flux.unit/other.flux.unit
        sp.funit_label = str(self.flux.unit/other.flux.unit)
        ref = ['date','observer','airmass','designation']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,'{} and {}'.format(getattr(self,r),getattr(other,r)))
        ref = ['source_key','data_key']
        for r in ref:
            setattr(sp,r,0)
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
# convert wavelength units
        other.toWaveUnit(self.wave.unit)

# make a copy and fill in wavelength to be overlapping
        sp = copy.deepcopy(self)
        sp.wave = self.wave.value[numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))]
        sp.wave=sp.wave*self.wave.unit

# generate interpolated axes
        f1 = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        f2 = interp1d(other.wave.value,other.flux.value,bounds_error=False,fill_value=0.)
        n1 = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
        n2 = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)

# divide
        sp.flux = numpy.divide(numpy.array(f1(sp.wave.value)),numpy.array(f2(sp.wave.value)))*self.flux.unit/other.flux.unit

# uncertainty
        sp.variance = numpy.multiply(sp.flux**2,((numpy.divide(n1(sp.wave.value),f1(sp.wave.value))**2)+(numpy.divide(n2(sp.wave.value),f2(sp.wave.value))**2)))
        sp.variance=sp.variance*((self.flux.unit/other.flux.unit)**2)
        sp.noise = sp.variance**0.5
        sp.snr = sp.computeSN()

# clean up infinities
        sp.flux = (numpy.where(numpy.absolute(sp.flux.value) == numpy.inf, numpy.nan, sp.flux.value))*self.flux.unit/other.flux.unit
        sp.cleanNoise()

# update information
        sp.name = self.name+' / '+other.name
        sp.funit = self.flux.unit/other.flux.unit
        sp.funit_label = str(self.flux.unit/other.flux.unit)
        ref = ['date','observer','airmass','designation']
        for r in ref:
            if r in self.__dict__.keys() and r in other.__dict__.keys():
                setattr(sp,r,'{} and {}'.format(getattr(self,r),getattr(other,r)))
        ref = ['source_key','data_key']
        for r in ref:
            setattr(sp,r,0)
        sp.history.append('Division of {} by {}'.format(self.name,other.name))
# reset original
        sp.original = copy.deepcopy(sp)
        return sp


    def cleanNoise(self,replace='median'):
        ns = numpy.array(self.noise.value)
        rep = 0.
        if replace=='median': rep = numpy.nanmedian(ns)
        ns[numpy.isnan(ns)==True] = rep
        ns[numpy.isinf(ns)==True] = rep
        var = ns**2
        self.noise = ns*self.noise.unit
        self.variance = (ns**2)*self.noise.unit*self.noise.unit
        return

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

    def addNoise(self,snr=0.):
        '''
        :Purpose: 
            Adds noise to a spectrum based on either the current uncertainties scaled by optional input S/N
        
        :Required Inputs: 
            None

        :Optional Inputs: 
            *snr*: Signal-to-noise ratio to use to scale uncertainty spectrum (default = 0. => use existing S/N)

        :Output: 
            None (spectrum object changed in place)

        :Example:
           >>> sp = splat.Spectrum(10001)
           >>> sp.addNoise(snr=10.)
           >>> sp.computeSN()
                9.5702358777108234
        '''
        if snr > 0.: 
            snr0 = self.computeSN()
            self.noise= [x*snr0/snr for x in self.noise.value]*self.noise.unit
            self.flux = self.flux+(numpy.random.normal(numpy.zeros(len(self.noise)),self.noise.value))*self.noise.unit
        return 

    def info(self):
        '''
        :Purpose: 
            Returns a summary of properties for the Spectrum object

        :Required Inputs: 
            None

        :Optional Inputs: 
            None

        :Output: 
            Text summary describing the Spectrum object

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.info()
                IRTF SpeX prism spectrum of SDSS J074201.41+205520.5
                Observed by nan on 2008 Jan 10
                Airmass = nan
                Source designation = J07420130+2055198
                Median S/N = 34.0
                SpeX Classification = T5.0
                Literature Classification = T5 from Burgasser, A. J. et al. (2006, ApJ, 637, 1067-1093)
                Spectrum key = 10931, Source key = 10526
                If you use these data, please cite:
                    Burgasser, A. J. et al. (2010, ApJ, 710, 1142-1169)
                    bibcode: 2010ApJ...710.1142B
                History:
                    SPEX_PRISM spectrum successfully loaded
        '''
        if self.ismodel == True:
            f = '\n{} for instrument {} with the following parmeters:'.format(self.modelset,self.instrument)
            for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
                if hasattr(self,ms): f+='\n\t{} = {} {}'.format(ms,getattr(self,ms),SPECTRAL_MODEL_PARAMETERS[ms]['unit'])
#            f+='\nSmoothed to slit width {} {}'.format(self.slit,SPECTRAL_MODEL_PARAMETERS['slit']['unit'])
            f+='\n\nIf you use this model, please cite {}'.format(spbib.shortRef(SPECTRAL_MODELS[self.modelset]['bibcode']))
            f+='\nbibcode = {}\n'.format(SPECTRAL_MODELS[self.modelset]['bibcode'])
        elif self.istransmission == True:
            f = '\n{} spectrum'.format(self.name)
#            f+='\nSmoothed to slit width {} {}'.format(self.slit,SPECTRAL_MODEL_PARAMETERS['slit']['unit'])
            f+='\n\nIf you use these data, please cite {}'.format(spbib.shortRef(self.bibcode))
            f+='\nbibcode = {}\n'.format(self.bibcode)
        else:
            f = '\n'
            if hasattr(self,'instrument_name'): f+='{} '.format(self.instrument_name)
            if hasattr(self,'name'): f+='spectrum of {}'.format(self.name)
            if hasattr(self,'observer') and hasattr(self,'date'): 
                if isinstance(self.observer,list):
                    for i in range(len(self.observer)):
                        f+='\nObserved by {} on {}'.format(self.observer[i],properDate(self.date[i],output='YYYY MMM DD'))
                else:
                    f+='\nObserved by {} on {}'.format(self.observer,properDate(self.date,output='YYYY MMM DD'))
            if hasattr(self,'airmass'): f+='\nAirmass = {}'.format(self.airmass)
            if hasattr(self,'designation'): f+='\nSource designation = {}'.format(self.designation)
            if hasattr(self,'median_snr'): f+='\nMedian S/N = {}'.format(self.median_snr)
            if hasattr(self,'spex_type'): f+='\nSpeX Classification = {}'.format(self.spex_type)
            if hasattr(self,'lit_type'): 
                if isinstance(self.lit_type,list):
                    for i in range(len(self.lit_type)):
                        f+='\nLiterature Classification = {} from {}'.format(self.lit_type[i],spbib.shortRef(self.lit_type_ref[i]))
                else:
                    f+='\nLiterature Classification = {} from {}'.format(self.lit_type,spbib.shortRef(self.lit_type_ref))
            if hasattr(self,'source_key') and hasattr(self,'data_key'): 
                if isinstance(self.source_key,list):
                    for i in range(len(self.source_key)):
                        f+='\nSpectrum key = {}, Source key = {}'.format(int(self.data_key[i]),int(self.source_key[i]))
                else:
                    f+='\nSpectrum key = {}, Source key = {}'.format(int(self.data_key),int(self.source_key))
            if self.published == 'Y':
                f+='\n\nIf you use these data, please cite:'
                if isinstance(self.data_reference,list):
                    for i in range(len(self.data_reference)):
                        f+='\n\t{}'.format(spbib.shortRef(self.data_reference[i]))
                        f+='\n\tbibcode: {}'.format(self.data_reference[i])
                else:
                    f+='\n\t{}'.format(spbib.shortRef(self.data_reference))
                    f+='\n\tbibcode: {}'.format(self.data_reference)
            else:
                f+='\n\nUNPUBLISHED DATA'

        f+='\n\nHistory:'
        for h in self.history:
            f+='\n\t{}'.format(h)
        print(f)
        return

    def export(self,filename='',clobber=True,csv=False,tab=True,delimiter='\t',save_header=True,save_noise=True,comment='#',*args,**kwargs):
        '''
        :Purpose: 
            Exports a Spectrum object to either a fits or ascii file, depending on file extension given.  
            If no filename is explicitly given, the Spectrum.filename attribute is used. 
            If the filename does not include the full path, the file is saved in the current directory.  
            Spectrum.export and `Spectrum.save()`_ function in the same manner.

        .. _`Spectrum.save()` : api.html#splat.core.Spectrum.save

        :Required Inputs: 
            None

        :Optional Inputs: 
            :param filename: String specifying the filename to save; filename can also be included as an argument; if not provided, Spectrum.filename is used; alternate keywords: `file`
            :param clobber: Set to True to overwrite file, or False to raise flag if file exists (default = True) 
            :param csv: Set to True to write a CSV (comma-delimited) file (default = False) 
            :param tab: Set to True to write a tab-delimited file (default = True) 
            :param delimiter: character or string to specify as delimiter between columns (default = '\t'); alternate keywords: `sep` 
            :param save_header: set to True to add header to ascii files (default = True) 
            :param save_noise: set to True to save the noise column (default = True) 
            :param comment: use to specify comment character (default = '#') 

        :Output: 
            An ascii or fits file with the data and header

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

# prep inputs
        if len(args) > 0:
            filename = args[0]
        filename = kwargs.get('file',filename)

        if filename == '' and 'filename' in list(self.__dict__.keys()):
            filename = self.filename

        if filename == '':
            print('\nWarning! no filename provided, data were not saved')
            return

# update source information
        self.filename = os.path.basename(filename)
        self.simplefilename = self.filename

# determine which type of file
        ftype = filename.split('.')[-1]

# fits file
        if (ftype == 'fit' or ftype == 'fits'):
#            try:
            data = numpy.vstack((self.wave.value,self.flux.value,self.noise.value))
            hdu = fits.PrimaryHDU(data)
            for k in list(self.header.keys()):
                if k.upper() not in ['HISTORY','COMMENT','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND'] and k.replace('#','') != '': # and k not in list(hdu.header.keys()):
# stupidy because of astropy's ridiculous unit issues
#                    if isUnit(self.header[k]):
#                        try:
 #                           hdu.header[k] = self.header[k].value
 #                       except:
 #                           hdu.header[k] = self.header[k].scale
 #                   else:
                    hdu.header[k] = str(self.header[k])
            for k in list(self.__dict__.keys()):
                if isinstance(self.__getattribute__(k),str) == True or (isinstance(self.__getattribute__(k),float) == True and numpy.isnan(self.__getattribute__(k)) == False) or isinstance(self.__getattribute__(k),int) == True or isinstance(self.__getattribute__(k),bool) == True:
                    hdu.header[k.upper()] = str(self.__getattribute__(k))
#            print(hdu.header)
            hdu.writeto(filename,clobber=clobber)
#            except:
#                raise NameError('Problem saving spectrum object to file {}'.format(filename))

# ascii file - by default tab delimited
        else:
            delimiter = kwargs.get('sep',delimiter)
            f = open(filename,'w')
            if save_header == True:
                for k in list(self.header.keys()):
                    if k.upper() not in ['HISTORY','COMMENT'] and k.replace('#','') != '':
                        f.write('{}{} = {}\n'.format(comment,k.upper(),self.header[k]))
                for k in list(self.__dict__.keys()):
                    if isinstance(self.__getattribute__(k),str) == True or (isinstance(self.__getattribute__(k),float) == True and numpy.isnan(self.__getattribute__(k)) == False) or isinstance(self.__getattribute__(k),int) == True or isinstance(self.__getattribute__(k),bool) == True:
                        f.write('{}{} = {}\n'.format(comment,k.upper(),self.__getattribute__(k)))
                f.write('{}WAVELENGTH{}FLUX{}UNCERTAINTY\n'.format(comment,delimiter,delimiter))
            if save_noise == True:
                for i in range(len(self.wave.value)): f.write('{}{}{}{}{}\n'.format(self.wave.value[i],delimiter,self.flux.value[i],delimiter,self.noise.value[i]))
            else:
                for i in range(len(self.wave.value)): f.write('{}{}{}\n'.format(self.wave.value[i],delimiter,self.flux.value[i]))
            f.close()

        self.history.append('Spectrum saved to {}'.format(filename))
        return


    def save(self,*args,**kwargs):
        '''
        :Purpose: Exports a Spectrum object to either a fits or ascii file, depending on file extension given.  If no filename is explicitly given, the Spectrum.filename attribute is used. If the filename does not include the full path, the file is saved in the current directory.  

        `Spectrum.export()`_ and `Spectrum.save` function in the same manner.

        .. _`Spectrum.export()` : api.html#splat.core.Spectrum.export
        '''
        self.export(*args,**kwargs)


    def flamToFnu(self):
        return self.toFnu()

    def toFnu(self):
        '''
        :Purpose: 
            Converts flux density F\_nu in units of Jy.  
            There is no change if the spectrum is already in F\_nu units.

        :Required Inputs:
            None
        
        :Optional Inputs:
            None
        
        :Outputs:
            None; Spectrum object is changed
        
        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.toFnu()
           >>> sp.flux.unit
            Unit("Jy")
        '''
        if self.fscale == 'Temperature' or self.flabel == r'${\lambda}F_{\lambda}$':
            self.reset()
        self.funit = u.Jy
        self.flabel = r'$F_{\nu}$'
        self.flux = self.flux.to(self.funit,equivalencies=u.spectral_density(self.wave))
        self.noise = self.noise.to(self.funit,equivalencies=u.spectral_density(self.wave))
        self.snr = self.computeSN()
        self.history.append('Converted to Fnu units of {}'.format(self.funit))
        return

    def fnuToFlam(self):
        return self.toFlam()

    def toFlam(self):
        '''
        :Purpose: 
            Converts flux density to F\_lambda in units of erg/s/cm\^2/Hz. 
            There is no change if the spectrum is already in F\_lambda units.
        
        :Required Inputs:
            None
        
        :Optional Inputs:
            None
        
        :Outputs:
            None; Spectrum object is changed
        
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
        if self.fscale == 'Temperature' or self.flabel == r'${\lambda}F_{\lambda}$':
            self.reset()
        self.funit = DEFAULT_FLUX_UNIT
        self.flabel = r'$F_{\lambda}$'
        self.flux = self.flux.to(self.funit,equivalencies=u.spectral_density(self.wave))
        self.noise = self.noise.to(self.funit,equivalencies=u.spectral_density(self.wave))
        self.variance = self.noise**2
        self.snr = self.computeSN()
        self.history.append('Converted to Flam units of {}'.format(self.funit))
        return

    def toSED(self):
        '''
        :Purpose: 
            Converts flux density in F\_lambda to lambda x F\_lambda with units of erg/s/cm\^2. 

        :Required Inputs:
            None
        
        :Optional Inputs:
            None
        
        :Outputs:
            None; Spectrum object is changed

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.toSED()
           >>> sp.flux.unit
            Unit("erg / (cm2 s)")
        '''
        if self.flabel == 'SED':
            return
# first convert to F_lambda
        self.toFlam()
# now convert to SED
        un = self.wave.unit*self.flux.unit
        self.flux = (self.wave*self.flux).to(DEFAULT_SED_UNIT)
        self.noise = (self.wave*self.noise).to(DEFAULT_SED_UNIT)
        self.variance = self.noise**2
        self.snr = self.computeSN()
        self.funit = DEFAULT_SED_UNIT
        self.flabel = r'${\lambda}F_{\lambda}$'
        self.history.append('Converted to SED units of {}'.format(self.funit))
        return

    def toAngstrom(self):
        '''
        :Purpose: 
            Converts wavelength to Angstrom
        
        :Required Inputs:
            None
        
        :Optional Inputs:
            None
        
        :Outputs:
            None; Spectrum object is changed

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
        :Purpose: 
            Converts wavelength to microns. 
        
        :Required Inputs:
            None
        
        :Optional Inputs:
            None
        
        :Outputs:
            None; Spectrum object is changed

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

    def toWaveUnit(self,wunit):
        '''
        :Purpose: 
            Converts wavelength to specified units. 
        
        :Required Inputs:
            None
        
        :Optional Inputs:
            None
        
        :Outputs:
            None; Spectrum object is changed

        :Example:
           >>> import splat
           >>> import astropy.units as u
           >>> sp = splat.Spectrum(file='somespectrum',wunit=u.Angstrom)
           >>> print(sp.wave.unit)
             Angstrom
           >>> sp.toWaveUnit(u.micron)
           >>> print(sp.wave.unit)
             micron
           >>> sp.toWaveUnit(u.s)
             Warning! failed to convert wavelength unit from micron to s; no change made
        '''
        try:
            self.wave = self.wave.to(wunit)
            self.wunit = wunit
            self.history.append('Converted wavelength to units of {}'.format(self.wunit))
        except:
            print('\nWarning! failed to convert wavelength unit from {} to {}; no change made'.format(self.wave.unit,wunit))            
        return


    def toFluxUnit(self,funit):
        '''
        :Purpose: 
            Converts flux and noise arrays to given flux units.
        
        :Required Inputs:
            None
        
        :Optional Inputs:
            None
        
        :Outputs:
            None; Spectrum object is changed

        :Example:
           >>> import splat
           >>> import astropy.units as u
           >>> sp = splat.Spectrum(file='somespectrum',wunit=u.Angstrom)
           >>> sp.flux.unit
            erg / (cm2 micron s)
           >>> sp.toFluxUnit(u.Watt/u.m/u.m)
           >>> sp.flux.unit
            W / m2
           >>> sp.toFluxUnit(u.erg/u.s)
           >>> sp.flux.unit
            Warning! failed to convert flux unit from W / m2 to erg / s; no change made
        '''
        try:
            self.flux = self.flux.to(funit,equivalencies=u.spectral_density(self.wave))
            self.noise = self.noise.to(funit,equivalencies=u.spectral_density(self.wave))
            self.variance = self.noise**2
            self.history.append('Converted to flux units of {}'.format(self.funit))
            self.snr = self.computeSN()
            self.funit = funit
        except:
            print('Warning! failed to convert flux unit from {} to {}; no change made'.format(self.flux.unit,funit))
        return


    def toWavelengths(self,wave):
        '''
        :Purpose: 
            Maps a spectrum onto a new wavelength grid via interpolation or integral resampling

        :Required Inputs:
            :param wave: wavelengths to map to

        :Optional Inputs:
            None
        
        :Outputs:
            None; Spectrum object is changed

        :Example:
           TBD
        '''
        
        if not isUnit(wave):
            wave = wave*self.wave.unit

        if numpy.nanmin(((self.wave).to(wave.unit)).value) > numpy.nanmin(wave.value) or \
            numpy.nanmax(((self.wave).to(wave.unit)).value) < numpy.nanmax(wave.value):
            print('\nWarning: input wavelength range {} to {} outside spectrum wavelength range {} to {}'.format(numpy.nanmin(wave.value),numpy.nanmax(wave.value),numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)))
        else:
            self.toWaveUnit(wave.unit)
            self.trim([numpy.nanmin(wave),numpy.nanmax(wave)])

# map onto wavelength grid; if spectrum has lower resolution, interpolate; otherwise integrate & resample
            funit = self.flux.unit
            if len(self.wave) <= len(wave):
                f = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
                n = interp1d(self.wave.value,self.noise.value,bounds_error=False,fill_value=0.)
                self.flux = f(wave.value)*funit
                self.noise = n(wave.value)*funit
            else:
                self.flux = integralResample(self.wave.value,self.flux.value,wave.value)*funit
                self.noise = integralResample(self.wave.value,self.noise.value,wave.value)*funit
            self.wave = wave
            self.variance = self.noise**2
            self.history.append('Resampled to new wavelength grid')
            self.snr = self.computeSN()

        return


    def toInstrument(self,instrument,pixel_resolution=3.,**kwargs):
        '''
        :Purpose: 
            Converts a spectrum to the parameters for a defined instrument

        :Required Inputs:
            :param instrument: the name of an instrument, must be defined in splat.INSTRUMENTS

        :Optional Inputs:
            None
        
        :Outputs:
            None; Spectrum object is changed

        :Example:
           TBD
        '''
        required_parameters = ['wave_range','resolution']
        method = kwargs.get('method','hamming')
        oversample = kwargs.get('oversample',5)
        overscan = kwargs.get('overscan',0.05)
        wunit = self.wave.unit
        funit = self.flux.unit

# set up trim and smoothing parameters
        instr = checkInstrument(instrument)
        if instr == False:
            for r in required_parameters:
                if kwargs.get(r,False) == False:
                    print('\nInstrument {} is not defined in SPLAT; no change made'.format(instrument))
            wave_range = kwargs['wave_range']
            resolution = kwargs['resolution']
            instrument = instrument.upper()
            instrument = instrument.replace(' ','_')
            self.resolution = resolution
            self.wave_range = wave_range
        else:
            instrument = instr
            wave_range = INSTRUMENTS[instr]['wave_range']
            resolution = INSTRUMENTS[instr]['resolution']
            for k in list(INSTRUMENTS[instr].keys()): setattr(self,k,INSTRUMENTS[instr][k])
        if not isUnit(wave_range):
            wave_range = wave_range*wunit
        wave_range.to(self.wave.unit)
        wave_range = wave_range.value
# limit to range of data
        wave_range[0] = numpy.nanmax([wave_range[0],numpy.nanmin(self.wave.value)])        
        wave_range[1] = numpy.nanmin([wave_range[1],numpy.nanmax(self.wave.value)])        

# generate output wave vector
        effres = resolution*pixel_resolution
        npix = numpy.floor(numpy.log(numpy.nanmax(wave_range)/numpy.nanmin(wave_range))/numpy.log(1.+1./effres))
        wave_out = numpy.array([numpy.nanmin(wave_range)*(1.+1./effres)**i for i in numpy.arange(npix)])
        wave_out = wave_out[wave_out>numpy.nanmin(self.wave.value)]
        wave_out = wave_out[wave_out<numpy.nanmax(self.wave.value)]

# generate smoothing wavelength vector
        a = numpy.linspace(0.,len(wave_out)-1,len(wave_out))
        b = numpy.linspace(0.,len(wave_out)-1,int(oversample*len(wave_out)))
        f = interp1d(a,wave_out)
        wave_oversample = f(b)

# trim relevant piece of spectrum 
        dw = overscan*(numpy.nanmax(wave_range)-numpy.nanmin(wave_range))
        wrng = [numpy.nanmax([numpy.nanmin(wave_range)-dw,numpy.nanmin(self.wave.value)])*self.wave.unit,\
                numpy.nanmin([numpy.nanmax(wave_range)+dw,numpy.nanmax(self.wave.value)])*self.wave.unit]
        self.trim(wrng)

# map onto oversampled grid and smooth; if model is lower resolution, interpolate; otherwise integrate & resample
        if len(self.wave) <= len(wave_oversample):
            f = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
            flux_oversample = f(wave_oversample)
            if numpy.isnan(self.variance.value[0]) == False:
                v = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
                var_oversample = v(wave_oversample)
            else:
                var_oversample = [numpy.nan for f in flux_oversample]
        else:
            flux_oversample = integralResample(self.wave.value,self.flux.value,wave_oversample)
            if numpy.isnan(self.variance.value[0]) == False:
                var_oversample = integralResample(self.wave.value,self.variance.value,wave_oversample)
            else:
                var_oversample = [numpy.nan for f in flux_oversample]
        self.wave = wave_oversample*wunit
        self.flux = flux_oversample*funit
        self.variance = var_oversample*(funit**2)
        self.noise = numpy.array([x**0.5 for x in self.variance.value])*funit

# smooth this in pixel space including oversample       
        self._smoothToSlitPixelWidth(pixel_resolution*oversample,method=method)

# resample down to final wavelength scale
        flux_smooth = integralResample(self.wave.value,self.flux.value,wave_out)
        if numpy.isnan(self.variance.value[0]) == False:
            var_smooth = integralResample(self.wave.value,self.variance.value,wave_out)
        else:
            var_smooth = [numpy.nan for f in flux_smooth]
        self.wave = wave_out*wunit
        self.flux = flux_smooth*funit
        self.variance = var_smooth*(funit**2)
        self.noise = numpy.array([x**0.5 for x in self.variance.value])*funit

# do final trim
        self.trim(wave_range)

        self.history.append('Converted to instrument {}: wave range = {} resolution = {}'.format(instrument,wave_range,resolution))
        self.snr = self.computeSN()
        return

    def rvShift(self,rv):
        '''
        :Purpose: Shifts the wavelength scale by a given radial velocity. This routine changes the underlying Spectrum object.
        
        :Example:
           >>> import splat
           >>> import astropy.units as u
           >>> sp.rvShift(15*u.km/u.s)
        '''
        if not isUnit(rv):
            rv=rv*(u.km/u.s)
        rv.to(u.km/u.s)
        self.wave = self.wave*(1.+(rv/const.c).to(u.m/u.m))
        self.history.append('Shifted spectrum by radial velocity {}'.format(rv))

        return

    def broaden(self,vbroad,kern=None,epsilon=0.6,method='rotation',verbose=False):
        '''
        :Purpose: 

            Broadens a spectrum in velocity space using a line spread function (LSF) either based on rotation or gaussian. 
            This routine changes the underlying Spectrum object.

        :Required Inputs:

            :param vbroad: broadening width, nominally in km/s
            
        :Optional Inputs:

            :param method: method of broadening, should be one of:

                - ``gaussian``: (default) Gaussian broadening, with vbroad equal to Gaussian sigma
                - ``rotation``: rotational broadening using splat.lsfRotation()
                - ``delta``: Delta function (no broadening)

            :param kern: input kernel, must be at least three elements wide (default = None)
            :param epsilon: epsilon parameter for limb darkening in rotational broadening (default = 0.6)
            :param verbose: provide extra feedback (default = False)

        :Outputs:

            None; Spectral flux is smoothed using the desired line spread function. No change is made to noise or other axes
            
        :Example:
           >>> import splat
           >>> sp = splat.Spectrum(10001)
           >>> sp.broaden(30.,method='rotation')
           >>> sp.info()
            History:
                SPEX_PRISM spectrum successfully loaded
                Rotationally broadened spectrum by 30.0 km/s
        '''
        report = ''
        if isUnit(vbroad):
            vbroad.to(u.km/u.s)
        else:
            vbroad=vbroad*(u.km/u.s)

# determine velocity sampling
        samp = numpy.nanmedian(numpy.absolute(self.wave.value-numpy.roll(self.wave.value,1)) / self.wave.value)
#        samp = numpy.abs((self.wave.value[numpy.floor(0.5*len(self.wave))+1]-self.wave.value[numpy.floor(0.5*len(self.wave))])/self.wave.value[numpy.floor(0.5*len(self.wave))])
        vsamp = (samp*const.c).to(u.km/u.s)

# velocity resolution is too low - use a delta function
        if kern != None:
            if len(kern) < 3:
                if verbose: print('\nWarning: input kernel {} must be at least three elements; setting to delta function'.format(kern))
                kern = None
                method = 'delta'

        if kern == None:
            if vsamp > vbroad:
                if verbose: print('\nWarning: velocity resolution {} is smaller than velocity broadening {}; setting to delta function'.format(vsamp,vbroad))
                method = 'delta'

# rotational broadening
            if 'rot' in method.lower():
                kern = lsfRotation(vbroad.value,vsamp.value,epsilon=epsilon)
                report = 'Rotationally broadened spectrum by {}'.format(vbroad)

# gaussian ±10 sigma
            elif 'gauss' in method.lower():
                n = numpy.ceil(20.*vbroad.value/vsamp.value)
                if n%2==0: n+=1
                x = numpy.arange(n)-0.5*(n-1.)
                kern = numpy.exp(-0.5*(x**2))
                report = 'Broadened spectrum using a Gaussian with velocity width {}'.format(vbroad)

# delta function (no smoothing)
            else:
                kern = numpy.zeros(5)
                kern[2] = 1.
                report = 'Applying delta line spread function (no broadening)'

        else:
                report = 'Broadened spectrum using a input line spread function'

# normalize kernel
        kern = kern/numpy.nansum(kern)

        funit = self.flux.unit
        a = (numpy.nanmax(self.wave.value)/numpy.nanmin(self.wave.value))**(1./len(self.wave))
        nwave = numpy.nanmin(self.wave.value)*(a**numpy.arange(len(self.wave)))
        nflux = self.flux.value*nwave
        ncflux = numpy.convolve(nflux, kern, 'same')
        self.flux = ncflux/nwave*funit
        self.history.append(report)

        return


    def rotate(self,vsini,epsilon=0.6,verbose=False):
        '''
        :Purpose:

            Rotationally broaden the lines of a spectrum; a shortcut call to `Spectrum.broaden()`_

        .. _`Spectrum.broaden()` : api.html#splat.core.Spectrum.broaden

        :Required Inputs:
            :param vsini: Rotational velocity in km/s

        :Optional Inputs:
            :param epsilon: epsilon parameter for limb darkening in rotational broadening (default = 0.6)
            :param verbose: provide extra feedback (default = False)

        :Outputs:

            None; Spectral flux is smoothed by rotational broadening
            
        :Example:
           >>> import splat
           >>> sp = splat.Spectrum(10001)
           >>> sp.vsini(30.)
           >>> sp.info()
            History:
                SPEX_PRISM spectrum successfully loaded
                Rotationally broadened spectrum by 30.0 km/s
        '''          
        self.broaden(vsini,method='rotation',epsilon=epsilon,verbose=verbose)
        return


    def fluxCalibrate(self,filt,mag,**kwargs):
        '''
        :Purpose: Flux calibrates a spectrum given a filter and a magnitude. The filter must be one of those listed in `splat.FILTERS.keys()`. It is possible to specifically set the magnitude to be absolute (by default it is apparent).  This function changes the Spectrum object's flux, noise and variance arrays.
        
        Required Inputs:

        :param filt: string specifiying the name of the filter
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

# check inputs
        if isinstance(mag,str) and isNumber(filt):
            m = mag
            mag = filt
            filt = m
        if not isNumber(mag) or not isinstance(filt,str):
            raise ValueError('\nSyntax for function is Spectrum.filterMag(filter,magnitude)')

        if self.fscale == 'Temperature' or self.fscale == 'SED':
            self.reset()
        if self.funit != DEFAULT_FLUX_UNIT:
            self.toFlam()
        absolute = kwargs.get('absolute',False)
        apparent = kwargs.get('apparent',not absolute)
        apmag,apmag_e = filterMag(self,filt,**kwargs)
# NOTE: NEED TO INCORPORATE UNCERTAINTY INTO SPECTRAL UNCERTAINTY
        if (~numpy.isnan(apmag)):
            self.scale(10.**(0.4*(apmag-mag)))
            if (absolute):
                self.fscale = 'Absolute'
                self.history.append('Flux calibrated with {} filter to an absolute magnitude of {}'.format(filt,mag))
            if (apparent):
                self.fscale = 'Apparent'
                self.history.append('Flux calibrated with {} filter to an apparent magnitude of {}'.format(filt,mag))
        self.snr = self.computeSN()

        return


    def filterMag(self,filt,**kwargs):
        '''
        :Purpose: 

            Wrapper for `filterMag()`_ function in splat.photometry

        .. _`filterMag()` : api.html#splat.photometry.filterMag
        
        Required Inputs:

            **filter**: string specifiying the name of the filter

        Optional Inputs:

            See `filterMag()`_

        Outputs:

            Returns tuple containing filter-based spectrophotometic magnitude and its uncertainty

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.fluxCalibrate('2MASS J',15.0)
           >>> sp.filterMag(sp,'2MASS J')
            (15.002545668628173, 0.017635234089677564)
        '''

        from .photometry import filterMag
        return filterMag(self,filt,**kwargs)


# determine maximum flux, by default in non telluric regions
    def fluxMax(self,maskTelluric=True,**kwargs):
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
        if len(self.flux) == 0:
            print('\nWarning: spectrum object has a flux vector of zero length - maybe empty?')
            return numpy.nan

        if maskTelluric == True:            
            try:
                fl = self.flux.value[numpy.where(numpy.logical_or(\
                    numpy.logical_and(self.wave.to(u.micron) > 0.9*u.micron,self.wave.to(u.micron) < 1.35*u.micron),
                    numpy.logical_and(self.wave.to(u.micron) > 1.42*u.micron,self.wave.to(u.micron) < 1.8*u.micron),
                    numpy.logical_and(self.wave.to(u.micron) > 1.92*u.micron,self.wave.to(u.micron) < 2.3*u.micron)))]
                if len(fl) > 0: return numpy.nanmax(fl)*self.funit
#                if isUnit(fl[0]): fl = [f.value for f in fl]
            except:
                pass
        
#        fl = self.flux.value[numpy.where(\
#                numpy.logical_and(self.wave > numpy.nanmin(self.wave.value)+0.1*(numpy.nanmax(self.wave)-numpy.nanmin(self.wave)),self.wave < numpy.nanmax(self.wave)-0.1*(numpy.nanmax(self.wave)-numpy.nanmin(self.wave))))]
#        if isUnit(fl[0]): fl = [f.value for f in fl]
        return numpy.nanmax(self.flux.value)*self.flux.unit


    def normalize(self,*args,**kwargs):
        '''
        :Purpose: 

            Normalize a spectrum to a maximum value of 1 (in its current units) either at a 
            particular wavelength or over a wavelength range

        :Required Inputs: 

            None

        :Optional Inputs: 

            :param wave_range: choose the wavelength range to normalize; can be a list specifying minimum and maximum or a single wavelength (default = None); alternate keywords: `wave_range`, `range`

        :Output: 

            None; spectrum is normalized

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
        rng = kwargs.get('wave_range',False)
        rng = kwargs.get('waverange',rng)
        rng = kwargs.get('range',rng)
        if len(args) > 0:
            rng = args[0]
        if rng is not False:
            if not isinstance(rng,list) and not isinstance(rng,numpy.ndarray):
                rng = [rng]
            if isUnit(rng[0]): rng = [r.to(self.wave.unit).value for r in rng]
            if isUnit(rng): rng = rng.to(self.wave.unit).value
            if numpy.nanmax(rng) > numpy.nanmax(self.wave.value) or numpy.nanmin(rng) < numpy.nanmin(self.wave.value):
                print('\nWarning: normalization range {} is outside range of spectrum wave array: {}'.format(rng,[numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)]))
            if len(rng) == 1:
                f = interp1d(self.wave.value,self.flux.value)
                scalefactor = f(rng[0])
            else:
                scalefactor = numpy.nanmax(self.flux.value[numpy.where(numpy.logical_and(self.wave.value > rng[0],self.wave.value < rng[1]))])
        else:
            scalefactor = self.fluxMax(**kwargs)
        if isUnit(scalefactor): scalefactor = scalefactor.value
        if scalefactor == 0.: print('\nWarning: normalize is attempting to divide by zero; ignoring')
        elif numpy.isnan(scalefactor) == True: print('\nWarning: normalize is attempting to divide by nan; ignoring')
        else: 
            self.scale(1./scalefactor)
            self.fscale = 'Normalized'
            self.history.append('Spectrum normalized')
            self.snr = self.computeSN()
        return

    def plot(self,**kwargs):
        '''
        :Purpose: 

            calls the `plotSpectrum()`_ function, by default showing the noise spectrum and zeropoints. 
            See the `plotSpectrum()`_ API listing for details.

        .. _`plotSpectrum()`: api.html#splat.plot.plotSpectrum

        :Output: A plot of the Spectrum object

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.plot()
        '''
        kwargs['legend'] = kwargs.get('legend',self.name)
        kwargs['showNoise'] = kwargs.get('showNoise',True)
        kwargs['showZero'] = kwargs.get('showZero',True)
#        kwargs['xrange'] = kwargs.get('xrange',[numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)])
        from .plot import plotSpectrum
        return plotSpectrum(self,**kwargs)


    def redden(self, av=0.0, rv=3.1, normalize=False, a=10., n=1.33, **kwargs):
        '''
        :Purpose:

            Redden a spectrum based on an either Mie theory or a standard interstellar profile
            using Cardelli, Clayton, and Mathis (1989 ApJ. 345, 245)

        :Required Inputs:

            None

        :Optional Inputs:

            :param av: Magnitude of reddening A_V (default = 0.)
            :param rv: Normalized extinction parameter, R_V = A(V)/E(B-V) (default = 3.1
            :param normalized: Set to True to normalize reddening function (default = False)

        :Outputs:

            None; spectral flux is changed

        :Example:

           >>> import splat
           >>> sp = splat.Spectrum(10001)                   # read in a source
           >>> spr = splat.redden(sp,av=5.,rv=3.2)          # redden to equivalent of AV=5

        **Note**
          This routine is still in beta form; only the CCM89 currently works

        '''
        w = self.wave.to(DEFAULT_WAVE_UNIT).value                           # assuming in microns!

        if kwargs.get('mie',False) == True:                 # NOT CURRENTLY FUNCTIONING
            x = 2*numpy.pi*a/w
            x0 = 2.*numpy.pi*a/0.55                 # for V-band
            qabs = -4.*x*((n**2-1)/(n**2+2)).imag
            qsca = (8./3.)*(x**4)*(((n**2-1)/(n**2+2))**2).real
    #        tau = numpy.pi*(a**2)*(qabs+qsca)
            tau = 1.5*(qabs+qsca)/a    # for constant mass
            qabs0 = -4.*x0*((n**2-1)/(n**2+2)).imag
            qsca0 = (8./3.)*(x0**4)*(((n**2-1)/(n**2+2))**2).real
    #        tau0 = numpy.pi*(a**2)*(qabs0+qsca0)
            tau0 = 1.5*(qabs0+qsca0)/a    # for constant mass
            scale = (10.**(-0.4*av))
            absfrac = scale*numpy.exp(numpy.max(tau)-tau)
            report = 'Reddened by Mie scattering using grain size {} and index of refraction {}'.format(a,n)
        else:
            x = 1./w
            a = 0.574*(x**1.61)
            b = -0.527*(x**1.61)
            absfrac = 10.**(-0.4*av*(a+b/rv))
            report = 'Reddened following Cardelli, Clayton, and Mathis (1989) using A_V = {} and R_V = {}'.format(av,rv)

        if normalize == True:
            absfrac = absfrac/numpy.median(absfrac)
            report = report+' and normalized'

        self.flux = numpy.array(self.flux.value)*numpy.array(absfrac)*self.flux.unit
        self.noise = numpy.array(self.noise.value)*numpy.array(absfrac)*self.noise.unit
        self.variance = self.noise**2
        self.history.append(report)

        return


    def remove(self,mask,others=[]):
        '''
        :Purpose: 

            Removes elements of wave, flux, noise specified by the True/False array

        :Required Input:

            :param mask: An array of booleans, ints, or floats, where True or 1 means remove the element 
        
        :Optional Input:

            :param others: list of other attributes that mask should be applied to(e.g., 'pixel') (default = [])

        :Output:

            Spectrum object has the flagged pixels removed from  wave, flux, noise arrays, and optional vectors

        :Example:
           >>> import splat, numpy
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> num = splat.numberList('1-10,18,30-50')
           >>> mask = [not(p in num) for p in numpy.arange(len(sp.wave))]
           >>> sp.remove(mask)
           >>> sp.showHistory()
                SPEX_PRISM spectrum successfully loaded
                Mask applied to remove 32 pixels
        '''
        msk = copy.deepcopy(mask)
        if len(msk) != len(self.flux):
            print('\nWarning: mask must be same length ({}) as wave/flux arrays ({}); not removing any pixels'.format(len(msk),len(self.wave)))
            return
        if isinstance(msk[0],float): msk = [int(x) for x in msk]
        if isinstance(msk[0],int): msk = [True if x==1 else False for x in msk]
        if not isinstance(msk[0],bool): print('\nWarning: cannot interpret mask {}; not removing any pixels'.format(mask))

# invert mask
        msk = numpy.array([not x for x in msk])
#        self.wave = (numpy.array(self.wave.value)[msk])*self.wave.unit
#        self.flux = (numpy.array(self.flux.value)[msk])*self.flux.unit
#        self.noise = (numpy.array(self.noise.value)[msk])*self.noise.unit
        self.wave = self.wave[msk]
        self.flux = self.flux[msk]
        self.noise = self.noise[msk]
        self.variance = self.noise**2
        self.snr = self.computeSN()

        if len(others) > 0:
            for k in others:
                if k in self.original.__dict__.keys():
                    try:
                        setattr(self,k,getattr(self,k)[msk])
                    except:
                        pass

        cnt = numpy.sum([1 if x == False else 0 for x in msk])
        self.history.append('Mask applied to remove {} pixels'.format(cnt))
        return


    def maskFlux(self,mask,others=[]):
        '''
        :Purpose: 

            Masks elements of flux (set to NaN) as specified by the True/False array

        :Required Inputs:

            :param **mask**: An array of booleans, ints or floats, where True or 1 means remove the element 
        
        :Optional Inputs:

            :param **others**: list of other attributes that mask should be applied (e.g., 'pixel') (default = [])

        :Output:

            Spectrum object has the flagged pixels set to nan in flux arrays

        :Example:
           >>> import splat,numpy
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> num = splat.numberList('1-10,18,30-50')
           >>> mask = [not(p in num) for p in numpy.arange(len(sp.wave))]
           >>> sp.remove(mask)
           >>> sp.showHistory()
                SPEX_PRISM spectrum successfully loaded
                Mask applied to remove 32 pixels
        '''
        msk = copy.deepcopy(mask)
        if len(msk) != len(self.flux):
            print('\nWarning: mask must be same length ({}) as flux arrays ({}); not removing any pixels'.format(len(msk),len(self.flux)))
            return
        if isinstance(msk[0],float): msk = [int(x) for x in msk]
        if isinstance(msk[0],int): msk = [True if x==1 else False for x in msk]
        if not isinstance(msk[0],bool): print('\nWarning: cannot interpret mask {}; not removing any pixels'.format(mask))

# mask pixels
        self.flux[msk] = numpy.nan

        if len(others) > 0:
            for k in others:
                if k in self.original.__dict__.keys():
                    try:
                        x = getattr(self,k)
                        x[msk] = numpy.nan
                        setattr(self,k,x)
                    except:
                        pass

        cnt = numpy.sum([1 if x == False else 0 for x in msk])
        self.history.append('Masking applied to {} pixels'.format(cnt))
        return

    def reset(self):
        '''
        :Purpose: 

            Restores a Spectrum to its original read-in state, removing scaling and smoothing. 

        :Required Inputs:

            None
        
        :Optional Inputs:

            None

        :Output:

            Spectrum object is restored to original parameters
        
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

        self.history.append('Returned to original state')
        self.original = copy.deepcopy(self)
        return


    def res(self,npixel=1.,method='median'):
        '''
        :Purpose: 

            Estimate the resolution of the data from the wavelength array

        :Required Inputs: 

            None

        :Optional Inputs: 

            :param: npixel = 1: number of pixels per resolution element
            :param: method = 'median': what statistic to report back; can be mean, average, median, min, or max

        :Outputs: 

            The resolution of the spectrum (single number)
        
        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.resolution()
        '''

        res = numpy.absolute(numpy.array([0.5*(self.wave[i+1]+self.wave[i])/(self.wave[i+1]-self.wave[i]) for i in range(len(self.wave)-1)])/float(npixel))

        if method=='median': return numpy.nanmedian(res)
        elif method=='mean' or method=='average': return numpy.nanmean(res)
        elif method=='min': return numpy.nanmin(res)
        elif method=='max': return numpy.nanmax(res)
        else: raise ValueError('\nCould not interpret method = {}'.format(method))


    def scale(self,factor):
        '''
        :Purpose: 

            Scales a Spectrum object's flux and noise values by a constant factor. 


        :Required Inputs:

            :param factor: A floating point number used to scale the Spectrum object

        :Optional Inputs:

            None

        :Output: 

            None; spectrum is scaled

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
        self.fscale = 'Scaled'
        self.history.append('Spectrum scaled by a factor of {}'.format(factor))
        return

    def setNoise(self,rng=[],floor=0.):
        '''
        :Purpose: 

            Sets the noise and variance array of the spectrum based on rough counting statistics; 
            this routine is useful if the input spectrum has no input noise array

        :Required Inputs:

            One of the following must be provided:

            * *rng*: 2-element list specifying the range
            * *floor* = 0.: floor of uncertainty in spectrum flux units

        :Optional Inputs:

            None

        :Output: 

            Updates the noise and variance elements of Spectrum class

        :Example:
           >>> import splat, numpy
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.normalize()
           >>> sp.computeSN()  # reports 25.9118
           >>> sp.noise = [numpy.nan for n in sp.noise]*sp.noise.unit
           >>> sp.computeSN()  # reports nan
           >>> sp.setNoise([floor=1.e-3) 
           >>> sp.computeSN()  # reports 25.5425
        '''
        if len(rng) != 2 and floor == 0.:
            print('\nMust specify a noise sampling range as a 2-element list or the floor flux density')
            return
        if floor == 0.:
            if not isUnit(rng):
                rng = rng*self.wave.unit
            rng = (rng.to(self.wave.unit)).value
            if rng[0] < numpy.nanmin(self.wave.value[1:]): rng[0] = numpy.nanmin(self.wave.value[1:])
            if rng[1] > numpy.nanmax(self.wave.value[:-1]): rng[1] = numpy.nanmax(self.wave.value[:-1])
            floor = numpy.nanmedian(numpy.array(self.flux.value)[numpy.where(numpy.logical_and(self.wave.value>rng[0],self.wave.value<rng[1]))])
        self.noise = (((self.flux.value/floor)+1.)**0.5)*floor*self.flux.unit
        self.variance = self.noise**2

        return

    def showHistory(self):
        '''
        :Purpose: 

            Report history of actions taken on a Spectrum object. 
            This can also be retrieved by printing the attribute Spectrum.history

        :Required Inputs:

            None

        :Optional Inputs:

            None

        :Output: 

            List of actions taken on spectrum

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
        :Purpose: 

            Smoothes a spectrum either by selecting a constant slit width 
            (smooth in spectral dispersion space), pixel width (smooth in pixel space) 
            or resolution (smooth in velocity space). One of these options must be selected 
            for any smoothing to happen. Changes spectrum directly.

        :Required Inputs:

            None

        :Optional Inputs:

            :param method: the type of smoothing window to use; see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details (default = 'hanning')
            :param resolution: Constant resolution to smooth to; see `_smoothToResolution()`_ (default = None)
            :param slitPixelWidth: Number of pixels to smooth in pixel space; see `_smoothToSlitPixelWidth()`_ (default = None)
            :param slitWidth: Number of pixels to smooth in angular space; see `_smoothToPixelWidth()`_ (default = None)

        .. _`_smoothToResolution()` : api.html#splat.core.Spectrum._smoothToResolution
        .. _`_smoothToPixelWidth()` : api.html#splat.core.Spectrum._smoothToPixelWidth
        .. _`_smoothToSlitPixelWidth()` : api.html#splat.core.Spectrum._smoothToSlitPixelWidth


        :Output: 

            None: spectrum is smoothed

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.smooth(resolution=30)
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

    def _smoothToResolution(self,res,oversample=10.,method='hamming',**kwargs):
        '''
        :Purpose: 

            Smoothes a spectrum to a constant or resolution (smooth in velocity space). 
            Note that no smoothing is done if requested resolution is greater than the current resolution

        :Required Inputs:

            :param res: desired resolution (lambda/delta_lambde)

        :Optional Inputs:

            :param overscale: number of samples in the smoothing window (default = 10)

            see other optional keywords in _`Spectrum.smooth()`

        .. _`Spectrum.smooth()` : api.html#splat.core.Spectrum.smooth

        :Outputs:
            None; spectrum is smoothed

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

        method = kwargs.get('method','hamming')
        kwargs['method'] = method

# add in resolution keyword if not present
        if 'resolution' not in list(self.__dict__.keys()):
            i = int(0.5*len(self.wave.value))
            self.resolution = self.wave.value[i]/numpy.abs(self.wave.value[i]-self.wave.value[i+1])

# do nothing if requested resolution is higher than current resolution
        if res <= self.resolution:
# sample onto a constant resolution grid at 5x current resolution
            r = res*oversample
            wave_range = self.waveRange()
            wave_range = [w.value for w in wave_range]
            npix = numpy.floor(numpy.log(numpy.nanmax(wave_range)/numpy.nanmin(wave_range))/numpy.log(1.+1./res))
            wave_out = numpy.array([numpy.nanmin(wave_range)*(1.+1./res)**i for i in numpy.arange(npix)])
            a = numpy.linspace(0.,len(wave_out)-1,len(wave_out))
            b = numpy.linspace(0.,len(wave_out)-1,int(oversample*len(wave_out)))
            f = interp1d(a,wave_out)
            wave_sample = f(b)

            if len(self.wave) <= len(wave_sample):
                f = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
                v = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
                flx_sample = f(wave_sample)
                var_sample = v(wave_sample)
            else:
                flx_sample = integralResample(self.wave.value,self.flux.value,wave_sample)
                var_sample = integralResample(self.wave.value,self.variance.value,wave_sample)

# now convolve a function to smooth resampled spectrum
            window = signal.get_window(method,2.*numpy.round(oversample))
            neff = numpy.sum(window)/numpy.nanmax(window)        # effective number of pixels
            flx_smooth = signal.convolve(flx_sample, window/numpy.sum(window), mode='same')
            var_smooth = signal.convolve(var_sample, window/numpy.sum(window), mode='same')/neff
# resample back to original wavelength grid
            wave_final = numpy.array(self.wave.value)
            wave_final = wave_final[wave_final <= numpy.max(wave_sample)]
            wave_final = wave_final[wave_final >= numpy.min(wave_sample)]
            if len(wave_final) >= len(wave_sample):
                f = interp1d(wave_sample,flx_smooth,bounds_error=False,fill_value=0.)
                v = interp1d(wave_sample,var_smooth,bounds_error=False,fill_value=0.)
                flx_final = f(wave_final)
                var_final = v(wave_final)
            else:
                flx_final = integralResample(wave_sample,flx_smooth,wave_final)
                var_final = integralResample(wave_sample,var_smooth,wave_final)

#            f = interp1d(wave_sample,flx_smooth,bounds_error=False,fill_value=0.)
#            v = interp1d(wave_sample,var_smooth,bounds_error=False,fill_value=0.)
#            self.flux = f(self.wave.value)*self.funit
#            self.variance = v(self.wave.value)*self.funit**2
            self.wave = wave_final*self.wave.unit
            self.flux = flx_final*self.flux.unit
            self.variance = var_final*(self.flux.unit**2)
            self.noise = numpy.array([ns**0.5 for ns in self.variance.value])*self.flux.unit
            self.snr = self.computeSN()
#            self.slitpixelwidth = self.slitpixelwidth*self.resolution/res
#            self.slitwidth = self.slitwidth*self.resolution/res
            self.resolution = res
            self.history.append('Smoothed to a constant resolution of {}'.format(self.resolution))
        else:
            print('\nTarget resolution {} greater than current resolution {}; no change made'.format(res,self.resolution))
        return

    def _smoothToSlitPixelWidth(self,width,**kwargs):
        '''
        :Purpose: 

            Smoothes a spectrum to a constant slit pixel width (smooth in pixel space). C
            Note that no smoothing is done if requested width is greater than the current slit width.

        :Required Inputs:

            :param width: smoothing scale in pixels

        :Optional Inputs:

            see other optional keywords in _`Spectrum.smooth()`

        .. _`Spectrum.smooth()` : api.html#splat.core.Spectrum.smooth

        :Outputs:
            None; spectrum is smoothed

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
        if width > 2.:
# convolve a function to smooth spectrum
            window = signal.get_window(method,int(numpy.round(width)))
            neff = numpy.sum(window)/numpy.nanmax(window)        # effective number of pixels
            self.flux = signal.convolve(self.flux.value, window/numpy.sum(window), mode='same')*self.funit
            self.variance = signal.convolve(self.variance.value, window/numpy.sum(window), mode='same')/neff*(self.funit**2)
            self.noise = [n**0.5 for n in self.variance.value]*self.funit
            self.snr = self.computeSN()
            self.resolution = self.resolution*self.slitpixelwidth/width
            self.slitwidth = self.slitwidth*width/self.slitpixelwidth
            self.slitpixelwidth = width
            self.history.append('Smoothed to slit pixel width of {}'.format(self.slitpixelwidth))
        else:
            print('\nTarget slit width {} is less than 2 pixels; no change made'.format(width))
        return

    def _smoothToSlitWidth(self,width,**kwargs):
        '''
        :Purpose: 

            Smoothes a spectrum to a constant slit angular width (smooth in dispersion space). 
            Note that no smoothing is done if requested width is greater than the current slit width.

        :Required Inputs:

            :param width: smoothing scale in arcseconds

        :Optional Inputs:

            see other optional keywords in _`Spectrum.smooth()`

        .. _`Spectrum.smooth()` : api.html#splat.core.Spectrum.smooth

        :Outputs:
            None; spectrum is smoothed

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
        if not isUnit(width): width=width*u.arcsec
        pwidth = self.slitpixelwidth*(width/self.slitwidth).value
        self._smoothToSlitPixelWidth(pwidth,**kwargs)
        return


    def toSurface(self,radius):
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
        if not isUnit(r): r=r*const.R_sun
        self.scale((((10.*u.pc/r).to(u.m/u.m)).value)**2)
        self.history.append('Converted to surface fluxes assuming a radius of {} solar radii'.format((r/const.R_sun).to(u.m/u.m)))
        self.fscale = 'Surface'
        return


    def toBrightnessTemperature(self,limbdarkening=False,limbdarkeningcoeff = 0.7):
        '''
        :Purpose: Convert to surface fluxes given a radius, assuming at absolute fluxes
        .. note: UNTESTED
        '''
        if self.fscale == 'Temperature':
            return
        if self.fscale != 'Surface':
            print('To convert to brightness temperature you must first scale spectrum to surface flux units')
            return
        if self.funit != DEFAULT_FLUX_UNIT: self.toFlam()
        fs = copy.deepcopy(self.flux).to(u.erg/u.s/u.cm**3)
        fse = copy.deepcopy(self.noise).to(u.erg/u.s/u.cm**3)
        w = copy.deepcopy(self.wave).to(u.cm)
        x = (2.*numpy.pi*const.h.to(u.erg*u.s)*(const.c.to(u.cm/u.s)**2)/(fs*(w**5))).to(u.m/u.m).value
        self.temperature = (const.h*const.c/(const.k_B*w)).to(u.K)/numpy.log(1.+x)
#        self.temperature_unc = self.flux*x/((1.+x)*numpy.log(1.+x))*((fs/fse).to(u.m/u.m))
        self.flux = (const.h*const.c/(const.k_B*w)).to(u.K)/numpy.log(1.+x)
        self.noise = self.flux*x/((1.+x)*numpy.log(1.+x))*((fs/fse).to(u.m/u.m))
        self.variance = self.noise**2
        self.snr = self.computeSN()
        self.history.append('Converted to brightness temperature')
        self.flabel = 'Temperature'
        self.fscale = ''
        self.funit = u.K
        return

    def toTemperature(self):
        self.brightnessTemperature()


    def trim(self,rng,**kwargs):
        '''
        :Purpose: 
            Trims a spectrum to be within a certain wavelength range or set of ranges. 
            Data outside of these ranges are excised from the wave, flux and noise arrays. 
            The full spectrum can be restored with the reset() procedure.

        :Required Inputs: 

            :param range: the range(s) over which the spectrum is retained - a series of nested 2-element arrays

        :Optional Inputs: 

            None

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

        mask = numpy.zeros(len(self.wave))

# some code to deal with various possibilities, ultimately leading to [ [r1a,r1b], [r2a,r2b], ...]
# convert a unit-ed quantity
        if isUnit(rng):
            try:
                rng.to(self.wave.unit).value
            except:
                raise ValueError('Could not convert trim range unit {} to spectrum wavelength unit {}'.format(rng.unit,self.wave.unit))

# single number = turn into small range
        if isinstance(rng,float):
            rng = [rng-0.01*(numpy.nanmax(self.wave.value)-numpy.nanmin(self.wave.value)),rng+0.01*(numpy.nanmax(self.wave.value)-numpy.nanmin(self.wave.value))]

        if isUnit(rng[0]):
            try:
                rng = [r.to(self.wave.unit).value for r in rng]
            except:
                raise ValueError('Could not convert trim range unit {} to spectrum wavelength unit {}'.format(rng[0].unit,self.wave.unit))

        if not isinstance(rng[0],list):
            rng = [rng]

        for r in rng:
#            if isUnit(r):
#                r=[(x*r.unit).to(self.wave.unit) for x in r]
#            if not isUnit(r[0]):
#                r = [x*self.wave.unit for x in r]
            if isUnit(r[0]):
                try:
                    r = [x.to(self.wave.unit).value for x in r]
                except:
                    raise ValueError('Could not convert trim range unit {} to spectrum wavelength unit {}'.format(r[0].unit,self.wave.unit))
            mask[numpy.where(numpy.logical_and(self.wave.value > r[0],self.wave.value < r[1]))] = 1
#        w = numpy.where(mask == 1)
        self.wave = self.wave[mask == 1]
        self.flux = self.flux[mask == 1]
        self.noise = self.noise[mask == 1]
        self.variance = self.variance[mask == 1]
#        self.flam = self.flux
#        self.nu = self.wave.to('Hz',equivalencies=u.spectral())
#        self.fnu = self.flux.to('Jy',equivalencies=u.spectral_density(self.wave))
#        self.noisenu = self.noise.to('Jy',equivalencies=u.spectral_density(self.wave))
        self.snr = self.computeSN()
        self.history.append('Spectrum trimmed to range {}'.format(rng))
        return

    def updateSourceInfo(self,verbose=False,radius=10.*u.arcsec,**kwargs):
        '''
        :Purpose: 
            Updates the source information for a spectrum object based on the SPLAT Source Database
            Uses either the spectrum object's or user supplied NAME, DESIGNATION, COORDINATE or RA & DEC
            If none of these are supplied, no search is done

        :Required Inputs:
            None 

        :Optional Inputs:
            *name*: Name of source
            *designation*: designation of source in format Jhhmmss[.]ss±ddmmss[.]ss
            *shortname*: shortname desigation of source in format Jhhmm±ddmm
            *coordinate*: coordinate of object, in astropy SkyCoord format or transferable via properCoordinates()
            *ra*, *dec*: RA and Declination of source in degrees
            *radius*: search radius (default = 10 arcseconds)
            *verbose*: set to True to provide feedback

        :Output:
            Spectrum object will have the Source Library keywords added or updated

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

        if len(list(kwargs.keys())) == 0:
            print('\nYou must provide some search terms for searchLibrary() to function')
            return

        spattr = self.__dict__.keys()
        spattr = [s.lower() for s in spattr]
        if 'designation' not in spattr: desig = False 
        else: desig = self.designation
        desig = kwargs.get('designation',desig)
        desig = kwargs.get('desig',desig)
        if 'shortname' not in spattr: sname = False 
        else: sname = self.shortname
        sname = kwargs.get('shortname',sname)
        sname = kwargs.get('sname',sname)
        sname = kwargs.get('short',sname)
        if 'coordinate' not in spattr: coord = False 
        else: coord = self.coordinate
        coord = kwargs.get('coordinate',coord)
        coord = kwargs.get('coord',coord)
        if 'name' not in spattr: name = False 
        else: name = self.name
        name = kwargs.get('name',name)
        if 'ra' not in spattr: ra = False 
        else: ra = self.ra
        ra = kwargs.get('ra',ra)
        if not isinstance(ra,float): ra = False
        if 'dec' not in spattr: dec = False 
        else: dec = self.dec
        dec = kwargs.get('dec',dec)
        if not isinstance(dec,float): dec = False

        s = pandas.DataFrame()

        if desig != False:
            s = searchLibrary(designation=desig)
            if len(s) == 0:
                s = searchLibrary(coordinate=desig,radius=radius)
        if name != False and len(s) == 0:
            s = searchLibrary(name=name)
        if sname != False and len(s) == 0:
            s = searchLibrary(shortname=sname)
        if coord != False and len(s) == 0:
            s = searchLibrary(coordinate=coord,radius=radius)
        if ra != False and dec != False and len(s) == 0:
            s = searchLibrary(coordinate=properCoordinates([ra,dec]),radius=radius)

        if len(s) == 0:
            if verbose: print('\nNo objects found in the SPLAT Source Database')
        else:
            if verbose: print('\nMatching to source {} designation {}'.format(s['NAME'].iloc[0],s['DESIGNATION'].iloc[0]))
            for k in list(DB_SOURCES.columns):
                if k in list(s.columns):
                    setattr(self,k.lower(),s[k].iloc[0])
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
def stitch(s1,s2,rng=[],verbose=False,scale=True,**kwargs):
    '''
    :Purpose: 

        Stitches together two spectra covering different wavelength scales.

    :Required Inputs: 

        :param s1: first spectrum object
        :param s2: second spectrum object (not necessarily in order)

    :Optional Inputs: 

        :param rng: range over which spectra are relatively scaled and combined (if desired); if neither spectrum cover this range, then no relative scaling is done and spectra are combined whereever they overlap (default: [])
        :param scale: set to True to relatively scale spectra over rng; if rng is not provided, or spectra to not overlap, this is automatially False (default: True)
        :param wunit: wavelength unit of final spectrum (default: wavelength unit of s1)
        :param funit: flux unit of final spectrum (default: flux unit of s1)
        :param verbose: set to True to provide feedback (default: False)

    :Output: 

        New spectrum object of stitched spectrum

    :Example:
       >>> import splat
       >>> spopt = splat.Spectrum(file='myopticalspectrum.fits')
       >>> spnir = splat.Spectrum(file='myopticalspectrum.fits')
       >>> sp = splat.stitch(spopt,spnir,rng=[0.8,0.9],trim=[0.35,2.4])
    '''

# parameters
    scaleflag = kwargs.get('scaleflag',scale)

# generate copies of spectrum objects
    sp1 = copy.deepcopy(s1)
    sp2 = copy.deepcopy(s2)
    
# assert common units
    wunit = kwargs.get('wunit',sp1.wave.unit)
    sp1.toWaveUnit(wunit)
    sp2.toWaveUnit(wunit)
    funit = kwargs.get('funit',sp1.flux.unit)
    sp1.toFluxUnit(funit)
    sp2.toFluxUnit(funit)

# reorder to shortest wavelength spectrum first?

# check range
    rngflag = True
    if len(rng) == 0:
        rng = [numpy.nanmax([numpy.nanmin(sp1.wave.value),numpy.nanmin(sp2.wave.value)]),\
            numpy.nanmin([numpy.nanmax(sp1.wave.value),numpy.nanmax(sp2.wave.value)])]*wunit
    if not isUnit(rng):
        rng=rng*wunit
    rng = (rng.to(wunit)).value
    if rng[0] < numpy.nanmin(sp1.wave.value[1:]): rng[0] = numpy.nanmin(sp1.wave.value[1:])
    if rng[0] < numpy.nanmin(sp2.wave.value[1:]): rng[0] = numpy.nanmin(sp2.wave.value[1:])
    if rng[1] > numpy.nanmax(sp1.wave.value[:-1]): rng[1] = numpy.nanmax(sp1.wave.value[:-1])
    if rng[1] > numpy.nanmax(sp2.wave.value[:-1]): rng[1] = numpy.nanmax(sp2.wave.value[:-1])
    if rng[0] >= rng[1]:
        if verbose== True: print('Stich region {} to {} does not overlap both spectra'.format(rng[0],rng[1]))
        rngflag = False
#        raise ValueError('Stich region {} to {} does not overlap both spectra'.format(rng[0],rng[1]))
#    print(rng,numpy.nanmin(sp1.wave.value),numpy.nanmax(sp1.wave.value),numpy.nanmin(sp2.wave.value),numpy.nanmax(sp2.wave.value))

# overlap and scaling
    if rngflag==True:

# interpolation of second spectrum
        f2r = interp1d(sp2.wave.value,sp2.flux.value)
        v2r = interp1d(sp2.wave.value,sp2.variance.value)

# find overlap region, assuming first spectrum sets the flux scale standard
# assume this minimizes chi^2 residuals
        w12 = numpy.where(numpy.logical_and(numpy.array(sp1.wave.value)>=rng[0],numpy.array(sp1.wave.value)<=rng[1]))
        wv12 = numpy.array(sp1.wave.value)[w12]
        flx12 = numpy.array(sp1.flux.value)[w12]
        var12 = numpy.array(sp1.variance.value)[w12]
        sp1mid = Spectrum(wave=wv12*wunit,flux=flx12*funit,noise=(var12**0.5)*funit)
        sp2mid = Spectrum(wave=wv12*wunit,flux=f2r(sp1mid.wave.value)*funit,noise=(v2r(sp1mid.wave.value)**0.5)*funit)

        if scaleflag == True:
            chi,scl = compareSpectra(sp1mid,sp2mid)

            vtot = numpy.zeros(len(wv12))
            vflag = True
            if not numpy.isnan(numpy.nanmedian(var12)): vtot=var12
            if not numpy.isnan(numpy.nanmedian(v2r(wv12))): vtot=vtot+v2r(wv12)

            if numpy.nanmedian(vtot) == 0.: vflag = False
            if vflag == True:
                scl = numpy.nansum(flx12*f2r(wv12)/vtot)/numpy.nansum(f2r(wv12)**2/vtot)
            else:
                scl = numpy.nansum(flx12*f2r(wv12))/numpy.nansum(f2r(wv12)**2)
            sp2.scale(scl)

        f2r = interp1d(sp2.wave.value,sp2.flux.value)
        v2r = interp1d(sp2.wave.value,sp2.variance.value)

# piece back together, starting with segments with minimum wavelength
    if numpy.nanmin(sp1.wave.value) < numpy.nanmin(sp2.wave.value):
        w1 = numpy.where(numpy.array(sp1.wave.value)<rng[0])
        wave = numpy.array(sp1.wave.value)[w1]
        flux = numpy.array(sp1.flux.value)[w1]
        variance = numpy.array(sp1.variance.value)[w1]
    else:
        w1 = numpy.where(numpy.array(sp2.wave.value)<rng[0])
        wave = numpy.array(sp2.wave.value)[w1]
        flux = numpy.array(sp2.flux.value)[w1]
        variance = numpy.array(sp2.variance.value)[w1]

# mixed region
    if rngflag == True:
        v1 = var12
        v2 = v2r(wv12)
        if numpy.isnan(numpy.nanmedian(v1)): v1 = v2
        if numpy.isnan(numpy.nanmedian(v2)): v2 = v1
        if vflag == 0:
            flxmid = (flx12/v1+f2r(wv12)/v2)/(1./v1+1./v2)
        else:
            flxmid = 0.5*(flx12+f2r(wv12))
        varmid = 1./(1./v1+1./v2)
        wave = numpy.append(wave,wv12)
        flux = numpy.append(flux,flxmid)
        variance = numpy.append(variance,varmid)

# second spectrum
    if numpy.nanmin(sp1.wave.value) < numpy.nanmin(sp2.wave.value):
        w2 = numpy.where(numpy.array(sp2.wave.value)>rng[1])
        wave = numpy.append(wave,numpy.array(sp2.wave.value)[w2])
        flux = numpy.append(flux,numpy.array(sp2.flux.value)[w2])
        variance = numpy.append(variance,numpy.array(sp2.variance.value)[w2])
    else:
        w2 = numpy.where(numpy.array(sp1.wave.value)>rng[1])
        wave = numpy.append(wave,numpy.array(sp1.wave.value)[w2])
        flux = numpy.append(flux,numpy.array(sp1.flux.value)[w2])
        variance = numpy.append(variance,numpy.array(sp1.variance.value)[w2])

# put back units
    wave = wave*wunit
    flux = flux*funit
    variance = variance*(funit**2)

# create new spectrum object containing combined spectrum
    sp = splat.Spectrum(wave=wave,flux=flux,noise=variance**0.5)
#    sp.wave = wave
#    sp.flux = flux
#    sp.variance = variance
#    sp.noise = sp.variance**0.5
#    sp.snr = sp.computeSN()
    sp.name = 'Stitched spectrum of {} and {}'.format(sp1.name,sp2.name)

# trim if desired - REMOVED
#    if len(trim) == 2:
#        if not isUnit(trim):
#            trim=trim*sp.wave.unit
#        trim.to(sp.wave.unit).value
#        sp.trim(trim)
        
    return sp


#####################################################
#################   DATA ACCESS   ###################
#####################################################


def getSpectrum(getList=False, limit=0, *args, **kwargs):
    '''
    :Purpose: 

        Gets a spectrum from the SPLAT library using various selection criteria. Calls searchLibrary_ to select spectra; if any found it routines an array of Spectrum objects, otherwise an empty array. 

    .. _searchLibrary : api.html#splat.core.searchLibrary

    :Output: 

        An array of Spectrum objects that satisfy the search criteria

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
        Retrieving 1 file
    >>> sparr = splat.getSpectrum(spt='M7')
        Retrieving 120 files
    >>> sparr = splat.getSpectrum(spt='T5',young=True)
        No files match search criteria
    '''

    if kwargs.get('lucky',False) == True: kwargs['published'] = True
    result = []
    kwargs['output'] = 'all'
    search = searchLibrary(*args, **kwargs)

    if len(search) > 0:
        files = []
        if len(search) == 1:
            files.append(search['DATA_FILE'].iloc[0])
        else:
            for i,x in enumerate(search['DATA_FILE']):
                files.append(search['DATA_FILE'].iloc[i])


# return just the filenames
        if getList == True:
            return files

        if len(files) == 1:
            if kwargs.get('lucky',False) == True:
                print('\nRetrieving 1 lucky file\n')
            else:
                print('\nRetrieving 1 file\n')
            skwargs = search.iloc[0].to_dict()
            result.append(Spectrum(files[0],**skwargs))
        else:
#            if (kwargs.get('lucky',False) == True):
#                print('\nRetrieving 1 lucky file\n')
#                ind = numpy.random.choice(numpy.arange(len(files)))
#                print(x)
#                result.append(Spectrum(files[ind],header=search[ind]))
#            else:
            if limit != 0 and limit < len(files):
                files = files[:limit]
                search = search.iloc[:limit]
            print('\nRetrieving {} files\n'.format(len(files)))
            for i,x in enumerate(files):
                skwargs = search.iloc[i].to_dict()
                result.append(Spectrum(x,header=dict(search.iloc[i]),**skwargs))

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
        stdtype = 'extreme subdwarf'
    elif kwargs.get('sd',False) or 'sd' in sptstr:
        stds = STDS_SD_SPEX
        kys = STDS_SD_SPEX_KEYS
        subclass = 'sd'
        stdtype = 'subdwarf'
    elif kwargs.get('vlg',False) or 'gamma' in sptstr:
        stds = STDS_VLG_SPEX
        kys = STDS_VLG_SPEX_KEYS
        subclass = 'gamma'
        stdtype = 'very low gravity dwarf'
    elif kwargs.get('intg',False) or 'beta' in sptstr:
        stds = STDS_INTG_SPEX
        kys = STDS_INTG_SPEX_KEYS
        subclass = 'beta'
        stdtype = 'intermediate gravity dwarf'
    else:
        stds = STDS_DWARF_SPEX
        kys = STDS_DWARF_SPEX_KEYS
        subclass = ''
        stdtype = 'dwarf'

    spt = typeToNum(spt,subclass=subclass)

# not a valid subtype
    if spt not in list(kys.keys()):
        print('Type {} is not in {} standards: try one of the following:'.format(spt,stdtype))
        print(sorted(list(kys.keys())))
        return Spectrum()

# not yet read in
    if spt not in list(stds.keys()):
        stds[spt] = Spectrum(kys[spt])
    
    return stds[spt]


def initializeStandards(**kwargs):
    return(initiateStandards())


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
    elif kwargs.get('dsd',False):
        stds = STDS_DSD_SPEX
        kys = STDS_DSD_SPEX_KEYS
    elif kwargs.get('esd',False):
        stds = STDS_ESD_SPEX
        kys = STDS_ESD_SPEX_KEYS
    elif kwargs.get('vlg',False):
        stds = STDS_VLG_SPEX
        kys = STDS_VLG_SPEX_KEYS
    elif kwargs.get('intg',False):
        stds = STDS_INTG_SPEX
        kys = STDS_INTG_SPEX_KEYS
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
#    sdb = copy.deepcopy(DB_SOURCES)
#    sdb['SELECT'] = [x in keys for x in sdb['SOURCE_KEY']]

    sdb = DB_SOURCES[[x in keys for x in DB_SOURCES['SOURCE_KEY']]]
#    if sum(sdb['SELECT']) == 0.:
    if len(sdb) == 0.:
        if kwargs.get('verbose',True) == True: print('No sources found with source key(s) = {}'.format(*keys))
        return False
    else:
#        db = sdb[:][numpy.where(sdb['SELECT']==1)]
        return sdb


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

# vectorize
    if isinstance(keys,list) == False:
        keys = [keys]

#    sdb = copy.deepcopy(DB_SPECTRA)
#    sdb['SELECT'] = [x in keys for x in sdb['DATA_KEY']]

#    if sum(sdb['SELECT']) == 0.:
#        if verbose: print('No spectra found with spectrum key {}'.format(keys[0]))
#        return False
#    else:
#        s2db = copy.deepcopy(DB_SOURCES)
#        db = join(sdb[:][numpy.where(sdb['SELECT']==1)],s2db,keys='SOURCE_KEY')
#        return db

    sdb = DB_SPECTRA[[x in keys for x in DB_SPECTRA['DATA_KEY']]]
    if len(sdb) == 0.:
        if kwargs.get('verbose',True) == True: print('No sources found with spectrum key(s) = {}'.format(*keys))
        return False
    else:
        return sdb


def searchLibrary(radius=10., instrument='SPEX-PRISM', *args, **kwargs):
    '''
    :Purpose: 

    Searches the SpeX database based on a series of keywords; 
    returns an astropy Table with the source and spectral information corresponding to the selected sources

    :Required Parameters:

    None

    :Optional Parameters:
        * :param name: search by source name (e.g., ``name = 'Gliese 570D'``); default = None
        * :param shortname: search be short name or list of short names (e.g. ``shortname = 'J1457-2124'``); default = None 
        * :param exclude_shortname: exclude a list of short names (e.g. ``excludeshortname = 'J1457-2124'``); default = None 
        * :param designation: search by full designation (e.g., ``designation = 'J11040127+1959217'``); default = None 
        * :param coordinate: search around a coordinate by a radius specified by radius keyword (e.g., ``coordinate = [180.,+30.], radius = 10.``); coordinate can be an astropy SkyCoord, list of [RA, DEC], or designation; default = None
        * :param source, source_key, id, id_key: search for a specific spectra based on single or list of id number or list of numbers (e.g., ``source_key = [10002,10005]``)
        * :param data_key: search for a specific spectra based on single or list of data key numbers (e.g., ``data_key = [10002,10005]``)
        * :param exclude_data_key: exclude specific spectra based on single or list of data key numbers (e.g., ``exclude_data_key = [10002,10005]``)
        * :param file: search by specific filename or list of filenames (e.g., ``file = 'myspectrum.fits``)
        * :param exclude_file: exclude by specific filename or list of filenames (e.g., ``exclude_file = 'myspectrum.fits``)
        * :param radius: search radius in arcseconds for coordinate search; default = 10.
        * :param spt: search by spectral type, which by default is the SpeX-based type; single value is exact, two-element array gives range (e.g., ``spt = 'M7'`` or ``spt = [24,39]``); can also specify:
            * :param spex_spt: same as ``spt``
            * :param opt_spt: same as ``spt`` for literature optical spectral types
            * :param nir_spt: same as ``spt`` for literature NIR spectral types
            * :param lit_spt: same as ``spt`` for literature spectral types from SIMBAD
        * :param jmag, hmag, kmag: select based on faint limit or range of J, H or Ks magnitudes (e.g., ``jmag = 11`` or ``jmag = [12,15]``); default = None
        * :param snr: search on minimum or range of S/N ratios (e.g., ``snr = 30.`` or ``snr = [50.,100.]``)
        * :param subdwarf, young, lowg, binary, sbinary, red, blue, giant, wd, standard, companion, peculiar: classes to search on or exclude (e.g., ``young = True``, ``giant = False``)
        * :param cluster: select sources in specfic clusters (e.g., ``cluster = 'TWA'``)
        * :param giant_class or luminosity_class: select sources based on luminosity class; this does not work for dwarfs (e.g., ``luminosity_class = 'I'``)
        * :param subdwarf_class or metallicity_class: select sources based on metallicty class; (e.g., ``metallicity_class = 'sd'``)
        * :param date: search by observation date (e.g., ``date = '20040322'``) or range of dates (e.g., ``date=[20040301,20040330]``)
        * :param instrument: search by instrument; must be designated instrument in splat.INSTRUMENTS variable (e.g., ``instrument = 'SPEX_PRISM'``)
        * :param reference: search by list of references (bibcodes) (e.g., ``reference = '2011ApJS..197...19K'``)
        * :param logic, combine: search logic, can be ``and`` (default) or ``or``
        * :param sort: by default returned table is sorted by designation; set this parameter to a column name to sort on a different parameter
        * :param reverse: set to True to do a reverse sort (default = False)
        * :param list: if True, return just a list of the data files (default = False)
        * :param lucky: if True, return one randomly selected spectrum from the selected sample (default = False)
        * :param output: returns desired column of selected results (default = 'all')

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
    object_classes = numpy.unique(numpy.sort(numpy.array([str(x) for x in splat.DB_SOURCES['OBJECT_TYPE']])))
    object_classes = object_classes[numpy.where(object_classes != '0')]  # eliminate masked element
    object_classes = object_classes[numpy.where(object_classes != 'nan')]  # eliminate masked element
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
    source_db['SELECT'] = numpy.zeros(len(source_db['SOURCE_KEY']))
    count = 0.

# search by source key
    idkey = kwargs.get('sourcekey',False)
    idkey = kwargs.get('source_key',idkey)
    idkey = kwargs.get('source',idkey)
    idkey = kwargs.get('idkey',idkey)
    idkey = kwargs.get('id_key',idkey)
    idkey = kwargs.get('id',idkey)
    if idkey != False:
        if not isinstance(idkey,list):
            idkey = [idkey]
        if isinstance(idkey[0],str):
            idkey = [int(i) for i in idkey]
        for s in idkey:
            source_db['SELECT'][source_db['SOURCE_KEY'] == s] += 1
        count+=1.

# search by name
    if kwargs.get('name',False) != False:
        nm = kwargs['name']
        if isinstance(nm,str):
            nm = [nm]
        if len(nm) > 0:
            source_db['NAMEGEN'] = [n.lower().replace(' ','') for n in source_db['NAME']]
            for n in nm:
                source_db['SELECT'][source_db['NAMEGEN'] == n.lower().replace(' ','')] += 1
            count+=1.
            del source_db['NAMEGEN']

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
            source_db['SELECT'][source_db['SHORTNAME'] == sn] += 1
        count+=1.

# exclude by shortname
    sname = kwargs.get('exclude_source',False)
    sname = kwargs.get('exclude_shortname',sname)
    if sname != False and len(sname) > 0:
        if isinstance(sname,str):
            sname = [sname]
        for sn in sname:
            if sn[0].lower() != 'j':
                sn = 'J'+sn
#            t = numpy.sum(source_db['SELECT'][numpy.where(source_db['SHORTNAME'] != sn)])
            source_db['SELECT'][source_db['SHORTNAME'] != sn] += 1
#            if numpy.sum(source_db['SELECT'][numpy.where(source_db['SHORTNAME'] != sn)]) > t:
#                print('rejected '+sn)
        count+=1.

# search by reference list
    if kwargs.get('discovery_reference',False) != False:
        refer = kwargs['discovery_reference']
        if isinstance(refer,str):
            refer = [refer]
        for r in refer:
            source_db['SELECT'][source_db['DISCOVERY_REFERENCE'] == r] += 1
        count+=1.

# search by designation
    desig = kwargs.get('designation',False)
    desig = kwargs.get('coordinate',desig)
    desig = kwargs.get('coord',desig)
    if desig != False:
        try:
            cc = properCoordinates(desig)
        except:
            print('\nWarning: {} is not a proper coordinate'.format(desig))
        else:

# make sure you can compare skycoords
            if 'RA' not in list(source_db.keys()) or 'DEC' not in list(source_db.keys()):
                if 'DESIGNATION' in list(source_db.keys()):
                    coords = [designationToCoordinate(d) for d in source_db['DESIGNATION']]
                    source_db['RA'] = [c.ra.degree for c in coords]
                    source_db['DEC'] = [c.dec.degree for c in coords]

            source_db['COORDFLAG'] = numpy.zeros(len(source_db))
            source_db['COORDSEPR'] = [numpy.abs((cc.ra.degree-r)*numpy.cos(cc.dec.degree*numpy.pi/180.)*3600.) for r in source_db['RA']]
            source_db['COORDSEPD'] = [numpy.abs((cc.dec.degree-d)*3600.) for d in source_db['DEC']]
            r = (source_db['COORDSEPR']**2+source_db['COORDSEPD']**2)**0.5
            print(radius,numpy.nanmin(r),source_db['DESIGNATION'].iloc[numpy.argmin(r)])
            chk = source_db[source_db['COORDSEPR'] <= radius]
            chk = chk[chk['COORDSEPD'] <= radius]
            if len(chk) > 0:
                for k in list(chk.index): source_db['SELECT'].iloc[k] += 1
            count+=1.
            del source_db['COORDFLAG'], source_db['COORDSEPR'], source_db['COORDSEPD']

#            s = []
#            for i in numpy.arange(len(source_db['RA'])):
#                try:        # to deal with a blank string
#                    s.append(SkyCoord(ra=float(source_db['RA'][i])*u.degree,dec=float(source_db['DEC'][i])*u.degree,frame='icrs'))
#                except:
#                    s.append(SkyCoord(ra=numpy.nan*u.degree,dec=numpy.nan*u.degree,frame='icrs'))
#                if numpy.mod(i,len(source_db['RA'])/10.) < 1 and i != 0:
#                    print('\b{:.0f}%...'.format(100*i/len(source_db['RA'])))
#            source_db['SKYCOORD'] = s
#        print('measuring separations')
#        source_db['SEPARATION'] = [cc.separation(source_db['SKYCOORDS'][i]).arcsecond for i in numpy.arange(len(source_db['SKYCOORDS']))]
#        source_db['SEPARATION'] = [cc.separation(c).arcsecond for c in source_db['SKYCOORD']]
#        print('done')
#        source_db['SELECT'][source_db['SEPARATION'] <= radius] += 1
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
        if not isinstance(spt_range,list) == True:        # one value = only this type
            spt_range = [spt_range,spt_range]
        if isinstance(spt_range[0],str) == True:          # convert to numerical spt
            spt_range = [typeToNum(spt_range[0]),typeToNum(spt_range[1])]
        source_db['SPTN'] = [typeToNum(x) for x in source_db[spt_type]]
        source_db['SELECT'][numpy.logical_and(source_db['SPTN'] >= spt_range[0],source_db['SPTN'] <= spt_range[1])] += 1
        count+=1.

# search by magnitude range
    if kwargs.get('jmag',False) != False:
        mag = kwargs['jmag']
        if not isinstance(mag,list): mag = [0,mag]
        source_db['SELECT'][numpy.logical_and(source_db['J_2MASS'] >= mag[0],source_db['J_2MASS'] <= mag[1])] += 1
        count+=1.
    if kwargs.get('hmag',False) != False:
        mag = kwargs['hmag']
        if not isinstance(mag,list): mag = [0,mag]
        source_db['SELECT'][numpy.logical_and(source_db['H_2MASS'] >= mag[0],source_db['H_2MASS'] <= mag[1])] += 1
        count+=1.
    if kwargs.get('kmag',False) != False:
        mag = kwargs['kmag']
        if not isinstance(mag,list): mag = [0,mag]
        source_db['SELECT'][numpy.logical_and(source_db['KS_2MASS'] >= mag[0],source_db['KS_2MASS'] <= mag[1])] += 1
        count+=1.

# low surface gravity
    if (kwargs.get('lowg','') != ''):
#        source_db['LOWG'] = [not numpy.ma.is_masked(i) for i in source_db['GRAVITY_CLASS_OPTICAL']] or [not numpy.ma.is_masked(i) for i in source_db['GRAVITY_CLASS_NIR']]
        source_db['LOWG'] = [source_db['GRAVITY_CLASS_OPTICAL'][i]=='alpha' or source_db['GRAVITY_CLASS_OPTICAL'][i]=='beta' or source_db['GRAVITY_CLASS_OPTICAL'][i]=='gamma' or source_db['GRAVITY_CLASS_OPTICAL'][i]=='delta' or source_db['GRAVITY_CLASS_NIR'][i]=='INT-G' or source_db['GRAVITY_CLASS_NIR'][i]=='VL-G' or source_db['GRAVITY_CLASS_NIR'][i]=='LOW-G' for i in range(len(source_db))]
        source_db['SELECT'][source_db['LOWG'] == kwargs.get('lowg')] += 1
        count+=1.
        del source_db['LOWG']

# specific gravity class
    flag = kwargs.get('gravity_class','')
    flag = kwargs.get('gravity',flag)
    if (flag != ''):
#        source_db['SELECT'][numpy.ma.filled(source_db['GRAVITY_CLASS_OPTICAL'],'') == flag] += 1
#        source_db['SELECT'][numpy.ma.filled(source_db['GRAVITY_CLASS_NIR'],'') == flag] += 1
        source_db['SELECT'][source_db['GRAVITY_CLASS_OPTICAL'] == flag] += 1
        source_db['SELECT'][source_db['GRAVITY_CLASS_NIR'] == flag] += 1
        count+=1.

# young => member of a young cluster
    if (kwargs.get('young','') != ''):
        source_db['INCLUSTER'] = [str(x).lower() != 'nan' for x in source_db['CLUSTER']]
        source_db['SELECT'][source_db['INCLUSTER'] == kwargs.get('young')] += 1
        count+=1.
        del source_db['INCLUSTER']

# young => member of a young cluster
    if (kwargs.get('cluster','') != '' and isinstance(kwargs.get('cluster'),bool)):
        source_db['INCLUSTER'] = [str(x).lower() != 'nan' for x in source_db['CLUSTER']]
        source_db['SELECT'][source_db['INCLUSTER'] == kwargs.get('cluster')] += 1
        count+=1.
        del source_db['INCLUSTER']

# specific cluster
    if (kwargs.get('cluster','') != '' and isinstance(kwargs.get('cluster'),str)):
        source_db['CLUSTER_FLAG'] = [str(i).lower() == kwargs.get('cluster').lower() for i in source_db['CLUSTER']]
        source_db['SELECT'][source_db['CLUSTER_FLAG'] == True] += 1
        count+=1.
        del source_db['CLUSTER_FLAG']

# select out object classes
    for oc in object_classes:
        if kwargs.get(oc.lower(),'') != '':
            if (kwargs[oc.lower()] == True):
                source_db['SELECT'][source_db['OBJECT_TYPE'] == oc.upper()] += 1
                count+=1.
            if (kwargs[oc.lower()] == False):
                source_db['SELECT'][source_db['OBJECT_TYPE'] != oc.upper()] += 1
                count+=1.

# giant
#    if (kwargs.get('giant','') != ''):
#        source_db['GIANT'] = [not numpy.ma.is_masked(i) for i in source_db['LUMINOSITY_CLASS']]
#        source_db['SELECT'][numpy.where(source_db['GIANT'] == kwargs.get('giant'))] += 1
#        count+=1.

# luminosity class - this is not quite right
    lclass = kwargs.get('giant_class','') 
    lclass = kwargs.get('luminosity_class',lclass) 
    if lclass != '':
#        if 'GIANT' not in source_db.keys():
#            source_db['GIANT'] = [not numpy.ma.is_masked(i) for i in source_db['LUMINOSITY_CLASS']]
        source_db['GIANT_FLAG'] = [i.lower() == lclass.lower() for i in numpy.ma.filled(source_db['LUMINOSITY_CLASS'],'')]
        source_db['SELECT'][source_db['GIANT_FLAG'] == True] += 1
        count+=1.
        del source_db['GIANT_FLAG']

# subdwarf
    if (kwargs.get('subdwarf','') != ''):
        source_db['SUBDWARF_FLAG'] = [str(i).lower() != 'nan' for i in source_db['METALLICITY_CLASS']]
        source_db['SELECT'][source_db['SUBDWARF_FLAG'] == kwargs.get('subdwarf')] += 1
        count+=1.
        del source_db['SUBDWARF_FLAG']

# metallicity class
    if (kwargs.get('subdwarf_class','') != ''):
        source_db['SUBDWARF_FLAG'] = [str(i).lower() == kwargs.get('subdwarf_class').lower() for i in source_db['METALLICITY_CLASS']]
        source_db['SELECT'][source_db['SUBDWARF_FLAG'] == True] += 1
        count+=1.
        del source_db['SUBDWARF_FLAG']

# red - THIS NEEDS TO BE CHANGED
    if (kwargs.get('red','') != ''):
        source_db['RED_FLAG'] = ['red' in str(i).lower() for i in source_db['LIBRARY']]
        source_db['SELECT'][source_db['RED_FLAG'] == kwargs.get('red')] += 1
        count+=1.
        del source_db['RED_FLAG']

# blue - THIS NEEDS TO BE CHANGED
    if (kwargs.get('blue','') != ''):
        source_db['BLUE_FLAG'] = ['blue' in str(i).lower() for i in source_db['LIBRARY']]
        source_db['SELECT'][source_db['BLUE_FLAG'] == kwargs.get('blue')] += 1
        count+=1.
        del source_db['BLUE_FLAG']

# binaries
    if (kwargs.get('binary','') != ''):
        source_db['BINARY_FLAG'] = [str(i).lower() != 'nan' for i in source_db['BINARY']]
        source_db['SELECT'][source_db['BINARY_FLAG'] == kwargs.get('binary')] += 1
        count+=1.
        del source_db['BINARY_FLAG']

# spectral binaries
    if (kwargs.get('sbinary','') != ''):
        source_db['SBINARY_FLAG'] = [str(i).lower() != 'nan' for i in source_db['SBINARY']]
        source_db['SELECT'][source_db['SBINARY_FLAG'] == kwargs.get('sbinary')] += 1
        count+=1.
        del source_db['SBINARY_FLAG']

# companions
    if (kwargs.get('companion','') != ''):
        source_db['COMPANION_FLAG'] = [str(i).lower() != 'nan' for i in source_db['COMPANION_NAME']]
        source_db['SELECT'][source_db['COMPANION_FLAG'] == kwargs.get('companion')] += 1
        count+=1.
        del source_db['COMPANION_FLAG']

# peculiars
    if (kwargs.get('peculiar','') != ''):
#        kwargs['vlm'] = False
        source_db['PECULIAR_FLAG'] = ['p' in str(i).lower() for i in source_db['LIT_TYPE']]
        source_db['SELECT'][source_db['PECULIAR_FLAG'] == kwargs.get('peculiar')] += 1
        count+=1.
        del source_db['PECULIAR_FLAG']

# select source keys
    if count > 0:
        if (logic == 'and'):
            source_db['SELECT'] = numpy.floor(source_db['SELECT']/count)
        elif (logic == 'or'):
            source_db['SELECT'] = numpy.ceil(source_db['SELECT']/count)

#        source_keys = source_db['SOURCE_KEY'][source_db['SELECT']==1]
        source_db = source_db[source_db['SELECT'] == 1]

# quit if there is nothing
    if len(source_db) == 0: 
        if verbose == True: print('\nNo sources found')
        return source_db

#        source_keys = list(source_db[source_db['SELECT']==1]['SOURCE_KEY'])
#        print(source_keys)
# no selection made on sources - choose everything
#    else:
#        source_keys = list(source_db['SOURCE_KEY'])


#    print(count,numpy.max(source_db['SELECT']),len(source_db[:][numpy.where(source_db['SELECT']==1)]),len(source_keys))


# read in spectral database
#    spectral_db = ascii.read(SPLAT_PATH+DB_FOLDER+SPECTRA_DB, delimiter='\t',fill_values='-99.',format='tab')
#    spectral_db = fetchDatabase(SPLAT_PATH+DB_FOLDER+SPECTRA_DB)
#    spectral_db = copy.deepcopy(DB_SPECTRA)

# merge with source_db selected sources
#    print(DB_SPECTRA['DATA_FILE'])
#    spectral_db = source_db.join(copy.deepcopy(DB_SPECTRA),on='SOURCE_KEY',rsuffix='_SP')
    spectral_db = source_db.join(DB_SPECTRA.set_index('SOURCE_KEY'),on='SOURCE_KEY',how='inner',rsuffix='_SP')
    spectral_db.reset_index(drop=True,inplace=True)

#    print(spectral_db['DATA_FILE'])

# having to force dtype here so SELECT remains an integer
    spectral_db['SELECT'] = numpy.zeros(len(spectral_db['DATA_KEY']))
    count = 0.

#    spectral_db['SOURCE_SELECT'] = [x in source_keys for x in spectral_db['SOURCE_KEY']]
#    print(spectral_db['SOURCE_KEY'][numpy.where(spectral_db['SOURCE_SELECT']==True)])

# search by data key
    datakey = kwargs.get('datakey',False)
    datakey = kwargs.get('data_key',datakey)
    if datakey != False:
        if not isinstance(datakey,list):
            datakey = [datakey]
        if isinstance(datakey[0],str):
            datakey = [int(i) for i in datakey]
        for s in datakey:
            spectral_db['SELECT'][spectral_db['DATA_KEY'] == s] += 1
        count+=1.

# exclude by data key
    if kwargs.get('exclude_data_key',False) != False:
        exkey = kwargs['exclude_data_key']
        if len(exkey) > 0:
            if isinstance(exkey,str):
                exkey = [exkey]
            for f in exkey:
                spectral_db['SELECT'][spectral_db['DATA_KEY'] != f] += 1
            count+=1.

# search by filename
    file = kwargs.get('file','')
    file = kwargs.get('filename',file)
    if (file != ''):
        if isinstance(file,str):
            file = [file]
        for f in file:
            spectral_db['SELECT'][spectral_db['DATA_FILE'] == f] += 1
        count+=1.


# exclude by filename
    if kwargs.get('exclude_file',False) != False:
        file = kwargs['exclude_file']
        if len(file) > 0:
            if isinstance(file,str):
                file = [file]
            for f in file:
                spectral_db['SELECT'][spectral_db['DATA_FILE'] != f] += 1
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
        spectral_db['SELECT'][numpy.logical_and(spectral_db['DATEN'] >= date[0],spectral_db['DATEN'] <= date[1])] += 1
        count+=1.
        del spectral_db['DATEN']

# search by S/N range
    if kwargs.get('snr',False) != False:
        snr = kwargs['snr']
        if not isinstance(snr,list):        # one value = minimum S/N
            snr = [float(snr),1.e9]
#        spectral_db['SNRN'] = [float('0'+str(x)) for x in spectral_db['MEDIAN_SNR']]
#        spectral_db['SELECT'][numpy.logical_and(numpy.logical_and(spectral_db['MEDIAN_SNR'] >= snr[0],spectral_db['MEDIAN_SNR'] <= snr[1]),not numpy.isnan(spectral_db['MEDIAN_SNR']))] += 1
        spectral_db['SELECT'][numpy.logical_and(spectral_db['MEDIAN_SNR'] >= snr[0],spectral_db['MEDIAN_SNR'] <= snr[1])] += 1
        count+=1.

# search by instrument
    if kwargs.get('instrument',False) != False:
        instr = checkInstrument(kwargs['instrument'])
        if instr == False and verbose==True: print('Warning: instrument {} is not a valid instrument name'.format(kwargs['instrument']))
        else:
            spectral_db['INST'] = [checkInstrument(i) for i in spectral_db['INSTRUMENT']]
            spectral_db['SELECT'][spectral_db['INST'] == instr] += 1
            count+=1.
            del spectral_db['INST']

# search by reference list
    drefer = kwargs.get('bibcode',False)
    drefer = kwargs.get('reference',drefer)
    drefer = kwargs.get('ref',drefer)
    if drefer != False:
        if isinstance(drefer,str):
            drefer = [drefer]
        for r in drefer:
            spectral_db['SELECT'][spectral_db['DATA_REFERENCE'] == r] += 1
        count+=1.

# search by spex type
    if spt_range != False and spt_type == 'SPEX_TYPE':
        if not isinstance(spt_range,list):        # one value = only this type
            spt_range = [spt_range,spt_range]
        if isinstance(spt_range[0],str):          # convert to numerical spt
            spt_range = [typeToNum(spt_range[0]),typeToNum(spt_range[1])]
        spectral_db['SPTN'] = [typeToNum(x) for x in spectral_db['SPEX_TYPE']]
        spectral_db['SELECT'][numpy.logical_and(spectral_db['SPTN'] >= spt_range[0],spectral_db['SPTN'] <= spt_range[1])] += 1
        count+=1.
        del spectral_db['SPTN']

# select by quality flag
    if kwargs.get('ok',False) != False or kwargs.get('quality','').upper() == 'OK':
        spectral_db['SELECT'][spectral_db['QUALITY_FLAG'] == 'OK'] += 1
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
        spectral_db['SELECT'][spectral_db['PUBLISHED'] != 'Y'] = 0.

#    print(spectral_db['SOURCE_KEY'][numpy.where(spectral_db['SELECT']==1)])
#    print(spectral_db['SOURCE_KEY'][numpy.where(spectral_db['SOURCE_SELECT']==True)])

# no matches
#    print(count,numpy.max(spectral_db['SOURCE_SELECT']),numpy.max(spectral_db['SELECT']))
#    print(len(spectral_db[:][numpy.where(spectral_db['SELECT']==1)]))
#    print(len(spectral_db[:][numpy.where(spectral_db['SOURCE_SELECT']==True)]))
#    print(len(spectral_db[:][numpy.where(numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==True))]))
    db = spectral_db[spectral_db['SELECT']==1]
#    sdb = spectral_db[spectral_db['SOURCE_SELECT']==True]
#    if len(spectral_db[numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==True)]) == 0:
    if len(db) == 0:
        if verbose: print('No spectra in the SPL database match the selection criteria')
        return db
    else:

# merge databases
#        db = join(spectral_db[:][numpy.where(numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==True))],source_db,keys='SOURCE_KEY')

# sort output - default is by designation
        sortkey = kwargs.get('sort','DESIGNATION')
        if sortkey.upper() == 'SNR': sortkey='MEDIAN_SNR'
        if sortkey.upper() in list(db.columns):
            db.sort_values(sortkey.upper(),ascending=(not kwargs.get('reverse',False)),inplace=True)

# select what to return
        if ref != 'all' and ref in list(db.columns):
            db = db[ref]

# reset index
        db.reset_index(drop=True,inplace=True)

# return only first or lucky or all
        if kwargs.get('first',False) == True:
            return db.iloc[[0]]
        elif kwargs.get('lucky',False) == True:
            return db.iloc[[numpy.random.choice(len(db))]]
        else:
            return db



def readSpectrum(verbose=False,*args,**kwargs):
    '''
    .. DOCS: will come back to this one
    '''
# keyword parameters
    folder = kwargs.get('folder','')
    catchSN = kwargs.get('catchSN',True)
    kwargs['model'] = False
    instrument = kwargs.get('instrument','')
    inst = checkInstrument(instrument)
    if inst != False: instrument = inst
#    local = kwargs.get('local',True)
    online = False
    dnldflag = False
#    local = not online
    url = kwargs.get('url',SPLAT_URL+DATA_FOLDER)

# filename
    file = kwargs.get('file','')
    file = kwargs.get('filename',file)
    if (len(args) > 0):
        file = args[0]
    kwargs['filename'] = file
    kwargs['model'] = False

# a filename must be passed
    if file == '':
        raise NameError('\nNo filename passed to readSpectrum')

# first pass: check if file is local
#    if online == False:
    if not os.path.exists(os.path.normpath(file)):
        nfile = folder+os.path.basename(kwargs['filename'])
        if not os.path.exists(os.path.normpath(nfile)):
            if verbose==True: print('Cannot find '+kwargs['filename']+' locally, trying online\n\n')
            online = checkAccess()
            if online == False:
                raise ValueError('\nCannot find files {} or {} locally, and you do not have remote access'.format(file,nfile))
        else:
            file=nfile
            kwargs['filename'] = file

# second pass: download file if necessary
#    online = not local
    if online == True:
        if checkOnline(url+file) == '':
            file = folder+os.path.basename(file)
            if checkOnline(url+file) == '':
                raise ValueError('\nCannot find file '+kwargs['filename']+' on SPLAT website\n\n')
# read in online file
#           file = kwargs['filename']
        try:
            if os.path.exists(os.path.normpath(os.path.basename(kwargs['filename']))):
                os.remove(os.path.normpath(os.path.basename(kwargs['filename'])))
            open(os.path.normpath(os.path.basename(kwargs['filename'])), 'wb').write(requests.get(url+file).content)
            dnldflag = True
        except:
            raise NameError('\nProblem reading in {} from SPLAT website'.format(kwargs['filename']))

# instrument specific reads
    if instrument.upper()=='APOGEE': output = _readAPOGEE(file,**kwargs)
    elif instrument.upper()=='BOSS': output = _readBOSS(file,**kwargs)
    elif instrument.upper()=='LDSS3': output = _readIRAF(file,**kwargs)
    elif instrument.upper()=='FIRE': output = _readFIRE(file,**kwargs)
    elif instrument.upper()=='MAGE': output = _readMAGE(file,**kwargs)
    elif instrument.upper()=='WFC3': output = _readWFC3(file,**kwargs)
    elif instrument.upper()=='KAST-RED' or instrument.upper()=='KAST-BLUE': output = _readKAST(file,**kwargs)
    else:

# determine which type of file
        ftype = file.split('.')[-1]
        zipflag = ''

# gzip compressed file - unzip and then rezip later
        if ftype == 'gz':
            zipflag = ftype
            file = file.replace('.'+ftype,'')
            with open(os.path.normpath(file), 'wb') as f_out, gzip.open(os.path.normpath(file+'.gz'), 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

# bz2 compressed file - unzip and then rezip later
        if ftype == 'bz2':
            zipflag = ftype
            file = file.replace('.'+ftype,'')
            with open(os.path.normpath(file), 'wb') as f_out, bz2.open(os.path.normpath(file+'.bz2'), 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

# fits file
        if (ftype == 'fit' or ftype == 'fits'):
#        df = fits.open(file)
#        with fits.open(file, ignore_missing_end=True) as data:
            with fits.open(os.path.normpath(file),ignore_missing_end=True) as data:
                data.verify('silentfix+ignore')
                if 'NAXIS3' in list(data[0].header.keys()):
                    d = numpy.copy(data[0].data[0,:,:])
                else:
                    d =  numpy.copy(data[0].data)
                header = data[0].header

# ascii file
        else:
            try:
                d = numpy.genfromtxt(os.path.normpath(file), comments='#', unpack=False, \
                    missing_values = ('NaN','nan'), filling_values = (numpy.nan)).transpose()
            except ValueError:
                d = numpy.genfromtxt(os.path.normpath(file), comments=';', unpack=False, \
                     missing_values = ('NaN','nan'), filling_values = (numpy.nan)).transpose()
            header = fits.Header()      # blank header

# delete file if this was an online read
#        if online and not local and os.path.exists(os.path.basename(file)):
#            os.remove(os.path.normpath(os.path.basename(file)))

# remove file if this was a zipped file
        if zipflag != '':
            os.remove(os.path.normpath(file))

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

        output = {'wave': wave, 'flux': flux, 'noise': noise, 'header': header}


# make sure arrays are numpy arrays
    output['wave'] = numpy.array(output['wave'])
    output['flux'] = numpy.array(output['flux'])
    output['noise'] = numpy.array(output['noise'])

# fix places where noise is claimed to be 0
    w = numpy.where(output['noise'] == 0.)
    output['noise'][w] = numpy.nan

# fix nans in flux
#    w = numpy.where(numpy.isnan(flux) == True)
#    flux[w] = 0.

# remove all parts of spectrum that are nans
    w = numpy.where(numpy.logical_and(numpy.isnan(output['wave']) == False,numpy.isnan(output['flux']) == False))
    output['wave'] = output['wave'][w]
    output['flux'] = output['flux'][w]
    output['noise'] = output['noise'][w]


# fix to catch badly formatted files where noise column is S/N
#    print(flux, numpy.median(flux))
    if catchSN == True:
          w = numpy.where(output['flux'] > numpy.nanmedian(output['flux']))
          if (numpy.nanmedian(output['flux'][w]/output['noise'][w]) < 1.):
              output['noise'] = output['flux']/output['noise']
              w = numpy.where(numpy.isnan(output['noise']))
              output['noise'][w] = numpy.nanmedian(output['noise'])

# add in instrument specific information
    if inst != False:
        for k in list(INSTRUMENTS[inst].keys()): output[k] = INSTRUMENTS[inst][k]  

# clean up
    if online==True and dnldflag == True:
        os.remove(os.path.normpath(os.path.basename(kwargs['filename'])))
    if 'wunit' not in list(output.keys()): output['wunit'] = kwargs.get('wunit',DEFAULT_WAVE_UNIT)
    if 'funit' not in list(output.keys()): output['funit'] = kwargs.get('funit',DEFAULT_FLUX_UNIT)
    return output


def _readAPOGEE(file,**kwargs):
    '''
    Reads APOGEE fits files following instructions from XXXX
    Assumes you have an ascap data file with data dimensions of [data],[uncertainty],[]
    '''

# make sure file is there
    if not os.path.exists(file):
        raise NameError('\nCould not find APOGEE file {}'.format(file))

# assess data model to use:
    model = kwargs.get('datamodel','apstar').lower()
    print(kwargs,model)
    if kwargs.get('apstar',False) == True: model='apstar'
    if kwargs.get('apvisit',False) == True: model='apvisit'
    if kwargs.get('aspcap',False) == True: model='aspcap'
    if kwargs.get('ap1d',False) == True: model='ap1d'

    hdulist = fits.open(file)
    header = hdulist[0].header
    output = {'wunit': u.Angstrom,'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom,'header': header}
# apstar data - combined fluxes and individual visits
# https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/TELESCOPE/LOCATION_ID/apStar.html
    if model=='apstar':
        output['wave'] = 10.**(numpy.linspace(hdulist[1].header['CRVAL1'],hdulist[1].header['CRVAL1']+hdulist[1].header['NAXIS1']*hdulist[1].header['CDELT1'],num=hdulist[1].header['NAXIS1'],endpoint=False))*u.Angstrom
        output['flux'] = numpy.array(hdulist[1].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['noise'] = numpy.array(hdulist[2].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['mask'] = numpy.array(hdulist[3].data[0])
        output['sky'] = numpy.array(hdulist[4].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['skynoise'] = numpy.array(hdulist[5].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['telluric'] = numpy.array(hdulist[6].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['telluricnoise'] = numpy.array(hdulist[7].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        if header['NVISITS'] > 1:
            output['flux_visits'] = [numpy.array(hdulist[1].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['noise_visits'] = [numpy.array(hdulist[2].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['mask_visits'] = [numpy.array(hdulist[3].data[i]) for i in range(2,header['NVISITS']+2,1)]
            output['sky_visits'] = [numpy.array(hdulist[4].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['skynoise_visits'] = [numpy.array(hdulist[5].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['telluric_visits'] = [numpy.array(hdulist[6].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['telluricnoise_visits'] = [numpy.array(hdulist[7].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]

# apvisit data - combined fluxes and individual visits
# https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/TELESCOPE/PLATE_ID/MJD5/apVisit.html 
    elif model=='apvisit':
        output['wave'] = numpy.array(hdulist[4].data[0])*u.Angstrom
        output['flux'] = numpy.array(hdulist[1].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['noise'] = numpy.array(hdulist[2].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['mask'] = numpy.array(hdulist[3].data[0])
        output['sky'] = numpy.array(hdulist[5].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['skynoise'] = numpy.array(hdulist[6].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['telluric'] = numpy.array(hdulist[7].data[0])
        output['telluricnoise'] = numpy.array(hdulist[8].data[0])
        if len(hdulist[1].data) > 1:
            output['wave_visits'] = [numpy.array(hdulist[4].data[i])*u.Angstrom for i in range(len(hdulist[1].data))]
            output['flux_visits'] = [numpy.array(hdulist[1].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(len(hdulist[1].data))]
            output['noise_visits'] = [numpy.array(hdulist[2].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(len(hdulist[1].data))]
            output['mask_visits'] = [numpy.array(hdulist[3].data[i]) for i in range(len(hdulist[1].data))]
            output['sky_visits'] = [numpy.array(hdulist[5].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(len(hdulist[1].data))]
            output['skynoise_visits'] = [numpy.array(hdulist[6].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(len(hdulist[1].data))]
            output['telluric_visits'] = [numpy.array(hdulist[7].data[i]) for i in range(len(hdulist[1].data))]
            output['telluricnoise_visits'] = [numpy.array(hdulist[8].data[i]) for i in range(len(hdulist[1].data))]

# aspcap data -  fluxes and individual visits
# https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/LOCATION_ID/aspcapStar.html
    elif model=='aspcap':
        output['wave'] = 10.**(numpy.linspace(hdulist[1].header['CRVAL1'],hdulist[1].header['CRVAL1']+hdulist[1].header['NAXIS1']*hdulist[1].header['CDELT1'],num=hdulist[1].header['NAXIS1'],endpoint=False))*u.Angstrom
        output['flux'] = numpy.array(hdulist[1].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['noise'] = numpy.array(hdulist[2].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['mask'] = numpy.array(hdulist[3].data[0])
        output['sky'] = numpy.array(hdulist[4].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['skynoise'] = numpy.array(hdulist[5].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['telluric'] = numpy.array(hdulist[6].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        output['telluricnoise'] = numpy.array(hdulist[7].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)
        if header['NVISITS'] > 1:
            output['flux_visits'] = [numpy.array(hdulist[1].data[i])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['noise_visits'] = [numpy.array(hdulist[2].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['mask_visits'] = [numpy.array(hdulist[3].data[0]) for i in range(2,header['NVISITS']+2,1)]
            output['sky_visits'] = [numpy.array(hdulist[4].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['skynoise_visits'] = [numpy.array(hdulist[5].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['telluric_visits'] = [numpy.array(hdulist[6].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]
            output['telluricnoise_visits'] = [numpy.array(hdulist[7].data[0])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) for i in range(2,header['NVISITS']+2,1)]

    elif model=='ap1d':
        raise ValueError('\nHave not implemented ap1d yet')

    else:
        raise ValueError('\nNeed to specify which APOGEE data model you are using (apstar, apvisit, aspcap, ap1d)')

    return output


def _readBOSS(file,**kwargs):
    '''
    Reads BOSS fits files following instructions from goo.gl/njCQp5
    '''
    if not os.path.exists(file):
        raise NameError('\nCould not find BOSS file {}'.format(file))

    hdulist=fits.open(file)
    header = hdulist[0].header
    wave = 10.**(hdulist[1].data['loglam'])*u.Angstrom  # log10(wavelength [Å]
    flux = numpy.array(hdulist[1].data['flux'])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)   # coadded calibrated flux in 10^-17 ergs/s/cm2/Å
    ivar = numpy.array(hdulist[1].data['ivar'])   # inverse variance of the flux
    noise = numpy.array([(i**(-0.5)) for i in ivar])*1.e-17 *(u.erg/u.s/u.cm/u.cm/u.Angstrom)
    sky = numpy.array(hdulist[1].data['sky'])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom)     # subtracted sky flux in 10^-17 ergs/s/cm2/Å
    model = numpy.array(hdulist[1].data['model'])*1.e-17*(u.erg/u.s/u.cm/u.cm/u.Angstrom) # pipeline best model fit used for classification and redshift
    
    return {'wave': wave,
          'flux': flux,
          'noise': noise,
          'header': header,
          'sky': sky,
          'model': model,
          'wunit': u.Angstrom,
          'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom}

def _readKAST(file,**kwargs):
    '''
    Reads Carl Melis's KAST fits files
    '''
    if not os.path.exists(file):
        raise NameError('\nCould not find KAST file {}'.format(file))

    hdulist=fits.open(file)
    header = hdulist[0].header
    flux = numpy.array(hdulist[0].data[0,:][0])
    noise = numpy.array(hdulist[0].data[3,:][0])
    if 'CRVAL1' in list(header.keys()) and 'CRDELT1' in list(header.keys()):
        wave = 10.**(numpy.linspace(header['CRVAL1'],header['CRVAL1']+len(flux)*header['CRDELT1'],len(flux)))
    elif 'CRVAL1' in list(header.keys()) and 'CD1_1' in list(header.keys()):
        wave = 10.**(numpy.linspace(header['CRVAL1'],header['CRVAL1']+len(flux)*header['CD1_1'],len(flux)))
    else:
        raise ValueError('\nCould not find appropriate header keywords to make wavelength axis')
    
    return {'wave': wave,
          'flux': flux,
          'noise': noise,
          'header': header,
          'wunit': u.Angstrom,
          'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom}
        

def _readIRAF(file, **kwargs):
    '''
    Reads IRAF data, including data reduced for LDSS3
    Wavelength data is described in the header.
    Does not assume that the data have a noise column
    '''
    if not os.path.exists(file):
        raise NameError('\nCould not find IRAF spectral file {}'.format(file))

    hdu = fits.open(file)
    header = hdu[0].header
    wave = (numpy.linspace(hdu[0].header['CRVAL1'],hdu[0].header['CRVAL1']+hdu[0].header['NAXIS1']*hdu[0].header['CDELT1'],
                                num=hdu[0].header['NAXIS1'],endpoint=False))*u.Angstrom

    flux=hdu[0].data*u.erg/u.s/u.cm/u.cm/u.Angstrom
    try:
        noise = hdu[1].data*u.erg/u.s/u.cm/u.cm/u.Angstrom
    except:
        noise = numpy.array([numpy.nan for x in flux]) #doesn't have noise
    return {'wave': wave,
            'flux': flux,
            'noise': noise ,
            'header': header,
            'wunit': u.Angstrom, 
            'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom}


def _readMAGE(file, **kwargs):
    '''
    Reads MagE spectral files as output from MASE.    
    This data don't have wavelength in the data, but the wavelength is described in the header.
    There may be some data from Mage that have 3 columns that include wavelength, flux, and noise.
    '''

    if not os.path.exists(file):
        raise NameError('\nCould not find MagE file {}'.format(file))
    wunit = u.Angstrom
    funit = u.erg/u.s/u.cm/u.cm/u.Angstrom

    hdulist = fits.open(file)
    header = hdulist[0].header
    if len(hdulist[0].data) == 2:
        wave = (10.**(numpy.linspace(hdulist[0].header['CRVAL1'],hdulist[0].header['CRVAL1']+hdulist[0].header['NAXIS1']*hdulist[0].header['CDELT1'],num=hdulist[0].header['NAXIS1'],endpoint=False)))*wunit
        flux = hdulist[0].data[0]*funit
        noise = hdulist[0].data[1]*funit
    elif len(hdulist[0].data) == 3:
        wave = hdulist[0].data[0]*wunit
        flux = hdulist[0].data[1]*funit
        noise = hdulist[0].data[2]*funit
    else:
        raise ValueError('\nWas expecting 2 or 3 data axes; instead found {}'.format(len(hdulist[0].data)))

    return {'wave': wave, 
            'flux': flux, 
            'noise': noise, 
            'header': header, 
            'wunit': wunit,
            'funit': funit}


def _readFIRE(file,**kwargs):
    '''
    Reads FIRE spectral files as output from FIREHOSE.
    
    '''
    if not os.path.exists(file):
        raise NameError('\nCould not find FIRE file {}'.format(file))
    wunit = u.micron
    funit = u.erg/u.s/u.cm/u.cm/u.micron

    hdulist = fits.open(file)
    if len(hdulist[0].data) != 3:
        raise ValueError('\n_readFIRE can only read data reduced via FIREHOSE, with wave, flux, noise in fits data channel')
    header = hdulist[0].header
    wave = hdulist[0].data[0]*wunit
    flux = hdulist[0].data[1]*funit
    noise = hdulist[0].data[2]*funit
    if kwargs.get('gagne',False) == True or numpy.median(wave.value) > 3.:
        wave = (hdulist[0].data[0]*u.Angstrom).to(wunit)
        flux = (hdulist[0].data[1]*u.erg/u.s/u.cm/u.cm/u.Angstrom).to(funit)
        noise = (hdulist[0].data[2]*u.erg/u.s/u.cm/u.cm/u.Angstrom).to(funit)
    
    return {'wave': wave, 
            'flux': flux, 
            'noise': noise, 
            'header': header, 
            'wunit': wunit,
            'funit': funit} 


def _readWFC3(file,**kwargs):
    '''
    Reads WFC3 spectral files reduced by the Axe software 
    http://www.stsci.edu/hst/wfc3/documents/WFC3_aXe_cookbook.pdf
    Handles both ascii and fits files from the WISP and HST-3D surveys
    '''
    if not os.path.exists(file):
        raise NameError('\nCould not find WFC3 file {}'.format(file))
    #for ascii files, (end with .dat, .ascii or .csv ?)
    if not file.endswith('.fits'):
        #sometimes columns are denoted by 1, 2, etc.
        try:
            data=ascii.read(file)
            wave= data['col1'] 
            flux=data['col2']
            noise=data['col3']
            contam=data['col4']
            output={'wave':wave,'flux':flux,'noise':noise,'contamination': contam,'wunit': u.Angstrom,'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom}
        except KeyError:
            try:
                output=splat.readSpectrum(file)
                output['wunit']=u.Angstrom
                output['funit']=u.erg/u.s/u.cm/u.cm/u.Angstrom
            except:
                raise ValueError('\nCould not read in WFC3 spectrum {}'.format(file))
    #for fits files
    else:
        hdu=fits.open(file)
        header=hdu[1].header
        wave=hdu[1].data['wave']
        flux=hdu[1].data['flux']
        noise=hdu[1].data['error']
        contam=hdu[1].data['contam']
        sens=hdu[1].data['sensitivity']
        output={'wave':wave,
                'flux':flux,
                'noise':noise,
                'contamination': contam,
                'header': header,
                'sensitivity': sens,
                'wunit': u.Angstrom,
                'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom}
   
    return output



#####################################################
###############   CLASSIFICATION   ##################
#####################################################


def classifyByIndex(sp,ref='burgasser',str_flag=True,rnd_flag=False,rem_flag=True,nsamples=100,nloop=5,verbose=False,indices={},sptoffset=0.,coeffs={},method='polynomial', **kwargs):

    '''
    :Purpose: 

    Determine the spectral type and uncertainty for a spectrum based on indices. Makes use of published index-SpT relations from `Reid et al. (2001) <http://adsabs.harvard.edu/abs/2001AJ....121.1710R>`_; `Testi et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...552L.147T>`_; `Allers et al. (2007) <http://adsabs.harvard.edu/abs/2007ApJ...657..511A>`_; and `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_. Returns 2-element tuple containing spectral type (numeric or string) and uncertainty.

    Required Inputs:

    :param sp: Spectrum class object, which should contain wave, flux and noise array elements.

    Optional Inputs:

    :param set: named set of indices to measure and compute spectral type

        - *'burgasser'*: H2O-J, CH4-J, H2O-H, CH4-H, CH4-K from `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_, applicable for types  (default)
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

    str_flag = kwargs.get('string', str_flag)
#    verbose = kwargs.get('verbose', False)
    rnd_flag = kwargs.get('round', rnd_flag)
    rem_flag = kwargs.get('remeasure', rem_flag)
#    nsamples = kwargs.get('nsamples', 100)
#    nloop = kwargs.get('nloop', 5)
    ref = kwargs.get('set',ref)
    kwargs['ref'] = ref
#    if (set.lower() not in allowed_sets):
#        print('\nWarning: index classification method {} not present; returning nan\n\n'.format(set))
#        return numpy.nan, numpy.nan

# Reid et al. (2001, AJ, 121, 1710)
    if (ref.lower() == 'reid'):
        if (rem_flag or len(args) == 0):
            indices = measureIndexSet(sp, **kwargs)
        sptoffset = 20.
        coeffs = { \
            'H2O-A': {'fitunc': 1.18, 'range': [18,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [-32.1, 23.4]}, \
            'H2O-B': {'fitunc': 1.02, 'range': [18,28], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [-24.9, 20.7]}}
        method='polynomial'

# Testi et al. (2001, ApJ, 522, L147)
    elif (ref.lower() == 'testi'):
        if (rem_flag or len(args) == 0):
            indices = measureIndexSet(sp, **kwargs)
        sptoffset = 10.
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
        method='polynomial'

# Burgasser (2007, ApJ, 659, 655) calibration
    elif (ref.lower() == 'burgasser'):
        if (rem_flag or len(args) == 0):
            indices = measureIndexSet(sp, **kwargs)
        sptoffset = 20.
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
        method='polynomial'

# Geballe et al. (2002, ApJ, 564, 466) calibration
    elif ref.lower() == 'geballe':
        if (rem_flag or len(args) == 0):
            i1 = measureIndexSet(sp, ref='geballe')
            i2 = measureIndexSet(sp, ref='martin')
            if sys.version_info.major == 2:
                indices = dict(i1.items() + i2.items())
            else:
                indices = dict(i1.items() | i2.items())
        spttypes = numpy.arange(20.,39.,1.)
        ranges = { \
            'PC3': [[2.4,2.6,20.],[2.6,2.86,21.],[2.85,3.25,22.],[3.25,4.25,23.],[4.25,6,24.]],\
            'Color-d2': [[4.5,5.5,20.],[5.5,6.5,21.],[6.5,7.5,22.],[7.5,10.,23.],[10,17,24.],[17.,23.,25.],[23.,25.,26.]],\
            'H2O-1.2': [[1.5,1.7,30.],[1.7,1.9,31.],[1.9,2.15,32.],[2.15,2.5,33.],[2.5,3.0,34.],[3.0,4.5,35.],[4.5,6.5,36.],[6.5,10.,37.],[10.,15.,38.]],\
            'H2O-1.5': [[1.2,1.27,20.],[1.27,1.35,21.],[1.35,1.43,22.],[1.43,1.5,23.],[1.5,1.55,24.],[1.55,1.6,25.],[1.6,1.65,26.],[1.65,1.7,27.],[1.7,1.8,28.],[1.8,1.95,29.],[1.95,2.2,30.],[2.2,2.5,31.],[2.5,3.0,32.],[3.0,3.5,33.],[3.5,4.5,34.],[4.5,5.5,35.],[5.5,7.,36.],[7.,9.,37.],[9.,12.,38.]],\
            'CH4-1.6': [[1.02,1.07,30.],[1.07,1.15,31.],[1.15,1.3,32.],[1.3,1.5,33.],[1.5,1.8,34.],[1.8,2.5,35.],[2.5,4,36.],[4.,6.,37.],[6.,9.,38.]],\
            'CH4-2.2': [[0.91,0.94,23.],[0.94,0.98,24.],[0.98,1.025,25.],[1.025,1.075,26.],[1.075,1.125,27.],[1.125,1.175,28.],[1.175,1.25,29.],[1.25,1.4,30.],[1.4,1.6,31.],[1.6,1.95,32.],[1.95,2.75,33.],[2.75,3.8,34.],[3.8,5.5,35.],[5.5,8.5,36.],[8.5,12.,37],[12.,18.,38.]],\
        }
        method='ranges'

# Allers et al. (2013, ApJ, 657, 511)
    elif (ref.lower() == 'allers'):
        if (rem_flag or len(args) == 0):
            i1 = measureIndexSet(sp, ref='mclean')
            i2 = measureIndexSet(sp, ref='slesnick')
            i3 = measureIndexSet(sp, ref='allers')
            if sys.version_info.major == 2:
                indices = dict(i1.items() + i2.items() + i3.items())
            else:
                indices = dict(i1.items() | i2.items() | i3.items())
        sptoffset = 10.
        coeffs = { \
            'H2O': {'fitunc': 0.390, 'range': [15,25], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [24.0476, -104.424, 169.388,-83.5437]}, \
            'H2O-1': {'fitunc': 1.097, 'range': [14,25], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [28.5982, -80.7404, 39.3513, 12.1927]}, \
            'H2OD': {'fitunc': 0.757, 'range': [20,28], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [-97.230, 229.884, -202.245, 79.4477]}, \
            'H2O-2': {'fitunc': 0.501, 'range': [14,22], 'spt': 0., 'sptunc': 99., 'mask': 1., \
            'coeff': [37.5013, -97.8144, 55.4580, 10.8822]}}
        method='polynomial'

    else:
        if len(indices.keys()) == 0: 
            print('Error: indices is an empty dictionary')
            return numpy.nan, numpy.nan
        if len(coeff.keys()) == 0: 
            print('Error: coeffs is an empty dictionary')
            return numpy.nan, numpy.nan

# check coefficient index names
    for i in list(coeffs.keys()):
        if i not in list(indices.keys()):
            print('Error: coefficient index {} is not one of the measured indices {}; remeasure indices'.format(i,list(indices.keys())))
            return numpy.nan, numpy.nan

# polynomial method
    if method=='polynomial':
        for index in coeffs.keys():
            coeffs[index]['spt'] = numpy.nan
            coeffs[index]['sptunc'] = numpy.nan
            coeffs[index]['mask'] = 0.
            if numpy.isfinite(indices[index][0]):
                coeffs[index]['spt'] = numpy.polyval(coeffs[index]['coeff'],indices[index][0])+sptoffset
                coeffs[index]['sptunc'] = coeffs[index]['fitunc']
                if (ref.lower() == 'testi'): coeffs[index]['spt'] = (coeffs[index]['spt']-10.)*10.+10.
# noise sim
                if numpy.isfinite(indices[index][1]):
                    vals = numpy.polyval(coeffs[index]['coeff'],numpy.random.normal(indices[index][0],indices[index][1],nsamples))
                    if (ref.lower() == 'testi'): vals = (vals-10.)*10.+10.
                    coeffs[index]['sptunc'] = (numpy.nanstd(vals)**2+coeffs[index]['fitunc']**2)**0.5
# unmask good values
            if coeffs[index]['spt'] >= numpy.nanmin(coeffs[index]['range']) and coeffs[index]['spt'] <= numpy.nanmax(coeffs[index]['range']) and numpy.isfinite(coeffs[index]['spt']): coeffs[index]['mask'] = 1.

# no good values
        mask = [coeffs[index]['mask'] for index in list(coeffs.keys())]
        if numpy.nansum(mask) == 0.:
            if verbose==True: print('\nNone of the indices in set {} returned viable values\n'.format(set))
            return numpy.nan, numpy.nan

# computed weighted mean with rejection, iterating to deal with indices outside ranges      
        for i in numpy.arange(nloop):
            wts = numpy.array([coeffs[index]['mask']/coeffs[index]['sptunc']**2 for index in list(coeffs.keys())])
            vals = numpy.array([coeffs[index]['mask']*coeffs[index]['spt']/coeffs[index]['sptunc']**2 for index in list(coeffs.keys())])
            w = numpy.where(numpy.isfinite(wts+vals))
            if len(w) == 0:
                if verbose==True: print('\nNone of the indices in set {} returned viable values\n'.format(set))
                return numpy.nan, numpy.nan
            sptn = numpy.nansum(vals[w])/numpy.nansum(wts[w])
            sptn_e = 1./numpy.nansum(wts[w])**0.5
            for index in coeffs.keys():
                if sptn < numpy.nanmin(coeffs[index]['range']) or sptn > numpy.nanmax(coeffs[index]['range']): coeffs[index]['mask'] = 0.

# report individual subtypes
        if verbose == True:
            for i in coeffs.keys():
                flg = '*'
                if coeffs[i]['mask'] == 0.: flg = ''
                print('{}{} = {:.3f}+/-{:.3f} = SpT = {}+/-{:.1f}'.format(flg,i,indices[i][0],indices[i][1],typeToNum(coeffs[i]['spt']),coeffs[i]['sptunc']))

# ranges method - NOT CURRENTLY CONSIDERING UNCERTAINTY
    if method=='ranges':
        spts = []
        spts_unc = []
        for index in list(ranges.keys()):
            spts.append(numpy.nan)
            spts_unc.append(numpy.nan)
            if indices[index][1] > 0.:
                for r in ranges[index]: 
                    if r[0] < indices[index][0] <= r[1]: spts[-1] = r[-1]
# PLACEHOLDER FOR UNCERTAINTY INCLUSION
            if verbose == True:
                flg = ''
                if numpy.isnan(spts[-1]): flg='*'
                print('{}{}: {:.3f}+/-{:.3f} => SpT = {}'.format(flg,index,indices[index][0],indices[index][1],typeToNum(spts[-1])))

        spts = numpy.array(spts)
        sptn = numpy.nanmean(spts)
        sptn_e = numpy.nanstd(spts)

# round off to nearest 0.5 subtypes if desired
    if (rnd_flag): sptn = 0.5*numpy.around(sptn*2.)

# change to string if desired
    if (str_flag): spt = typeToNum(sptn,uncertainty=sptn_e)
    else: spt = sptn

    if kwargs.get('allmeasures',False) == True and method == 'polynomial':
        output = {}
        for k in coeffs.keys():
            output[k] = {'spt': coeffs[k]['spt'], 'spt_e': coeffs[k]['sptunc'], 'index': indices[k][0], 'index_e': indices[k][0]}
        output['result'] = (spt,sptn_e)
        return output
    else:
        return spt, sptn_e



def classifyByStandard(sp, std_class='dwarf', *args, **kwargs):
    '''
    :Purpose: 
        Determine the spectral type and uncertainty for a
        spectrum by direct comparison to defined spectral standards.  
        Dwarf standards span M0-T9 and include the standards listed in
        `Burgasser et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...637.1067B>`_, `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_ and `Cushing et al. (2011) <http://adsabs.harvard.edu/abs/2011ApJ...743...50C>`_. 
        Comparison to subdwarf and extreme subdwarf standards may also be done.
        Returns the best match or an F-test weighted mean and uncertainty. There is an option
        to follow the procedure of `Kirkpatrick et al. (2010)
        <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_, fitting only in
        the 0.9-1.4 micron region.

    :Required Parameters: 

        * :param sp: Spectrum class object, which should contain wave, flux and noise array elements (required)

    :Optional Parameters: 

        * :param sptrange: set to the spectral type range over which comparisons should be made, can be a two-element array of strings or numbers (optional, default = ['M0','T9'])
        * :param statistic: string defining which statistic to use in comparison; available options are:

                - *'chisqr'*: (DEFAULT) compare by computing chi squared value (requires spectra with noise values)
                - *'stddev'*: compare by computing standard deviation
                - *'stddev_norm'*: compare by computing normalized standard deviation
                - *'absdev'*: compare by computing absolute deviation

        * :param std_class: the type of standard to compare to; allowed options are 'dwarf' (default), 'subdwarf', 'sd', 'dsd', 'esd', 'lowg', 'vlg', and 'intg'. These can also be set using the following keywords

            - :param dwarf: set to True to compare to dwarf standards
            - :param sd: set to True to compare to subdwarf standards (`subdwarf` does the same)
            - :param dsd: set to True to compare to dwarf/subdwarf standards 
            - :param esd: set to True to compare to extreme subdwarf standards 
            - :param vlg: set to True to compare to very low gravity standards (`lowg` does the same)
            - :param intg: set to True to compare to intermediate gravity standards

        * :param all: compare to standards across all metallicity and gravity types
        * :param method: set to ``'kirkpatrick'`` to follow the `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_ method, fitting only to the 0.9-1.4 micron band (optional, default = '')
        * :param best: set to True to return the best fit standard type (optional, default = True)
        * :param average: set to True to return an chi-square weighted type only (optional, default = False)
        * :param compareto: set to the single standard (string or number) you want to compare to (optional, default = None)
        * :param plot: set to True to generate a plot comparing best fit template to source; can also set keywords associated with plotSpectrum_ routine (optional, default = False)
        * :param string: set to True to return spectral type as a string (optional, default = True)
        * :param return_standard: set to True to return a Spectrum class of the standard properly scaled (optional, default = False)
        * :param return_statistic: set to True to return a best match spectral type and statistic (instead of uncertainty) (optional, default = False)
        * :param verbose: set to True to give extra feedback (optional, default = False)

    Users can also set keyword parameters defined in plotSpectrum_ and compareSpectra_ routine.

    .. _compareSpectra : api.html#splat.core.compareSpectra
    .. _plotSpectrum : api.html#splat.plot.plotSpectrum

    :Output: 

        A tuple listing the best match standard and uncertainty based on F-test weighting and systematic uncertainty of 0.5 subtypes

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
    allowed_classes = ['dwarf','subdwarf','sd','esd','dsd','lowg','vlg','intg','all']
    for a in allowed_classes:
        if kwargs.get(a,False) == True: std_class = a
    std_class = std_class.lower()
    if std_class not in allowed_classes:
        if verbose==True: print('\nStandard class {} unknown; defaulting to dwarf'.format(std_class))
        std_class = 'dwarf'

    if std_class == 'dwarf':
        initiateStandards()
        stds = STDS_DWARF_SPEX
        subclass = ''
        stdtype = 'Dwarf'
        if verbose==True: print('Using dwarf standards')
    elif std_class == 'sd' or std_class == 'subdwarf':
        initiateStandards(sd=True)
        stds = STDS_SD_SPEX
        subclass = 'sd'
        stdtype = 'Subdwarf'
        if verbose==True: print('Using subdwarf standards')
    elif std_class == 'dsd':
        initiateStandards(dsd=True)
        stds = STDS_DSD_SPEX
        subclass = ''
        stdtype = 'Mild subdwarf'
        if verbose == True: print('Using dwarf standards')
    elif std_class == 'esd':
        initiateStandards(esd=True)
        stds = STDS_ESD_SPEX
        subclass = 'esd'
        stdtype = 'Extreme Subdwarf'
        if verbose == True: print('Using extreme subdwarf standards')
    elif std_class == 'vlg' or std_class == 'lowg':
        initiateStandards(vlg=True)
        stds = STDS_VLG_SPEX
        subclass = ''
        stdtype = 'Very Low Gravity'
        if verbose == True: print('Using very low gravity standards')
    elif std_class == 'intg':
        initiateStandards(intg=True)
        stds = STDS_INTG_SPEX
        subclass = ''
        stdtype = 'Intermediate Gravity'
        if verbose == True: print('Using intermediate low gravity standards')
    elif std_class == 'all':
        initiateStandards()
        initiateStandards(sd=True)
        initiateStandards(dsd=True)
        initiateStandards(esd=True)
        initiateStandards(vlg=True)
        initiateStandards(intg=True)
        stds = STDS_DWARF_SPEX.copy()
        stds.update(STDS_VLG_SPEX)
        stds.update(STDS_INTG_SPEX)
        stds.update(STDS_SD_SPEX)
        stds.update(STDS_DSD_SPEX)
        stds.update(STDS_ESD_SPEX)
        subclass = ''
        stdtype = 'Mixed'
        if verbose == True: print('Using all of the standards')
    else:
        raise ValueError('\nUnknown class type {}'.format(std_class))

# select desired spectral range
    std_spt = numpy.array(list(stds.keys()))
    std_sptn = numpy.array([typeToNum(s) for s in list(stds.keys())])
    spt_sample = std_spt[numpy.where(numpy.logical_and(std_sptn >= sptrange[0],std_sptn<=sptrange[1]))]
    spt_sample_n = std_sptn[numpy.where(numpy.logical_and(std_sptn >= sptrange[0],std_sptn<=sptrange[1]))]
#    spt_allowed = numpy.array([typeToNum(s) for s in stds.keys()])
#    spt_sample = spt_allowed[numpy.where(spt_allowed >= sptrange[0])]
#    spt_sample = spt_sample[numpy.where(spt_sample <= sptrange[1])]

# determine comparison range based on method
    if (kwargs.get('method','').lower() == 'kirkpatrick'):
        fit_ranges = [[0.9,1.4]]         # as prescribed in Kirkpatrick et al. 2010, ApJS,
    elif (kwargs.get('method','').lower() == ''):
        fit_ranges = [[0.7,2.45]]       # by default, compare whole spectrum
    else:
        print('\nWarning: do not recognize method = {}'.format(kwargs['method']))
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
#        chisq,scale = compareSpectra(sp,stds[typeToNum(t,subclass=subclass)],fit_ranges=fit_ranges,statistic=statistic,novar2=True)
        chisq,scale = compareSpectra(sp,stds[t],fit_ranges=fit_ranges,statistic=statistic,novar2=True)
        stat.append(chisq)
        sspt.append(t)
        if (verbose):
            print('Type {}: statistic = {}, scale = {}'.format(t, chisq, scale))

# list of sorted standard files and spectral types
    sorted_stdsptnum = [x for (y,x) in sorted(zip(stat,sspt))]

# select either best match or an ftest-weighted average
    if (best_flag or len(stat) == 1):
        spt = sorted_stdsptnum[0]
        sptn = typeToNum(spt)
        sptn_e = unc_sys
    else:
        ssptn = [typeToNum(s) for s in sspt]
        try:
            st = stat.value
        except:
            st = stat
        if numpy.isnan(numpy.median(sp.noise)):
            mean,var = weightedMeanVar(ssptn,st)
        else:
            mean,var = weightedMeanVar(ssptn,st,method='ftest',dof=sp.dof)
        if (var**0.5 < 1.):
            sptn = numpy.round(mean*2)*0.5
        else:
            sptn = numpy.round(mean)
        sptn_e = (unc_sys**2+var)**0.5
        spt = typeToNum(sptn,uncertainty=sptn_e,subclass=subclass)

# string or not?
    if (kwargs.get('string', True) == True):
        output_spt = str(spt)
    else:
        output_spt = sptn

    if verbose:
        print('\nBest match to {} standard'.format(sorted_stdsptnum[0]))
        print('Best spectral type = {}+/-{}'.format(output_spt,sptn_e))

# plot spectrum compared to best spectrum
    if (kwargs.get('plot',False) != False):
#        spstd = Spectrum(file=sorted_stdfiles[0])
#        print(typeToNum(sorted_stdsptnum[0],subclass=subclass))
        spstd = stds[sorted_stdsptnum[0]]
#        getStandard(typeToNum(sorted_stdsptnum[0],subclass=subclass))
        chisq,scale = compareSpectra(sp,spstd,fit_ranges=fit_ranges,statistic=statistic)
        spstd.scale(scale)
        if kwargs.get('colors',False) == False:
            kwargs['colors'] = ['k','r','b']
        if kwargs.get('labels',False) == False:
            kwargs['labels'] = [sp.name,'{} Standard'.format(sorted_stdsptnum[0]),'Difference']
        from .plot import plotSpectrum
        if kwargs.get('difference',True):
            kwargs['labels'].append('Difference')
            pl = plotSpectrum(sp,spstd,sp-spstd,**kwargs)
        else:
            pl = plotSpectrum(sp,spstd,**kwargs)

    if kwargs.get('return_standard',False) == True: 
        spstd = stds[sorted_stdsptnum]
        chisq,scale = compareSpectra(sp,spstd,fit_ranges=fit_ranges,statistic=statistic)
        spstd.scale(scale)
        return spstd
    elif kwargs.get('return_statistic',False) == True: 
        return output_spt, numpy.nanmin(stat)
    else:
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

    .. _classifyByStandard : api.html#splat.core.classifyByStandard
    .. _searchLibrary : api.html#splat.core.searchLibrary
    .. _plotSpectrum : api.html#splat.plot.plotSpectrum

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

    if verbose == True:
        print('\nComparing to {} templates\n'.format(len(files)))
        if len(files) > 100:
            print('This may take some time!\n\n'.format(len(files)))

# do comparison
    stat = []
    scl = []
    for i,d in enumerate(dkey):

# INSERT TRY STATEMNT HERE?

        s = Spectrum(d)

        mkwargs = copy.deepcopy(kwargs)
        if 'plot' in (mkwargs.keys()): mkwargs['plot'] = False

        stt,scale = compareSpectra(sp,s,fit_ranges=fit_ranges,statistic=statistic,novar2=True,**mkwargs)
        stat.append(stt)
        scl.append(scale)
        if verbose == True: print(s)

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
        s = Spectrum(sorted_dkey[0])
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

    if verbose == True:
        s = Spectrum(sorted_dkey[0])
        print('\nBest match = {} with spectral type {}'.format(s,typeToNum(sorted_spt[0])))
        print('Mean spectral type = {}+/-{}'.format(output_spt,sptn_e))

# return dictionary of results
    return {'result': (output_spt,sptn_e), \
        'statistic': sorted(stat)[0:nbest], 'spt': sorted_spt[0:nbest], \
        'scale': sorted_scale[0:nbest], \
        'spectra': [Spectrum(d) for d in sorted_dkey[0:nbest]]}



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



def compareSpectra(s1, s2, statistic='chisqr',scale=True, novar2=True, plot=False, verbose=False, **kwargs):
    '''
    :Purpose: 

        Compare two spectra against each other using a pre-selected statistic. 
        Returns the value of the desired statistic as well as the optimal scale factor. 

    :Required Parameters: 

        :param sp1: First spectrum class object, which sets the wavelength scale
        :param sp2: Second spectrum class object, interpolated onto the wavelength scale of sp1

    :Optional Parameters: 

        :param: statistic = 'chisqr': string defining which statistic to use in comparison; available options are (also stat):

            - *'chisqr'* or *'chi'*: compare by computing chi squared value (requires spectra with noise values)
            - *'stddev'*: compare by computing standard deviation
            - *'stddev_norm'*: compare by computing normalized standard deviation
            - *'absdev'*: compare by computing absolute deviation

        :param: scale = True: If True, finds the best scale factor to minimize the statistic
        :param: fit_ranges = [range of first spectrum]: 2-element array or nested array of 2-element arrays specifying the wavelength ranges to be used for the fit, assumed to be measured in microns; this is effectively the opposite of mask_ranges (also fit_range, fitrange, fitrng, comprange, comprng)
        :param: mask = numpy.zeros(): Array specifiying which wavelengths to mask; must be an array with length equal to the wavelength scale of ``sp1`` with only 0 (OK) or 1 (mask).
        :param: mask_ranges = None: Multi-vector array setting wavelength boundaries for masking data, assumed to be in microns
        :param: mask_telluric = False: Set to True to mask pre-defined telluric absorption regions
        :param: mask_standard = False: Like ``mask_telluric``, with a slightly tighter cut of 0.80-2.35 micron
        :param: weights = numpy.ones(): Array specifying the weights for individual wavelengths; must be an array with length equal to the wavelength scale of ``sp1``; need not be normalized
        :param: novar2 = True: Set to True to compute statistic without considering variance of ``sp2``
        :param: plot = False: Set to True to plot ``sp1`` with scaled ``sp2`` and difference spectrum overlaid
        :param: verbose = False: Set to True to report things as you're going along

    :Output: 
        statistic and optimal scale factor for the comparison

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

    sp1 = copy.deepcopy(s1)
    sp2 = copy.deepcopy(s2)
    
# make sure spectra are on the same wavelength and flux unit scales
    sp2.toWaveUnit(sp1.wave.unit)
    sp2.toFluxUnit(sp1.flux.unit)

    fit_ranges = kwargs.get('fit_ranges',[[numpy.nanmin(sp1.wave),numpy.nanmax(sp1.wave)]])
    fit_ranges = kwargs.get('fit_range',fit_ranges)
    fit_ranges = kwargs.get('fitrange',fit_ranges)
    fit_ranges = kwargs.get('fitrng',fit_ranges)
    fit_ranges = kwargs.get('comprange',fit_ranges)
    fit_ranges = kwargs.get('comprng',fit_ranges)
#    mask_ranges = kwargs.get('mask_ranges',[])
#    mask_standard = kwargs.get('mask_standard',False)
#    mask_telluric = kwargs.get('mask_telluric',mask_standard)
    var_flag = novar2
    if numpy.isnan(numpy.max(sp1.variance.value)) == True: var_flag = False
    statistic = kwargs.get('stat',statistic)
    minreturn = 1.e-60
    scale_factor = 1.

# create interpolation function for second spectrum
    f = interp1d(sp2.wave.value,sp2.flux.value,bounds_error=False,fill_value=0.)
    if var_flag:
        v = interp1d(sp2.wave.value,[numpy.nan for s in sp2.wave],bounds_error=False,fill_value=numpy.nan)
    else:
        v = interp1d(sp2.wave.value,sp2.variance.value,bounds_error=False,fill_value=numpy.nan)
# total variance - funny form to cover for nans
    vtot = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)],axis=0)

# manage fit ranges and generate fit mask
    if len(fit_ranges) == 0: fit_ranges = [[numpy.nanmin(sp1.wave),numpy.nanmax(sp1.wave)]]
    if len(fit_ranges) == 2 and not isinstance(fit_ranges[0],list): fit_ranges = [fit_ranges]
    for i,m in enumerate(fit_ranges):
        if not isUnit(m):
            fit_ranges[i] = (m*u.micron).to(sp1.wave.unit)
    fit_mask = kwargs.get('fit_mask',1.-generateMask(sp1.wave,mask_ranges=fit_ranges))

# generate masking array and combine with fit mask
    reject_mask = numpy.array(kwargs.get('mask',generateMask(sp1.wave,**kwargs)))
# mask flux < 0
    reject_mask[numpy.where(numpy.logical_or(sp1.flux < 0,f(sp1.wave) < 0))] = 1
    mask = numpy.clip(fit_mask+reject_mask,0,1)

# set the weights
    weights = kwargs.get('weights',numpy.ones(len(sp1.wave)))
    weights = weights*(1.-mask)

# comparison statistics
# switch to standard deviation if no uncertainty
    if numpy.isnan(numpy.nanmax(vtot)):
        statistic = 'stddev'
        if verbose==True:
            print('No uncertainties provided; using the {} statistic by default'.format(statistic))
    else:
        if verbose==True:
            print('Comparing spectra using the {} statistic'.format(statistic))

# chi^2
    if (statistic == 'chisqr' or statistic == 'chisq' or statistic == 'chi'):
# compute scale factor
        if scale == True:
            scale_factor = numpy.nansum(weights*sp1.flux.value*f(sp1.wave.value)/vtot)/ \
                numpy.nansum(weights*f(sp1.wave.value)*f(sp1.wave.value)/vtot)
# correct variance
        vtot = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor**2)],axis=0)
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave.value)*scale_factor)**2/vtot)
        unit = sp1.funit/sp1.funit

# normalized standard deviation
    elif (statistic == 'stddev_norm' or statistic == 'stdev_norm'):
# compute scale factor
        if scale == True:
            scale_factor = numpy.nansum(weights*sp1.flux.value)/ \
                numpy.nansum(weights*f(sp1.wave.value))
# correct variance
        vtot = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor**2)],axis=0)
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave.value)*scale_factor)**2)/ \
            numpy.median(sp1.flux.value)**2
        unit = sp1.funit/sp1.funit

# standard deviation
    elif (statistic == 'stddev' or statistic == 'stdev'):
# compute scale factor
        if scale == True:
            scale_factor = numpy.nansum(weights*sp1.flux.value*f(sp1.wave.value))/ \
                numpy.nansum(weights*f(sp1.wave.value)*f(sp1.wave.value))
# correct variance
        vtot = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor**2)],axis=0)
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave.value)*scale_factor)**2)
        unit = sp1.funit**2

# absolute deviation
    elif (statistic == 'absdev'):
# compute scale factor
        if scale == True:
            scale_factor = numpy.nansum(weights*sp1.flux.value)/ \
                numpy.nansum(weights*f(sp1.wave.value))
# correct variance
        vtot = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor**2)],axis=0)
        stat = numpy.nansum(weights*abs(sp1.flux.value-f(sp1.wave.value)*scale_factor))
        unit = sp1.funit

# error
    else:
        print('Error: statistic {} for compareSpectra not available'.format(statistic))
        return numpy.nan, numpy.nan

# plot spectrum compared to best spectrum
    if plot == True:
        spcomp = sp2.copy()
        spcomp.scale(scale_factor)
        kwargs['colors'] = kwargs.get('colors',['k','r','b'])
        kwargs['title'] = kwargs.get('title',sp1.name+' vs '+sp2.name)
        from .plot import plotSpectrum
        plotSpectrum(sp1,spcomp,sp1-spcomp,labels=[sp1.name,sp2.name,'{} = {}'.format(statistic,stat)],**kwargs)


    return numpy.nanmax([stat,minreturn])*unit, scale_factor


def generateMask(wv,mask=[],mask_range=[-99.,-99.],mask_telluric=False,mask_standard=False,**kwargs):
    '''
    :Purpose: Generates a mask array based on wavelength vector and optional inputs on what to mask.

    :Output: A mask array, where 0 = OK and 1 = ignore

    :Example:
    '''

# parameter check
    wave = copy.deepcopy(wv)
    if isinstance(wv,splat.Spectrum): wave = wv.wave
    if not isUnit(wv): wave = wave*DEFAULT_WAVE_UNIT
    if not isinstance(wave.value,list) and not isinstance(wave.value,numpy.ndarray):
        raise ValueError('\nInput parameter should be an array of wavelengths; you passed {}'.format(wv))

# generate initial mask 
    if len(mask) == 0:
        mask = numpy.zeros(len(wave))
    mask = numpy.array(mask)
    if len(mask) != len(wave):
        raise ValueError('\nInitial mask of length {} is not the same as wave array of len {}'.format(len(mask),len(wave)))


# a standard masking
    mask_ranges = kwargs.get('mask_ranges',[mask_range])
    if mask_standard == True:
        mask_telluric = True
        mask_ranges.append([0.,0.8]*DEFAULT_WAVE_UNIT)        # standard short cut

# mask telluric bands
    if mask_telluric == True:
#        mask_ranges.append([0.,0.65]*DEFAULT_WAVE_UNIT)        # meant to clear out short wavelengths
        mask_ranges.append([1.35,1.42]*DEFAULT_WAVE_UNIT)
        mask_ranges.append([1.8,1.92]*DEFAULT_WAVE_UNIT)
        mask_ranges.append([2.45,99.]*DEFAULT_WAVE_UNIT)        # meant to clear out long wavelengths

# make sure quantities are all correct
    mask_ranges_apply = []
    for i,m in enumerate(mask_ranges):
        if not isUnit(m): m = m*wave.unit
        m.to(wave.unit)
        mask_ranges_apply.append(m)

# apply to mask
    for m in mask_ranges_apply:
        mask[numpy.where(numpy.logical_and(wave.value >= numpy.nanmin(m.value),wave.value <= numpy.nanmax(m.value)))]= 1

    return mask



def measureEW(sp,line,width=[0.,0.],continuum=[0.,0.],plot=False,file='',output_unit=u.Angstrom,nsamp=100,verbose=True,recenter=True,absorption=True,name=''):

# input checks
    if isUnit(line):
        line_center = line.to(sp.wave.unit).value
    else: line_center = copy.deepcopy(line)
    if numpy.nanmin(sp.wave.value) > line_center or numpy.nanmax(sp.wave.value) < line_center:
        raise ValueError('\nLine {} is outside spectral data limits of {} to {}'.format(line_center*sp.wave.unit,numpy.nanmin(sp.wave),numpy.nanmax(sp.wave)))
    if isinstance(width,int) or isinstance(width,float):
        line_width = [numpy.abs(width),numpy.abs(width)]
    else: line_width = copy.deepcopy(width) 
    if line_width[0] == 0.:
        ic = numpy.nanargmin(numpy.array([numpy.abs(a-line_center) for a in sp.wave.value]))
        width = numpy.nanmax([line_center/6.e2,numpy.array(sp.wave.value)[ic+5]-line_center])
        line_width = [numpy.abs(width),numpy.abs(width)]
    if isinstance(continuum,int) or isinstance(continuum,float):
        continuum_width = [numpy.abs(continuum),2.*numpy.abs(continuum)]
    else: continuum_width = copy.deepcopy(continuum) 
    if continuum_width[0] == 0.:
        continuum_width = [line_width[0]*2.,line_width[0]*3.]
        
# preset fail condition
    ew = numpy.nan
    ew_unc = numpy.nan
    line_center_measure = numpy.nan
    line_center_measure_unc = numpy.nan
    rv = numpy.nan
    rv_unc = numpy.nan
    
# first compute value
    if numpy.nanmin(sp.wave.value) <= line_center-continuum_width[1] and numpy.nanmax(sp.wave.value) >= line_center+continuum_width[1]:

# refine line centering
        line_center_measure = line_center
        if recenter == True:    
            for i in range(5):
                wc = numpy.where(numpy.logical_and(sp.wave.value >= line_center_measure-line_width[0],sp.wave.value <= line_center_measure+line_width[1]))
                wv = numpy.array(sp.wave.value[wc])
                fl = numpy.array(sp.flux.value[wc])
                if absorption == True: line_center_measure = wv[numpy.nanargmin(fl)]
                else: line_center_measure = wv[numpy.nanargmax(fl)]
        rv = ((line_center_measure-line_center)/line_center)*const.c.to(u.km/u.s)

        w = numpy.where(numpy.logical_and(sp.wave.value >= line_center_measure-continuum_width[1],sp.wave.value <= line_center_measure+continuum_width[1]))
        if len(w[0]) > 0:
            f = interp1d(sp.wave.value[w],sp.flux.value[w],bounds_error=False,fill_value=0.)
            wline = numpy.linspace(line_center_measure-line_width[0],line_center_measure+line_width[1],nsamp)
# either fit over a range or fit specific wavelengths
            if len(continuum_width) == 2:
                wcont = numpy.append(numpy.linspace(line_center_measure-continuum_width[1],line_center_measure-continuum_width[0],nsamp),numpy.linspace(line_center_measure+continuum_width[0],line_center_measure+continuum_width[1],nsamp))        
            else:
                wcont = continuum_width
            fline = f(wline)
            fcont = f(wcont)
            pcont = numpy.poly1d(numpy.polyfit(wcont,fcont,1))
            fcontfit = pcont(wline)
#            print(wline,fline)
#            print(wcont,fcont)
            ew = (trapz((numpy.ones(len(wline))-(fline/fcontfit)), wline)*sp.wave.unit).to(output_unit)
            if plot == True:
                plt.clf()
                plt.plot(sp.wave,sp.flux,'k-')
                plt.plot(wline,fline,'r-')
                plt.plot(wline,fcontfit,'b-')
                plt.plot([line_center_measure,line_center_measure],[numpy.nanmin(sp.flux.value[w]),numpy.nanmax(sp.flux.value[w])],'k--')
                plt.xlim([numpy.nanmin(wcont),numpy.nanmax(wcont)])
                plt.ylim([numpy.nanmin(sp.flux.value[w]),numpy.nanmax(sp.flux.value[w])])
                plt.title(name+' {:.4f}'.format(line_center_measure))
                if file != '': plt.savefig(file)
            
            if sp.noise.value[0] != numpy.nan:
# default uncertainty = 2 x median uncertainty
                ew_unc = (2.*numpy.nanmedian(sp.noise.value[w])/pcont(line_center_measure)*(numpy.nanmax(wline)-numpy.nanmin(wline))*sp.wave.unit).to(output_unit)
# MC for errors
                ews = []
                lns = []
                rvs = []
                spvar = copy.deepcopy(sp)
                for i in range(nsamp):
                    spvar.flux = numpy.random.normal(sp.flux.value,sp.noise.value)*sp.flux.unit
                    line_center_measure_var = line_center
                    if recenter == True:    
                        for i in range(5):
                            wc = numpy.where(numpy.logical_and(sp.wave.value >= line_center_measure_var-line_width[0],sp.wave.value <= line_center_measure_var+line_width[1]))
                            wv = numpy.array(sp.wave.value[wc])
                            fl = numpy.array(spvar.flux.value[wc])
                            if absorption == True: line_center_measure_var = wv[numpy.nanargmin(fl)]
                            else: line_center_measure_var = wv[numpy.nanargmax(fl)]
                    lns.append(line_center_measure_var)
                    rvs.append((((line_center_measure_var-line_center)/line_center)*const.c.to(u.km/u.s)).value)

                    f = interp1d(sp.wave.value[w],spvar.flux.value[w],bounds_error=False,fill_value=0.)
                    wline = numpy.linspace(line_center_measure_var-line_width[0],line_center_measure_var+line_width[1],nsamp)
                    if len(continuum_width) == 2:
                        wcont = numpy.append(numpy.linspace(line_center_measure_var-continuum_width[1],line_center_measure_var-continuum_width[0],nsamp),numpy.linspace(line_center_measure_var+continuum_width[0],line_center_measure_var+continuum_width[1],nsamp))        
                    else:
                        wcont = continuum_width
                    fline = f(wline)
                    fcont = f(wcont)
                    pcont = numpy.poly1d(numpy.polyfit(wcont,fcont,1))
                    fcontfit = pcont(wline)
                    ews.append(((trapz((numpy.ones(len(wline))-(fline/fcontfit)), wline)*sp.wave.unit).to(output_unit)).value)
                rv_unc = numpy.std(rvs)*u.km/u.s
                ew_unc = numpy.sqrt(ew_unc**2+(numpy.std(ews)*output_unit)**2)
                line_center_measure_unc = numpy.std(lns)
                            
    return {'ew': ew, 
            'ew_unc': ew_unc,
            'line_center': line_center_measure*sp.wave.unit,
            'line_center_unc': line_center_measure_unc*sp.wave.unit,
            'rv': rv,
            'rv_unc': rv_unc}


def measureEWElement(sp,element,wave_range=[0.,0.],getNist=False,**kwargs):
    if wave_range[0] == 0.: wave_range = [numpy.nanmin(sp.wave),numpy.nanmax(sp.wave)]
    if not isUnit(wave_range[0]): wave_range = [w*sp.wave.unit for w in wave_range]  

    if getNist==True:
        from .database import queryNist
        t = queryNist(element,wave_range,**kwargs)
        if len(t) == 0: return {}
        lines = [(a*u.Angstrom).to(DEFAULT_WAVE_UNIT) for a in list(t['Observed'])]

    else:
        el = element.lower().strip()
        if el == 'na i':
            lines=[0.81832556,0.81947905,0.81948237,0.864992,0.865089,1.138145,1.140378,2.206242,2.208969]
        elif el == 'k i':
            lines=[1.1690219,1.1769637,1.1772838,1.2432274,1.2522141]
        elif el == 'cs i':
            lines = [0.85211316530,0.894347423876,1.3588293,1.46949087]
        elif el == 'rb i':
            lines = [0.8271410,0.8868512,1.475241]
        elif el == 'mg i':
            lines = [1.1828171,1.2083649,1.4877608,1.5024997,1.5040246,1.5740706,1.7108631]
        elif el == 'ca i':
            lines = [0.941697,1.977679,2.261410,2.263110,2.265741]
        elif el == 'ca ii':
            lines = [0.849802,0.854209,0.866214]
        else:
            print('\nHave not curated lines set for {}'.format(element))
            return {}
        lines = [l*DEFAULT_WAVE_UNIT for l in lines]
            
# measure lines and output to dictionary        
    output = {}
    for l in lines:
#        try:
        res = measureEW(sp,l,**kwargs)
        output['{:.4f}'.format(l.value)] = res
#        except:
#            pass
    return output
        


def measureEWSet(sp,ref='rojas',*args,**kwargs):
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
    ref = kwargs.get('set',ref)

# determine combine method
    if ('rojas' in ref.lower()):
        reference = 'EW measures from Rojas-Ayala et al. (2012)'
        names = ['Na I 2.206/2.209','Ca I 2.26']
        ews = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        tmp = measureEW(sp,2.206,width=[2.2020, 2.2120],continuum=[2.1965, 2.2125, 2.2175],recenter=False,**kwargs)
        ews[0] = tmp['ew']
        errs[0] = tmp['ew_unc']
        tmp = measureEW(sp,2.2635,width=[2.2580, 2.2690],continuum=[2.2510, 2.2580, 2.2705, 2.2760],recenter=False,**kwargs)
        ews[0] = tmp['ew']
        errs[0] = tmp['ew_unc']
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


def measureIndex(sp,ranges,method='ratio',sample='integrate',nsamples=100,noiseFlag=False,plot=False,verbose=False,**pkwargs):
    '''
    :Purpose: Measure an index on a spectrum based on defined methodology
                measure method can be mean, median, integrate
                index method can be ratio = 1/2, valley = 1-2/3, OTHERS
                output is index value and uncertainty
    .. will also come back to this one

# NOTE: NOISE FLAG ERROR IS BACKWARDS HERE    '''

# error checking on number of arguments provided
    if len(ranges) == 0:
        print('\nmeasureIndex needs at least 1 sample region')
        return numpy.nan, numpy.nan
    if len(ranges) == 1:
        method = 'single'
    if not isinstance(ranges[0],list) and not isinstance(ranges[0],numpy.ndarray) and not isinstance(ranges[0],tuple) and not isinstance(ranges[0],set):
        print('\nmeasureIndex needs a list of wavelength ranges, you entered {} which has types {}'.format(ranges,type(ranges[0])))
        return numpy.nan, numpy.nan
    for r in ranges:
        if len(r) != 2:
            print('\nProblem with range {} in input ranges {}: must be 2-element array'.format(r,ranges))
            return numpy.nan, numpy.nan
    if (len(ranges) < 2 and (method == 'ratio' or method == 'change')):
        print('\nIndex method {} needs at least 2 sample regions'.format(method))
        return numpy.nan, numpy.nan
    if (len(ranges) < 3 and (method == 'line' or method == 'allers' or method == 'inverse_line'  or method == 'sumnum' or method == 'sumdenom')):
        print('\nIndex method {} needs at least 3 sample regions'.format(method))
        return numpy.nan, numpy.nan


# define the sample vectors
    value = numpy.zeros(len(ranges))
    value_sim = numpy.zeros((len(ranges),nsamples))

# loop over all sampling regions
    for i,waveRng in enumerate(ranges):

# convert units
        if isUnit(waveRng):
            waveRng = (waveRng.to(sp.wave.unit)).value
        elif isUnit(waveRng[0]):
            waveRng = [(w.to(sp.wave.unit)).value for w in waveRng]
        else:
            waveRng = ((waveRng*DEFAULT_WAVE_UNIT).to(sp.wave.unit)).value
        xNum = (numpy.arange(0,nsamples+1.0)/nsamples)* \
            (numpy.nanmax(waveRng)-numpy.nanmin(waveRng))+numpy.nanmin(waveRng)

# identify measureable regions
        w = numpy.where(numpy.logical_and(\
            numpy.logical_and(numpy.isnan(sp.flux.value) == False,numpy.isnan(sp.noise.value) == False),\
            numpy.logical_and(sp.wave.value >= numpy.nanmin(waveRng),sp.wave.value <= numpy.nanmax(waveRng))))
        if len(w[0]) == 0:
            noiseFlag = True
            w = numpy.where(numpy.logical_and(\
                numpy.isnan(sp.flux) == False,\
                numpy.logical_and(sp.wave.value >= numpy.nanmin(waveRng),sp.wave.value <= numpy.nanmax(waveRng))))
        if len(w[0]) == 0:
            if verbose: print('Warning: no data in the wavelength range {}'.format(waveRng))
            return numpy.nan,numpy.nan

# compute intepolated flux and noise
#        print(waveRng,len(w),numpy.min(w),numpy.max(w))
        w = padWhereArray(w,len(sp.wave))
#        print(waveRng,len(w),numpy.min(w),numpy.max(w))
        f = interp1d(sp.wave.value[w],sp.flux.value[w],bounds_error=False,fill_value=numpy.nan)
        yNum = f(xNum)
        if noiseFlag == False:
            s = interp1d(sp.wave.value[w],sp.noise.value[w],bounds_error=False,fill_value=0.)
            yNum_e = s(xNum)

# first compute the actual value
        if (sample == 'integrate'):
            value[i] = trapz(yNum,xNum)
        elif (sample == 'average'):
            value[i] = numpy.nanmean(yNum)
        elif (sample == 'sum'):
            value[i,j] = numpy.nansum(yNum)
        elif (sample == 'median'):
            value[i] = numpy.median(yNum)
        elif (sample == 'maximum'):
            value[i] = numpy.nanmax(yNum)
        elif (sample == 'minimum'):
            value[i] = numpy.nanmin(yNum)
        else:
            value[i] = numpy.nanmean(yNum)

# now do Monte Carlo measurement of value and uncertainty
# THERE IS A PROBLEM HERE - ALL VALUE_SIM = NAN
        if noiseFlag == False: 
            for j in numpy.arange(0,nsamples):

# sample variance
# METHOD 1 - THIS SEEMS RIGHT BUT VERY SMALL UNCERTAINTIES
#                yVar = yNum+numpy.random.normal(0.,1.,size=len(yNum))*yNum_e
                yVar = yNum+numpy.random.normal(0.,3.,size=len(yNum))*yNum_e
# METHOD 1 - THIS IS THE SAME AS 1
#                yVar = numpy.random.normal(yNum,yNum_e)
# METHOD 3 - THIS SEEMS TOO COARSE
#                yVar = yNum+numpy.random.normal(0.,1.)*yNum_e
#                print(numpy.std((yVar-yNum)/yNum_e))

# choose function for measuring indices
                if (sample == 'integrate'):
                    value_sim[i,j] = trapz(yVar,xNum)
                elif (sample == 'average'):
                    value_sim[i,j] = numpy.nanmean(yVar)
                elif (sample == 'sum'):
                    value_sim[i,j] = numpy.nansum(yVar)
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
        if (method == 'single'):
            val = value[0]
            vals = value_sim[0,:]
        elif (method == 'ratio'):
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
        elif (method == 'sumnum'):
            val = (value[0]+value[1])/value[2]
            vals = (value_sim[0,:]+value_sim[1,:])/value_sim[2,:]
        elif (method == 'sumdenom'):
            val = value[0]/(value[1]+value[2])
            vals = value_sim[0,:]/(value_sim[1,:]+value_sim[2,:])
        elif (method == 'allers'):
            val = (((numpy.mean(ranges[0])-numpy.mean(ranges[1]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value[2] \
                + ((numpy.mean(ranges[2])-numpy.mean(ranges[0]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value[1]) \
                /value[0]
            vals = (((numpy.mean(ranges[0])-numpy.mean(ranges[1]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value_sim[2,:] \
                + ((numpy.mean(ranges[2])-numpy.mean(ranges[0]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value_sim[1,:]) \
                /value_sim[0,:]
        else:
            val = value[0]/value[1]
            vals = value_sim[0,:]/value_sim[1,:]

# PLOTTING/VISUALIZATION?
        if plot == True:
            from splat.plot import visualizeIndices, plotSpectrum
            bands = []
            for r in ranges: bands.append(r)
            inddict = {pkwargs.get('name','Index'): {'ranges': bands, 'value': value}}
    #        visualizeIndices(sp,inddict,**kwargs)        
            plotSpectrum(sp,bands=bands,bandlabels=[pkwargs.get('name','') for b in bands],**pkwargs)        

# output mean, standard deviation

    if noiseFlag == True: return val, numpy.nan
    return val, numpy.nanstd(vals)





def measureIndexSet(sp,ref='burgasser',verbose=False,**kwargs):
    '''
    :Purpose: Measures indices of ``sp`` from specified sets. Returns dictionary of indices.
    :param sp: Spectrum class object, which should contain wave, flux and noise array elements
    :param set: string defining which indices set you want to use; options include:

            - *bardalez*: H2O-J, CH4-J, H2O-H, CH4-H, H2O-K, CH4-K, K-J, H-dip, K-slope, J-slope, J-curve, H-bump, H2O-Y from `Bardalez Gagliuffi et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...794..143B>`_
            - *burgasser*: H2O-J, CH4-J, H2O-H, CH4-H, H2O-K, CH4-K, K-J from `Burgasser et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...637.1067B>`_
            - *tokunaga*: K1, K2 from `Tokunaga & Kobayashi (1999) <http://adsabs.harvard.edu/abs/1999AJ....117.1010T>`_
            - *reid*: H2O-A, H2O-B from `Reid et al. (2001) <http://adsabs.harvard.edu/abs/2001AJ....121.1710R>`_
            - *geballe*: H2O-1.2, H2O-1.5, CH4-2.2, Color-d2 from `Geballe et al. (2002) <http://adsabs.harvard.edu/abs/2002ApJ...564..466G>`_
            - *allers*: H2O, FeH-z, VO-z, FeH-J, KI-J, H-cont from `Allers et al. (2007) <http://adsabs.harvard.edu/abs/2007ApJ...657..511A>`_, `Allers & Liu (2013) <http://adsabs.harvard.edu/abs/2013ApJ...772...79A>`_
            - *testi*: sHJ, sKJ, sH2O-J, sH2O-H1, sH2O-H2, sH2O-K from `Testi et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...552L.147T>`_
            - *slesnick*: H2O-1, H2O-2, FeH from `Slesnick et al. (2004) <http://adsabs.harvard.edu/abs/2004ApJ...610.1045S>`_
            - *mclean*: H2OD from `McLean et al. (2003) <http://adsabs.harvard.edu/abs/2003ApJ...596..561M>`_
            - *rojas*: H2O-K2 from `Rojas-Ayala et al.(2012) <http://adsabs.harvard.edu/abs/2012ApJ...748...93R>`_

    :type ref: optional, default = 'burgasser'

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1555+0954')[0]
    >>> print splat.measureIndexSet(sp, ref = 'reid')
        {'H2O-B': (1.0531856077273236, 0.0045092074790538221), 'H2O-A': (0.89673318593633422, 0.0031278302105038594)}
    '''
# keyword parameters

    ref = kwargs.get('set',ref)

    if ('allers' in ref.lower()):
        reference = 'Indices from Allers et al. (2007), Allers & Liu (2013)'
        refcode = '2013ApJ...772...79A'
        names = ['H2O','FeH-z','VO-z','FeH-J','KI-J','H-cont']
        ranges = [([1.55,1.56],[1.492,1.502]),\
            ([0.99135,1.00465],[0.97335,0.98665],[1.01535,1.02865]),\
            ([1.05095,1.06505],[1.02795,1.04205],[1.07995,1.09405]),\
            ([1.19880,1.20120],[1.19080,1.19320],[1.20680,1.20920]),\
            ([1.23570,1.25230],[1.21170,1.22830],[1.26170,1.27830]),\
            ([1.54960,1.57040],[1.45960,1.48040],[1.65960,1.68040])]
        methods = ['ratio','allers','allers','allers','allers','allers']
        samples = ['average']*len(names)
    elif ('bardalez' in ref.lower()):
        reference = 'Indices from Bardalez Gagliuffi et al. (2014)'
        refcode = '2014ApJ...794..143B'
        names = ['H2O-J','CH4-J','H2O-H','CH4-H','H2O-K','CH4-K','K-J','H-dip','K-slope','J-slope','J-curve','H-bump','H2O-Y']
        ranges = [([1.14,1.165],[1.26,1.285]),\
            ([1.315,1.335],[1.26,1.285]),\
            ([1.48,1.52],[1.56,1.60]),\
            ([1.635,1.675],[1.56,1.60]),\
            ([1.975,1.995],[2.08,2.12]),\
            ([2.215,2.255],[2.08,2.12]),\
            ([2.06,2.10],[1.25,1.29]),\
            ([1.61,1.64],[1.56,1.59],[1.66,1.69]),\
            ([2.06,2.10],[2.10,2.14]),\
            ([1.27,1.30],[1.30,1.33]),\
            ([1.04,1.07],[1.26,1.29],[1.14,1.17]),\
            ([1.54,1.57],[1.66,1.69]),\
            ([1.04,1.07],[1.14,1.17])]
        methods = ['ratio','ratio','ratio','ratio','ratio','ratio','ratio','inverse_line','ratio','ratio','line','ratio','ratio',]
        samples = ['integrate']*len(names)
    elif ('burgasser' in ref.lower()):
        reference = 'Indices from Burgasser et al. (2006)'
        refcode = '2006ApJ...637.1067B'
        names = ['H2O-J','CH4-J','H2O-H','CH4-H','H2O-K','CH4-K','K-J']
        ranges = [([1.14,1.165],[1.26,1.285]),\
            ([1.315,1.335],[1.26,1.285]),\
            ([1.48,1.52],[1.56,1.60]),\
            ([1.635,1.675],[1.56,1.60]),\
            ([1.975,1.995],[2.08,2.12]),\
            ([2.215,2.255],[2.08,2.12]),\
            ([2.06,2.10],[1.25,1.29])]
        methods = ['ratio']*len(names)
        samples = ['integrate']*len(names)
    elif ('geballe' in ref.lower()):
        reference = 'Indices from Geballe et al. (2002)'
        refcode = '2002ApJ...564..466G'
        names = ['Color-d2','Cont-1.0','H2O-1.2','H2O-1.5','H2O-2.2','CH4-1.6','CH4-2.2']
        ranges = [\
            ([0.96,0.98],[0.735,0.755]),\
            ([1.04,1.05],[0.875,0.885]),\
            ([1.26,1.29],[1.13,1.16]),\
            ([1.57,1.59],[1.46,1.48]),\
            ([2.09,2.11],[1.975,1.995]),\
            ([1.56,1.6],[1.635,1.675]),\
            ([2.08,2.12],[2.215,2.255]),\
            ]
        methods = ['ratio']*len(names)
        samples = ['integrate']*len(names)
    elif ('kirkpatrick' in ref.lower()):
        reference = 'Indices from Kirkpatrick et al. (1999)'
        refcode = '1999ApJ...519..802K'
        names = ['Rb-a','Rb-b','Na-a','Na-b','Cs-a','Cs-b','TiO-a','TiO-b','VO-a','VO-b','CrH-a','CrH-b','FeH-a','FeH-b','Color-a','Color-b','Color-c','Color-d']
        ranges = [\
            ([.77752,.77852],[.78152,.78252],[.77952,.78052]),\
            ([.79226,.79326],[.79626,.79726],[.79426,.79526]),\
            ([.81533,.81633],[.81783,.81883]),\
            ([.81533,.81633],[.81898,.81998]),\
            ([.84961,.85061],[.85361,.85461],[.85161,.85261]),\
            ([.89185,.89285],[.89583,.89683],[.89385,.89485]),\
            ([.7033,.7048],[.7058,.7073]),\
            ([.8400,.8415],[.8435,.8470]),\
            ([.7350,.7370],[.7550,.7570],[.7430,.7470]),\
            ([.7860,.7880],[.8080,.8100],[.7960,.8000]),\
            ([.8580,.8600],[.8621,.8641]),\
            ([.9940,.9960],[.9970,.9990]),\
            ([.8660,.8680],[.8700,.8720]),\
            ([.9863,.9883],[.9908,.9928]),\
            ([.9800,.9850],[.7300,.7350]),\
            ([.9800,.9850],[.7000,.7050]),\
            ([.9800,.9850],[.8100,.8150]),\
            ([.9675,.9850],[.7350,.7550]),\
        ]
        methods = ['line','line','ratio','ratio','line','line','ratio','ratio','sumnum','sumnum','ratio','ratio','ratio','ratio','ratio','ratio','ratio','ratio']
        samples = ['sum']*len(names)
    elif ('martin' in ref.lower()):
        reference = 'Indices from Martin et al. (1999)'
        refcode = '1999AJ....118.2466M'
        names = ['PC3','PC6','CrH1','CrH2','FeH1','FeH2','H2O1','TiO1','TiO2','VO1','VO2']
        ranges = [\
            ([.823,.827],[.754,.758]),\
            ([.909,.913],[.650,.654]),\
            ([.856,.860],[.861,.865]),\
            ([.984,.988],[.997,1.001]),\
            ([.856,.860],[.8685,.8725]),\
            ([.984,.988],[.990,.994]),\
            ([.919,.923],[.928,.932]),\
            ([.700,.704],[.706,.710]),\
            ([.838,.842],[.844,.848]),\
            ([.754,.758],[.742,.746]),\
            ([.799,.803],[.790,.794]),\
        ]
        methods = ['ratio']*len(names)
        samples = ['integrate']*len(names)
    elif ('mclean' in ref.lower()):
        reference = 'Indices from McLean et al. (2003)'
        refcode = '2003ApJ...596..561M'
        names = ['H2OD']
        ranges = [([1.951,1.977],[2.062,2.088])]
        methods = ['ratio']
        samples = ['average']
    elif ('reid' in ref.lower()):
        reference = 'Indices from Reid et al. (2001)'
        refcode = '2001AJ....121.1710R'
        names = ['H2O-A','H2O-B']
        ranges = [([1.33,1.35],[1.28,1.30]),([1.47,1.49],[1.59,1.61])]
        methods = ['ratio']*len(names)
        samples = ['average']*len(names)
    elif ('rojas' in ref.lower()):
        reference = 'Indices from Rojas-Ayala et al.(2012)'
        refcode = '2012ApJ...748...93R'
        names = ['H2O-K2num','H2O-K2den']
        ranges = [([2.070,2.090],[2.235,2.255]),\
            ([2.235,2.255],[2.360,2.380])]
        methods = ['ratio']*len(names)
        samples = ['average']*len(names)
    elif ('slesnick' in ref.lower()):
        reference = 'Indices from Slesnick et al. (2004)'
        refcode = '2004ApJ...610.1045S'
        names = ['H2O-1','H2O-2','FeH']
        ranges = [([1.335,1.345],[1.295,1.305]),\
            ([2.035,2.045],[2.145,2.155]),\
            ([1.1935,1.2065],[1.2235,1.2365])]
        methods = ['ratio']*len(names)
        samples = ['average']*len(names)
    elif ('testi' in ref.lower()):
        reference = 'Indices from Testi et al. (2001)'
        refcode = '2001ApJ...552L.147T'
        names = ['sHJ','sKJ','sH2O_J','sH2O_H1','sH2O_H2','sH2O_K']
        ranges = [([1.265,1.305],[1.6,1.7]),\
            ([1.265,1.305],[2.12,2.16]),\
            ([1.265,1.305],[1.09,1.13]),\
            ([1.60,1.70],[1.45,1.48]),\
            ([1.60,1.70],[1.77,1.81]),\
            ([2.12,2.16],[1.96,1.99])]
        methods = ['change']*len(names)
        samples = ['average']*len(names)
    elif ('tokunaga' in ref.lower()):
        reference = 'Indices from Tokunaga & Kobayashi (1999)'
        refcode = '1999AJ....117.1010T'
        names = ['K1','K2']
        ranges = [([2.1,2.18],[1.96,2.04]),([2.2,2.28],[2.1,2.18])]
        methods = ['change']*len(names)
        samples = ['average']*len(names)
#        inds = numpy.zeros(len(names))
#        errs = numpy.zeros(len(names))
#        inds[0],errs[0] = measureIndex(sp,[2.1,2.18],[1.96,2.04],method='change',sample='average',**kwargs)
#        inds[1],errs[1] = measureIndex(sp,[2.2,2.28],[2.1,2.18],method='change',sample='average',**kwargs)
    else:
        print('{} is not one of the sets used for measureIndexSet'.format(ref))
        return numpy.nan

    result = {'ref': ref, 'bibcode': refcode}
    pkwargs = copy.deepcopy(kwargs) 
    if 'method' in list(pkwargs.keys()): del pkwargs['method']
    if 'sample' in list(pkwargs.keys()): del pkwargs['sample']
    for i,n in enumerate(names):
        ind,err = measureIndex(sp,ranges[i],method=methods[i],sample=samples[i],verbose=verbose,**pkwargs)
        result[n] = (ind,err)
        if verbose == True: print('Index {}: {:.3f}+/-{:.3f}'.format(n,ind,err))

# some exceptions
    if 'bardalez' in ref.lower():
        result['H-dip'] = tuple([i*0.5 for i in result['H-dip']])
        result['J-curve'] = tuple([i*2. for i in result['J-curve']])
    elif 'rojas' in ref.lower():
        ind = result['H2O-K2num'][0]/result['H2O-K2den'][0]
        err = [ind*numpy.sqrt((result['H2O-K2num'][1]/result['H2O-K2num'][0])**2+(result['H2O-K2den'][1]/result['H2O-K2den'][0])**2)]
        result['H2O-K2'] = (ind,err)
        del result['H2O-K2num'], result['H2O-K2den']


# output dictionary of indices
#    result = {names[i]: (inds[i],errs[i]) for i in numpy.arange(len(names))}
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

    h2ok2,h2ok2_e = measureIndexSet(sp, ref='rojas')['H2O-K2']
    ew = measureEWSet(sp,ref='rojas')
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



def addUserData(folders=[],default_info={},verbose=True):
    '''
    :Purpose:

        Reads in list of folders with properly processed model sets, checks them, and adds them to the SPECTRAL_MODELS global variable

    :Required Inputs:

        None

    :Optional Inputs:

        * :param folders = []: By default model folders are set in the .splat_spectral_models file; 
        alternately (or in addition) folders of models can be included as an input list.
        * :param default_info = {}: default parameter set to use for models; superceded by 'info.txt' file if present in model folder 
        * :param verbose = False: provide verbose feedback

    :Outputs:
        
        None, simply adds new model sets to SPECTRAL_MODELS global variable

    '''
# default information dictionary
    if len(default_info.keys()) == 0:
        default_info = {'folder': '', 'name': '', 'citation': '', 'bibcode': '', 'altnames': [], 'rawfolder': '', 'default': {'teff': 1500, 'logg': 5.0, 'z': 0.}}

# read in folders specified in .splat_spectral_models
    if os.path.exists(HOME_FOLDER+'/'+EXTERNAL_SPECTRAL_MODELS_FILE):
        with open(HOME_FOLDER+'/'+EXTERNAL_SPECTRAL_MODELS_FILE, 'r') as frd: x = frd.read()
        folders.extend(x.split('\n'))
        if '' in folders: folders.remove('')

# check and read in the new folders in the SPECTRAL_MODELS dictionary
    if len(folders) > 0:
        for i,f in enumerate(folders):
            flag = 0
            minfo = copy.deepcopy(default_info)
            if minfo['folder'] == '': minfo['folder'] = f
            if minfo['name'] == '': minfo['name'] = os.path.normpath(f).split('/')[-1]
            subfiles = os.listdir(minfo['folder'])
# no duplicate models (for now)
            if minfo['name'] in list(SPECTRAL_MODELS.keys()):
                print('\nWarning: spectral model set {} already exists in SPECTRAL_MODELS library; ignoring this one'.format(minfo['name']))
                flag = 1
# make sure RAW directory exists (indicates models have been processed)
            if 'RAW' not in subfiles:
                print('\nWarning: did not find a RAW directory in {}; please process this model set using splat.model._processModels()'.format(minfo['folder']))
                flag = 1
# check for additional information file
            if 'info.txt' not in subfiles:
                print('\nWarning: did not find info.txt file in {}; using default values for model information'.format(minfo['folder']))
            else:
#                try:
                f = minfo['folder']
                with open(f+'/info.txt', 'r') as frd: x = frd.read()
                lines = x.split('\n')
                if '' in lines: lines.remove('')
                lines = [x.split('\t') for x in lines]
                minfo = dict(lines)
                minfo['folder'] = f
                for k in list(default_info.keys()):
                    if k not in list(minfo.keys()): minfo[k] = default_info[k]
                for k in list(SPECTRAL_MODEL_PARAMETERS.keys()):
                    if k in list(minfo.keys()): minfo['default'][k] = minfo[k]
                    if 'default_'+k in list(minfo.keys()): minfo['default'][k] = minfo['default_'+k]
                minfo['altnames'] = minfo['altnames'].split(',')
#                except:
#                    print('\nWarning: problem reading info.txt file in {}; using default values for model information'.format(minfo['folder']))
            if flag == 0:
                if verbose == True: print('\nAdding {} models to SPLAT model set'.format(minfo['name']))
                SPECTRAL_MODELS[minfo['name']] = copy.deepcopy(minfo)
                del minfo
    return




# run test program if calling from command line
if __name__ == '__main__':
    pass
#    test_gooddata()
#    splat.test()
#    test_info()
