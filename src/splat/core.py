# -*- coding: utf-8 -*-
from __future__ import print_function

# imports - internal
import copy
import bz2
import glob
import gzip
import os
import shutil
import sys
import warnings

# imports - external
import matplotlib
#try: matplotlib.use('agg')
#except: pass
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
from scipy.integrate import trapezoid as trapz        # for numerical integration
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


# databases - using the .txt files for now, will need to change to SQL at a future date
# these are now done using pandas
#DB_SOURCES = ascii.read(DB_SOURCES_FILE)
#DB_SPECTRA = ascii.read(DB_SPECTRA_FILE)
# DB_SOURCES = pandas.read_csv(os.path.normpath(DB_SOURCES_FILE),delimiter='\t')
# DB_SPECTRA = pandas.read_csv(os.path.normpath(DB_SPECTRA_FILE),delimiter='\t')

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

#DB_ALL_SPECTRA = 
#for k in list(DB_SPECTRA_DEFAULT_PARAMETERS.keys()): DB_ALL_SPECTRA[k] = []

def _processNewData(verbose=True,**kwargs):
    pass

def _addData(folder,database=pandas.DataFrame(),input_file=DB_SPECTRA_INPUT_FILE,default_parameters=DB_SPECTRA_DEFAULT_PARAMETERS,allowed_file_extensions=SPECTRA_FILES_EXTENSIONS,verbose=False):
    '''
    '''
    sfile = os.path.join(folder,input_file)
    if os.path.exists(sfile) == True:
        try: pd = pandas.read_csv(sfile)
        except: 
            if verbose==True: print('Error reading in {} in folder {}; ignoring this folder'.format(input_file,folder))
            return database
        for k in list(pd.columns):
            tmp = checkDict(k,default_parameters)
            if tmp!=False: pd.rename({k: tmp},axis='columns',inplace=True)
    else:
        pd = pandas.DataFrame()
        files = []
        for k in allowed_file_extensions:
            files.extend(glob.glob('{}/*.{}'.format(folder,k)))
        if len(files) == 0:
            if verbose==True: print('Warning: no spectral data files in folder {}; ignoring this folder'.format(folder))
            return database
        files = list(set(files))
        pd['FILENAME'] = [os.path.basename(f) for f in files]
    for k in list(default_parameters.keys()):
        if k not in list(pd.columns):
            if default_parameters[k]['initialize'] == 'required':
                print('Warning: required keyword {} is not in database file {} in folder {}; ignoring this folder'.format(k,input_file,folder))
                return database
            elif default_parameters[k]['initialize'] == 'create':
                if k=='FOLDER': pd[k] = [folder]*len(pd)
                else: pd[k] = [default_parameters[k]['default']]*len(pd)
            elif default_parameters[k]['initialize'] == 'default':
                pd[k] = [default_parameters[k]['default']]*len(pd)
            else:
                raise ValueError('Unexpected initialize command {} for parameter {}'.format(default_parameters[k]['initialize'],k))
    if os.path.exists(sfile) == False and len(pd) > 0: 
#            shutil.copy2(sfile,sfile.replace('.csv','_old.csv'))
        pd.to_csv(sfile,index=False)
        if verbose==True: print('Created database file {} in folder {}; you may want to update with relevant data'.format(input_file,folder))
    if len(pd) > 0:
        if verbose==True: print('Adding {} sources from {} to spectral database'.format(len(pd),folder))
        if len(database) == 0: database = copy.deepcopy(pd)
        else: database = pandas.concat([database,pd],axis=0,ignore_index=True,sort=False)
    return database


def _initializeAllData(database=pandas.DataFrame(),override=False,verbose=False,**kwargs):
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

# default data provided with package    
    data_folders = directoryTree(SPECTRAL_DATA_FOLDER)
#    for d in data_folders: print(d)

# as specified in .splat_spectral_data
    if os.path.exists(EXTERNAL_DATA_FILE):
        with open(EXTERNAL_DATA_FILE, 'r') as frd: x = frd.read()
        mfolders = x.split('\n')
        if len(mfolders) > 0:
            for m in mfolders: data_folders.extend(directoryTree(m))
#    for d in data_folders: print(d)

# specified in environmental variable SPLAT_DATA as folder1:folder2:...
    if os.environ.get('SPLAT_DATA') != None:
        mfolders = str(os.environ['SPLAT_DATA']).split(':')
        if len(mfolders) > 0:
            for m in mfolders: data_folders.extend(directoryTree(m))
#    for d in data_folders: print(d)

# check that the folders exist
    data_folders = list(set(data_folders))
    for x in ['',' ']:
        if x in data_folders: data_folders.remove(x)
    rm = []
    for m in data_folders:
        if os.path.exists(m) == False: rm.append(m)
    if len(rm) > 0:
        print('could not find {}'.format(rm))
        for m in rm: data_folders.remove(m)
    if len(data_folders) == 0:
        if verbose == True: print('\nNo folders containing spectral data were found to be present')
        return
#    for d in data_folders: print(d)

# go through each folder and load in the relavant spectra.csv files, or create them
    for d in data_folders: database = _addData(d,database=database,verbose=verbose)

# drop duplicate file names
    if len(database)>0:
        ln = len(database)
        database.drop_duplicates(subset='FILENAME',keep='first',inplace=True)
        if len(database) != ln and verbose==True:
            print('Dropped {} duplicates; {} remaining'.format(ln-len(database),len(database)))

    return data_folders,database

DATA_FOLDERS,DB_ALL_SPECTRA = _initializeAllData()


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
    :param wave_label: label of wavelength (default = 'Wavelength')
    :param wave_unit: unit in which wavelength is measured (default = u.micron)
    :param wave_unit_label: label of the unit of wavelength (default = 'micron')
    :param flux_label: label of flux density (default = r'F\\_lambda')
    :param flux_unit: unit in which flux density is measured (default = u.erg/(u.cm**2 * u.s * u.micron)
    :param flux_unit_label: label of the unit of flux density (default = r'erg cm\\^-2 s\\^-1 micron\\^-1')
    :param resolution: Resolution of spectrum (default = median lam/lam step/2.)
    :param slitpixelwidth: Width of the slit measured in subpixel values (default = 3.33)
    :param slitwidth: Actual width of the slit, measured in arcseconds. Default value is the ``slitpixelwidth`` multiplied an assumed (for SpeX) spectrograph pixel scale of 0.15 arcseconds 
    :param header: header info of the spectrum (default = Table())
    :param filename: a string containing the spectrum's filename (default = '')
    :param file: same as filename (default = '')
    :param idkey: spectrum key of the desired spectrum (default = False)

    :Example:
       >>> import splat
       >>> sp = Spectrum(filename='myspectrum.fits')      # read in a file
       >>> sp = Spectrum('myspectrum.fits')               # same
       >>> sp = Spectrum(10002)                           # read in spectrum with data_key = 10002
       >>> sp = Spectrum(wave=wavearray,flux=fluxarray)   # create objects with wavelength & flux arrays
    '''

    def __init__(self, *args, **kwargs):
# some presets - DO WE NEED ALL OF THESE?
        sdb = pandas.DataFrame()
#        self.ismodel = kwargs.get('ismodel',False)
#        self.istransmission = kwargs.get('istransmission',False)
        self.wave_unit = kwargs.get('wave_unit',DEFAULT_WAVE_UNIT)
        self.wave_label = kwargs.get('wave_label',r'Wavelength')
#        self.wave_unit_label = kwargs.get('wave_unit_label',self.wave_unit)
        self.flux_unit = kwargs.get('flux_unit',DEFAULT_FLUX_UNIT)
        self.flux_label = kwargs.get('flux_label',r'F$_{\lambda}$')
        if kwargs.get('surface',False) == True: self.flux_label = 'Surface {}'.format(self.flux_label)
        if kwargs.get('apparent',False) == True: self.flux_label = 'Apparent {}'.format(self.flux_label)
        if kwargs.get('absolute',False) == True: self.flux_label = 'Absolute {}'.format(self.flux_label)
        if kwargs.get('normalized',False) == True: self.flux_label = 'Normalized {}'.format(self.flux_label)
        if kwargs.get('dimensionless',False)==True: self.flux_unit = u.dimensionless_unscaled
#        self.header = kwargs.get('header',fits.PrimaryHDU())
        self.filename = kwargs.get('file','')
        self.filename = kwargs.get('filename',self.filename)
        self.name = kwargs.get('name','')
        self.data_key = kwargs.get('idkey',False)
        self.data_key = kwargs.get('data_key',self.data_key)
        self.instrument = kwargs.get('instrument','')
        inst = checkInstrument(self.instrument) 
        if inst != False: 
            for k in list(INSTRUMENTS[inst].keys()): setattr(self,k,kwargs.get(k,INSTRUMENTS[inst][k]))
#        self.instrument_mode = kwargs.get('instrument_mode','')
#        self.runfast = kwargs.get('runfast',True)
        self.published = kwargs.get('published',False)
        self.bibcode = kwargs.get('bibcode','')
        self.history = kwargs.get('history',[])
        self.wave = []
        self.flux = []
        self.noise = []
        self.variance = []
        verbose = kwargs.get('verbose',True)

# return empty Spectrum object if something goes wrong
        empty = kwargs.get('empty',False)
        if empty==True: return(self)

# process arguments
# option 1: a filename is given
        if len(args) == 1:
            if isinstance(args[0],str): 
                self.filename = args[0]
                empty=True

# option 2: a spectrum data key is given
            elif isNumber(args[0])==True:
                self.data_key = int(args[0])
#                try:
                sdb = keySpectrum(self.data_key)
                if len(sdb) > 0: self.filename = sdb['DATA_FILE'].iloc[0]
                else: 
                    if verbose==True: print('Could not find data key {} in SPLAT database'.format(args[0]))
                empty=True
#                except:
#                    print('Warning: problem reading in spectral database')

# option 3: a dictionary is given - check for the appropriate keys
            elif isinstance(args[0],dict)==True:
                if 'wave' in list(args[0].keys()) and 'flux' in list(args[0].keys()):
                    for k in list(args[0].keys()): setattr(self,k,args[0][k])
                else: 
                    if verbose==True: print('Passed a dictionary that is missing wave and/or flux keys')
                    empty=True


# option 4: a pandas dataframe is given  
            elif isinstance(args[0],pandas.core.frame.DataFrame)==True:
                if 'wave' in list(args[0].columns) and 'flux' in list(args[0].columns):
                    for k in list(args[0].columns): setattr(self,k,args[0][k])
                else: 
                    if verbose==True: print('Passed a pandas array that is missing wave and/or flux columns')
                    empty=True

# option 5: a list or numpy array or arrays is given:  
# interpret as wave, flux, and optionally noise
            elif isinstance(args[0],list)==True or isinstance(args[0],numpy.ndarray)==True:
                if len(args[0]) >=2 and (isinstance(args[0][0],list)==True or isinstance(args[0][0],numpy.ndarray)==True):
                    self.wave = args[0][0]
                    self.flux = args[0][1]
                    if len(args[0]) >=3: self.noise = args[0][2]
                else: 
                    if verbose==True: print('Input array needs at least two arrays for wave and flux, but only contains {}'.format(len(args[0])))
                    empty=True
            else:
                if verbose==True: print('Do not recognize input as filename, database key, list, dict, numpy array, or pandas table')
                empty=True

# option 6: multiple lists or numpy arrays are given
# interpret as wave, flux, and optionally noise
        elif len(args) > 1:
            if (isinstance(args[0],list) or isinstance(args[0],numpy.ndarray)) and \
                (isinstance(args[1],list) or isinstance(args[1],numpy.ndarray)):
                self.wave = args[0]
                self.flux = args[1]
            else:
                if verbose==True: print('Multiple inputs need to be lists or numpy arrays')
                empty=True

            if len(args) > 2:
                if isinstance(args[2],list) or isinstance(args[2],numpy.ndarray):
                    kwargs['noise'] = kwargs.get('noise',args[2])

# option 7: wave, flux and optionally noise included in keyword arguments
# allow for this not being an option
        else: 
            for k in ['wave','flux','noise']:
                if len(kwargs.get(k,[])) > 0:
                    if isinstance(kwargs[k],list) or isinstance(kwargs[k],numpy.ndarray): 
                        setattr(self,k,kwargs[k])
                    else: 
                        if verbose==True: print('{} keyword must be a list or numpy array; you passed {}'.format(k,type(kwargs[k])))
                        empty=True
                else:
                    if k in 'wave flux':
#                        if verbose==True: print('{} keyword not included or passed empty list'.format(k))
                        empty=True


# read in file if container is currently empty but there is a filename
        if empty==True and self.filename != '':
            self.simplefilename = os.path.basename(self.filename)
            if self.name=='': self.name=self.simplefilename

# set up parameters
            mkwargs = {}
#            mkwargs['filename']=self.filename
            mkwargs['instrument']=self.instrument
            mkwargs['folder'] = kwargs.get('folder','./')
            mkwargs['wave_unit'] = self.wave_unit
            mkwargs['flux_unit'] = self.flux_unit
            mkwargs['delimiter']=kwargs.get('delimiter','')
            mkwargs['comment']=kwargs.get('comment','#')
            mkwargs['file_type']=kwargs.get('file_type','')
            mkwargs['verbose']=kwargs.get('verbose',False)
#            self.file = self.filename

# is this in the SPLAT database? if so use the default folder
            if self.filename in list(DB_SPECTRA['DATA_FILE']): 
                mkwargs['folder'] = SPECTRAL_DATA_FOLDER                
                sdb = searchLibrary(file=self.filename)

# return prior spectrum - THIS IS NOT WORKING SO COMMENTED OUT
#            if self.filename in list(SPECTRA_READIN.keys()) and self.runfast == True:
#                self = SPECTRA_READIN[self.filename]
#                return

# read in spectrum, being careful not to overwrite specifically assigned quantities
            rs = readSpectrum(self.filename,**mkwargs)
            if 'wave' in rs.keys():
                for k in list(rs.keys()): 
                    if k not in list(kwargs.keys()): setattr(self,k.lower(),rs[k])
                empty=False

# some extras
            others = ['pixel','mask','flag','flags','model','background']
            for o in others:
                if len(kwargs.get(o,'')) > 0:
                    setattr(self,o,kwargs[o])


# None of this worked; create an empty Spectrum object (can be used for copying)
        if empty==True:
            print('Warning: Creating an empty Spectrum object')
            return

# process spectral data
# convert to numpy arrays
        if not isinstance(self.wave,numpy.ndarray): self.wave = numpy.array(self.wave)
        if not isinstance(self.flux,numpy.ndarray): self.flux = numpy.array(self.flux)
        if len(self.noise)==0: self.noise = [numpy.nan]*len(self.wave)
        if not isinstance(self.noise,numpy.ndarray): self.noise = numpy.array(self.noise)
# enforce positivity of noise array
        if kwargs.get('no_negative_noise',True)==True: self.noise = numpy.absolute(self.noise)

# assure wave, flux, noise have units
        if not isUnit(self.wave): self.wave = numpy.array(self.wave)*self.wave_unit
        if not isUnit(self.flux): self.flux = numpy.array(self.flux)*self.flux_unit
        if not isUnit(self.noise): self.noise = numpy.array(self.noise)*self.flux_unit

# some conversions - move off to Spectrum methods
        # self.flam = self.flux
        # try: self.nu = self.wave.to(u.Hz,equivalencies=u.spectral())
        # except: pass
        # try:
        #     self.fnu = self.flux.to(u.Jansky,equivalencies=u.spectral_density(self.wave))
        #     self.fnu_unit = u.Jansky
        # except: pass
        # try: self.noisenu = self.noise.to(u.Jansky,equivalencies=u.spectral_density(self.wave))
        # except: pass
        # self.temperature = numpy.zeros(len(self.flux))
# calculate variance & S/N
#            self.variance = numpy.array([n**2 for n in self.noise.value])*self.noise.unit*self.noise.unit
        self.variance = self.noise**2
#            self.snr = self.computeSN()
        self.snr = numpy.nanmedian(self.flux/self.noise)

# estimate resolution - by default central lam/lam spacing/3
        self.resolution = kwargs.get('resolution',
            numpy.nanmedian(self.wave.value)/numpy.absolute(numpy.nanmedian(self.wave.value-numpy.roll(self.wave.value,3))))

# populate information on source and spectrum from database
# COULD POSSIBLY MOVE THIS TO A SEPARATE FUNCTION
#        print(sdb)
#        sdb = searchLibrary(**mkwargs)
        if len(sdb) > 0:
            for k in list(sdb.columns): setattr(self,k.lower(),str(sdb[k].iloc[0]))

# convert some data into numbers
        # kconv = ['ra','dec','julian_date','median_snr','resolution','airmass',\
        # 'jmag','jmag_error','hmag','hmag_error','kmag','kmag_error','source_key']
        # for k in kconv:
        #     try: setattr(self,k,float(getattr(self,k)))
        #     except: pass
# this is to make sure the database resolution is the default value
# IS THIS NECESSARY?
        # try:
        #     if kwargs.get('resolution',False) == False or kwargs.get('instrument',False) == False: 
        #         kwargs['resolution'] = self.resolution
        # except: pass

# instrument specific information
# THIS HAS BEEN REMOVED
#         hkys = list(self.header.keys())

# # automated stuff for spex data
#         if 'INSTRUME' in hkys:
#             if 'spex' in self.header['INSTRUME'].lower() and 'GRAT' in hkys:
#                 if 'lowres15' in self.header['GRAT'].lower() or 'prism' in self.header['GRAT'].lower(): self.instrument = 'SPEX-PRISM'
#                 if 'shortxd' in self.header['GRAT'].lower() or 'sxd' in self.header['GRAT'].lower(): self.instrument = 'SPEX-SXD'
#         if 'INSTR' in hkys:
#             if 'spex' in self.header['INSTR'].lower() and 'GRAT' in hkys:
#                 if 'lowres15' in self.header['GRAT'].lower() or 'prism' in self.header['GRAT'].lower(): self.instrument = 'SPEX-PRISM'
#                 if 'shortxd' in self.header['GRAT'].lower() or 'sxd' in self.header['GRAT'].lower(): self.instrument = 'SPEX-SXD'
#         if 'spex' in self.instrument.lower():
#             if 'observation_date' in list(self.__dict__.keys()): dt = self.observation_date
#             elif 'OBS-DATE' in hkys: dt = self.header['OBS-DATE'].replace('-','')
#             elif 'OBS_DATE' in hkys: dt = self.header['OBS_DATE'].replace('-','')
#             elif 'DATE_OBS' in hkys: dt = self.header['DATE_OBS'].replace('-','')
#             elif 'DATE-OBS' in hkys: dt = self.header['DATE-OBS'].replace('-','')
#             else: dt = '20000101'
# #            if int(dt) > 20140800:
# #                self.instrument = self.instrument.replace('SPEX','USPEX')

# # populate defaults
#         inst = checkInstrument(self.instrument)
#         if inst == False: 
# #            print(self.instrument)
#             inst = 'UNKNOWN'
#         self.instrument = inst
#         for k in list(INSTRUMENTS[self.instrument].keys()): setattr(self,k,INSTRUMENTS[self.instrument][k])
#         self.slitpixelwidth = (self.slitwidth/self.pixelscale).value
#         if kwargs.get('resolution',False) != False:
#             if kwargs['resolution'] > 0.:
#                 self.slitwidth = self.slitwidth*self.resolution/kwargs['resolution']
#                 self.slitpixelwidth = self.slitpixelwidth*self.resolution/kwargs['resolution']
#                 self.resolution = kwargs['resolution']
#         if kwargs.get('slitwidth',False) != False:
#             sl = kwargs['slitwidth']
#             if isUnit(sl): sl = sl.value
#             if s1 > 0.:
#                 self.resolution = self.resolution*(self.slitwidth.value)/sl
#                 self.slitpixelwidth = self.slitpixelwidth*s1/(self.slitwidth.value)
#                 self.slitwidth = kwargs['slitwidth']
#                 if not isUnit(self.slitwidth):
#                     self.slitwidth = self.slitwidth*u.arcsec
#         if kwargs.get('slitpixelwidth',False) != False:
#             if kwargs['slitpixelwidth'] > 0.:
#                 self.resolution = self.resolution*self.slitpixelwidth/kwargs['slitpixelwidth']
#                 self.slitwidth = self.slitwidth*kwargs['slitpixelwidth']/self.slitpixelwidth
#                 self.slitpixelwidth = kwargs['slitpixelwidth']
# #        else:
# #            kys = list(splat.INSTRUMENTS.keys())
# #            for k in list(INSTRUMENTS[kys[0]].keys()): setattr(self,k,'UNKNOWN')
# #            self.instrument_name = self.instrument
# #            self.slitpixelwidth = 3.
# #            self.resolution = (self.wave[int(0.5*len(self.wave))]/(self.wave[int(0.5*len(self.wave))+2]-self.wave[int(0.5*len(self.wave))-1])).value
# #            self.slitwidth = 1.0*u.arcsec # not sure I should set this
# #            print('Warning: {} is not one of the instruments defined for SPLAT'.format(self.instrument))
# #            print('Setting slit width to {}, slit pixel width to {} and resolution to {}'.format(self.slitwidth,self.slitpixelwidth,self.resolution))

#         self.dof = numpy.round(len(self.wave)/self.slitpixelwidth)


# information on model
# MOVE THIS TO READ MODEL FUNCTION
#         if self.ismodel == True:
#             self.teff = kwargs.get('teff',numpy.nan)
#             self.logg = kwargs.get('logg',numpy.nan)
#             self.z = kwargs.get('z',numpy.nan)
#             self.fsed = kwargs.get('fsed',numpy.nan)
#             self.cld = kwargs.get('cld',numpy.nan)
#             self.kzz = kwargs.get('kzz',numpy.nan)
#             self.slit = kwargs.get('slit',numpy.nan)
#             self.model = kwargs.get('model','')
#             mset = checkSpectralModelName(self.model)
# #            print(self.model,mset)
#             if mset != False:
#                 self.model = mset
#                 self.modelset = mset
#                 for k in list(SPECTRAL_MODELS[mset].keys()):
# #                    print(k,SPECTRAL_MODELS[mset][k])
#                     setattr(self,k.lower(),SPECTRAL_MODELS[mset][k])
# #            self.name = self.model+' model'
# #            self.shortname = self.name
#             self.published = 'Y'
# #            print(self.name)

# information on transmission/filter spectrum
# MOVE THIS TO READ TRANSMISSION FUNCTION
        # elif self.istransmission == True:
        #     self.shortname = self.name
        #     self.published = 'Y'
# populate header
# THIS HAS BEEN REMOVED            
        # else:
        #     kconv = {'designation': 'DESIG','name': 'NAME','shortname': 'SNAME','ra': 'RA_DEC','dec': 'DEC_DEC','slitwidth': 'SLTW_ARC','source_key': 'SRC_KEY','data_key': 'DATA_KEY','observer': 'OBSERVER', 'bibcode': 'BIB_DATA','program_pi': 'PI','program_number': 'PROGRAM','airmass': 'AIRMASS','reduction_spextool_version': 'VERSION','reduction_person': 'RED_PERS','reduction_date': 'RED_DATE','observation_date': 'OBSDATE','julian_date': 'JDATE','median_snr': 'SNR','resolution': 'RES', 'instrument': 'INSTRUME','wave_unit': 'XUNITS','flux_unit': 'YUNITS', 'wave_label': 'XTITLE', 'flux_label': 'YTITLE', 'opt_type': 'SPT_OPT', 'lit_type': 'SPT_LIT', 'nir_type': 'SPT_NIR', 'spex_type': 'SPT_SPEX', 'gravity_class_nir': 'GRAV_NIR', 'gravity_class_opt': 'GRAV_OPT', 'metallicity_class': 'ZCLASS', 'luminosity_class': 'LUMCLASS', 'color_extremity': 'COLOREX','mu': 'MU','mu_e': 'E_MU', 'mu_ra': 'MU_RA', 'mu_dec': 'MU_DEC', 'parallax': 'PARALLAX', 'parallax_e': 'E_PARALL', 'vtan': 'VTAN','vtan_e': 'E_VTAN','rv': 'RV','rv_e': 'E_RV', 'vsini': 'VSINI', 'vsini_e': 'E_VSINI','distance': 'DISTANCE', 'distance_e': 'E_DISTAN','j_2mass': 'J_2MASS', 'h_2mass': 'H_2MASS', 'ks_2mass': 'K_2MASS', 'j_2mass_e': 'E_J_2MAS', 'h_2mass_e': 'E_H_2MAS', 'ks_2mass_e': 'E_K_2MAS', 'object_type': 'OBJ_TYPE', 'binary': 'BINARY','sbinary': 'SPBINARY', 'companion_name': 'COMPNAME', 'cluster': 'CLUSTER' }
        #     for k in list(kconv.keys()):
        #         if kconv[k].upper() not in list(self.header.keys()):
        #             try:
        #                 self.header[kconv[k]] = getattr(self,k)
        #             except:
        #                 self.header[kconv[k]] = ''
        #     if 'DATE_OBS' not in list(self.header.keys()) and 'observation_date' in list(self.__dict__.keys()):
        #         self.header['DATE_OBS'] = '{}-{}-{}'.format(self.observation_date[:4],self.observation_date[4:6],self.observation_date[6:])
        #     if 'TIME_OBS' not in list(self.header.keys()) and 'observation_time' in list(self.__dict__.keys()):
        #         self.header['TIME_OBS'] = self.observation_time.replace(' ',':')


        if 'history' in list(self.__dict__.keys()):
            self.history.append('{} spectrum successfully loaded'.format(self.instrument))

# create a copy to store as the original
        self.original = copy.deepcopy(self)

# add to previous read spectra
        if self.filename != '': SPECTRA_READIN[self.filename] = self

        return


    def mapTo(self,other,overhang=0.1):
        '''
        Purpose: maps spectrum onto the wavelength scale of another spectrum

        '''
        self.toWaveUnit(other.wave.unit)
        flux_unit = self.flux.unit
        trng = [numpy.nanmin(other.wave.value),numpy.nanmax(other.wave.value)]
        dt = numpy.abs(trng[1]-trng[0])
        trng = [trng[0]-overhang*dt,trng[1]+overhang*dt]*other.wave.unit
        self.trim(trng)
        self.flux = reMap(self.wave.value,self.flux.value,other.wave.value)*flux_unit
        self.noise = reMap(self.wave.value,self.noise.value,other.wave.value)*flux_unit
        self.wave = other.wave
        self.variance = [x**2 for x in self.noise.value]*flux_unit**2
        if 'history' in list(self.__dict__.keys()):
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
        other.toWaveUnit(self.wave_unit)
        other.toFluxUnit(self.flux_unit)

# make a copy and fill in wavelength to be overlapping
        sp = copy.deepcopy(self)
        sp.wave = self.wave.value[numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))]
        sp.wave=sp.wave*self.wave_unit

# generate interpolated axes
        f1 = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        f2 = interp1d(other.wave.value,other.flux.value,bounds_error=False,fill_value=0.)
        n1 = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
        n2 = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)

# add
        sp.flux = (f1(sp.wave.value)+f2(sp.wave.value))*self.flux_unit

# uncertainty
        sp.variance = (n1(sp.wave.value)+n2(sp.wave.value))*(self.flux_unit**2)
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
        sp.flux_unit = self.flux.unit*other.flux.unit
        sp.flux_unit_label = str(self.flux.unit*other.flux.unit)
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
        sp.flux = (numpy.where(numpy.absolute(sp.flux.value) == numpy.inf, numpy.nan, sp.flux.value))*self.flux_unit/other.flux.unit
        sp.cleanNoise()

# update information
        sp.name = self.name+' / '+other.name
        sp.flux_unit = self.flux.unit/other.flux.unit
        sp.fflux_unit_label = str(self.flux.unit/other.flux.unit)
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
        sp.flux_unit = self.flux.unit/other.flux.unit
        sp.flux_unit_label = str(self.flux.unit/other.flux.unit)
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

    def computeSN(self,rng=[],statistic='median',offset=0.1):
        '''
        :Purpose: Compute a representative S/N value as the median value of S/N among the top 50% of flux values
        
        :Output: the S/N value

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.computeSN()
           115.96374031163553
        '''
        if isinstance(rng,int)==True or isinstance(rng,float)==True: rng = [rng]
        if len(rng)==0:
            w = numpy.where(self.flux.value > numpy.median(self.flux.value))
        if len(rng) == 1:
            rng = [rng[0]-offset,rng[0]+offset]
        if len(rng) >= 2:
            if isUnit(rng[0]): rng = [x.to(u.self.wave.unit).value for x in rng]
            w = numpy.where(numpy.logical_and(self.wave.value >=numpy.nanmin(rng),self.wave.value <= numpy.nanmax(rng)))
        if 'mean' in statistic: 
            return numpy.nanmean(self.flux.value[w]/self.noise.value[w])
        elif 'max' in statistic: 
            return numpy.nanmax(self.flux.value[w]/self.noise.value[w])
        elif 'min' in statistic: 
            return numpy.nanmin(self.flux.value[w]/self.noise.value[w])
        else:
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
           >>> sp = Spectrum(10001)
           >>> sp.addNoise(snr=10.)
           >>> sp.computeSN()
                9.5702358777108234
        '''
        if snr > 0.: 
            snr0 = self.computeSN()
            self.noise= [x*snr0/snr for x in self.noise.value]*self.noise.unit
            self.flux = self.flux+(numpy.random.normal(numpy.zeros(len(self.noise)),self.noise.value))*self.noise.unit
        return 

    def info(self,printout=True):
        '''
        Purpose
        -------

        Returns a summary of properties for the Spectrum object

        Parameters
        ----------

        None

        Outputs
        -------

        Text summary describing the Spectrum object

        Example
        -------

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

        Dependencies
        ------------
        None
        '''
#         if self.ismodel == True:
#             f = '\n{} for instrument {} with the following parmeters:'.format(self.modelset,self.instrument)
#             for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
#                 if hasattr(self,ms): f+='\n\t{} = {} {}'.format(ms,getattr(self,ms),SPECTRAL_MODEL_PARAMETERS[ms]['unit'])
# #            f+='\nSmoothed to slit width {} {}'.format(self.slit,SPECTRAL_MODEL_PARAMETERS['slit']['unit'])
#             f+='\n\nIf you use this model, please cite {}'.format(spbib.shortRef(SPECTRAL_MODELS[self.modelset]['bibcode']))
#             f+='\nbibcode = {}\n'.format(SPECTRAL_MODELS[self.modelset]['bibcode'])
#         elif self.istransmission == True:
#             f = '\n{} spectrum'.format(self.name)
# #            f+='\nSmoothed to slit width {} {}'.format(self.slit,SPECTRAL_MODEL_PARAMETERS['slit']['unit'])
#             f+='\n\nIf you use these data, please cite {}'.format(spbib.shortRef(self.bibcode))
#             f+='\nbibcode = {}\n'.format(self.bibcode)
#         else:
        f = '\n'
        if hasattr(self,'instrument'): f+='{} '.format(self.instrument)
        if hasattr(self,'name'): f+='spectrum of {}'.format(self.name)
        if hasattr(self,'observer') and hasattr(self,'date'): 
            if isinstance(self.observer,list):
                for i in range(len(self.observer)):
                    f+='\nObserved by {} on {}'.format(self.observer[i],properDate(self.date[i],output='YYYY MMM DD'))
            else:
                f+='\nObserved by {} on {}'.format(self.observer,properDate(self.date,output='YYYY MMM DD'))
        if hasattr(self,'airmass'): f+='\nAirmass = {:.2f}'.format(float(self.airmass))
        if hasattr(self,'designation'): f+='\nSource designation = {}'.format(self.designation)
        if hasattr(self,'median_snr'): f+='\nMedian S/N = {:.0f}'.format(float(self.median_snr))
        if hasattr(self,'spex_type'): f+='\nSpeX Classification = {}'.format(self.spex_type)
# these lines are currently broken
        # if hasattr(self,'lit_type'): 
        #     if isinstance(self.lit_type,list):
        #         for i in range(len(self.lit_type)):
        #             f+='\nLiterature Classification = {} from {}'.format(self.lit_type[i],spbib.shortRef(self.lit_type_ref[i]))
        #     else:
        #         f+='\nLiterature Classification = {} from {}'.format(self.lit_type,spbib.shortRef(self.lit_type_ref))
        if hasattr(self,'source_key') and hasattr(self,'data_key'): 
            # if isinstance(self.source_key,list):
            #     for i in range(len(self.source_key)):
            #         f+='\nSpectrum key = {}, Source key = {}'.format(int(self.data_key[i]),int(self.source_key[i]))
            # else:
            f+='\nSpectrum key = {}, Source key = {}'.format(int(self.data_key),int(self.source_key))
        if hasattr(self,'published') and hasattr(self,'bibcode'):
            if self.published == True or self.bibcode != '':
                f+='\n\nIf you use these data, please cite:'
            # if isinstance(self.data_reference,list):
            #     for i in range(len(self.data_reference)):
            #         f+='\n\t{}'.format(spbib.shortRef(self.data_reference[i]))
            #         f+='\n\tbibcode: {}'.format(self.data_reference[i])
            # else:
                f+='\n\t{}'.format(spbib.shortRef(self.bibcode))
                f+='\n\tbibcode: {}'.format(self.bibcode)
        else: f+='\n\nUNPUBLISHED DATA'

        f+='\n\nHistory:'
        if 'history' in list(self.__dict__.keys()):
            for h in self.history: f+='\n\t{}'.format(h)
        if printout==True: 
            print(f)
            return
        return f


    def export(self,filename='',clobber=True,csv=False,tab=True,delimiter='\t',save_header=True,save_noise=True,comment='#',file_type='',*args,**kwargs):
        '''
        :Purpose: 
            Exports a Spectrum object to either a fits or ascii file, depending on file extension given.  
            If no filename is explicitly given, the Spectrum.filename attribute is used. 
            If the filename does not include the full path, the file is saved in the current directory.  
            Spectrum.export and `Spectrum.save()`_ function in the same manner.

        .. _`Spectrum.save()` : api.html#Spectrum.save

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
        if file_type == '': file_type = filename.split('.')[-1]

# fits file
        if 'fit' in file_type:
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
            hdu.writeto(filename,overwrite=clobber)
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
            if save_noise == True:
                f.write('{}WAVELENGTH{}FLUX{}UNCERTAINTY\n'.format(comment,delimiter,delimiter))
                for i in range(len(self.wave.value)): f.write('{}{}{}{}{}\n'.format(self.wave.value[i],delimiter,self.flux.value[i],delimiter,self.noise.value[i]))
            else:
                f.write('{}WAVELENGTH{}FLUX\n'.format(comment,delimiter))
                for i in range(len(self.wave.value)): f.write('{}{}{}\n'.format(self.wave.value[i],delimiter,self.flux.value[i]))
            f.close()

        if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum saved to {}'.format(filename))
        return


    def save(self,*args,**kwargs):
        '''
        :Purpose: Exports a Spectrum object to either a fits or ascii file, depending on file extension given.  If no filename is explicitly given, the Spectrum.filename attribute is used. If the filename does not include the full path, the file is saved in the current directory.  

        `Spectrum.export()`_ and `Spectrum.save` function in the same manner.

        .. _`Spectrum.export()` : api.html#Spectrum.export
        '''
        self.export(*args,**kwargs)


    def flamToFnu(self):
        return self.toFnu()

    def toFnu(self):
        '''
        :Purpose: 
            Converts flux density r'F\\_nu' in units of Jy.  
            There is no change if the spectrum is already in r'F\\_nu' units.

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
        if self.flux_label == 'Temperature' or self.flux_label == r'${\lambda}F_{\lambda}$':
            self.reset()
        self.flux_unit = u.Jy
        self.flux_label = r'$F_{\nu}$'
        self.flux = self.flux.to(self.flux_unit,equivalencies=u.spectral_density(self.wave))
        self.noise = self.noise.to(self.flux_unit,equivalencies=u.spectral_density(self.wave))
        self.snr = self.computeSN()
        if 'history' in list(self.__dict__.keys()): self.history.append('Converted to Fnu units of {}'.format(self.flux_unit))
        return

    def fnuToFlam(self):
        return self.toFlam()

    def toFlam(self):
        '''
        :Purpose: 
            Converts flux density to r'F\\_lambda' in units of r'erg/s/cm\\^2/Hz'. 
            There is no change if the spectrum is already in r'F\\_lambda' units.
        
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
        if self.flux_label == 'Temperature' or self.flux_label == r'${\lambda}F_{\lambda}$':
            self.reset()
        self.flux_unit = DEFAULT_FLUX_UNIT
        self.flux_label = r'$F_{\lambda}$'
        self.flux = self.flux.to(self.flux_unit,equivalencies=u.spectral_density(self.wave))
        self.noise = self.noise.to(self.flux_unit,equivalencies=u.spectral_density(self.wave))
        self.variance = self.noise**2
        self.snr = self.computeSN()
        if 'history' in list(self.__dict__.keys()): self.history.append('Converted to Flam units of {}'.format(self.flux_unit))
        return

    def toSED(self):
        '''
        :Purpose: 
            Converts flux density in r'F\\_lambda' to lambda x r'F\\_lambda' with units of r'erg/s/cm\\^2'. 

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
        if self.flux_label == 'SED':
            return
# first convert to F_lambda
        self.toFlam()
# now convert to SED
        un = self.wave.unit*self.flux.unit
        self.flux = (self.wave*self.flux).to(DEFAULT_SED_UNIT)
        self.noise = (self.wave*self.noise).to(DEFAULT_SED_UNIT)
        self.variance = self.noise**2
        self.snr = self.computeSN()
        self.flux_unit = DEFAULT_SED_UNIT
        self.flux_label = r'${\lambda}F_{\lambda}$'
        if 'history' in list(self.__dict__.keys()): self.history.append('Converted to SED units of {}'.format(self.flux_unit))
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
        self.wave_unit = u.Angstrom
        self.wave = self.wave.to(self.wave_unit)
        if 'history' in list(self.__dict__.keys()): self.history.append('Converted wavelength to units of {}'.format(self.wave_unit))
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
           >>> sp = Spectrum(file='somespectrum',wave_unit=u.Angstrom)
           >>> sp.wave.unit
            Unit("Angstrom")
           >>> sp.toMicron()
           >>> sp.wave.unit
            Unit("micron")
        '''
        self.wave_unit = u.micron
        self.wave = self.wave.to(self.wave_unit)
        if 'history' in list(self.__dict__.keys()): self.history.append('Converted wavelength to units of {}'.format(self.wave_unit))
        return

    def toWaveUnit(self,wave_unit):
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
           >>> sp = Spectrum(file='somespectrum',wave_unit=u.Angstrom)
           >>> print(sp.wave.unit)
             Angstrom
           >>> sp.toWaveUnit(u.micron)
           >>> print(sp.wave.unit)
             micron
           >>> sp.toWaveUnit(u.s)
             Warning! failed to convert wavelength unit from micron to s; no change made
        '''
        try:
            self.wave = self.wave.to(wave_unit)
            self.wave_unit = wave_unit
            if 'history' in list(self.__dict__.keys()): self.history.append('Converted wavelength to units of {}'.format(self.wave_unit))
        except:
            print('\nWarning! failed to convert wavelength unit from {} to {}; no change made'.format(self.wave.unit,wave_unit))            
        return


    def toFluxUnit(self,flux_unit):
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
           >>> sp = Spectrum(file='somespectrum',wave_unit=u.Angstrom)
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
            self.flux = self.flux.to(flux_unit,equivalencies=u.spectral_density(self.wave))
            self.noise = self.noise.to(flux_unit,equivalencies=u.spectral_density(self.wave))
            self.variance = self.noise**2
            if 'history' in list(self.__dict__.keys()): self.history.append('Converted to flux units of {}'.format(self.flux_unit))
            self.snr = self.computeSN()
            self.flux_unit = flux_unit
        except:
            print('Warning! failed to convert flux unit from {} to {}; no change made'.format(self.flux.unit,flux_unit))
        return


    def toWavelengths(self,wave,force=True,verbose=False):
        '''
        :Purpose: 
            Maps a spectrum onto a new wavelength grid via interpolation or integral resampling

        :Required Inputs:
            :param wave: wavelengths to map to

        :Optional Inputs:
            :param force = True: proceed with conversion even if wavelength ranges are not perfectly in range
            :param verbose = False: provide verbose feedback
        
        :Outputs:
            None; Spectrum object is changed

        :Example:
           TBD
        '''
        
        if not isUnit(wave):
            wave = wave*self.wave.unit

        flag = 0
        if numpy.nanmin(((self.wave).to(wave.unit)).value) > numpy.nanmin(wave.value) or \
            numpy.nanmax(((self.wave).to(wave.unit)).value) < numpy.nanmax(wave.value):
            if force is not True: 
                print('\nWarning: input wavelength range {} to {} is partially outside spectrum wavelength range {} to {}'.format(numpy.nanmin(wave.value),numpy.nanmax(wave.value),numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)))
                flag = 1
        elif numpy.nanmin(((self.wave).to(wave.unit)).value) > numpy.nanmax(wave.value) or \
            numpy.nanmax(((self.wave).to(wave.unit)).value) < numpy.nanmin(wave.value):
            print('\nWarning: input wavelength range {} to {} is completely outside spectrum wavelength range {} to {}'.format(numpy.nanmin(wave.value),numpy.nanmax(wave.value),numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)))
            flag = 1

        if flag==0:
            self.toWaveUnit(wave.unit)
            self.trim([numpy.nanmin(wave)-0.05*(numpy.nanmax(wave)-numpy.nanmin(wave)),numpy.nanmax(wave)+0.05*(numpy.nanmax(wave)-numpy.nanmin(wave))])

# map onto wavelength grid; if spectrum has lower resolution, interpolate; otherwise integrate & resample
            flux_unit = self.flux.unit
            if len(self.wave) <= len(wave):
                f = interp1d(self.wave.value,self.flux.value,kind='slinear',bounds_error=False,fill_value='extrapolate')
                self.flux = f(wave.value)*flux_unit
                if numpy.isfinite(numpy.nanmin(self.noise))==True:
                    n = interp1d(self.wave.value,self.noise.value,bounds_error=False,fill_value='extrapolate')
                    self.noise = n(wave.value)*flux_unit
                else: self.noise = self.flux*numpy.nan
            else:
                self.flux = integralResample(self.wave.value,self.flux.value,wave.value,method='splat')*flux_unit
                if numpy.isfinite(numpy.nanmin(self.noise))==True:
                    self.noise = integralResample(self.wave.value,self.noise.value,wave.value,method='splat')*flux_unit
                else: self.noise = self.flux*numpy.nan
            self.wave = wave
            self.variance = self.noise**2
            if 'history' in list(self.__dict__.keys()): self.history.append('Resampled to new wavelength grid')
            self.snr = self.computeSN()
            if verbose: print('\nSpectrum resampled onto input wavelength grid')
        else: 
            print('\nNo change made to spectrum')
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
        wave_unit = self.wave.unit
        flux_unit = self.flux.unit

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
            wave_range = wave_range*wave_unit
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
        self.wave = wave_oversample*wave_unit
        self.flux = flux_oversample*flux_unit
        self.variance = var_oversample*(flux_unit**2)
        self.noise = numpy.array([x**0.5 for x in self.variance.value])*flux_unit

# smooth this in pixel space including oversample       
        self._smoothToSlitPixelWidth(pixel_resolution*oversample,method=method)

# resample down to final wavelength scale
        flux_smooth = integralResample(self.wave.value,self.flux.value,wave_out)
        if numpy.isnan(self.variance.value[0]) == False:
            var_smooth = integralResample(self.wave.value,self.variance.value,wave_out)
        else:
            var_smooth = [numpy.nan for f in flux_smooth]
        self.wave = wave_out*wave_unit
        self.flux = flux_smooth*flux_unit
        self.variance = var_smooth*(flux_unit**2)
        self.noise = numpy.array([x**0.5 for x in self.variance.value])*flux_unit

# do final trim
        self.trim(wave_range)

        if 'history' in list(self.__dict__.keys()): self.history.append('Converted to instrument {}: wave range = {} resolution = {}'.format(instrument,wave_range,resolution))
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
        if 'history' in list(self.__dict__.keys()): self.history.append('Shifted spectrum by radial velocity {}'.format(rv))

        return

    def redshift(self,z):
        '''
        :Purpose: Shifts the wavelength scale by a given (positive) redshift. A negative redshift is interpreted as shifting back to the restframe.
        
        :Example:
           >>> import splat
           >>> sp.redshift(0.56)
        '''
        shift = 1.+z
        if z<0: shift = 1./(1.-z)
        self.wave = self.wave*shift
        if 'history' in list(self.__dict__.keys()): 
            if z<0: self.history.append('Shifted spectrum by redshift {}'.format(z))
            else: self.history.append('Shifted spectrum back to restframe with redshift {}'.format(z))

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
           >>> sp = Spectrum(10001)
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

# NOTE: THIS IS CURRENTLY NOT FUNCTIONAL
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

        flux_unit = self.flux.unit
        a = (numpy.nanmax(self.wave.value)/numpy.nanmin(self.wave.value))**(1./len(self.wave))
        nwave = numpy.nanmin(self.wave.value)*(a**numpy.arange(len(self.wave)))
        nflux = self.flux.value*nwave
        ncflux = numpy.convolve(nflux, kern, 'same')
        self.flux = (ncflux/nwave)*flux_unit
        if 'history' in list(self.__dict__.keys()): self.history.append(report)
        if verbose==True: print(report)

        return


    def rotate(self,vsini,epsilon=0.6,verbose=False):
        '''
        :Purpose:

            Rotationally broaden the lines of a spectrum; a shortcut call to `Spectrum.broaden()`_

        .. _`Spectrum.broaden()` : api.html#Spectrum.broaden

        :Required Inputs:
            :param vsini: Rotational velocity in km/s

        :Optional Inputs:
            :param epsilon: epsilon parameter for limb darkening in rotational broadening (default = 0.6)
            :param verbose: provide extra feedback (default = False)

        :Outputs:

            None; Spectral flux is smoothed by rotational broadening
            
        :Example:
           >>> import splat
           >>> sp = Spectrum(10001)
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

        :param absolute: set to True to specify that the given magnitude is an absolute magnitude, which sets the ``flux_label`` keyword in the Spectrum object to 'Absolute Flux' (default = False)
        :param apparent: set to True to specify that the given magnitude is an apparent magnitude, which sets the ``flux_label`` flag in the Spectrum object to 'Apparent Flux' (default = False)

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

        if self.flux_label == 'Temperature' or self.flux_label == 'SED':
            self.reset()
        # if self.flux_unit != DEFAULT_FLUX_UNIT:
        #     self.toFlam()
        absolute = kwargs.get('absolute',False)
        apparent = kwargs.get('apparent',not absolute)
        apmag,apmag_e = filterMag(self,filt,**kwargs)
# NOTE: NEED TO INCORPORATE UNCERTAINTY INTO SPECTRAL UNCERTAINTY
        if (~numpy.isnan(apmag)):
            self.scale(10.**(0.4*(apmag-mag)))
            if absolute == True:
                self.flux_label = 'Absolute {}'.format(self.flux_label.split(' ')[-1])
                if 'history' in list(self.__dict__.keys()): self.history.append('Flux calibrated with {} filter to an absolute magnitude of {}'.format(filt,mag))
            if apparent == True:
                self.flux_label = 'Apparent {}'.format(self.flux_label.split(' ')[-1])
                if 'history' in list(self.__dict__.keys()): self.history.append('Flux calibrated with {} filter to an apparent magnitude of {}'.format(filt,mag))
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

#        from .photometry import filterMag
        return filterMag(self,filt,**kwargs)


# determine maximum flux, by default in non telluric regions
    def fluxMax(self,wrange=[],maskTelluric=True,**kwargs):
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

        if len(wrange) < 2:
            wrange = [numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)]*self.wave.unit
        if isUnit(wrange): wrange = wrange.to(u.micron).value
        if isUnit(wrange[0]): wrange = [x.to(u.micron).value for x in wrange]

        fl = self.flux.value[numpy.logical_and(self.wave.to(u.micron).value > numpy.nanmin(wrange),self.wave.to(u.micron).value < numpy.nanmax(wrange))]
        wv = self.wave.to(u.micron).value[numpy.logical_and(self.wave.to(u.micron).value > numpy.nanmin(wrange),self.wave.to(u.micron).value < numpy.nanmax(wrange))]

        if maskTelluric == True:            
            try:
                msk = numpy.zeros(len(fl))
                msk[numpy.where(numpy.logical_or(\
                    numpy.logical_and(wv > 1.35,wv < 1.42),
                    numpy.logical_and(wv > 1.8,wv < 1.95)))] = 1
                fl = fl[msk == 0]
                if len(fl) > 0: return numpy.nanmax(fl)*self.flux_unit
#                if isUnit(fl[0]): fl = [f.value for f in fl]
            except:
                pass
        
#        fl = self.flux.value[numpy.where(\
#                numpy.logical_and(self.wave > numpy.nanmin(self.wave.value)+0.1*(numpy.nanmax(self.wave)-numpy.nanmin(self.wave)),self.wave < numpy.nanmax(self.wave)-0.1*(numpy.nanmax(self.wave)-numpy.nanmin(self.wave))))]
#        if isUnit(fl[0]): fl = [f.value for f in fl]
        return numpy.nanmax(fl)*self.flux.unit


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
        med = kwargs.get('median',False)
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
                if med:                
                    scalefactor = numpy.nanmedian(self.flux.value[numpy.where(numpy.logical_and(self.wave.value > rng[0],self.wave.value < rng[1]))])
                else:
                    scalefactor = numpy.nanmax(self.flux.value[numpy.where(numpy.logical_and(self.wave.value > rng[0],self.wave.value < rng[1]))])
        else:
            if med:
                scalefactor = numpy.nanmedian(self.flux.value)
            else:
                scalefactor = self.fluxMax(**kwargs)
        if isUnit(scalefactor): scalefactor = scalefactor.value
        if scalefactor == 0.: print('\nWarning: normalize is attempting to divide by zero; ignoring')
        elif numpy.isnan(scalefactor) == True: print('\nWarning: normalize is attempting to divide by nan; ignoring')
        else: 
            self.scale(1./scalefactor)
            self.flux_label = 'Normalized {}'.format(self.flux_label.split(' ')[-1])
            if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum normalized')
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
           >>> sp = Spectrum(10001)                   # read in a source
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
        if 'history' in list(self.__dict__.keys()): self.history.append(report)

        return


    def remove(self,mask,others=[]):
        '''
        :Purpose: 

            Removes elements of wave, flux, noise specified by the True/False array

        :Required Input:

            :param mask: Either a mask array (an array of booleans, ints, or floats, where True or 1 removes the element)
                or a 2-element array that defines the wavelegnth range to mask 
        
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

# wavelength range given
        if len(msk) == 2:
            if not isUnit(msk): msk = msk*self.wave.unit
            msk.to(self.wave.unit)
            msk = generateMask(self.wave,mask_range=msk)

        if len(msk) != len(self.flux):
            print('\nWarning: mask must be same length ({}) as wave/flux arrays ({}); not removing any pixels'.format(len(msk),len(self.wave)))
            return
        if isinstance(msk[0],float): msk = [int(x) for x in msk]
        if isinstance(msk[0],int): msk = [True if x==1 else False for x in msk]
        if not isinstance(msk[0],bool) and not isinstance(msk[0],numpy.bool_): print('\nWarning: cannot interpret mask {}; not removing any pixels'.format(mask))

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
        if 'history' in list(self.__dict__.keys()): self.history.append('Mask applied to remove {} pixels'.format(cnt))
        return


    def replace(self,mask,replace_value,replace_noise=True,replace_flux=True):
        '''
        :Purpose: 

            Replaces flux and noise values using a mask and specified value

        :Required Inputs:

            :param mask: Either a mask array (an array of booleans, ints, or floats, where True or 1 removes the element)
                or a 2-element array that defines the wavelength range to replace 
            :param replace_value: value with which the masked elements should be replaced
        
        :Optional Inputs:

            :param replace_flux = True: replace elements in the noise array
            :param replace_noise = True: replace elements in the flux array

        :Output:

            Spectrum object has the flagged pixels replaced in flux, noise arrays, and optional vectors

        :Example:
           >>> import splat, numpy
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> num = splat.numberList('1-10,18,30-50')
           >>> mask = [not(p in num) for p in numpy.arange(len(sp.wave))]
           >>> sp.replace(mask,numpy.nan)
           >>> sp.showHistory()
                SPEX_PRISM spectrum successfully loaded
                Mask applied to replace 32 pixels with value nan
        '''
        msk = copy.deepcopy(mask)

# wavelength range given
        if len(msk) == 2:
            if not isUnit(msk): msk = msk*self.wave.unit
            msk.to(self.wave.unit)
            msk = generateMask(self.wave,mask_range=msk)

        if len(msk) != len(self.flux):
            print('\nWarning: mask must be same length ({}) as wave/flux arrays ({}); not removing any pixels'.format(len(msk),len(self.wave)))
            return
        if isinstance(msk[0],float): msk = [int(x) for x in msk]
        if isinstance(msk[0],int): msk = [True if x==1 else False for x in msk]
        if not isinstance(msk[0],bool): print('\nWarning: cannot interpret mask {}; not removing any pixels'.format(mask))

# check units of replacement value
        if not isUnit(replace_value) and not numpy.isnan(replace_value): 
            replace_value = replace_value*self.flux.unit
        if isUnit(replace_value):
            if replace_value.unit != self.flux.unit:
                try:
                    replace_value = replace_value.to(self.flux.unit)
                except:
                    replace_value = replace_value.value*self.flux.unit

# invert mask
        msk = numpy.array([not x for x in msk])
#        self.wave = (numpy.array(self.wave.value)[msk])*self.wave.unit
#        self.flux = (numpy.array(self.flux.value)[msk])*self.flux.unit
#        self.noise = (numpy.array(self.noise.value)[msk])*self.noise.unit
        if replace_flux: self.flux[msk] = replace_value
        if replace_noise: self.noise[msk] = replace_value
        self.variance = self.noise**2
        self.snr = self.computeSN()

        cnt = numpy.sum([1 if x == False else 0 for x in msk])
        if 'history' in list(self.__dict__.keys()): self.history.append('Mask applied to replace {} pixels'.format(cnt))
        return



    def clean(self,action='remove',replace_value=0.):
        '''
        :Purpose: 

            Cleans a spectrum by either removing or replacing nan values

        :Required Inputs:

            None
        
        :Optional Inputs:

            :param action = 'remove': specify as either 'remove' or 'replace'
            :param replace_value = 0.: for replace, what value to replace with

        :Output:

            Spectrum object is modified to have nan pixels "cleaned"

        :Example:
           >>> import splat,numpy
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.flux[0] = numpy.nan
           >>> sp.clean()
           >>> sp.remove(mask)
           >>> sp.showHistory()
                SPEX_PRISM spectrum successfully loaded
                Mask applied to remove 1 pixels
        '''
# clean out data points with nans in flux
        msk = [numpy.isnan(x) for x in self.flux.value]
# clean out data points with nans in noise
        msk2 = [numpy.isnan(x) for x in self.noise.value]
        msk = msk or msk2
        if action=='remove': self.remove(msk)
        elif action=='replace': self.replace(msk,replace_value)
        else: print('\nWarning: ambiguous action {} for clean; no action taken'.format(action))
        return

    def maskWave(self,rng,apply=True,**kwargs):
        '''
        Purpose
        -------

        Generate a mask array based on min/max range of wavelengths

        Parameters
        ----------

        rng : list of numpy.ndarray
            2-element array that specifies the minimum and maximum wavelength range to mask
            elements can be unitted values or pure numbers (interpretted to be same units as spectrum)

        apply = True : bool [optional]
            If True, apply the mask using maskFlux()


        Outputs
        -------

        None: changes the Spectrum object by:
            * creates and/or sets the mask keyword within the Spectrum object to an array of booleans
            * if apply==True, masks the relevant range of data

        Example
        -------
        
        TBD

        Dependencies
        ------------
            `maskFlux()`_
            `isUnit()`_
            numpy

        .. _`isUnit()` : api.html#splat.utilities.isUnit
        .. _`maskFlux()` : api.html#splat.core.maskFlux

        '''
        if 'mask' not in list(self.__dict__.keys()): self.mask = numpy.array([False]*len(self.wave))
        if isinstance(rng,list) == False and isinstance(rng,numpy.ndarray) == False: 
            raise ValueError('Trim range should be 2-element list; you passed {}'.format(rng))
        if isUnit(rng[0]): rng = [r.to(self.wave.unit).value for r in rng]
        w = numpy.where(numpy.logical_and(self.wave.value>=numpy.nanmin(rng),self.wave.value<=numpy.nanmax(rng)))
        if len(w[0])>0: self.mask[w] = True 
        if 'history' in list(self.__dict__.keys()): self.history.append('Masked {} pixels in wavelength range {} to {}'.format(len(numpy.where(self.mask==True)),numpy.nanmin(rng),numpy.nanmax(rng)))
        if apply==True: self.maskFlux(**kwargs)
        return


    def maskSN(self,limit,apply=True,**kwargs):
        '''
        Purpose
        -------

        Apply a mask based on S/N upper limit

        Parameters
        ----------

        limit : float
            S/N upper limit; all fluxes with values below this S/N are masked

        apply = True : bool [optional]
            If True, apply the mask using maskFlux()

        Outputs
        -------

        None: changes the Spectrum object by:
            * creates and/or sets the mask keyword within the Spectrum object to an array of booleans
            * if apply==True, masks the relevant range of data

        Example
        -------
        
        TBD

        Dependencies
        ------------
            `maskFlux()`_
            `isUnit()`_
            numpy

        .. _`isUnit()` : api.html#splat.utilities.isUnit
        .. _`maskFlux()` : api.html#splat.core.maskFlux

        '''
        if 'mask' not in list(self.__dict__.keys()): self.mask = numpy.array([False]*len(self.wave))
        if isUnit(limit): limit  = limit.value
        w = numpy.where(numpy.logical_or(self.flux.value/self.noise.value<limit,numpy.isnan(self.flux.value/self.noise.value)==True))
        if len(w[0])>0: 
            if len(w[0])==len(self.wave): 
                if verbose==True: print('Warning: all pixels would be masked with this S/N limit; not applying mask')
            else: 
                self.mask[w] = True 
                if apply==True: self.maskFlux(**kwargs)
                if 'history' in list(self.__dict__.keys()): self.history.append('Masked {} pixels with S/N < {}'.format(len(w[0]),limit))
        return


    def maskFlux(self,mask=[],replace_noise=True,replace_value=numpy.nan,others=[]):
        '''
        Purpose
        -------

        Applies a mask array to the elements of a Spectrum object

        Parameters
        ----------

        mask = [] : list or numpy.ndarray [optional]
            array with same length as Spectrum object wavelength array that specifies pixels to be masked
            array values can be floating point, int, or bool, with 1 = True (mask) and 0 = False (don't mask)

        replace_noise = True : bool [optional]
            If True, apply the mask to the uncertainty array

        replace_value = numpy.nan : float [optional]
            Value with which to replace masked elements of flux

        others = []] : list [optional]
            Other elements of the Spectrum object to apply mask (e.g., wave, flags, mask, etc.)
            These elements should be present and must be the same length as mask

        Outputs
        -------

        None: changes the Spectrum object by:
            * masking the flux array with numpy.nan (or whatever is replace_value)
            * optionally masks other arrays (noise, etc.)

        Example
        -------
        
           >>> import splat,numpy
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> num = splat.numberList('1-10,18,30-50')
           >>> mask = [not(p in num) for p in numpy.arange(len(sp.wave))]
           >>> sp.remove(mask)
           >>> sp.showHistory()
                SPEX_PRISM spectrum successfully loaded
                Mask applied to remove 32 pixels


        Dependencies
        ------------
            copy
            numpy

        '''

        msk = copy.deepcopy(mask)
        if len(msk) == 0:
            if 'mask' in list(self.__dict__.keys()): 
                msk = copy.deepcopy(self.mask)
            else: 
                print('Empty mask provided and no mask in Spectrum object; not removing any pixels')
                return

# wavelength range given
        # if len(msk) == 2:
        #     if not isUnit(msk): msk = msk*self.wave.unit
        #     msk.to(self.wave.unit)
        #     msk = generateMask(self.wave,mask_range=msk)

        if len(msk) != len(self.flux):
            print('\nWarning: mask must be same length ({}) as flux arrays ({}); not removing any pixels'.format(len(msk),len(self.flux)))
            return
        if isinstance(msk[0],float): msk = [int(x) for x in msk]
        if isinstance(msk[0],int): msk = [True if x==1 else False for x in msk]
        if not isinstance(msk[0],bool) and not isinstance(msk[0],numpy.bool_): 
            print('\nWarning: cannot interpret mask {}; not removing any pixels'.format(msk))
            return
#            print(type(msk[0]),isinstance(msk[0],bool))

# mask pixels
        self.flux[msk] = replace_value
        if replace_noise == True: self.noise[msk] = replace_value
        self.variance = self.noise**2

        if len(others) > 0:
            for k in others:
                if k in list(self.original.__dict__.keys()):
                    try:
                        x = getattr(self,k)
                        x[msk] = replace_value
                        setattr(self,k,x)
                    except:
                        pass

        cnt = numpy.sum([1 if x == False else 0 for x in msk])
        if 'history' in list(self.__dict__.keys()): self.history.append('Masking applied to {} pixels'.format(cnt))
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
        for k in list(self.original.__dict__.keys()):
            if k != 'history':
                try:
                    setattr(self,k,getattr(self.original,k))
                except:
                    pass

        if 'history' in list(self.__dict__.keys()): self.history.append('Returned to original state')
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


    def scale(self,factor,noiseonly=False):
        '''
        :Purpose: 

            Scales a Spectrum object's flux and noise values by a constant factor. 


        :Required Inputs:

            :param factor: A floating point number used to scale the Spectrum object

        :Optional Inputs:

            :param noiseonly = False: scale only the noise and variance, useful when uncertainty is under/over estimated

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
        if noiseonly == False: 
            self.flux = self.flux*factor
            self.flux_label = 'Scaled {}'.format(self.flux_label.split(' ')[-1])
        self.noise = self.noise*factor
        self.variance = self.noise**2
        if 'history' in list(self.__dict__.keys()): 
            if noiseonly == True: self.history.append('Spectrum noise scaled by a factor of {}'.format(factor))
            else: self.history.append('Spectrum scaled by a factor of {}'.format(factor))
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
        if 'history' in list(self.__dict__.keys()): 
            for h in self.history: print(h)
        return

    def smooth(self,smv,resolution=False,slitwidth=False,**kwargs):
        '''
        :Purpose: 

            Smoothes a spectrum either by selecting a constant slit width 
            (smooth in spectral dispersion space), pixel width (smooth in pixel space) 
            or resolution (smooth in velocity space). One of these options must be selected 
            for any smoothing to happen. Changes spectrum directly.

        :Required Inputs:

            None

        :Optional Inputs:

            :param method: the type of smoothing window to use; see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details (default = 'hamming')
            :param resolution: Constant resolution to smooth to; see `_smoothToResolution()`_ (default = None)
            :param slitPixelWidth: Number of pixels to smooth in pixel space; see `_smoothToSlitPixelWidth()`_ (default = None)
            :param slitWidth: Number of pixels to smooth in angular space; see `_smoothToPixelWidth()`_ (default = None)

        .. _`_smoothToResolution()` : api.html#Spectrum._smoothToResolution
        .. _`_smoothToPixelWidth()` : api.html#Spectrum._smoothToPixelWidth
        .. _`_smoothToSlitPixelWidth()` : api.html#Spectrum._smoothToSlitPixelWidth


        :Output: 

            None: spectrum is smoothed

        :Example:
           >>> import splat
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.smooth(resolution=30)
        '''
        method = kwargs.get('method','hamming')
        kwargs['method'] = method
        swargs = copy.deepcopy(kwargs)

# smooth by resolution        
        if isinstance(resolution,bool):
            if resolution == True: 
                self._smoothToResolution(smv,**swargs)
                return
        else:
            self._smoothToResolution(resolution,**swargs)
            return

# smooth by slit width (arcseconds)        
        if isinstance(slitwidth,bool):
            if slitwidth == True: 
                self._smoothToSlitWidth(smv,**swargs)
                return
        else:
            self._smoothToSlitWidth(slitwidth,**swargs)
            return

# smooth by slit pixel width (default)        
        if isinstance(slitwidth,bool):
            if slitwidth == True: 
                self._smoothToSlitWidth(smv,**swargs)
                return
        else:
            self._smoothToSlitWidth(slitwidth,**swargs)
            return
            

        smv = kwargs.get('slitPixelWidth',smv)
        smv = kwargs.get('slitpixelwidth',smv)
        smv = kwargs.get('pixelwidth',smv)
        smv = kwargs.get('pixel',smv)
        self._smoothToSlitPixelWidth(smv,**swargs)
        return

# old
        # smv = kwargs.get('resolution',None)
        # smv = kwargs.get('res',smv)
        # if smv != None:

        # smv = kwargs.get('slitWidth',None)
        # smv = kwargs.get('slitwidth',smv)
        # smv = kwargs.get('slit',smv)
        # if smv != None:
        #     self._smoothToSlitWidth(smv,**swargs)
        #     return
        # return


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

        .. _`Spectrum.smooth()` : api.html#Spectrum.smooth

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
                flx_sample = f(wave_sample)
                if numpy.isfinite(self.variance.value).any() == True:
                    v = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
                    var_sample = v(wave_sample)
            else:
                flx_sample = integralResample(self.wave.value,self.flux.value,wave_sample)
                if numpy.isfinite(self.variance.value).any() == True:
                    var_sample = integralResample(self.wave.value,self.variance.value,wave_sample)

# now convolve a function to smooth resampled spectrum
            window = signal.get_window(method,int(2.*numpy.round(oversample)))
            neff = numpy.sum(window)/numpy.nanmax(window)        # effective number of pixels
            flx_smooth = signal.convolve(flx_sample, window/numpy.sum(window), mode='same',method='direct')
            if numpy.isfinite(self.variance.value).any() == True:
                var_smooth = signal.convolve(var_sample, window/numpy.sum(window), mode='same',method='direct')/neff
# resample back to original wavelength grid
            wave_final = numpy.array(self.wave.value)
            wave_final = wave_final[wave_final <= numpy.max(wave_sample)]
            wave_final = wave_final[wave_final >= numpy.min(wave_sample)]
            if len(wave_final) >= len(wave_sample):
                f = interp1d(wave_sample,flx_smooth,bounds_error=False,fill_value=0.)
                flx_final = f(wave_final)
                if numpy.isfinite(self.variance.value).any() == True:
                    v = interp1d(wave_sample,var_smooth,bounds_error=False,fill_value=0.)
                    var_final = v(wave_final)
                else:
                    var_final = numpy.ones(len(wave_final))*numpy.nan
            else:
                flx_final = integralResample(wave_sample,flx_smooth,wave_final)
                if numpy.isfinite(self.variance.value).any() == True:
                    var_final = integralResample(wave_sample,var_smooth,wave_final)
                else:
                    var_final = numpy.ones(len(wave_final))*numpy.nan

#            f = interp1d(wave_sample,flx_smooth,bounds_error=False,fill_value=0.)
#            v = interp1d(wave_sample,var_smooth,bounds_error=False,fill_value=0.)
#            self.flux = f(self.wave.value)*self.flux_unit
#            self.variance = v(self.wave.value)*self.flux_unit**2
            self.wave = wave_final*self.wave.unit
            self.flux = flx_final*self.flux.unit
            self.variance = var_final*(self.flux.unit**2)
            self.noise = numpy.array([ns**0.5 for ns in self.variance.value])*self.flux.unit
            self.snr = self.computeSN()
#            self.slitpixelwidth = self.slitpixelwidth*self.resolution/res
#            self.slitwidth = self.slitwidth*self.resolution/res
            self.resolution = res
            if 'history' in list(self.__dict__.keys()): self.history.append('Smoothed to a constant resolution of {}'.format(self.resolution))
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

        .. _`Spectrum.smooth()` : api.html#Spectrum.smooth

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
        method = kwargs.get('method','hamming')
        kwargs['method'] = method

# add in resolution keyword if not present
        if 'resolution' not in list(self.__dict__.keys()):
            i = int(0.5*len(self.wave.value))
            self.resolution = self.wave.value[i]/numpy.abs(self.wave.value[i]-self.wave.value[i+1])

# do nothing if requested resolution is higher than current resolution
        if width > 2.:
# convolve a function to smooth spectrum
            window = signal.get_window(method,int(numpy.round(width)))
            neff = numpy.nansum(window)/numpy.nanmax(window)        # effective number of pixels
            self.flux = signal.convolve(self.flux.value, window/numpy.nansum(window), mode='same',method='direct')*self.flux_unit
            self.variance = signal.convolve(self.variance.value, window/numpy.nansum(window), mode='same',method='direct')/neff*(self.flux_unit**2)
            self.noise = [n**0.5 for n in self.variance.value]*self.flux_unit
            self.snr = self.computeSN()
            self.resolution = self.resolution/width
#            self.slitwidth = self.slitwidth*width/self.slitpixelwidth
#            self.slitpixelwidth = width
            if 'history' in list(self.__dict__.keys()): self.history.append('Smoothed to pixel width of {}'.format(width))
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

        .. _`Spectrum.smooth()` : api.html#Spectrum.smooth

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
        method = kwargs.get('method','hamming')
        kwargs['method'] = method
        if not isUnit(width): width=width*u.arcsec

# add in resolution keyword if not present
        if 'resolution' not in list(self.__dict__.keys()):
            i = int(0.5*len(self.wave.value))
            self.resolution = self.wave.value[i]/numpy.abs(self.wave.value[i]-self.wave.value[i+1])
# add in resolution keyword if not present
        if 'slitpixelwidth' not in list(self.__dict__.keys()): self.slitpixelwidth = 1

# if slitwidth not present, do nothing
        if 'slidwidth' not in list(self.__dict__.keys()): return

        pwidth = self.slitpixelwidth*(width/self.slitwidth).value
        self._smoothToSlitPixelWidth(pwidth,**kwargs)
        return


    def toSurface(self,radius):
        '''
        :Purpose: Convert to surface fluxes given a radius, assuming at absolute fluxes
        .. note:: UNTESTED, NEED TO ADD IN UNCERTAINTIES
        '''
        if 'Surface' in self.flux_label: return
        if 'Absolute' not in self.flux_label: 
            print('To convert to surface fluxes you must first scale spectrum to absolute (10 pc) flux units')
            return
        r = copy.deepcopy(radius)
        if not isUnit(r): r=r*const.R_sun
        self.scale((((10.*u.pc/r).to(u.m/u.m)).value)**2)
        if 'history' in list(self.__dict__.keys()): self.history.append('Converted to surface fluxes assuming a radius of {} solar radii'.format((r/const.R_sun).to(u.m/u.m)))
        self.flux_label = 'Surface {}'.format(self.flux_label.split(' ')[-1])
        return


    def toBrightnessTemperature(self,limbdarkening=False,limbdarkeningcoeff = 0.7):
        '''
        :Purpose: Convert to surface fluxes given a radius, assuming at absolute fluxes
        .. note: UNTESTED
        '''
        if 'Temperature' in self.flux_label: return
        if 'Surface' not in self.flux_label: 
            print('To convert to brightness temperature you must first scale spectrum to surface flux units')
            return
        if self.flux_unit != DEFAULT_FLUX_UNIT: self.toFlam()
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
        if 'history' in list(self.__dict__.keys()): self.history.append('Converted to brightness temperature')
        self.flux_label = 'Temperature'
        self.flux_unit = u.K
        return

    def toTemperature(self):
        self.brightnessTemperature()


    def sample(self,rng,method='median',verbose=ERROR_CHECKING):
        '''
        :Purpose: 
            Obtains a sample of spectrum over specified wavelength range

        :Required Inputs: 

            :param range: the range(s) over which the spectrum is sampled
            a single 2-element array or array of 2-element arrays

        :Optional Inputs: 

            None

        :Example:
            TBD
        '''

# single number = turn into small range
        if isinstance(rng,float):
            rng = [rng-0.01*(numpy.nanmax(self.wave.value)-numpy.nanmin(self.wave.value)),rng+0.01*(numpy.nanmax(self.wave.value)-numpy.nanmin(self.wave.value))]

        if not isinstance(rng,list): rng = list(rng)

        if isUnit(rng[0]):
            try: rng = [r.to(self.wave.unit).value for r in rng]
            except: raise ValueError('Could not convert trim range unit {} to spectrum wavelength unit {}'.format(rng.unit,self.wave.unit))

        w = numpy.where(numpy.logical_and(self.wave.value >= rng[0],self.wave.value <= rng[1]))
        if len(w[0])>0:
            if method.lower() in ['median','med']: val = numpy.nanmedian(self.flux.value[w])
            elif method.lower() in ['mean','average','ave']: val = numpy.nanmean(self.flux.value[w])
            elif method.lower() in ['max','maximum']: val = numpy.nanmax(self.flux.value[w])
            elif method.lower() in ['min','minimum']: val = numpy.nanmin(self.flux.value[w])
            elif method.lower() in ['std','stddev','stdev','rms']: val = numpy.nanstd(self.flux.value[w])
            elif method.lower() in ['unc','uncertainty','noise','error']: val = numpy.nanmedian(self.noise.value[w])
            elif method.lower() in ['sn','snr','signal-to-noise','s/n']: val = numpy.nanmedian(self.flux.value[w]/self.noise.value[w])
            else: raise ValueError('Did not recongize sampling method {}'.format(method))
            return val
        else:
            if verbose==True: print('Sampling range {} outside wavelength range of data'.format(rng))
            return numpy.nan


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
        if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum trimmed to range {}'.format(rng))
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



class NewSpectrum(object):
    '''
    :Description: 
        (Updated) class for spectral data in SPLAT.
        This is a temporary structure until astropy.specutils is completed

    WORK IN PROGRESS    
    '''
    def __init__(self,*args,**kwargs):
        core_attributes = {
            'instrument': {'default': DEFAULT_INSTRUMENT, 'altname': ['INST','INSTR']},
            'name': {'default': 'Unknown source', 'altname': ['SOURCE','SOURCE_NAME']},
            'wave': {'default': [], 'altname': ['WAVELENGTH','W','LAM','LAMBDA']},
            'flux': {'default': [], 'altname': ['FLX','FLUX_DENSITY']},
            'noise': {'default': [], 'altname': ['UNCERTAINTY','UNC','ERROR','ERR','NS']},
            'wave_unit': {'default': DEFAULT_WAVE_UNIT, 'altname': ['WUNIT','WAVELENGTH_UNIT','LAM_UNIT','LAMBDA_UNIT']},
            'flux_unit': {'default': DEFAULT_FLUX_UNIT, 'altname': ['FUNIT','FLX_UNIT','FLUX_DENSITY_UNIT']},
            'header': {'default': {}, 'altname': ['HD','HEAD']},
            'filename': {'default': '', 'altname': ['FILE','FILE_NAME']},
            'folder': {'default': '', 'altname': ['FILE','FILE_NAME']},
            }

# set all inputs
        for k in list(core_attributes.keys()): 
            setattr(self,k,core_attributes[k]['default'])
            if k in list(kwargs.keys()): setattr(self,k,kwargs[k])
            if k.lower() in list(kwargs.keys()): setattr(self,k,kwargs[k.lower()])
            if k.upper() in list(kwargs.keys()): setattr(self,k,kwargs[k.upper()])
            for kk in core_attributes[k]['altname']:
                if kk in list(kwargs.keys()): setattr(self,k,kwargs[kk])
                if kk.lower() in list(kwargs.keys()): setattr(self,k,kwargs[kk.lower()])
                if kk.upper() in list(kwargs.keys()): setattr(self,k,kwargs[kk.upper()])
        for k in list(kwargs.keys()): setattr(self,k.lower(),kwargs[k])

# check if filename is given
        if len(args) > 0:
            if isinstance(args[0],str): self.filename = args[0]

# read in file
        if self.filename != '':
            if os.path.exists(self.filename)==False:
                if os.path.exists(os.path.join(self.folder,self.filename))==False:
                    raise ValueError('Could not find file {} or {}'.format(self.filename,os.path.join(self.folder,self.filename)))
                else:
                    self.filename = os.path.join(self.folder,self.filename)
# Have we already read this in? if so, just copy
# NOTE: THIS IS NOT WORKING SO COMMENTED OUT
#            if self.filename in list(SPECTRA_READIN.keys()):
#                self = copy.deepcopy(SPECTRA_READIN[self.filename])
#                return

# read in spectrum
# NOTE: NEED TO PUT IN SOME CONTEXT HERE FOR READING DIFFERENT TYPES OF FILES
# THIS WILL BE CONTAINED IN READSPECTRUM CALL
            if 'file_type' not in list(self.__dict__.keys()): self.file_type=''
            rs = readSpectrum(self.filename,file_type=self.file_type)
            for k in list(rs.keys()): 
                setattr(self,k.lower(),rs[k])
                for kk in list(core_attributes.keys()):
                    if k in core_attributes[kk]['altname']: setattr(self,kk,rs[k])
                    if k.lower() in core_attributes[kk]['altname']: setattr(self,kk,rs[k])
                    if k.upper() in core_attributes[kk]['altname']: setattr(self,kk,rs[k])

# let's make sure we have what we need (error checking)
        for k in ['wave','flux','noise']: 
            if not isinstance(getattr(self,k),numpy.ndarray): setattr(self,k,numpy.array(getattr(self,k)))
        if len(self.wave) == 0: raise ValueError('Spectrum object must be initiated with a wave array or this must be in read in file')
        if len(self.flux) == 0: raise ValueError('Spectrum object must be initiated with a flux array or this must be in read in file')
        if len(self.noise) == 0: self.noise = numpy.array([numpy.nan]*len(self.flux))
        if len(self.wave) != len(self.flux): raise ValueError('Input error: wavelength array has length {} while flux array has length {}'.format(len(self.wave),len(self.flux)))
        if len(self.wave) != len(self.noise): raise ValueError('Input error: wavelength array has length {} while unc array has length {}'.format(len(self.wave),len(self.noise)))

# set units
        if not isUnit(self.wave_unit): self.wave_unit = DEFAULT_WAVE_UNIT
        if not isUnit(self.flux_unit): self.flux_unit = DEFAULT_FLUX_UNIT
        tmp = checkDict(self.instrument,INSTRUMENTS)
        if tmp!=False:
            self.instrument=tmp
            for k in list(INSTRUMENTS[tmp].keys()): setattr(self,k,INSTRUMENTS[tmp][k])
        if isUnit(self.wave)==False: self.wave = self.wave*self.wave_unit
        try: self.wave=self.wave.to(self.wave_unit)
        except: pass
        if isUnit(self.flux)==False: self.flux = self.flux*self.flux_unit
        try: self.flux=self.flux.to(self.flux_unit)
        except: pass
        if isUnit(self.noise)==False: self.noise = self.noise*self.flux_unit
        try: self.noise=self.noise.to(self.flux_unit)
        except: pass

# clean up
        self.variance = self.noise**2
        self.original = copy.deepcopy(self)
        if 'history' in list(self.__dict__.keys()): self.history = ['{} spectrum of {} successfully loaded'.format(self.instrument,self.name)]
        return

    def setbase(self):
        '''
        :Purpose: Sets the current state of spectrum as the default, eliminates prior original
        '''
        self.original = copy.deepcopy(self)
        return

    def reset(self):
        '''
        :Purpose: Reset a spectrum to its original state
        '''
        for k in list(self.original.__dict__.keys()):
            if k != 'history':
                try:
                    setattr(self,k,getattr(self.original,k))
                except:
                    pass

        if 'history' in list(self.__dict__.keys()): self.history.append('Returned to original state')
#        self.original = copy.deepcopy(self)
        return

    def clean(self):
        '''
        :Purpose: Cleans up spectrum elements to make sure they are properly configured
        '''
# set up units
        try: self.wave_unit = self.wave.unit
        except: self.wave_unit = DEFAULT_WAVE_UNIT
        try: self.flux_unit = self.flux.unit
        except: self.flux_unit = DEFAULT_FLUX_UNIT

# clean wavelength vector
        if isUnit(self.wave): self.wave = numpy.array(self.wave.value)*self.wave_unit
        else: self.wave = numpy.array(self.wave)*self.wave_unit

# clean flux vector
        for k in ['flux','noise']:
            if isUnit(getattr(self,k)): setattr(self,k,numpy.array(getattr(self,k).value)*self.flux_unit)
            else: setattr(self,k,numpy.array(getattr(self,k))*self.flux_unit)
# set variance
        self.variance = self.noise**2
# need to: 
        if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum cleaned')
        return

    def __copy__(self):
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
        :Purpose: 

            A representation of addition for Spectrum objects which correctly interpolates 
            as a function of wavelength and combines variances

        :Output: 

            A new Spectrum object equal to the spectral sum of the inputs
        '''          
        try: other.wave = other.wave.to(self.wave.unit)
        except: raise ValueError('Cannot add spectra with wave units {} and {}'.format(self.wave.unit,other.wave.unit))
        try: other.flux = other.flux.to(self.flux.unit)
        except: raise ValueError('Cannot add spectra with flux units {} and {}'.format(self.flux.unit,other.flux.unit))

# establish the baseline wavelength grid
        out = copy.deepcopy(self)
        wave = numpy.array(copy.deepcopy(self.wave.value))
        wself = numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))
        out.wave = wave[wself]
        out.wave = out.wave*self.wave.unit
#       wother = numpy.where(numpy.logical_and(other.wave.value <= numpy.nanmax(out.wave.value),other.wave.value >= numpy.nanmin(out.wave.value)))

# do the math
        fself = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        fother = interp1d(other.wave.value,(other.flux.to(self.flux.unit)).value,bounds_error=False,fill_value=0.)
        out.flux = (fself(out.wave.value)+fother(out.wave.value))*self.flux.unit
# special for variance
        fself = interp1d(self.wave.value,(self.noise.value)**2,bounds_error=False,fill_value=0.)
        fother = interp1d(other.wave.value,((other.noise.to(self.flux.unit)).value)**2,bounds_error=False,fill_value=0.)
        if numpy.random.choice(numpy.isfinite(self.noise.value))==True and numpy.random.choice(numpy.isfinite(other.noise.value))==True: 
            out.variance = (fself(out.wave.value)+fother(out.wave.value))*(self.flux.unit**2)
        elif numpy.random.choice(numpy.isfinite(self.noise.value))==True:
            out.variance = fself(out.wave.value)*(self.flux.unit**2)
        elif numpy.random.choice(numpy.isfinite(other.noise.value))==True:
            out.variance = fother(out.wave.value)*(self.flux.unit**2)
        else:
            out.variance = numpy.array([numpy.nan]*len(out.wave))
        out.noise = out.variance**0.5

# update other information
        out.name = self.name+' + '+other.name
        out.history.append('Sum of {} and {}'.format(self.name,other.name))
        out.original = copy.deepcopy(out)
        return out
    
    def __sub__(self,other):
        '''
        :Purpose: 

            A representation of subtraction for Spectrum objects which correctly interpolates 
            as a function of wavelength and combines variances

        :Output: 

            A new Spectrum object equal to the spectral difference of the inputs
        '''        
        try: other.wave = other.wave.to(self.wave.unit)
        except: raise ValueError('Cannot subtract spectra with wave units {} and {}'.format(self.wave.unit,other.wave.unit))
        try: other.flux = other.flux.to(self.flux.unit)
        except: raise ValueError('Cannot subtract spectra with flux units {} and {}'.format(self.flux.unit,other.flux.unit))

# establish the baseline wavelength grid
        out = copy.deepcopy(self)
        wave = numpy.array(copy.deepcopy(self.wave.value))
        wself = numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))
        out.wave = wave[wself]
        out.wave = out.wave*self.wave.unit
#       wother = numpy.where(numpy.logical_and(other.wave.value <= numpy.nanmax(out.wave.value),other.wave.value >= numpy.nanmin(out.wave.value)))

# do the math
        fself = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        fother = interp1d(other.wave.value,(other.flux.to(self.flux.unit)).value,bounds_error=False,fill_value=0.)
        out.flux = (fself(out.wave.value)-fother(out.wave.value))*self.flux.unit
# special for variance
        fself = interp1d(self.wave.value,(self.noise.value)**2,bounds_error=False,fill_value=0.)
        fother = interp1d(other.wave.value,((other.noise.to(self.flux.unit)).value)**2,bounds_error=False,fill_value=0.)
        if numpy.random.choice(numpy.isfinite(self.noise.value))==True and numpy.random.choice(numpy.isfinite(other.noise.value))==True: 
            out.variance = (fself(out.wave.value)+fother(out.wave.value))*(self.flux.unit**2)
        elif numpy.random.choice(numpy.isfinite(self.noise.value))==True:
            out.variance = fself(out.wave.value)*(self.flux.unit**2)
        elif numpy.random.choice(numpy.isfinite(other.noise.value))==True:
            out.variance = fother(out.wave.value)*(self.flux.unit**2)
        else:
            out.variance = numpy.array([numpy.nan]*len(out.wave))
        out.noise = out.variance**0.5

# update other information
        out.name = self.name+' - '+other.name
        out.history.append('Difference of {} and {}'.format(self.name,other.name))
        out.original = copy.deepcopy(out)
        return out


    def __mul__(self,other):
        '''
        :Purpose: 

            A representation of multiplication for Spectrum objects which correctly interpolates 
            as a function of wavelength and combines variances

        :Output: 

            A new Spectrum object equal to the spectral product of the inputs
        '''
        try: other.wave = other.wave.to(self.wave.unit)
        except: raise ValueError('Cannot multiply spectra with wave units {} and {}'.format(self.wave.unit,other.wave.unit))

# establish the baseline wavelength grid
        out = copy.deepcopy(self)
        wave = numpy.array(copy.deepcopy(self.wave.value))
        wself = numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))
        out.wave = wave[wself]
        out.wave = out.wave*self.wave.unit
#       wother = numpy.where(numpy.logical_and(other.wave.value <= numpy.nanmax(out.wave.value),other.wave.value >= numpy.nanmin(out.wave.value)))

# do the math
        fself = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        fother = interp1d(other.wave.value,other.flux.value,bounds_error=False,fill_value=0.)
        out.flux = numpy.multiply(fself(out.wave.value),fother(out.wave.value))*self.flux.unit*other.flux.unit
        out.flux = out.flux.decompose()
# special for variance
        fselfv = interp1d(self.wave.value,(self.noise.value)**2,bounds_error=False,fill_value=0.)
        fotherv = interp1d(other.wave.value,(other.noise.value)**2,bounds_error=False,fill_value=0.)
        if numpy.random.choice(numpy.isfinite(self.noise.value))==True and numpy.random.choice(numpy.isfinite(other.noise.value))==True: 
            out.variance = numpy.multiply(out.flux.value**2,((numpy.divide(fselfv(out.wave.value),fself(out.wave.value))**2)+(numpy.divide(fotherv(out.wave.value),fother(out.wave.value))**2)))*(self.flux.unit**2)*(other.flux.unit**2)
        elif numpy.random.choice(numpy.isfinite(self.noise.value))==True:
            out.variance = numpy.multiply(out.flux.value**2,(numpy.multiply(fselfv(out.wave.value),fself(out.wave.value))**2))*(self.flux.unit**2)*(other.flux.unit**2)
        elif numpy.random.choice(numpy.isfinite(other.noise.value))==True:
            out.variance = numpy.multiply(out.flux.value**2,(numpy.multiply(flxs(out.wave.value),fother(out.wave.value))**2))*(self.flux.unit**2)*(other.flux.unit**2)
        else:
            out.variance = numpy.array([numpy.nan]*len(out.wave))
        out.variance = out.variance.decompose()
        out.noise = out.variance**0.5

# update other information
        out.name = self.name+' x '+other.name
        out.history.append('Product of {} and {}'.format(self.name,other.name))
        out.original = copy.deepcopy(out)
        return out
    

    def __div__(self,other):
        '''
        :Purpose: 

            A representation of division for Spectrum objects which correctly interpolates 
            as a function of wavelength and combines variances

        :Output: 

            A new Spectrum object equal to the spectral division of the inputs
        '''
        try: other.wave = other.wave.to(self.wave.unit)
        except: raise ValueError('Cannot divide spectra with wave units {} and {}'.format(self.wave.unit,other.wave.unit))

# establish the baseline wavelength grid
        out = copy.deepcopy(self)
        wave = numpy.array(copy.deepcopy(self.wave.value))
        wself = numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))
        out.wave = wave[wself]
        out.wave = out.wave*self.wave.unit
#       wother = numpy.where(numpy.logical_and(other.wave.value <= numpy.nanmax(out.wave.value),other.wave.value >= numpy.nanmin(out.wave.value)))

# do the math
        fself = interp1d(self.wave.value,self.flux.value,bounds_error=False,fill_value=0.)
        fother = interp1d(other.wave.value,other.flux.value,bounds_error=False,fill_value=0.)
        out.flux = numpy.divide(fself(out.wave.value),fother(out.wave.value))*self.flux.unit/other.flux.unit
        out.flux = out.flux.decompose()
# special for variance
        fselfv = interp1d(self.wave.value,(self.noise.value)**2,bounds_error=False,fill_value=0.)
        fotherv = interp1d(other.wave.value,(other.noise.value)**2,bounds_error=False,fill_value=0.)
        if numpy.random.choice(numpy.isfinite(self.noise.value))==True and numpy.random.choice(numpy.isfinite(other.noise.value))==True: 
            out.variance = numpy.multiply(out.flux.value**2,((numpy.divide(fselfv(out.wave.value),fself(out.wave.value))**2)+(numpy.divide(fotherv(out.wave.value),fother(out.wave.value))**2)))*(self.flux.unit**2)/(other.flux.unit**2)
        elif numpy.random.choice(numpy.isfinite(self.noise.value))==True:
            out.variance = numpy.multiply(out.flux.value**2,(numpy.multiply(fselfv(out.wave.value),fself(out.wave.value))**2))*(self.flux.unit**2)/(other.flux.unit**2)
        elif numpy.random.choice(numpy.isfinite(other.noise.value))==True:
            out.variance = numpy.multiply(out.flux.value**2,(numpy.multiply(flxs(out.wave.value),fother(out.wave.value))**2))*(self.flux.unit**2)/(other.flux.unit**2)
        else:
            out.variance = numpy.array([numpy.nan]*len(out.wave))
        out.variance = out.variance.decompose()
        out.noise = out.variance**0.5


# update other information
        out.name = self.name+' / '+other.name
        out.history.append('Division of {} by {}'.format(self.name,other.name))
        out.original = copy.deepcopy(out)
        return out

    def __truediv__(self,other):
        return self.__div__(other)
    
    def scale(self,fact):
        '''
        Scale spectrum by a float value
        '''
        for k in ['flux','noise']:
            if k in list(self.__dict__.keys()): setattr(self,k,getattr(self,k)*fact)
        self.variance = self.noise**2
        if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum scaled by factor {}'.format(fact))
        return

    def normalize(self,rng=[]):
        '''
        Normalize spectrum
        '''
        if len(rng)==0: 
            rng=[numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)]
        if isUnit(rng[0]): rng = [r.to(self.wave.unit).value for r in rng]
        if isUnit(rng): rng = rng.to(self.wave.unit).value
        if numpy.nanmin(rng) > numpy.nanmax(self.wave.value) or numpy.nanmax(rng) < numpy.nanmin(self.wave.value):
            print('Warning: normalization range {} is outside range of spectrum wave array ({}); ignoring'.format(rng,[numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)]))
            return
        if numpy.nanmax(rng) > numpy.nanmax(self.wave.value): rng[1] = numpy.nanmax(self.wave.value)
        if numpy.nanmin(rng) < numpy.nanmin(self.wave.value): rng[0] = numpy.nanmin(self.wave.value)

        w = numpy.where(numpy.logical_and(self.wave.value>=numpy.nanmin(rng),self.wave.value<=numpy.nanmax(rng)))
        factor = numpy.nanmax(self.flux.value[w])
        if factor == 0.: 
            print('\nWarning: normalize is attempting to divide by zero; ignoring')
        elif numpy.isnan(factor) == True: 
            print('\nWarning: normalize is attempting to divide by nan; ignoring')
        else: 
            self.scale(1./factor)
            if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum normalized')
        return

    def smooth(self,scale,method='median'):
        '''
        Do a boxcar smooth
        '''
        xsamp = numpy.arange(0,len(self.wave)-scale,scale)
        self.wave = numpy.array([self.wave.value[x+int(0.5*scale)] for x in xsamp])*self.wave.unit
        for k in ['flux','noise']:
            setattr(self,k,numpy.array([numpy.nanmedian(getattr(self,k).value[x:x+scale]) for x in xsamp])*getattr(self,k).unit)
        self.variance = self.noise**2
        if 'history' in list(self.__dict__.keys()): self.history.append('Smoothed by a scale of {} pixels'.format(scale))
    

    def convertWave(self,wave_unit):
        '''
        Convert the wavelength to a new unit
        '''
        if not isUnit(wave_unit): raise ValueError('{} is not a unit'.format(wave_unit))
        try:
            self.wave = self.wave.to(wave_unit)
            if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum wavelength converted to units of {}'.format(wave_unit))
        except:
            print('Cannot convert spectrum with wave units {} to {}'.format(self.wave.unit,wave_unit))
        return

    def convertFlux(self,flux_unit):
        '''
        Convert the flux to a new unit
        '''
        if not isUnit(flux_unit): raise ValueError('{} is not a unit'.format(flux_unit))
        for k in ['flux','noise']:
            if k in list(self.__dict__.keys()):
                try:
                    setattr(self,k,(getattr(self,k).to(flux_unit)))
                except:
                    print('Cannot convert spectrum element {} into flux units {}'.format(k,flux_unit))
        self.variance = self.noise**2
        if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum flux converted to units of {}'.format(flux_unit))
        return

    def clean(self,positive=True):
        '''
        Clean spectrum of nan, zero flux, and mask
        '''
        if 'mask' in list(self.__dict__.keys()): mask = copy.deepcopy(self.mask)
        else: mask = numpy.zeros(len(self.wave))
        mask[numpy.isnan(self.flux.value)==True] = 1
        mask[numpy.isnan(self.noise.value)==True] = 1
        if postive==True: mask[self.flux.value<=0] = 1
        for k in ['wave','flux','noise']:
            if k in list(self.__dict__.keys()):
                setattr(self,k,(getattr(self,k).value[mask==0])*getattr(self,k).unit)
        self.variance = self.noise**2
        if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum cleaned of {:.0f} pixels'.format(numpy.total(mask)))
        return

    def plot(self,**kwargs):
        '''
        Plot the spectrum
        '''
        f,fig = plt.subplots(figsize=kwargs.get('figsize',[6,4]))
        fig.plot(self.wave.value,self.flux.value,c=kwargs.get('color','k'),ls=kwargs.get('ls','-'),alpha=kwargs.get('alpha',1))
        leg = [kwargs.get('label',self.name)]
        fig.legend(leg)
        if kwargs.get('plot_uncertainty',True)==True:
            fig.plot(self.wave.value,self.noise.value,c=kwargs.get('unc_color','grey'),ls=kwargs.get('unc_ls','-'),alpha=kwargs.get('unc_alpha',1))
        if kwargs.get('plot_zero',True)==True:
            fig.plot(self.wave.value,numpy.zeros(len(self.wave)),c=kwargs.get('zero_color','k'),ls=kwargs.get('zero_ls','--'),alpha=kwargs.get('zero_alpha',1))
        fig.set_xlim(kwargs.get('xlim',[numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)]))
        fig.set_ylim(kwargs.get('ylim',[-1.1*numpy.absolute(numpy.quantile(self.flux.value,0.05)),1.5*numpy.quantile(self.flux.value,0.95)]))
        fig.set_xlabel(kwargs.get('xlabel',r'Wavelength ({})'.format(self.wave.unit)),fontsize=kwargs.get('fontsize',16))
        fig.set_ylabel(kwargs.get('xlabel',r'Flux ({})'.format(self.flux.unit)),fontsize=kwargs.get('fontsize',16))
        return fig

    def toFile(self,file,clobber=True,csv=False,delimiter='\t',save_header=True,save_noise=False,save_background=False,save_mask=False,comment='#',file_type='',**kwargs):
        '''
        Exports a spectrum to a file
        '''
# what are we saving?
        output = [self.wave.value,self.flux.value,self.noise.value]
        labels = ['Wave ({})'.format(self.wave.unit),'Flux ({})'.format(self.flux.unit),'Uncertainty ({})'.format(self.noise.unit),]
# determine which type of file
        if file_type=='': file_type = file.split('.')[-1]
# fits file
        if 'fit' in file_type:
            output = tuple(output)
            data = numpy.vstack(output)
            hdu = fits.PrimaryHDU(data)
            for k in list(self.header.keys()):
                if k.upper() not in ['HISTORY','COMMENT','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND'] and k.replace('#','') != '': # and k not in list(hdu.header.keys()):
                    hdu.header[k] = str(self.header[k])
            for k in list(self.__dict__.keys()):
                if isinstance(getattr(self,k),str) == True or isinstance(getattr(self,k),int) == True or isinstance(getattr(self,k),bool) == True or (isinstance(getattr(self,k),float) == True and numpy.isnan(getattr(self,k)) == False):
                    hdu.header[k.upper()] = str(getattr(self,k))
            hdu.writeto(file,overwrite=clobber)

# ascii file - by default tab delimited
        else:
            f = open(file,'w')
            if save_header == True:
                for k in list(self.header.keys()):
                    if k.upper() not in ['HISTORY','COMMENT'] and k.replace('#','') != '':
                        f.write('{}{} = {}\n'.format(comment,k.upper(),self.header[k]))
                for k in list(self.__dict__.keys()):
                    if isinstance(getattr(self,k),str) == True or isinstance(getattr(self,k),int) == True or isinstance(getattr(self,k),bool) == True or (isinstance(getattr(self,k),float) == True and numpy.isnan(getattr(self,k)) == False):
                        f.write('{}{} = {}\n'.format(comment,k.upper(),getattr(self,k)))
                lhead = '{}{}'.format(comment,labels[0])
                for l in labels[1:]: lhead=lhead+'{}{}'.format(delimiter,l)
                f.write('{}\n'.format(lhead))
            for i in range(len(self.wave.value)): 
                ln = '{}'.format(output[0][i])
                for j in range(1,len(labels)): ln=ln+'{}{}'.format(delimiter,output[j][i])
                f.write('{}\n'.format(ln))
            f.close()

        if 'history' in list(self.__dict__.keys()): self.history.append('Spectrum saved to {}'.format(file))
        return

    def write(self,file,**kwargs): 
        self.toFile(file,**kwargs)
        return
        
    def save(self,file,**kwargs): 
        self.toFile(file,**kwargs)
        return

    def remove(self,mask,others=[]):
        '''
        :Purpose: 

            Removes elements of wave, flux, unc arrays as specified by a  mask array

        :Required Input:

            :param mask: Either a mask array (an array of booleans, ints, or floats, 
                where True or 1 removes the element) or a 2-element array that defines 
                the wavelegnth range to mask 
        
        :Optional Input:

            :param others=[]: list of other attributes that mask should be applied to (e.g., 'pixel') 

        :Output:

            Spectrum object has the flagged pixels removed from wave, flux, unc arrays, 
            and optional arrays

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

# wavelength range given
        if len(msk) == 2 and len(self.flux) != 2:
            if not isUnit(msk): msk = msk*self.wave.unit
            msk.to(self.wave.unit)
            msk = generateMask(self.wave,mask_range=msk)

        if len(msk) != len(self.flux):
            print('\nWarning: mask must be same length ({}) as wave/flux arrays ({}); not removing any pixels'.format(len(msk),len(self.wave)))
            return
        if isinstance(msk[0],float): msk = [int(x) for x in msk]
        if isinstance(msk[0],int): msk = [True if x==1 else False for x in msk]
        if not isinstance(msk[0],bool) and not isinstance(msk[0],numpy.bool_): print('\nWarning: cannot interpret mask {}; not removing any pixels'.format(mask))

# invert and apply mask
        msk = numpy.array([not x for x in msk])
        self.wave = self.wave[msk]
        self.flux = self.flux[msk]
        self.noise = self.noise[msk]
        self.variance = self.noise**2

        if len(others) > 0:
            for k in others:
                if k in self.__dict__.keys():
                    if len(getattr(self,k)) == len(msk):
                        try: setattr(self,k,getattr(self,k)[msk])
                        except: pass

        cnt = numpy.sum([1 if x == False else 0 for x in msk])
        if 'history' in list(self.__dict__.keys()): self.history.append('Mask applied to remove {} pixels'.format(cnt))
        return


    def replace(self,mask,replace_value,replace_flux=True,replace_noise=False,others=[]):
        '''
        :Purpose: 

            Replaces flux and noise values using a mask and specified value

        :Required Inputs:

            :param mask: Either a mask array (an array of booleans, ints, or floats, where True or 1 removes the element)
                or a 2-element array that defines the wavelength range to replace 
            :param replace_value: value with which the masked elements should be replaced
        
        :Optional Inputs:

            :param replace_flux = True: replace elements in the noise array
            :param replace_noise = True: replace elements in the flux array

        :Output:

            Spectrum object has the flagged pixels replaced in flux, noise arrays, and optional vectors

        :Example:
           >>> import splat, numpy
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> num = splat.numberList('1-10,18,30-50')
           >>> mask = [not(p in num) for p in numpy.arange(len(sp.wave))]
           >>> sp.replace(mask,numpy.nan)
           >>> sp.showHistory()
                SPEX_PRISM spectrum successfully loaded
                Mask applied to replace 32 pixels with value nan
        '''
        msk = copy.deepcopy(mask)

# wavelength range given
        if len(msk) == 2 and len(self.flux) != 2:
            if not isUnit(msk): msk = msk*self.wave.unit
            msk.to(self.wave.unit)
            msk = generateMask(self.wave,mask_range=msk)

        if len(msk) != len(self.flux):
            print('\nWarning: mask must be same length ({}) as wave/flux arrays ({}); not removing any pixels'.format(len(msk),len(self.wave)))
            return
        if isinstance(msk[0],float): msk = [int(x) for x in msk]
        if isinstance(msk[0],int): msk = [True if x==1 else False for x in msk]
        if not isinstance(msk[0],bool): print('\nWarning: cannot interpret mask {}; not removing any pixels'.format(mask))

# check units of replacement value
        if not isUnit(replace_value) and not numpy.isnan(replace_value): 
            replace_value = replace_value*self.flux.unit
        if isUnit(replace_value):
            try:
                replace_value = replace_value.to(self.flux.unit)
            except:
                raise ValueError('replacement value {} does not have the same unit as flux or unc array ({})'.format(replace_value,self.flux.unit))

# invert and apply mask
        msk = numpy.array([not x for x in msk])
        if replace_flux==True: self.flux[msk] = replace_value
        if replace_noise==True: self.noise[msk] = replace_value
        if len(others) > 0:
            for k in others:
                if k in self.__dict__.keys():
                    if len(getattr(self,k)) == len(msk):
                        try: 
                            tmp = getattr(self,k)
                            tmp[msk] = replace_value
                            setattr(self,k,tmp)
                        except: pass
        self.variance = self.noise**2
        self.snr = self.computeSN()

        cnt = numpy.sum([1 if x == False else 0 for x in msk])
        if 'history' in list(self.__dict__.keys()): self.history.append('Mask applied to replace {} pixels'.format(cnt))
        return



    def clean(self,action='remove',replace_value=0.):
        '''
        :Purpose: 

            Cleans a spectrum by either removing or replacing nan values

        :Required Inputs:

            None
        
        :Optional Inputs:

            :param action = 'remove': specify as either 'remove' or 'replace'
            :param replace_value = 0.: for replace, what value to replace with

        :Output:

            Spectrum object is modified to have nan pixels "cleaned"

        :Example:
           >>> import splat,numpy
           >>> sp = splat.getSpectrum(lucky=True)[0]
           >>> sp.flux[0] = numpy.nan
           >>> sp.clean()
           >>> sp.remove(mask)
           >>> sp.showHistory()
                SPEX_PRISM spectrum successfully loaded
                Mask applied to remove 1 pixels
        '''
# clean out data points with nans in flux
        msk = [numpy.isnan(x) for x in self.flux.value]
# clean out data points with nans in noise
        msk2 = [numpy.isnan(x) for x in self.noise.value]
        msk = msk or msk2
        if action=='remove': self.remove(msk)
        elif action=='replace': self.replace(msk,replace_value)
        else: print('\nWarning: ambiguous action {} for clean; no action taken'.format(action))
        return


    def loadSourceInfo(self,verbose=True,radius=10.*u.arcsec,source_library=DB_SOURCES,**kwargs):
        '''
        WORK IN PROGRESS

        :Purpose: 
            Loads in source information for a spectrum object based on the SPLAT source database
            It uses either the spectrum object's or user-supplied NAME, DESIGNATION, COORDINATE or RA & DEC
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
            Spectrum object will have source information keywords added or updated

        THIS IS A WORK IN PROGRESS

        :Example:

        '''
        search_parameters = {
            'name': {'altname': ['SOURCE_NAME','SOURCE'], 'value': None},
            'shortname': {'altname': ['SNAME','SHORT_NAME','SHNAME'], 'value': None},
            'designation': {'altname': ['DESIG'], 'value': None},
            'coordinate': {'altname': ['COORD','COO','POSITION'], 'value': None},
            'ra': {'altname': ['RIGHT_ASCENSION','RIGHT ASCENSION','R'], 'value': None},
            'dec': {'altname': ['DECLINATION','DECL','D'], 'value': None},
        }

# feed in possible search parameters
        spattr = self.__dict__.keys()
        spattr = [k.lower() for k in spattr]
        for k in list(search_parameters.keys()):
            if k in spattr: search_parameters[k]['value'] = getattr(self,k)
            elif k in list(kwargs.keys()): search_parameters[k]['value'] = kwargs[k]
            elif k.lower() in list(kwargs.keys()): search_parameters[k]['value'] = kwargs[k.lower()]
            elif k.upper() in list(kwargs.keys()): search_parameters[k]['value'] = kwargs[k.upper()]
            else:
                for kk in search_parameters[k]['altname']:
                    if kk.lower() in spattr: search_parameters[k]['value'] = getattr(self,kk.lower())
                    elif kk.lower() in list(kwargs.keys()): search_parameters[k]['value'] = kwargs[kk.lower()]
                    elif kk.upper() in list(kwargs.keys()): search_parameters[k]['value'] = kwargs[kk.upper()]
                    else: pass

# some checking
        for k in ['name','shortname','designation']:
            if not isinstance(search_parameters[k]['value'],str): search_parameters[k]['value'] = None
        for k in ['ra','dec']:
            if not isinstance(search_parameters[k]['value'],float) and not isinstance(search_parameters[k]['value'],numpy.float64): search_parameters[k]['value'] = None
        if not isinstance(search_parameters['coordinate']['value'],SkyCoord): 
            try: search_parameters['coordinate']['value'] = properCoordinates(search_parameters['coordinate']['value'])
            except: pass
        if not isinstance(search_parameters['coordinate']['value'],SkyCoord): search_parameters['coordinate']['value'] = None

# now do search
        s = pandas.DataFrame()
        cnt = 0
        if search_parameters['name']['value'] != None:
            if verbose==True: print('Searching name = {}'.format(search_parameters['name']['value']))
            cnt+=1
            s = searchLibrary(name=search_parameters['name']['value'])
        if search_parameters['shortname']['value'] != None and len(s) == 0:
            if verbose==True: print('Searching shortname = {}'.format(search_parameters['shortname']['value']))
            cnt+=1
            s = searchLibrary(shortname=search_parameters['shortname']['value'])
        if search_parameters['designation']['value'] != None and len(s) == 0:
            if verbose==True: print('Searching designation = {}'.format(search_parameters['designation']['value']))
            cnt+=1
            s = searchLibrary(designation=search_parameters['designation']['value'])
            if len(s) == 0:
                s = searchLibrary(coordinate=search_parameters['designation']['value'],radius=radius)
        if search_parameters['coordinate']['value'] != None and len(s) == 0:
            if verbose==True: print('Searching coordinate = {}'.format(search_parameters['coordinate']['value']))
            cnt+=1
            s = searchLibrary(coordinate=search_parameters['coordinate']['value'],radius=radius)
        if search_parameters['ra']['value'] != None != False and search_parameters['dec']['value'] != None and len(s) == 0:
            if verbose==True: print('Searching ra,dec = {},{}'.format(search_parameters['ra']['value'],search_parameters['dec']['value']))
            cnt+=1
            s = searchLibrary(coordinate=properCoordinates([search_parameters['ra']['value'],search_parameters['dec']['value']]),radius=radius)
        if cnt==0:
            print('\nThere are no search parameters (name, shortname, desigation, coordinate, ra/dec) provided to search on; no source information uploaded')
            return

        if len(s) == 0:
            if verbose: print('\nNo objects found in the SPLAT source database')
        else:
            if verbose: print('\nMatched to source {}'.format(s['NAME'].iloc[0]))
            for k in list(source_library.columns):
                if k in list(s.columns):
# NOTE: CURRENTLY RETURNING FIRST OBJECT IN LIST
                    setattr(self,k.lower(),s[k].iloc[0])
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
           >>> sp = Spectrum(10001)                   # read in a source
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
        if 'history' in list(self.__dict__.keys()): self.history.append(report)

        return

    def mapTo(self,newwave,overhang=0.1):
        '''
        Purpose: maps spectrum onto the wavelength scale of another spectrum
        THIS NEEDS TO BE UPDATED

        '''
        if not isUnit(newwave): newwave = newwave*self.wave.unit
        wave_unit = self.wave.unit
        try: self.convertWave(newwave.unit)
        except: raise ValueError('Attempted to map spectrum with wavelength unit {} to wave grid with unit {}'.format(self.wave.unit,newwave.unit))
        if numpy.nanmin(newwave.value) > numpy.nanmax(self.wave.value) or numpy.nanmax(newwave.value) < numpy.nanmin(self.wave.value):
            self.convertWave(wave_unit)
            raise ValueError('New wave range {} to {}{} is outside range of spectrum {} to {}{}'.format(numpy.nanmin(newwave.value),numpy.nanmax(newwave.value),newwave.unit,numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value),self.wave.unit))

        flux_unit = self.flux.unit
        self.flux = reMap(self.wave.value,self.flux.value,newwave.value)*flux_unit
        self.noise = reMap(self.wave.value,self.noise.value,newwave.value)*flux_unit
        self.wave = newwave
        self.variance = self.noise**2
        if 'history' in list(self.__dict__.keys()): self.history.append('Mapped onto wavelength grid of {}'.format(other))
        return
        


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
        :param wave_unit: wavelength unit of final spectrum (default: wavelength unit of s1)
        :param flux_unit: flux unit of final spectrum (default: flux unit of s1)
        :param verbose: set to True to provide feedback (default: False)

    :Output: 

        New spectrum object of stitched spectrum

    :Example:
       >>> import splat
       >>> spopt = Spectrum(file='myopticalspectrum.fits')
       >>> spnir = Spectrum(file='myopticalspectrum.fits')
       >>> sp = splat.stitch(spopt,spnir,rng=[0.8,0.9],trim=[0.35,2.4])
    '''

# parameters
    scaleflag = kwargs.get('scaleflag',scale)
    vflag = kwargs.get('vflag',scale)

# generate copies of spectrum objects
    sp1 = copy.deepcopy(s1)
    sp2 = copy.deepcopy(s2)
    
# assert common units
    wave_unit = kwargs.get('wave_unit',sp1.wave.unit)
    sp1.toWaveUnit(wave_unit)
    sp2.toWaveUnit(wave_unit)
    flux_unit = kwargs.get('flux_unit',sp1.flux.unit)
    sp1.toFluxUnit(flux_unit)
    sp2.toFluxUnit(flux_unit)

# reorder to shortest wavelength spectrum first?

# check range
    rngflag = True
    if len(rng) == 0:
        rng = [numpy.nanmax([numpy.nanmin(sp1.wave.value),numpy.nanmin(sp2.wave.value)]),\
            numpy.nanmin([numpy.nanmax(sp1.wave.value),numpy.nanmax(sp2.wave.value)])]*wave_unit
    if not isUnit(rng):
        rng=rng*wave_unit
    rng = (rng.to(wave_unit)).value
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
        if verbose==True: print('Scaling over the region {} to {}'.format(rng[0],rng[1]))
        f2r = interp1d(sp2.wave.value,sp2.flux.value,bounds_error=False,fill_value=numpy.nan)
        v2r = interp1d(sp2.wave.value,sp2.variance.value,bounds_error=False,fill_value=numpy.nan)

# find overlap region, assuming first spectrum sets the flux scale standard
# assume this minimizes chi^2 residuals
        w12 = numpy.where(numpy.logical_and(numpy.array(sp1.wave.value)>=rng[0],numpy.array(sp1.wave.value)<=rng[1]))
        wv12 = numpy.array(sp1.wave.value)[w12]
        flx12 = numpy.array(sp1.flux.value)[w12]
        var12 = numpy.array(sp1.variance.value)[w12]
        sp1mid = Spectrum(wave=wv12*wave_unit,flux=flx12*flux_unit,noise=(var12**0.5)*flux_unit)
        sp2mid = Spectrum(wave=wv12*wave_unit,flux=f2r(sp1mid.wave.value)*flux_unit,noise=(v2r(sp1mid.wave.value)**0.5)*flux_unit)

        if scaleflag == True:
            chi,scl = compareSpectra(sp1mid,sp2mid)

            vtot = numpy.zeros(len(wv12))
            if not numpy.isnan(numpy.nanmedian(var12)): vtot=var12
            if not numpy.isnan(numpy.nanmedian(v2r(wv12))): vtot=vtot+v2r(wv12)

            if numpy.nanmedian(vtot) == 0.: vflag = False
            if numpy.isnan(numpy.nanmedian(vtot)) == True: vflag = False
            if vflag == True:
                scl = numpy.nansum(flx12*f2r(wv12)/vtot)/numpy.nansum(f2r(wv12)**2/vtot)
            else:
                scl = numpy.nansum(flx12*f2r(wv12))/numpy.nansum(f2r(wv12)**2)
            sp2.scale(scl)

        f2r = interp1d(sp2.wave.value,sp2.flux.value,bounds_error=False,fill_value=numpy.nan)
        v2r = interp1d(sp2.wave.value,sp2.variance.value,bounds_error=False,fill_value=numpy.nan)

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
        if numpy.isnan(numpy.nanmedian(v1))==True and numpy.isnan(numpy.nanmedian(v2))==True:
            flxmid = 0.5*(flx12+f2r(wv12))
        else:
            flxmid = (flx12/v1+f2r(wv12)/v2)/(1./v1+1./v2)
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
    wave = wave*wave_unit
    flux = flux*flux_unit
    variance = variance*(flux_unit**2)

# create new spectrum object containing combined spectrum
    sp = Spectrum(wave=wave,flux=flux,noise=variance**0.5)
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


def getSpectrum(output='spectra', limit=20, verbose=True, key_name='DATA_KEY',*args, **kwargs):
    '''
    Purpose
    -------

        Gets a spectrum from the SPLAT library using various selection criteria. 
        Calls searchLibrary_ to select spectra; if any found it routines an array of 
        Spectrum objects; otherwise, it returns an empty array. 

    Parameters
    ----------

    output = 'spectra' : str [optional]
        Determines what is returned to user; options are:
            * 'spectra': return array of Spectrum objects [default]
            * 'files': return a list of filenames
            * 'table': return pandas table of spectral information (output of `searchLibrary()`_)

    limit = 0 : int [optional]
        Sets limit for the number of returns. Set to < 0 to return unlimited number

    key = 'DATA_KEY' : str [optional]
        Keyword for the data key in the database file

    verbose = True : bool [optional]
        Set to True to provide feedback

    **kwargs : additional keywords
        Include additional keywords for `searchLibrary()`_


    Outputs
    -------

    result : array or pandas Dataframe
        An array of Spectrum objects or filenames or pandas Dataframe for sources
        that satisfy the search criteria


    Example
    -------
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
        Retrieving 1 file
    >>> sparr = splat.getSpectrum(spt='M7')
        Retrieving 120 files
    >>> sparr = splat.getSpectrum(spt='T5',young=True)
        No files match search criteria


    Dependencies
    ------------
        `searchLibrary()`_

    .. _searchLibrary : api.html#splat.core.searchLibrary

    '''
# set default search parameters - NOT NECESSARY
    # if kwargs.get('lucky',False) == True: kwargs['published'] = True
    # if kwargs.get('output','')=='': kwargs['output'] = 'all'

# search library
    search = searchLibrary(*args, **kwargs)

    if len(search) == 0: 
        if verbose==True: print('\nNo files match search criteria\n\n')
        return []
    # else:
    #     files = list(search['DATA_FILE'])
        # if len(search) == 1: files.append(search['DATA_FILE'].iloc[0])
        # else:
        #     for i,x in enumerate(search['DATA_FILE']):
        #         files.append(search['DATA_FILE'].iloc[i])

# return search table
    if 'table' in output: return search 

# return just the filenames
    if 'file' in output: return list(search['DATA_FILE']) 

# return Spectrum objects
# limit number
    if limit > 0 and limit < len(search): search = search.iloc[:limit]
    result=[]
# string on feedback
    if verbose==True:
        txt = ' file'
        if kwargs.get('lucky',False) == True: txt=' lucky'+txt
        if len(search)>1: txt = txt+'s'
        print('\nreturning {:.0f}{}\n'.format(len(search),txt))
    for i,k in enumerate(search[key_name]):
        skwargs = search.iloc[i].to_dict()
        result.append(Spectrum(k))
#         if len(files) == 1:
#             if kwargs.get('lucky',False) == True:
#                 print('\nRetrieving 1 lucky file\n')
#             else:
#                 print('\nRetrieving 1 file\n')
#             skwargs = search.iloc[0].to_dict()
#             result.append(Spectrum(files[0],**skwargs))
#         else:
# #            if (kwargs.get('lucky',False) == True):
#                print('\nRetrieving 1 lucky file\n')
#                ind = numpy.random.choice(numpy.arange(len(files)))
#                print(x)
#                result.append(Spectrum(files[ind],header=search[ind]))
#            else:
            # if limit != 0 and limit < len(files):
            #     files = files[:limit]
            #     search = search.iloc[:limit]
            # print('\nRetrieving {} files\n'.format(len(files)))
            # for i,x in enumerate(files):
            #     skwargs = search.iloc[i].to_dict()
            #     result.append(Spectrum(x,header=dict(search.iloc[i]),**skwargs))

    # else:
    #     if checkAccess() == False:
    #         sys.stderr.write('\nNo published files match search criteria\n\n')
    #     else:
    #         sys.stderr.write('\nNo files match search criteria\n\n')

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
    elif kwargs.get('dsd',False) or 'd/sd' in sptstr:
        stds = STDS_DSD_SPEX
        kys = STDS_DSD_SPEX_KEYS
        subclass = 'd/sd'
        stdtype = 'mild subdwarf'
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


def initializeStandards(*args,**kwargs):
    '''
    See initiateStandards()
    '''
    return initiateStandards(*args,**kwargs)


def initiateStandards(spt=[],dwarf=True,subdwarf=False,young=False,allstds=False,dsd=False,sd=False,esd=False,intg=False,vlg=False,reset=True,output=False,verbose=False):
    '''
    Purpose
    -------

    Initiates the spectral standards in the SpeX library. By default this loads the dwarfs standards, 
    but you can also load in subdwarf and low gravity standards as well. Once loaded, these standards remain in memory.

    Parameters
    ----------

    spt = [] : None, int, float, str, or list [optional]
        spectral type range to limit standards; if empty list or None, reads in all relevant standards

    dwarf = True : boolean [optional]
        by default, read in dwarf standards

    subdwarf = False : boolean [optional]
        read in d/sd, sd and esd standards

    young = False : boolean [optional]
        read in intg and vlg standards

    allstds = False : boolean [optional]
        read in dwarf, subdwarf, and young standards

    dsd = False : boolean [optional]
        read in mild subdwarf standards

    sd = False : boolean [optional]
        read in subdwarf standards

    esd = False : boolean [optional]
        read in extreme subdwarf standards

    int = False : boolean [optional]
        read in intermediate low gravity standards

    vlg = False : boolean [optional]
        read in very low gravity standards

    reset = False : boolean [optional]
        set to True to re-read in spectral standards (replaces those in place)

    output = False : boolean [optional]
        set to True to return the standards dictionary as output

    verbose = False : boolean [optional]
        set to True to have program return verbose output

    Outputs
    -------

    By default, code updates the relevant code dictionary variable with key (e.g., STDS_DWARF_SPEX) with value pairs
    set to spectral type: Spectrum object. If keyword output = True, will also return dictionary as output

    Example
    -------

        >>> import splat
        >>> splat.initiateStandards(spt=['L0','L5'])
        >>> splat.STDS_DWARF_SPEX

        {'L0.0': SPEX-PRISM spectrum of 2MASP J0345432+254023 (L0.0 Std),
         'L1.0': SPEX-PRISM spectrum of 2MASSW J2130446-084520 (L1.0 Std),
         'L2.0': SPEX-PRISM spectrum of Kelu-1 (L2.0 Std),
         'L3.0': SPEX-PRISM spectrum of 2MASSW J1506544+132106 (L3.0 Std),
         'L4.0': SPEX-PRISM spectrum of 2MASS J21580457-1550098 (L4.0 Std),
         'L5.0': SPEX-PRISM spectrum of SDSS J083506.16+195304.4 (L5.0 Std)}

        >>> splat.initiateStandards(spt=['L0','L5'],intg=True)
        >>> splat.STDS_INTG_SPEX

        {'L0.0beta': SPEX-PRISM spectrum of 2MASSW J1552591+294849 (L0.0beta Std),
         'L1.0beta': SPEX-PRISM spectrum of 2MASS J02271036-1624479 (L1.0beta Std),
         'L2.0beta': SPEX-PRISM spectrum of LSR 0602+3910 (L2.0beta Std),
         'L3.0beta': SPEX-PRISM spectrum of 2MASSI J1726000+153819 (L3.0beta Std)}

    Dependencies
    ------------
        
        `typeToNum()`_

    .. _`typeToNum()` : api.html#splat.core.typeToNum

    '''


# choose what kind of standards desired - d, sd, esd
# and read in standards into dictionary if they haven't been already
    if allstds==True:
        if verbose==True: print('Reading in all standards')
        initiateStandards(spt=spt,reset=reset,verbose=verbose,output=False)
        initiateStandards(spt=spt,subdwarf=True,reset=reset,verbose=verbose,output=False)
        initiateStandards(spt=spt,young=True,reset=reset,verbose=verbose,output=False)
        if output==True:
            stds = copy.deepcopy(STDS_DWARF_SPEX)
            for x in [STDS_SD_SPEX_KEYS,STDS_DSD_SPEX_KEYS,STDS_ESD_SPEX_KEYS,STDS_VLG_SPEX_KEYS,STDS_INTG_SPEX_KEYS]:
                stds.update(x)
            return stds
        else: return
    if subdwarf==True:
        if verbose==True: print('Reading in all metal-poor standards')
        initiateStandards(spt=spt,sd=True,reset=reset,verbose=verbose,output=False)
        initiateStandards(spt=spt,dsd=True,reset=reset,verbose=verbose,output=False)
        initiateStandards(spt=spt,esd=True,reset=reset,verbose=verbose,output=False)
        if output==True:
            stds = copy.deepcopy(STDS_SD_SPEX_KEYS)
            for x in [STDS_DSD_SPEX_KEYS,STDS_ESD_SPEX_KEYS]:
                stds.update(x)
            return stds
        else: return
    if young==True:
        if verbose==True: print('Reading in all young standards')
        initiateStandards(spt=spt,intg=True,reset=reset,verbose=verbose,output=False)
        initiateStandards(spt=spt,vlg=True,reset=reset,verbose=verbose,output=False)
        if output==True:
            stds = copy.deepcopy(STDS_VLG_SPEX_KEYS)
            stds.update(STDS_INTG_SPEX_KEYS)
            return stds
        else: return

# assign dictionary to fill
    if sd==True:
        if verbose==True: print('Reading in subdwarf (sd) standards')
        stds = STDS_SD_SPEX
        kys = copy.deepcopy(STDS_SD_SPEX_KEYS)
    elif dsd==True:
        if verbose==True: print('Reading in mild subdwarf (d/sd) standards')
        stds = STDS_DSD_SPEX
        kys = copy.deepcopy(STDS_DSD_SPEX_KEYS)
    elif esd==True:
        if verbose==True: print('Reading in extreme subdwarf (esd) standards')
        stds = STDS_ESD_SPEX
        kys = copy.deepcopy(STDS_ESD_SPEX_KEYS)
    elif vlg==True:
        if verbose==True: print('Reading in very low gravity (vlg) standards')
        stds = STDS_VLG_SPEX
        kys = copy.deepcopy(STDS_VLG_SPEX_KEYS)
    elif intg==True:
        if verbose==True: print('Reading in intermediate low gravity (intg) standards')
        stds = STDS_INTG_SPEX
        kys = copy.deepcopy(STDS_INTG_SPEX_KEYS)
    elif dwarf==True:
        if verbose==True: print('Reading in dwarf standards')
        stds = STDS_DWARF_SPEX
        kys = copy.deepcopy(STDS_DWARF_SPEX_KEYS)
    else:
        print('WARNING: now spectral classes indicated to read in, returning empty dictionary')
        return

# read in spectral types as desired by comparing numerical 
    if not isinstance(spt,type(None)) and (isinstance(spt,list) and len(spt)>0):
        newkys = {}
        sptset = copy.deepcopy(spt)
        if not isinstance(sptset,list): sptset = [sptset]
        if isinstance(sptset[0],str): sptset = [typeToNum(x) for x in sptset]
# single type
        if len(sptset)==1:
            for k in kys:
                if typeToNum(k)==sptset[0]: newkys[k] = kys[k]
# range min to max inclusive
        elif len(sptset)==2:
            for k in kys:
                if sptset[0]<=typeToNum(k)<=sptset[1]: newkys[k] = kys[k]
# individual types
        else:
            for t in sptset:
                for k in keys:
                    if typeToNum(k)==t: newkys[k] = kys[k]
        kys = copy.deepcopy(newkys)

# now read in missing keys as needed
    for t in list(kys.keys()):
        if t not in list(stds.keys()) or reset==True:
            stds[t] = Spectrum(kys[t])
            stds[t].normalize()
#            stds[t].name += ' ({} Std)'.format(t)
            if verbose==True: print('\tReading in {}'.format(stds[t].name))

    if output==True: return stds
    return


def keySource(keys,verbose=True):
    '''
    Purpose
    -------

    Searches DB_SOURCE for source matching input source key(s) in column 'SOURCE_KEY'

    Parameters
    ----------

    keys : int or list
        integer or list of integers corresponding to source keys

    verbose = True : boolean [optional]
        set to True to have program return verbose output

    Outputs
    -------

    pandas DataFrame containing the rows in DB_SOURCE that match input keys, 
    or empty pandas DataFrame

    Example
    -------

        TBD

    Dependencies
    ------------
        
        `keySearch()`_

    .. _`keySearch()` : api.html#splat.utilities.keySearch

    '''
    return keySearch(keys,key_name='SOURCE_KEY',db=DB_SOURCE,verbose=verbose)


def keySpectrum(keys,verbose=True):
    '''
    Purpose
    -------

    Searches DB_SPECTRA for source matching input source key(s) in column 'DATA_KEY'

    Parameters
    ----------

    keys : int or list
        integer or list of integers corresponding to source keys

    verbose = True : boolean [optional]
        set to True to have program return verbose output

    Outputs
    -------

    pandas DataFrame containing the rows in DB_SPECTRA that match input keys, 
    or empty pandas DataFrame

    Example
    -------

        TBD

    Dependencies
    ------------
        
        `keySearch()`_

    .. _`keySearch()` : api.html#splat.utilities.keySearch

    '''
    return keySearch(keys,key_name='DATA_KEY',db=DB_SPECTRA,verbose=verbose)

# def keySource(keys, verbose=True db=DB_SOURCES, **kwargs):
#     '''
#     :Purpose: Takes a source key and returns a table with the source information
#     :param keys: source key or a list of source keys
#     :Example:
#     >>> import splat
#     >>> print spl.keySource(10001)
#         SOURCE_KEY           NAME              DESIGNATION    ... NOTE SELECT
#         ---------- ------------------------ ----------------- ... ---- ------
#              10001 SDSS J000013.54+255418.6 J00001354+2554180 ...        True
#     >>> print spl.keySource([10105, 10623])
#         SOURCE_KEY          NAME             DESIGNATION    ... NOTE SELECT
#         ---------- ---------------------- ----------------- ... ---- ------
#              10105 2MASSI J0103320+193536 J01033203+1935361 ...        True
#              10623 SDSS J09002368+2539343 J09002368+2539343 ...        True
#     >>> print spl.keySource(1000001)
#         No sources found with source key 1000001
#         False
#     '''

# # vectorize
#     if isinstance(keys,list) == False: keys = [keys]

# #    sdb = ascii.read(DB_FOLDER+SOURCES_DB, delimiter='\t',fill_values='-99.',format='tab')
# #    sdb = fetchDatabase(DB_FOLDER+SOURCES_DB)
# #    sdb = copy.deepcopy(DB_SOURCES)
# #    sdb['SELECT'] = [x in keys for x in sdb['SOURCE_KEY']]

#     sdb = db[[x in keys for x in db[key_name]]]
#     if len(sdb) == 0.:
#         if kwargs.get('verbose',True) == True: print('No sources found with source key(s) = {}'.format(*keys))
#         return False
#     else:
# #        db = sdb[:][numpy.where(sdb['SELECT']==1)]
#         return sdb


# def keySpectrum(keys, **kwargs):
#     '''
#     :Purpose: Takes a spectrum key and returns a table with the spectrum and source information
#     :param keys: spectrum key or a list of source keys
#     :Example:
#     >>> import splat
#     >>> print spl.keySpectrum(10001)
#         DATA_KEY SOURCE_KEY    DATA_FILE     ... COMPANION COMPANION_NAME NOTE_2
#         -------- ---------- ---------------- ... --------- -------------- ------
#            10001      10443 10001_10443.fits ...
#     >>> print spl.keySpectrum([10123, 11298])
#         DATA_KEY SOURCE_KEY    DATA_FILE     ... COMPANION COMPANION_NAME NOTE_2
#         -------- ---------- ---------------- ... --------- -------------- ------
#            11298      10118 11298_10118.fits ...
#            10123      10145 10123_10145.fits ...
#     >>> print spl.keySpectrum(1000001)
#         No spectra found with spectrum key 1000001
#         False
#     '''

# # vectorize
#     if isinstance(keys,list) == False:
#         keys = [keys]

# #    sdb = copy.deepcopy(DB_SPECTRA)
# #    sdb['SELECT'] = [x in keys for x in sdb['DATA_KEY']]

# #    if sum(sdb['SELECT']) == 0.:
# #        if verbose: print('No spectra found with spectrum key {}'.format(keys[0]))
# #        return False
# #    else:
# #        s2db = copy.deepcopy(DB_SOURCES)
# #        db = join(sdb[:][numpy.where(sdb['SELECT']==1)],s2db,keys='SOURCE_KEY')
# #        return db

#     sdb = DB_SPECTRA[[x in keys for x in DB_SPECTRA['DATA_KEY']]]
#     if len(sdb) == 0.:
#         if kwargs.get('verbose',True) == True: print('No sources found with spectrum key(s) = {}'.format(*keys))
#         return False
#     else:
#         return sdb


def searchLibrary(radius=10., instrument='SPEX-PRISM', source_database=DB_SOURCES, spectra_database=DB_SPECTRA, *args, **kwargs):
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

# check the input catalog
#    source_db = ascii.read(DB_FOLDER+SOURCES_DB, delimiter='\t', fill_values='-99.', format='tab')
#    source_db = fetchDatabase(SOURCES_DB)
    source_db = copy.deepcopy(source_database)
    spectra_db = copy.deepcopy(spectra_database)
    required_keywords = ['SOURCE_KEY']
    for k in required_keywords:
        if k not in list(source_db.keys()): raise ValueError('Source database must include keyword {}'.format(k))

# program parameters
    ref = kwargs.get('output','all')
    verbose = kwargs.get('verbose',False)

# prep object classes
    object_classes = numpy.array([])
    if 'OBJECT_TYPE' in list(source_db.keys()):
        object_classes = numpy.unique(numpy.sort(numpy.array([str(x) for x in DB_SOURCES['OBJECT_TYPE']])))
    object_classes = object_classes[numpy.where(object_classes != '0')]  # eliminate masked element
    object_classes = object_classes[numpy.where(object_classes != 'nan')]  # eliminate masked element

# logic of search
    logic = 'and'         # default combination
    logic = kwargs.get('combine',logic).lower()
    logic = kwargs.get('logic',logic).lower()
    if (logic != 'and' and logic != 'or'):
        raise ValueError('\nLogical operator '+logic+' not supported\n\n')

# first search by source parameters
    source_db['SELECT'] = numpy.zeros(len(source_db['SOURCE_KEY']))
    count = 0.

# search by source key
    idkey = False
    for k in ['source_key','source','idkey','id_key','id']: idkey = kwargs.get(k,idkey)
    if idkey != False:
        if not isinstance(idkey,list): idkey = [idkey]
        if isinstance(idkey[0],str): idkey = [int(i) for i in idkey]
        for s in idkey: source_db['SELECT'][source_db['SOURCE_KEY'] == s] += 1
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
    sname = kwargs.get('reject_shortname',sname)
    if sname != False and len(sname) > 0:
        if 'SHORTNAME' not in source_db.keys():
            source_db['SHORTNAME'] = [designationToShortName(x) for x in source_db['DESIGNATION']]
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
#            print(radius,numpy.nanmin(r),source_db['DESIGNATION'].iloc[numpy.argmin(r)])
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
#    spectral_db = ascii.read(DB_FOLDER+SPECTRA_DB, delimiter='\t',fill_values='-99.',format='tab')
#    spectral_db = fetchDatabase(DB_FOLDER+SPECTRA_DB)
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
            spectral_db['SELECT'][spectral_db['BIBCODE'] == r] += 1
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
#    if (not checkAccess() or kwargs.get('published',False) or kwargs.get('public',False)):
#        spectral_db['SELECT'][spectral_db['PUBLISHED'] != 'Y'] = 0.

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



def readSpectrum(file,folder='',file_type='',wave_unit=DEFAULT_WAVE_UNIT,
    flux_unit=DEFAULT_FLUX_UNIT,dimensionless=False,comment='#',delimiter='',
    crval1='CRVAL1',cdelt1='CDELT1',catch_sn=True,remove_nans=True,no_zero_noise=True,
    use_instrument_reader=True,instrument=DEFAULT_INSTRUMENT,instrument_param={},verbose=False):
    '''
    Purpose
    -------

    Reads in spectral data from a variety of formats

    Parameters
    ----------

    file : string
        filename of data to be read in; if full path not provided, routine will search in local directory
        or in folder indicated by folder keyword; can also pass a URL to a remote file

    folder = '' : string [optional]
        full path to folder containing file

    file_type = '' : string [optional]
        a string indicating some flags that specify what kind of file this is; some examples include:
            * csv = comma separated file (sets delimiter = ',')
            * tab or tsv = tab delimited file (sets delimiter = '\t')
            * pipe = pipe delimited file (sets delimiter = '|')
            * latex = latex-style table data (sets delimiter = ' & ')
            * waveheader = wavelength solution is contained within the fits header
            * wavelog = wavelength solution is logarithmic (common for echelle data)
            * sdss = sets both waveheader and wavelog

    wave_unit = DEFAULT_WAVE_UNIT : astropy.unit [optional]
        units of wavelength axis, by default specified by the DEFAULT_WAVE_UNIT global parameter

    flux_unit = DEFAULT_FLUX_UNIT : astropy.unit [optional]
        units of flux and uncertainty axes, by default specified by the DEFAULT_FLUX_UNIT global parameter
        note that you can specify a unitless number 

    dimensionless = False : boolean [optional]
        set to True to set the flux units to a dimensionless quantity (for transmission, reflectance)

    comment = '#' : string [optional]
        for ascii files, character that indicates the file line is a comment (to be ignored)

    delimiter = '' : string [optional]
        for ascii files, character that separates columns of values

    crval1,cdelt1 = 'CRVAL1','CDELT1' : string [optional]
        for fits files for which the wavelength solution is embedded in header, 
        these are the keywords containing the zeroth wavelength and linear change coefficient

    catch_sn = True : boolean [optional]
        set to True to check if uncertainty axis is actually signal-to-noise, 
        by checking if median(flux/uncertainty) < 1 [NOTE: NO LONGER IMPLEMENTED]

    remove_nans = True : boolean [optional]
        set to True to remove all wave/flux/noise values for which wave or flux are nan

    no_zero_noise = True : boolean [optional]
        set to True to set all elements of noise array that are zero to numpy.nan;
        this helps in later computations of S/N or fit statistics

    use_instrument_reader = True : bool [optional]
        set to True to use the default instrument read in the DEFAULT-INSTRUMENT global parameter
        checked against INSTRUMENTS dictionary

    instrument = DEFAULT_INSTRUMENT : string [optional]
        instrument by which data was acquired, by default the DEFAULT-INSTRUMENT global parameter
        checked against INSTRUMENTS dictionary

    instrument_param = {} : dict [optional]
        instrument-specific parameters to pass to instrument reader

    verbose = False : boolean [optional]
        set to True to have program return verbose output

    Outputs
    -------

    output : dictionary
        dictionary containing the following elements:
            * wave : numpy.array = array of wavelength values (not unitted)
            * flux : numpy.array = array of flux values (not unitted)
            * noise : numpy.array = array of uncertainty values (not unitted)
            * wave_unit : astropy.unit = units of wavelength axis
            * flux_unit : astropy.unit = units of flux and uncertainty axes
            * header : dictionary = header of file (if fits)

    Example
    -------

       Case 1: An ascii file named 'spectrum.csv' with three comma-delimited columns of wavelength in Angstroms 
       and flux and uncertainty in erg/s/cm2/Angstrom

           >>> import splat
           >>> import astropy.units as u
           >>> splat.readSpectrum('spectrum.csv',wave_unit=u.Angstrom,flux_unit=u.erg/u.s/u.cm/u.cm/u.Angstrom)

            {'wave': <Quantity [0.45, 0.5 , 0.55, 0.6 , 0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 ,
                        0.95, 1.  , 1.05, 1.1 , 1.15, 1.2 , 1.25, 1.3 , 1.35, 1.4 ,
                        1.45, 1.5 , 1.55, 1.6 , 1.65, 1.7 , 1.75, 1.8 , 1.85, 1.9 ,
                        1.95, 2.  , 2.05, 2.1 , 2.15, 2.2 , 2.25, 2.3 , 2.35, 2.4 ,
                        2.45] Angstrom>,
             'flux': <Quantity [0.965, 0.985, 1.   , 1.007, 1.012, 1.014, 1.014, 1.012, 1.009,
                        1.005, 1.003, 1.004, 1.004, 1.005, 1.008, 1.012, 1.015, 1.019,
                        1.026, 1.034, 1.043, 1.052, 1.062, 1.071, 1.081, 1.09 , 1.101,
                        1.109, 1.12 , 1.13 , 1.139, 1.148, 1.158, 1.17 , 1.18 , 1.188,
                        1.196, 1.202, 1.208, 1.214, 1.22 ] erg / (Angstrom cm2 s)>,
             'noise': <Quantity [0.022, 0.01 ,   nan, 0.01 , 0.015, 0.018, 0.019, 0.02 , 0.024,
                        0.025, 0.028, 0.032, 0.036, 0.04 , 0.044, 0.048, 0.05 , 0.053,
                        0.055, 0.057, 0.061, 0.064, 0.068, 0.07 , 0.074, 0.077, 0.081,
                        0.083, 0.086, 0.09 , 0.092, 0.095, 0.099, 0.102, 0.107, 0.109,
                        0.112, 0.115, 0.116, 0.118, 0.12 ] erg / (Angstrom cm2 s)>,
             'header': ,
             'wave_unit': Unit("Angstrom"),
             'flux_unit': Unit("erg / (Angstrom cm2 s)")}

       Case 2: A fits file named 'spectrum.fits' that is a asteroid spectrum, wavelength in micron 
       and flux and uncertainty in dimensionless reflectance

           >>> import splat
           >>> import astropy.units as u
           >>> splat.readSpectrum('asteroid.fits',wave_unit=u.micron,dimensionless=True)

        {'wave': <Quantity [0.65777063, 0.65988046, 0.66198677, 0.66409034, 0.6661918 ,
                    0.6682921 , 0.67039174, 0.6724916 , 0.6745923 , 0.6766945 ,
                    0.6787989 , 0.6809061 , 0.68301684, 0.6851316 , 0.68725103,
                    0.68937576, 0.69150627, 0.6936432 , 0.6957871 , 0.6979384 ,
                    0.7000977 , 0.70226544, 0.70444214, 0.7066283 , 0.7088244 ,
                    ...
                    2.5862362 , 2.5901842 , 2.594177  ] micron>,
         'flux': <Quantity [0.2733175 , 0.35201955, 0.35244268, 0.23883328, 0.30472383,
                    0.3025949 , 0.21126787, 0.22968623, 0.3400901 , 0.18981682,
                    0.36001053, 0.23199332, 0.32443064, 0.24086563, 0.18240955,
                    0.25307003, 0.20817572, 0.20775121, 0.06874909, 0.9939198 ,
                    0.23347081, 0.22479972, 0.2291573 , 0.2468396 , 0.22524261,
                    ...
                    0.28856662, 0.2951897 , 0.29580346]>,
         'noise': <Quantity [1.94668144e-01, 3.19979459e-01, 2.37098530e-01, 1.75096914e-01,
                    2.45749980e-01, 2.58316606e-01, 1.30246907e-01, 2.03696743e-01,
                    1.72521144e-01, 8.74038711e-02, 1.51685730e-01, 1.09243065e-01,
                    8.36683586e-02, 7.38894045e-02, 4.08149362e-02, 4.65220883e-02,
                    3.61229926e-02, 4.17911708e-02, 7.78270513e-02, 1.64021879e-01,
                    ...
                    6.11881772e-03, 6.60514692e-03, 8.59303772e-03]>,
         'header': SIMPLE  =                    T / Written by IDL:  Tue Jun 29 15:51:22 2021      
         BITPIX  =                  -32 / IEEE single precision floating point           
         NAXIS   =                    2 /                                                
         NAXIS1  =                  564 /                                                
         NAXIS2  =                    4 /                                                
         ...                                                                                
         HISTORY com_spec_2122-2129.fits. Therefore the spectrum is unitless.            ,
         'wave_unit': Unit("micron"),
         'flux_unit': Unit(dimensionless)}

    Note that `readSpectrum()` is called by the Spectrum class; you can also create a Spectrum object directly
    from the output of `readSpectrum()`:

        >>> d = splat.readSpectrum('asteroid.fits',wave_unit=u.micron,dimensionless=True)
        >>> sp = Spectrum(d,name='Asteroid',instrument='SPEX-PRISM',flux_label='Reflectance')
        >>> sp.info()

        IRTF SpeX prism spectrum of Asteroid

        UNPUBLISHED DATA

        History:
            SPEX-PRISM spectrum successfully loaded

    Dependencies
    ------------
        
        `checkAccess()`_
        `checkInstrument()`_
        `checkOnline()`_
        `isUnit()`_
        `_readAPOGEE()`_
        `_readBOSS()`_
        `_readIRAF()`_
        `_readFIRE()`_
        `_readMAGE()`_
        `_readWFC3()`_
        `_readKAST()`_
        astropy.io.fits
        gzip
        os
        shutil

    .. _`checkInstrument()` : api.html#splat.utilities.checkInstrument
    .. _`checkAccess()` : api.html#splat.utilities.checkAccess
    .. _`checkOnline()` : api.html#splat.utilities.checkOnline
    .. _`isUnit()` : api.html#splat.utilities.isUnit
    .. _`_readAPOGEE()` : api.html#splat.core._readAPOGEE
    .. _`_readBOSS()` : api.html#splat.core._readBOSS
    .. _`_readIRAF()` : api.html#splat.core._readIRAF
    .. _`_readFIRE()` : api.html#splat.core._readFIRE
    .. _`_readMAGE()` : api.html#splat.core._readMAGE
    .. _`_readWFC3()` : api.html#splat.core._readWFC3
    .. _`_readKAST()` : api.html#splat.core._readKAST

    '''
# check inputs and keyword parameters
    # for k in ['file','filename','data_file','datafile']:
    #     if k in list(kwargs.keys()): file = kwargs[k]
    if file == '': raise NameError('\nNo filename passed to readSpectrum')
    if not(isUnit(wave_unit)):
        if verbose==True: print('Warning: wave_unit {} is not an astropy unit; using default {}'.format(wave_unit,DEFAULT_WAVE_UNIT))
        wave_unit = DEFAULT_WAVE_UNIT
    if not(isUnit(flux_unit)):
        if verbose==True: print('Warning: flux_unit {} is not an astropy unit; using default {}'.format(flux_unit,DEFAULT_FLUX_UNIT))
        flux_unit = DEFAULT_FLUX_UNIT
    if dimensionless==True: flux_unit = u.dimensionless_unscaled

# if a url, make sure it exists
    if file[:4]=='http':
        if requests.get(file).status_code!=requests.codes.ok:
            raise ValueError('Cannot find remote file {}; check URL or your online status'.format(file))

# if a local file, make sure it exists
    else:
        f0 = copy.deepcopy(file)
        if os.path.exists(os.path.normpath(file)) == False:
            file = '{}.{}'.format(file,'gz') 
        if os.path.exists(os.path.normpath(file)) == False:
            file = os.path.join(folder,os.path.basename(f0))
        if os.path.exists(os.path.normpath(file)) == False:
            file = '{}.{}'.format(file,'gz') 
        if os.path.exists(os.path.normpath(file)) == False:
            raise ValueError('Cannot find {} locally or in folder {}\n\n'.format(f0,folder))

# determine which type of file
    file_type = '{} {}'.format(file.split('.')[-1],file_type)

# gzipped file - extract root
#    for k in ['gz','bz2','zip']:
    gzflag = False
    if 'gz' in file_type:
        file_type = '{} {}'.format((file.replace('.gz','')).split('.')[-1],file_type)
        gzflag = True

# instrument specific read shortcut - not working?
    readin = False
    inst = checkInstrument(instrument)
#    if inst != False and use_instrument_reader==True: 
        # this mode is not working
        # instrument = inst
        # if INSTRUMENTS[instrument]['reader'] != '': 
        #     output = INSTRUMENTS[instrument]['reader'](file,verbose=verbose)
        #     readin = True

# other reads
    if readin==False:

# fits - can be done with fits.open as local or online and w/ or w/o gzip/bzip2/pkzip
        if 'fit' in file_type:
            with fits.open(os.path.normpath(file),ignore_missing_end=True,ignore_missing_simple=True,do_not_scale_image_data=True) as hdu:
                hdu.verify('silentfix+ignore')
                header = hdu[0].header
                if 'NAXIS3' in list(header.keys()): d = numpy.copy(hdu[0].data[0,:,:])
                else: d =  numpy.copy(hdu[0].data)
# make sure file is oriented correctly
            if numpy.shape(d)[0]>numpy.shape(d)[1]: d = numpy.transpose(d)

# some specific formatting cases
            if 'sdss' in file_type and 'fit' in file_type: file_type='waveheader wavelog {}'.format(file_type)
# wavelength is in header 
            if 'waveheader' in file_type and 'fit' in file_type and len(d[:,0])<3:
                flux = d[0,:]
                if crval1 in list(header.keys()) and cdelt1 in list(header.keys()):
                    wave = numpy.polyval([float(header[cdelt1]),float(header[crval1])],numpy.arange(len(flux)))
                else: 
                    raise ValueError('\nCannot find {} and {} keywords in header of fits file {}'.format(crval1,cdelt1,file))
# wavelength is explicitly in data array 
            else:
                wave = d[0,:]
                flux = d[1,:]
            if len(d[:,0]) > 2: noise = d[2,:]
            else: noise = [numpy.nan]*len(wave)

# ascii - can be done with pandas as local or online and w/ or w/o gzip/bzip2/pkzip
        else:
            if 'csv' in file_type and delimiter=='': delimiter = ','
            elif ('tsv' in file_type or 'tab' in file_type or 'txt' in file_type) and delimiter=='': delimiter = '\t'
            elif 'pipe' in file_type and delimiter=='': delimiter = '|'
            elif 'latex' in file_type and delimiter=='': delimiter = '&'

# initial read
            dp = pandas.read_csv(file,delimiter=delimiter,comment=comment,header=0)
# if numbers in first row, replace with header
            if isNumber(dp.columns[0])==True:
                cnames = ['wave','flux']
                if len(dp.columns)>2: cnames.append('noise')
                if len(dp.columns)>3: 
                    for i in range(len(dp.columns))-3: cnames.append('c{}'.format(i))
                dp = pandas.read_csv(file,delimiter=delimiter,comment=comment,names=cnames)
# assume order wave, flux, noise
# for now ignoring any other columns
            wave = numpy.array(dp[dp.columns[0]])
            flux = numpy.array(dp[dp.columns[1]])
            if len(dp.columns)>2: noise = numpy.array(dp[dp.columns[2]])
            else: noise = [numpy.nan]*len(dp)
# placeholder header
            header = fits.Header()      # blank header

#  wavelength scale is logarithmic
        if 'wavelog' in file_type: wave = 10.**wave

# final output dictionary
        output = {'wave': wave, 
                  'flux': flux, 
                  'noise': noise, 
                  'header': header, 
                  'wave_unit': wave_unit, 
                  'flux_unit': flux_unit}

# make sure arrays are numpy arrays
    output['wave'] = numpy.array(output['wave'])
    output['flux'] = numpy.array(output['flux'])
    output['noise'] = numpy.array(output['noise'])

# make sure arrays have units
    if not isUnit(output['wave']): output['wave'] = output['wave']*wave_unit
    output['wave'].to(wave_unit)
    if not isUnit(output['flux']): output['flux'] = output['flux']*flux_unit
    output['flux'].to(flux_unit)
    if not isUnit(output['noise']): output['noise'] = output['noise']*flux_unit
    output['noise'].to(flux_unit)

# remove all parts of spectrum that are nans
    if remove_nans==True:
        w = numpy.where(numpy.logical_and(numpy.isnan(output['wave']) == False,numpy.isnan(output['flux']) == False))
        output['wave'] = output['wave'][w]
        output['flux'] = output['flux'][w]
        output['noise'] = output['noise'][w]

# force places where noise is zero to be NaNs
    if no_zero_noise==True:
        output['noise'][numpy.where(output['noise'] == 0.)] = numpy.nan


# add in instrument specific information
# this includes default wavelength and flux units
# LEAVING THIS OUT FOR THE MOMENT
#    if inst != False:
#        for k in list(INSTRUMENTS[inst].keys()): output[k] = INSTRUMENTS[inst][k]  

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
#    print(kwargs,model)
    if kwargs.get('apstar',False) == True: model='apstar'
    if kwargs.get('apvisit',False) == True: model='apvisit'
    if kwargs.get('aspcap',False) == True: model='aspcap'
    if kwargs.get('ap1d',False) == True: model='ap1d'

    hdulist = fits.open(file)
    header = hdulist[0].header
    output = {'wave_unit': u.Angstrom,'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom,'header': header}
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
          'wave_unit': u.Angstrom,
          'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom}

def _readKAST(file,**kwargs):
    '''
    Reads kastredux KAST fits files
    '''
    wave_unit = INSTRUMENTS['KAST-RED']['wave_unit']
    flux_unit = INSTRUMENTS['KAST-RED']['flux_unit']

    if not os.path.exists(file):
        raise NameError('\nCould not find KAST file {}'.format(file))

    with fits.open(os.path.normpath(file),ignore_missing_end=True,ignore_missing_simple=True,do_not_scale_image_data=True) as data:
        data.verify('silentfix+ignore')
        if 'NAXIS3' in list(data[0].header.keys()):
            d = numpy.copy(data[0].data[0,:,:])
        else:
            d =  numpy.copy(data[0].data)
        header = data[0].header
    if header['NAXIS']==1:
        w0,w1 = None,None
        for x in ['CRVAL1']:
            if x in list(header.keys()): w0 = header[x]
        for x in ['CDELT1','CD1_1']:
            if x in list(header.keys()): w1 = header[x]
        if w0==None or w1==None: raise ValueError('Could not find the necessary keywords in fits header to generate wavelength scale')
        if 'linear' in header['CTYPE1'].lower(): wave = numpy.arange(len(d))*w1+w0
        else: raise ValueError('Cannot handle wavelength scale type {}'.format(header['CTYPE1']))
        wave=wave*wave_unit
        flux = d*flux_unit
    else:
        wave=d[0,:]*wave_unit
        flux=d[1,:]*flux_unit
        if len(d[:,0])>2: noise=d[2,:]*flux_unit
    
    return {'wave': wave,
          'flux': flux,
          'noise': noise,
          'header': header,
          'wave_unit': wave_unit,
          'flux_unit': flux_unit}
INSTRUMENTS['KAST-RED']['reader'] = _readKAST
INSTRUMENTS['KAST-BLUE']['reader'] = _readKAST


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
            'wave_unit': u.Angstrom, 
            'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom}


def _readMAGE(file, **kwargs):
    '''
    Reads MagE spectral files as output from MASE.    
    This data don't have wavelength in the data, but the wavelength is described in the header.
    There may be some data from Mage that have 3 columns that include wavelength, flux, and noise.
    '''

    if not os.path.exists(file):
        raise NameError('\nCould not find MagE file {}'.format(file))
    wave_unit = u.Angstrom
    flux_unit = u.erg/u.s/u.cm/u.cm/u.Angstrom

    hdulist = fits.open(file)
    header = hdulist[0].header
    if len(hdulist[0].data) == 2:
        wave = (10.**(numpy.linspace(hdulist[0].header['CRVAL1'],hdulist[0].header['CRVAL1']+hdulist[0].header['NAXIS1']*hdulist[0].header['CDELT1'],num=hdulist[0].header['NAXIS1'],endpoint=False)))*wave_unit
        flux = hdulist[0].data[0]*flux_unit
        noise = hdulist[0].data[1]*flux_unit
    elif len(hdulist[0].data) == 3:
        wave = hdulist[0].data[0]*wave_unit
        flux = hdulist[0].data[1]*flux_unit
        noise = hdulist[0].data[2]*flux_unit
    else:
        raise ValueError('\nWas expecting 2 or 3 data axes; instead found {}'.format(len(hdulist[0].data)))

    return {'wave': wave, 
            'flux': flux, 
            'noise': noise, 
            'header': header, 
            'wave_unit': wave_unit,
            'flux_unit': flux_unit}


def _readFIRE(file,**kwargs):
    '''
    Reads FIRE spectral files as output from FIREHOSE.
    
    '''
    if not os.path.exists(file):
        raise NameError('\nCould not find FIRE file {}'.format(file))
    wave_unit = u.micron
    flux_unit = u.erg/u.s/u.cm/u.cm/u.micron

    hdulist = fits.open(file)
    if len(hdulist[0].data) != 3:
        raise ValueError('\n_readFIRE can only read data reduced via FIREHOSE, with wave, flux, noise in fits data channel')
    header = hdulist[0].header
    wave = hdulist[0].data[0]*wave_unit
    flux = hdulist[0].data[1]*flux_unit
    noise = hdulist[0].data[2]*flux_unit
    if kwargs.get('gagne',False) == True or numpy.median(wave.value) > 3.:
        wave = (hdulist[0].data[0]*u.Angstrom).to(wave_unit)
        flux = (hdulist[0].data[1]*u.erg/u.s/u.cm/u.cm/u.Angstrom).to(flux_unit)
        noise = (hdulist[0].data[2]*u.erg/u.s/u.cm/u.cm/u.Angstrom).to(flux_unit)
    
    return {'wave': wave, 
            'flux': flux, 
            'noise': noise, 
            'header': header, 
            'wave_unit': wave_unit,
            'flux_unit': flux_unit} 


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
            output={'wave':wave,'flux':flux,'noise':noise,'contamination': contam,'wave_unit': u.Angstrom,'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom}
        except KeyError:
            try:
                output=readSpectrum(file)
                output['wave_unit']=u.Angstrom
                output['flux_unit']=u.erg/u.s/u.cm/u.cm/u.Angstrom
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
                'wave_unit': u.Angstrom,
                'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom}
   
    return output



#####################################################
###############   CLASSIFICATION   ##################
#####################################################


def classifyByIndex(sp,ref='burgasser',string_flag=True,round_flag=False,remeasure=True,nsamples=100,nloop=5,verbose=False,indices={},param={}, allmeasures=False,**kwargs):

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

    for k in ['set','reference','indexset']: ref = kwargs.get(k,ref)
    for k in ['str_flag','string']: string_flag = kwargs.get(k,string_flag)
    for k in ['round','rnd_flag','roundoff']: round_flag = kwargs.get(k,round_flag)
    for k in ['rem_flag']: remeasure = kwargs.get(k,remeasure)
    for k in ['parameters']: param = kwargs.get(k,param)
#    kwargs['ref'] = ref
#    if (set.lower() not in allowed_sets):
#        print('\nWarning: index classification method {} not present; returning nan\n\n'.format(set))
#        return numpy.nan, numpy.nan

# # Reid et al. (2001, AJ, 121, 1710)
#     if (ref.lower() == 'reid'):
#         if (rem_flag or len(args) == 0):
#             indices = measureIndexSet(sp, **kwargs)
#         sptoffset = 20.
#         coeffs = { \
#             'H2O-A': {'fitunc': 1.18, 'range': [18,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [-32.1, 23.4]}, \
#             'H2O-B': {'fitunc': 1.02, 'range': [18,28], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [-24.9, 20.7]}}
#         method='polynomial'

# # Testi et al. (2001, ApJ, 522, L147)
#     elif (ref.lower() == 'testi'):
#         if (rem_flag or len(args) == 0):
#             indices = measureIndexSet(sp, **kwargs)
#         sptoffset = 10.
#         coeffs = { \
#             'sHJ': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [-1.87, 1.67]}, \
#             'sKJ': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [-1.20, 2.01]}, \
#             'sH2O_J': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [1.54, 0.98]}, \
#             'sH2O_H1': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [1.27, 0.76]}, \
#             'sH2O_H2': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [2.11, 0.29]}, \
#             'sH2O_K': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [2.36, 0.60]}}
#         method='polynomial'

# # Burgasser (2007, ApJ, 659, 655) calibration
#     elif (ref.lower() == 'burgasser'):
#         if (rem_flag or len(args) == 0):
#             indices = measureIndexSet(sp, **kwargs)
#         sptoffset = 20.
#         coeffs = { \
#             'H2O-J': {'fitunc': 0.8, 'range': [20,39], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [1.038e2, -2.156e2,  1.312e2, -3.919e1, 1.949e1]}, \
#             'H2O-H': {'fitunc': 1.0, 'range': [20,39], 'spt': 0., 'sptunc': 99., 'mask': 1.,  \
#             'coeff': [9.087e-1, -3.221e1, 2.527e1, -1.978e1, 2.098e1]}, \
#             'CH4-J': {'fitunc': 0.7, 'range': [30,39], 'spt': 0., 'sptunc': 99., 'mask': 1.,  \
#             'coeff': [1.491e2, -3.381e2, 2.424e2, -8.450e1, 2.708e1]}, \
#             'CH4-H': {'fitunc': 0.3, 'range': [31,39], 'spt': 0., 'sptunc': 99., 'mask': 1.,  \
#             'coeff': [2.084e1, -5.068e1, 4.361e1, -2.291e1, 2.013e1]}, \
#             'CH4-K': {'fitunc': 1.1, 'range': [20,37], 'spt': 0., 'sptunc': 99., 'mask': 1.,  \
#             'coeff': [-1.259e1, -4.734e0, 2.534e1, -2.246e1, 1.885e1]}}
#         method='polynomial'

# # Geballe et al. (2002, ApJ, 564, 466) calibration
#     elif ref.lower() == 'geballe':
#         if (rem_flag or len(args) == 0):
#             i1 = measureIndexSet(sp, ref='geballe')
#             i2 = measureIndexSet(sp, ref='martin')
#             if sys.version_info.major == 2:
#                 indices = dict(i1.items() + i2.items())
#             else:
#                 indices = dict(i1.items() | i2.items())
#         spttypes = numpy.arange(20.,39.,1.)
#         ranges = { \
#             'PC3': [[2.4,2.6,20.],[2.6,2.86,21.],[2.85,3.25,22.],[3.25,4.25,23.],[4.25,6,24.]],\
#             'Color-d2': [[4.5,5.5,20.],[5.5,6.5,21.],[6.5,7.5,22.],[7.5,10.,23.],[10,17,24.],[17.,23.,25.],[23.,25.,26.]],\
#             'H2O-1.2': [[1.5,1.7,30.],[1.7,1.9,31.],[1.9,2.15,32.],[2.15,2.5,33.],[2.5,3.0,34.],[3.0,4.5,35.],[4.5,6.5,36.],[6.5,10.,37.],[10.,15.,38.]],\
#             'H2O-1.5': [[1.2,1.27,20.],[1.27,1.35,21.],[1.35,1.43,22.],[1.43,1.5,23.],[1.5,1.55,24.],[1.55,1.6,25.],[1.6,1.65,26.],[1.65,1.7,27.],[1.7,1.8,28.],[1.8,1.95,29.],[1.95,2.2,30.],[2.2,2.5,31.],[2.5,3.0,32.],[3.0,3.5,33.],[3.5,4.5,34.],[4.5,5.5,35.],[5.5,7.,36.],[7.,9.,37.],[9.,12.,38.]],\
#             'CH4-1.6': [[1.02,1.07,30.],[1.07,1.15,31.],[1.15,1.3,32.],[1.3,1.5,33.],[1.5,1.8,34.],[1.8,2.5,35.],[2.5,4,36.],[4.,6.,37.],[6.,9.,38.]],\
#             'CH4-2.2': [[0.91,0.94,23.],[0.94,0.98,24.],[0.98,1.025,25.],[1.025,1.075,26.],[1.075,1.125,27.],[1.125,1.175,28.],[1.175,1.25,29.],[1.25,1.4,30.],[1.4,1.6,31.],[1.6,1.95,32.],[1.95,2.75,33.],[2.75,3.8,34.],[3.8,5.5,35.],[5.5,8.5,36.],[8.5,12.,37],[12.,18.,38.]],\
#         }
#         method='ranges'

# # Allers et al. (2013, ApJ, 657, 511)
#     elif (ref.lower() == 'allers'):
#         if (rem_flag or len(args) == 0):
#             i1 = measureIndexSet(sp, ref='mclean')
#             i2 = measureIndexSet(sp, ref='slesnick')
#             i3 = measureIndexSet(sp, ref='allers')
#             if sys.version_info.major == 2:
#                 indices = dict(i1.items() + i2.items() + i3.items())
#             else:
#                 indices = dict(i1.items() | i2.items() | i3.items())
#         sptoffset = 10.
#         coeffs = { \
#             'H2O': {'fitunc': 0.390, 'range': [15,25], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [24.0476, -104.424, 169.388,-83.5437]}, \
#             'H2O-1': {'fitunc': 1.097, 'range': [14,25], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [28.5982, -80.7404, 39.3513, 12.1927]}, \
#             'H2OD': {'fitunc': 0.757, 'range': [20,28], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [-97.230, 229.884, -202.245, 79.4477]}, \
#             'H2O-2': {'fitunc': 0.501, 'range': [14,22], 'spt': 0., 'sptunc': 99., 'mask': 1., \
#             'coeff': [37.5013, -97.8144, 55.4580, 10.8822]}}
#         method='polynomial'

#     else:
#         if len(indices.keys()) == 0: 
#             print('Error: indices is an empty dictionary')
#             return numpy.nan, numpy.nan
#         if len(coeff.keys()) == 0: 
#             print('Error: coeffs is an empty dictionary')
#             return numpy.nan, numpy.nan

# check if using a pre-defined classification relation
# use that if classification parameters not provided
    tmp = checkDict(ref,INDEX_CLASSIFICATION_RELATIONS)
    if tmp==False and len(param) == 0: 
        print('Error: input citation {} is not one of the predefined classification sets: {}'.format(ref,list(INDEX_CLASSIFICATION_RELATIONS.keys())))
        return numpy.nan, numpy.nan
    elif len(param) == 0: 
        param = INDEX_CLASSIFICATION_RELATIONS[tmp]
        if verbose: print('Using {} index classification method'.format(tmp))
    else: pass

# check parameter dictionary has correct keys
    for k in ['indices','method','sets']:
        if k not in list(param.keys()): raise ValueError('Index classification parameter set must include the parameter {}; contains {}'.format(k,list(param.keys())))
    if 'decimal' not in list(param.keys()): param['decimal'] = False
    if 'sptoffset' not in list(param.keys()): param['sptoffset'] = 0
    if 'min_indices' not in list(param.keys()): param['min_indices'] = 1

# check coefficient index names are present in indices input; if not remeasure
    for i in list(param['indices'].keys()):
        if i not in list(indices.keys()):
            if remeasure==False: print('Error: coefficient index {} is not one of the measured indices provided in input indices; remeasuring'.format(i))
            remeasure = True

# remeasure indices            
    if remeasure==True or False in [k in indices.keys() for k in list(param['indices'].keys())]:
        for s in param['sets']:
            ind0 = copy.deepcopy(indices)
            tmp = checkDict(s,INDEX_SETS)
            if tmp==False: raise ValueError('Index set {} is not one of predefined sets: {}; correct parameter dictionary'.format(s,list(INDEX_SETS.keys())))
            ind = measureIndexSet(sp,ref=s,verbose=verbose,**kwargs)
            if sys.version_info.major == 2: indices = ind0.update(ind)
            else: indices = {**ind0,**ind}
# determine classifications from polynomials
    if param['method']=='polynomial':
        coeffs = {}
        for index in list(param['indices'].keys()):
#            if index in ['jtype', 'ktype']: continue
            coeffs[index] = {}
            coeffs[index]['spt'] = numpy.nan
            coeffs[index]['sptunc'] = numpy.nan
            coeffs[index]['mask'] = 0.
# check dictionary for each index has the right keys; if not, skip it
            skip = False
            for k in ['coeff','fitunc','fitunc','range']:
                if k not in list(param['indices'][index].keys()): 
                    if verbose==True: print('Parameter dictionary is missing key {}'.format(k))
                    skip = True
            if skip == False and numpy.isfinite(indices[index][0]):
                coeffs[index]['spt'] = numpy.polyval(param['indices'][index]['coeff'],indices[index][0])+param['sptoffset']
                coeffs[index]['sptunc'] = param['indices'][index]['fitunc']
                if param['decimal']==True: coeffs[index]['spt'] = (coeffs[index]['spt']-param['sptoffset'])*10.+param['sptoffset']
# noise sim
                if numpy.isfinite(indices[index][1]):
                    vals = numpy.polyval(param['indices'][index]['coeff'],numpy.random.normal(indices[index][0],indices[index][1],nsamples))
                    if param['decimal']==True: vals = (vals-param['sptoffset'])*10.+param['sptoffset']
                    coeffs[index]['sptunc'] = (numpy.nanstd(vals)**2+param['indices'][index]['fitunc']**2)**0.5

#                if verbose==True: print('{}: index = {:.3f}+/-{:.3f}, spt = {:.2f}+/-{:.2f}'.format(index,indices[index][0],indices[index][1],coeffs[index]['spt'],coeffs[index]['sptunc']))            

# unmask good values
            if coeffs[index]['spt'] >= numpy.nanmin(param['indices'][index]['range']) and \
                coeffs[index]['spt'] <= numpy.nanmax(param['indices'][index]['range']) and \
                numpy.isfinite(coeffs[index]['spt']): coeffs[index]['mask'] = 1.

# not enough good values
#         mask = [coeffs[index]['mask'] for index in list(coeffs.keys())]
# #        print(numpy.nansum(mask), param['min_indices'])
#         if numpy.nansum(mask) < param['min_indices']:
#             if verbose==True: print('\nNot of enough indices in set {} returned viable values\n'.format(ref))
#             return numpy.nan, numpy.nan

# if Allers method, need to compute the J and K spectral types
# TEMPORARILY REMOVING THIS AT IT BREAKS THE CODE
        #print('METHOD', ref.lower())
        # if ref.lower() in ['allers','all13','allers13','allers2013']:
        #     coeffs['jtype'] = {}
        #     coeffs['ktype'] = {}
        #     # Calculate The J and K spectral types for Allers method
        #     jtype = classifyByStandard(sp, std_class='dwarf', fit_ranges=[1.07,1.40]) 
        #     #print('JType', jtype)
        #     #print(typeToNum(jtype[0]))
        #     coeffs['jtype']['spt']    = typeToNum(jtype[0])
        #     coeffs['jtype']['sptunc'] = 1.
        #     coeffs['jtype']['mask']   = 1.
        #     ktype = classifyByStandard(sp, std_class='dwarf', fit_ranges=[1.90,2.20]) 
        #     #print('KType', ktype)
        #     #print(typeToNum(ktype[0]))
        #     coeffs['ktype']['spt']    = typeToNum(ktype[0])
        #     coeffs['ktype']['sptunc'] = 1.
        #     coeffs['ktype']['mask']   = 1.

# computed weighted mean with rejection, iterating to deal with indices outside ranges   
        #print(coeffs)   
        #sys.exit()
        for i in numpy.arange(nloop):
            wts = numpy.array([coeffs[index]['mask']/coeffs[index]['sptunc']**2 for index in list(coeffs.keys())])
            vals = numpy.array([coeffs[index]['mask']*coeffs[index]['spt']/coeffs[index]['sptunc']**2 for index in list(coeffs.keys())])
            w = numpy.where(numpy.isfinite(wts+vals))
            if len(w) == 0:
                if verbose==True: print('\nNone of the indices in set {} returned viable values\n'.format(ref))
                return numpy.nan, numpy.nan
            sptn = numpy.nansum(vals[w])/numpy.nansum(wts[w])
            sptn_e = 1./numpy.nansum(wts[w])**0.5
            for index in coeffs.keys():
                if sptn < numpy.nanmin(param['indices'][index]['range']) or sptn > numpy.nanmax(param['indices'][index]['range']): coeffs[index]['mask'] = 0.

# report individual subtypes
        if verbose == True:
            for i in coeffs.keys():
                flg = ''
                if coeffs[i]['mask'] == 0.: flg = '*'
#                if i in ['jtype', 'ktype']: continue
                print('{}{} = {:.3f}+/-{:.3f} = SpT = {}+/-{:.1f}'.format(flg,i,indices[i][0],indices[i][1],typeToNum(coeffs[i]['spt']),coeffs[i]['sptunc']))

# not enough good values
        if verbose == True: print('* = not included in classification')
        mask = [coeffs[index]['mask'] for index in list(coeffs.keys())]
        if numpy.nansum(mask) < param['min_indices']:
            if verbose==True: print('\nNot of enough indices in set {} returned viable values\n'.format(ref))
            return numpy.nan, numpy.nan

# ranges method - NOT CURRENTLY CONSIDERING UNCERTAINTY
    elif param['method']=='ranges':
        spts = []
        spts_unc = []
        for index in list(param['indices'].keys()):
            spts.append(numpy.nan)
            spts_unc.append(numpy.nan)
            if indices[index][1] > 0.:
                for r in param['indices'][index]: 
                    if r[0] < indices[index][0] <= r[1]: spts[-1] = r[-1]
#
#
# PLACEHOLDER FOR UNCERTAINTY INCLUSION
#
#
            if verbose == True:
                flg = ''
                if numpy.isnan(spts[-1]): flg='*'
                print('{}{}: {:.3f}+/-{:.3f} => SpT = {}'.format(flg,index,indices[index][0],indices[index][1],typeToNum(spts[-1])))

        if verbose == True: print('* = not included in classification')
        spts = numpy.array(spts)
        sptn = numpy.nanmean(spts)
        sptn_e = numpy.nanstd(spts)

    else: raise ValueError('Do not recognize index classification method {}'.format(param['method']))

# round off to nearest 0.5 subtypes if desired
    if round_flag==True: sptn = 0.5*numpy.around(sptn*2.)

# change to string if desired
    if string_flag==True: spt = typeToNum(sptn,uncertainty=sptn_e)
    else: spt = sptn

    if 'allmeasures' == True and param['method'] == 'polynomial':
        output = {}
        for k in coeffs.keys():
            output[k] = {'spt': coeffs[k]['spt'], 'spt_e': coeffs[k]['sptunc'], 'index': indices[k][0], 'index_e': indices[k][1]}
        output['result'] = (spt,sptn_e)
        return output
    else:
#        return spt, sptn_e, (jtype, ktype) # NOT SURE WHERE THIS CAME IN?
        return spt, sptn_e



def classifyByStandard(sp, std_class='dwarf',dof=-1, verbose=False,**kwargs):
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

#    verbose = kwargs.get('verbose',False)
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
    if dof<0: dof=len(sp.wave)


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
    allowed_classes = ['dwarf','sd','esd','dsd','vlg','intg','subdwarf','young','all']
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
    elif std_class == 'sd':
        initiateStandards(sd=True)
        stds = STDS_SD_SPEX
        subclass = 'sd'
        stdtype = 'subdwarf'
        if verbose==True: print('Using subdwarf standards')
    elif std_class == 'dsd':
        initiateStandards(dsd=True)
        stds = STDS_DSD_SPEX
        subclass = ''
        stdtype = 'mild subdwarf'
        if verbose == True: print('Using dwarf standards')
    elif std_class == 'esd':
        initiateStandards(esd=True)
        stds = STDS_ESD_SPEX
        subclass = 'esd'
        stdtype = 'extreme subdwarf'
        if verbose == True: print('Using extreme subdwarf standards')
    elif std_class == 'vlg':
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
    elif std_class == 'subdwarf':
        initiateStandards(sd=True)
        initiateStandards(dsd=True)
        initiateStandards(esd=True)
        stds = STDS_SD_SPEX.copy()
        stds.update(STDS_DSD_SPEX)
        stds.update(STDS_ESD_SPEX)
        stdtype = 'subdwarf'
        if verbose == True: print('Using d/sd, sd, and esd standards')
    elif std_class == 'young':
        initiateStandards(vlg=True)
        initiateStandards(intg=True)
        stds = STDS_VLG_SPEX.copy()
        stds.update(STDS_INTG_SPEX)
        stdtype = 'young'
        if verbose == True: print('Using INTG and VLG standards')
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
    elif (kwargs.get('method','').lower() == 'cruz'):
        fit_ranges = [[[0.87,1.39]],[[1.41,1.89]],[[1.91,2.39]]]    # as prescribed in Cruz et al. 2018, AJ,
    elif (kwargs.get('method','').lower() == ''):
        fit_ranges = [[0.7,2.45]]       # by default, compare whole spectrum
    else:
        print('\nWarning: do not recognize method = {}'.format(kwargs['method']))
        fit_ranges = [[0.7,2.45]]       # by default, compare whole spectrum
    fit_ranges = kwargs.get('fit_ranges',fit_ranges)
    fit_ranges = kwargs.get('fitrange',fit_ranges)
    fit_ranges = kwargs.get('fitrng',fit_ranges)
    fit_ranges = kwargs.get('comprange',fit_ranges)
    fit_ranges = kwargs.get('comprng',fit_ranges)
    if not isinstance(fit_ranges[0],list):
        fit_ranges = [fit_ranges]
    if verbose==True: print('Fitting over range {}'.format(fit_ranges))

# compute fitting statistics
    stat = []
    sspt = []

    for t in spt_sample:
#        chisq,scale = compareSpectra(sp,stds[typeToNum(t,subclass=subclass)],fit_ranges=fit_ranges,statistic=statistic,novar2=True)
        if (kwargs.get('method','').lower() == 'cruz'):
            chisq,scale = compareSpectraCruz(sp,stds[t],fit_ranges=fit_ranges,statistic=statistic,novar2=True,verbose=verbose)
        else: 
            chisq,scale = compareSpectra(sp,stds[t],fit_ranges=fit_ranges,statistic=statistic,novar2=True,verbose=verbose)
        stat.append(chisq)
        sspt.append(t)
        if verbose==True: print('Type {}: statistic = {}, scale = {}'.format(t, chisq, scale))

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
            mean,var = weightedMeanVar(ssptn,st,method='ftest',dof=dof)
        if (var**0.5 < 1.):
            sptn = numpy.round(mean*2)*0.5
        else:
            sptn = numpy.round(mean)
        sptn_e = (unc_sys**2+var**2)**0.5
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
        spstd = copy.deepcopy(stds[sorted_stdsptnum[0]])
#        getStandard(typeToNum(sorted_stdsptnum[0],subclass=subclass))
        if (kwargs.get('method','').lower() == 'cruz'):
            chisq,scale = compareSpectraCruz(sp,spstd,fit_ranges=fit_ranges,statistic=statistic,verbose=verbose)
            spstd.flux[numpy.where( (spstd.wave.value >= numpy.min(fit_ranges[0])) & ((spstd.wave.value <= numpy.max(fit_ranges[0])) ) )] *= scale[0]
            spstd.flux[numpy.where( (spstd.wave.value >= numpy.min(fit_ranges[1])) & ((spstd.wave.value <= numpy.max(fit_ranges[1])) ) )] *= scale[1]
            spstd.flux[numpy.where( (spstd.wave.value >= numpy.min(fit_ranges[2])) & ((spstd.wave.value <= numpy.max(fit_ranges[2])) ) )] *= scale[2]
            spstd.flux[numpy.where( ( (spstd.wave.value < numpy.min(fit_ranges[0])) | (spstd.wave.value > numpy.max(fit_ranges[0])) ) &
                                    ( (spstd.wave.value < numpy.min(fit_ranges[1])) | (spstd.wave.value > numpy.max(fit_ranges[1])) ) &
                                    ( (spstd.wave.value < numpy.min(fit_ranges[2])) | (spstd.wave.value > numpy.max(fit_ranges[2])) ) 
                                   )] = numpy.nan
        else:
            chisq,scale = compareSpectra(sp,spstd,fit_ranges=fit_ranges,statistic=statistic,verbose=verbose)
            spstd.scale(scale)
        if kwargs.get('colors',False) == False:
            kwargs['colors'] = ['k','r','b']
        if kwargs.get('labels',False) == False:
            kwargs['labels'] = [sp.name,'{} Standard'.format(sorted_stdsptnum[0]),'Difference']
        from .plot import plotSpectrum
        if kwargs.get('difference',True):
            kwargs['labels'].append('Difference')
            pl = plotSpectrum([sp,spstd,sp-spstd],**kwargs)
        else:
            pl = plotSpectrum([sp,spstd],**kwargs)

    if verbose==True: print('Fit-ranges:', fit_ranges)

    if kwargs.get('return_standard',False) == True: 
        spstd = copy.deepcopy(stds[sorted_stdsptnum[0]])
        if (kwargs.get('method','').lower() == 'cruz'):
            chisq,scale = compareSpectraCruz(sp,spstd,fit_ranges=fit_ranges,statistic=statistic,verbose=verbose)
        else:
            chisq,scale = compareSpectra(sp,spstd,fit_ranges=fit_ranges,statistic=statistic,verbose=verbose)
        spstd.scale(scale)
        return spstd
    elif kwargs.get('return_statistic',False) == True: 
        return output_spt, numpy.nanmin(stat)
    else:
        return output_spt, sptn_e




def classifyByTemplate(sp, templates, output='best', nbest=1, maxtemplates=100, plot=False, verbose=ERROR_CHECKING, dof=-1, **kwargs):
    '''
    Purpose
    -------

    Compares a Spectrum object to a dictionary of templates, identifiying the best fit template

    Parameters
    ----------

    sp : Spectrum class
        spectrum to copmare to templates

    templates : dict
        dictionary of spectrum objects, where each key corresponds to the template identifier (e.g., spectral type)
        and points to a Spectrum class object of the template

    output = 'best' : str [optional]
        signifies what should be returned by call to classifyByTemplate()

    nbest = 1 : int [optional]
        indicate what number of best-fit templates to return, in order of lowest fit statistic 

    maxtemplates = 100 : int [optional]
        maximum number of templates to compare to if time is an issue; -1 is unlimited 

    plot = False : bool [optional]
        plot comparison between source and best-fit template 

    verbose = False : bool [optional]
        return verbose feedback during execution 


    Outputs
    -------

    Depends on value of output keyword:
        * output = 'best' : returns the single or array (nbest > 1) of identifiers of the best-fit template(s)
        * output = 'template' : returns the single or array (nbest > 1) of Spectrum class objects of the best-fit template(s)
        * output = 'stats' : returns the single or array (nbest > 1) of tuples containing the identifier, fit statistic, and scaling factor for best-fit template(s)
        * output = 'all' : returns a dictionary where each key is a template identifier, which point to dictionaries 
        that contain the template, statistic, and scale factor


    Example
    -------
        TBD

    Dependencies
    ------------
        `splat.compareSpectra()`_
        copy

    '''

# check inputs
    if isinstance(sp,Spectrum)==False: raise ValueError('You must input a Spectrum object for comparison source, you passed {}'.format(type(sp)))
    if isinstance(templates,dict)==False: raise ValueError('You must pass a dictionary of templates; you passed {}'.format(type(templates)))
    if isinstance(list(templates.items())[0][1],Spectrum)==False: raise ValueError('Your templates must be Spectrum objects; you passed {}'.format(type(list(templates.items())[0][1])))
    if maxtemplates>0 and len(list(templates.items()))>maxtemplates and verbose==True:
        print('Warning: maxtemplates is set to {} and your template dictionary has {} templates; change maxtemplates to fit full list'.format(maxtemplates,len(list(templates.items()))))
    if nbest<1: nbest==1
    nbest = int(nbest)
    if verbose==True: print('Comparing spectrum of {} to {} templates'.format(sp.name,len(list(templates.items()))))

# compare spectra
    stats,scales,ids = [],[],[]
    for k in list(templates.keys()):

        mkwargs = copy.deepcopy(kwargs)
        mkwargs['plot'] = False

        st,scl = compareSpectra(sp,templates[k],**mkwargs)
        stats.append(st)
        scales.append(scl)
        ids.append(k)
        if verbose == True: print('\tcomparison to {}: stat = {}, scale = {}'.format(k,st,scl))

# sort by fit statistic (assumes smallest values is best)
    sorted_ids = [x for (y,x) in sorted(zip(stats,ids))]
    sorted_scales = [x for (y,x) in sorted(zip(stats,scales))]
    stats.sort()
    if verbose==True: print('Best fit template = {}, statistic = {}, scale = {}'.format(sorted_ids[0],stats[0],sorted_scales[0]))

# plot spectrum compared to best-fit spectrum
    if plot==True:
        mkwargs = copy.deepcopy(kwargs)
        mkwargs['plot'] = True

        st,scl = compareSpectra(sp,templates[sorted_ids[0]],**mkwargs)

# return based on output
    if output=='best':
        if nbest==1: return sorted_ids[0]
        else: return sorted_ids[:nbest]
    elif output=='template':
        if nbest==1: 
            tmp = copy.deepcopy(templates[sorted_ids[0]])
            tmp.scale(sorted_scales[0])
            return tmp
        else: 
            out = []
            for i in range(nbest):
                tmp = copy.deepcopy(templates[sorted_ids[i]])
                tmp.scale(sorted_scales[i])
            return out
    elif output=='stats':
        if nbest==1: return (sorted_ids[0],stats[0],sorted_scales[0])
        else: 
            out = []
            for i in range(nbest): out.append((sorted_ids[i],stats[i],sorted_scales[i]))
            return out
    elif output=='all':
        out = {}
        for i,k in enumerate(sorted_ids):
            tmp = copy.deepcopy(templates[k])
            tmp.scale(sorted_scales[i])
            out[k] = {'spectrum': tmp, 'statistic': stats[i], 'scale': sorted_scales[i]}
        return out
    else:
        raise ValueError('Do not recognize output format {}'.format(output))




# # select either best match or an ftest-weighted average
#     if (best_flag or len(stat) == 1):
#         spt = sorted_stdsptnum[0]
#         sptn = typeToNum(spt)
#         sptn_e = unc_sys
#     else:
#         ssptn = [typeToNum(s) for s in sspt]
#         try:
#             st = stat.value
#         except:
#             st = stat
#         if numpy.isnan(numpy.median(sp.noise)):
#             mean,var = weightedMeanVar(ssptn,st)
#         else:
#             mean,var = weightedMeanVar(ssptn,st,method='ftest',dof=dof)
#         if (var**0.5 < 1.):
#             sptn = numpy.round(mean*2)*0.5
#         else:
#             sptn = numpy.round(mean)
#         sptn_e = (unc_sys**2+var**2)**0.5
#         spt = typeToNum(sptn,uncertainty=sptn_e,subclass=subclass)

#     spt_type = kwargs.get('spt_type','literature')
#     spt = kwargs.get('spt',[10.,39.9])
#     spt = kwargs.get('spt_range',spt)
#     nbest = kwargs.get('nbest',1)
#     verbose = kwargs.get('verbose',False)
#     published = kwargs.get('published','')
#     published = kwargs.get('public',published)
#     statistic = kwargs.get('statistic','chisqr')
#     statistic = kwargs.get('stat',statistic)
#     force = kwargs.get('force',False)
#     maxtemplates = kwargs.get('maxtemplates',100)
#     select = kwargs.get('select','')
#     select = kwargs.get('set',select)
# #   placeholder for a systematic uncertainty term
#     unc_sys = 0.
#     if (kwargs.get('method','') == 'kirkpatrick'):
#         fit_ranges = [[0.9,1.4]]         # as prescribed in Kirkpatrick et al. 2010, ApJS,
#     else:
#         fit_ranges = [[0.7,2.45]]       # by default, compare whole spectrum
#     fit_ranges = kwargs.get('fitrange',fit_ranges)
#     fit_ranges = kwargs.get('fitrng',fit_ranges)
#     fit_ranges = kwargs.get('comprange',fit_ranges)
#     fit_ranges = kwargs.get('comprng',fit_ranges)
#     if not isinstance(fit_ranges[0],list):
#         fit_ranges = [fit_ranges]
#     if dof<0: dof = len(sp.wave)

# #  canned searches
# #  constrain spectral types
#     if ('lit' in spt_type.lower()):
#         spt_type = 'LIT_TYPE'
#     elif ('opt' in spt_type.lower() or 'optical' in select):
#         spt_type = 'OPT_TYPE'
#     elif ('nir' in spt_type.lower()):
#         spt_type = 'NIR_TYPE'
#     else:
#         spt_type = 'LIT_TYPE'

#     if ('m dwarf' in select.lower() or kwargs.get('mdwarf',False)):
#         spt = [numpy.max([10,spt[0]]),numpy.min([19.9,spt[-1]])]
#     if ('l dwarf' in select.lower() or kwargs.get('ldwarf',False)):
#         spt = [numpy.max([20,spt[0]]),numpy.min([29.9,spt[-1]])]
#     if ('t dwarf' in select.lower() or kwargs.get('tdwarf',False)):
#         spt = [numpy.max([30,spt[0]]),numpy.min([39.9,spt[-1]])]
#     if ('vlm' in select.lower() or kwargs.get('vlm',False)):
#         spt = [numpy.max([17,spt_range[0]]),numpy.min([39.9,spt_range[-1]])]

# #  constrain S/N
#     snr = 0.
#     if ('high sn' in select.lower()):
#         snr = 100.
#     snr = kwargs.get('snr',snr)

# #  don't compare to same spectrum
#     try:
#         excludefile = [sp.filename]
#     except:
#         excludefile = []
#     if kwargs.get('excludefile',False) != False:
#         e = kwargs.get('excludefile')
#         if isinstance(e,list):
#             excludefile.extend(e)
#         else:
#             excludefile.append(e)
#     try:
#         excludekey = [sp.data_key]
#     except:
#         excludekey = []
#     if kwargs.get('excludekey',False) != False:
#         e = kwargs.get('excludekey')
#         if isinstance(e,list):
#             excludekey.extend(e)
#         else:
#             excludekey.append(e)
#     try:
#         excludeshortname = [sp.shortname]
#     except:
#         excludeshortname = []
#     if kwargs.get('excludeshortname',False) != False:
#         e = kwargs.get('excludeshortname')
#         if isinstance(e,list):
#             excludeshortname.extend(e)
#         else:
#             excludeshortname.append(e)
# #    print(excludefile, excludekey, excludeshortname)

# # other classes
#     giant = ''
#     if 'giant' in select.lower() or kwargs.get('giant',False):
#         giant = True
#     if 'not giant' in select.lower():
#         giant = False
#     companion = ''
#     if 'companion' in select.lower() or kwargs.get('companion',False):
#         companion = True
#     if 'not companion' in select.lower():
#         companion = False
#     young = ''
#     if 'young' in select.lower() or kwargs.get('young',False):
#         young = True
#     if 'not young' in select.lower():
#         young = False
#     binary = ''
#     if 'binary' in select.lower() or kwargs.get('binary',False):
#         binary = True
#     if 'not binary' in select.lower():
#         binary = False
#     spbinary = ''
#     if 'spectral binary' in select.lower() or kwargs.get('sbinary',False):
#         spbinary = True
#     if 'not spectral binary' in select.lower():
#         spbinary = False

# # REARRANGE THIS - SEND IN KWARGS WITH OUTPUT, LOGIC SET, AND THE REST ARE UP TO USER?

#     lib = searchLibrary(excludefile=excludefile,excludekey=excludekey,excludeshortname=excludeshortname, \
#         snr=snr,spt_type=spt_type,spt=spt,published=published, \
#         giant=giant,companion=companion,young=young,binary=binary,spbinary=spbinary,output='all',logic='and')

# # first search for the spectra desired - parameters are set by user
#     if len(lib) == 0:
#         print('\nNo templates available for comparison\n\n')
#         return numpy.nan, numpy.nan

#     if len(lib) > maxtemplates and force == False:
#         print('\nToo many templates ({}) for classifyByTemplate; set force=True to override this\n\n'.format(len(lib)))
#         return numpy.nan, numpy.nan

#     files = lib['DATA_FILE']
#     dkey = lib['DATA_KEY']
#     sspt = [typeToNum(s) for s in lib[spt_type]]

#     if verbose == True:
#         print('\nComparing to {} templates\n'.format(len(files)))
#         if len(files) > 100:
#             print('This may take some time!\n\n'.format(len(files)))

# # do comparison
#     stat = []
#     scl = []
#     for i,d in enumerate(dkey):

# # INSERT TRY STATEMNT HERE?

#         s = Spectrum(d)

#         mkwargs = copy.deepcopy(kwargs)
#         if 'plot' in (mkwargs.keys()): mkwargs['plot'] = False

#         stt,scale = compareSpectra(sp,s,fit_ranges=fit_ranges,statistic=statistic,novar2=True,**mkwargs)
#         stat.append(stt)
#         scl.append(scale)
#         if verbose == True: print(s)

# # list of sorted standard files and spectral types
#     sorted_dkey = [x for (y,x) in sorted(zip(stat,dkey))]
#     sorted_spt = [x for (y,x) in sorted(zip(stat,sspt))]
#     sorted_scale = [x for (y,x) in sorted(zip(stat,scl))]

# # select either best match or an ftest-weighted average
#     if (kwargs.get('best',False) or len(stat) == 1):
#         sptn = sorted_spt[0]
#         sptn_e = unc_sys
#     else:
#         mean,var = weightedMeanVar(sspt,stat,method='ftest',dof=dof)
# # allow 1/2 subtypes if uncertainty is less than 1.0
#         if (var**0.5 < 1.):
#             sptn = numpy.round(mean*2.)*0.5
#         else:
#             sptn = numpy.round(mean)
#         sptn_e = (unc_sys**2+var)**0.5

# # plot spectrum compared to best spectrum
#     if (kwargs.get('plot',False) != False):
#         s = Spectrum(sorted_dkey[0])
# #        chisq,scale = compareSpectra(s,sp,fit_ranges=[comprng],stat='chisqr',novar2=True)
#         s.scale(sorted_scale[0])
#         kwargs['legend'] = [sp.name,s.name]
#         kwargs['colors'] = ['k','r','b']
#         from .plot import plotSpectrum
#         plotSpectrum(sp,s,sp-s,**kwargs)

# # string or not?
#     if (kwargs.get('string', True) == True):
#         output_spt = typeToNum(sptn,uncertainty=sptn_e)
#     else:
#         output_spt = sptn

#     if verbose == True:
#         s = Spectrum(sorted_dkey[0])
#         print('\nBest match = {} with spectral type {}'.format(s,typeToNum(sorted_spt[0])))
#         print('Mean spectral type = {}+/-{}'.format(output_spt,sptn_e))

# # return dictionary of results
#     return {'result': (output_spt,sptn_e), \
#         'statistic': sorted(stat)[0:nbest], 'spt': sorted_spt[0:nbest], \
#         'scale': sorted_scale[0:nbest], \
#         'spectra': [Spectrum(d) for d in sorted_dkey[0:nbest]]}



def classifyGravity(sp, output='classification',verbose=ERROR_CHECKING, **kwargs):
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

# NB: WHOEVER EDITED THIS BROKE THE ROUTINE - PLEASE TEST BEFORE UPDATING
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
        'H-cont': {'M5.0': [numpy.nan,numpy.nan], 'M6.0': [.988,.994], 'M7.0': [.981,.990],'M8.0': [.963,.984],'M9.0': [.949,.979],'L0.0': [.935,.972],'L1.0': [.914,.968],'L2.0': [.906,.964],'L3.0': [.898,.960],'L4.0': [.885,.954],'L5.0': [.869,.949],'L6.0': [.874,.950],'L7.0': [0.888,0.952]},\
        }

    if verbose==True: print('\nGravity Classification')
    gravscore = {}
    gravscore['gravity_class'] = 'UNKNOWN'

# Calculate Allers indices and their uncertainties
    ind = kwargs.get('indices',False)
    if ind == False:
        ind = measureIndexSet(sp,set='allers')
    gravscore['indices'] = ind

# Determine the object's NIR spectral type and its uncertainty
    sptn = kwargs.get('spt',False)
    if sptn == False:
#        sptn, spt_e, spts2 = classifyByIndex(sp,string=False,ref='allers2013')
        sptn, spt_e = classifyByIndex(sp,string=False,ref='allers2013')
        if numpy.isnan(sptn):
            if verbose==True: print('Spectral type could not be determined from indices; try entering with spt keyword')
            if output=='allmeasures': return gravscore
            else: return gravscore['gravity_class']

    if isinstance(sptn,str): sptn = typeToNum(sptn)
    Spt = typeToNum(numpy.round(sptn))
    gravscore['spt'] = Spt
    # if verbose==True: 
    #     print('\tSpT = {} (J-type = {}; K-type = {})'.format(Spt, spts2[0][0], spts2[1][0]))

#Check whether the NIR SpT is within gravity sensitive range values
    if ((sptn < 16.0) or (sptn > 27.0)):
        if verbose==True: print('Spectral type '+typeToNum(sptn)+' outside range for gravity classification')
        if output=='allmeasures': return gravscore
        else: return gravscore['gravity_class']

#Creates an empty array with dimensions 4x1 to fill in later with 5 gravscore values
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
            if verbose==True: print('\t{}: {:.3f}+/-{:.3f} => {}'.format(k,ind[k][0], ind[k][1], val))
        if k == 'FeH-z' or k=='KI-J':
            if numpy.isnan(grav[k][Spt][0]):
                val = numpy.nan
            if ind[k][0] <= grav[k][Spt][0]:
                val = 1.0
            if ind[k][0] <= grav[k][Spt][1]:
                val = 2.0
            if verbose==True: print('\t{}: {:.3f}+/-{:.3f} => {}'.format(k,ind[k][0], ind[k][1], val))
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
    if verbose==True: print('\tGravity Class = {}\n'.format(gravscore['gravity_class']))


# plot spectrum against standard
    if (kwargs.get('plot',False) != False):
        spt,unc = classifyByStandard(sp,compareto=Spt,method='kirkpatrick',**kwargs)

# return gravity class or entire dictionary
    if output=='allmeasures': return gravscore
    else: return gravscore['gravity_class']



def classifySB(sp,ref='burgasser2010',output='classification',spt='',indices=None,classify='standard',classify_parameters=None,index_definitions=SPECTRAL_BINARY_INDICES,verbose=False,**kwargs):
    """
    Determines whether a spectrum of an ultracool dwarf is that of a spectral binary using the index-based methods
    of Burgasser et al. (2010; 2010ApJ...710.1142B) and Bardalez Gagliuffi et al. (2014; 2014ApJ...794..143B)


    Parameters
    ----------
    sp : Spectrum class object


    ref : string, default = 'burgasser2010'
        reference for the index/calibration set to use for determining spectral binary status, one of the following:
            * 'burgasser2010': indices from Burgasser et al. (2010; 2010ApJ...710.1142B) appropriate for spectral types L7--T5
            * 'bardalez2014': indices from Bardalez Gagliuffi et al. (2014; 2014ApJ...794..143B) appropriate for spectral types M7-L8

    output : string, default = 'classification'
        what to return from function call, one of the following
            * 'classification': returns a string indicating the spectral binary classification, either 'strong','weak', or 'non-binary'
            * 'criteria': return an dictionary with the index measures as key point to whether the criteria were satisfied, and an overall classification
            * 'allmeasures': return an dictionary containing index measurements, criteria outcomes, and overall classification

    spt : string, default = ''
        string or number indicating input spectral type; if not provided, spectrum is classified

    indices : dict, default = {}
        dictionary containing index measurements, with keys corresponding to index names required for spectral binary classification;
        if empty or partially complete, indices will be remeasured

    classify : string, default = 'indices'
        if spectrum is to be reclassified, determines method for classification; options are
            * 'standard': use classifyByStandard()
            * 'template': use classifyByTemplate(); classify_parameters should include template array
            * 'indices': use classifyByIndex(); classify_parameters should include classification reference

    classify_parameters : dict, default = {}
        if spectrum is to be reclassified, optional parameters to pass onto classification functions
    
    index_definitions : dict, default = SPECTRAL_BINARY_INDICES
        dictionary containing the relevant spectral binary index parameters; this should be kept default to avoid code issues
    
    verbose : bool, default = False
        If True, return verbose output

    Returns
    -------
    string or dict
        Output depends on wht is set by 'output' parameter:
            * 'classification': returns a string indicating the spectral binary classification, either 'strong','weak', or 'non-binary'
            * 'criteria': return an dictionary with the index measures as key point to whether the criteria were satisfied, and an overall classification
            * 'allmeasures': return an dictionary containing index measurements, criteria outcomes, and overall classification

    Examples
    --------
    TBD

    See Also
    --------
    classifyByIndex : classify using a pre-defined set of indices
    classifyByStandard : classify using a pre-defined set of spectral standards
    classifyByTemplate : classify using an input set of spectra
    classifyGravity : determine the gravity class of a spectrum

    """

# keyword replacement
    for k in ['set','reference']: ref = kwargs.get(k,ref)
    for k in ['out','return']: output = kwargs.get(k,output)
    for k in ['type','spectral_type','spectral type']: spt = kwargs.get(k,spt)
    for k in ['index','measures']: indices = kwargs.get(k,indices)
    for k in ['classify_method']: classify = kwargs.get(k,classify)

# input checking
    if isinstance(sp,Spectrum)==False: 
        raise ValueError('Input variable should be a single Spectrum object; you passed {}'.format(type(sp)))
    tmp = checkDict(ref,index_definitions)
    if tmp==False: raise ValueError('Reference set {} is not one of the options for the spectral binary indices: {}'.format(ref,list(index_definitions,keys())))
    index_data = copy.deepcopy(index_definitions[tmp])
    if indices==None: indices = {}
    if classify_parameters==None: classify_parameters = {}

# is spectrum classified? if not, conduct classification
    if isinstance(spt,str)==True and index_data['spt']==True:
        if spt=='' and 'SPT' not in list(indices.keys()):
            if classify.lower() in ['template','templates']: 
                if len(classify_parameters.keys())==0: 
                    if verbose==True: print('classifyByTemplate requires a template list to be passed in classify_parmaeters dictionary, which is empty; running classifyByStandard instead')
                    classify = 'standard'
                else: 
                    spt,spt_e = classifyByTemplate(sp,**classify_parameters)
                    if verbose==True: print('\nClassified spectrum by template: result = {}'.format(spt))
            if classify.lower() in ['standard','std','standards']: 
                if len(classify_parameters.keys())==0: 
                    classify_parameters['method']='kirkpatrick'
                    classify_parameters['telluric']=True
                spt,spt_e = classifyByStandard(sp,**classify_parameters)
                if verbose==True: print('\nClassified spectrum by standards: result = {}'.format(spt))
            if classify.lower() in ['index','indices','ind']: 
                if len(classify_parameters.keys())==0: classify_parameters['ref']='burgasser2007'
                spt,spt_e = classifyByIndex(sp,**classify_parameters)
                if verbose==True: print('\nClassified spectrum by indices: result = {}'.format(spt))
            if spt=='' : 
                if verbose==True: print('classify parameter {} is not one of the allowed values; classifying by indices'.format(classify))
#    if not isNumber(spt): spt = typeToNum(spt)

# check indices - first build up set of required indices
    remeasure=False
    if len(list(indices.keys()))==0: remeasure=True
    else:
        index_required = []
        for k in index_data['relations']: index_required.extend(k['indices'])
        index_required = list(set(index_required))
        index_required.remove('SPT')
        for k in index_required:
            if k not in list(indices.keys()): remeasure=True
    if remeasure==True:
        if verbose==True: print('\nMeasuring indices using the {} definitions'.format(index_data['index_set']))
        indices_measured = measureIndexSet(sp,ref=index_data['index_set'])
        for k in list(indices_measured.keys()):
            if k not in list(indices.keys()): 
                indices[k] = indices_measured[k]
                if verbose==True and isinstance(indices[k],tuple)==True: 
#                    print(k,indices[k])
                    print('\t{} = {:.2f}+/-{:.2f}'.format(k,indices[k][0],indices[k][1]))

    if index_data['spt']==True:
        if 'SPT' not in list(indices.keys()): indices['SPT'] = spt
        if isinstance(indices['SPT'],list): indices['SPT']=indices['SPT'][0]
        if isinstance(indices['SPT'],str): indices['SPT']=typeToNum(indices['SPT'])
        indices['SPT']=indices['SPT']+index_data['sptoffset']

    indices_applied = copy.deepcopy(indices)
    for k in list(indices.keys()):
        if isinstance(indices_applied[k],list): indices_applied[k] = indices_applied[k][0]
        if isinstance(indices_applied[k],tuple): indices_applied[k] = indices_applied[k][0]

# now apply index methods
    criteria=[]
    class_dict = {}

    if verbose==True: print('\nUsing the {} spectral binary classification indices'.format(tmp))
    for k in index_data['relations']:
        crit = 1
# extract specific values of x, y
        if k['measure']=='ratio1':
            x = indices_applied[k['indices'][0]]/indices_applied[k['indices'][1]]
            key = '{}/{} vs {}'.format(k['indices'][0],k['indices'][1],k['indices'][2])
        else: x = indices_applied[k['indices'][0]]
        if k['measure']=='ratio2':
            y = indices_applied[k['indices'][-2]]/indices_applied[k['indices'][-1]]
            key = '{} vs {}/{}'.format(k['indices'][0],k['indices'][1],k['indices'][2])
        else: y = indices_applied[k['indices'][-1]]
        if k['measure']!='ratio1' and k['measure']!='ratio2':
            key = '{} vs {}'.format(k['indices'][0],k['indices'][1])
# apply test
        ycrit = numpy.polyval(k['coeff'],x)
        if k['direction']=='low' and y>ycrit: crit = 0
        if k['direction']=='high' and y<ycrit: crit = 0
# addition xlim, ylim constraints
        if 'xlim' in list(k.keys()):
            if x < k['xlim'][0] or x > k['xlim'][1]: crit = 0
        if 'ylim' in list(k.keys()):
            if y < k['ylim'][0] or y > k['ylim'][1]: crit = 0
        criteria.append(crit)
        class_dict[key] = True if crit == 1 else False
        if verbose==True: print('\t{} test: {}'.format(key,class_dict[key]))


# classify
    classification = 'fail'
    if numpy.sum(numpy.array(criteria))>=index_data['weak']: classification = 'weak'
    if numpy.sum(numpy.array(criteria))>=index_data['strong']: classification = 'strong'
    if indices['SPT'] < index_data['spt_range'][0] or indices['SPT'] > index_data['spt_range'][1]: 
        if verbose==True: print('\nSpectral classification {} is outside range of {} relation: {}'.format(typeToNum(indices['SPT']),ref,[typeToNum(x) for x in index_data['spt_range']]))
        classification='fail'
    if verbose==True: print('\nSpectral binary classification: {}'.format(classification))

# return output
    if output in ['class','classification']: return classification    
    if output in ['criteria','cases','sets','individual','detail']: return class_dict
    else: 
        outval = {}
        outval['indices'] = indices
        outval['criteria'] = class_dict
        outval['classification'] = classification
        return outval
    return



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
    if sp1.wave.unit != sp2.wave.unit: sp2.toWaveUnit(sp1.wave.unit)
    if sp1.flux.unit != sp2.flux.unit: sp2.toFluxUnit(sp1.flux.unit)

    fit_ranges = [[numpy.nanmin(sp1.wave),numpy.nanmax(sp1.wave)]]
    for k in ['fit_ranges','fit_range','fitrange','fitrng','comprange','comprng']:
        if k in list(kwargs.keys()): fit_ranges = kwargs[k]

#    mask_ranges = kwargs.get('mask_ranges',[])
#    mask_standard = kwargs.get('mask_standard',False)
#    mask_telluric = kwargs.get('mask_telluric',mask_standard)
    var_flag = novar2
    if numpy.isnan(numpy.max(sp2.variance.value)) == True: var_flag = True
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
    vtot[numpy.isnan(vtot)==True] = numpy.nanmedian(vtot)

# manage fit ranges and generate fit mask
    if len(fit_ranges) == 0: fit_ranges = [[numpy.nanmin(sp1.wave),numpy.nanmax(sp1.wave)]]
    if len(fit_ranges) == 2 and not isinstance(fit_ranges[0],list): fit_ranges = [fit_ranges]
    for i,m in enumerate(fit_ranges):
        if not isUnit(m):
            fit_ranges[i] = (m*u.micron).to(sp1.wave.unit)
    fit_mask = kwargs.get('fit_mask',1.-generateMask(sp1.wave,mask_ranges=fit_ranges))

# generate masking array and combine with fit mask
#    reject_mask = numpy.array(kwargs.get('mask',generateMask(sp1.wave,**kwargs)))
    reject_mask = numpy.array(kwargs.get('mask',numpy.zeros(len(sp1.wave))))
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
        vtot[numpy.isnan(vtot)==True] = numpy.nanmedian(vtot)
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave.value)*scale_factor)**2/vtot)
        unit = sp1.flux_unit/sp1.flux_unit

# normalized standard deviation
    elif (statistic == 'stddev_norm' or statistic == 'stdev_norm'):
# compute scale factor
        if scale == True:
            scale_factor = numpy.nansum(weights*sp1.flux.value)/ \
                numpy.nansum(weights*f(sp1.wave.value))
# correct variance
        vtot = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor**2)],axis=0)
        vtot[numpy.isnan(vtot)==True] = numpy.nanmedian(vtot)
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave.value)*scale_factor)**2)/ \
            numpy.nanmedian(sp1.flux.value)**2
        unit = sp1.flux_unit/sp1.flux_unit

# standard deviation
    elif (statistic == 'stddev' or statistic == 'stdev'):
# compute scale factor
        if scale == True:
            scale_factor = numpy.nansum(weights*sp1.flux.value*f(sp1.wave.value))/ \
                numpy.nansum(weights*f(sp1.wave.value)*f(sp1.wave.value))
# correct variance
        vtot = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor**2)],axis=0)
        vtot[numpy.isnan(vtot)==True] = numpy.nanmedian(vtot)
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave.value)*scale_factor)**2)
        unit = sp1.flux_unit**2

# absolute deviation
    elif (statistic == 'absdev'):
# compute scale factor
        if scale == True:
            scale_factor = numpy.nansum(weights*sp1.flux.value)/ \
                numpy.nansum(weights*f(sp1.wave.value))
# correct variance
        vtot = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor**2)],axis=0)
        vtot[numpy.isnan(vtot)==True] = numpy.nanmedian(vtot)
        stat = numpy.nansum(weights*abs(sp1.flux.value-f(sp1.wave.value)*scale_factor))
        unit = sp1.flux_unit

# error
    else:
        print('Error: statistic {} for compareSpectra not available'.format(statistic))
        return numpy.nan, numpy.nan

# plot spectrum compared to best spectrum
    if plot == True:
        spcomp = copy.deepcopy(sp2)
        spcomp.scale(scale_factor)
        kwargs['colors'] = kwargs.get('colors',['k','m','b'])
        kwargs['title'] = kwargs.get('title',sp1.name+' vs '+sp2.name)
        from .plot import plotSpectrum
        plotSpectrum([sp1,spcomp,sp1-spcomp],labels=[sp1.name,sp2.name,'{} = {}'.format(statistic,stat)],**kwargs)


    return numpy.nanmax([stat,minreturn])*unit, scale_factor


def compareSpectraCruz(s1, s2, statistic='chisqr',scale=True, novar2=True, plot=False, verbose=False, **kwargs):
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
    if sp1.wave.unit != sp2.wave.unit: sp2.toWaveUnit(sp1.wave.unit)
    if sp1.flux.unit != sp2.flux.unit: sp2.toFluxUnit(sp1.flux.unit)

    fit_ranges = [[[0.87,1.39]],[[1.41,1.89]],[[1.91,2.39]]]    # as prescribed in Cruz et al. 2018, AJ

#    mask_ranges = kwargs.get('mask_ranges',[])
#    mask_standard = kwargs.get('mask_standard',False)
#    mask_telluric = kwargs.get('mask_telluric',mask_standard)
    var_flag = novar2
    if numpy.isnan(numpy.max(sp2.variance.value)) == True: var_flag = True
    if numpy.isnan(numpy.max(sp1.variance.value)) == True: var_flag = False
    statistic = kwargs.get('stat',statistic)
    minreturn = 1.e-60
    scale_factor = 1.
    stats, scale_factors = [], []

# create interpolation function for second spectrum
    f = interp1d(sp2.wave.value,sp2.flux.value,bounds_error=False,fill_value=0.)
    if var_flag:
        v = interp1d(sp2.wave.value,[numpy.nan for s in sp2.wave],bounds_error=False,fill_value=numpy.nan)
    else:
        v = interp1d(sp2.wave.value,sp2.variance.value,bounds_error=False,fill_value=numpy.nan)
# total variance - funny form to cover for nans
    vtot = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)],axis=0)
    vtot[numpy.isnan(vtot)==True] = numpy.nanmedian(vtot)

# manage fit ranges and generate fit mask
    if not isUnit(fit_ranges[0]):
        fit_ranges[0] = (fit_ranges[0]*u.micron).to(sp1.wave.unit)
    if not isUnit(fit_ranges[1]):
        fit_ranges[1] = (fit_ranges[1]*u.micron).to(sp1.wave.unit)
    if not isUnit(fit_ranges[2]):
        fit_ranges[2] = (fit_ranges[2]*u.micron).to(sp1.wave.unit)
    fit_mask1 = kwargs.get('fit_mask',1.-generateMask(sp1.wave,mask_ranges=fit_ranges[0]))
    fit_mask2 = kwargs.get('fit_mask',1.-generateMask(sp1.wave,mask_ranges=fit_ranges[1]))
    fit_mask3 = kwargs.get('fit_mask',1.-generateMask(sp1.wave,mask_ranges=fit_ranges[2]))

# generate masking array and combine with fit mask
#    reject_mask = numpy.array(kwargs.get('mask',generateMask(sp1.wave,**kwargs)))
    reject_mask = numpy.array(kwargs.get('mask',numpy.zeros(len(sp1.wave))))
# mask flux < 0
    reject_mask[numpy.where(numpy.logical_or(sp1.flux < 0,f(sp1.wave) < 0))] = 1
    mask1 = numpy.clip(fit_mask1+reject_mask,0,1)
    mask2 = numpy.clip(fit_mask2+reject_mask,0,1)
    mask3 = numpy.clip(fit_mask3+reject_mask,0,1)

# set the weights
    weights1 = kwargs.get('weights',numpy.ones(len(sp1.wave)))
    weights2 = kwargs.get('weights',numpy.ones(len(sp1.wave)))
    weights3 = kwargs.get('weights',numpy.ones(len(sp1.wave)))
    weights1 = weights1*(1.-mask1)
    weights2 = weights2*(1.-mask2)
    weights3 = weights3*(1.-mask3)

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
            scale_factor1 = numpy.nansum(weights1*sp1.flux.value*f(sp1.wave.value)/vtot)/ \
                numpy.nansum(weights1*f(sp1.wave.value)*f(sp1.wave.value)/vtot)
            scale_factor2 = numpy.nansum(weights2*sp1.flux.value*f(sp1.wave.value)/vtot)/ \
                numpy.nansum(weights2*f(sp1.wave.value)*f(sp1.wave.value)/vtot)
            scale_factor3 = numpy.nansum(weights3*sp1.flux.value*f(sp1.wave.value)/vtot)/ \
                numpy.nansum(weights3*f(sp1.wave.value)*f(sp1.wave.value)/vtot)

# correct variance
        vtot1 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor1**2)],axis=0)
        vtot2 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor2**2)],axis=0)
        vtot3 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor3**2)],axis=0)
        vtot1[numpy.isnan(vtot1)==True] = numpy.nanmedian(vtot1)
        vtot2[numpy.isnan(vtot2)==True] = numpy.nanmedian(vtot2)
        vtot3[numpy.isnan(vtot3)==True] = numpy.nanmedian(vtot3)
        stat1 = numpy.nansum(weights1*(sp1.flux.value-f(sp1.wave.value)*scale_factor1)**2/vtot1)
        stat2 = numpy.nansum(weights2*(sp1.flux.value-f(sp1.wave.value)*scale_factor2)**2/vtot2)
        stat3 = numpy.nansum(weights3*(sp1.flux.value-f(sp1.wave.value)*scale_factor3)**2/vtot3)
        stat = stat1+stat2+stat3
        unit = sp1.flux_unit/sp1.flux_unit

# normalized standard deviation (NEED TO DO ALL OF THESE TOO!)
    elif (statistic == 'stddev_norm' or statistic == 'stdev_norm'):
# compute scale factor
        if scale == True:
            scale_factor1 = numpy.nansum(weights1*sp1.flux.value)/ \
                numpy.nansum(weights1*f(sp1.wave.value))
            scale_factor2 = numpy.nansum(weights2*sp1.flux.value)/ \
                numpy.nansum(weights2*f(sp1.wave.value))
            scale_factor3 = numpy.nansum(weights3*sp1.flux.value)/ \
                numpy.nansum(weights3*f(sp1.wave.value))
# correct variance
        vtot1 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor1**2)],axis=0)
        vtot2 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor2**2)],axis=0)
        vtot3 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor3**2)],axis=0)
        vtot1[numpy.isnan(vtot1)==True] = numpy.nanmedian(vtot1)
        vtot2[numpy.isnan(vtot2)==True] = numpy.nanmedian(vtot2)
        vtot3[numpy.isnan(vtot3)==True] = numpy.nanmedian(vtot3)
        stat1 = numpy.nansum(weights1*(sp1.flux.value-f(sp1.wave.value)*scale_factor1)**2)/ \
            numpy.nanmedian(sp1.flux.value)**2
        stat2 = numpy.nansum(weights2*(sp1.flux.value-f(sp1.wave.value)*scale_factor2)**2)/ \
            numpy.nanmedian(sp1.flux.value)**2
        stat3 = numpy.nansum(weights3*(sp1.flux.value-f(sp1.wave.value)*scale_factor3)**2)/ \
            numpy.nanmedian(sp1.flux.value)**2
        stat = stat1+stat2+stat3
        unit = sp1.flux_unit/sp1.flux_unit

# standard deviation
    elif (statistic == 'stddev' or statistic == 'stdev'):
# compute scale factor
        if scale == True:
            scale_factor1 = numpy.nansum(weights1*sp1.flux.value*f(sp1.wave.value))/ \
                numpy.nansum(weights1*f(sp1.wave.value)*f(sp1.wave.value))
            scale_factor2 = numpy.nansum(weights2*sp1.flux.value*f(sp1.wave.value))/ \
                numpy.nansum(weights2*f(sp1.wave.value)*f(sp1.wave.value))
            scale_factor3 = numpy.nansum(weights3*sp1.flux.value*f(sp1.wave.value))/ \
                numpy.nansum(weights3*f(sp1.wave.value)*f(sp1.wave.value))
# correct variance
        vtot1 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor1**2)],axis=0)
        vtot2 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor2**2)],axis=0)
        vtot3 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor3**2)],axis=0)
        vtot1[numpy.isnan(vtot1)==True] = numpy.nanmedian(vtot1)
        vtot2[numpy.isnan(vtot2)==True] = numpy.nanmedian(vtot2)
        vtot3[numpy.isnan(vtot3)==True] = numpy.nanmedian(vtot3)
        stat1 = numpy.nansum(weights1*(sp1.flux.value-f(sp1.wave.value)*scale_factor1)**2)
        stat2 = numpy.nansum(weights2*(sp1.flux.value-f(sp1.wave.value)*scale_factor2)**2)
        stat3 = numpy.nansum(weights3*(sp1.flux.value-f(sp1.wave.value)*scale_factor3)**2)
        stat = stat1+stat2+stat3
        unit = sp1.flux_unit**2

# absolute deviation
    elif (statistic == 'absdev'):
# compute scale factor
        if scale == True:
            scale_factor1 = numpy.nansum(weights1*sp1.flux.value)/ \
                numpy.nansum(weights1*f(sp1.wave.value))
            scale_factor2 = numpy.nansum(weights2*sp1.flux.value)/ \
                numpy.nansum(weights2*f(sp1.wave.value))
            scale_factor3 = numpy.nansum(weights3*sp1.flux.value)/ \
                numpy.nansum(weights3*f(sp1.wave.value))
# correct variance
        vtot1 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor1**2)],axis=0)
        vtot2 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor2**2)],axis=0)
        vtot3 = numpy.nanmax([sp1.variance.value,v(sp1.wave.value)*(scale_factor3**2)],axis=0)
        vtot1[numpy.isnan(vtot1)==True] = numpy.nanmedian(vtot1)
        vtot2[numpy.isnan(vtot2)==True] = numpy.nanmedian(vtot2)
        vtot3[numpy.isnan(vtot3)==True] = numpy.nanmedian(vtot3)
        stat1 = numpy.nansum(weights1*abs(sp1.flux.value-f(sp1.wave.value)*scale_factor1))
        stat2 = numpy.nansum(weights2*abs(sp1.flux.value-f(sp1.wave.value)*scale_factor2))
        stat3 = numpy.nansum(weights3*abs(sp1.flux.value-f(sp1.wave.value)*scale_factor3))
        stat = stat1+stat2+stat3
        unit = sp1.flux_unit

# error
    else:
        print('Error: statistic {} for compareSpectra not available'.format(statistic))
        return numpy.nan, numpy.nan

# plot spectrum compared to best spectrum
    if plot == True:
        spcomp = copy.deepcopy(sp2)
        spcomp.flux[numpy.where( (spcomp.wave >= numpy.min(fit_ranges[0])) & ((spcomp.wave <= numpy.max(fit_ranges[0])) ) )] *= scale_factor1
        spcomp.flux[numpy.where( (spcomp.wave >= numpy.min(fit_ranges[1])) & ((spcomp.wave <= numpy.max(fit_ranges[1])) ) )] *= scale_factor2
        spcomp.flux[numpy.where( (spcomp.wave >= numpy.min(fit_ranges[2])) & ((spcomp.wave <= numpy.max(fit_ranges[2])) ) )] *= scale_factor3
        spcomp.flux[numpy.where( ( (spcomp.wave < numpy.min(fit_ranges[0])) | (spcomp.wave > numpy.max(fit_ranges[0])) ) &
                                 ( (spcomp.wave < numpy.min(fit_ranges[1])) | (spcomp.wave > numpy.max(fit_ranges[1])) ) &
                                 ( (spcomp.wave < numpy.min(fit_ranges[2])) | (spcomp.wave > numpy.max(fit_ranges[2])) ) 
                                )] = numpy.nan
        kwargs['colors'] = kwargs.get('colors',['k','m','b'])
        kwargs['title'] = kwargs.get('title',sp1.name+' vs '+sp2.name)
        from .plot import plotSpectrum
        plotSpectrum([sp1,spcomp,sp1-spcomp],labels=[sp1.name,sp2.name,'{} = {}'.format(statistic,stat)],**kwargs)

    scale_factors = [scale_factor1, scale_factor2, scale_factor3]

    return numpy.nanmax([stat,minreturn])*unit, scale_factors




def generateMask(wv,mask=[],mask_range=[-99.,-99.],mask_telluric=False,mask_standard=False,**kwargs):
    '''
    :Purpose: Generates a mask array based on wavelength vector and optional inputs on what to mask.

    :Output: A mask array, where 0 = OK and 1 = ignore

    :Example:
    '''

# parameter check
    wave = copy.deepcopy(wv)
    if isinstance(wv,Spectrum): wave = wv.wave
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








def measureEW(sp,lc,width=0.,continuum=[0.,0.],plot=False,file='',continuum_width=False,output_unit=u.Angstrom,nsamp=100,nmc=100,verbose=True,recenter=True,absorption=True,name='',continuum_fit_order=1,debug=False):

# input checks
    if isUnit(lc):
        line_center = lc.to(sp.wave.unit).value
    else: line_center = copy.deepcopy(lc)
    if not isinstance(line_center,int) and not isinstance(line_center,float):
        raise ValueError('Second input value should be a single float number; you entered {}'.format(lc))        
    if numpy.nanmin(sp.wave.value) > line_center or numpy.nanmax(sp.wave.value) < line_center:
        raise ValueError('\nLine {} is outside spectral data limits of {} to {}'.format(line_center*sp.wave.unit,numpy.nanmin(sp.wave),numpy.nanmax(sp.wave)))

    if isUnit(width):
        line_width = width.to(sp.wave.unit).value
    else: line_width = copy.deepcopy(width)
    if not isinstance(line_width,int) and not isinstance(line_width,float):
        raise ValueError('Second input value should be a single float number; you entered {}'.format(width))        
    if float(line_width) == 0.:
        ic = numpy.nanargmin(numpy.array([numpy.abs(a-line_center) for a in sp.wave.value]))
        line_width = numpy.absolute(sp.wave.value[ic+2]-sp.wave.value[ic-2])

    if isUnit(continuum):
        cont = continuum.to(sp.wave.unit).value
    else: cont = copy.deepcopy(continuum)
    if isinstance(cont,int) or isinstance(cont,float):
        cont = [line_width,line_width+numpy.abs(cont)]
    if cont[0] == cont[1]:
        cont = [line_width,line_width*2.]
    if numpy.nanmax(cont) < line_center:
        cont = [c+line_center for c in cont]
    if len(cont) < 4:
        if len(cont) == 2:
            if continuum_width != False:
                cont = [cont[0]-continuum_width/2, cont[0]+continuum_width/2, cont[1]-continuum_width/2,cont[1]+continuum_width/2]
        else:
            cont = [2*line_center-cont[-1],2*line_center-cont[-2],cont[-1],cont[-2]]

    if debug==True: print('Line center = {}, Line width = {}, Continuum = {}'.format(line_center,line_width,cont))

# preset fail condition
    ew                      = numpy.nan
    ew_unc                  = numpy.nan
    line_center_measure     = numpy.nan
    line_center_measure_unc = numpy.nan
    rv                      = numpy.nan
    rv_unc                  = numpy.nan
    
# first compute value
    samplerng = [numpy.nanmin(cont)-0.1*(numpy.nanmax(cont)-numpy.nanmin(cont)),numpy.nanmax(cont)+0.1*(numpy.nanmax(cont)-numpy.nanmin(cont))]
    if numpy.nanmin(sp.wave.value) <= samplerng[0] and numpy.nanmax(sp.wave.value) >= samplerng[1]:

# refine line centering
        line_center_measure = line_center
        if recenter == True:    
            for i in range(5):
                wc = numpy.where(numpy.logical_and(sp.wave.value >= line_center_measure-line_width,sp.wave.value <= line_center_measure+line_width))
                wv = numpy.array(sp.wave.value[wc])
                fl = numpy.array(sp.flux.value[wc])
                if absorption == True: line_center_measure = wv[numpy.nanargmin(fl)]
                else: line_center_measure = wv[numpy.nanargmax(fl)]
        rv = ((line_center_measure-line_center)/line_center)*const.c.to(u.km/u.s)

        w = numpy.where(numpy.logical_and(sp.wave.value >= samplerng[0],sp.wave.value <= samplerng[1]) )
        if len(w[0]) > 0:
            f = interp1d(sp.wave.value[w],sp.flux.value[w],bounds_error=False,fill_value=numpy.nan)
            wline = numpy.linspace(line_center_measure-line_width,line_center_measure+line_width,nsamp)
            wcont = numpy.append(numpy.linspace(cont[0],cont[1],nsamp),numpy.linspace(cont[-2],cont[-1],nsamp))        
            fline = f(wline)
            fcont = f(wcont)

            goodind  = numpy.where(~numpy.isnan(fcont))
            pcont    = numpy.poly1d(numpy.polyfit(wcont[goodind],fcont[goodind],continuum_fit_order))
            fcontfit = pcont(wline)

            ew = (trapz((numpy.ones(len(wline))-(fline/fcontfit)), wline)*sp.wave.unit).to(output_unit)
            if plot == True:
                plt.clf()
                plt.plot(sp.wave.value,sp.flux.value,'k-')
                plt.fill_between(sp.wave.value,sp.flux.value-sp.noise.value,sp.flux.value+sp.noise.value,color='grey',alpha=0.2)
                plt.plot(wline,fline,'m-')
                plt.plot(wcont[wcont<line_center_measure],fcont[wcont<line_center_measure],'b-')
                plt.plot(wcont[wcont>line_center_measure],fcont[wcont>line_center_measure],'b-')
                plt.plot(wline,fcontfit,'b-')
                plt.plot([line_center_measure,line_center_measure],[numpy.nanmin(sp.flux.value[w]),numpy.nanmax(sp.flux.value[w])],'k--')
                plt.xlim(samplerng)
                plt.ylim([numpy.nanmin(sp.flux.value[w]),numpy.nanmax(sp.flux.value[w])])
                plt.title(name+' {:.4f}'.format(line_center_measure))
                plt.xlabel('Wavelength ({})'.format(sp.wave.unit))
                plt.ylabel('Flux Density ({})'.format(sp.flux.unit))
                if file != '': plt.savefig(file)


# MC for errors            
            if numpy.median(sp.noise.value[w]) != numpy.nan and numpy.median(sp.noise.value[w]) != 0.:
# default uncertainty = 2 x median uncertainty --> NOT USING
                ew_unc = (2.*numpy.nanmedian(sp.noise.value[w])/pcont(line_center_measure)*(numpy.nanmax(wline)-numpy.nanmin(wline))*sp.wave.unit).to(output_unit)
                ews,lns,rvs = [],[],[]
                for i in range(nmc):
                    wvv = sp.wave.value[w]
                    flxv = numpy.random.normal(sp.flux.value[w],sp.noise.value[w])
#                    spvar.flux = numpy.random.normal(sp.flux.value,sp.noise.value)*sp.flux.unit
                    line_center_measure_var = line_center
                    if recenter == True:    
                        for i in range(5):
                            wc = numpy.where(numpy.logical_and(wvv >= line_center_measure_var-line_width,wvv <= line_center_measure_var+line_width))
                            wc = numpy.where(numpy.logical_and(wvv >= line_center_measure_var-line_width,wvv <= line_center_measure_var+line_width))
                            wv = numpy.array(wvv[wc])
                            fl = numpy.array(flxv[wc])
                            if absorption == True: line_center_measure_var = wv[numpy.nanargmin(fl)]
                            else: line_center_measure_var = wv[numpy.nanargmax(fl)]
                    lns.append(line_center_measure_var)
                    rvs.append((((line_center_measure_var-line_center)/line_center)*const.c.to(u.km/u.s)).value)

                    f = interp1d(wvv,flxv,bounds_error=False,fill_value=0.)
                    wline = numpy.linspace(line_center_measure_var-line_width,line_center_measure_var+line_width,nsamp)
                    wcont = numpy.append(numpy.linspace(cont[0],cont[1],nsamp),numpy.linspace(cont[-2],cont[-1],nsamp))        
                    fline = f(wline)
                    fcont = f(wcont)
                    pcont = numpy.poly1d(numpy.polyfit(wcont,fcont,continuum_fit_order))
                    fcontfit = pcont(wline)
                    ews.append(((trapz((numpy.ones(len(wline))-(fline/fcontfit)), wline)*sp.wave.unit).to(output_unit)).value)
                rv_unc = numpy.std(rvs)*u.km/u.s
#                ew_unc = numpy.sqrt(ew_unc**2+(numpy.std(ews)*output_unit)**2)
                ew_unc = (numpy.std(ews)*output_unit)
                line_center_measure_unc = numpy.std(lns)

    else:
        print('Spectral range {} to {} does not cover measurement region for EW {} to {}'.format(numpy.nanmin(sp.wave.value),numpy.nanmax(sp.wave.value),samplerng[0],samplerng[1]))
                            
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
        elif el == 'fe i':
            lines = [1.0340886]
        elif el == 'al i':
            lines = [1.3123416,1.3150753,1.6718974,1.6750520]
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
        


def measureEWSet(sp,ref='rojas',**kwargs):
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
    for alts in ['set','ref']: ref = kwargs.get(alts,ref)

# if index information is not passed, then check with INDEX_SET
    tmp = checkDict(ref,EW_SETS)
    if tmp==False: raise ValueError('Index set {} is not currently available'.format(ref))
    ref = copy.deepcopy(tmp)
    info = copy.deepcopy(EW_SETS[ref])

    result = {'reference': info['reference'],'bibcode': info['bibcode']}
    for ftr in list(info['features'].keys()):
        tmp = measureEW(sp,info['features'][ftr]['linecenter'],width=info['features'][ftr]['width'],continuum=info['features'][ftr]['continuum'],continuum_fit_order=info['continuum_fit_order'],recenter=info['features'][ftr]['recenter'],**kwargs)
        result[ftr] = tmp


# # determine combine method
#     if ('rojas' in ref.lower()):
#         reference = 'EW measures from Rojas-Ayala et al. (2012)'
#         names = ['Na I 2.206/2.209','Ca I 2.26']
#         ews = numpy.zeros(len(names))
#         errs = numpy.zeros(len(names))
#         tmp = measureEW(sp,2.206,width=[2.2020, 2.2120],continuum=[2.1965, 2.2125, 2.2175],recenter=False,**kwargs)
#         ews[0] = tmp['ew']
#         errs[0] = tmp['ew_unc']
#         tmp = measureEW(sp,2.2635,width=[2.2580, 2.2690],continuum=[2.2510, 2.2580, 2.2705, 2.2760],recenter=False,**kwargs)
#         ews[0] = tmp['ew']
#         errs[0] = tmp['ew_unc']
#     else:
#         print('{} is not one of the sets used for measureIndexSet'.format(set))
#         return numpy.nan

# # output dictionary of indices
#     result = {names[i]: (ews[i],errs[i]) for i in numpy.arange(len(names))}
#     result['reference'] = reference
#     result['names'] = names
# #    result['reference'] = reference
# #    return inds,errs,names

    return result


def measureIndex(sp,ranges,method='ratio',sample='integrate',weights=[],fitorder=1,nsamples=100,noiseFlag=False,plot=False,verbose=False,**pkwargs):
    '''
    :Purpose: Measure an index on a spectrum based on defined methodology
                measure method can be mean, median, integrate
                index method can be ratio = 1/2, valley = 1-2/3, OTHERS
                output is index value and uncertainty
    .. will also come back to this one

# NOTE: NOISE FLAG ERROR IS BACKWARDS HERE
    '''

# error checking on number of arguments provided
    if len(ranges) == 0:
        raise ValueError('measureIndex needs at least 1 sample region; zero given')
    if len(ranges) == 1:
        method = 'single'
    if not isinstance(ranges[0],list) and not isinstance(ranges[0],numpy.ndarray) and not isinstance(ranges[0],tuple) and not isinstance(ranges[0],set):
        raise ValueError('measureIndex needs a list of wavelength ranges, you entered {} which has type {}'.format(ranges,type(ranges[0])))
    for r in ranges:
        if len(r) != 2:
            raise ValueError('Problem with range {} in input ranges {}: must be 2-element array'.format(r,ranges))
    if (len(ranges) < 2 and (method == 'ratio' or method == 'change')):
        raise ValueError('Index method {} needs at least 2 sample regions'.format(method))
    if (len(ranges) < 3 and (method == 'line' or method == 'allers' or method == 'inverse_line'  or method == 'sumnum' or method == 'sumdenom')):
        raise ValueError('Index method {} needs at least 3 sample regions'.format(method))

# weights
    while len(weights)<len(ranges): weights.append(1.)

# define the sample vectors
    value = numpy.zeros(len(ranges))
    value_sim = numpy.zeros((len(ranges),nsamples))
    xNums,yNums,yNums_e = [],[],[]

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
        xNums.append(xNum)

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
        w = padWhereArray(w,len(sp.wave)-1)
#        print(waveRng,len(w),numpy.min(w),numpy.max(w))
        f = interp1d(sp.wave.value[w],sp.flux.value[w],bounds_error=False,fill_value=numpy.nan)
        yNum = f(xNum)
        yNums.append(yNum)
        if noiseFlag == False:
            s = interp1d(sp.wave.value[w],sp.noise.value[w],bounds_error=False,fill_value=0.)
            yNum_e = s(xNum)
            yNums_e.append(yNum_e)

# first compute the actual value
        if (sample == 'integrate'):
            value[i] = trapz(yNum,xNum)
        elif (sample == 'average'):
            value[i] = numpy.nanmean(yNum)
        elif (sample == 'sum'):
            value[i] = numpy.nansum(yNum)
        elif (sample == 'median'):
            value[i] = numpy.nanmedian(yNum)
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
                    value_sim[i,j] = numpy.nanmedian(yVar)
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
        val = (weights[0]*value[0])/(weights[1]*value[1])
        vals = (weights[0]*value_sim[0,:])/(weights[1]*value_sim[1,:])
    elif (method == 'line'):
        val = 0.5*(weights[0]*value[0]+weights[1]*value[1])/(weights[2]*value[2])
        vals = 0.5*(weights[0]*value_sim[0,:]+weights[1]*value_sim[1,:])/(weights[2]*value_sim[2,:])
    elif (method == 'inverse_line'):
        val = 2.*weights[0]*value[0]/(weights[1]*value[1]+weights[2]*value[2])
        vals = 2.*weights[0]*value_sim[0,:]/(weights[1]*value_sim[1,:]+weights[2]*value_sim[2,:])
    elif (method == 'change'):
        val = 2.*(weights[0]*value[0]-weights[1]*value[1])/(weights[0]*value[0]+weights[1]*value[1])
        vals = 2.*(weights[0]*value_sim[0,:]-weights[1]*value_sim[1,:])/(weights[0]*value_sim[0,:]+weights[1]*value_sim[1,:])
    elif (method == 'sumnum'):
        val = (weights[0]*value[0]+weights[1]*value[1])/(weights[2]*value[2])
        vals = (weights[0]*value_sim[0,:]+weights[1]*value_sim[1,:])/(weights[2]*value_sim[2,:])
    elif (method == 'sumdenom'):
        val = weights[0]*value[0]/(weights[1]*value[1]+weights[2]*value[2])
        vals = weights[0]*value_sim[0,:]/(weights[1]*value_sim[1,:]+weights[2]*value_sim[2,:])
# continuum fitting methods
    elif (method == 'fitnum'):
        cfit = numpy.polyfit(numpy.array([xNums[0],xNums[1]]).ravel(),numpy.array([yNums[0],yNums[1]]).ravel(),fitorder)
        val = weights[0]*numpy.polyval(cfit,numpy.nanmean(xNums[-1]))/(weights[-1]*value[-1])
        if noiseFlag == False: 
            vals = []
            for j in numpy.arange(0,nsamples):
                yVar1 = yNums[0]+numpy.random.normal(0.,3.,size=len(yNums[0]))*yNums_e[0]
                yVar2 = yNums[1]+numpy.random.normal(0.,3.,size=len(yNums[1]))*yNums_e[1]
                cfit = numpy.polyfit(numpy.array([xNums[0],xNums[1]]).ravel(),numpy.array([yVar1,yVar2]).ravel(),fitorder)
                vals.append(weights[0]*numpy.polyval(cfit,numpy.nanmean(xNums[-1]))/(weights[-1]*value_sim[-1,j]))
    elif (method == 'fitdenom'):
        cfit = numpy.polyfit(numpy.array([xNums[1],xNums[2]]).ravel(),numpy.array([yNums[1],yNums[2]]).ravel(),fitorder)
        val = weights[0]*value[0]/(weights[-1]*numpy.polyval(cfit,xNums[-1]))
        if noiseFlag == False: 
            vals = []
            for j in numpy.arange(0,nsamples):
                yVar1 = yNums[0]+numpy.random.normal(0.,3.,size=len(yNums[0]))*yNums_e[0]
                yVar2 = yNums[1]+numpy.random.normal(0.,3.,size=len(yNums[1]))*yNums_e[1]
                cfit = numpy.polyfit(numpy.array([xNums[0],xNums[1]]).ravel(),numpy.array([yVar1,yVar2]).ravel(),fitorder)
                vals.append(weights[0]*value_sim[0,j]/(weights[-1]*numpy.polyval(cfit,xNums[-1])))
    elif (method == 'doubleratio'):
        val = ((weights[0]/weights[1])/(weights[1]/weights[2]))*((value[0]/value[1])/(value[1]/value[2]))
        vals = ((weights[0]/weights[1])/(weights[1]/weights[2]))*((value_sim[0,:]/value_sim[1,:])/(value_sim[1,:]/value_sim[2,:]))
    elif (method == 'allers'):
        val = (((numpy.mean(ranges[0])-numpy.mean(ranges[1]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value[2] \
            + ((numpy.mean(ranges[2])-numpy.mean(ranges[0]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value[1]) \
            /value[0]
        vals = (((numpy.mean(ranges[0])-numpy.mean(ranges[1]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value_sim[2,:] \
            + ((numpy.mean(ranges[2])-numpy.mean(ranges[0]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value_sim[1,:]) \
            /value_sim[0,:]
    else:
        val = (weights[0]*value[0])/(weights[1]*value[1])
        vals = (weights[0]*value_sim[0,:])/(weights[1]*value_sim[1,:])

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





def measureIndexSet(sp,ref='burgasser',index_info={},info=False,verbose=False,indices_keyword='indices',range_keyword='ranges',method_keyword='method',sample_keyword='sample',**kwargs):
    '''
    :Purpose: 

        Measures indices of ``sp`` from specified sets. Returns dictionary of indices.

    :Required Inputs:
    
        * :param sp: Spectrum class object, which should contain wave, flux and noise array elements

    :Optional Inputs:
    
        * :param ref='burgasser': string identifying the index set you want to measure, as contained in the variable splat.INDEX_SETS
        * :param index_info={}: a dictionary defining the index set to measure; this dictionary should contain within it:

            * 'indices': points to a dictionary containing the index definitions, each of which are also a dictionary containing the 
            keywords 'range' (the wavelength range),'sample' (spectrum sampling), and 'method'(index calculation method)
            * alternately, a set of individual dictionaries containing 'range','sample', and 'method' keywords 
            * 'bibcode': optional bibcode citation, if 'indices' keyword is present

        * :param info=False: set to True to give the list of avaialble index sets in SPLAT
        * :param indices_keyword='indices': alternate keyword name convection for indices 
        * :param range_keyword='ranges': alternate keyword name convection for index wavelength ranges
        * :param method_keyword='method': alternate keyword name convection for spectral sampling 
        * :param sample_keyword='sample': alternate keyword name convection for index measurement method
        * :param verbose=False: set to True to give extra feedback

        Additional keyword options for measureIndex()_ can also be supplied

    :Output:

        A dictionary containing the names of each index pointed to tuple of (measurement,uncertainty)

    :Example:

        >>> import splat
        >>> sp = splat.getSpectrum(shortname='1555+0954')[0]
        >>> splat.measureIndexSet(sp, ref = 'reid')
             {'H2O-A': (0.896763035762751, 0.0009482437569842745),
             'H2O-B': (1.0534404728276665, 0.001329359465953495)}

    .. _measureIndex() : api.html#splat.core.measureIndex

    '''

# keyword parameters
    for k in ['reference','set']:
        ref = kwargs.get(k,ref)
    for k in ['indices','index','index_dict']:
        index_info = kwargs.get(k,index_info)
    for k in ['option','options','information','available']:
        info = kwargs.get(k,info)

# just return information on available sets
# NOTE THAT THIS HAS BEEN MOVED TO EMPIRICAL INFO FUNCTION
    if info==True:
        tmp = checkDict(ref,INDEX_SETS)
        if tmp==False: 
            print('Reference {} IS NOT in the available set of indices, which are as follows:\n'.format(ref))
            for k in list(INDEX_SETS.keys()):
                inds = list(INDEX_SETS[k]['indices'].keys())
                bibref = ''
                if 'bibcode' in list(INDEX_SETS[k].keys()): bibref = INDEX_SETS[k]['bibcode']
                r = inds[0]
                for i in inds[1:]: r=r+','+i
                print('{} ({}): indices: {}'.format(k,bibref,r))
        else:
            print('Reference {} IS in the available set of indices'.format(ref))
            if 'bibcode' in list(INDEX_SETS[tmp].keys()): print('Bibcode: {}'.format(INDEX_SETS[tmp]['bibcode']))
            inds = list(INDEX_SETS[tmp]['indices'].keys())
            r = inds[0]
            for i in inds[1:]: r=r+','+i
            print('Indices: {}'.format(r))
        return

# if index information is not passed, then check with INDEX_SET
    if len(index_info) == 0:
        tmp = checkDict(ref,INDEX_SETS)
        if tmp==False: raise ValueError('Index set {} is not currently available'.format(ref))
        ref = copy.deepcopy(tmp)
        index_info = copy.deepcopy(INDEX_SETS[ref])

# check that the relevant keywords are in the index dictionary
    dkeys = list(index_info.keys())
    if indices_keyword not in dkeys:
        for n in dkeys:
            for k in [range_keyword,method_keyword,sample_keyword]:
                if k not in list(index_info[indices_keyword][n].keys()):
                    raise ValueError('Keyword "{}" is missing from top level of dictionary, or keyword "{}" is missing for index {}'.format(indices_keyword,n))
        tmp = {}
        tmp[indices_keyword] = index_info
        index_info = copy.deepcopy(tmp)

    names = list(index_info[indices_keyword])
    for n in names:
        for k in [range_keyword,method_keyword,sample_keyword]:
            if k not in list(index_info[indices_keyword][n].keys()):
                raise ValueError('Keyword "{}" must be present for each defined index in input index dictionary, but is missing for index {}'.format(k,n))

# check on weights

# measure indices
    result = {'reference': ref,'bibcode': INDEX_SETS[ref]['bibcode']}
    for n in names:
        fitorder,weights=1,[]
        if 'weights' in list(index_info[indices_keyword][n].keys()): weights=index_info[indices_keyword][n]['weights']
        if 'fitorder' in list(index_info[indices_keyword][n].keys()): fitorder=index_info[indices_keyword][n]['fitorder']
        ind,err = measureIndex(sp,index_info[indices_keyword][n][range_keyword],method=index_info[indices_keyword][n][method_keyword],sample=index_info[indices_keyword][n][sample_keyword],weights=weights,fitorder=fitorder,verbose=verbose,**kwargs)
        result[n] = (ind,err)
        if verbose == True: print('Index {}: {:.3f}+/-{:.3f}'.format(n,ind,err))

# some modifications
    if ref=='bardalez':
        result['H-dip'] = tuple([i*0.5 for i in result['H-dip']])
        result['J-curve'] = tuple([i*2. for i in result['J-curve']])
    # elif 'rojas' in ref.lower():
    #     ind = result['H2O-K2num'][0]/result['H2O-K2den'][0]
    #     err = [ind*numpy.sqrt((result['H2O-K2num'][1]/result['H2O-K2num'][0])**2+(result['H2O-K2den'][1]/result['H2O-K2den'][0])**2)]
    #     result['H2O-K2'] = (ind,err)
    #     del result['H2O-K2num'], result['H2O-K2den']

    return result




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
        default_info = {'folder': '', 'name': '', 'citation': '', 'bibcode': '', 'altname': [], 'rawfolder': '', 'default': {'teff': 1500, 'logg': 5.0, 'z': 0.}}

# read in folders specified in .splat_spectral_models
    if os.path.exists(os.path.join(HOME_FOLDER,EXTERNAL_SPECTRAL_MODELS_FILE)):
        with open(os.path.join(HOME_FOLDER,EXTERNAL_SPECTRAL_MODELS_FILE), 'r') as frd: x = frd.read()
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
                minfo['altname'] = minfo['altname'].split(',')
#                except:
#                    print('\nWarning: problem reading info.txt file in {}; using default values for model information'.format(minfo['folder']))
            if flag == 0:
                if verbose == True: print('\nAdding {} models to SPLAT model set'.format(minfo['name']))
                SPECTRAL_MODELS[minfo['name']] = copy.deepcopy(minfo)
                del minfo
    return


def theWorks(sp,measure_classification=True,classification_classes=[],classification_plot='',\
    measure_index_set=True,index_sets=list(INDEX_SETS.keys()),\
    measure_index_classification=True,index_class_sets=list(INDEX_CLASSIFICATION_RELATIONS.keys()),class_string=True,\
    measure_gravity=True,\
    measure_spectral_binary=True,spectral_binary_sets=list(SPECTRAL_BINARY_INDICES.keys()),\
    measure_spectrophotometry=True,spectrophotometry_bands=list(FILTERS.keys()),\
    nsamples=100,verbose=ERROR_CHECKING):
    """
    Conducts a full analysis of a spectrum

    Parameters
    ----------

    TBD

    Returns
    -------

    TBD

    Examples
    --------

    TBD

    See Also
    --------

    TBD

    """     

    output = {}

# classification by template 
    if measure_classification==True:
        cout = {}
        if verbose==True: print('\nClassification by standards:')
        if len(classification_classes) == 0: classification_classes = ['all','d','sd','esd','usd','intg','vlg']
        initializeStandards(all=True)
        for cl in classification_classes:
            if len(classification_plot.split('.'))>1:
                suffix = classification_plot.split('.')[-1]
                cfile = classification_plot.replace('.'+suffix,'')+'_'+cl+'.'+suffix
                cout[cl] = classifyByStandard(sp,std_class=cl,telluric=True,method='kirkpatrick',plot=True,output=cfile,colors=['k','m','b'],string=class_string)
            else:
                cout[cl] = classifyByStandard(sp,std_class=cl,telluric=True,method='kirkpatrick',plot=False,colors=['k','m','b'],string=class_string)
            if verbose==True: print('\t{}: {}+/-{:.1f}'.format(cl,cout[cl][0],cout[cl][1]))
        output['spt-template'] = cout

# indices 
    if measure_index_set==True:
        cout = {}
        if verbose==True: print('\nIndices:')
        for iset in index_sets:
            cout[iset] = measureIndexSet(sp,iset)
            if verbose==True: 
                print('\t{}'.format(iset))
                for k in list(cout[iset].keys()): print('\t\t{}: {:.3f}+/-{:.3f}'.format(k,cout[iset][k][0],cout[iset][k][1]))
        output['indices'] = cout

# index-based classifications
    if measure_index_classification==True:
        cout = {}
        if verbose==True: print('\nIndex based classifications:')
        for iset in index_class_sets:
            cout[iset] = classifyByIndex(sp,iset)
            if verbose==True: print('\t{}: {}+/-{:.1f}'.format(iset,cout[iset][0],cout[iset][1]))
        output['spt-indices'] = cout

# gravity
    if measure_gravity==True:
        output['gravity'] = classifyGravity(sp,output='allmeasures',verbose=verbose)

# spectral binary
    if measure_spectral_binary==True:
        cout = {}
        if verbose==True: print('\nSpectral binary:')
        for iset in spectral_binary_sets:
            cout[iset] = classifySB(sp,ref=iset,output='allmeasures',verbose=verbose)
        output['spectral_binary'] = cout
        
# spectrophotometry
    if measure_spectrophotometry==True:
        cout = {}
        if verbose==True: print('\nSpectrophotometry:')
        for band in spectrophotometry_bands:
            cout[band] = sp.filterMag(band)
            if verbose==True: print('\t{}: {:.2f}+/-{:.2f}'.format(band,cout[band][0],cout[band][1]))
        output['spectrophotometry'] = cout
        
    return output








# run test program if calling from command line
if __name__ == '__main__':
    pass
#    test_gooddata()
#    splat.test()
#    test_info()
