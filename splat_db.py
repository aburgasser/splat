"""
.. note::
         These are the database functions for SPLAT 
"""

import astropy
import copy
from datetime import datetime
import os
import re
import requests
import sys
from scipy import stats
import numpy
from astropy.io import ascii            # for reading in spreadsheet
from astropy.table import Table
import splat

SPLAT_DB_FOLDER = '/db/'
SOURCE_DB_FILE = SPLAT_DB_FOLDER+'source_data.txt'
SPECTRAL_DB_FILE = SPLAT_DB_FOLDER+'spectral_data.txt'
PHOTOMETRY_DB_FILE = SPLAT_DB_FOLDER+'photometry_data.txt'
TMPFILENAME = 'splattmpfile'

# change the command prompt
sys.ps1 = 'splat db> '




def assignSourceID(coordinate,**kwargs)
    '''
    Purpose
        Assigns a source ID based on current source database and possible overlap.

    :Required parameters:
        :param coordinate: astropy Coordinate object with coordinates (RA, Dec) of source, to check against current database

    :Optional parameters:
        :param new: Force the return of a new ID number regardless of search
        :type logical: optional, default = False
        :param validate: Ask for validation from user that a source is new
        :type logical: optional, default = False

    :Output:
        - An id key (long integer) for the source

    '''
    pass


def assignSpectralID(**kwargs)
    '''
    Purpose
        Assigns a spectral ID based on current spectral database

    :Required parameters:
        None

    :Optional parameters:
        None

    :Output:
        - An id key (long integer) for the spectrum

    '''
    pass


def bibtext(bibcode,**kwargs)
    '''
    Purpose
        Takes a bibcode and returns a dictionary containing the bibtex information; looks either in internal SPLAT
            or user-supplied bibfile, or seeks online. If nothing found, gives a soft warning and returns False

    :Required parameters:
        :param bibcode: Bibcode string to look up (e.g., '2014ApJ...787..126L')

    :Optional parameters:
        :param bibfile: Filename for bibfile to use in place of SPLAT internal one
        :type string: optional, default = ''
        :param online: If True, go directly online; if False, do not try to go online 
        :type logical: optional, default = null

    :Output:
        - A dictionary containing the bibtext fields, or False if not found

    '''
    pass


def getPhotometry(coordinate,**kwargs)
    '''
    Purpose
        Downloads photometry for a source using astroquery (?)

    :Required parameters:
        :param coordinate: astropy Coordinate object with coordinates (RA, Dec) of source to be searched

    :Optional parameters:
 
        :param radius: radius in arcseconds for matching
        :type float: optional, default = 5.
        :param validate: Ask for validation from user for matched photometry
        :type logical: optional, default = False
        :param 2MASS: Download 2MASS photometry
        :type logical: optional, default = True
        :param SDSS: Download 2MASS photometry
        :type logical: optional, default = True
        :param UKIDSS: Download 2MASS photometry
        :type logical: optional, default = True
        :param WISE: Download 2MASS photometry
        :type logical: optional, default = True
        :param DENIS: Download 2MASS photometry
        :type logical: optional, default = True

    :Output:
        - A table? filename? if new photometry

    '''
    pass




def importSpectra(input,**kwargs)
    '''
    Purpose
        imports a set of spectra into the SPLAT library; requires manager access.

    :Required parameters:
        :param input: Can be one of the following:
            - A string or array of strings containing spectra filenames, which must either be fits files
                or an ascii files with three columns (tab- or comma-delimited) containing wavelength, flux and noise.
            - A string filename for a spreadsheet (ascii, tab- or comma-delimited) listing the input spectra, one per row.
                At least one column must be ``filename``; the remaining may be one of the optional inputs listed below

    :Optional parameters:
        :param source_key: Source identification key if a source is already in the SPLAT database; use this for a *new* observation of a previously observed source
        :type long: optional, default = null
        :param spectrum_key: Spectrum identification if a spectrum is already in the SPLAT database; use this for a *replacing* a previous observation
        :type long: optional, default = null
        :param name: Source name
        :type string: optional, default = ''
        :param ra: Source Right Ascension (used for searching for catalog information)
        :type double: optional, default = null
        :param dec: Source Declination (used for searching for catalog information)
        :type double: optional, default = null
        :param designation: Source designation (e.g., 'J1234567+0123456'; used for searching for catalog information)
        :type string: optional, default = ''
        :param opt_type: Source optical spectral type (e.g., 'M7')
        :type string: optional, default = ''
        :param opt_type_reference: Bibcode for optical type reference (e.g., '2014ApJ...787..126L')
        :type string: optional, default = ''
        :param nir_type: Source near-infrared spectral type (e.g., 'M7')
        :type string: optional, default = ''
        :param nir_type_reference: Bibcode for optical type reference (e.g., '2014ApJ...787..126L')
        :type string: optional, default = ''
        :param publication_reference: Publication reference bibcode (e.g., '2014ApJ...787..126L')
        :type string: optional, default = ''
        :param observer: Name of observer
        :type string: optional, default = ''
        :param date_observed: Observation date (e.g., '20100425', '2010 Apr 25')
        :type string: optional, default = ''

    :Output:
        - Source DB update file: spreadsheet containing update to source_data.txt, saved locally as UPDATE_source_data.txt
        - Spectral DB update file: spreadsheet containing update to spectral_data.txt, saved locally as UPDATE_spectral_data.txt
        - Photometry DB update file: spreadsheet containing update to photometry_data.txt, saved locally as UPDATE_photometry_data.txt

    '''
    pass




# main testing of program
if __name__ == '__main__':
    basefolder = '/Users/adam/projects/splat/exercises/ex9/'
    sp = splat.getSpectrum(shortname='1047+2124')[0]        # T6.5 radio emitter
    spt,spt_e = splat.classifyByStandard(sp,spt=['T2','T8'])
    teff,teff_e = splat.typeToTeff(spt)
    sp.fluxCalibrate('MKO J',splat.typeToMag(spt,'MKO J')[0],absolute=True)
    table = modelFitMCMC(sp, mask_standard=True, initial_guess=[teff, 5.3, 0.], zstep=0.1, nsamples=100,savestep=0,filebase=basefolder+'fit1047',verbose=True)

