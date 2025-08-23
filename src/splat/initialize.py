# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
.. note::
         These commands initialize the SPLAT code 
"""

import os
import sys
import numpy
from astropy import units as u
import pandas


# reference parameters
VERSION = '2025.08.23'
__version__ = VERSION
SPLAT_URL = 'http://splat.physics.ucsd.edu/splat/'
DOCUMENTATION_URL = 'http://splat.physics.ucsd.edu/splat/'
GITHUB_URL = 'https://github.com/aburgasser/splat/'
ZENODO_URL = 'https://github.com/aburgasser/splat/'
DOI = '10.48550/arXiv.1707.00062'
BIBCODE = '2017ASInC..14....7B'
CITATION = 'Burgasser et al. (2017, Astro. Soc. India Conf. Series 14, p. 7)'
EMAIL = 'aburgasser@gmail.com'
TMPFILENAME = 'splattmpfile'
MONTHS = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ERROR_CHECKING = False

# SPLAT authors
AUTHORS = [
    'Christian Aganze',
    'Daniella Bardalez Gagliuffi',
    'Jessica Birky',
    'Adam Burgasser (PI)',
    'Caleb Choban',
    'Andrew Davis',
    'Ivanna Escala',
    'Joshua Hazlett',
    'Carolina Herrara Hernandez',
    'Aishwarya Iyer',
    'Yuhui Jin',
    'Michael Lopez',
    'Dorsa Majidi',
    'Diego Octavio Talavera Maya',
    'Alex Mendez',
    'Gretel Mercado',
    'Niana Mohammed',
    'Elizabeth Hilario Moreno',
    'Jonathan Parra',
    'Maitrayee Sahi',
    'Adrian Suarez',
    'Melisa Tallis',
    'Chris Theissen',
    'Tomoki Tamiya',
    'Steven Truong',
    'Russell Van Linge',
    'Joman Wong',
]

# folders
CODE_PATH = os.path.dirname(os.path.abspath(__file__))
HOME_FOLDER = os.path.expanduser('~')
RESOURCES_FOLDER = os.path.join(CODE_PATH,'resources/')
DB_FOLDER = os.path.join(CODE_PATH,'db/')
#DOCS_FOLDER = os.path.join(CODE_PATH,'../../docs/')
FILTER_FOLDER = os.path.join(RESOURCES_FOLDER,'Filters/')
SPECTRAL_DATA_FOLDER = os.path.join(RESOURCES_FOLDER,'Spectra/Public/SPEX-PRISM')
SPECTRAL_MODEL_FOLDER = os.path.join(RESOURCES_FOLDER,'SpectralModels/')
EVOLUTIONARY_MODEL_FOLDER = os.path.join(RESOURCES_FOLDER,'EvolutionaryModels/')
TELLURIC_MODEL_FOLDER = os.path.join(RESOURCES_FOLDER,'Telluric/')
CITATION_RESOURCES_FOLDER = os.path.join(RESOURCES_FOLDER,'Citations/')
EXTERNAL_DATA_FILE = os.path.join(HOME_FOLDER,'.splat_data')
EXTERNAL_SPECTRAL_MODELS_FILE = os.path.join(HOME_FOLDER,'.splat_spectral_models')
EXTERNAL_EVOLUTIONARY_MODELS_FILE = os.path.join(HOME_FOLDER,'.splat_evolutionary_models')
#ACCESS_FILE = '.splat_access'

# set up external files if not present
for x in [EXTERNAL_DATA_FILE,EXTERNAL_SPECTRAL_MODELS_FILE,EXTERNAL_EVOLUTIONARY_MODELS_FILE]:
    if os.path.exists(x)==False:
        try:
            with open(x, 'w') as fp:
                pass
        except: print('Warning: cannot create backup folder {}'.format(x))

#set user SPLAT model path from environmental variable
# SPLAT_USER_MODELS = './'
# if os.environ.get('SPLAT_MODELS') != None:
#     SPLAT_USER_MODELS = os.environ['SPLAT_MODELS']

# Unit standards
DEFAULT_WAVE_UNIT = u.micron
DEFAULT_FLUX_UNIT = u.erg/u.s/u.cm/u.cm/u.micron
DEFAULT_SED_UNIT = u.erg/u.s/u.cm/u.cm
DEFAULT_CROSSSECTION_UNIT = u.cm**2
DEFAULT_INSTRUMENT = 'SPEX-PRISM'

# OLD SPECTRAL LIBRARY PARAMETERS
#DATA_FOLDER = os.path.join(SPECTRAL_DATA_FOLDER,DEFAULT_INSTRUMENT,'')
#DEFAULT_DATA_FOLDER = os.path.join(SPECTRAL_DATA_FOLDER,/')
DB_SPECTRA_FILE = os.path.join(DB_FOLDER,'spectral_data.txt')
DB_SOURCES_FILE = os.path.join(DB_FOLDER,'source_data.txt')
# read these in
DB_SOURCES = pandas.read_csv(os.path.normpath(DB_SOURCES_FILE),delimiter='\t')
DB_SPECTRA = pandas.read_csv(os.path.normpath(DB_SPECTRA_FILE),delimiter='\t')

# NEW SPECTRAL LIBRARY PARAMETERS
DB_SPECTRA_INPUT_FILE = 'spectra.csv'
DB_SPECTRA_DEFAULT_PARAMETERS = {
    'FILENAME': {'altname': ['FILE','DATAFILE','DATA_FILE','FILE_NAME'], 'default': '', 'type': str, 'initialize': 'required'},
#    'DATA_KEY': {'altname': [], 'type': int, 'required': False},
    'FOLDER': {'altname': ['PATH','FOLD','DIR','DIRECTORY'], 'default': '', 'type': str, 'initialize': 'create'},
    'INSTRUMENT': {'altname': [], 'default': DEFAULT_INSTRUMENT, 'type': str, 'initialize': 'default'},
    'BIBCODE': {'altname': ['REFERENCE','DATA_BIBCODE','REF','BIB','DATA_REFERENCE','DATA_REF'], 'default': '', 'type': str, 'initialize': 'default'},
    'OBSERVATION_DATE': {'altname': ['OBSDATE','OBS_DATE','DATE'], 'default': '', 'type': str, 'initialize': 'default'},
    'OBSERVATION_TIME': {'altname': ['OBSTIME','OBS_TIME'], 'default': numpy.nan, 'type': float, 'initialize': 'default'},
    'OBSERVATION_MJD': {'altname': ['MJD','OBS_MJD','JULIAN_DATE'], 'default': numpy.nan, 'type': float, 'initialize': 'default'},
    'PROGRAM': {'altname': ['PROGRAM_ID','PROGRAM_NUMBER'], 'default': '', 'type': str, 'initialize': 'default'},
    'OBSERVER': {'altname': ['PI'], 'default': '', 'type': str, 'initialize': 'default'},
    'AIRMASS': {'altname': ['Z'], 'default': numpy.nan, 'type': float, 'initialize': 'default'},
    'INTEGRATION': {'altname': ['TINT','TIME','INT_TIME','INTEGRATION_TIME','EXPOSURE','EXPTIME','EXP_TIME','ITIME'], 'default': numpy.nan, 'type': float, 'initialize': 'default'},
    'CONDITIONS': {'altname': ['WEATHER'], 'default': '', 'type': str, 'initialize': 'default'},
    'NAME': {'altname': ['SOURCE','SOURCE_NAME'], 'default': '', 'type': str, 'initialize': 'default'},
}
DB_SPECTRA_INPUT_FILE_KEY = 'FILENAME'
SPECTRA_FILES_EXTENSIONS = ['txt','asc','fit','fits','dat','gz','bz2','zip']

DB_SOURCES_INPUT_FILE = 'sources.csv'
DB_SOURCES_DEFAULT_PARAMETERS = {
    'NAME': {'altname': ['SOURCE_KEY','SOURCE_NAME'], 'type': int, 'required': False},
    'DESIGNATION': {'altname': ['DESIG'], 'type': str, 'required': False},
    'RA': {'altname': ['RIGHT_ASCENSION'], 'type': float, 'required': False},
    'DEC': {'altname': ['DECLINATION'], 'type': float, 'required': False},
}
DB_SOURCES_INPUT_FILE_KEY = 'NAME'

# Library constants
# #LIBRARY_DATA_FOLDERS = [LIBRARY_PUBLIC_FOLDER]
# LIBRARY_DEFAULT_COLUMNS = {
#     'FILENAME': {'altname': ['FILE','FILES','FILENAMES','PATH','PATHS'], 'type': 'required', 'default': '', 'error': False},
#     'INSTRUMENT': {'altname': ['INST','INSTR','INSTRUMENTS'], 'type': 'add', 'default': DEFAULT_INSTRUMENT, 'error': False},
#     'LIBRARY': {'altname': ['LIB','LIBRARIES','CATALOG','CAT','CATALOGS'], 'type': 'add', 'default': '', 'error': False},
#     'DESIGNATION': {'altname': ['DES','DESIG','DESIGNATIONS'], 'type': 'optional', 'default': '', 'error': False},
#     'COORDINATE': {'altname': ['COORD','COORDINATES'], 'type': 'optional', 'default': '', 'error': False},
#     'RA': {'altname': ['RAS','R','RIGHT_ASCENSION'], 'type': 'optional', 'default': 0., 'error': True},
#     'DEC': {'altname': ['DECS','D','DECLINATION'], 'type': 'optional', 'default': 0., 'error': True},
#     'NAME': {'altname': ['NAMES'], 'type': 'optional', 'default': '', 'error': False},
#     'SHORTNAME': {'altname': ['SHORTNAMES','SNAME','SNAMES'], 'type': 'optional', 'default': '', 'error': False},
#     'SPECTRUM': {'altname': ['SPEC','SPECTRA'], 'type': 'add', 'default': None, 'error': False},
#     'BIBCODE': {'altname': ['DATA_BIBCODE','REFERENCE','REF','BIBCODE','BIB','DATA_REFERENCE','DATA_REF','CITE','CITATION'], 'type': 'optional', 'default': '', 'error': False},
#     'OBSERVATION_DATE': {'altname': ['OBSDATE','OBS_DATE','DATE'], 'type': 'optional', 'default': '', 'error': False},
#     'OBSERVATION_MJD': {'altname': ['JD','OBS_JD','JULIAN_DATE','JULIAN','MJD','MODIFIED_JULIAN_DATE','MODIFIED_JD','OBS_MJD','JDATE','MJDATE'], 'type': 'optional', 'default': 0., 'error': False},
#     'SOURCE_ID': {'altname': ['ID','SRC','SRC_ID','SOURCE'], 'type': 'optional', 'default': -1, 'error': False},
#     'SPT': {'altname': ['SPECTRAL_TYPE','TYPE','SPECTYPE','STYPE'], 'type': 'optional', 'default': 0., 'error': True},
#     'OSPT': {'altname': ['OPTICAL_SPECTRAL_TYPE','OPTICAL_TYPE','OPT_TYPE','OPT_SPT','TYPE_OPT','SPT_OPT','SPECTRAL_TYPE_OPT','SPECTRAL_TYPE_OPTICAL'], 'type': 'optional', 'default': 0., 'error': True},
#     'NSPT': {'altname': ['NIR_SPECTRAL_TYPE','NIR_TYPE','NIR_TYPE','NIR_SPT','TYPE_NIR','SPT_NIR','SPECTRAL_TYPE_NIR','IR_SPECTRAL_TYPE','IR_TYPE','IR_TYPE','IR_SPT','TYPE_IR','SPT_IR','SPECTRAL_TYPE_IR','SPECTRAL_TYPE_IR','IRSPT','INFRARED_SPECTRAL_TYPE','INFRARED_TYPE','INFRARED_TYPE','INFRARED_SPT','TYPE_INFRARED','SPT_INFRARED','SPECTRAL_TYPE_INFRARED','SPECTRAL_TYPE_INFRARED','INFRAREDSPT'], 'type': 'optional', 'default': 0., 'error': True},
#     'SSPT': {'altname': ['SIMBAD_SPECTRAL_TYPE','SIMBAD_TYPE','SIM_TYPE','SIM_SPT','TYPE_SIM','SPT_SIM','SPECTRAL_TYPE_SIM','SPECTRAL_TYPE_SIMBAD'], 'type': 'optional', 'default': 0., 'error': True},
#     # 'RV': {'altname': ['RADIAL_VELOCITY','VRAD','VR','RADVEL','RAD_VEL','RAD_V'], 'type': 'optional', 'default': 0., 'error': True},
#     # 'VSINI': {'altname': ['ROTATIONAL_VELOCITY','VROT','ROTVEL','ROT_VEL','ROT_V'], 'type': 'optional', 'default': 0., 'error': True},
#     # 'PLX': {'altname': ['PARALLAX','PI'], 'type': 'optional', 'default': 0., 'error': True},
#     # 'MU': {'altname': ['PROPER_MOTION','PM'], 'type': 'optional', 'default': 0., 'error': True},
#     # 'MU_RA': {'altname': ['RA_PROPER_MOTION','PROPER_MOTION_RA','MU_RA','MUR','MU_R'], 'type': 'optional', 'default': 0., 'error': True},
#     # 'MU_DEC': {'altname': ['DEC_PROPER_MOTION','PROPER_MOTION_DEC','MU_DEC','MUD','MU_D'], 'type': 'optional', 'default': 0., 'error': True},
# }
# # also magnitude convention: MAG_[FILTER] where filter is preferably one listed below
# LIBRARY_DEFAULT_DATABASE_FILE = 'spectra.csv'
# LIBRARY_DEFAULT_DATA_FOLDER = 'spectra/'
# LIBRARY = {}
# SOURCES_DEFAULT_DATABASE_FILE = 'sources.csv'
# SOURCES = {}


# dwarf spectral standards
STDS_DWARF_SPEX_KEYS = { \
    'M0.0': 11335, #'11335_10505.fits',\
    'M1.0': 11364, #'11364_10806.fits',\
    'M2.0': 11181, #'11181_10187.fits',\
    'M3.0': 10823, #'10823_11422.fits',\
    'M4.0': 12004, #'12004_10444.fits',\
    'M5.0': 10829, #'10829_10104.fits',\
    'M6.0': 11182, #'11182_10188.fits',\
    'M7.0': 10822, #'10822_11283.fits',\
    'M8.0': 10824, #'10824_11423.fits',\
    'M9.0': 10821, #'10821_11058.fits',\
    'L0.0': 10107, #'10107_10315.fits',\
    'L1.0': 11072, #'11072_11527.fits',\
    'L2.0': 10600, #'10600_10957.fits',\
    'L3.0': 10592, #'10592_11111.fits',\
    'L4.0': 10675, #'10675_11572.fits',\
    'L5.0': 10351, #'10351_10583.fits',\
    'L6.0': 10375, #'10375_10696.fits',\
    'L7.0': 10678, #'10678_10105.fits',\
    'L8.0': 10115, #'10115_11254.fits',\
    'L9.0': 10268, #'10268_10237.fits',\
    'T0.0': 10771, #'10771_10871.fits',\
    'T1.0': 10767, #'10767_10591.fits',\
    'T2.0': 10017, #'10017_10945.fits',\
    'T3.0': 10034, #'10034_10874.fits',\
    'T4.0': 10143, #'10143_11632.fits',\
    'T5.0': 10021, #'10021_11106.fits',\
    'T6.0': 10200, #'10200_11236.fits',\
    'T7.0': 10159, #'10159_10513.fits',\
    'T8.0': 10126, #'10126_10349.fits',\
    'T9.0': 11536, #'11536_10509.fits'}
    }

# subdwarf spectral standards
STDS_DSD_SPEX_KEYS = { \
    'd/sdM4.0': 10523, # LSPM J0713+2151
    'd/sdM5.0': 10232, # 2MASS J2059203+175223
    'd/sdM6.0': 10198, # LSR 1610-0040
    'd/sdM7.0': 10863, # NLTT 57956
    'd/sdM8.0': 10040, # 2MASS J15561873+1300527
#    'd/sdM9.0': 10367, # SSSPM 1444-2019
    'd/sdL0.0': 10146, # 2MASS J00412179+3547133
    'd/sdL1.0': 10506, # 2MASS J17561080+2815238 Greco+2019 sdL1
    'd/sdL3.0': 12230, # WISEA J003338.45+282732.4
    'd/sdL6.0': 11509, # SDSS J133148.92-011651.4 Zhang+2018 peculiar L6
    'd/sdL7.0': 10552, # 2MASS J11582077+0435014 Burgasser+2025 d/sdL7 Greco+2019 sdL7
    'd/sdT8.0': 10185, # 2MASS J09393548-2448279 Burgasser+2025 d/sdT8
    }

STDS_SD_SPEX_KEYS = { \
    'sdM2.0': 10223, # LHS 3181 Greco+2019 sdM2
    'sdM4.0': 10528, # LSPM J0949+1746
    'sdM5.0': 10221, # LHS 407 Greco+2019 sdM5
    'sdM5.5': 11670, # APMPM 1523-0245
    'sdM6.0': 10265, # LHS 1074 Greco+2019 sdM6
    'sdM7.0': 10197, # LHS 377 Greco+2019 sdM7
    'sdM8.5': 10123, # 2MASS J01423153+0523285 Greco+2019 sdM8.5
#    'sdM9.5': 10188, # SSSPM 1013-1356 
    'sdL0.0': 11973, # WISE J04592121+1540592
    'sdL1.0': 10367, # SSSPM 1444-2019 Zhang+2017 esdL1
    'sdL5.0': 11240, # SDSS J1416+1348A Zhang+2017 sdL7
    'sdT6.0': 10171, # 2MASS J0937+2931 Burgasser+2025 sdT6
    'sdT7.0': 11377, # ULAS J1416+1348B Burgasser+2025 sdT7
    }

# extreme subdwarf spectral standards
STDS_ESD_SPEX_KEYS = { \
#    'esdM0.0': 10763, # LHS 217 Greco+2019 esdM0
    'esdM4.0': 10366, # LHS 375 Greco+2019 esdM4
    'esdM5.0': 10229, # LP 589-7 Greco+2019 esdM5
    'esdM6.0': 10359, # LHS 2023 Greco+2019 esdM6
    'esdM7.0': 10521, # APMPM J0559-2903  Greco+2019 esdM7
    'esdM8.5': 10278, # LEHPM 2-59  Greco+2019 esdM8
    'esdL0.0': 10188, # SSSPM J1013-1356 Zhang+2017 usdL0 Greco+2019 sdM9.5
    'esdL3.0': 10364, # SDSS J125637.16-022452.2 Zhang+2017 usdL3 Greco+2019 sdL3.5
    'esdL4.0': 10203, # 2MASS J16262034+3925190 Zhang+2017 usdL4 Greco+2019 sdL4
    }

# young spectral standards from Allers & Liu (2013) and Cruz et al. (2017)
STDS_VLG_SPEX_KEYS = { \
#    'M5.0': 10228,  # TWA11C - don't have this
#    'M6.0': 10228,  # TWA8B - don't have this
    'M7.0gamma': 11199,  # 0335+2342
    'M8.0gamma': 10799,  # TWA27A
    'M9.0gamma': 10797,  # TWA26
    'L0.0gamma': 10228,  # 0141-4633 - proposed by Cruz et al. 2017
    'L1.0gamma': 11178,  # 0518-2756
    'L2.0gamma': 11696,  # 0536-1920
    'L3.0gamma': 11198,  # 2208+2921
    'L4.0gamma': 11177,  # 0501-0010 - proposed by Cruz et al. 2017
    'L6.0gamma': 10455,  # 2244+2043
    }

STDS_INTG_SPEX_KEYS = { \
    'M8.0beta': 12572,  # 0019+4614
    'L0.0beta': 11113,  # 1552+2948
    'L1.0beta': 10845,  # 0227-1624 - proposed by Cruz et al. 2017
    'L2.0beta': 11304,  # 0602+3910
    'L3.0beta': 11070,  # 1726+1538
    'L6.0beta': 10678,  # 0103+1935
    }

# holding arrays
SPECTRA_READIN = {}
STDS_DWARF_SPEX = {}
STDS_SD_SPEX = {}
STDS_DSD_SPEX = {}
STDS_ESD_SPEX = {}
STDS_VLG_SPEX = {}
STDS_INTG_SPEX = {}

# filters
# this information is from the SVO filter profile service: http://svo2.cab.inta-csic.es/svo/theory/fps/
FILTERS = { \
    '2MASS_J': {'file': 'j_2mass.txt.gz', 'description': '2MASS J-band', 'zeropoint': 1594.0, 'method': 'vega', 'rsr': True, 'altname': []}, \
    '2MASS_H': {'file': 'h_2mass.txt.gz', 'description': '2MASS H-band', 'zeropoint': 1024.0, 'method': 'vega', 'rsr': True, 'altname': []}, \
    '2MASS_KS': {'file': 'ks_2mass.txt.gz', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega', 'rsr': True, 'altname': ['2MASS_K']}, \
#    '2MASS_K': {'file': 'ks_2mass.txt.gz', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega'}, \
#    '2MASS_Ks': {'file': 'ks_2mass.txt.gz', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega'}, \
    'BESSEL_U': {'file': 'Bessel_U.txt.gz', 'description': 'Bessel U-band', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': ['U']}, \
    'BESSEL_B': {'file': 'Bessel_B.txt.gz', 'description': 'Bessel B-band', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': ['B']}, \
    'BESSEL_V': {'file': 'Bessel_V.txt.gz', 'description': 'Bessel V-band', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': ['V']}, \
    'BESSEL_R': {'file': 'Bessel_R.txt.gz', 'description': 'Bessel R-band', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': ['R']}, \
    'BESSEL_I': {'file': 'Bessel_I.txt.gz', 'description': 'Bessel I-band', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': ['I']}, \
    'COUSINS_I': {'file': 'i_cousins.txt.gz', 'description': 'Cousins I-band', 'zeropoint': 2405.3, 'method': 'vega', 'rsr': False, 'altname': ['IC']}, \
    'DECAM_U': {'file': 'DECam_u.txt.gz', 'description': 'DECam u-band', 'zeropoint': 1568.5, 'method': 'vega', 'rsr': False, 'altname': ['DECCAM_U','DEC_U']}, \
    'DECAM_G': {'file': 'DECam_g.txt.gz', 'description': 'DECam g-band', 'zeropoint': 3909.11, 'method': 'vega', 'rsr': False, 'altname': ['DECCAM_G','DEC_G']}, \
    'DECAM_R': {'file': 'DECam_r.txt.gz', 'description': 'DECam r-band', 'zeropoint': 3151.44, 'method': 'vega', 'rsr': False, 'altname': ['DECCAM_R','DEC_R']}, \
    'DECAM_I': {'file': 'DECam_i.txt.gz', 'description': 'DECam i-band', 'zeropoint': 2584.6, 'method': 'vega', 'rsr': False, 'altname': ['DECCAM_I','DEC_I']}, \
    'DECAM_Z': {'file': 'DECam_z.txt.gz', 'description': 'DECam z-band', 'zeropoint': 2273.09, 'method': 'vega', 'rsr': False, 'altname': ['DECCAM_Z','DEC_Z']}, \
    'DECAM_Y': {'file': 'DECam_y.txt.gz', 'description': 'DECam y-band', 'zeropoint': 2205.95, 'method': 'vega', 'rsr': False, 'altname': ['DECCAM_Y','DEC_Y']}, \
    'DECAM_VR': {'file': 'DECam_vr.txt.gz', 'description': 'DECam z-band', 'zeropoint': 4000., 'method': 'vega', 'rsr': False, 'altname': ['DECCAM_VR','DEC_VR']}, \
    'DES_U': {'file': 'DES_u.txt.gz', 'description': 'DES u-band (filter + atm)', 'zeropoint': 1568.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'DES_G': {'file': 'DES_g.txt.gz', 'description': 'DES g-band (filter + atm)', 'zeropoint': 3909.11, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'DES_R': {'file': 'DES_r.txt.gz', 'description': 'DES r-band (filter + atm)', 'zeropoint': 3151.44, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'DES_I': {'file': 'DES_i.txt.gz', 'description': 'DES i-band (filter + atm)', 'zeropoint': 2584.6, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'DES_Z': {'file': 'DES_z.txt.gz', 'description': 'DES z-band (filter + atm)', 'zeropoint': 2273.09, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'DES_Y': {'file': 'DES_y.txt.gz', 'description': 'DES y-band (filter + atm)', 'zeropoint': 2205.95, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'EUCLID_VIS': {'file': 'Euclid_VIS.txt.gz', 'description': 'Euclid VIS-band', 'zeropoint': 2666.99, 'method': 'ab', 'rsr': False, 'altname': ['EVIS','E_VIS','EUC_VIS']}, \
    'EUCLID_Y': {'file': 'Euclid_Y.txt.gz', 'description': 'Euclid Y-band', 'zeropoint': 1898.78, 'method': 'ab', 'rsr': False, 'altname': ['EY','E_Y','EUC_Y']}, \
    'EUCLID_J': {'file': 'Euclid_J.txt.gz', 'description': 'Euclid H-band', 'zeropoint': 1353.80, 'method': 'ab', 'rsr': False, 'altname': ['EJ','E_J','EUC_J']}, \
    'EUCLID_H': {'file': 'Euclid_H.txt.gz', 'description': 'Euclid K-band', 'zeropoint': 921.81, 'method': 'ab', 'rsr': False, 'altname': ['EH','E_H','EUC_H']}, \
    'FOURSTAR_J': {'file': 'fourstar-j.txt.gz', 'description': 'FOURSTAR J-band', 'zeropoint': 1581.2, 'method': 'vega', 'rsr': False, 'altname': ['4star j']}, \
    'FOURSTAR_J1': {'file': 'fourstar-j1.txt.gz', 'description': 'FOURSTAR J1-band', 'zeropoint': 1978.7, 'method': 'vega', 'rsr': False, 'altname': ['4star j1']}, \
    'FOURSTAR_J2': {'file': 'fourstar-j2.txt.gz', 'description': 'FOURSTAR J2-band', 'zeropoint': 1774.5, 'method': 'vega', 'rsr': False, 'altname': ['4star j2']}, \
    'FOURSTAR_J3': {'file': 'fourstar-j3.txt.gz', 'description': 'FOURSTAR J3-band', 'zeropoint': 1488.8, 'method': 'vega', 'rsr': False, 'altname': ['4star j3']}, \
    'FOURSTAR_H': {'file': 'fourstar-h.txt.gz', 'description': 'FOURSTAR H-band', 'zeropoint': 1054.9, 'method': 'vega', 'rsr': False, 'altname': ['4star h']}, \
    'FOURSTAR_H_SHORT': {'file': 'fourstar-hshort.txt.gz', 'description': 'FOURSTAR H short', 'zeropoint': 1119.1, 'method': 'vega', 'rsr': False, 'altname': ['4star h short','4star h-short','4star hs','fourstar hs','fourstar h1']}, \
    'FOURSTAR_H_LONG': {'file': 'fourstar-hlong.txt.gz', 'description': 'FOURSTAR H long', 'zeropoint': 980.7, 'method': 'vega', 'rsr': False, 'altname': ['4star h long','4star h-long','4star hl','fourstar hl','fourstar h2']}, \
    'FOURSTAR_KS': {'file': 'fourstar-ks.txt.gz', 'description': 'FOURSTAR Ks-band', 'zeropoint': 675.7, 'method': 'vega', 'rsr': False, 'altname': ['4star k','4star ks','fourstar k']}, \
    'FOURSTAR_1.18': {'file': 'fourstar-118.txt.gz', 'description': 'FOURSTAR 1.18 micron narrow band', 'zeropoint': 675.7, 'method': 'vega', 'rsr': False, 'altname': ['4star 1.18','4star 118','fourstar 118']}, \
    'FOURSTAR_2.09': {'file': 'fourstar-209.txt.gz', 'description': 'FOURSTAR 2.09 micron narrow band', 'zeropoint': 675.7, 'method': 'vega', 'rsr': False, 'altname': ['4star 2.09','4star 209','fourstar 209']}, \
    'GAIA_G': {'file': 'GAIA_G.txt.gz', 'description': 'GAIA G-band', 'zeropoint': 3534.7, 'method': 'vega', 'rsr': False, 'altname': ['gaia']}, \
    'GAIA_B': {'file': 'GAIA_Bp.txt.gz', 'description': 'GAIA Bp-band', 'zeropoint': 3296.2, 'method': 'vega', 'rsr': False, 'altname': ['gaia-bp']}, \
    'GAIA_R': {'file': 'GAIA_Rp.txt.gz', 'description': 'GAIA Rp-band', 'zeropoint': 2620.3, 'method': 'vega', 'rsr': False, 'altname': ['gaia-rp']}, \
    'HAWK_Y': {'file': 'hawk-y.txt.gz', 'description': 'HAWK Y-band', 'zeropoint': 2092.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_J': {'file': 'hawk-j.txt.gz', 'description': 'HAWK J-band', 'zeropoint': 1543.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_H': {'file': 'hawk-h.txt.gz', 'description': 'HAWK H-band', 'zeropoint': 1053.6, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_H2': {'file': 'hawk-h2.txt.gz', 'description': 'HAWK H2-band', 'zeropoint': 688.8, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_CH4': {'file': 'hawk-ch4.txt.gz', 'description': 'HAWK CH4-band', 'zeropoint': 1093.4, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_KS': {'file': 'hawk-ks.txt.gz', 'description': 'HAWK Ks-band', 'zeropoint': 675.3, 'method': 'vega', 'rsr': False, 'altname': ['hawk k']}, \
    'HAWK_BRG': {'file': 'hawk-brg.txt.gz', 'description': 'HAWK Brackett Gamma', 'zeropoint': 638.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_NB1060': {'file': 'hawk-nb1060.txt.gz', 'description': 'HAWK Narrow Band 1060', 'zeropoint': 2003.27, 'method': 'vega', 'rsr': False, 'altname': ['hawk 1060']}, \
    'HAWK_NB1190': {'file': 'hawk-nb1190.txt.gz', 'description': 'HAWK Narrow Band 1190', 'zeropoint': 1697.50, 'method': 'vega', 'rsr': False, 'altname': ['hawk 1190']}, \
    'HAWK_NB2090': {'file': 'hawk-nb2090.txt.gz', 'description': 'HAWK Narrow Band 2090', 'zeropoint': 706.68, 'method': 'vega', 'rsr': False, 'altname': ['hawk 2090']}, \
    'ACS_F435W': {'file': 'HST-ACS_F435W.txt.gz', 'description': 'HST ACS WFC F435W', 'zeropoint': 4036.38, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'ACS_F606W': {'file': 'HST-ACS_F606W.txt.gz', 'description': 'HST ACS WFC F606W', 'zeropoint': 3351.09, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'ACS_F814W': {'file': 'HST-ACS_F814W.txt.gz', 'description': 'HST ACS WFC F814W', 'zeropoint': 2440.74, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'IRAC_CH1': {'file': 'irac1.txt.gz', 'description': 'IRAC Channel 1 (3.6 micron)', 'zeropoint': 280.9, 'method': 'vega', 'rsr': True, 'altname': ['irac 1','irac 3.6','[3.6]']}, \
    'IRAC_CH2': {'file': 'irac2.txt.gz', 'description': 'IRAC Channel 2 (4.5 micron)', 'zeropoint': 179.7, 'method': 'vega', 'rsr': True, 'altname': ['irac 2','irac 4.5','[4.5]']}, \
    'IRAC_CH3': {'file': 'irac3.txt.gz', 'description': 'IRAC Channel 3 (5.8 micron)', 'zeropoint': 115.0, 'method': 'vega', 'rsr': True, 'altname': ['irac 3','irac 5.8','[5.8]']}, \
    'IRAC_CH4': {'file': 'irac4.txt.gz', 'description': 'IRAC Channel 4 (8.0 micron)', 'zeropoint': 64.13, 'method': 'vega', 'rsr': True, 'altname': ['irac 4','irac 8.0','[8.0]']}, \
    'KEPLER': {'file': 'Kepler.txt.gz', 'description': 'Kepler bandpass', 'zeropoint': 3033.1, 'method': 'vega', 'rsr': False, 'altname': ['kep','kepler k','kp']}, \
    'LSST_U': {'file': 'LSST_u.txt.gz', 'description': 'LSST u', 'zeropoint': 2038.22, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'LSST_G': {'file': 'LSST_g.txt.gz', 'description': 'LSST g', 'zeropoint': 3990.48, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'LSST_R': {'file': 'LSST_r.txt.gz', 'description': 'LSST r', 'zeropoint': 3163.86, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'LSST_I': {'file': 'LSST_i.txt.gz', 'description': 'LSST i', 'zeropoint': 2576.07, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'LSST_Z': {'file': 'LSST_z.txt.gz', 'description': 'LSST z', 'zeropoint': 2261.40, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'LSST_Y': {'file': 'LSST_y.txt.gz', 'description': 'LSST y', 'zeropoint': 2165.23, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F0560W': {'file': 'JWST-MIRI_F560W.txt.gz', 'description': 'JWST MIRI F0560W', 'zeropoint': 115.29, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F0770W': {'file': 'JWST-MIRI_F770W.txt.gz', 'description': 'JWST MIRI F0770W', 'zeropoint': 65.08, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F1000W': {'file': 'JWST-MIRI_F1000W.txt.gz', 'description': 'JWST MIRI F1000W', 'zeropoint': 38.51, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F1065C': {'file': 'JWST-MIRI_F1065C.txt.gz', 'description': 'JWST MIRI F1065C', 'zeropoint': 33.89, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F1130W': {'file': 'JWST-MIRI_F1130W.txt.gz', 'description': 'JWST MIRI F1130W', 'zeropoint': 29.63, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F1140C': {'file': 'JWST-MIRI_F1140C.txt.gz', 'description': 'JWST MIRI F1140C', 'zeropoint': 29.66, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F1280W': {'file': 'JWST-MIRI_F1280W.txt.gz', 'description': 'JWST MIRI F1280W', 'zeropoint': 23.52, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F1500W': {'file': 'JWST-MIRI_F1500W.txt.gz', 'description': 'JWST MIRI F1500W', 'zeropoint': 17.12, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F1550C': {'file': 'JWST-MIRI_F1550C.txt.gz', 'description': 'JWST MIRI F1550C', 'zeropoint': 15.93, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F1800W': {'file': 'JWST-MIRI_F1800W.txt.gz', 'description': 'JWST MIRI F1800W', 'zeropoint': 11.99, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F2100W': {'file': 'JWST-MIRI_F2100W.txt.gz', 'description': 'JWST MIRI F2100W', 'zeropoint': 9.06, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F2300C': {'file': 'JWST-MIRI_F2300C.txt.gz', 'description': 'JWST MIRI F2300C', 'zeropoint': 7.62, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MIRI_F2550W': {'file': 'JWST-MIRI_F2550W.txt.gz', 'description': 'JWST MIRI F2550W', 'zeropoint': 6.07, 'method': 'vega', 'rsr': True, 'altname': []}, \
    'MKO_J_ATM': {'file': 'j_atm_mko.txt.gz', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_H_ATM': {'file': 'h_atm_mko.txt.gz', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_K_ATM': {'file': 'k_atm_mko.txt.gz', 'description': 'MKO K-band + atmosphere', 'zeropoint': 647.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
# using WFCAM Y for MKO Y
    'MKO_Y': {'file': 'wfcam-y.txt.gz', 'description': 'MKO Y-band + atmosphere', 'zeropoint': 1562.3, 'method': 'vega', 'rsr': False, 'altname': ['Y']}, \
    'MKO_J': {'file': 'mko_j.txt.gz', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_H': {'file': 'mko_h.txt.gz', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_K': {'file': 'mko_ks.txt.gz', 'description': 'MKO K-band', 'zeropoint': 647.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_KP': {'file': 'mko_kp.txt.gz', 'description': 'MKO Kp-band', 'zeropoint': 693.7, 'method': 'vega', 'rsr': False, 'altname': ['mko k prime']}, \
    'MKO_LP': {'file': 'mko_lp.txt.gz', 'description': 'MKO Lp-band', 'zeropoint': 248.3, 'method': 'vega', 'rsr': False, 'altname': ['mko l','mko l prime']}, \
    'MKO_MP': {'file': 'mko_mp.txt.gz', 'description': 'MKO Mp-band', 'zeropoint': 164.7, 'method': 'vega', 'rsr': False, 'altname': ['mko m','mko m prime']}, \
    'NICMOS_F090M': {'file': 'nic1_f090m.txt.gz', 'description': 'NICMOS F090M', 'zeropoint': 2255.0, 'method': 'vega', 'rsr': False, 'altname': ['F090M']}, \
    'NICMOS_F095N': {'file': 'nic1_f095n.txt.gz', 'description': 'NICMOS F095N', 'zeropoint': 2044.6, 'method': 'vega', 'rsr': False, 'altname': ['F095N']}, \
    'NICMOS_F097N': {'file': 'nic1_f097n.txt.gz', 'description': 'NICMOS F097N', 'zeropoint': 2275.4, 'method': 'vega', 'rsr': False, 'altname': ['F097N']}, \
    'NICMOS_F108N': {'file': 'nic1_f108n.txt.gz', 'description': 'NICMOS F108N', 'zeropoint': 1937.3, 'method': 'vega', 'rsr': False, 'altname': ['F108N']}, \
    'NICMOS_F110M': {'file': 'nic1_f110m.txt.gz', 'description': 'NICMOS F110M', 'zeropoint': 1871.8, 'method': 'vega', 'rsr': False, 'altname': ['F110M']}, \
    'NICMOS_F110W': {'file': 'nic1_f110w.txt.gz', 'description': 'NICMOS F110W', 'zeropoint': 1768.5, 'method': 'vega', 'rsr': False, 'altname': ['']}, \
    'NICMOS_F113N': {'file': 'nic1_f113n.txt.gz', 'description': 'NICMOS F113N', 'zeropoint': 1821.0, 'method': 'vega', 'rsr': False, 'altname': ['F113N']}, \
    'NICMOS_F140W': {'file': 'nic1_f140w.txt.gz', 'description': 'NICMOS F140W', 'zeropoint': 1277.1, 'method': 'vega', 'rsr': False, 'altname': ['']}, \
    'NICMOS_F145M': {'file': 'nic1_f145m.txt.gz', 'description': 'NICMOS F145M', 'zeropoint': 1242.0, 'method': 'vega', 'rsr': False, 'altname': ['F145M']}, \
    'NICMOS_F160W': {'file': 'nic1_f160w.txt.gz', 'description': 'NICMOS F160W', 'zeropoint': 1071.7, 'method': 'vega', 'rsr': False, 'altname': ['']}, \
    'NICMOS_F164N': {'file': 'nic1_f164n.txt.gz', 'description': 'NICMOS F164N', 'zeropoint': 1003.0, 'method': 'vega', 'rsr': False, 'altname': ['']}, \
    'NICMOS_F165M': {'file': 'nic1_f165m.txt.gz', 'description': 'NICMOS F165M', 'zeropoint': 1023.6, 'method': 'vega', 'rsr': False, 'altname': ['F165M']}, \
    'NICMOS_F166N': {'file': 'nic1_f166n.txt.gz', 'description': 'NICMOS F166N', 'zeropoint': 1047.7, 'method': 'vega', 'rsr': False, 'altname': ['F166N']}, \
    'NICMOS_F170M': {'file': 'nic1_f170m.txt.gz', 'description': 'NICMOS F170M', 'zeropoint': 979.1, 'method': 'vega', 'rsr': False, 'altname': ['F170M']}, \
    'NICMOS_F187N': {'file': 'nic1_f187n.txt.gz', 'description': 'NICMOS F187N', 'zeropoint': 803.7, 'method': 'vega', 'rsr': False, 'altname': ['F187N']}, \
    'NICMOS_F190N': {'file': 'nic1_f190n.txt.gz', 'description': 'NICMOS F190N', 'zeropoint': 836.5, 'method': 'vega', 'rsr': False, 'altname': ['F190N']}, \
    'NIRC2_J': {'file': 'nirc2-j.txt.gz', 'description': 'NIRC2 J-band', 'zeropoint': 1562.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_H': {'file': 'nirc2-h.txt.gz', 'description': 'NIRC2 H-band', 'zeropoint': 1075.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_HCONT': {'file': 'nirc2-hcont.txt.gz', 'description': 'NIRC2 H-continuum band', 'zeropoint': 1044.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_K': {'file': 'nirc2-k.txt.gz', 'description': 'NIRC2 K-band', 'zeropoint': 648.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_KP': {'file': 'nirc2-kp.txt.gz', 'description': 'NIRC2 Kp-band', 'zeropoint': 689.3, 'method': 'vega', 'rsr': False, 'altname': ['nirc2 k prime']}, \
    'NIRC2_KS': {'file': 'nirc2-ks.txt.gz', 'description': 'NIRC2 Ks-band', 'zeropoint': 676.2, 'method': 'vega', 'rsr': False, 'altname': ['nirc2 k short']}, \
    'NIRC2_KCONT': {'file': 'nirc2-kcont.txt.gz', 'description': 'NIRC2 K continuum-band', 'zeropoint': 605.9, 'method': 'vega', 'rsr': False, 'altname': ['nirc2 k continuum']}, \
    'NIRC2_FE2': {'file': 'nirc2-fe2.txt.gz', 'description': 'NIRC2 Fe II', 'zeropoint': 1019.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_LP': {'file': 'nirc2-lp.txt.gz', 'description': 'NIRC2 LP', 'zeropoint': 248.0, 'method': 'vega', 'rsr': False, 'altname': ['nirc2 l prime','nirc2 l']}, \
    'NIRC2_M': {'file': 'nirc2-ms.txt.gz', 'description': 'NIRC2 M', 'zeropoint': 165.8, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRCAM_F070W': {'file': 'JWST-NIRCAM_F070W.txt.gz', 'description': 'JWST NIRCAM F070W (wide 0.70 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F090W': {'file': 'JWST-NIRCAM_F090W.txt.gz', 'description': 'JWST NIRCAM F090W (wide 0.90 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F115W': {'file': 'JWST-NIRCAM_F115W.txt.gz', 'description': 'JWST NIRCAM F115W (wide 1.15 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F140M': {'file': 'JWST-NIRCAM_F140M.txt.gz', 'description': 'JWST NIRCAM F140M (medium 1.40 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F150W': {'file': 'JWST-NIRCAM_F150W.txt.gz', 'description': 'JWST NIRCAM F150W (wide 1.50 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F150W2': {'file': 'JWST-NIRCAM_F150W2.txt.gz', 'description': 'JWST NIRCAM F150W2 (wide 1.50 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F162M': {'file': 'JWST-NIRCAM_F162M.txt.gz', 'description': 'JWST NIRCAM F162M (medium 1.62 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F164N': {'file': 'JWST-NIRCAM_F164N.txt.gz', 'description': 'JWST NIRCAM F164N (narrow 1.64 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F182M': {'file': 'JWST-NIRCAM_F182M.txt.gz', 'description': 'JWST NIRCAM F182M (medium 1.82 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F187N': {'file': 'JWST-NIRCAM_F187N.txt.gz', 'description': 'JWST NIRCAM F187N (narrow 1.87 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F200W': {'file': 'JWST-NIRCAM_F200W.txt.gz', 'description': 'JWST NIRCAM F200W (wide 2.00 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F210M': {'file': 'JWST-NIRCAM_F210M.txt.gz', 'description': 'JWST NIRCAM F210M (medium 2.10 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F212N': {'file': 'JWST-NIRCAM_F212N.txt.gz', 'description': 'JWST NIRCAM F212N (narrow 2.12 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F250M': {'file': 'JWST-NIRCAM_F250M.txt.gz', 'description': 'JWST NIRCAM F250M (medium 2.50 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F277W': {'file': 'JWST-NIRCAM_F277W.txt.gz', 'description': 'JWST NIRCAM F277W (wide 2.77 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F300M': {'file': 'JWST-NIRCAM_F300M.txt.gz', 'description': 'JWST NIRCAM F300M (medium 3.00 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F322W2': {'file': 'JWST-NIRCAM_F322W2.txt.gz', 'description': 'JWST NIRCAM F322W2 (wide 3.22 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F323N': {'file': 'JWST-NIRCAM_F323N.txt.gz', 'description': 'JWST NIRCAM F323N (narrow 3.23 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F335M': {'file': 'JWST-NIRCAM_F335M.txt.gz', 'description': 'JWST NIRCAM F335M (medium 3.35 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F356W': {'file': 'JWST-NIRCAM_F356W.txt.gz', 'description': 'JWST NIRCAM F356W (wide 3.56 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F360M': {'file': 'JWST-NIRCAM_F360M.txt.gz', 'description': 'JWST NIRCAM F360M (medium 3.60 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F405N': {'file': 'JWST-NIRCAM_F405N.txt.gz', 'description': 'JWST NIRCAM F405N (narrow 4.05 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F410M': {'file': 'JWST-NIRCAM_F410M.txt.gz', 'description': 'JWST NIRCAM F410M (medium 4.10 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F430M': {'file': 'JWST-NIRCAM_F430M.txt.gz', 'description': 'JWST NIRCAM F430M (medium 4.30 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F444W': {'file': 'JWST-NIRCAM_F444W.txt.gz', 'description': 'JWST NIRCAM F444W (wide 4.44 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': ['F444W']}, \
    'NIRCAM_F460M': {'file': 'JWST-NIRCAM_F460M.txt.gz', 'description': 'JWST NIRCAM F460M (medium 4.60 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F466N': {'file': 'JWST-NIRCAM_F466N.txt.gz', 'description': 'JWST NIRCAM F466N (narrow 4.66 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F470N': {'file': 'JWST-NIRCAM_F470N.txt.gz', 'description': 'JWST NIRCAM F470N (narrow 4.70 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F480M': {'file': 'JWST-NIRCAM_F480M.txt.gz', 'description': 'JWST NIRCAM F480M (medium 4.80 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'PANSTARRS_G': {'file': 'panstarrs-g.txt.gz', 'description': 'PANSTARRS g-band', 'zeropoint': 3909.11, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_R': {'file': 'panstarrs-r.txt.gz', 'description': 'PANSTARRS r-band', 'zeropoint': 3151.44, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_W': {'file': 'panstarrs-w.txt.gz', 'description': 'PANSTARRS w-band', 'zeropoint': 3024.76, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_I': {'file': 'panstarrs-i.txt.gz', 'description': 'PANSTARRS i-band', 'zeropoint': 2584.6, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_Z': {'file': 'panstarrs-z.txt.gz', 'description': 'PANSTARRS z-band', 'zeropoint': 2273.09, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_Y': {'file': 'panstarrs-y.txt.gz', 'description': 'PANSTARRS y-band', 'zeropoint': 2205.95, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'ROMAN_F062': {'file': 'Roman-WFI_F062.txt.gz', 'description': 'Roman WFI F062', 'zeropoint': 3174.18, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'ROMAN_F087': {'file': 'Roman-WFI_F087.txt.gz', 'description': 'Roman WFI F087', 'zeropoint': 2294.04, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'ROMAN_F106': {'file': 'Roman-WFI_F106.txt.gz', 'description': 'Roman WFI F106', 'zeropoint': 1967.08, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'ROMAN_F129': {'file': 'Roman-WFI_F129.txt.gz', 'description': 'Roman WFI F129', 'zeropoint': 1483.60, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'ROMAN_F146': {'file': 'Roman-WFI_F146.txt.gz', 'description': 'Roman WFI F146', 'zeropoint': 1683.10, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'ROMAN_F158': {'file': 'Roman-WFI_F158.txt.gz', 'description': 'Roman WFI F158', 'zeropoint': 1355.34, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'ROMAN_F184': {'file': 'Roman-WFI_F184.txt.gz', 'description': 'Roman WFI F184', 'zeropoint': 1396.74, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'ROMAN_F213': {'file': 'Roman-WFI_F213.txt.gz', 'description': 'Roman WFI F213', 'zeropoint': 1091.96, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'ROMAN_GRISM': {'file': 'Roman-WFI_Grism.txt.gz', 'description': 'Roman WFI Grism passband', 'zeropoint': 854.12, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'ROMAN_PRISM': {'file': 'Roman-WFI_Prism.txt.gz', 'description': 'Roman WFI Prism passband', 'zeropoint': 675.72, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'SDSS_U': {'file': 'sdss-u.txt.gz', 'description': 'SDSS u-band', 'zeropoint': 1568.5, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'SDSS_G': {'file': 'sdss-g.txt.gz', 'description': 'SDSS g-band', 'zeropoint': 3965.9, 'method': 'ab', 'rsr': False, 'altname': ['g']}, \
    'SDSS_R': {'file': 'sdss-r.txt.gz', 'description': 'SDSS r-band', 'zeropoint': 3162.0, 'method': 'ab', 'rsr': False, 'altname': ['r']}, \
    'SDSS_I': {'file': 'sdss-i.txt.gz', 'description': 'SDSS i-band', 'zeropoint': 2602.0, 'method': 'ab', 'rsr': False, 'altname': ['i']}, \
    'SDSS_Z': {'file': 'sdss-z.txt.gz', 'description': 'SDSS z-band', 'zeropoint': 2244.7, 'method': 'ab', 'rsr': False, 'altname': ['z']}, \
    'SKYMAPPER_U': {'file': 'skymapper-u.txt.gz', 'description': 'SkyMapper u-band', 'zeropoint': 1320.1, 'method': 'ab', 'rsr': False, 'altname': ['skymapper u'], 'citation': '2011PASP..123..789B'}, \
    'SKYMAPPER_V': {'file': 'skymapper-v.txt.gz', 'description': 'SkyMapper v-band', 'zeropoint': 2771.8, 'method': 'ab', 'rsr': False, 'altname': ['skymapper v'], 'citation': '2011PASP..123..789B'}, \
    'SKYMAPPER_G': {'file': 'skymapper-g.txt.gz', 'description': 'SkyMapper g-band', 'zeropoint': 3728.2, 'method': 'ab', 'rsr': False, 'altname': ['skymapper g'], 'citation': '2011PASP..123..789B'}, \
    'SKYMAPPER_R': {'file': 'skymapper-r.txt.gz', 'description': 'SkyMapper r-band', 'zeropoint': 3186.0, 'method': 'ab', 'rsr': False, 'altname': ['skymapper r'], 'citation': '2011PASP..123..789B'}, \
    'SKYMAPPER_I': {'file': 'skymapper-i.txt.gz', 'description': 'SkyMapper i-band', 'zeropoint': 2495.7, 'method': 'ab', 'rsr': False, 'altname': ['skymapper i'], 'citation': '2011PASP..123..789B'}, \
    'SKYMAPPER_Z': {'file': 'skymapper-z.txt.gz', 'description': 'SkyMapper z-band', 'zeropoint': 2227.6, 'method': 'ab', 'rsr': False, 'altname': ['skymapper z'], 'citation': '2011PASP..123..789B'}, \
    'SPECULOOS': {'file': 'SPECULOOS_iz.txt.gz', 'description': 'SPECULOOS iz bandpass', 'zeropoint': 2317.40, 'method': 'vega', 'rsr': False, 'altname': ['speculoos_iz','iz']}, \
    'TESS': {'file': 'TESS.txt.gz', 'description': 'TESS bandpass', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_Z': {'file': 'ukidss-z.txt.gz', 'description': 'UKIDSS Z-band', 'zeropoint': 2261.4, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_Y': {'file': 'ukidss-y.txt.gz', 'description': 'UKIDSS Y-band', 'zeropoint': 2057.2, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_J': {'file': 'ukidss-j.txt.gz', 'description': 'UKIDSS J-band', 'zeropoint': 1556.8, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_H': {'file': 'ukidss-h.txt.gz', 'description': 'UKIDSS H-band', 'zeropoint': 1038.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_K': {'file': 'ukidss-k.txt.gz', 'description': 'UKIDSS K-band', 'zeropoint': 644.1, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_Z': {'file': 'vista_z.txt.gz', 'description': 'VISTA Z-band', 'zeropoint': 2263.81, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_Y': {'file': 'vista_y.txt.gz', 'description': 'VISTA Y-band', 'zeropoint': 2087.32, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_J': {'file': 'vista_j.txt.gz', 'description': 'VISTA J-band', 'zeropoint': 1554.03, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_H': {'file': 'vista_h.txt.gz', 'description': 'VISTA H-band', 'zeropoint': 1030.40, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_KS': {'file': 'vista_ks.txt.gz', 'description': 'VISTA Ks-band', 'zeropoint': 674.83, 'method': 'vega', 'rsr': False, 'altname': []}, \
# need to add in WFC3 optical filters    
    'WFC3_F098M': {'file': 'HST-WFC3_F098M.txt.gz', 'description': 'WFC3 F098M', 'zeropoint': 2154.5, 'method': 'vega', 'rsr': False, 'altname': ['F098M']}, \
    'WFC3_F105W': {'file': 'HST-WFC3_F105W.txt.gz', 'description': 'WFC3 F105W', 'zeropoint': 1975.2, 'method': 'vega', 'rsr': False, 'altname': ['WFC3-Y','F105W']}, \
    'WFC3_F110W': {'file': 'HST-WFC3_F110W.txt.gz', 'description': 'WFC3 F110W', 'zeropoint': 1738.4, 'method': 'vega', 'rsr': False, 'altname': ['WFC3-YJ','F110W']}, \
    'WFC3_F125W': {'file': 'HST-WFC3_F125W.txt.gz', 'description': 'WFC3 F125W', 'zeropoint': 1564.3, 'method': 'vega', 'rsr': False, 'altname': ['WFC3-J','F125W']}, \
    'WFC3_F126N': {'file': 'HST-WFC3_F126N.txt.gz', 'description': 'WFC3 F126N', 'zeropoint': 1552.5, 'method': 'vega', 'rsr': False, 'altname': ['F126N']}, \
    'WFC3_F127M': {'file': 'HST-WFC3_F127M.txt.gz', 'description': 'WFC3 F127M', 'zeropoint': 1496.5, 'method': 'vega', 'rsr': False, 'altname': ['F127M']}, \
    'WFC3_F128N': {'file': 'HST-WFC3_F128N.txt.gz', 'description': 'WFC3 F128N', 'zeropoint': 1392.6, 'method': 'vega', 'rsr': False, 'altname': ['F128N']}, \
    'WFC3_F130N': {'file': 'HST-WFC3_F130N.txt.gz', 'description': 'WFC3 F130N', 'zeropoint': 1475.9, 'method': 'vega', 'rsr': False, 'altname': ['F130N']}, \
    'WFC3_F132N': {'file': 'HST-WFC3_F132N.txt.gz', 'description': 'WFC3 F132N', 'zeropoint': 1466.6, 'method': 'vega', 'rsr': False, 'altname': ['F132N']}, \
    'WFC3_F139M': {'file': 'HST-WFC3_F139M.txt.gz', 'description': 'WFC3 F139M', 'zeropoint': 1342.8, 'method': 'vega', 'rsr': False, 'altname': ['F139M']}, \
    'WFC3_F140W': {'file': 'HST-WFC3_F140W.txt.gz', 'description': 'WFC3 F140W', 'zeropoint': 1324.8, 'method': 'vega', 'rsr': False, 'altname': ['F140W']}, \
    'WFC3_F153M': {'file': 'HST-WFC3_F153M.txt.gz', 'description': 'WFC3 F153M', 'zeropoint': 1142.0, 'method': 'vega', 'rsr': False, 'altname': ['F153M']}, \
    'WFC3_F160W': {'file': 'HST-WFC3_F160W.txt.gz', 'description': 'WFC3 F160W', 'zeropoint': 1138.1, 'method': 'vega', 'rsr': False, 'altname': ['WFC3-H','F160W']}, \
    'WFC3_F164N': {'file': 'HST-WFC3_F164N.txt.gz', 'description': 'WFC3 F164N', 'zeropoint': 1005.5, 'method': 'vega', 'rsr': False, 'altname': ['F164N']}, \
    'WFC3_F167N': {'file': 'HST-WFC3_F167N.txt.gz', 'description': 'WFC3 F167N', 'zeropoint': 1030.0, 'method': 'vega', 'rsr': False, 'altname': ['F167N']}, \
    'WFC3_G102': {'file': 'HST-WFC3_G102.txt.gz', 'description': 'WFC3 G102', 'zeropoint': 2074.54, 'method': 'vega', 'rsr': False, 'altname': ['G102']}, \
    'WFC3_G141': {'file': 'HST-WFC3_G141.txt.gz', 'description': 'WFC3 G141', 'zeropoint': 1355.32, 'method': 'vega', 'rsr': False, 'altname': ['G141']}, \
#    'WFC3_F127M': {'file': 'wfc3_F127M.txt.gz', 'description': 'WFC3 F127M', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
#    'WFC3_F139M': {'file': 'wfc3_F139M.txt.gz', 'description': 'WFC3 F139M', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
#    'WFC3_F164N': {'file': 'wfc3_F164N.txt.gz', 'description': 'WFC3 F164N', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
#    'WFC3_F167N': {'file': 'wfc3_F167N.txt.gz', 'description': 'WFC3 F167N', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_Z': {'file': 'wfcam-z.txt.gz', 'description': 'UKIRT WFCAM Z', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_Y': {'file': 'wfcam-y.txt.gz', 'description': 'UKIRT WFCAM Y', 'zeropoint': 2040.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_J': {'file': 'wfcam-j.txt.gz', 'description': 'UKIRT WFCAM J', 'zeropoint': 1548.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_H': {'file': 'wfcam-h.txt.gz', 'description': 'UKIRT WFCAM H', 'zeropoint': 1027.1, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_H2': {'file': 'wfcam-h2.txt.gz', 'description': 'UKIRT WFCAM H2', 'zeropoint': 677.1, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_BRG': {'file': 'wfcam-brg.txt.gz', 'description': 'UKIRT WFCAM Brackett Gamma', 'zeropoint': 645.5, 'method': 'vega', 'rsr': False, 'altname': ['wfcam brackett gamma']}, \
    'WFCAM_K': {'file': 'wfcam-k.txt.gz', 'description': 'UKIRT WFCAM K', 'zeropoint': 630.0, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_J': {'file': 'wirc_jcont.txt.gz', 'description': 'WIRC J-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_H': {'file': 'wirc_hcont.txt.gz', 'description': 'WIRC H-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_K': {'file': 'wirc_kcont.txt.gz', 'description': 'WIRC K-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_CO': {'file': 'wirc_co.txt.gz', 'description': 'WIRC CO', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_CH4S': {'file': 'wirc_ch4s.txt.gz', 'description': 'WIRC CH4S', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_CH4L': {'file': 'wirc_ch4l.txt.gz', 'description': 'WIRC CH4L', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_FE2': {'file': 'wirc_feii.txt.gz', 'description': 'WIRC Fe II', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_BRGAMMA': {'file': 'wirc_brgamma.txt.gz', 'description': 'WIRC H I Brackett Gamma', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': ['wirc brackett gamma']}, \
    'WIRC_PABETA': {'file': 'wirc_pabeta.txt.gz', 'description': 'WIRC H I Paschen Beta', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': ['wirc paschen beta']}, \
    'WIRCAM_Y': {'file': 'wircam-cfht-y.txt.gz', 'description': 'CFHT WIRCAM Y', 'zeropoint': 2073.32, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_J': {'file': 'wircam-cfht-j.txt.gz', 'description': 'CFHT WIRCAM J', 'zeropoint': 1551.01, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_H': {'file': 'wircam-cfht-h.txt.gz', 'description': 'CFHT WIRCAM H', 'zeropoint': 1044.35, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_KS': {'file': 'wircam-cfht-ks.txt.gz', 'description': 'CFHT WIRCAM Ks', 'zeropoint': 674.62, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_KCONT': {'file': 'wircam-cfht-kcont.txt.gz', 'description': 'CFHT WIRCAM K-cont', 'zeropoint': 636.17, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_CH4_OFF': {'file': 'wircam-cfht-ch4s.txt.gz', 'description': 'CFHT WIRCAM CH4-off', 'zeropoint': 987.39, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_CH4_ON': {'file': 'wircam-cfht-ch4l.txt.gz', 'description': 'CFHT WIRCAM CH4-on', 'zeropoint': 1076.31, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WISE_W1': {'file': 'wise_w1.txt.gz', 'description': 'WISE W1 (3.5 micron)', 'zeropoint': 309.54, 'method': 'vega', 'rsr': True, 'altname': ['W1']}, \
    'WISE_W2': {'file': 'wise_w2.txt.gz', 'description': 'WISE W2 (4.6 micron)', 'zeropoint': 171.79, 'method': 'vega', 'rsr': True, 'altname': ['W2']}, \
    'WISE_W3': {'file': 'wise_w3.txt.gz', 'description': 'WISE W3 (13 micron)', 'zeropoint': 31.67, 'method': 'vega', 'rsr': True, 'altname': ['W3']}, \
    'WISE_W4': {'file': 'wise_w4.txt.gz', 'description': 'WISE W4 (22 micron)', 'zeropoint': 8.363, 'method': 'vega', 'rsr': True, 'altname': ['W4']} \
    }
VEGAFILE = os.path.join(FILTER_FOLDER,'vega_kurucz.txt.gz')

# instrument defaults
INSTRUMENT_DEFINITION_FILE = 'instrument.txt'
INSTRUMENT_DEFAULT_PARAMETERS = {
    'instrument_name': {'altname': ['name','instrument','inst'], 'type': str, 'default': 'UNKNOWN'},
    'altname': {'altname': [], 'type': list, 'default': []},
    'instrument_bibcode': {'altname': ['bib','cite','citation'], 'type': str, 'default': ''},
    'resolution': {'altname': ['res'], 'type': float, 'default': numpy.nan},
    'wunit': {'altname': ['waveunit','wavelength_unit'], 'type': u.quantity.Quantity, 'default': u.micron},
    'funit': {'altname': ['fluxunit','flux_unit'], 'type': u.quantity.Quantity, 'default': u.erg/u.s/u.cm/u.cm/u.micron}
    }
INSTRUMENTS = {
#	'SPEX': {'instrument_name': 'SpeX prism', 'pixelscale': 0.15*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': ['']},
    'UNKNOWN': {'instrument_name': 'UNKNOWN', 'pixelscale': 1.*u.arcsec, 'slitwidth': 1.*u.arcsec, 'altname': ['UNK'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'reader': ''},
    'RAW': {'instrument_name': 'RAW', 'pixelscale': 1.*u.arcsec, 'slitwidth': 1.*u.arcsec, 'altname': [], 'wave_unit': u.micron, 'reader': ''},
    'SED': {'instrument_name': 'SED', 'pixelscale': 1.*u.arcsec, 'wave_range': [0.1,100]*u.micron, 'slitwidth': 1.*u.arcsec, 'altname': ['SPECTRAL_ENERGY_DISTRIBUTION'], 'resolution': 100, 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'reader': ''},
    'APOGEE': {'instrument_name': 'SDSS APOGEE', 'pixelscale': 2./3.*u.arcsec, 'wave_range': [1.51,1.70]*u.micron, 'slitwidth': 2.*u.arcsec, 'resolution': 22500, 'norders': 1, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['APO'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2010SPIE.7735E..1CW', 'reader': '_readAPOGEE'},
    'BOSS': {'instrument_name': 'SDSS BOSS', 'pixelscale': 2./3.*u.arcsec, 'wave_range': [3700,10400]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 2000, 'norders': 1, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['SDSS','BOSS','EBOSS'], 'wave_unit': u.Angstrom, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'instrument_bibcode': '2013AJ....146...32S', 'reader': '_readBOSS'},
    'DEIMOS': {'instrument_name': 'Keck DEIMOS', 'pixelscale': 0.1185*u.arcsec, 'wave_range': [5000,9700]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'disperser': '600 l/mm', 'resolution': 2000, 'norders': 1, 'readnoise': 2.5, 'darkcurrent': 0., 'gain': 1.2, 'altname': ['DEIMOS'], 'instrument_bibcode': '2003spie.4841.1657F', 'wave_unit': u.Angstrom, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'reader': ''},
    'FIRE': {'instrument_name': 'Magellan FIRE', 'pixelscale': 0.18*u.arcsec, 'wave_range': [0.82,2.51]*u.micron, 'slitwidth': 0.6*u.arcsec, 'resolution': 6000, 'norders': 21, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['FIRE'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2013PASP..125..270S', 'reader': '_readFIRE'},
    'IRS-SL': {'instrument_name': 'Spitzer IRS Short-Low', 'pixelscale': 1.8*u.arcsec, 'wave_range': [5.2,14.5]*u.micron, 'slitwidth': 1.8*u.arcsec, 'resolution': 100, 'norders': 1, 'readnoise': 30., 'darkcurrent': 10., 'gain': 4.6, 'altname': ['IRS','IRS Short Low'], 'wave_unit': u.micron, 'flux_unit': u.Jy, 'instrument_bibcode': '2004ApJS..154...18H', 'reader': ''},
    'IRS-LL': {'instrument_name': 'Spitzer IRS Short-Low', 'pixelscale': 5.1*u.arcsec, 'wave_range': [14,38.]*u.micron, 'slitwidth': 5.1*u.arcsec, 'resolution': 90, 'norders': 1, 'readnoise': 30., 'darkcurrent': 40., 'gain': 4.6, 'altname': ['IRS Long Low'], 'wave_unit': u.micron, 'flux_unit': u.Jy, 'instrument_bibcode': '2004ApJS..154...18H', 'reader': ''},
    'JWST-NIRSPEC-PRISM': {'instrument_name': 'JWST NIRSpec (prism mode)', 'pixelscale': 0.43*u.arcsec, 'disperser': 'prism', 'wave_range': [0.6,6.0]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 100, 'norders': 1, 'readnoise': 3.8, 'darkcurrent': 0., 'gain': 1.9, 'altname': ['JWST-NIRSPEC'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'reader': ''},
    'JWST-NIRCAM-GRISM': {'instrument_name': 'JWST NIRCam (grism mode)', 'pixelscale': 0.43*u.arcsec, 'disperser': 'prism', 'wave_range': [2.4,5.2]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 100, 'norders': 1, 'readnoise': 3.8, 'darkcurrent': 0., 'gain': 1.9, 'altname': ['JWST-NIRCAM'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'reader': ''},
    'KAST-RED': {'instrument_name': 'Lick KAST Red Channel', 'pixelscale': 0.43*u.arcsec, 'disperser': '600/7500', 'wave_range': [5700,9200]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 1200, 'norders': 1, 'readnoise': 3.8, 'darkcurrent': 0., 'gain': 1.9, 'altname': ['KAST','KAST-R'], 'wave_unit': u.Angstrom, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'reader': '_readKAST'},
    'KAST-BLUE': {'instrument_name': 'Lick KAST Blue Channel', 'pixelscale': 0.43*u.arcsec, 'disperser': '600/4310', 'wave_range': [3300,5520]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 950, 'norders': 1, 'readnoise': 3.7, 'darkcurrent': 0., 'gain': 1.2, 'altname': ['KAST-B'], 'wave_unit': u.Angstrom, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'reader': '_readKAST'},
    'LDSS-3': {'instrument_name': 'Magellan LDSS-3', 'pixelscale': 0.189*u.arcsec, 'disperser': 'VPH-RED', 'wave_range': [6000,10000]*u.Angstrom, 'slitwidth': 0.75*u.arcsec, 'resolution': 1810, 'norders': 1, 'readnoise': 4.07, 'darkcurrent': 0., 'gain': 1, 'altname': ['LDSS3'], 'wave_unit': u.Angstrom, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'instrument_bibcode': '1994PASP..106..983A', 'reader': '_readIRAF'},
    'LRIS-RED': {'altname': ['LRIS','LRISR'], 'instrument_name': 'Keck LRIS red channel longslit', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4.6, 'darkcurrent': 0., 'gain': 1.2, 'disperser': '400/8500', 'wave_range': [6300,10100]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1200, 'norders': 1, 'wave_unit': u.Angstrom, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'instrument_bibcode': '1995PASP..107..375O', 'reader': ''},
    'LRIS-BLUE': {'altname': ['LRISB'], 'instrument_name': 'Keck LRIS blue channel longslit', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4., 'darkcurrent': 0., 'gain': 1.6, 'disperser': '400/3400', 'wave_range': [1270,5740]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 530, 'norders': 1, 'wave_unit': u.Angstrom, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'instrument_bibcode': '1995PASP..107..375O', 'reader': ''},
    'MAGE': {'instrument_name': 'Magellan MAGE', 'pixelscale': 0.3*u.arcsec, 'wave_range': [3100,10000]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 4100, 'norders': 13, 'readnoise': 2.9, 'darkcurrent': 1.0, 'gain': 1.02, 'altname': ['MagE'], 'wave_unit': u.Angstrom, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'instrument_bibcode': '2008SPIE.7014E..54M', 'reader': '_readMAGE'},
    'NIRI-J': {'instrument_name': 'Gemini NIRI J-band (G5202)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [1.05,1.41]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 610, 'norders': 1, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altname': ['NIRI','NIRI G5202'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2003PASP..115.1388H', 'reader': ''},
    'NIRI-H': {'instrument_name': 'Gemini NIRI J-band (G5203)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [1.43,1.96]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 825, 'norders': 1, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altname': ['NIRI G5203'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2003PASP..115.1388H', 'reader': ''},
    'NIRI-K': {'instrument_name': 'Gemini NIRI J-band (G5204)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [1.90,2.49]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 780, 'norders': 1, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altname': ['NIRI G5204'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2003PASP..115.1388H', 'reader': ''},
    'NIRI-L': {'instrument_name': 'Gemini NIRI J-band (G5205)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [2.99,4.15]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 690, 'norders': 1, 'readnoise': 50., 'darkcurrent': 0.25, 'gain': 12.3, 'altname': ['NIRI G5205'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2003PASP..115.1388H', 'reader': ''},
    'NIRES': {'instrument_name': 'Keck NIRES', 'pixelscale': 0.123*u.arcsec, 'wave_range': [0.94,2.45]*u.micron, 'slitwidth': 0.55*u.arcsec, 'resolution': 2700, 'norders': 5, 'orders': [3,4,5,6,7], 'order_wave_range': [[0.94,1.06],[0.95,1.23],[1.13,1.48],[1.42,1.85],[1.88,2.46]],'readnoise': 15., 'darkcurrent': 0.13, 'gain': 3.8, 'altname': [], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2000SPIE.4008.1048M', 'reader': ''},
    'NIRSPEC': {'instrument_name': 'Keck NIRSPEC', 'pixelscale': 0.43/3.*u.arcsec, 'wave_range': [0.95,5.5]*u.micron, 'slitwidth': 0.43*u.arcsec, 'resolution': 25000, 'norders': 8, 'readnoise': 23., 'darkcurrent': 0.8, 'gain': 5.8, 'altname': ['NIRSPAO'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2000SPIE.4008.1048M', 'reader': ''},
    'SPEX-PRISM': {'instrument_name': 'IRTF SpeX prism', 'pixelscale': 0.15*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': ['SPEX','PRISM'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2003PASP..115..362R', 'reader': ''},
    'SPEX-SXD': {'instrument_name': 'IRTF SpeX SXD', 'pixelscale': 0.15*u.arcsec, 'wave_range': [0.8,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': ['SXD'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2003PASP..115..362R', 'reader': ''},
    'SPEX-LXD1.9': {'instrument_name': 'IRTF SpeX LXD 1.9 micron', 'pixelscale': 0.15*u.arcsec, 'wave_range': [1.95,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': ['SPEX LXD','LXD'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2003PASP..115..362R', 'reader': ''},
    'SPEX-LXD2.1': {'instrument_name': 'IRTF SpeX LXD 2.1 micron', 'pixelscale': 0.15*u.arcsec, 'wave_range': [2.15,5.0]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': [], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2003PASP..115..362R', 'reader': ''},
    'SPEX-LXD2.3': {'instrument_name': 'IRTF SpeX LXD 2.3 micron', 'pixelscale': 0.15*u.arcsec, 'wave_range': [2.25,5.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': [], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'instrument_bibcode': '2003PASP..115..362R', 'reader': ''},
#   'USPEX': {'instrument_name': 'Updated SpeX prism', 'pixelscale': 0.10*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['']},
#   'USPEX_PRISM': {'instrument_name': 'IRTF Updated SpeX prism', 'pixelscale': 0.10*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['USPEX','UPRISM'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron},
#    'USPEX_SXD': {'instrument_name': 'IRTF Updated SpeX SXD', 'pixelscale': 0.10*u.arcsec, 'wave_range': [0.7,2.55]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['USXD'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron},
#    'USPEX_LXD_SHORT': {'instrument_name': 'IRTF Updated SpeX LXD short', 'pixelscale': 0.10*u.arcsec, 'wave_range': [1.67,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500., 'norders': 8, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['ULXD','LXD_SHORT','LXDS'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron},
#    'USPEX_LXD_LONG': {'instrument_name': 'IRTF Updated SpeX LXD long', 'pixelscale': 0.10*u.arcsec, 'wave_range': [1.98,5.3]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500., 'norders': 7, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['LXD_LONG','LXDL'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron},
    'WFC3-G102': {'instrument_name': 'HST WFC3 IR G102', 'pixelscale': 0.128*u.arcsec, 'wave_range': [0.8,1.15]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 210., 'norders': 7, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altname': ['G102'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'reader': '_readWFC3'},
    'WFC3-G141': {'instrument_name': 'HST WFC3 IR G141', 'pixelscale': 0.128*u.arcsec, 'wave_range': [1.075,1.70]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 130., 'norders': 7, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altname': ['WFC3','HST WFC3','HST WFC3 IR','WFC3 IR','G141'], 'wave_unit': u.micron, 'flux_unit': u.erg/u.s/u.cm/u.cm/u.micron, 'reader': '_readWFC3'},
}

#experimental
INSTRUMENTS_ALT = { 
    'SED': {'instrument_name': 'SED', 'pixelscale': 0.*u.arcsec, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['SPECTRAL_ENERGY_DISTRIBUTION'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'DEFAULT': {'wave_range': [0.1,100.]*u.micron, 'slitwidth': 0.*u.arcsec, 'resolution': 100, 'norders': 1, 'default': True}, 
    }}, 
    'APOGEE': {'instrument_name': 'SDSS APOGEE', 'pixelscale': 2./3.*u.arcsec, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['APO'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'DEFAULT': {'wave_range': [1.51,1.70]*u.micron, 'slitwidth': 2.*u.arcsec, 'resolution': 22500, 'norders': 1, 'default': True}, 
    }}, 
    'BOSS': {'instrument_name': 'SDSS BOSS', 'pixelscale': 2./3.*u.arcsec, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['BOSS','EBOSS'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': { 
        'DEFAULT': {'wave_range': [3700,10400]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': True}, 
    }}, 
    'FIRE': {'instrument_name': 'Magellan FIRE', 'pixelscale': 0.18*u.arcsec, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['FIRE'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': {
        'ECHELLE': {'wave_range': [0.82,2.51]*u.micron, 'slitwidth': 0.6*u.arcsec, 'resolution': 6000, 'norders': 21, 'default': True},
        'PRISM': {'wave_range': [0.82,2.51]*u.micron, 'slitwidth': 0.6*u.arcsec, 'resolution': 450, 'norders': 1, 'default': False},
    }},
    'IRS': {'instrument_name': 'Spitzer IRS', 'readnoise': 30., 'gain': 4.6, 'altname': ['IRS'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': {
        'SL': {'wave_range': [5.2,14.5]*u.micron, 'slitwidth': 1.8*u.arcsec, 'pixelscale': 1.8*u.arcsec, 'darkcurrent': 10., 'resolution': 100, 'norders': 1, 'default': True},
        'SH': {'wave_range': [9.9,19.6]*u.micron, 'slitwidth': 2.3*u.arcsec, 'pixelscale': 2.3*u.arcsec, 'darkcurrent': 10., 'resolution': 600, 'norders': 1, 'default': False},
        'LL': {'wave_range': [14,38.]*u.micron, 'slitwidth': 5.1*u.arcsec, 'pixelscale': 5.1*u.arcsec, 'darkcurrent': 40., 'resolution': 90, 'norders': 1, 'default': False},
        'LH': {'wave_range': [18.7,37.2]*u.micron, 'slitwidth': 4.5*u.arcsec, 'pixelscale': 4.5*u.arcsec, 'darkcurrent': 40., 'resolution': 90, 'norders': 1, 'default': False},
    }}, 
    'LDSS-3': {'instrument_name': 'Magellan LDSS-3', 'pixelscale': 0.189*u.arcsec, 'readnoise': 4.07, 'darkcurrent': 0., 'gain': 1, 'altname': ['LDSS3'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': {
        'VPH-RED': {'wave_range': [6000,10000]*u.Angstrom, 'slitwidth': 0.75*u.arcsec, 'resolution': 1810, 'norders': 1, 'default': True}, 
        'VPH-BLUE': {'wave_range': [3800,6200]*u.Angstrom, 'slitwidth': 0.75*u.arcsec, 'resolution': 1900, 'norders': 1, 'default': False}, 
        'VPH-ALL': {'wave_range': [4250,10000]*u.Angstrom, 'slitwidth': 0.75*u.arcsec, 'resolution': 860, 'norders': 1, 'default': False}, 
    }}, 
    'LRIS-BLUE': {'instrument_name': 'Keck LRIS blue channel longslit', 'altname': ['LRISB'], 'instrument_reference': '1998SPIE.3355...81M', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4., 'darkcurrent': 0., 'gain': 1.6, 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': { 
        '300/5000': {'default': False, 'wave_range': [1570,7420]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 570, 'norders': 1}, 
        '400/3400': {'default': True, 'wave_range': [1270,5740]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 500, 'norders': 1}, 
        '600/4000': {'default': False, 'wave_range': [3010,5600]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1000, 'norders': 1}, 
        '1200/3400': {'default': False, 'wave_range': [2910,3890]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 2200, 'norders': 1}, 
    }}, 
    'LRIS-RED': {'instrument_name': 'Keck LRIS red channel longslit', 'altname': ['LRIS','LRISR'], 'instrument_reference': '1995PASP..107..375O', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4.6, 'darkcurrent': 0., 'gain': 1.2, 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': { 
        '300/5000': {'default': False, 'wave_range': [3700,8270]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 540, 'norders': 1}, 
        '400/8500': {'default': True, 'wave_range': [6200,10900]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1200, 'norders': 1}, 
        '600/5000': {'default': False, 'wave_range': [3700,6640]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1050, 'norders': 1}, 
        '600/7500': {'default': False, 'wave_range': [6200,9140]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1600, 'norders': 1}, 
        '600/10000': {'default': False, 'wave_range': [8700,10900]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 2130, 'norders': 1}, 
    }}, 
    'MAGE': {'instrument_name': 'Magellan MAGE', 'pixelscale': 0.3*u.arcsec, 'readnoise': 2.9, 'darkcurrent': 1.0, 'gain': 1.02, 'altname': ['MagE'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': { 
        'DEFAULT': {'wave_range': [3100,10000]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 4100, 'norders': 13, 'default': True}, 
    }}, 
    'NIRI': {'instrument_name': 'Gemini NIRI', 'pixelscale': 0.47/4.*u.arcsec, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altname': ['GEMINI NIRI'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': {
        'J': {'default': True, 'wave_range': [1.05,1.41]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 610, 'norders': 1},
        'H': {'default': False, 'wave_range': [1.43,1.96]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 825, 'norders': 1},
        'K': {'default': False, 'wave_range': [1.90,2.49]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 780, 'norders': 1},
        'L': {'default': False, 'wave_range': [2.99,4.15]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 690, 'norders': 1},
        'M': {'default': False, 'wave_range': [4.45,5.45]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 770, 'norders': 1},
    }}, 
    'NIRSPEC': {'instrument_name': 'Keck NIRSPEC', 'pixelscale': 0.43/3.*u.arcsec, 'readnoise': 23., 'darkcurrent': 0.8, 'gain': 5.8, 'altname': ['NIRSPAO'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'LOWRES_N1': {'wave_range': [0.947,1.121]*u.micron, 'slitwidth': 0.38*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': False}, 
        'LOWRES_N2': {'wave_range': [1.089,1.293]*u.micron, 'slitwidth': 0.38*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': False}, 
        'LOWRES_N3': {'wave_range': [1.143,1.375]*u.micron, 'slitwidth': 0.38*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': False}, 
        'LOWRES_N4': {'wave_range': [1.241,1.593]*u.micron, 'slitwidth': 0.38*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': False}, 
        'LOWRES_N5': {'wave_range': [1.431,1.808]*u.micron, 'slitwidth': 0.38*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': False}, 
        'LOWRES_N6A': {'wave_range': [1.558,2.000]*u.micron, 'slitwidth': 0.38*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': False}, 
        'LOWRES_N6B': {'wave_range': [1.937,2.315]*u.micron, 'slitwidth': 0.38*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': False}, 
        'LOWRES_N7': {'wave_range': [1.997,2.428]*u.micron, 'slitwidth': 0.38*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': False}, \
        'HIGHRES_N3': {'wave_range': [1.18,1.30]*u.micron, 'slitwidth': 0.43*u.arcsec, 'resolution': 25000, 'norders': 8, 'default': True}, 
    }}, 
    'SPEX': {'instrument_name': 'IRTF SpeX', 'pixelscale': 0.15*u.arcsec, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': ['OLD SPEX'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'PRISM': {'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'default': True}, 
        'SXD': {'wave_range': [0.8,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'default': False}, 
        'LXD1.9': {'wave_range': [1.95,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'default': False}, 
        'LXD2.1': {'wave_range': [2.15,5.0]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'default': False}, 
        'LXD2.3': {'wave_range': [2.25,5.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'default': False}, 
    }}, 
    'USPEX': {'instrument_name': 'IRTF SpeX (upgraded)', 'pixelscale': 0.10*u.arcsec, 'readnoise': 5., 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['SPEX'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'PRISM': {'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'default': True}, 
        'SXD': {'wave_range': [0.8,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'default': False}, 
        'LXD_SHORT': {'wave_range': [1.67,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 8, 'default': False}, 
        'LXD_LONG': {'wave_range': [1.98,5.3]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'default': False}, 
    }}, 
    'WFC3-UVIS': {'instrument_name': 'HST WFC3 UVIS', 'pixelscale': 0.0395*u.arcsec, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altname': ['WFC3-OPT','WFC3-VIS','WFC3-UV'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'DEFAULT': {'wave_range': [1900,4500]*u.Angstrom, 'slitwidth': 0.3*u.arcsec, 'resolution': 70, 'norders': 1, 'default': True}, 
    }}, 
    'WFC3-IR': {'instrument_name': 'HST WFC3 IR', 'pixelscale': 0.128*u.arcsec, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altname': ['WFC3'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'G102': {'wave_range': [0.8,1.15]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 210, 'norders': 1, 'default': False}, 
        'G141': {'wave_range': [1.075,1.70]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 130, 'norders': 1, 'default': True}, 
    }}, 
    }

# Spectral model defaults
SPECTRAL_MODEL_PARAMETERS_INORDER = ['teff','logg','z','kzz','fsed','co','cld','y','enrich','zc','zo','zn','broad','logpmin','logpmax','ad','radius']
SPECTRAL_MODEL_PARAMETERS = {\
    'teff': {'name': 'temperature', 'prefix': 't', 'unit': u.K, 'default': 1000.0, 'title': r'$T_{eff}$ (K)', 'type': 'continuous'}, \
    'logg': {'name': 'gravity', 'prefix': 'g', 'unit': u.dex, 'default': 5.0, 'title': r'$\log{g}$ (cgs)', 'type': 'continuous'}, \
    'z': {'name': 'metallicity', 'prefix': 'z', 'unit': u.dex, 'default': 0., 'title': '$Z$', 'type': 'continuous'}, \
    'fsed': {'name': 'rainout', 'prefix': 'f', 'unit': u.m/u.m, 'default': 'nc', 'title': r'$f_{sed}$', 'type': 'discrete'}, \
    'cld': {'name': 'cloud', 'prefix': 'cld', 'unit': u.m/u.m, 'default': 'nc', 'title': 'Cloud or Condensation Treatment', 'type': 'discrete'}, \
    'kzz': {'name': 'mixing', 'prefix': 'k', 'unit': u.m/u.m, 'default': 'eq', 'title': r'$log\ \kappa_{zz}$ (cgs)', 'type': 'continuous'},\
    'ad': {'name': 'adiabat', 'prefix': 'ad', 'unit': u.m/u.m, 'default': 1., 'title': 'Adiabatic Index', 'type': 'continuous'},\
    'y': {'name': 'He abundance', 'prefix': 'y', 'unit': u.dex, 'default': 0.27, 'title': '$Y$', 'type': 'continuous'}, \
    'enrich': {'name': 'alpha enrichment', 'prefix': 'en', 'unit': u.dex, 'default': 0., 'title': 'Alpha Element Enrichment', 'type': 'continuous'},\
    'zc': {'name': 'C enrichment', 'prefix': 'ca', 'unit': u.dex, 'default': 0., 'title': 'Carbon Enrichment', 'type': 'continuous'},\
    'zo': {'name': 'O enrichment', 'prefix': 'ox', 'unit': u.dex, 'default': 0., 'title': 'Oxygen Enrichment', 'type': 'continuous'},\
    'zn': {'name': 'N enrichment', 'prefix': 'ni', 'unit': u.dex, 'default': 0., 'title': 'Nitrogen Enrichment', 'type': 'continuous'},\
    'co': {'name': 'C/O ratio', 'prefix': 'co', 'unit': u.dex, 'default': 0.54, 'title': 'C/O ratio', 'type': 'continuous'},\
    'broad': {'name': 'broadening', 'prefix': 'br', 'unit': u.m/u.m, 'default': 'A', 'title': 'Alkali Line Broadening Prescription', 'type': 'discrete'},\
    'logpmin': {'name': 'log pressure top', 'prefix': 'pt', 'unit': u.dex, 'default': -8., 'title': 'log Minimum Pressure (bar)', 'type': 'continuous'},\
    'logpmax': {'name': 'log pressure bottom', 'prefix': 'pb', 'unit': u.dex, 'default': 4., 'title': 'log Maximum Pressure (bar)', 'type': 'continuous'},\
    'radius': {'name': 'radius', 'prefix': 'rad', 'unit': u.Rsun, 'default': 0., 'title': r'Radius (R$_{\odot}$)', 'type': 'continuous'},\
}
SPECTRAL_MODELS = {\
#    'gaia': {'folder': SPECTRAL_MODEL_FOLDER+'/gaia/', 'name': 'AMES GAIA', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altname': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'rawfolder': HOME_FOLDER+'/models/phoenix/nextgen/fullres/', 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
#    'alvarado24': {'instruments': {}, 'name': 'Alvarado 2024', 'citation': 'Alvarado et al. (2024)', 'bibcode': '', 'altname': ['alvarado','alv24','sand','sandy'], 'default': {'teff': 1500., 'logg': 6.0, 'z': 0.10, 'enrich': 0.0}}, \
    'atmo20': {'instruments': {}, 'name': 'ATMO2020', 'citation': 'Phillips et al. (2020)', 'bibcode': '2020A%26A...637A..38P', 'altname': ['atmo','phillips','phi20','atmos2020','atmos20','atmo2020','atmo20'], 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.0,'kzz': 0.0,'cld': 'LC','broad': 'A','ad': 1.0,'logpmin': -8, 'logpmax': 4}}, \
    'atmo20pp': {'instruments': {}, 'name': 'ATMO2020++', 'citation': 'Meisner et al. (2023)', 'bibcode': '2023AJ....166...57M', 'altname': ['atmo++','meisner23','mei23','atmo2020++','atmo20++','atmos2020++','atmos20++'], 'default': {'teff': 1200., 'logg': 5.0, 'z': 0.0,'kzz': 4.0}}, \
    'btcond': {'instruments': {}, 'name': 'BT Cond', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['dusty-cond','bt-cond','btc'], 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.0, 'enrich': 0.0}}, \
    'btdusty16': {'instruments': {}, 'name': 'BT Dusty 2016', 'citation': 'TBD', 'bibcode': '', 'altname': ['btdusty2016','dusty16','dusty2016','dusty-bt','bt-dusty','bt-dusty2016','btdusty','bt-dusty16','btd'], 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0, 'enrich': 0.0}}, \
    'btnextgen': {'instruments': {}, 'name': 'BT NextGen', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['nextgen-bt','btnextgen','btn'], 'default': {'teff': 3000., 'logg': 5.0, 'z': 0.0, 'enrich': 0.}}, \
    'btsettl08': {'instruments': {}, 'name': 'BT Settl 2008', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['allard','allard12','allard2012','btsettl','btsettled','btsettl08','btsettl2008','BTSettl2008','bts','bts08'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'enrich': 0.}}, \
    'btsettl15': {'instruments': {}, 'name': 'BT Settl 2015', 'citation': 'Allard et al. (2015)', 'bibcode': '2015A&A...577A..42B', 'altname': ['allard15','allard2015','btsettl015','btsettl2015','BTSettl2015','bts15'],  'default': {'teff': 1500., 'logg': 5.0, 'z': 0.}}, \
    'burrows06': {'instruments': {}, 'name': 'Burrows et al. (2006)', 'citation': 'Burrows et al. (2006)', 'bibcode': '2006ApJ...640.1063B', 'altname': ['burrows','burrows2006'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'cld': 'nc'}}, \
    'cond01': {'instruments': {}, 'name': 'AMES Cond', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altname': ['cond','cond-ames','amescond'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'dusty01': {'instruments': {}, 'name': 'AMES Dusty', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altname': ['dusty','dusty-ames','amesdusty'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'dback24': {'instruments': {}, 'name': 'Sonora Diamondback', 'citation': 'Morley et al. (2024)', 'bibcode': '2024arXiv240200758M', 'altname': ['diamondback','sonora-diamondback','sonora-dback','dback24','diamondback24','morley24','mor24'], 'default': {'teff': 1200., 'logg': 5.0, 'z': 0., 'fsed': 'f2'}}, \
    'drift': {'instruments': {}, 'name': 'Drift', 'citation': 'Witte et al. (2011)', 'bibcode': '2011A&A...529A..44W', 'altname': ['witte','witte11','witte2011','helling'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
    'dusty01': {'instruments': {}, 'name': 'AMES Dusty', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altname': ['dusty','dusty-ames','amesdusty'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'gerasimov20': {'instruments': {}, 'name': 'Gerasimov et al. 2020', 'citation': 'Gerasimov et al. (2020)', 'bibcode': '2020RNAAS...4..214G', 'altname': ['phxlowz','ger20'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
    'helios': {'instruments': {}, 'name': 'Helios', 'citation': 'Malik et al. (2017)', 'bibcode': '2017AJ....153...56M', 'altname': ['hel','malik2017','malik17','mal17'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'co': 0.5}}, \
    'elfowl24': {'instruments': {}, 'name': 'Sonora Elfowl', 'citation': 'Mukherjee et al. (2024)', 'bibcode': '2024ApJ...963...73M', 'altname': ['elfowl','sonora-elfowl','elfowl24','mukherjee','mukherjee24','muk24'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'co': 1, 'kzz': 2.0}}, \
    'elfowl24-ph3': {'instruments': {}, 'name': 'Modified Sonora Elfowl', 'citation': 'Beiler et al. (2024)', 'bibcode': '2024ApJ...973...60B', 'altname': ['ph3','elfowl24-mod','beiler','beiler24','bei24'], 'default': {'teff': 500., 'logg': 5.0, 'z': 0., 'co': 1, 'kzz': 2.0}}, \
    'karalidi21': {'instruments': {}, 'name': 'Sonora Cholla', 'citation': 'Karalidi et al. (2021)', 'bibcode': '2021ApJ...923..269K', 'altname': ['karalidi2021','karalidi','sonora-cholla','cholla'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'kzz': 4.0}}, \
    'lacy23': {'instruments': {}, 'name': 'Lacy & Burrows (2023)', 'citation': 'Lacy & Burrows (2023)', 'bibcode': '2023ApJ...950....8L', 'altname': ['lacy2023','lac23','lacy'], 'default': {'teff': 400., 'logg': 4.0, 'z': 0., 'cld': 'nc', 'kzz': 0.}}, \
    'lowz': {'instruments': {}, 'name': 'LowZ models', 'citation': 'Meisner et al. (2021)', 'bibcode': '2021ApJ...915..120M', 'altname': ['meisner2021','mei21','line21','line2021'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'kzz': 2.0, 'co': 0.85}}, \
    'lowz-cold': {'instruments': {}, 'name': 'LowZ Cold Models', 'citation': 'M. Line (priv. comm.)', 'bibcode': '', 'altname': ['line24','lin24'], 'default': {'teff': 500., 'logg': 4.5, 'z': 0., 'co': 0.91}}, \
    'madhu11': {'instruments': {}, 'name': 'Madhusudhan et al. (2011)', 'citation': 'Madhusudhan et al. (2011)', 'bibcode': '2011ApJ...737...34M', 'altname': ['madhu','madhusudhan','madhusudhan11','madhu2011','madhusudhan2011'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.,'cld': 'ae60', 'kzz': 0.0,'fsed': -1.0}}, \
    'meisner23': {'instruments': {}, 'name': 'Meisner et al. (2023)', 'citation': 'Meisner et al. (2023)', 'bibcode': '2023AJ....166...57M', 'altname': ['meisner','mei23'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'kzz': 5.0, 'ad': 3.95}}, \
    'morley12': {'instruments': {}, 'name': 'Morley et al. (2012)', 'citation': 'Morley et al. (2012)', 'bibcode': '2012ApJ...756..172M', 'altname': ['morley','morley2012'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'fsed': 'f5'}}, \
    'morley14': {'instruments': {}, 'name': 'Morley et al. (2014)', 'citation': 'Morley et al. (2014)', 'bibcode': '2014ApJ...787...78M', 'altname': ['morley2014'], 'default': {'teff': 300., 'logg': 5.0, 'z': 0., 'fsed': 'f5', 'cld': 'h50'}}, \
#    'nextgen99': {'instruments': {}, 'name': 'Phoenix NextGen', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altname': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
    'phillips24': {'instruments': {}, 'name': 'Phillips et al. (2024)', 'citation': 'Phillips et al. (2024)', 'bibcode': '2024ApJ...961..210P', 'altname': ['phi24'], 'default': {'teff': 800., 'logg': 5.0, 'z': 0., 'kzz': 4.0, 'zo': 0.0}}, \
    'petrus22': {'instruments': {}, 'name': 'Petrus et al. (2022)', 'citation': 'Petrus et al. (2022)', 'bibcode': '', 'altname': ['pet22','petrus'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'kzz': 5.0, 'co': 1.0, 'ad': 1.03}}, \
    'sand24': {'instruments': {}, 'name': 'Alvarado et al. 2024', 'citation': 'Alvarado et al. 2024', 'bibcode': '2024RNAAS...8..134A', 'altname': ['sand','san24','sand2024'], 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.1, 'enrich': 0.0}}, \
    'saumon08': {'instruments': {}, 'name': 'Saumon & Marley 2008', 'citation': 'Saumon & Marley 2008', 'bibcode': '2008ApJ...689.1327S', 'altname': ['sau08','saumon2008'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
    'saumon12': {'instruments': {}, 'name': 'Saumon et al. 2012', 'citation': 'Saumon & Marley 2008', 'bibcode': '2012ApJ...750...74S', 'altname': ['saumon','sau12','saumon2012'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
    'sonora18': {'instruments': {}, 'name': 'Sonora Alpha', 'citation': 'Marley et al. (2018)', 'bibcode': 'marley_mark_2018_1309035', 'altname': ['marley','marley18','marley2018','sonora2018'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'cld': 'nc'}}, \
    'sonora21': {'instruments': {}, 'name': 'Sonora Bobcat', 'citation': 'Marley et al. (2021)', 'bibcode': '2021ApJ...920...85M', 'altname': ['marley2021','sonora','sonora2021','bobcat','sonora-bobcat'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'co': 1}}, \
#    'gerasimov23': {'instruments': {}, 'name': 'Gerasimov 2023', 'citation': 'Gerasimov et al. (2023)', 'bibcode': '', 'altname': ['gerasimov','ger23'], 'default': {'teff': 1500., 'logg': 5.0, 'z': -0.5, 'enrich': 0.30}}, \
    'tremblin15': {'instruments': {}, 'name': 'Tremblin et al. 2015', 'citation': 'Tremblin et al. 2015', 'bibcode': '2015ApJ...804L..17T', 'altname': ['tremblin','tre15','tremblin2015'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0, 'kzz': 8.0, 'ad': 1.20}}, \
#    'tremblin16': {'instruments': {}, 'name': 'Tremblin et al. 2016', 'citation': 'Tremblin et al. 2016', 'bibcode': '2016ApJ...817L..19T', 'altname': ['tremblin','tre16','tremblin2016'], 'default': {'teff': 1300., 'logg': 5.0, 'z': 0.1, 'kzz': 6.0, 'ad': 1.05}}, \
    'veyette': {'instruments': {}, 'name': 'Veyette et al. 2017', 'citation': 'Veyette et al. 2017', 'bibcode': '2017ApJ...851...26V', 'altname': ['veyette17','veyette2017'], 'default': {'teff': 3000., 'logg': 5.0, 'z': 0.0, 'enrich': 0.0, 'carbon': 0.0, 'oxygen': 0.0}}, \
}


# Evolutionary model defaults
EVOLUTIONARY_MODEL_PARAMETERS = {\
    'mass': {'unit': u.solMass, 'default': 0.05, 'title': '$M$', 'altname': ['m','masses']},\
    'age': {'unit': u.Gyr, 'default': 5., 'title': '$\tau$', 'altname': ['a','time','ages']},\
    'temperature': {'unit': u.K, 'default': 1000.0, 'title': r'$T_{eff}$', 'altname': ['t','teff','temp','temperatures']},\
    'gravity': {'unit': u.dex(u.cm / u.s**2), 'default': 5.0, 'title': r'$\log{g}', 'altname': ['g','logg','gravities','grav']},\
    'luminosity': {'unit': u.dex(u.solLum), 'default': -5., 'title': r'$\log{L_{bol}/L_{\odot}}$', 'altname': ['l','lbol','lum','luminosities']},\
    'radius': {'unit': u.solRad, 'default': 0.1, 'title': r'$R_{\odot}$', 'altname': ['r','rad','radii','radiuses']},
    'metallicity': {'unit': u.dex(), 'default': 0., 'title': '[M/H]', 'altname': ['z','metal','abundance']},
    'y': {'unit': u.dex(), 'default': 0., 'title': 'Y', 'altname': ['he']},
    'l_mix': {'unit': u.cm/u.cm, 'default': 1.0, 'title': 'Mixing Length', 'altname': ['alpha']},
    'cloud': {'unit': u.dex(), 'default': '', 'title': 'Cloud parameter', 'altname': ['cld']},
    'Mv': {'unit': u.dex(), 'default': '', 'title': 'Absolute V', 'altname': ['mv']},
}
EVOLUTIONARY_MODELS = {\
    'atmo2020-ceq': {'file': 'atmo2020-ceq.csv','citation': 'Phillips et al. (2020) CEQ', 'bibcode': '2020A&A...637A..38P', 'altname': ['phi20','phillips20','phillips2020','phi20-ceq','phillips20-ceq','phillips2020-ceq','atmo','atmo20','atmo2020'], 'default': {'metallicity': 0.}},\
    'atmo2020-neq-strong': {'file': 'atmo2020-neq-strong.csv','citation': 'Phillips et al. (2020) NEQ (strong)', 'bibcode': '2020A&A...637A..38P', 'altname': ['phi20-neq','phillips20-neq','phillips2020-neq','phi20-neq-strong','phillips20-neq-strong','phillips2020-neq-strong','atmo-neq','atmo20-neq','atmo2020-neq'], 'default': {'metallicity': 0.}},\
    'atmo2020-neq-weak': {'file': 'atmo2020-neq-weak.csv','citation': 'Phillips et al. (2020) NEQ (weak)', 'bibcode': '2020A&A...637A..38P', 'altname': ['phi20-neq-weak','phillips20-neq-weak','phillips2020-neq-weak'], 'default': {'metallicity': 0.}},\
    'atmo2020-eos-ceq': {'file': 'atmo2020-eos-ceq.csv','citation': 'Phillips et al. (2020) CEQ new EOS', 'bibcode': '2020A&A...637A..38P', 'altname': ['phi20-eos','phillips20-eos','phillips2020-eos','phi20-eos-ceq','phillips20-eos-ceq','phillips2020-eos-ceq','atmo-eos','atmo20-eos','atmo2020-eos'], 'default': {'metallicity': 0.}},\
    'atmo2020-eos-neq-strong': {'file': 'atmo2020-eos-neq-strong.csv','citation': 'Phillips et al. (2020) NEQ (strong) new EOS', 'bibcode': '2020A&A...637A..38P', 'altname': ['phi20-eos-neq','phillips20-eos-neq','phillips2020-eos-neq','phi20-eos-neq-strong','phillips20-eos-neq-strong','phillips2020-eos-neq-strong','atmo-eos-neq','atmo20-eos-neq','atmo2020-eos-neq'], 'default': {'metallicity': 0.}},\
    'atmo2020-eos-neq-weak': {'file': 'atmo2020-eos-neq-weak.csv','citation': 'Phillips et al. (2020) NEQ (weak) new EOS', 'bibcode': '2020A&A...637A..38P', 'altname': ['phi20-eos-neq-weak','phillips20-eos-neq-weak','phillips2020-eos-neq-weak'], 'default': {'metallicity': 0.}},\
    'baraffe1997': {'file': 'baraffe1997.csv','citation': 'Baraffe et al. (1997)', 'bibcode': '1997A&A...327.1054B', 'altname': ['bar97','baraffe97'], 'default': {'metallicity': 0.}},\
    'baraffe1998': {'file': 'baraffe1998.csv','citation': 'Baraffe et al. (1998)', 'bibcode': '1998A&A...337..403B', 'altname': ['bar98','baraffe98'], 'default': {'metallicity': 0., 'y': 0.275, 'l_mix': 1.}},\
    'baraffe2001-cond': {'file': 'baraffe2001-cond.csv','citation': 'Baraffe et al. (2001) COND', 'bibcode': '2001ASPC..243..571B', 'altname': ['bar01','baraffe01','baraffe2001','bar01-cond','baraffe01-cond','cond','cond01','cond2001'], 'default': {'metallicity': 0.}},\
    'baraffe2001-dusty': {'file': 'baraffe2001-dusty.csv','citation': 'Baraffe et al. (2001) DUSTY', 'bibcode': '2001ASPC..243..571B', 'altname': ['bar01-dusty','baraffe01-dusty','dusty01','dusty2001'], 'default': {'metallicity': 0.}},\
    'baraffe2003': {'file': 'baraffe2003.csv','citation': 'Baraffe et al. (2003)', 'bibcode': '2003A&A...402..701B', 'altname': ['bar03','baraffe03'], 'default': {}},\
    'baraffe2015': {'file': 'baraffe2015.csv','citation': 'Baraffe et al. (2015)', 'bibcode': '2015A&A...577A..42B', 'altname': ['baraffe','bar15','baraffe15'], 'default': {}},\
    'burrows2001': {'file': 'burrows2001.csv','citation': 'Burrows et al. (2001)', 'bibcode': '2001RvMP...73..719B', 'altname': ['burrows','bur01','burrows01'], 'default': {}},\
    'chabrier1997': {'file': 'chabrier1997.csv','citation': 'Chabrier & Baraffe (1997)', 'bibcode': '1997A&A...328...83C', 'altname': ['chabrier','cha97','chabrier97'], 'default': {'metallicity': 0., 'y': 0.275}},\
    'chabrier2000': {'file': 'chabrier2000.csv','citation': 'Chabrier et al. (2000)', 'bibcode': '2000ApJ...542..464C', 'altname': ['dusty','cha0','chabrier00','dusty00','dusty2000'], 'default': {'metallicity': 0.}},\
    'saumon2008': {'file': 'saumon2008.csv','citation': 'Saumon et al. (2008)', 'bibcode': '2008ApJ...689.1327S', 'altname': ['saumon','sau08','saumon08'], 'default': {'metallicity': 0., 'cloud': 'hybrid'}},\
    'saumon2008-nc': {'file': 'saumon2008-nc.csv','citation': 'Saumon et al. (2008) no clouds', 'bibcode': '2008ApJ...689.1327S', 'altname': ['saumon','sau08','saumon08'], 'default': {'metallicity': 0., 'cloud': 'nc'}},\
    'saumon2008-hybrid': {'file': 'saumon2008-hybrid.csv','citation': 'Saumon et al. (2008) hybrid', 'bibcode': '2008ApJ...689.1327S', 'altname': ['saumon','sau08','saumon08'], 'default': {'metallicity': 0., 'cloud': 'hybrid'}},\
    'saumon2008-f2': {'file': 'saumon2008-f2.csv','citation': 'Saumon et al. (2008) fsed = 2', 'bibcode': '2008ApJ...689.1327S', 'altname': ['saumon','sau08','saumon08'], 'default': {'metallicity': 0., 'cloud': 'f2'}},\
    'marley2019': {'file': 'marley2019.csv','citation': 'Marley et al. (2019)', 'bibcode': 'marley2019', 'url': 'https://zenodo.org/record/1405206#.W4bNhmQzrq0', 'altname': ['marley19','mar19'], 'default': {}},\
    'marley2021': {'file': 'marley2021.csv','citation': 'Marley et al. (2021)', 'bibcode': '2021ApJ...920...85M', 'url': 'https://zenodo.org/record/5063476', 'altname': ['marley','sonora','marley21','mar21','bobcat'], 'default': {'metallicity': 0.}},\
}

# labels for plotting features
FEATURE_LABELS = { \
    'h2o': {'altname': [], 'label': r'H$_2$O', 'type': 'band', 'wavelengths': [[0.925,0.95]*u.micron,[1.08,1.20]*u.micron,[1.325,1.550]*u.micron,[1.72,2.14]*u.micron]}, \
    'ch4': {'altname': [], 'label': r'CH$_4$', 'type': 'band', 'wavelengths': [[1.1,1.24]*u.micron,[1.28,1.44]*u.micron,[1.6,1.76]*u.micron,[2.2,2.35]*u.micron]}, \
    'co': {'altname': [], 'label': r'CO', 'type': 'band', 'wavelengths': [[2.29,2.39]*u.micron]}, \
    'tio': {'altname': [], 'label': r'TiO', 'type': 'band', 'wavelengths': [[0.6569,0.6852]*u.micron,[0.705,0.727]*u.micron,[0.76,0.80]*u.micron,[0.825,0.831]*u.micron,[0.845,0.86]*u.micron]}, \
#        'tio': {'label': r'TiO', 'type': 'band', 'wavelengths': [[0.76,0.80]*u.micron,[0.825,0.831]*u.micron]}, \
    'vo': {'altname': [], 'label': r'VO', 'type': 'band', 'wavelengths': [[1.04,1.08]*u.micron]}, \
    'young vo': {'altname': [], 'label': r'VO', 'type': 'band', 'wavelengths': [[1.17,1.20]*u.micron]}, \
#        'feh': {'label': r'FeH', 'type': 'band', 'wavelengths': [[0.86,0.90]*u.micron,[0.98,1.03],*u.micron[1.19,1.25]*u.micron,[1.57,1.64]*u.micron]}, \
    'cah': {'altname': [], 'label': r'CaH', 'type': 'band', 'wavelengths': [[0.6346,0.639]*u.micron,[0.675,0.705]*u.micron]}, \
    'crh': {'altname': [], 'label': r'CrH', 'type': 'band', 'wavelengths': [[0.8611,0.8681]*u.micron]}, \
    'feh': {'altname': [], 'label': r'FeH', 'type': 'band', 'wavelengths': [[0.8692,0.875]*u.micron,[0.98,1.03]*u.micron,[1.19,1.25]*u.micron,[1.57,1.64]*u.micron]}, \
    'h2': {'altname': ['cia h2'], 'label': r'H$_2$', 'type': 'band', 'wavelengths': [[1.5,2.4]*u.micron]}, \
    'sb': {'altname': ['binary','lt binary','spectral binary'], 'label': r'*', 'type': 'band', 'wavelengths': [[1.6,1.64]*u.micron]}, \
    'h': {'altname': ['hi','h1'], 'label': r'H I', 'type': 'line', 'wavelengths': [[1.004,1.005]*u.micron,[1.093,1.094]*u.micron,[1.281,1.282]*u.micron,[1.944,1.945]*u.micron,[2.166,2.166]*u.micron]},\
    'na': {'altname': ['nai','na1'], 'label': r'Na I', 'type': 'line', 'wavelengths': [[0.8186,0.8195]*u.micron,[1.136,1.137]*u.micron,[2.206,2.209]*u.micron]}, \
    'cs': {'altname': ['csi','cs1'], 'label': r'Cs I', 'type': 'line', 'wavelengths': [[0.8521,0.8521]*u.micron,[0.8943,0.8943]*u.micron]}, \
    'rb': {'altname': ['rbi','rb1'], 'label': r'Rb I', 'type': 'line', 'wavelengths': [[0.78,0.78]*u.micron,[0.7948,0.7948]*u.micron]}, \
    'mg': {'altname': ['mgi','mg1'], 'label': r'Mg I', 'type': 'line', 'wavelengths': [[1.7113336,1.7113336]*u.micron,[1.5745017,1.5770150]*u.micron,[1.4881595,1.4881847,1.5029098,1.5044356]*u.micron,[1.1831422,1.2086969]*u.micron]}, \
    'ca': {'altname': ['cai','ca1'], 'label': r'Ca I', 'type': 'line', 'wavelengths': [[0.6573,0.6573]*u.micron,[2.263110,2.265741]*u.micron,[1.978219,1.985852,1.986764]*u.micron,[1.931447,1.945830,1.951105]*u.micron]}, \
    'caii': {'altname': ['ca2'], 'label': r'Ca II', 'type': 'line', 'wavelengths': [[1.184224,1.195301]*u.micron,[0.985746,0.993409]*u.micron]}, \
    'al': {'altname': ['ali','al1'], 'label': r'Al I', 'type': 'line', 'wavelengths': [[1.672351,1.675511]*u.micron,[1.3127006,1.3154345]*u.micron]}, \
    'fe': {'altname': ['fei','fe1'], 'label': r'Fe I', 'type': 'line', 'wavelengths': [[1.5081407,1.5494570]*u.micron,[1.25604314,1.28832892]*u.micron,[1.14254467,1.15967616,1.16107501,1.16414462,1.16931726,1.18860965,1.18873357,1.19763233]*u.micron]}, \
    'k': {'altname': ['ki','k1'], 'label': r'K I', 'type': 'line', 'wavelengths': [[0.7699,0.7665]*u.micron,[1.169,1.177]*u.micron,[1.244,1.252]*u.micron]}, \
}

EW_SETS = {
    'rojas2012': {'altname': ['rojas','rojas12','roj12'], 'reference': 'Rojas et al. (2012)','bibcode': '2012ApJ...748...93R', 'continuum_fit_order': 1, 'features': {\
        'nai': {'linecenter': 2.207*u.micron,'width': 0.005*u.micron,'continuum': [2.1965,2.201,2.2125,2.2175]*u.micron, 'recenter': False},\
        'cai': {'linecenter': 2.2635*u.micron,'width': 0.0055*u.micron,'continuum': [2.251,2.258,2.2705,2.2760]*u.micron, 'recenter': False},\
    }},
    'terrien2012': {'altname': ['terrien','terrien12','ter12'], 'reference': 'Terrien et al. (2012)','bibcode': '2012ApJ...747L..38T', 'continuum_fit_order': 1, 'features': {\
        'ki-h': {'linecenter': 1.5171*u.micron,'width': 0.0012*u.micron,'continuum': [1.4644,1.4710,1.5206,1.5368]*u.micron, 'recenter': False},\
        'cai-h1': {'linecenter': 1.6159*u.micron,'width': 0.00125*u.micron,'continuum': [1.593,1.6056,1.6248,1.6421]*u.micron, 'recenter': False},\
        'cai-h2': {'linecenter': 1.6203*u.micron,'width': 0.00165*u.micron,'continuum': [1.593,1.6056,1.6248,1.6421]*u.micron, 'recenter': False},\
        'nai-k': {'linecenter': 2.2074*u.micron,'width': 0.00355*u.micron,'continuum': [2.1835,2.1875,2.2335,2.2529]*u.micron, 'recenter': False},\
        'cai-k': {'linecenter': 2.2638*u.micron,'width': 0.0037*u.micron,'continuum': [2.2335,2.2529,2.2717,2.2781]*u.micron, 'recenter': False},\
    }},
    'mann2013': {'altname': ['mann','mann13','man13'], 'reference': 'Mann et al. (2013)', 'bibcode': '2013AJ....145...52M', 'continuum_fit_order': 1, 'features': {\
        'f09': {'linecenter': 1.1396*u.micron,'width': 0.0013*u.micron, 'recenter': False,'continuum': [1.126,1.13,1.153,1.158]*u.micron},\
        'f10': {'linecenter': 1.2698*u.micron,'width': 0.0049*u.micron,'recenter': False,'continuum': [1.255,1.2634,1.27,1.273]*u.micron},\
        'f11': {'linecenter': 1.2908*u.micron,'width': 0.0010*u.micron, 'recenter': False,'continuum': [1.27,1.273,1.295,1.297]*u.micron},\
        'f12': {'linecenter': 1.3148*u.micron,'width': 0.0025*u.micron, 'recenter': False,'continuum': [1.304,1.307,1.3214,1.327]*u.micron},\
        'f13': {'linecenter': 1.3344*u.micron,'width': 0.00115*u.micron, 'recenter': False,'continuum': [1.3214,1.327,1.409,1.415]*u.micron},\
        'f14': {'linecenter': 1.4766*u.micron,'width': 0.00205*u.micron, 'recenter': False,'continuum': [1.4644,1.471,1.4921,1.4965]*u.micron},\
        'f15': {'linecenter': 1.4836*u.micron,'width': 0.00115*u.micron, 'recenter': False,'continuum': [1.4644,1.471,1.4921,1.4965]*u.micron},\
        'f16': {'linecenter': 1.5172*u.micron,'width': 0.00165*u.micron, 'recenter': False,'continuum': [1.506,1.509,1.519,1.522]*u.micron},\
        'f17': {'linecenter': 1.6158*u.micron,'width': 0.00125*u.micron, 'recenter': False,'continuum': [1.592,1.596,1.623,1.631]*u.micron},\
        'f18': {'linecenter': 1.7261*u.micron,'width': 0.0016*u.micron, 'recenter': False,'continuum': [1.6935,1.698,1.753,1.757]*u.micron},\
        'f19': {'linecenter': 2.2079*u.micron,'width': 0.0034*u.micron, 'recenter': False,'continuum': [2.194,2.1985,2.213,2.219]*u.micron},\
        'f20': {'linecenter': 2.3242*u.micron,'width': 0.0019*u.micron, 'recenter': False,'continuum': [2.305,2.3105,2.36,2.364]*u.micron},\
        'f21': {'linecenter': 2.3342*u.micron,'width': 0.00175*u.micron, 'recenter': False,'continuum': [2.305,2.3105,2.36,2.364]*u.micron},\
        'f22': {'linecenter': 2.3844*u.micron,'width': 0.00175*u.micron, 'recenter': False,'continuum': [2.371,2.376,2.395,2.405]*u.micron},\
    }},
    'mann2014': {'altname': ['mann14','man14'], 'reference': 'Mann et al. (2014)','bibcode': '2014AJ....147..160M', 'continuum_fit_order': 1, 'features': {\
        'nai': {'linecenter': 2.2079*u.micron,'width': 0.0034*u.micron,'continuum': [2.194,2.1985,2.213,2.219]*u.micron, 'recenter': False},\
        'cai': {'linecenter': 2.264*u.micron,'width': 0.0030*u.micron,'continuum': [2.245,2.252,2.2717,2.2781]*u.micron, 'recenter': False},\
    }},
    'newton2015': {'altname': ['newton15','new15'], 'reference': 'Newton et al. (2015)','bibcode': '2015ApJ...800...85N', 'continuum_fit_order': 1, 'features': {\
        'mg1': {'linecenter': 1.4882*u.micron,'width': 0.002/2*u.micron,'continuum': [1.4790,1.4850,1.4900,1.4950]*u.micron, 'recenter': False},\
        'mg2': {'linecenter': 1.5040*u.micron,'width': 0.004/2*u.micron,'continuum': [1.4965,1.5000,1.5070,1.5120]*u.micron, 'recenter': False},\
        'k'  : {'linecenter': 1.5170*u.micron,'width': 0.002/2*u.micron,'continuum': [1.5105,1.5135,1.5185,1.5215]*u.micron, 'recenter': False},\
        'mg3': {'linecenter': 1.5760*u.micron,'width': 0.004/2*u.micron,'continuum': [1.5640,1.5680,1.5785,1.5815]*u.micron, 'recenter': False},\
        'si' : {'linecenter': 1.59025*u.micron,'width': 0.0045/2*u.micron,'continuum': [1.5845,1.5875,1.5925,1.5955]*u.micron, 'recenter': False},\
        'co1' : {'linecenter': 1.6205*u.micron,'width': 0.003/2*u.micron,'continuum': [1.6120,1.6150,1.6265,1.6295]*u.micron, 'recenter': False},\
        'co2' : {'linecenter': 1.6255*u.micron,'width': 0.002/2*u.micron,'continuum': [1.6120,1.6150,1.6265,1.6295]*u.micron, 'recenter': False},\
        'al-a' : {'linecenter': 1.6725*u.micron,'width': 0.002/2*u.micron,'continuum': [1.6550,1.6650,1.6780,1.6820]*u.micron, 'recenter': False},\
        'al-b' : {'linecenter': 1.676*u.micron,'width': 0.003/2*u.micron,'continuum': [1.6550,1.6650,1.6780,1.6820]*u.micron, 'recenter': False},\
        'mg' : {'linecenter': 1.71125*u.micron,'width': 0.0025/2*u.micron,'continuum': [1.7025,1.7055,1.7130,1.7160]*u.micron, 'recenter': False},\
    }},
}

# all Mann continuua:
            # 'continuum': [[0.979,0.989],[1.061,1.065],[1.126,1.13],[1.153,1.158],[1.189,1.193],[1.214,1.218],[1.225,1.230],[1.255,1.2634],[1.27,1.273],[1.295,1.297],[1.304,1.307],[1.3214,1.327],[1.409,1.415]]*u.micron},\
            # 'continuum': [[1.444,1.448],[1.4644,1.471],[1.4921,1.4965],[1.506,1.509],[1.519,1.522],[1.592,1.596],[1.623,1.631],[1.6935,1.698],[1.753,1.757]]*u.micron},\
            # 'continuum': [[1.886,1.89],[1.935,1.94],[1.961,1.97],[2.05,2.054],[2.08,2.087],[2.133,2.1351],[2.153,2.159],[2.167,2.172],[2.194,2.1985],[2.213,2.219],[2.245,2.252],[2.2717,2.2781],[2.285,2.29],[2.305,2.3105],[2.36,2.364],[2.371,2.376],[2.395,2.405]]*u.micron},\


# index definitions
# NOTE: SOME OF THESE APPEAR TO BE INCOMPLETE, AND ROJAS IS ODDLY DEFINED
# NEED TO ADD: LEPINE, GIZIS
INDEX_SETS = {
    'allers2013': {'altname': ['allers','all13','allers13'], 'bibcode': '2013ApJ...772...79A', 'indices': {\
        'H2O': {'ranges': ([1.55,1.56]*u.micron,[1.492,1.502]*u.micron), 'method': 'ratio', 'sample': 'average'},\
        'FeH-z': {'ranges': ([0.99135,1.00465]*u.micron,[0.97335,0.98665]*u.micron,[1.01535,1.02865]*u.micron), 'method': 'allers', 'sample': 'average'},\
        'VO-z': {'ranges': ([1.05095,1.06505]*u.micron,[1.02795,1.04205]*u.micron,[1.07995,1.09405]*u.micron), 'method': 'allers', 'sample': 'average'},\
        'FeH-J': {'ranges': ([1.19880,1.20120]*u.micron,[1.19080,1.19320]*u.micron,[1.20680,1.20920]*u.micron), 'method': 'allers', 'sample': 'average'},\
        'KI-J': {'ranges': ([1.23570,1.25230]*u.micron,[1.21170,1.22830]*u.micron,[1.26170,1.27830]*u.micron), 'method': 'allers', 'sample': 'average'},\
        'H-cont': {'ranges': ([1.54960,1.57040]*u.micron,[1.45960,1.48040]*u.micron,[1.65960,1.68040]*u.micron), 'method': 'allers', 'sample': 'average'},\
    }},
    'burgasser2006': {'altname': ['burgasser','burgasser06','bur06'], 'bibcode': '2006ApJ...637.1067B', 'indices': {\
        'CH4-J': {'ranges': ([1.315,1.335]*u.micron,[1.26,1.285]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'CH4-H': {'ranges': ([1.635,1.675]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'CH4-K': {'ranges': ([2.215,2.255]*u.micron,[2.08,2.12]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O-J': {'ranges': ([1.14,1.165]*u.micron,[1.26,1.285]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O-H': {'ranges': ([1.48,1.52]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O-K': {'ranges': ([1.975,1.995]*u.micron,[2.08,2.10]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'K/J': {'ranges': ([2.06,2.10]*u.micron,[1.25,1.29]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H-dip': {'ranges': ([1.61,1.64]*u.micron,[1.56,1.59]*u.micron,[1.66,1.69]*u.micron), 'method': 'inverse_line', 'sample': 'integrate'},\
    }},
    'burgasser2023': {'altname': ['burgasser23','bur23'], 'bibcode': '', 'indices': {\
        'CH4-J': {'ranges': ([1.315,1.335]*u.micron,[1.26,1.285]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'CH4-H': {'ranges': ([1.635,1.675]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'CH4-K': {'ranges': ([2.215,2.255]*u.micron,[2.08,2.12]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'H2O-J': {'ranges': ([1.14,1.165]*u.micron,[1.26,1.285]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'H2O-H': {'ranges': ([1.48,1.52]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'H2O-K': {'ranges': ([1.975,1.995]*u.micron,[2.08,2.10]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'Y/J': {'ranges': ([1.005,1.045]*u.micron,[1.25,1.29]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'K/J': {'ranges': ([2.06,2.10]*u.micron,[1.25,1.29]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'K/H': {'ranges': ([2.06,2.10]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'H-dip': {'ranges': ([1.61,1.64]*u.micron,[1.56,1.59]*u.micron,[1.66,1.69]*u.micron), 'method': 'inverse_line', 'sample': 'median'},\
    }},
    'burgasser2025': {'altname': ['burgasser23','bur23'], 'bibcode': '', 'indices': {\
        'CH4-1.3': {'ranges': ([1.315,1.335]*u.micron,[1.26,1.285]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'CH4-1.6': {'ranges': ([1.635,1.675]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'CH4-3.3': {'ranges': ([3.22,3.42]*u.micron,[3.93,4.13]*u.micron), 'method': 'ratio', 'sample': 'average'},
        'CO-2.3': {'ranges': ([2.33,2.37]*u.micron,[2.12,2.16]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'CO-4.6': {'ranges': ([4.5,4.6]*u.micron,[4.33,4.43]*u.micron), 'method': 'ratio', 'sample': 'average'},
        'CO2-4.2': {'ranges': ([4.22,4.32]*u.micron,[4.05,4.15]*u.micron), 'method': 'ratio', 'sample': 'average'},
        'H2O-1.1': {'ranges': ([1.14,1.165]*u.micron,[1.26,1.285]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'H2O-1.5': {'ranges': ([1.48,1.52]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'H2O-2.0': {'ranges': ([1.975,1.995]*u.micron,[2.08,2.10]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'H2O-5.8': {'ranges': ([5.65,5.95]*u.micron,[6.6,6.9]*u.micron,[6.1,6.4]*u.micron), 'weights':[1,0.562,0.474],'method': 'sumnum', 'sample': 'average'},\
        'CH4-7.6': {'ranges': ([7.35,7.95]*u.micron,[9.7,10.3]*u.micron), 'method': 'ratio', 'sample': 'average'},\
        'NH3-11': {'ranges': ([10.5,11.1]*u.micron,[9.7,10.3]*u.micron), 'method': 'ratio', 'sample': 'average'},\
        'SIH4-4.6': {'ranges': ([4.55,4.65]*u.micron,[4.36,4.]*u.micron), 'method': 'ratio', 'sample': 'average'},
        'SO2-4.0': {'ranges': ([3.95,4.05]*u.micron,[3.8,3.9]*u.micron), 'method': 'ratio', 'sample': 'average'},
        'Silicate': {'ranges': ([8.7,9.3]*u.micron,[7.2,7.8]*u.micron,[11.2,11.8]*u.micron), 'method': 'fitdenom', 'sample': 'average', 'fitorder': 1},\
        'PH3': {'ranges': ([4.285,4.335]*u.micron,[4.22,4.27]*u.micron,[4.36,4.41]*u.micron), 'method': 'avegdenom', 'sample': 'average'},
        'Y/J': {'ranges': ([1.005,1.045]*u.micron,[1.25,1.29]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'K/J': {'ranges': ([2.06,2.10]*u.micron,[1.25,1.29]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'K/H': {'ranges': ([2.06,2.10]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'median'},\
        'M/J': {'ranges': ([4.1,4.15]*u.micron,[1.24,1.29]*u.micron,), 'method': 'ratio', 'sample': 'average'},
        'M/K': {'ranges': ([4.1,4.15]*u.micron,[2.06,2.11]*u.micron,), 'method': 'ratio', 'sample': 'average'},
        'M/L': {'ranges': ([4.1,4.15]*u.micron,[3,3.05]*u.micron,), 'method': 'ratio', 'sample': 'average'},
        'H-dip': {'ranges': ([1.61,1.64]*u.micron,[1.56,1.59]*u.micron,[1.66,1.69]*u.micron), 'method': 'inverse_line', 'sample': 'median'},\
    }},
    'bardalez2014': {'altname': ['bardalez','bardalez14','bar14'], 'bibcode': '2014ApJ...794..143B', 'indices': {\
        'H2O-J': {'ranges': ([1.14,1.165]*u.micron,[1.26,1.285]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'CH4-J': {'ranges': ([1.315,1.335]*u.micron,[1.26,1.285]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O-H': {'ranges': ([1.48,1.52]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'CH4-H': {'ranges': ([1.635,1.675]*u.micron,[1.56,1.60]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O-K': {'ranges': ([1.975,1.995]*u.micron,[2.08,2.12]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'CH4-K': {'ranges': ([2.215,2.255]*u.micron,[2.08,2.12]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'K/J': {'ranges': ([2.06,2.10]*u.micron,[1.25,1.29]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H-dip': {'ranges': ([1.61,1.64]*u.micron,[1.56,1.59]*u.micron,[1.66,1.69]*u.micron), 'method': 'inverse_line', 'sample': 'integrate'},\
        'K-slope': {'ranges': ([2.06,2.10]*u.micron,[2.10,2.14]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'J-slope': {'ranges': ([1.27,1.30]*u.micron,[1.30,1.33]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'J-curve': {'ranges': ([1.04,1.07]*u.micron,[1.26,1.29]*u.micron,[1.14,1.17]*u.micron), 'method': 'line', 'sample': 'integrate'},\
        'H-bump': {'ranges': ([1.54,1.57]*u.micron,[1.66,1.69]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O-Y': {'ranges': ([1.04,1.07]*u.micron,[1.14,1.17]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
    }},
    'canty2013': {'altname': ['canty','canty13','can13'], 'bibcode': '2010ApJ...722..971C', 'indices': {\
        'H2(K)': {'ranges': ([2.16,2.18]*u.micron,[2.23,2.25]*u.micron), 'method': 'ratio', 'sample': 'median'},\
    }},
    'covey2010': {'altname': ['covey','covey10','cov10'], 'bibcode': '2010ApJ...722..971C', 'indices': {\
        'H2O-H': {'ranges': ([1.595,1.615]*u.micron,[1.68,1.70]*u.micron,[1.76,1.78]*u.micron), 'method': 'doubleratio', 'sample': 'average'},\
        'H2O-K': {'ranges': ([2.18,2.2]*u.micron,[2.27,2.29]*u.micron,[2.36,2.38]*u.micron), 'method': 'doubleratio', 'sample': 'average'},\
    }},
    'cushing2006': {'altname': ['cus06','cushing06'], 'bibcode': '2006ApJ...648..614C', 'indices': {\
        'IRS-H2O': {'ranges': ([6.175,6.325]*u.micron,[5.725,5.875]*u.micron,[6.675,6.825]*u.micron), 'weights':[1,0.562,0.474],'method': 'sumdenom', 'sample': 'average'},\
        'IRS-CH4': {'ranges': ([9.85,10.15]*u.micron,[8.35,8.65]*u.micron), 'method': 'ratio', 'sample': 'average'},\
        'IRS-NH3': {'ranges': ([9.85,10.15]*u.micron,[10.65,10.95]*u.micron), 'method': 'ratio', 'sample': 'average'},\
    }},
    'geballe2002': {'altname': ['geballe','geballe02','geb02'], 'bibcode': '2002ApJ...564..466G', 'indices': {\
        'Color-d2': {'ranges': ([0.96,0.98]*u.micron,[0.735,0.755]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'Cont-1.0': {'ranges': ([1.04,1.05]*u.micron,[0.875,0.885]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O-1.2': {'ranges': ([1.26,1.29]*u.micron,[1.13,1.16]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O-1.5': {'ranges': ([1.57,1.59]*u.micron,[1.46,1.48]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O-2.2': {'ranges': ([2.09,2.11]*u.micron,[1.975,1.995]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'CH4-1.6': {'ranges': ([1.56,1.6]*u.micron,[1.635,1.675]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'CH4-2.2': {'ranges': ([2.08,2.12]*u.micron,[2.215,2.255]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
    }},
    'kirkpatrick1999': {'altname': ['kirkpatrick','kirkpatrick19','kir99'], 'bibcode': '1999ApJ...519..802K', 'indices': {\
        'Rb-a': {'ranges': ([.77752,.77852]*u.micron,[.78152,.78252]*u.micron,[.77952,.78052]*u.micron), 'method': 'line', 'sample': 'sum'},\
        'Rb-b': {'ranges': ([.79226,.79326]*u.micron,[.79626,.79726]*u.micron,[.79426,.79526]*u.micron), 'method': 'line', 'sample': 'sum'},\
        'Na-a': {'ranges': ([.81533,.81633]*u.micron,[.81783,.81883]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'Na-b': {'ranges': ([.81533,.81633]*u.micron,[.81898,.81998]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'Cs-a': {'ranges': ([.84961,.85061]*u.micron,[.85361,.85461]*u.micron,[.85161,.85261]*u.micron), 'method': 'line', 'sample': 'sum'},\
        'Cs-b': {'ranges': ([.89185,.89285]*u.micron,[.89583,.89683]*u.micron,[.89385,.89485]*u.micron), 'method': 'line', 'sample': 'sum'},\
        'TiO-a': {'ranges': ([.7033,.7048]*u.micron,[.7058,.7073]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'TiO-b': {'ranges': ([.8400,.8415]*u.micron,[.8435,.8470]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'VO-a': {'ranges': ([.7350,.7370]*u.micron,[.7550,.7570]*u.micron,[.7430,.7470]*u.micron), 'method': 'sumnum', 'sample': 'sum'},\
        'VO-b': {'ranges': ([.7860,.7880]*u.micron,[.8080,.8100]*u.micron,[.7960,.8000]*u.micron), 'method': 'sumnum', 'sample': 'sum'},\
        'CrH-a': {'ranges': ([.8580,.8600]*u.micron,[.8621,.8641]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'CrH-b': {'ranges': ([.9940,.9960]*u.micron,[.9970,.9990]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'FeH-a': {'ranges': ([.8660,.8680]*u.micron,[.8700,.8720]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'FeH-b': {'ranges': ([.9863,.9883]*u.micron,[.9908,.9928]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'Color-a': {'ranges': ([.9800,.9850]*u.micron,[.7300,.7350]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'Color-b': {'ranges': ([.9800,.9850]*u.micron,[.7000,.7050]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'Color-c': {'ranges': ([.9800,.9850]*u.micron,[.8100,.8150]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
        'Color-d': {'ranges': ([.9675,.9850]*u.micron,[.7350,.7550]*u.micron), 'method': 'ratio', 'sample': 'sum'},\
    }},
    'mann2013': {'altname': ['mann','mann13','man13'], 'bibcode': '2013AJ....145...52M', 'indices': {\
        'H2O-J': {'ranges': ([1.210,1.230]*u.micron,[1.313,1.333]*u.micron,[1.311,1.351]*u.micron), 'method': 'doubleratio', 'sample': 'median'},\
    }},
    'martin1999': {'altname': ['martin','martin99','mar99'], 'bibcode': '1999AJ....118.2466M', 'indices': {\
        'PC3': {'ranges': ([.823,.827]*u.micron,[.754,.758]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'PC6': {'ranges': ([.909,.913]*u.micron,[.650,.654]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'CrH1': {'ranges': ([.856,.860]*u.micron,[.861,.865]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'CrH2': {'ranges': ([.984,.988]*u.micron,[.997,1.001]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'FeH1': {'ranges': ([.856,.860]*u.micron,[.8685,.8725]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'FeH2': {'ranges': ([.984,.988]*u.micron,[.990,.994]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'H2O1': {'ranges': ([.919,.923]*u.micron,[.928,.932]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'TiO1': {'ranges': ([.700,.704]*u.micron,[.706,.710]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'TiO2': {'ranges': ([.838,.842]*u.micron,[.844,.848]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'VO1': {'ranges': ([.754,.758]*u.micron,[.742,.746]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
        'VO2': {'ranges': ([.799,.803]*u.micron,[.790,.794]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
    }},
    'mclean2003': {'altname': ['mclean','mclean03','mcl03'], 'bibcode': '2003ApJ...596..561M', 'indices': {\
        'H2OD': {'ranges': ([1.951,1.977]*u.micron,[2.062,2.088]*u.micron), 'method': 'ratio', 'sample': 'average'},\
    }},
    'reid2001': {'altname': ['reid','reid01','rei01'], 'bibcode': '2001AJ....121.1710R', 'indices': {\
        'H2O-A': {'ranges': ([1.33,1.35]*u.micron,[1.28,1.30]*u.micron), 'method': 'ratio', 'sample': 'average'},\
        'H2O-B': {'ranges': ([1.47,1.49]*u.micron,[1.59,1.61]*u.micron), 'method': 'ratio', 'sample': 'average'},\
    }},
    'rojas2012': {'altname': ['rojas','rojas12','roj12'], 'bibcode': '2012ApJ...748...93R', 'indices': {\
        'H2O-K2': {'ranges': ([2.070,2.090]*u.micron,[2.235,2.255]*u.micron,[2.360,2.380]*u.micron), 'method': 'doubleratio', 'sample': 'median'},\
    }},
    'slesnick2004': {'altname': ['slesnick','slesnick04','sle04'], 'bibcode': '2004ApJ...610.1045S', 'indices': {\
        'H2O-1': {'ranges': ([1.335,1.345]*u.micron,[1.295,1.305]*u.micron), 'method': 'ratio', 'sample': 'average'},\
        'H2O-2': {'ranges': ([2.035,2.045]*u.micron,[2.145,2.155]*u.micron), 'method': 'ratio', 'sample': 'average'},\
        'FeH': {'ranges': ([1.1935,1.2065]*u.micron,[1.2235,1.2365]*u.micron), 'method': 'ratio', 'sample': 'average'},\
    }},
    'suarez2022': {'altname': ['sua22','suarez22'], 'bibcode': '2022MNRAS.513.5701S', 'indices': {\
        'H2O': {'ranges': ([6.1,6.4]*u.micron,[5.65,5.95]*u.micron,[6.6,6.9]*u.micron), 'weights':[1,0.562,0.474],'method': 'sumdenom', 'sample': 'average'},\
        'CH4': {'ranges': ([9.7,10.3]*u.micron,[7.35,7.95]*u.micron), 'method': 'ratio', 'sample': 'average'},\
        'NH3': {'ranges': ([9.7,10.3]*u.micron,[10.5,11.1]*u.micron), 'method': 'ratio', 'sample': 'average'},\
        'Silicate': {'ranges': ([7.2,7.8]*u.micron,[11.2,11.8]*u.micron,[8.7,9.3]*u.micron), 'method': 'fitnum', 'sample': 'average', 'fitorder': 1},\
    }},
    'testi2001': {'altname': ['testi','testi01','tes01'], 'bibcode': '2001ApJ...552L.147T', 'indices': {\
        'sHJ': {'ranges': ([1.265,1.305]*u.micron,[1.6,1.7]*u.micron), 'method': 'change', 'sample': 'average'},\
        'sKJ': {'ranges': ([1.265,1.305]*u.micron,[2.12,2.16]*u.micron), 'method': 'change', 'sample': 'average'},\
        'sH2O_J': {'ranges': ([1.265,1.305]*u.micron,[1.09,1.13]*u.micron), 'method': 'change', 'sample': 'average'},\
        'sH2O_H1': {'ranges': ([1.60,1.70]*u.micron,[1.45,1.48]*u.micron), 'method': 'change', 'sample': 'average'},\
        'sH2O_H2': {'ranges': ([1.60,1.70]*u.micron,[1.77,1.81]*u.micron), 'method': 'change', 'sample': 'average'},\
        'sH2O_K': {'ranges': ([2.12,2.16]*u.micron,[1.96,1.99]*u.micron), 'method': 'change', 'sample': 'average'},\
    }},
    'tokunaga1999': {'altname': ['tokunaga','tokunaga99','tok99'], 'bibcode': '1999AJ....117.1010T', 'indices': {\
        'K1': {'ranges': ([2.1,2.18]*u.micron,[1.96,2.04]*u.micron), 'method': 'change', 'sample': 'average'},\
        'K2': {'ranges': ([2.2,2.28]*u.micron,[2.1,2.18]*u.micron), 'method': 'change', 'sample': 'average'},\
    }},
}


# classification indices
INDEX_CLASSIFICATION_RELATIONS = {
    'reid2001': {'altname': ['reid','reid01','rei01'], 'bibcode': '2001AJ....121.1710R', 'method': 'polynomial', 'sptoffset': 20., 'decimal': False, 'min_indices': 1, 'sets': ['reid2001'], 'indices': { \
            'H2O-A': {'fitunc': 1.18, 'range': [18,26], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [-32.1, 23.4]}, \
            'H2O-B': {'fitunc': 1.02, 'range': [18,28], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [-24.9, 20.7]},
    }},
    'testi2001': {'altname': ['testi','testi01','tes01'], 'bibcode': '2001ApJ...552L.147T', 'method': 'polynomial', 'sptoffset': 10., 'decimal': True, 'min_indices': 2, 'sets': ['testi2001'], 'indices': { \
            'sHJ': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [-1.87, 1.67]}, \
            'sKJ': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [-1.20, 2.01]}, \
            'sH2O_J': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [1.54, 0.98]}, \
            'sH2O_H1': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [1.27, 0.76]}, \
            'sH2O_H2': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [2.11, 0.29]}, \
            'sH2O_K': {'fitunc': 0.5, 'range': [20,26], 'spt': 0., 'sptunc': 99., 'mask': 1.,  'coeff': [2.36, 0.60]},
    }},
    'slesnick2004': {'altname': ['slesnick','slesnick04','sle04'], 'bibcode': '2004ApJ...610.1045S', 'method': 'polynomial', 'sptoffset': 10., 'decimal': False, 'min_indices': 1, 'sets': ['slesnick2004'], 'indices': { \
            'H2O-1': {'fitunc': 1.2, 'range': [10,30], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [-35.35,33.71]}, \
            'H2O-2': {'fitunc': 0.53, 'range': [12,23], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [-27.10,34.13]}, \
            'FeH': {'fitunc': 0.66, 'range': [13,23], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [-31.59,36.31]}, \
    }},
    'burgasser2007': {'altname': ['burgasser','burgasser07','bur07'], 'bibcode': '2007ApJ...659..655B', 'method': 'polynomial', 'sptoffset': 20., 'decimal': False, 'min_indices': 2, 'sets': ['burgasser2006'], 'indices': { \
            'H2O-J': {'fitunc': 0.8, 'range': [20,39], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [1.038e2, -2.156e2,  1.312e2, -3.919e1, 1.949e1]}, \
            'H2O-H': {'fitunc': 1.0, 'range': [20,39], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [9.087e-1, -3.221e1, 2.527e1, -1.978e1, 2.098e1]}, \
            'CH4-J': {'fitunc': 0.7, 'range': [30,39], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [1.491e2, -3.381e2, 2.424e2, -8.450e1, 2.708e1]}, \
            'CH4-H': {'fitunc': 0.3, 'range': [31,39], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [2.084e1, -5.068e1, 4.361e1, -2.291e1, 2.013e1]}, \
            'CH4-K': {'fitunc': 1.1, 'range': [20,37], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [-1.259e1, -4.734e0, 2.534e1, -2.246e1, 1.885e1]}, \
    }},
    'geballe2002': {'altname': ['geballe','geballe02','geb02'], 'bibcode': '2002ApJ...564..466G', 'method': 'ranges', 'sptoffset': 20., 'decimal': False, 'min_indices': 2, 'sets': ['geballe2002','martin1999'], 'indices': { \
            'PC3': [[2.4,2.6,20.],[2.6,2.86,21.],[2.85,3.25,22.],[3.25,4.25,23.],[4.25,6,24.]],\
            'Color-d2': [[4.5,5.5,20.],[5.5,6.5,21.],[6.5,7.5,22.],[7.5,10.,23.],[10,17,24.],[17.,23.,25.],[23.,25.,26.]],\
            'H2O-1.2': [[1.5,1.7,30.],[1.7,1.9,31.],[1.9,2.15,32.],[2.15,2.5,33.],[2.5,3.0,34.],[3.0,4.5,35.],[4.5,6.5,36.],[6.5,10.,37.],[10.,15.,38.]],\
            'H2O-1.5': [[1.2,1.27,20.],[1.27,1.35,21.],[1.35,1.43,22.],[1.43,1.5,23.],[1.5,1.55,24.],[1.55,1.6,25.],[1.6,1.65,26.],[1.65,1.7,27.],[1.7,1.8,28.],[1.8,1.95,29.],[1.95,2.2,30.],[2.2,2.5,31.],[2.5,3.0,32.],[3.0,3.5,33.],[3.5,4.5,34.],[4.5,5.5,35.],[5.5,7.,36.],[7.,9.,37.],[9.,12.,38.]],\
            'CH4-1.6': [[1.02,1.07,30.],[1.07,1.15,31.],[1.15,1.3,32.],[1.3,1.5,33.],[1.5,1.8,34.],[1.8,2.5,35.],[2.5,4,36.],[4.,6.,37.],[6.,9.,38.]],\
            'CH4-2.2': [[0.91,0.94,23.],[0.94,0.98,24.],[0.98,1.025,25.],[1.025,1.075,26.],[1.075,1.125,27.],[1.125,1.175,28.],[1.175,1.25,29.],[1.25,1.4,30.],[1.4,1.6,31.],[1.6,1.95,32.],[1.95,2.75,33.],[2.75,3.8,34.],[3.8,5.5,35.],[5.5,8.5,36.],[8.5,12.,37],[12.,18.,38.]]
    }},
    'allers2013': {'altname': ['allers','allers13','all13'], 'bibcode': '2013ApJ...657..511A', 'method': 'polynomial', 'sptoffset': 10., 'decimal': False, 'min_indices': 2, 'sets': ['allers2013','mclean2003','slesnick2004'], 'indices': { \
             'H2O': {'fitunc': 0.390, 'range': [15,25], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [24.0476, -104.424, 169.388,-83.5437]}, \
             'H2O-1': {'fitunc': 1.097, 'range': [14,25], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [28.5982, -80.7404, 39.3513, 12.1927]}, \
             'H2OD': {'fitunc': 0.757, 'range': [20,28], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [-97.230, 229.884, -202.245, 79.4477]}, \
             'H2O-2': {'fitunc': 0.501, 'range': [14,22], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [37.5013, -97.8144, 55.4580, 10.8822]},\
             # 'jtype': {'fitunc': 0., 'range': [14,28], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [37.5013, -97.8144, 55.4580, 10.8822]},\ # NOTE THESE ADDITIONS BROKE THE CODE
             # 'ktype': {'fitunc': 0., 'range': [14,28], 'spt': 0., 'sptunc': 99., 'mask': 1., 'coeff': [37.5013, -97.8144, 55.4580, 10.8822]},\
    }},
}

# indices for spectral binary identification
SPECTRAL_BINARY_INDICES = {
    'burgasser2010': {'altname': ['burgasser','burgasser10','bur10'], 'bibcode': '2010ApJ...710.1142B', 'index_set': 'burgasser2006', 'spt': True, 'sptoffset': 0, 'spt_range':[27,35],'strong':3,'weak':2, 'relations': [ \
            {'indices': ['H2O-J','H2O-K'], 'coeff': [0.2/0.325,0.3], 'xlim':[0.325,0.65],'measure': 'individual', 'direction': 'high'},\
            {'indices': ['CH4-H','CH4-K'], 'coeff': [0.425/0.4,-0.2875],'xlim':[0.6,1.0], 'measure': 'individual', 'direction': 'high'},\
            {'indices': ['CH4-H','K/J'], 'coeff': [0.125/0.35,0.01786],'xlim':[0.65,1.0], 'measure': 'individual', 'direction': 'high'},\
            {'indices': ['H2O-H','H-dip'], 'coeff': [0.49],'xlim':[0.5,0.875],'measure': 'individual', 'direction': 'low'},\
            {'indices': ['SPT','H2O-J','H2O-H'], 'coeff': [-0.075/2,2.10625],'xlim':[28.5,33.5],'ylim':[0.18,0.925], 'measure': 'ratio2', 'direction': 'low'},\
            {'indices': ['SPT','H2O-J','CH4-K'], 'coeff': [0.2/6.,-0.325],'xlim':[28.5,34.5], 'measure': 'ratio2', 'direction': 'low'},\
            # {'indices': ['H2O-J','H2O-K'], 'points': ((0.325,0.5),(0.65,0.7)), 'method': 'individual', 'direction': 'y-high'},\
            # {'indices': ['CH4-H','CH4-K'], 'points': ((0.6,0.35),(1.0,0.775)), 'method': 'individual', 'direction': 'y-high'},\
            # {'indices': ['CH4-H','K/J'], 'points': ((0.65,0.25),(1.0,0.375)), 'method': 'individual', 'direction': 'y-high'},\
            # {'indices': ['H2O-H','H-dip'], 'points': ((0.5,0.49),0.875,0.49), 'method': 'individual', 'direction': 'y-low'},\
            # {'indices': ['SPT','H2O-J','H2O-H'], 'points': ((28.5,0.925),(31.5,0.925),(33.5,0.85)), 'method': 'ratio2', 'direction': 'y-low'},\
            # {'indices': ['SPT','H2O-J','CH4-K'], 'points': ((28.5,0.625),(34.5,0.825)), 'method': 'ratio2', 'direction': 'y-low'},\
    ]},
    'bardalez2014': {'altname': ['bardalez','bardalez14','bar14'], 'bibcode': '2014ApJ...794..143B', 'index_set': 'bardalez2014', 'spt': True, 'sptoffset': 0, 'spt_range':[17,28], 'strong':8,'weak':4,'relations': [ \
            {'indices': ['SPT','CH4-H'], 'coeff': [-0.00043,0.0253,0.7178-0.0354], 'xlim': [17,29],'measure': 'individual', 'direction': 'low'},\
            {'indices': ['H2O-J','CH4-H'], 'coeff': [-0.08,1.09], 'xlim': [-10,0.9],'measure': 'individual', 'direction': 'low'},\
            {'indices': ['H2O-J','H-bump'], 'coeff': [0.16,0.806], 'xlim': [-10,0.9],'measure': 'individual', 'direction': 'high'},\
            {'indices': ['CH4-J','CH4-H'], 'coeff': [-0.56,1.41], 'xlim': [-10,0.74],'ylim':[-10,1.04],'measure': 'individual', 'direction': 'low'},\
            {'indices': ['CH4-J','H-bump'], 'coeff': [1.0,0.24], 'xlim': [-10,0.74],'ylim':[0.91,10.],'measure': 'individual', 'direction': 'high'},\
            {'indices': ['CH4-H','J-slope'], 'coeff': [1.250,0.-207], 'xlim': [-10,1.03],'ylim':[1.03,10.],'measure': 'individual', 'direction': 'high'},\
            {'indices': ['CH4-H','J-curve'], 'coeff': [1.245,-1.565,2.224+0.088], 'xlim': [-10,1.03],'ylim':[1.97,10.],'measure': 'individual', 'direction': 'high'},\
            {'indices': ['CH4-H','H-bump'], 'coeff': [1.36,-4.26,3.89-0.013], 'ylim':[-10,1.03],'measure': 'individual', 'direction': 'low'},\
            {'indices': ['J-slope','H-dip'], 'coeff': [0.20,0.27], 'xlim':[1.03,10.],'measure': 'individual', 'direction': 'low'},\
            {'indices': ['J-slope','H-bump'], 'coeff': [-2.75,3.84], 'ylim':[0.91,10.],'measure': 'individual', 'direction': 'high'},\
            {'indices': ['K-slope','H2O-Y'], 'coeff': [12.036,-20,8.973+0.064], 'xlim':[0.93,0.96],'measure': 'individual', 'direction': 'high'},\
            {'indices': ['J-curve','H-bump'], 'coeff': [0.269,-1.326,2.479+0.048], 'xlim':[2.0,10.],'ylim':[0.92,10],'measure': 'individual','direction': 'high'},\
    ]},
}
# Empirical relations - SpT to Teff
SPT_TEFF_RELATIONS = {
    'dupuy_saumon': {'altname': ['dupuy','dupuy17','dupuy2017','dup17','dupuy17_saumon','dupuy2017_saumon','dup17-saumon'], 'reference': 'Dupuy et al. (2017)','bibcode': '2017ApJS..231...15D',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [6.001,-284.52,4544.3], 'range': [21.5,35.],'fitunc': 80.},
    'dupuy_lyon': {'altname': ['dupuy17_lyon','dupuy2017_lyon','dup17_lyon'], 'reference': 'Dupuy et al. (2017)','bibcode': '2017ApJS..231...15D',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [4.582,-238.03,4251.0], 'range': [17.,35.],'fitunc': 90.},
    'faherty': {'altname': ['faherty16','faherty2016','fah16'], 'reference': 'Faherty et al. (2016)','bibcode': '2016ApJS..225...10F',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [1.546e-4,-1.606e-2,6.318e-1,-1.191e1,1.155e2,-7.005e2,4.747e3], 'range': [17.,38.],'fitunc': 113.},
    'faherty_young': {'altname': ['faherty16_young','faherty2016_young','fah16yng'], 'reference': 'Faherty et al. (2016)','bibcode': '2016ApJS..225...10F',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [1.330,-6.68637e1,1.23542e3,-1.00688e4,3.27664e4], 'range': [17.,27.],'fitunc': 180.},
    'faherty_young2': {'altname': ['faherty16_young2','faherty2016_young2','fah16yng2'], 'reference': 'Faherty et al. (2016)','bibcode': '2016ApJS..225...10F',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [9.106e-4,-1.016e-1,4.578,-1.066e2,1.360e3,-9.183e3,2.795e4], 'range': [17.,27.],'fitunc': 198.},
    'faherty_group': {'altname': ['faherty16_group','faherty2016_group','fah16grp'], 'reference': 'Faherty et al. (2016)','bibcode': '2016ApJS..225...10F',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [7.383e0,-3.44522e2,4.87986e3], 'range': [17.,27.],'fitunc': 172.},
    'filippazzo': {'altname': ['filippazzo15','filippazzo2015','fil15'], 'reference': 'Filippazzo et al. (2015)','bibcode': '2015ApJ...810..158F',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [1.546e-4, -1.606e-2, 6.318e-1, -1.191e1, 1.155e2, -7.005e2, 4.747e3], 'range': [16.,39.],'fitunc': 113.},
    'golimowski': {'altname': ['golimowski04','golimowski2004','gol04'], 'reference': 'Golimowski et al. (2004)','bibcode': '2004AJ....127.3516G',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [9.5373e-4,-9.8598e-2,4.0323,-8.3099e1,9.0951e2,-5.1287e3,1.4322e4],'range': [16.,38.],'fitunc': 124.},
    'gonzales2018': {'altname': ['gonzales','gonzales18','gonzalez2018','gonzalez18','gonzalez','gon18','subdwarf'],'reference': 'Gonzales et al. (2018)','bibcode': '2018ApJ...864..100G', 
        'method': 'polynomial', 'sptoffset': 10., 'fitunc': 108, 'range' : [17., 27.], 'coeff': [-117,3721.]},
    'kirkpatrick2019': {'altname': ['kirkpatrick','kirkpatrick19','kirk19','kir19'], 'reference': 'Kirkpatrick et al. (2019)','bibcode': '2019ApJS..240...19K',
        'method': 'polynomial', 'sptoffset': 30, 'coeff': [9.89240,-286.401,2335.64],'range': [36.,44.],'fitunc': 68.},
    'kirkpatrick2021_l': {'altname': ['kirkpatrick_l','kirkpatrick21_l','kirk21_l','kir21_l'], 'reference': 'Kirkpatrick et al. (2021)','bibcode': '2021ApJS..253....7K',
        'method': 'polynomial', 'sptoffset': 20, 'coeff': [4.0301e0,-1.4496e2,2.2375e3],'range': [20,28.75],'fitunc': 134.},
    'kirkpatrick2021_lt': {'altname': ['kirkpatrick_lt','kirkpatrick21_lt','kirk21_lt','kir21_lt'], 'reference': 'Kirkpatrick et al. (2021)','bibcode': '2021ApJS..253....7K',
        'method': 'polynomial', 'sptoffset': 20, 'coeff': [-1.8309e1,1.4379e3],'range': [28.75,34.75],'fitunc': 79.},
    'kirkpatrick2021_ty': {'altname': ['kirkpatrick_ty','kirkpatrick21_ty','kirk21_ty','kir21_ty'], 'reference': 'Kirkpatrick et al. (2021)','bibcode': '2021ApJS..253....7K',
        'method': 'polynomial', 'sptoffset': 20, 'coeff': [6.3701e0,-3.6865e2,5.1413e3],'range': [34.75,42],'fitunc': 79.},
    'looper': {'altname': ['looper08','looper2008','lop08'], 'reference': 'Looper et al. (2008)','bibcode': '2008ApJ...685.1183L',
        'method': 'polynomial', 'sptoffset': 20.,'coeff': [9.084e-4,-4.255e-2,6.414e-1,-3.101,1.950,-108.094,2319.92],'range': [20.,38.],'fitunc': 87.},
    'pecaut': {'altname': ['mamajek','pecaut13','pecaut2013'], 'reference': 'Pecaut & Mamajek (2013)','bibcode':'2013ApJS..208....9P', 'url': 'http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt',
        'method': 'interpolate','range': [0.,29.],'fitunc': 108.,
        'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,34.5,35.,35.5,36.,37.,37.5,38.,38.5,39.,39.5,40.,40.5,41.,41.5,42.], \
        'values': [5280.,5240.,5170.,5140.,5040.,4990.,4830.,4700.,4600.,4540.,4410.,4330.,4230.,4190.,4070.,4000.,3940.,3870.,3800.,3700.,3650.,3550.,3500.,3410.,3250.,3200.,3100.,3030.,3000.,2850.,2710.,2650.,2600.,2500.,2440.,2400.,2320.,2250.,2100.,1960.,1830.,1700.,1590.,1490.,1410.,1350.,1300.,1260.,1230.,1200.,1160.,1120.,1090.,1050.,1010.,960.,840.,770.,700.,610.,530.,475.,420.,390.,350.,325.,250.],\
        'rms': numpy.zeros(67)+108.},
    'marocco': {'altname': ['marocco13','marocco2013','mar13'], 'reference': 'Marocco et al. (2013)','bibcode': '2013AJ....146..161M',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [7.4211e-5,-8.43736e-3,3.90319e-1,-9.46896,129.141,-975.953,3561.47,-1613.82],'range': [17.,38.],'fitunc': 140.},
    'sanghi2023': {'altname': ['sanghi23','san23','sanghi'], 'reference': 'Sanghi et al. (2023)','bibcode': '2023ApJ...959...63S',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [-5.75357642e-02,3.71753828e+00,-7.99773680e+01,5.52980543e+02,1.58703841e+03],'range': [16.,39.],'fitunc': 175.},
    'sanghi2023_opt': {'altname': ['sanghi23_opt','san23_opt','sanghi_opt'], 'reference': 'Sanghi et al. (2023)','bibcode': '2023ApJ...959...63S',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [-7.58233882e-02,5.08274865e+00,-1.15315067e+02,9.22014212e+02,2.89211204e+02],'range': [16.,39.],'fitunc': 183.},
    'sanghi2023_nir': {'altname': ['sanghi23_nir','san23_nir','sanghi_nir'], 'reference': 'Sanghi et al. (2023)','bibcode': '2023ApJ...959...63S',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [-4.05522698e-02,2.54129573e+00,-5.16856405e+01,2.80510088e+02,2.44117451e+03],'range': [16.,39.],'fitunc': 180.},
    'stephens': {'altname': ['stephens09','stephens2009','ste09'], 'reference': 'Stephens et al. (2009)','bibcode': '2009ApJ...702..154S',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [-0.0025492,0.17667,-4.4727,54.67,-467.26,4400.],'range': [16.,38.],'fitunc': 100.},
    'stephens_alt': {'altname': ['stephens09_alt','stephens2009_alt','ste09alt'], 'reference': 'Stephens et al. (2009)','bibcode': '2009ApJ...702..154S',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [-0.011997,1.2315,-50.472,1031.9,-10560.,44898.],'range': [23.,38.],'fitunc': 100.},
    'zhang2017': {'altname': ['zhang','zha17'], 'reference': 'Zhang et al. (2017)','bibcode': '2017MNRAS.464.3040Z',
        'method': 'polynomial', 'sptoffset': 10.,'coeff': [-0.1606,1.686,-107.8,3706],'range': [15.5,27],'fitunc': 32.5},
}


# Empirical relations - SpT to Color
SPT_COLORS_RELATIONS = {
    'best2018': {'altname': ['best','best18'], 'reference': 'Best et al. (2018)','bibcode':'2018ApJS..234....1B', 'sptoffset': 0, 'method': 'interpolate', 'scatter': 0.05, 'filters': {'PANSTARRS_G','PANSTARRS_R','PANSTARRS_I','PANSTARRS_Z','PANSTARRS_Y','2MASS_J','2MASS_H','2MASS_KS','WISE_W1','WISE_W2','WISE_W3'},
        'colors': {
            'PANSTARRS_G-PANSTARRS_R': {
                'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23], \
                'values': [1.19,1.22,1.21,1.21,1.23,1.31,1.33,1.40,1.53,1.79,1.85,2.00,2.30,2.70]},
            'PANSTARRS_R-PANSTARRS_I': {
                'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28], \
                'values': [0.67,0.85,1.02,1.22,1.46,1.88,2.13,2.55,2.69,2.58,2.35,2.35,2.27,2.23,2.11,2.00,1.91,2.13]},
            'PANSTARRS_I-PANSTARRS_Z': {
                'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35], \
                'values': [0.31,0.39,0.46,0.55,0.67,0.87,0.98,1.21,1.38,1.44,1.47,1.48,1.45,1.49,1.63,1.77,1.89,2.10,2.45,2.45,2.38,2.30,3.53,2.85,3.49]},
            'PANSTARRS_Z-PANSTARRS_Y': {
                'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
                'values': [0.17,0.20,0.23,0.27,0.32,0.44,0.51,0.67,0.81,0.92,0.95,0.97,0.97,1.01,1.02,1.05,1.05,1.02,1.05,1.16,1.21,1.39,1.45,1.55,1.71,1.74,1.74,1.80,1.76,1.37]},
            'PANSTARRS_Y-2MASS_J': {
                'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
                'values': [1.12,1.14,1.16,1.20,1.25,1.34,1.40,1.54,1.66,1.77,1.82,1.94,2.00,2.12,2.15,2.18,2.09,2.19,2.14,2.11,2.21,2.19,2.22,2.24,2.43,2.49,2.54,2.56,2.57,2.61]},
            '2MASS_J-2MASS_H': {
                'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
                'values': [0.66,0.64,0.62,0.60,0.59,0.59,0.60,0.63,0.68,0.72,0.76,0.80,0.91,0.95,1.06,1.06,0.99,1.12,1.13,1.07,0.93,0.86,0.80,0.59,0.32,0.18,0.05,0.08,0.09,0.24]},
            '2MASS_H-2MASS_KS': {
                'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37], \
                'values': [0.18,0.21,0.22,0.24,0.26,0.31,0.33,0.39,0.43,0.48,0.48,0.51,0.59,0.64,0.64,0.67,0.65,0.73,0.65,0.58,0.56,0.33,0.25,0.17,0.11,0.03,0.06,0.01]},
            '2MASS_KS-WISE_W1': {
                'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37], \
                'values': [0.10,0.11,0.12,0.14,0.16,0.19,0.20,0.22,0.26,0.31,0.32,0.35,0.42,0.52,0.61,0.67,0.71,0.73,0.85,0.78,0.78,0.63,0.57,0.56,0.19,0.15,0.37,0.33]},
            'WISE_W1-WISE_W2': {
                'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
                'values': [0.02,0.07,0.12,0.16,0.18,0.21,0.22,0.23,0.23,0.26,0.27,0.26,0.29,0.30,0.32,0.33,0.38,0.42,0.53,0.57,0.58,0.72,0.91,1.07,1.35,1.85,2.08,2.26,2.79,2.95]},
            'WISE_W2-WISE_W3': {
                'spt': [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,36,37,38,39], \
                'values': [0.21,0.25,0.32,0.35,0.44,0.51,0.46,0.55,0.50,0.64,0.65,0.76,0.96,1.12,0.88,0.92,1.41,1.44,1.30,1.21,1.35,1.40,1.82]},
            },
        },
    'kirkpatrick2019': {'altname': ['kirkpatrick19','kirk19','kir19'], 'bibcode': '2019ApJS..240...19K', 'sptoffset': 30, 'method': 'polynomial', 
        'colors': {
            '2MASS_H-WISE_W2': {'fitunc': 0.51, 'range': [36,44], 'coeff': [-0.0264937,0.806149,-6.76528,19.7444]},\
            'WISE_W1-WISE_W2': {'fitunc': 0.39, 'range': [36,44], 'coeff': [0.0290939,-0.0609403,0.812400,-1.29441]},\
            '2MASS_H-IRAC_CH1': {'fitunc': 0.51, 'range': [36,44], 'coeff': [-0.0291438,0.876955,-7.36894,21.3744]},\
            'IRAC_CH1-IRAC_CH2': {'fitunc': 0.19, 'range': [36,44], 'coeff': [-0.00245289,0.0645336,-0.229115,0.769877]},\
            },
        },
    'kirkpatrick2021': {'altname': ['kirkpatrick','kirkpatrick121','kirk21','kir21'], 'bibcode': '2020arXiv201111616K', 'sptoffset': 20, 'method': 'polynomial', 
        'colors': {
            'MKO_J-IRAC_CH2': {'fitunc': 0.44, 'range': [30,42.], 'coeff': [3.6779e-7,-2.7663e-5,8.1335e-4,-1.1351e-2,7.2989e-2,-1.7678e-1,1.8527e-1,1.8153e0]},\
            '2MASS_H-IRAC_CH2': {'fitunc': 0.42, 'range': [30,42.], 'coeff': [3.211e-7,-2.1740e-5,5.7149e-4,-7.1625e-3,4.0283e-2,-7.3996e-2,6.2704e-2,1.1150e0]},\
            'IRAC_CH1-IRAC_CH2': {'fitunc': 0.19, 'range': [30,42.], 'coeff': [3.2520e-4,8.1897e-4,-2.6015e-2,2.6662e-2]},\
            'WISE_W1-WISE_W2': {'fitunc': 0.28, 'range': [30,42.], 'coeff': [5.7825e-4,-4.6379e-3,2.9069e-2,2.2668e-2]},\
            },
        },
    'leggett2017': {'altname': ['leggett','ydwarf','leggett17'], 'reference': 'Leggett et al. (2017)','bibcode': '2017ApJ...842..118L','sptoffset': 30,'method': 'polynomial', 'range' : [37.5,42], 'fitunc': 0.1, 
        'colors': {
            'MKO_Y-WFC3_F105W': {'fitunc' : 0.1, 'range' : [37.5,42], 'coeff': [0.027813,-0.601128,2.42719]}, \
            'MKO_J-WFC3_F125W': {'fitunc' : 0.1, 'range' : [37.5,42], 'coeff': [0.020710,-0.398120,1.228187]}, \
            'MKO_J-WFC3_F127M': {'fitunc' : 0.1, 'range' : [37.5,42], 'coeff': [0.014005,-0.1183161,1.092170]}, \
            'MKO_H-WFC3_F160W': {'fitunc' : 0.1, 'range' : [37.5,42], 'coeff': [0.001134,-0.008060,-0.21659]}, \
            'MKO_H-WIRC_CH4S': {'fitunc' : 0.1, 'range' : [37.5,42], 'coeff': [-0.006130,0.204865,-0.679045]}
            },
        },
    'pecaut2013': {'altname': ['mamajek','pecaut','pecaut13'], 'reference': 'Pecaut & Mamajek (2013)','bibcode':'2013ApJS..208....9P', 'url': 'http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt','sptoffset': 0, 'method': 'interpolate', 'scatter': 0.05, 'filters': {'BESSEL_U','BESSEL_B','BESSEL_V','BESSEL_R','BESSEL_I','GAIA_G','SDSS_I','SDSS_Z','UKIDSS_Y','2MASS_J','2MASS_H','2MASS_KS','WISE_W1','WISE_W2','WISE_W3','WISE_W4'},
        'colors': {
            'BESSEL_U-BESSEL_B': {
                'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16], \
                'values': [0.436, 0.456, 0.491, 0.528, 0.581, 0.691, 0.776, 0.917, 1.004, 1.028, 1.081, 1.144, 1.184, 1.213, 1.221, 1.216, 1.21, 1.204, 1.184, 1.172, 1.17, 1.17, 1.175, 1.181, 1.2, 1.222, 1.23, 1.24, 1.3, 1.3]},
            'BESSEL_B-BESSEL_V': {
                'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18], \
                'values': [0.816, 0.825, 0.842, 0.86, 0.884, 0.938, 0.986, 1.05, 1.1, 1.116, 1.15, 1.2, 1.24, 1.31, 1.33, 1.38, 1.42, 1.408, 1.441, 1.475, 1.486, 1.5, 1.522, 1.544, 1.602, 1.661, 1.72, 1.874, 1.91, 2, 2.06, 2.06, 2.17, 2.2]},
            'GAIA_G-BESSEL_V': {
                'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19.,19.5,20.], \
                'values': [-0.221, -0.225, -0.232, -0.241, -0.254, -0.288, -0.322, -0.371, -0.412, -0.425, -0.454, -0.495, -0.528, -0.581, -0.595, -0.628, -0.69, -0.65, -0.74, -0.82, -0.85, -0.92, -1.02, -1.09, -1.29, -1.41, -1.55, -1.74, -2.01, -2.14, -2.64, -2.98, -3.03, -3.08, -3.04, -3, -3.1, -3.2]},
            'BESSEL_V-BESSEL_R': {
                'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19.,19.5], \
                'values': [0.443, 0.448, 0.457, 0.467, 0.482, 0.513, 0.537, 0.592, 0.64, 0.654, 0.685, 0.728, 0.759, 0.796, 0.806, 0.843, 0.866, 0.889, 0.924, 0.959, 0.978, 1.001, 1.041, 1.079, 1.178, 1.241, 1.345, 1.446, 1.656, 1.95, 2.003, 2.18, 2.16, 2.15, 1.967, 1.89, 2.51]},
            'BESSEL_V-BESSEL_I': {
                'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19.,19.5,20.,21.,22.,23.,24.,25.], \
                'values': [0.853, 0.862, 0.879, 0.895, 0.92, 0.974, 1.013, 1.108, 1.19, 1.216, 1.272, 1.357, 1.42, 1.505, 1.529, 1.632, 1.699, 1.766, 1.886, 2.019, 2.089, 2.173, 2.306, 2.42, 2.68, 2.831, 3.073, 3.277, 3.664, 4.13, 4.31, 4.45, 4.56, 4.64, 4.71, 4.75, 4.79, 4.82, 4.91, 5.05, 5.29, 5.57, 6.28]},
            'BESSEL_V-2MASS_KS': {
                'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19.,19.5,20.,21.,22.,23.,24.,25.], \
                'values': [1.953, 1.977, 2.021, 2.066, 2.132, 2.274, 2.38, 2.582, 2.733, 2.781, 2.883, 3.034, 3.143, 3.288, 3.33, 3.487, 3.584, 3.68, 3.84, 4.02, 4.12, 4.24, 4.43, 4.6, 5, 5.25, 5.64, 5.94, 6.5, 7.3, 7.6, 8.05, 8.45, 8.73, 8.92, 9, 9.3, 9.45, 9.7, 10, 10.4, 10.9, 11.4]},
            '2MASS_J-2MASS_H': {
                'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,34.5,35.,35.5,36.,37.,37.5,38.,38.5,39.], \
                'values': [0.387, 0.393, 0.402, 0.412, 0.427, 0.459, 0.483, 0.52, 0.544, 0.552, 0.568, 0.591, 0.601, 0.613, 0.617, 0.623, 0.625, 0.626, 0.62, 0.613, 0.607, 0.6, 0.589, 0.579, 0.558, 0.557, 0.564, 0.58, 0.588, 0.605, 0.609, 0.613, 0.65, 0.67, 0.685, 0.749, 0.826, 0.79, 0.8, 0.87, 1, 1.14, 1.13, 1.08, 1.1, 1.14, 1.1, 1.02, 1.02, 0.86, 0.68, 0.35, 0.2, 0.2, 0.2, 0.1, 0, 0.2, 0.2, 0.2, 0.1]},
            '2MASS_H-2MASS_KS': {
                'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,34.5,35.,35.5,36.,37.,37.5,38.,39.,40.,40.5,41.5], \
                'values': [0.091, 0.092, 0.094, 0.096, 0.098, 0.104, 0.109, 0.118, 0.125, 0.127, 0.132, 0.139, 0.148, 0.159, 0.162, 0.176, 0.184, 0.193, 0.208, 0.225, 0.228, 0.234, 0.244, 0.252, 0.269, 0.282, 0.301, 0.311, 0.329, 0.352, 0.364, 0.386, 0.422, 0.45, 0.47, 0.48, 0.505, 0.5, 0.54, 0.57, 0.63, 0.63, 0.65, 0.64, 0.62, 0.63, 0.63, 0.54, 0.45, 0.27, 0.08, -0.19, -0.06, -0.08, -0.1, -0.03, 0, -0.05, -0.05,-0.2,-0.5,-0.6,-0.8]},
            'WISE_W1-WISE_W2': {
                'spt': [15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,34.5,35.,35.5,36.,37.,37.5,38.,40.,40.5,41.,42.], \
                'values': [0.17, 0.19, 0.21, 0.22, 0.24, 0.25, 0.26, 0.265, 0.27, 0.27, 0.27, 0.28, 0.28, 0.29, 0.3, 0.32, 0.36, 0.41, 0.48, 0.57, 0.68, 0.82, 0.99, 1.19, 1.43, 1.57, 1.7, 1.86, 2.02, 2.38, 2.59, 2.79, 2.4, 2.5, 2.6, 3.4]},
            'SDSS_I-SDSS_Z': {
                'spt': [10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,34.5,35.,35.5,36.,37.,37.5,38.], \
                'values': [0.33, 0.37, 0.41, 0.47, 0.53, 0.57, 0.61, 0.66, 0.71, 0.81, 0.91, 1.13, 1.45, 1.58, 1.77, 1.85, 1.93, 1.96, 1.99, 2, 2.01, 2.02, 2.04, 2.1, 2.2, 2.33, 2.51, 2.71, 2.93, 3.15, 3.36, 3.55, 3.7, 3.82, 3.9, 3.93, 3.95, 3.97, 3.98, 4.01, 4.05, 4.08]},
            'SDSS_Z-UKIDSS_Y': {
                'spt': [15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,34.5,35.,35.5,36.,37.,37.5,38.], \
                'values': [0.47, 0.52, 0.6, 0.64, 0.7, 0.74, 0.77, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.97, 1, 1.04, 1.09, 1.16, 1.23, 1.33, 1.43, 1.55, 1.68, 1.75, 1.81, 1.89, 1.96, 2.11, 2.19, 2.26]},
            },
        },
    'schmidt2015': {'altname': ['schmidt','schmidt15'], 'reference': 'Schmidt et al. (2015)','bibcode':'2015AJ....149..158S', 'sptoffset': 0, 'method': 'interpolate', 'scatter': 0.1, 'filters': {'SDSS_I','SDSS_Z','2MASS_J','2MASS_H','2MASS_KS','WISE_W1','WISE_W2','WISE_W3'},
        'colors': {
            'SDSS_I-SDSS_Z': {
                'spt': [17,18,19,20,21,22,23,25], \
                'values': [1.17,1.48,1.66,1.82,1.81,1.86,2.13]},
            'SDSS_I-2MASS_J': {
                'spt': [17,18,19,20,21,22,23,25], \
                'values': [2.89,3.41,3.79,4.22,4.41,4.45,4.51,4.88]},
            'SDSS_I-2MASS_J': {
                'spt': [17,18,19,20,21,22,23,25], \
                'values': [3.86,4.42,4.91,5.44,5.73,5.83,5.95,6.4]},
            'SDSS_Z-2MASS_J': {
                'spt': [17,18,19,20,21,22,23,25,26], \
                'values': [1.72,1.95,2.14,2.37,2.5,2.59,2.67,2.82,2.76]},
            '2MASS_J-2MASS_H': {
                'spt': [17,18,19,20,21,22,24,23,25,26], \
                'values': [0.61,0.64,0.68,0.74,0.80,0.89,0.91,0.91,0.95,0.96]},
            '2MASS_J-2MASS_KS': {
                'spt': [17,18,19,20,21,22,24,23,25,26], \
                'values': [0.96,1.03,1.12,1.20,1.31,1.45,1.52,1.47,1.53,1.54]},
            '2MASS_H-2MASS_KS': {
                'spt': [17,18,19,20,21,22,24,23,25,26], \
                'values': [0.34,0.39,0.43,0.46,0.52,0.56,0.63,0.58,0.60,0.59]},
            '2MASS_KS-WISE_W1': {
                'spt': [17,18,19,20,21,22,24,23,25,26], \
                'values': [0.18,0.20,0.25,0.32,0.37,0.44,0.41,0.55,0.61,0.75]},
            'WISE_W1-WISE_W2': {
                'spt': [17,18,19,20,21,22,24,23,25,26], \
                'values': [0.2,0.22,0.24,0.27,0.26,0.27,0.31,0.28,0.30,0.34]},
            'WISE_W2-WISE_W3': {
                'spt': [17,18,19,20,23], \
                'values': [0.29,0.34,0.42,0.44,0.78]},
            },
        },
    'skrzypek2015': {'altname': ['skrzypek','skrzypek15'], 'reference': 'Skrzypek et al. (2015)','bibcode': '2015A%26A...574A..78S','range': [15,38], 'method': 'interpolate', 'filters': ['i','z','y','j','h','k','w1','w2'],'scatter': 0.07,
        'colors': { 
            'SDSS_I-SDSS_Z': {'fitunc' : 0.07, 'spt' : numpy.arange(15,38.1), 'values': [0.91,1.45,1.77,1.93,1.99,2.01,2.02,2.04,2.1,2.2,2.33,2.51,2.71,2.93,3.15,3.36,3.55,3.7,3.82,3.9,3.95,3.98,4.01,4.08]}, \
            'SDSS_Z-UKIDSS_Y': {'fitunc' : 0.07, 'spt' : numpy.arange(15,38.1), 'values': [0.47,0.6,0.7,0.77,0.82,0.86,0.88,0.9,0.92,0.94,0.97,1.0,1.04,1.09,1.16,1.23,1.33,1.43,1.55,1.68,1.81,1.96,2.11,2.26]}, \
            'UKIDSS_Y-MKO_J': {'fitunc' : 0.07, 'spt' : numpy.arange(15,38.1), 'values': [0.55,0.67,0.78,0.87,0.96,1.04,1.11,1.18,1.23,1.27,1.31,1.33,1.35,1.21,1.2,1.19,1.19,1.18,1.18,1.17,1.16,1.16,1.15,1.15]}, \
            'MKO_J-MKO_H': {'fitunc' : 0.07, 'spt' : numpy.arange(15,38.1), 'values': [0.45,0.53,0.56,0.58,0.6,0.63,0.67,0.73,0.79,0.86,0.91,0.96,0.97,0.96,0.9,0.8,0.65,0.46,0.25,0.02,-0.19,-0.35,-0.43,-0.36]}, \
            'MKO_H-MKO_K': {'fitunc' : 0.07, 'spt' : numpy.arange(15,38.1), 'values': [0.32,0.39,0.44,0.47,0.51,0.54,0.58,0.63,0.67,0.71,0.74,0.75,0.75,0.71,0.65,0.56,0.45,0.31,0.16,0.01,-0.11,-0.19,-0.2,-0.09]}, \
            'MKO_K-WISE_W1': {'fitunc' : 0.07, 'spt' : numpy.arange(15,38.1), 'values': [0.11,0.22,0.25,0.26,0.27,0.29,0.33,0.4,0.48,0.56,0.65,0.72,0.77,0.79,0.79,0.76,0.71,0.65,0.59,0.55,0.54,0.59,0.7,0.9]}, \
            'WISE_W1-WISE_W2': {'fitunc' : 0.07, 'spt' : numpy.arange(15,38.1), 'values': [0.17,0.21,0.24,0.26,0.27,0.27,0.28,0.28,0.29,0.3,0.32,0.36,0.41,0.48,0.57,0.68,0.82,0.99,1.19,1.43,1.7,2.02,2.38,2.79]}
            },
        },
    'theissen2017': {'altname': ['theissen','theissen17','the17'], 'reference': 'Theissen et al. (2017)','bibcode':'2017AJ....153...92T', 'sptoffset': 0, 'method': 'interpolate', 'scatter': 0.1, 'filters': {'SDSS_I','SDSS_Z'},
        'colors': {
            'SDSS_I-SDSS_Z': {
                'spt': [15,16,17,18,19,20,21,22,23,24,25,26,27,28], \
                'values': [0.86,1.06,1.17,1.48,1.66,1.85,1.86,1.85,1.88,2.13,2.16,2.43,2.59,2.81]},
            },
        },
    }


# Empirical relations - SpT to BC
SPT_BC_RELATIONS = {
    'dupuy2013': {'altname': ['dupuy13','dupuy'], 'bibcode': '2013Sci...341.1492D', 'sptoffset': 10, 'method': 'interpolate', 'scatter': 0.6, 'filters': {
        'MKO_Y': {'spt': [38,38.5,39,40.,40.5], 'values': [1.6,1.1,1.0,0.8,-0.7], 'rms': [0.6,0.6,0.6,0.6,0.6]}, 
        'MKO_J': {'spt': [38,38.5,39,39.5,40.,40.5], 'values': [2.4,2.0,2.0,1.9,0.9,-0.4], 'rms': [0.6,0.6,0.6,0.6,0.6,0.6]}, 
        'MKO_H': {'spt': [38,38.5,39,39.5,40.,40.5], 'values': [2.1,1.6,1.7,1.5,0.4,-1.5], 'rms': [0.6,0.6,0.6,0.6,0.6,0.6]},
        }}, 
    'filippazzo2015': {'altname': ['filippazzo','filippazzo-field','filippazzo15','filippazzo15-field','fillippazzo','fillippazzo-field','filipazo','filipazo-field','filippazo','filippazo-field'], 'bibcode': '2015ApJ...810..158F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.163, 'range' : [16,39], 'coeff': [-1.555e-5,1.200e-3,-3.437e-2,4.566e-1,-2.862e0,8.842e0]},
        '2MASS_KS': {'fitunc' : 0.243, 'range' : [16,38], 'coeff': [-5.474e-6,4.549e-4,-1.462e-2,2.204e-1,-1.508e0,6.815e0]},
        }},
    'filippazzo2015-young': {'altname': ['filippazzo-young','filippazzo15-young','fillippazzo-young','filipazo-young','filippazo-young'], 'bibcode': '2015ApJ...810..158F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.189, 'range' : [17,38], 'coeff': [-1.283e-6,-1.118e-4,1.165e-2,-3.020e-1,2.899e0,-7.352e0]},
        '2MASS_KS': {'fitunc' : 0.126, 'range' : [17,38], 'coeff': [2.742e-6,-2.633e-4,9.711e-3,-1.759e-1,1.565e0,-2.174e0]},
        }},
    'liu2010': {'altname': ['liu10','liu'], 'bibcode': '2010ApJ...722..311L', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.14, 'range' : [16,38.5], 'coeff': [1.462266e-6,-1.558986e-4,6.540717e-3,-1.367605e-1,1.491738e0,-8.053993e0,1.890448e1]},
        'MKO_H': {'fitunc' : 0.07, 'range' : [16,38.5], 'coeff': [1.148133e-6,-1.171595e-4,4.733874e-3,-9.618535e-2,1.027185e0,-5.426683e0,1.366709e1]},
        'MKO_K': {'fitunc' : 0.08, 'range' : [16,38.5], 'coeff': [3.159780e-7,-3.177629e-5,1.282657e-3,-2.647464e-2,2.868014e-1,-1.471358e0,5.795848e0]},
        }},
    'schmidt2014': {'altname': ['schmidt','schmidt14'], 'reference': 'Schmidt et al. (2014)','bibcode':'2014PASP..126..642S', 'sptoffset': 0, 'method': 'interpolate', 'scatter': 0.1, 'filters': {
        'SDSS_Z': {
            'spt': [17,18,19,20,21,22,23,24,25,26,27,28], \
            'values': [0.18,-0.08,-0.30,-0.57,-0.58,-0.72,-0.79,-0.64,-1.03,-1.06,-1.52,-1.03]},
        '2MASS_J': {
            'spt': [17,18,19,20,21,22,23,24,25,26,27,28], \
            'values': [1.97,1.97,1.96,1.96,1.94,1.85,1.81,1.97,1.58,1.70,1.24,1.70]},
        '2MASS_KS': {
            'spt': [17,18,19,20,21,22,23,24,25,26,27,28], \
            'values': [2.92,3.00,3.13,3.19,3.24,3.36,3.28,3.36,3.38,3.31,3.31,3.35]},
        }},
    'pecaut2013': {'altname': ['mamajek','pecaut','pecaut13'], 'reference': 'Pecaut & Mamajek (2013)','bibcode':'2013ApJS..208....9P', 'url': 'http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt','sptoffset': 0, 'method': 'interpolate', 'scatter': 0.05, 'filters': {
        'BESSEL_V': {
            'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.], \
            'values': [-0.23, -0.26, -0.27, -0.29, -0.32, -0.41, -0.45, -0.56, -0.6, -0.68, -0.75, -0.81, -0.92, -0.95, -1.03, -1.1, -1.16, -1.38, -1.44, -1.57, -1.65, -1.76, -1.97, -2.27, -2.59, -3.05, -3.28, -3.8, -4.36, -4.6, -5.06, -5.46, -5.7, -5.8, -5.9, -6.13, -6.33, -6.51, -6.56, -6.74, -7.28, -7.76]},
        }},
    'sanghi2023': {'altname': ['sanghi23','san23','sanghi'], 'bibcode': '2023ApJ...959...63S', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'SDSS_Z': {'fitunc' : 0.208, 'range' : [16,39], 'coeff': [-3.67327336e-05,2.00785031e-03,-3.01625334e-02,-3.28812471e-02,8.82748312e-01]},
        '2MASS_J': {'fitunc' : 0.157, 'range' : [16,39], 'coeff': [-4.49787817e-06,2.91250402e-04,-6.23850138e-03,4.91437247e-02,-1.00428244e-01,1.83270671e+00]},
        'MKO_J': {'fitunc' : 0.141, 'range' : [16,39], 'coeff': [-7.12860976e-06,5.03938151e-04,-1.27017514e-02,1.41154739e-01,-7.05942778e-01,3.34867561e+00]},
        '2MASS_H': {'fitunc' : 0.114, 'range' : [16,39], 'coeff': [-1.63230829e-05,9.23336345e-04,-1.90813876e-02,1.78079042e-01,2.06836504e+00]},
        'MKO_H': {'fitunc' : 0.091, 'range' : [16,39], 'coeff': [-2.89674748e-05,1.68671698e-03,-3.51960013e-02,3.16405418e-01,1.60830399e+00]},
        '2MASS_KS': {'fitunc' : 0.169, 'range' : [16,39], 'coeff': [-3.28963542e-04,9.29897400e-03,-3.73942004e-02,2.96143878e+00]},
        'MKO_K': {'fitunc' : 0.148, 'range' : [16,39], 'coeff': [3.72518239e-05,-2.73558092e-03,6.31017798e-02,-5.30973644e-01,4.56079546e+00]},
        'WISE_W1': {'fitunc' : 0.193, 'range' : [16,39], 'coeff': [-2.74443149e-07,3.99592136e-05,-2.05225834e-03,4.92274677e-02,-5.96094136e-01,3.60659733e+00,-5.36228544e+00]},
        'WISE_W2': {'fitunc' : 0.198, 'range' : [16,39], 'coeff': [-8.53309439e-07,9.68864595e-05,-4.27523017e-03,9.35831952e-02,-1.07147549e+00,6.19843221e+00,-1.07410524e+01]},
        }},
    }

# Empirical relations - SpT to Lbol
SPT_LBOL_RELATIONS = {
    'filippazzo2015': {'altname': ['filippazzo','filippazzo15','fillippazzo','filipazo','filippazo'], 'bibcode': '2015ApJ...810..158F', 
        'sptoffset': 10., 'method': 'polynomial', 'fitunc' : 0.133, 'range' : [16,39], 'coeff': [2.736e-7,-3.220e-5,1.449e-3,-3.207e-2,3.727e-1,-2.310e0,2.787e0]},
    'gonzales2018': {'altname': ['gonzales','gonzales18','gonzalez2018','gonzalez18','gonzalez','gon18','subdwarf'],'bibcode': '2018ApJ...864..100G', 
        'sptoffset': 10., 'method': 'polynomial', 'fitunc': 0.144, 'range' : [17., 27.], 'coeff': [-0.125,-2.10]},
    'pecaut2013': {'altname': ['mamajek','pecaut','pecaut13'], 'reference': 'Pecaut & Mamajek (2013)','bibcode':'2013ApJS..208....9P', 'url': 'http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt','sptoffset': 0, 'method': 'interpolate', 'scatter': 0.05,
        'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,34.5,35.,35.5,36.,37.,37.5,38.,38.5,39.,40.],
        'values': [-0.33, -0.36, -0.39, -0.47, -0.51, -0.58, -0.64, -0.67, -0.68, -0.78, -0.84, -0.9, -0.96, -0.98, -1.1, -1.18, -1.2, -1.27, -1.4, -1.47, -1.57, -1.68, -1.78, -2.07, -2.2, -2.31, -2.52, -2.79, -3.02, -3.09, -3.21, -3.29, -3.36, -3.44, -3.57, -3.55, -3.57, -3.7, -3.84, -3.98, -4.11, -4.24, -4.36, -4.46, -4.55, -4.61, -4.66, -4.69, -4.73, -4.77, -4.84, -4.9, -4.95, -5.04, -5.12, -5.37, -5.54, -5.71, -5.93, -6.15, -6.52]
        },
    'sanghi2023': {'altname': ['sanghi23','san23','sanghi'],'bibcode': '2023ApJ...959...63S',
        'sptoffset': 10., 'method': 'polynomial', 'fitunc': 0.157, 'range' : [16., 39.], 'coeff': [3.85072198e-07,-4.22643837e-05,1.79418353e-03,-3.75870300e-02,4.10403429e-01,-2.34765187e+00,2.38575339e+00]},
    }

# Empirical relations - Absolute magnitude to Lbol
ABSMAG_LBOL_RELATIONS = {
    'dupuy2017': {'altname': ['dupuy17','dupuy'], 'bibcode': '2017ApJS..231...15D',  'method': 'polynomial', 'filters': {
        'MKO_H': {'fitunc' : 0.023, 'range' : [9.6,13.3], 'coeff': [1.06200e-2,-3.51721e-1,3.46876e1,-1.3282e1]},
        'MKO_K': {'fitunc' : 0.05, 'range' : [9.1,17.8], 'coeff': [4.54547e-4,-3.068824e-2,8.162709e-1,-1.0671188e1,6.811147e1,-1.72188e2]},
        '2MASS_H': {'fitunc' : 0.023, 'range' : [9.2,13.3], 'coeff': [8.9222e-3,-2.90797e-1,2.74259e0,-1.0426e1]},
        '2MASS_KS': {'fitunc' : 0.05, 'range' : [8.8,16.6], 'coeff': [9.8821e-5,-7.85837e-3,2.357643e-1,-3.364101e0,2.258776e1,-5.9877e1]}},
    },
}

# Empirical relations - SpT to Absolute magnitude
SPT_ABSMAG_RELATIONS = {
    'best2018': {'altname': ['best','best18','bes18'], 'reference': 'Best et al. (2018)','bibcode':'2018ApJS..234....1B', 'sptoffset': 0, 'method': 'interpolate', 'scatter': 0.05, 'filters': {
        'PANSTARRS_G': {\
            'spt': [16,17,18,19,20,21,22,23], \
            'values': [16.71,18.11,19.19,19.95,20.33,20.85,21.24,22.51],\
            'rms': [0.44,0.56,0.59,0.43,0.40,0.78,0.5,0.5]},
        'PANSTARRS_R': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,28], \
            'values': [15.37,16.76,17.74,18.14,18.37,18.74,19.02,19.61,20.60,20.74,21.21,22.88],\
            'rms': [0.43,0.50,0.51,0.37,0.31,0.28,0.29,0.39,0.56,0.37,0.78,0.5]},
        'PANSTARRS_I': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,35], \
            'values': [13.25,14.18,15.00,15.62,16.00,16.41,16.73,17.40,18.35,18.71,19.27,20.09,20.38,20.09,20.22,21.10,21.97,22.69],\
            'rms': [0.34,0.39,0.48,0.39,0.26,0.25,0.26,0.34,0.38,0.33,0.65,0.36,0.79,0.5,1.14,0.5,0.5,0.5]},
        'PANSTARRS_Z': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
            'values': [12.27,12.98,13.61,14.17,14.52,14.93,15.30,15.88,16.68,16.94,17.35,18.18,18.10,17.69,17.98,18.84,18.26,18.08,18.02,19.20,19.82,21.17,21.52,21.82],\
            'rms': [0.32,0.34,0.45,0.37,0.25,0.23,0.24,0.21,0.33,0.29,0.62,0.26,0.22,0.73,0.5,0.21,0.23,0.25,0.39,0.22,0.32,0.78,0.52,0.5]},
        'PANSTARRS_Y': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
            'values': [11.76,12.31,12.81,13.23,13.58,13.97,14.33,14.89,15.66,15.87,16.27,17.13,17.04,16.57,16.77,17.45,16.75,16.50,16.32,17.43,18.06,19.34,19.75,20.37],\
            'rms': [0.30,0.31,0.43,0.36,0.23,0.21,0.24,0.29,0.32,0.28,0.61,0.25,0.21,0.72,0.5,0.16,0.13,0.22,0.38,0.18,0.32,0.81,0.59,0.5]},
        '2MASS_J': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
            'values': [10.36,10.77,11.15,11.46,11.76,12.03,12.32,12.77,13.51,13.69,14.18,14.94,14.90,14.46,14.56,15.25,14.54,14.26,13.89,14.94,15.53,16.78,17.18,17.75],\
            'rms': [0.30,0.30,0.42,0.34,0.18,0.15,0.21,0.24,0.28,0.25,0.60,0.20,0.13,0.71,0.5,0.12,0.06,0.16,0.36,0.12,0.27,0.76,0.51,0.5]},
        '2MASS_H': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
            'values': [9.76,10.14,10.47,10.74,11.00,11.23,11.41,11.82,12.45,12.63,13.19,13.82,13.77,13.39,13.62,14.39,13.73,13.67,13.57,14.76,15.48,16.70,17.09,17.51],\
            'rms': [0.30,0.31,0.43,0.35,0.23,0.21,0.25,0.29,0.3,0.30,0.62,0.31,0.20,0.73,0.5,0.18,0.15,0.24,0.40,0.24,0.37,0.78,0.5,0.5]},
        '2MASS_KS': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37], \
            'values': [9.43,9.75,10.03,10.27,10.52,10.71,10.83,11.19,11.85,11.94,12.58,13.12,13.04,12.79,13.08,14.07,13.48,13.44,13.49,14.77,15.37,16.70],\
            'rms': [0.31,0.31,0.43,0.37,0.24,0.24,0.29,0.32,0.49,0.34,0.68,0.32,0.22,0.74,0.5,0.33,0.26,0.32,0.39,0.23,0.41,0.93]},
        'WISE_W1': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
            'values': [9.22,9.52,9.78,9.95,10.22,10.36,10.43,10.66,11.18,11.26,11.85,12.42,12.20,11.94,12.43,13.39,12.94,12.84,13.31,14.52,15.08,16.23,16.58,16.70],\
            'rms': [0.30,0.31,0.43,0.37,0.26,0.26,0.32,0.36,0.43,0.34,0.71,0.48,0.28,0.74,0.5,0.38,0.39,0.61,0.41,0.25,0.35,0.84,0.84,0.5]},
        'WISE_W2': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39], \
            'values': [9.01,9.31,9.55,9.68,9.94,10.10,10.14,10.39,10.84,10.91,11.49,11.99,11.66,11.45,11.94,12.68,12.06,11.72,11.90,12.69,13.02,14.11,14.03,13.81],\
            'rms': [0.30,0.31,0.44,0.37,0.26,0.26,0.33,0.38,0.46,0.38,0.75,0.57,0.30,0.77,0.5,0.33,0.33,0.45,0.39,0.24,0.42,0.95,0.92,0.5]},
        'WISE_W3': {\
            'spt': [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39], \
            'values': [8.72,8.96,9.18,9.24,9.35,9.66,9.62,9.91,10.05,10.13,10.53,10.87,10.66,10.48,11.02,11.26,10.61,10.70,11.37,11.61,12.26,12.27,11.97],\
            'rms': [0.34,0.33,0.47,0.42,0.32,0.31,0.37,0.37,0.44,0.51,0.87,1.08,0.36,0.75,0.5,0.76,0.40,0.5,0.5,0.37,0.76,0.64,0.5]},
    }},
    'burgasser2007': {'altname': ['burgasser','burgasser07','bur07'], 'bibcode': '2007ApJ...659..655B', 'sptoffset': 20, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.30, 'range' : [20., 38.], 'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, 
        'MKO_H': {'fitunc' : 0.27, 'range' : [20., 38.], 'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, 
        'MKO_K': {'fitunc' : 0.26, 'range' : [20., 38.], 'coeff': [.0000001051, -.000006985, .0001807, -.002271, .01414, -.04024, .05129, .2322, 10.45]},
    }},
    'cruz2003': {'altname': ['cruz','cruz03','cru03'], 'bibcode': '2003AJ....126.2421C', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.30, 'range' : [16.,28.], 'coeff': [-6.892e-4,3.453e-2,-6.193e-1,5.043,-4.410]},
    }},
    'dahn2002': {'altname': ['dahn','dahn02','dah02'], 'bibcode': '2002AJ....124.1170D', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.25, 'range' : [17.,28.], 'coeff': [0.341,8.38]},
    }},
    'dupuy2012': {'altname': ['dupuy','dupuy12','dup12'], 'bibcode': '2012ApJS..201...19D', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc': 0.40, 'range' : [16., 39.], 'coeff': [-.00000252638, .000285027, -.0126151, .279438, -3.26895, 19.5444, -35.1560]}, 
        'MKO_J': {'fitunc' : 0.39, 'range' : [16., 39.], 'coeff' : [-.00000194920, .000227641, -.0103332, .232771, -2.74405, 16.3986, -28.3129]}, 
        'MKO_H': {'fitunc': 0.38, 'range' : [16., 39.], 'coeff': [-.00000224083, .000251601, -.0110960, .245209, -2.85705, 16.9138, -29.7306]}, 
        'MKO_K': {'fitunc': 0.40, 'range' : [16., 39.], 'coeff': [-.00000104935, .000125731, -.00584342, .135177, -1.63930, 10.1248, -15.2200]}, 
        'MKO_LP': {'fitunc': 0.28, 'range': [16., 39.], 'coeff': [0.00000, 0.00000, .0000546366, -.00293191, .0530581,  -.196584, 8.89928]}, 
        '2MASS_J': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [-.000000784614, .000100820, -.00482973, .111715, -1.33053, 8.16362, -9.67994]}, 
        '2MASS_H': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [-.00000111499, .000129363, -.00580847, .129202, -1.50370, 9.00279, -11.7526]}, 
        '2MASS_KS': {'fitunc': 0.43, 'range':[16., 39.], 'coeff': [1.06693e-4, -6.42118e-3, 1.34163e-1, -8.67471e-1, 1.10114e1]}, 
        'IRAC_CH1': {'fitunc': 0.29, 'range':[16., 39.], 'coeff': [6.50191e-5,-3.60108e-3,6.91081e-2,-3.35222e-1,9.34220e0]}, 
        'IRAC_CH2': {'fitunc': 0.22, 'range':[16., 39.], 'coeff': [5.82107e-5,-3.63435e-3,7.65343e-2,-4.39968e-1,9.73946e0]}, 
        'IRAC_CH3': {'fitunc': 0.32, 'range':[16., 39.], 'coeff': [1.03507e-4,-6.22795e-3,1.29019e-1,-9.01820e-1,1.10834e1]}, 
        'IRAC_CH4': {'fitunc': 0.27, 'range':[16., 39.], 'coeff': [6.89733e-5,-4.12294e-3,8.43465e-2,-5.29595e-1,9.97853e0]}, 
        'WISE_W1': {'fitunc': 0.39, 'range':[16., 39.], 'coeff': [1.58040e-5, -3.33944e-4, -4.38105e-3, 3.55395e-1, 7.14765]}, 
        'WISE_W2': {'fitunc': 0.35, 'range':[16., 39.], 'coeff': [1.78555e-5, -8.81973e-4, 1.14325e-2, 1.92354e-1, 7.46564]},
        'WISE_W3': {'fitunc': 0.43, 'range':[16., 39.], 'coeff': [2.37656e-5,-1.28563e-3,2.01740e-2,6.64242e-2,7.81181e0]},
        'WISE_W4': {'fitunc': 0.76, 'range':[16., 39.], 'coeff': [-2.16042e-3,1.14630e-1,7.78974e0]},
    }},
    'dupuy2013': {'altname': ['dupuy13','dup13'], 'bibcode': '2013Sci...341.1492D', 'sptoffset': 0, 'method': 'interpolate', 'filters': {
        'MKO_Y': {'spt': [38,38.5,39,40.], 'values': [17.4,18.81,19.26,20.24], 'rms': [0.25,0.51,0.88,0.17]}, 
        'MKO_J': {'spt': [38,38.5,39,39.5,40.], 'values': [16.43,17.87,18.39,17.68,20.09], 'rms': [0.46,0.44,0.95,0.37,0.25]}, 
        'MKO_H': {'spt': [38,38.5,39,39.5,40.], 'values': [16.82,18.2,18.77,18.08,20.6], 'rms': [0.43,0.45,1.08,0.39,0.25]}, 
        'MKO_K': {'spt': [38,38.5,39,40.], 'values': [16.93,18.27,18.89,20.7], 'rms': [0.8,0.4,0.57,0.18]}, 
        'IRAC_CH1': {'spt': [38,38.5,39,39.5,40.], 'values': [15.11,15.83,16.17,15.58,16.99], 'rms': [0.15,0.22,0.23,0.41,0.21]}, 
        'IRAC_CH2': {'spt': [38,38.5,39,39.5,40.], 'values': [13.4,13.79,14.09,13.51,14.66], 'rms': [0.21,0.12,0.2,0.43,0.28]},
    }}, 
    'faherty2012': {'altname': ['faherty12','fah12'],'bibcode': '2012ApJ...752...56F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.30, 'range' : [20., 38.], 'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, 
        'MKO_H': {'fitunc' : 0.27, 'range' : [20., 38.], 'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, 
        'MKO_K': {'fitunc' : 0.28, 'range' : [20., 38.], 'coeff' : [.0000816516, -.00469032, .0940816, -.485519, 9.76100]},
    }},
    'faherty2016': {'altname': ['faherty','faherty2016','faherty-field','fah16'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.402, 'range' : [16., 39.], 'coeff': [3.478e-5,-2.684e-3,7.771e-2,-1.058,7.157,-8.350]},
        '2MASS_H': {'fitunc' : 0.389, 'range' : [16., 39.], 'coeff' : [2.841e-5,-2.217e-3,6.551e-2,-9.174e-1,6.406,-7.496]},
        '2MASS_KS': {'fitunc' : 0.537, 'range' : [16., 39.], 'coeff' : [2.540e-5,-1.997e-3,5.978e-2,-8.481e-1,5.970,-6.704]},
        'WISE_W1': {'fitunc' : 0.365, 'range' : [16., 39.], 'coeff' : [8.337e-6,-6.897e-4,2.258e-2,-3.603e-1,2.991,-1.664e-1]},
        'WISE_W2': {'fitunc' : 0.398, 'range' : [16., 39.], 'coeff' : [8.190e-6,-6.938e-4,2.283e-2,-3.655e-1,3.032,-5.043e-1]},
        'WISE_W3': {'fitunc' : 0.446, 'range' : [16., 39.], 'coeff' : [-1.024e-6,9.477e-5,-2.573e-3,1.520e-2,3.365e-1,6.462]},
    }},
    'faherty2016-young': {'altname': ['faherty-young','fah16y'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.647, 'range' : [17., 27.], 'coeff': [4.032e-3,-1.416e-1,2.097,8.478e-1]},
        '2MASS_H': {'fitunc' : 0.634, 'range' : [17., 27.], 'coeff' : [2.642e-3,-1.049e-1,1.753,1.207]},
        '2MASS_KS': {'fitunc' : 0.640, 'range' : [17., 27.], 'coeff' : [-1.585e-2,7.338e-1,4.537]},
        'WISE_W1': {'fitunc' : 0.648, 'range' : [17., 27.], 'coeff' : [-1.397e-2,5.955e-1,5.247]},
        'WISE_W2': {'fitunc' : 0.694, 'range' : [17., 27.], 'coeff' : [-1.507e-2,5.944e-1,5.061]},
        'WISE_W3': {'fitunc' : 0.717, 'range' : [17., 27.], 'coeff' : [-1.003e-4,-1.670e-3,2.023e-1,7.529]},
    }},
    'faherty2016-group': {'altname': ['faherty-group','fah16g'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.660, 'range' : [17., 27.], 'coeff': [-3.825e-3,1.370e-1,-9.279e-1,10.141]},
        '2MASS_H': {'fitunc' : 0.603, 'range' : [17., 27.], 'coeff' : [-3.909e-3,1.346e-1,-9.347e-1,9.728]},
        '2MASS_KS': {'fitunc' : 0.556, 'range' : [17., 27.], 'coeff' : [-4.006e-3,1.378e-1,-1.031,9.916]},
        'WISE_W1': {'fitunc' : 0.551, 'range' : [17., 27.], 'coeff' : [-4.483e-3,1.505e-1,-1.208,10.403]},
        'WISE_W2': {'fitunc' : 0.616, 'range' : [17., 27.], 'coeff' : [-6.821e-3,2.322e-1,-2.133,13.322]},
        'WISE_W3': {'fitunc' : 0.427, 'range' : [17., 27.], 'coeff' : [-5.684e-3,1.993e-1,-1.987,13.972]},
    }},
    'filippazzo2015': {'altname': ['filippazzo','filippazzo15','fillippazzo','filipazo','filippazo','fil15'],'bibcode': '2015ApJ...810..158F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [3.478e-5, -2.684e-3, 7.771e-2, -1.058, 7.157, -8.350]}, 
        'WISE_W2': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [8.190e-6, -6.938e-4, 2.283e-2, -3.655e-1, 3.032, -5.043e-1]},
    }},
    'hawley2002': {'altname': ['hawley02','hawley','haw02'], 'bibcode': '2002AJ....123.3409H', 'sptoffset': 0, 'method': 'interpolate', 'filters': {
        '2MASS_J': {'spt': numpy.arange(10.,36.1,1.), 'values': [6.45,6.72,6.98,7.24,8.34,9.44,10.18,10.92,11.14,11.43,11.72,12.,12.29,12.58,12.87,13.16,14.31,14.45,14.58,14.72,14.86,14.99,15.13,15.27,15.4,15.54,15.68], 'rms': numpy.zeros(27)+0.4},
    }}, 
    'gonzales2018': {'altname': ['gonzales','gonzales18','gonzalez2018','gonzalez18','gonzalez','gon18','subdwarf'],'bibcode': '2018ApJ...864..100G', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.26, 'range' : [17., 27.], 'coeff': [0.263,8.49]},
        '2MASS_H': {'fitunc' : 0.24, 'range' : [17., 27.], 'coeff': [0.304,7.77]},
        '2MASS_KS': {'fitunc' : 0.21, 'range' : [17., 27.], 'coeff': [0.344,7.29]},
        'WISE_W1': {'fitunc' : 0.27, 'range' : [17., 27.], 'coeff': [0.241,7.72]},
        'WISE_W2': {'fitunc' : 0.24, 'range' : [17., 27.], 'coeff': [0.228,7.62]},
    }},
    'knapp2004': {'altname': ['knapp','knapp04','kna04'], 'bibcode': '2004AJ....127.3553K', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.40, 'range' : [21., 39.], 'coeff' : [-7.923e-5,3.986e-3,-6.848e-2,4.500e-1,-6.278e-1,12.03]}, 
        'MKO_K': {'fitunc': 0.30, 'range' : [21., 39.], 'coeff': [-7.351e-5,3.524e-3,-5.819e-2,3.876e-1,-6.485e-1,10.93]}, 
    }},
    'kiman2019': {'altname': ['kiman','kiman19','kim19'], 'bibcode': '2019AJ....157..231K', 'sptoffset': 0, 'method': 'interpolate', 'filters': {
        'GAIA_G': {'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26], 'values': [8.13,8.78,9.35,9.97,10.77,11.86,12.92,13.54,14.60,15.26,16.11,16.82,17.11,17.89,18.51,18.92], 'rms': [0.76,0.83,0.69,0.74,0.77,0.85,0.53,0.56,0.55,0.53,0.44,0.31,0.40,0.59,0.27,0.17]}, 
        'GAIA_Rp': {'spt': [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26], 'values': [7.20,7.76,8.26,8.81,9.55,10.54,11.50,12.07,13.02,13.62,14.45,15.14,15.42,16.21,16.80,16.98], 'rms': [0.74,0.82,0.67,0.72,0.75,0.81,0.50,0.52,0.51,0.50,0.42,0.31,0.39,0.58,0.25,0.20]}, 
        'GAIA_Bp': {'spt': [10,11,12,13,14,15,16,17,18,19], 'values': [8.99,9.78,10.51,11.27,12.25,13.38,14.67,15.48,17.04,17.92], 'rms': [0.78,0.84,0.70,0.75,0.80,0.87,0.61,0.70,0.97,0.60]}, 
    }},
    'kirkpatrick2019': {'altname': ['kirkpatrick19','kirk19','kir19'], 'bibcode': '2019ApJS..240...19K', 'sptoffset': 30, 'method': 'polynomial', 'filters': {
        '2MASS_H': {'fitunc' : 0.67, 'range' : [36,44.], 'coeff' : [-0.0344809,1.05122,-8.66856,36.9714]}, 
        'WISE_W1': {'fitunc' : 0.50, 'range' : [36,44.], 'coeff' : [-0.00124337,0.0844831,-0.276439,13.8175]}, 
        'WISE_W2': {'fitunc' : 0.29, 'range' : [36,44.], 'coeff' : [-0.00682090,0.211474,-1.59820,16.3585]}, 
        'WISE_W3': {'fitunc' : 0.41, 'range' : [36,44.], 'coeff' : [0.0159426,0.0319455,10.7315]}, 
        'IRAC_CH1': {'fitunc' : 0.39, 'range' : [36,44.], 'coeff' : [-0.00878279,0.267153,-1.78453,17.0849]}, 
        'IRAC_CH2': {'fitunc' : 0.28, 'range' : [36,44.], 'coeff' : [-0.00635074,0.203183,-1.56047,16.3304]}, 
    }},
    'kirkpatrick2021': {'altname': ['kirkpatrick','kirkpatrick21','kirk21','kir21'], 'bibcode': '2020arXiv201111616K', 'sptoffset': 20, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.60, 'range' : [30.,42.], 'coeff' : [2.1526e-7,-1.8672e-5,6.3147e-4,-9.98929e-3,7.1759e-2,-1.9013e-1,3.3790e-1,1.1808e1]}, 
        '2MASS_H': {'fitunc' : 0.57, 'range' : [30.,42.], 'coeff' : [5.4808e-7,-4.0259e-5,1.1719e-3,-1.6688e-2,1.1696e-1,-3.5467e-1,6.0330e-1,1.0966e1]}, 
        'IRAC_CH1': {'fitunc' : 0.38, 'range' : [36,42.], 'coeff' : [-9.9981e-8,4.7761e-6,-5.0487e-5,-7.1402e-4,1.7669e-2,-1.0725e-1,3.4919e-1,9.9424e0]}, 
        'IRAC_CH2': {'fitunc' : 0.21, 'range' : [36,42.], 'coeff' : [2.9971e-8,-3.3136e-6,1.3230e-4,-2.3844e-3,1.9711e-2,-6.2186e-2,1.8897e-1,1.0071e1]}, 
    }},
    'liu2016-ir': {'altname': ['liu','liu-ir','liu2016','liu2016-ir','liu16-ir','liu16-field','liu16-field-ir','liu-field','liu-field-ir'],'bibcode': '2016ApJ...833...96L', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.35, 'range' : [16., 28.], 'coeff': [0.402,8.437]},
        'MKO_J': {'fitunc' : 0.37, 'range' : [16., 28.], 'coeff': [0.368,8.062]},
        'MKO_H': {'fitunc' : 0.33, 'range' : [16., 28.], 'coeff': [0.348,7.516]},
        'MKO_K': {'fitunc' : 0.31, 'range' : [16., 28.], 'coeff': [0.319,7.258]},
        '2MASS_J': {'fitunc' : 0.38, 'range' : [16., 28.], 'coeff': [0.391,7.848]},
        '2MASS_H': {'fitunc' : 0.36, 'range' : [16., 28.], 'coeff': [0.340,7.555]},
        '2MASS_KS': {'fitunc' : 0.33, 'range' : [16., 28.], 'coeff': [0.317,7.292]},
        'WISE_W1': {'fitunc' : 0.26, 'range' : [16., 28.], 'coeff': [0.243,7.698]},
        'WISE_W2': {'fitunc' : 0.24, 'range' : [16., 28.], 'coeff': [0.236,7.466]},
    }},
    'liu2016-optical': {'altname': ['liu-optical','liu16-optical','liu16-field-optical','liu-field-optical'],'bibcode': '2016ApJ...833...96L', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.29, 'range' : [16., 28.], 'coeff': [0.404,8.677]},
        'MKO_J': {'fitunc' : 0.26, 'range' : [16., 28.], 'coeff': [0.364,8.131]},
        'MKO_H': {'fitunc' : 0.24, 'range' : [16., 28.], 'coeff': [0.331,7.820]},
        'MKO_K': {'fitunc' : 0.25, 'range' : [16., 28.], 'coeff': [0.312,7.433]},
        '2MASS_J': {'fitunc' : 0.30, 'range' : [16., 28.], 'coeff': [0.380,8.021]},
        '2MASS_H': {'fitunc' : 0.30, 'range' : [16., 28.], 'coeff': [0.332,7.723]},
        '2MASS_KS': {'fitunc' : 0.29, 'range' : [16., 28.], 'coeff': [0.308,7.491]},
        'WISE_W1': {'fitunc' : 0.29, 'range' : [16., 28.], 'coeff': [0.255,7.610]},
        'WISE_W2': {'fitunc' : 0.29, 'range' : [16., 28.], 'coeff': [0.231,7.559]},
    }},
    'liu2016-vlg': {'altname': ['liu-vlg','liu16-vlg'],'bibcode': '2016ApJ...833...96L', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.77, 'range' : [16., 27.], 'coeff': [0.799,3.689]},
        'MKO_J': {'fitunc' : 0.72, 'range' : [16., 27.], 'coeff': [0.731,3.475]},
        'MKO_H': {'fitunc' : 0.64, 'range' : [16., 27.], 'coeff': [0.633,3.778]},
        'MKO_K': {'fitunc' : 0.56, 'range' : [16., 27.], 'coeff': [0.553,3.898]},
        '2MASS_J': {'fitunc' : 0.72, 'range' : [16., 27.], 'coeff': [0.745,3.406]},
        '2MASS_H': {'fitunc' : 0.64, 'range' : [16., 27.], 'coeff': [0.632,3.734]},
        '2MASS_KS': {'fitunc' : 0.55, 'range' : [16., 27.], 'coeff': [0.573,3.699]},
        'WISE_W1': {'fitunc' : 0.46, 'range' : [16., 27.], 'coeff': [0.457,4.392]},
        'WISE_W2': {'fitunc' : 0.48, 'range' : [16., 27.], 'coeff': [0.451,4.039]},
    }},
    'liu2016-intg': {'altname': ['liu-intg','liu16-intg'],'bibcode': '2016ApJ...833...96L', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.36, 'range' : [20., 27.], 'coeff': [0.562,6.778]},
        'MKO_J': {'fitunc' : 0.50, 'range' : [20., 27.], 'coeff': [0.537,5.941]},
        'MKO_H': {'fitunc' : 0.43, 'range' : [20., 27.], 'coeff': [0.415,6.557]},
        'MKO_K': {'fitunc' : 0.39, 'range' : [20., 27.], 'coeff': [0.355,6.555]},
        '2MASS_J': {'fitunc' : 0.50, 'range' : [20., 27.], 'coeff': [0.545,5.914]},
        '2MASS_H': {'fitunc' : 0.47, 'range' : [20., 27.], 'coeff': [0.411,6.552]},
        '2MASS_KS': {'fitunc' : 0.36, 'range' : [20., 27.], 'coeff': [0.356,6.573]},
        'WISE_W1': {'fitunc' : 0.13, 'range' : [20., 27.], 'coeff': [0.237,7.436]},
        'WISE_W2': {'fitunc' : 0.18, 'range' : [20., 27.], 'coeff': [0.193,7.561]},
    }},
    'looper2008': {'altname': ['looper','looper08'],'bibcode': '2008ApJ...685.1183L', 'sptoffset': 20, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.29, 'range' : [20., 38.], 'coeff': [-5.462e-6,2.595e-4,-3.915e-3,1.663e-2,3.690e-2,1.255e-1,11.817]},
        '2MASS_H': {'fitunc' : 0.29, 'range' : [20., 38.], 'coeff' : [-4.218e-6,1.987e-4,-2.970e-3,1.261e-2,3.032e-2,1.125e-1,11.010]},
        '2MASS_KS': {'fitunc' : 0.33, 'range' : [20., 38.], 'coeff' : [-4.104e-6,1.911e-4,-2.864e-3,1.299e-2,2.565e-2,7.369e-2,10.521]},
    }},
    'pecaut2013': {'altname': ['mamajek','pecaut','pecaut13'], 'reference': 'Pecaut & Mamajek (2013)','bibcode':'2013ApJS..208....9P', 'url': 'http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt','sptoffset': 0, 'method': 'interpolate', 'filters': {
        'BESSEL_V': {\
            'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.], \
            'values': [5.76, 5.8, 5.89, 5.97, 6.19, 6.34, 6.57, 6.79, 6.98, 7.04, 7.36, 7.6, 7.8, 8.01, 8.15, 8.47, 8.69, 8.91, 9.2, 9.69, 9.97, 10.3, 10.7, 11.14, 12.19, 12.8, 13.57, 14.3, 15.51, 16.62, 17.07, 17.81, 18.42, 18.84, 19.14, 19.36, 19.75, 20, 20.5, 20.9, 21.7, 22.3, 23.1],\
            'rms': numpy.zeros(43)+0.05},
        'GAIA_G': {\
            'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.], \
            'values': [5.53, 5.56, 5.63, 5.69, 5.89, 6.04, 6.26, 6.42, 6.55, 6.62, 6.95, 7.16, 7.32, 7.49, 7.58, 7.84, 8, 8.26, 8.46, 8.87, 9.12, 9.38, 9.68, 10.05, 10.9, 11.39, 12.02, 12.56, 13.5, 14.48, 14.32, 14.83, 15.33, 15.73, 16.05, 16.29, 16.33, 16.36, 16.83, 17.24, 17.76, 18.32, 18.86, 19.25, 19.3, 20],\
            'rms': numpy.zeros(46)+0.05},
        '2MASS_J': {\
            'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,34.5,35.,35.5,36.,37.,37.5,38.,38.5,39.,39.5,40.,40.5,41.,41.5,42.],\
            'values': [4.29, 4.31, 4.37, 4.41, 4.57, 4.63, 4.76, 4.85, 4.92, 4.94, 5.18, 5.3, 5.41, 5.56, 5.6, 5.78, 5.92, 6.04, 6.19, 6.51, 6.69, 6.89, 7.01, 7.4, 8.02, 8.39, 8.79, 9.25, 9.93, 10.28, 10.47, 10.76, 10.68, 11.23, 11.45, 11.53, 11.78, 11.84, 12.14, 12.34, 12.93, 13.17, 13.6, 13.99, 14.34, 14.47, 14.47, 14.46, 14.34, 14.32, 14.45, 14.67, 14.8, 15.02, 15.28, 15.61, 15.77, 16.51, 17.24, 18.45, 19.67, 19.52, 21, 22, 23, 25, 28.3],\
            'rms': numpy.zeros(67)+0.05},
        '2MASS_KS': {\
            'spt': [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,34.5,35.,35.5,36.,37.,37.5,38.,39.,39.5,40.5,41.,41.5,42.],\
            'values': [3.81, 3.82, 3.87, 3.9, 4.04, 4.07, 4.16, 4.21, 4.25, 4.26, 4.48, 4.57, 4.66, 4.78, 4.82, 4.98, 5.11, 5.22, 5.36, 5.67, 5.85, 6.06, 6.27, 6.54, 7.19, 7.55, 7.93, 8.36, 9.01, 9.32, 9.47, 9.76, 9.97, 10.11, 10.22, 10.3, 10.45, 10.55, 10.8, 10.9, 11.3, 11.4, 11.82, 12.27, 12.62, 12.7, 12.74, 12.9, 12.87, 13.19, 13.69, 14.51, 14.66, 14.9, 15.18, 15.54, 16.36, 16.85, 17.43, 18.48, 19.32, 21.5, 23, 23.5, 24],\
            'rms': numpy.zeros(65)+0.05},
    }},
    'tinney2003': {'altname': ['tinney','tinney03'],'bibcode': '2003AJ....126..975T', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'COUSINS_I': {'fitunc' : 0.37, 'range' : [20., 37.5], 'coeff': [-2.49821e-6,1.04398e-3,-6.49719e-2,1.56038,-1.58296e1,7.22089e1]},
        'UKIDSS_Z': {'fitunc' : 0.29, 'range' : [20., 37.5], 'coeff': [-9.97226e-7,1.05950e-4,-4.57019e-3,1.02898e-1,-1.29357e0,8.96822e0,-3.08010e1,4.99447e1]},
        'UKIDSS_K': {'fitunc' : 0.40, 'range' : [20., 37.5], 'coeff': [1.14139e-5,-8.86885e-4,2.68071e-2,-3.89554e-1,2.95440e0,8.14626e-1]},
        '2MASS_KS': {'fitunc' : 0.38, 'range' : [20., 37.5], 'coeff': [-1.25074e-5,1.63124e-3,-7.42418e-2,1.54509,-1.47407e1,6.27861e1]},
        'UKIRT_J': {'fitunc' : 0.30, 'range' : [20., 37.5], 'coeff': [-9.91110e-7,1.05811e-4,-4.58399e-3,1.03572e-1,-1.30526e0,9.06701e0,-3.13411e1,5.04642e1]},
        '2MASS_J': {'fitunc' : 0.36, 'range' : [20., 37.5], 'coeff': [-2.80824e-6,3.41146e-4,-1.73848e-2,4.82120e-1,-7.86911,7.57222e1,-3.98105e2,8.94012e2]},
    }},
    'tinney2014': {'altname': ['tinney14'],'bibcode': '2014ApJ...796...39T', 'sptoffset': 0, 'method': 'interpolate', 'filters': {
        'MKO_J': {'spt': [36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,42], 'values': [15.22,15.49,16.39,16.66,17.9,18.35,19.08,20.32,22.39,22.18,25.76], 'rms': [0.31,0.37,0.72,0.36,0.46,0.9,0.97,1.25,1.,0.76,3.52]}, 
        'WISE_W2': {'spt': [36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,42], 'values': [12.86,13.28,13.39,13.44,13.75,13.92,14.28,14.65,15.2,14.78,15.76], 'rms': [0.17,0.48,0.27,0.23,0.22,0.24,0.46,0.35,1.,0.77,2.15]},
    }},
    'zhang2013': {'altname': ['zhang','zhang13','zha13','zhang_dwarf'],'bibcode': '2013MNRAS.434.2664Z', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'SDSS_R': {'fitunc' : 0.71, 'range' : [11., 28.], 'coeff': [-1.6751e-4,8.4503e-3,-1.5169e-1,1.1111,-1.8452,9.8326]},
        'SDSS_I': {'fitunc' : 0.70, 'range' : [11., 28.], 'coeff': [-1.0634e-4,5.3482e-3,-9.6139e-2,7.1114e-1,-1.0400,8.6772]},
        'SDSS_Z': {'fitunc' : 0.60, 'range' : [11., 28.], 'coeff': [-6.2600e-5,2.9652e-3,-4.9888e-2,3.3379e-1,-4.7257e-2,7.8303]},
        '2MASS_J': {'fitunc' : 0.52, 'range' : [11., 29.], 'coeff': [-5.3057e-5,2.7283e-3,-4.9840e-2,3.6508e-1,-2.9785e-1,6.2826]},
        '2MASS_H': {'fitunc' : 0.54, 'range' : [11., 29.], 'coeff': [-5.8608e-5,3.0343e-3,-5.5512e-2,4.0397e-1,-3.9566e-1,5.7798]},
        '2MASS_KS': {'fitunc' : 0.53, 'range' : [11., 29.], 'coeff': [-5.5258e-5,2.8626e-3,-5.2136e-2,3.7299e-1,-3.0303e-1,5.4490]},
    }},
    'zhang2017': {'altname': ['zhang17','zha17'],'bibcode': '2017MNRAS.464.3040Z', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.40, 'range' : [10., 27.], 'coeff': [8.53625e-4,-1.76459e-2,3.17384e-1,8.64788]},
        'MKO_H': {'fitunc' : 0.40, 'range' : [10., 27.], 'coeff': [2.90020e-4,-4.54248e-3,2.71013e-1,8.19731]},
        '2MASS_J': {'fitunc' : 0.40, 'range' : [10., 27.], 'coeff': [8.48172e-4,-1.75984e-2,3.16187e-1,8.68342]},
        '2MASS_H': {'fitunc' : 0.41, 'range' : [10., 27.], 'coeff': [4.32261e-4,-7.53663e-3,2.81607e-1,8.18494]},
    }},
    'zhang2019lsd': {'altname': ['zhang19lsd','zha19lsd','zhang_lsubdwarf','zhang_lsd'],'bibcode': '2019MNRAS.486.1260Z', 'sptoffset': 20., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.1063, 'range' : [20., 27.], 'coeff': [0.23129,12.1660]},
        'WISE_W1': {'fitunc' : 0.1735, 'range' : [20., 27.], 'coeff': [0.1787,10.3005]},
        'WISE_W2': {'fitunc' : 0.1647, 'range' : [20., 27.], 'coeff': [0.1830,9.9392]},
    }},
    'zhang2019tsd': {'altname': ['zhang19tsd','zha19tsd','zhang_tsubdwarf','zhang_tsd'],'bibcode': '2019MNRAS.486.1260Z', 'sptoffset': 30., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.6426, 'range' : [35.5, 38.], 'coeff': [0.8742,11.2935]},
        'MKO_J': {'fitunc' : 0.8135, 'range' : [35.5, 39.], 'coeff': [1.2294,7.9159]},
        'MKO_H': {'fitunc' : 0.7675, 'range' : [35.5, 39.], 'coeff': [1.0964,9.2711]},
        'MKO_K': {'fitunc' : 0.7876, 'range' : [35.5, 38.], 'coeff': [1.2030,9.1042]},
        'MKO_KS': {'fitunc' : 0.7865, 'range' : [35.5, 38.], 'coeff': [1.1996,8.9844]},
        'WISE_W1': {'fitunc' : 0.4584, 'range' : [35.5, 39.], 'coeff': [0.7678,10.5810]},
        'WISE_W2': {'fitunc' : 0.2887, 'range' : [35.5, 39.], 'coeff': [0.4232,10.2525]},
    }},
    'zhang2019esd': {'altname': ['zhang19esd','zha19esd','zhang_esubdwarf','zhang_esd','zhang_usd'],'bibcode': '2019MNRAS.486.1260Z', 'sptoffset': 20., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.122, 'range' : [20., 23.], 'coeff': [0.2443,11.5748]},
        'WISE_W1': {'fitunc' : 0.1084, 'range' : [20., 27.], 'coeff': [0.2610,10.0081]},
        'WISE_W2': {'fitunc' : 0.1213, 'range' : [20., 27.], 'coeff': [0.2250,9.7670]},
    }},
}

# Empirical relations - SpT to chi value (Halpha EW to LHa/Lbol)
SPT_CHI_RELATIONS = {
    'schmidt2014': {'altname': ['schmidt','schmidt14'], 'reference': 'Schmidt et al. (2014)','bibcode':'2014PASP..126..642S', 'sptoffset': 0, 'method': 'interpolate', 'scale': 1.e-6,
        'spt': [17,18,19,20,21,22,23,24,25,26,27], \
        'values': [10.28,4.26,2.52,1.98,2.25,2.11,1.67,1.16,1.46,1.23,0.73],\
        'scatter': [3.13,1.18,0.58,0.27,0.11,0.36,0.22,0.3,0.28,0.3,0.3],\
        },
    'douglas2014': {'altname': ['douglas','douglas14'], 'reference': 'Douglas et al. (2014)','bibcode':'2014ApJ...795..161D', 'sptoffset': 0, 'method': 'interpolate', 'scale': 1.e-5,
        'spt': [10,11,12,13,14,15,16,17,18,19], \
        'values': [6.6453,6.0334,5.2658,4.4872,3.5926,2.4768,1.7363,1.2057,0.6122,0.3522],\
        'scatter': [0.6207,0.5326,0.5963,0.4967,0.5297,0.4860,0.3475,0.3267,0.2053,0.1432],\
        },
}

# telscope locations
TELESCOPES = {
    'KECK': {'lat': 19.8263*u.deg, 'lon': -155.4783*u.deg, 'height': 4160*u.m, 'altname': ['MAUNA_KEA','SUBARU','IRTF','CHFT','UKIRT','GEMINI_N','GEMINI_NORTH']},
    'VLT': {'lat': -24.6275*u.deg, 'lon': -70.4044*u.deg, 'height': 2636*u.m, 'altname': ['PARANAL','CERRO_PARANAL','VERY_LARGE_TELESCOPE','ESA']},
}

# Bibliography defaults
BIBFILE = os.path.join(CITATION_RESOURCES_FOLDER,'splat_bibs.bib')
JOURNALS_LONGNAMES = {\
    'aap': 'Astronomy & Astrophysics',\
    'actaa': 'Acta Astronomica',\
    'aj': 'Astronomical Journal',\
    'apj': 'Astrophysical Journal',\
    'apjl': 'Astrophysical Journal Letters',\
    'apjs': 'Astrophysical Journal Supplement',\
    'araa': 'Annual Reviews of Astronomy & Astrophysics',\
    'icarus': 'Icarus',\
    'mnras': 'Monthly Notices of the Royal Astronomical Society',\
    'nat': 'Nature',\
    'natas': 'Nature Astronomy',\
    'pasp': 'Proceedings of the Astronomical Society of the Pacific',\
    'pnas': 'Proceedings of the National Academy of Sciences',\
    'rnaas': 'Research Notes of the American Astronomical Society',\
    'sci': 'Science',\
    'solphys': 'Solar Physics',\
}
CITATION_URL_BASE = 'https://ui.adsabs.harvard.edu/abs/'

#######################################################
#######################################################
###############   DISPLAY ON LOAD IN  #################
#######################################################
#######################################################

print('\n\nWelcome to the Spex Prism Library Analysis Toolkit (SPLAT)!')
print('You are currently using version {}\n'.format(VERSION))
print('If you make use of any features of this toolkit for your research, please remember to cite the SPLAT paper:')
print('\n{}; Bibcode: {}\n'.format(CITATION,BIBCODE))
print('If you make use of any spectra or models in this toolkit, please remember to cite the original source.')
print('Please report any errors are feature requests to our github page, {}\n\n'.format(GITHUB_URL))
if ERROR_CHECKING==True: print('Currently running in error checking mode')

