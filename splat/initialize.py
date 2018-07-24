# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
.. note::
         These commands initialize the SPLAT code 
"""

import os
import sys
import numpy
import pandas
from astropy import units as u


# program constants
VERSION = '2018.07.11'
__version__ = VERSION
SPLAT_URL = 'http://splat.physics.ucsd.edu/splat/'
DOCUMENTATION_URL = 'http://pono.ucsd.edu/~adam/splat/'
GITHUB_URL = 'https://github.com/aburgasser/splat/'
BIBCODE = '2017arXiv170700062B'
EMAIL = 'aburgasser@gmail.com'
INSTRUMENT_DEFINITION_FILE = 'instrument.txt'
# no longer used
#DB_PHOTOMETRY_FILE = 'photometry_data.txt'
BIBFILE = 'splat_bibs.bib'
TMPFILENAME = 'splattmpfile'
HOME_FOLDER = os.path.expanduser('~')
FILTER_FOLDER = '/resources/Filters/'
SPECTRAL_MODEL_FOLDER = '/resources/SpectralModels/'
EVOLUTIONARY_MODEL_FOLDER = '/resources/EvolutionaryModels/'
DOCS_FOLDER = '/docs/'
DOCS_INDEX_HTML = '/docs/_build/html/index.html'
WEB_HTML_BASE = '/docs/_templates/'
DB_FOLDER = '/db/'
ACCESS_FILE = '.splat_access'
EXTERNAL_SPECTRAL_MODELS_FILE = '.splat_spectral_models'
EXTERNAL_EVOLUTIONARY_MODELS_FILE = '.splat_evolutionary_models'
MONTHS = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

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
    'Russell van Linge',
]

# SET VARIOUS ENVIRONMENTAL VARIABLES
#set the SPLAT PATH, either from set environment variable or from sys.path
SPLAT_PATH = './'
if os.environ.get('SPLAT_PATH') != None:
    SPLAT_PATH = os.environ['SPLAT_PATH']
# get from PYTHONPATH
if os.environ.get('PYTHONPATH') != None and SPLAT_PATH == './':
    path = os.environ['PYTHONPATH']
    for i in path.split(':'):
        if 'splat' in i:
            SPLAT_PATH = i
# get from system path
if SPLAT_PATH == './':
    checkpath = ['splat' in r for r in sys.path]
    if max(checkpath):
        SPLAT_PATH = sys.path[checkpath.index(max(checkpath))]

#set user SPLAT model path from environmental variable
SPLAT_USER_MODELS = './'
if os.environ.get('SPLAT_MODELS') != None:
    SPLAT_USER_MODELS = os.environ['SPLAT_MODELS']


# Unit standards
DEFAULT_WAVE_UNIT = u.micron
DEFAULT_FLUX_UNIT = u.erg/u.s/u.cm/u.cm/u.micron
DEFAULT_SED_UNIT = u.erg/u.s/u.cm/u.cm


# Spectral data parameters
DATA_FOLDER = '/resources/Spectra/SPEX-PRISM/'
EXTERNAL_DATA_FILE = '.splat_data'
DB_SOURCES_FILE = 'source_data.txt'
DB_SPECTRA_FILE = 'spectral_data.txt'
DB_SPECTRA_DEFAULT_PARAMETERS = {
    'DATA_FILE': {'altname': ['FILE','FILENAME'], 'type': str},
    'DATA_KEY': {'altname': [], 'type': int},
    'INSTRUMENT': {'altname': [], 'type': str},
    'BIBCODE': {'altname': ['DATA_BIBCODE','REFERENCE','REF','BIBCODE','BIB','DATA_REFERENCE','DATA_REF'], 'type': str},
    'OBSERVATION_DATE': {'altname': ['OBSDATE','OBS_DATE','DATE'], 'type': str},
    'OBSERVATION_MJD': {'altname': ['MJD','OBS_MJD','JULIAN_DATE'], 'type': float},
    'PROGRAM': {'altname': [], 'type': str},
    'OBSERVER': {'altname': [], 'type': str},
    'AIRMASS': {'altname': ['Z'], 'type': float},
    'INTEGRATION': {'altname': ['TINT','TIME','INT_TIME','INTEGRATION_TIME'], 'type': float},
    'CONDITIONS': {'altname': ['WEATHER'], 'type': float},
    'SOURCE_KEY': {'altname': [], 'type': int},
    'NAME': {'altname': ['SOURCE_NAME'], 'type': str},
    'DESIGNATION': {'altname': ['DESIG'], 'type': str},
    'RA': {'altname': ['RIGHT_ASCENSION'], 'type': float},
    'DEC': {'altname': ['DECLINATION'], 'type': float},
}

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
    'T9.0': 11536} #'11536_10509.fits'}

# subdwarf spectral standards
STDS_DSD_SPEX_KEYS = { \
    'd/sdM4.0': 10523, #'11670_11134.fits',\  10501
    'd/sdM5.0': 10232, #'10265_10045.fits',\  11498
    'd/sdM6.0': 10198, #'10265_10045.fits',\  11208
    'd/sdM7.0': 10863, #'10197_11074.fits',\  11704
    'd/sdM8.0': 10040, #'10197_11074.fits',\  11190
    'd/sdM9.0': 10367, #'10197_11074.fits',\  11086
    'd/sdL0.0': 10146, #'11972_10248.fits',\
    'd/sdL1.0': 10506, #'11972_10248.fits',\
    'd/sdL7.0': 10552, #'11972_10248.fits',\
    }

STDS_SD_SPEX_KEYS = { \
    'sdM2.0': 10223, #'11670_11134.fits',\
    'sdM4.0': 10528, #'11670_11134.fits',\
    'sdM5.0': 10221, #'11670_11134.fits',\
    'sdM5.5': 11670, #'11670_11134.fits',\
    'sdM6.0': 10265, #'10265_10045.fits',\
    'sdM7.0': 10197, #'10197_11074.fits',\
    'sdM8.0': 10123, #'10123_10145.fits',\
    'sdM9.5': 10188, #'10188_10700.fits',\
    'sdL0.0': 11973, #'11972_10248.fits',\
    'sdL3.5': 10364, #'10364_10946.fits',\
    'sdL4.0': 10203} #'10203_11241.fits'}

# extreme subdwarf spectral standards
STDS_ESD_SPEX_KEYS = { \
    'esdM0.0': 10763, #'10229_10163.fits',\
    'esdM4.0': 10366, #'10229_10163.fits',\
    'esdM5.0': 10229, #'10229_10163.fits',\
#    'esdM6.5': '_10579.fits',\
    'esdM6.5': 10359, #'10521_10458.fits',\
    'esdM7.5': 10521, #'10521_10458.fits',\
    'esdM8.5': 10278} #'10278_10400.fits'}
# EMPTY DICTIONARY

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
    'L6.0gamma': 10455}  # 2244+2043

STDS_INTG_SPEX_KEYS = { \
    'M8.0beta': 12572,  # 0019+4614
    'L0.0beta': 11113,  # 1552+2948
    'L1.0beta': 10845,  # 0227-1624 - proposed by Cruz et al. 2017
    'L2.0beta': 11304,  # 0602+3910
    'L3.0beta': 11070,  # 1726+1538
    'L6.0beta': 10678}  # 0103+1935


# containters for spectra
#STDS_DWARF_SPEX = {}
#STDS_SD_SPEX = {}
#STDS_ESD_SPEX = {}

# filters
# this information is from the SVO filter profile service: http://svo2.cab.inta-csic.es/svo/theory/fps/
FILTERS = { \
    '2MASS_J': {'file': 'j_2mass.txt', 'description': '2MASS J-band', 'zeropoint': 1594.0, 'method': 'vega', 'rsr': True, 'altname': []}, \
    '2MASS_H': {'file': 'h_2mass.txt', 'description': '2MASS H-band', 'zeropoint': 1024.0, 'method': 'vega', 'rsr': True, 'altname': []}, \
    '2MASS_KS': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega', 'rsr': True, 'altname': ['2MASS_K']}, \
#    '2MASS_K': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega'}, \
#    '2MASS_Ks': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega'}, \
    'BESSEL_I': {'file': 'bessel_i.txt', 'description': 'Bessel I-band', 'zeropoint': 2405.3, 'method': 'vega', 'rsr': False, 'altname': ['I']}, \
    'COUSINS_I': {'file': 'i_cousins.txt', 'description': 'Cousins I-band', 'zeropoint': 2405.3, 'method': 'vega', 'rsr': False, 'altname': ['IC']}, \
    'FOURSTAR_J': {'file': 'fourstar-j.txt', 'description': 'FOURSTAR J-band', 'zeropoint': 1581.2, 'method': 'vega', 'rsr': False, 'altname': ['4star j']}, \
    'FOURSTAR_J1': {'file': 'fourstar-j1.txt', 'description': 'FOURSTAR J1-band', 'zeropoint': 1978.7, 'method': 'vega', 'rsr': False, 'altname': ['4star j1']}, \
    'FOURSTAR_J2': {'file': 'fourstar-j2.txt', 'description': 'FOURSTAR J2-band', 'zeropoint': 1774.5, 'method': 'vega', 'rsr': False, 'altname': ['4star j2']}, \
    'FOURSTAR_J3': {'file': 'fourstar-j3.txt', 'description': 'FOURSTAR J3-band', 'zeropoint': 1488.8, 'method': 'vega', 'rsr': False, 'altname': ['4star j3']}, \
    'FOURSTAR_H': {'file': 'fourstar-h.txt', 'description': 'FOURSTAR H-band', 'zeropoint': 1054.9, 'method': 'vega', 'rsr': False, 'altname': ['4star h']}, \
    'FOURSTAR_H_SHORT': {'file': 'fourstar-hshort.txt', 'description': 'FOURSTAR H short', 'zeropoint': 1119.1, 'method': 'vega', 'rsr': False, 'altname': ['4star h short','4star h-short','4star hs','fourstar hs','fourstar h1']}, \
    'FOURSTAR_H_LONG': {'file': 'fourstar-hlong.txt', 'description': 'FOURSTAR H long', 'zeropoint': 980.7, 'method': 'vega', 'rsr': False, 'altname': ['4star h long','4star h-long','4star hl','fourstar hl','fourstar h2']}, \
    'FOURSTAR_KS': {'file': 'fourstar-ks.txt', 'description': 'FOURSTAR Ks-band', 'zeropoint': 675.7, 'method': 'vega', 'rsr': False, 'altname': ['4star k','4star ks','fourstar k']}, \
    'FOURSTAR_1.18': {'file': 'fourstar-118.txt', 'description': 'FOURSTAR 1.18 micron narrow band', 'zeropoint': 675.7, 'method': 'vega', 'rsr': False, 'altname': ['4star 1.18','4star 118','fourstar 118']}, \
    'FOURSTAR_2.09': {'file': 'fourstar-209.txt', 'description': 'FOURSTAR 2.09 micron narrow band', 'zeropoint': 675.7, 'method': 'vega', 'rsr': False, 'altname': ['4star 2.09','4star 209','fourstar 209']}, \
    'HAWK_Y': {'file': 'hawk-y.txt', 'description': 'HAWK Y-band', 'zeropoint': 2092.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_J': {'file': 'hawk-j.txt', 'description': 'HAWK J-band', 'zeropoint': 1543.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_H': {'file': 'hawk-h.txt', 'description': 'HAWK H-band', 'zeropoint': 1053.6, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_H2': {'file': 'hawk-h2.txt', 'description': 'HAWK H2-band', 'zeropoint': 688.8, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_CH4': {'file': 'hawk-ch4.txt', 'description': 'HAWK CH4-band', 'zeropoint': 1093.4, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_KS': {'file': 'hawk-ks.txt', 'description': 'HAWK Ks-band', 'zeropoint': 675.3, 'method': 'vega', 'rsr': False, 'altname': ['hawk k']}, \
    'HAWK_BRG': {'file': 'hawk-brg.txt', 'description': 'HAWK Brackett Gamma', 'zeropoint': 638.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'HAWK_NB1060': {'file': 'hawk-nb1060.txt', 'description': 'HAWK Narrow Band 1060', 'zeropoint': 2003.27, 'method': 'vega', 'rsr': False, 'altname': ['hawk 1060']}, \
    'HAWK_NB1190': {'file': 'hawk-nb1190.txt', 'description': 'HAWK Narrow Band 1190', 'zeropoint': 1697.50, 'method': 'vega', 'rsr': False, 'altname': ['hawk 1190']}, \
    'HAWK_NB2090': {'file': 'hawk-nb2090.txt', 'description': 'HAWK Narrow Band 2090', 'zeropoint': 706.68, 'method': 'vega', 'rsr': False, 'altname': ['hawk 2090']}, \
    'IRAC_CH1': {'file': 'irac1.txt', 'description': 'IRAC Channel 1 (3.6 micron)', 'zeropoint': 280.9, 'method': 'vega', 'rsr': True, 'altname': ['irac 1','irac 3.6','[3.6]']}, \
    'IRAC_CH2': {'file': 'irac2.txt', 'description': 'IRAC Channel 2 (4.5 micron)', 'zeropoint': 179.7, 'method': 'vega', 'rsr': True, 'altname': ['irac 2','irac 4.5','[4.5]']}, \
    'IRAC_CH3': {'file': 'irac3.txt', 'description': 'IRAC Channel 3 (5.8 micron)', 'zeropoint': 115.0, 'method': 'vega', 'rsr': True, 'altname': ['irac 3','irac 5.8','[5.8]']}, \
    'IRAC_CH4': {'file': 'irac4.txt', 'description': 'IRAC Channel 4 (8.0 micron)', 'zeropoint': 64.13, 'method': 'vega', 'rsr': True, 'altname': ['irac 4','irac 8.0','[8.0]']}, \
    'KEPLER': {'file': 'Kepler.txt', 'description': 'Kepler bandpass', 'zeropoint': 3033.1, 'method': 'vega', 'rsr': False, 'altname': ['kep','kepler k','kp']}, \
    'MKO_J_ATM': {'file': 'j_atm_mko.txt', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_H_ATM': {'file': 'h_atm_mko.txt', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_K_ATM': {'file': 'k_atm_mko.txt', 'description': 'MKO K-band + atmosphere', 'zeropoint': 647.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_J': {'file': 'mko_j.txt', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_H': {'file': 'mko_h.txt', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_K': {'file': 'mko_ks.txt', 'description': 'MKO K-band', 'zeropoint': 647.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'MKO_KP': {'file': 'mko_kp.txt', 'description': 'MKO Kp-band', 'zeropoint': 693.7, 'method': 'vega', 'rsr': False, 'altname': ['mko k prime']}, \
    'MKO_LP': {'file': 'mko_lp.txt', 'description': 'MKO Lp-band', 'zeropoint': 248.3, 'method': 'vega', 'rsr': False, 'altname': ['mko l','mko l prime']}, \
    'MKO_MP': {'file': 'mko_mp.txt', 'description': 'MKO Mp-band', 'zeropoint': 164.7, 'method': 'vega', 'rsr': False, 'altname': ['mko m','mko m prime']}, \
    'NICMOS_F090M': {'file': 'nic1_f090m.txt', 'description': 'NICMOS F090M', 'zeropoint': 2255.0, 'method': 'vega', 'rsr': False, 'altname': ['F090M']}, \
    'NICMOS_F095N': {'file': 'nic1_f095n.txt', 'description': 'NICMOS F095N', 'zeropoint': 2044.6, 'method': 'vega', 'rsr': False, 'altname': ['F095N']}, \
    'NICMOS_F097N': {'file': 'nic1_f097n.txt', 'description': 'NICMOS F097N', 'zeropoint': 2275.4, 'method': 'vega', 'rsr': False, 'altname': ['F097N']}, \
    'NICMOS_F108N': {'file': 'nic1_f108n.txt', 'description': 'NICMOS F108N', 'zeropoint': 1937.3, 'method': 'vega', 'rsr': False, 'altname': ['F108N']}, \
    'NICMOS_F110M': {'file': 'nic1_f110m.txt', 'description': 'NICMOS F110M', 'zeropoint': 1871.8, 'method': 'vega', 'rsr': False, 'altname': ['F110M']}, \
    'NICMOS_F110W': {'file': 'nic1_f110w.txt', 'description': 'NICMOS F110W', 'zeropoint': 1768.5, 'method': 'vega', 'rsr': False, 'altname': ['F110W']}, \
    'NICMOS_F113N': {'file': 'nic1_f113n.txt', 'description': 'NICMOS F113N', 'zeropoint': 1821.0, 'method': 'vega', 'rsr': False, 'altname': ['F113N']}, \
    'NICMOS_F140W': {'file': 'nic1_f140w.txt', 'description': 'NICMOS F140W', 'zeropoint': 1277.1, 'method': 'vega', 'rsr': False, 'altname': ['F140W']}, \
    'NICMOS_F145M': {'file': 'nic1_f145m.txt', 'description': 'NICMOS F145M', 'zeropoint': 1242.0, 'method': 'vega', 'rsr': False, 'altname': ['F145M']}, \
    'NICMOS_F160W': {'file': 'nic1_f160w.txt', 'description': 'NICMOS F160W', 'zeropoint': 1071.7, 'method': 'vega', 'rsr': False, 'altname': ['F160W']}, \
    'NICMOS_F164N': {'file': 'nic1_f164n.txt', 'description': 'NICMOS F164N', 'zeropoint': 1003.0, 'method': 'vega', 'rsr': False, 'altname': ['F164N']}, \
    'NICMOS_F165M': {'file': 'nic1_f165m.txt', 'description': 'NICMOS F165M', 'zeropoint': 1023.6, 'method': 'vega', 'rsr': False, 'altname': ['F165M']}, \
    'NICMOS_F166N': {'file': 'nic1_f166n.txt', 'description': 'NICMOS F166N', 'zeropoint': 1047.7, 'method': 'vega', 'rsr': False, 'altname': ['F166N']}, \
    'NICMOS_F170M': {'file': 'nic1_f170m.txt', 'description': 'NICMOS F170M', 'zeropoint': 979.1, 'method': 'vega', 'rsr': False, 'altname': ['F170M']}, \
    'NICMOS_F187N': {'file': 'nic1_f187n.txt', 'description': 'NICMOS F187N', 'zeropoint': 803.7, 'method': 'vega', 'rsr': False, 'altname': ['F187N']}, \
    'NICMOS_F190N': {'file': 'nic1_f190n.txt', 'description': 'NICMOS F190N', 'zeropoint': 836.5, 'method': 'vega', 'rsr': False, 'altname': ['F190N']}, \
    'NIRC2_J': {'file': 'nirc2-j.txt', 'description': 'NIRC2 J-band', 'zeropoint': 1562.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_H': {'file': 'nirc2-h.txt', 'description': 'NIRC2 H-band', 'zeropoint': 1075.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_HCONT': {'file': 'nirc2-hcont.txt', 'description': 'NIRC2 H-continuum band', 'zeropoint': 1044.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_K': {'file': 'nirc2-k.txt', 'description': 'NIRC2 K-band', 'zeropoint': 648.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_KP': {'file': 'nirc2-kp.txt', 'description': 'NIRC2 Kp-band', 'zeropoint': 689.3, 'method': 'vega', 'rsr': False, 'altname': ['nirc2 k prime']}, \
    'NIRC2_KS': {'file': 'nirc2-ks.txt', 'description': 'NIRC2 Ks-band', 'zeropoint': 676.2, 'method': 'vega', 'rsr': False, 'altname': ['nirc2 k short']}, \
    'NIRC2_KCONT': {'file': 'nirc2-kcont.txt', 'description': 'NIRC2 K continuum-band', 'zeropoint': 605.9, 'method': 'vega', 'rsr': False, 'altname': ['nirc2 k continuum']}, \
    'NIRC2_FE2': {'file': 'nirc2-fe2.txt', 'description': 'NIRC2 Fe II', 'zeropoint': 1019.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRC2_LP': {'file': 'nirc2-lp.txt', 'description': 'NIRC2 LP', 'zeropoint': 248.0, 'method': 'vega', 'rsr': False, 'altname': ['nirc2 l prime','nirc2 l']}, \
    'NIRC2_M': {'file': 'nirc2-ms.txt', 'description': 'NIRC2 M', 'zeropoint': 165.8, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'NIRCAM_F070W': {'file': 'jwst-nircam-F070W.txt', 'description': 'JWST NIRCAM F070W (wide 0.70 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F090W': {'file': 'jwst-nircam-F090W.txt', 'description': 'JWST NIRCAM F090W (wide 0.90 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F115W': {'file': 'jwst-nircam-F115W.txt', 'description': 'JWST NIRCAM F115W (wide 1.15 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F140M': {'file': 'jwst-nircam-F140M.txt', 'description': 'JWST NIRCAM F140M (medium 1.40 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F150W': {'file': 'jwst-nircam-F150W.txt', 'description': 'JWST NIRCAM F150W (wide 1.50 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F150W2': {'file': 'jwst-nircam-F150W2.txt', 'description': 'JWST NIRCAM F150W2 (wide 1.50 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F162M': {'file': 'jwst-nircam-F162M.txt', 'description': 'JWST NIRCAM F162M (medium 1.62 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F164N': {'file': 'jwst-nircam-F164N.txt', 'description': 'JWST NIRCAM F164N (narrow 1.64 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F182M': {'file': 'jwst-nircam-F182M.txt', 'description': 'JWST NIRCAM F182M (medium 1.82 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F187N': {'file': 'jwst-nircam-F187N.txt', 'description': 'JWST NIRCAM F187N (narrow 1.87 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F200W': {'file': 'jwst-nircam-F200W.txt', 'description': 'JWST NIRCAM F200W (wide 2.00 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F210W': {'file': 'jwst-nircam-F210M.txt', 'description': 'JWST NIRCAM F210M (medium 2.10 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F212N': {'file': 'jwst-nircam-F212N.txt', 'description': 'JWST NIRCAM F212N (narrow 2.12 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F250M': {'file': 'jwst-nircam-F250M.txt', 'description': 'JWST NIRCAM F250M (medium 2.50 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F277W': {'file': 'jwst-nircam-F277W.txt', 'description': 'JWST NIRCAM F277W (wide 2.77 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F300M': {'file': 'jwst-nircam-F300M.txt', 'description': 'JWST NIRCAM F300M (medium 3.00 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F322W2': {'file': 'jwst-nircam-F322W2.txt', 'description': 'JWST NIRCAM F322W2 (wide 3.22 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F323N': {'file': 'jwst-nircam-F323N.txt', 'description': 'JWST NIRCAM F323N (narrow 3.23 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F335M': {'file': 'jwst-nircam-F335M.txt', 'description': 'JWST NIRCAM F335M (medium 3.35 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F356W': {'file': 'jwst-nircam-F356W.txt', 'description': 'JWST NIRCAM F356W (wide 3.56 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F360M': {'file': 'jwst-nircam-F360M.txt', 'description': 'JWST NIRCAM F360M (medium 3.60 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F405N': {'file': 'jwst-nircam-F405N.txt', 'description': 'JWST NIRCAM F405N (narrow 4.05 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F410M': {'file': 'jwst-nircam-F410M.txt', 'description': 'JWST NIRCAM F410M (medium 4.10 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F430M': {'file': 'jwst-nircam-F430M.txt', 'description': 'JWST NIRCAM F430M (medium 4.30 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F444W': {'file': 'jwst-nircam-F444W.txt', 'description': 'JWST NIRCAM F444W (wide 4.44 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F460M': {'file': 'jwst-nircam-F460M.txt', 'description': 'JWST NIRCAM F460M (medium 4.60 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F466N': {'file': 'jwst-nircam-F466N.txt', 'description': 'JWST NIRCAM F466N (narrow 4.66 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F470N': {'file': 'jwst-nircam-F470N.txt', 'description': 'JWST NIRCAM F470N (narrow 4.70 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'NIRCAM_F480M': {'file': 'jwst-nircam-F480M.txt', 'description': 'JWST NIRCAM F480M (medium 4.80 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altname': []}, \
    'PANSTARRS_G': {'file': 'panstarrs-g.txt', 'description': 'PANSTARRS g-band', 'zeropoint': 3909.11, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_R': {'file': 'panstarrs-r.txt', 'description': 'PANSTARRS r-band', 'zeropoint': 3151.44, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_W': {'file': 'panstarrs-w.txt', 'description': 'PANSTARRS w-band', 'zeropoint': 3024.76, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_I': {'file': 'panstarrs-i.txt', 'description': 'PANSTARRS i-band', 'zeropoint': 2584.6, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_Z': {'file': 'panstarrs-z.txt', 'description': 'PANSTARRS z-band', 'zeropoint': 2273.09, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'PANSTARRS_Y': {'file': 'panstarrs-y.txt', 'description': 'PANSTARRS y-band', 'zeropoint': 2205.95, 'method': 'ab', 'rsr': False, 'altname': []}, \
    'SDSS_U': {'file': 'sdss-u.txt', 'description': 'SDSS u-band', 'zeropoint': 1568.5, 'method': 'ab', 'rsr': False, 'altname': ['u']}, \
    'SDSS_G': {'file': 'sdss-g.txt', 'description': 'SDSS g-band', 'zeropoint': 3965.9, 'method': 'ab', 'rsr': False, 'altname': ['g']}, \
    'SDSS_R': {'file': 'sdss-r.txt', 'description': 'SDSS r-band', 'zeropoint': 3162.0, 'method': 'ab', 'rsr': False, 'altname': ['r']}, \
    'SDSS_I': {'file': 'sdss-i.txt', 'description': 'SDSS i-band', 'zeropoint': 2602.0, 'method': 'ab', 'rsr': False, 'altname': ['i']}, \
    'SDSS_Z': {'file': 'sdss-z.txt', 'description': 'SDSS z-band', 'zeropoint': 2244.7, 'method': 'ab', 'rsr': False, 'altname': ['z']}, \
    'TESS': {'file': 'TESS.txt', 'description': 'TESS bandpass', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_Z': {'file': 'ukidss-z.txt', 'description': 'UKIDSS Z-band', 'zeropoint': 2261.4, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_Y': {'file': 'ukidss-y.txt', 'description': 'UKIDSS Y-band', 'zeropoint': 2057.2, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_J': {'file': 'ukidss-j.txt', 'description': 'UKIDSS J-band', 'zeropoint': 1556.8, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_H': {'file': 'ukidss-h.txt', 'description': 'UKIDSS H-band', 'zeropoint': 1038.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'UKIDSS_K': {'file': 'ukidss-k.txt', 'description': 'UKIDSS K-band', 'zeropoint': 644.1, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_Z': {'file': 'vista_z.txt', 'description': 'VISTA Z-band', 'zeropoint': 2263.81, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_Y': {'file': 'vista_y.txt', 'description': 'VISTA Y-band', 'zeropoint': 2087.32, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_J': {'file': 'vista_j.txt', 'description': 'VISTA J-band', 'zeropoint': 1554.03, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_H': {'file': 'vista_h.txt', 'description': 'VISTA H-band', 'zeropoint': 1030.40, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'VISTA_KS': {'file': 'vista_ks.txt', 'description': 'VISTA Ks-band', 'zeropoint': 674.83, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F098M': {'file': 'HST-WFC3_IR_F098M.txt', 'description': 'WFC3 F098M', 'zeropoint': 2154.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F105W': {'file': 'HST-WFC3_IR_F105W.txt', 'description': 'WFC3 F105W', 'zeropoint': 1975.2, 'method': 'vega', 'rsr': False, 'altname': ['wfc3 y']}, \
    'WFC3_F110W': {'file': 'HST-WFC3_IR_F110W.txt', 'description': 'WFC3 F110W', 'zeropoint': 1738.4, 'method': 'vega', 'rsr': False, 'altname': ['wfc3 yj']}, \
    'WFC3_F125W': {'file': 'HST-WFC3_IR_F125W.txt', 'description': 'WFC3 F125W', 'zeropoint': 1564.3, 'method': 'vega', 'rsr': False, 'altname': ['wfc3 j']}, \
    'WFC3_F126N': {'file': 'HST-WFC3_IR_F126N.txt', 'description': 'WFC3 F126N', 'zeropoint': 1552.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F127M': {'file': 'HST-WFC3_IR_F127M.txt', 'description': 'WFC3 F127M', 'zeropoint': 1496.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F128N': {'file': 'HST-WFC3_IR_F128N.txt', 'description': 'WFC3 F128N', 'zeropoint': 1392.6, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F130N': {'file': 'HST-WFC3_IR_F130N.txt', 'description': 'WFC3 F130N', 'zeropoint': 1475.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F132N': {'file': 'HST-WFC3_IR_F132N.txt', 'description': 'WFC3 F132N', 'zeropoint': 1466.6, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F139M': {'file': 'HST-WFC3_IR_F139M.txt', 'description': 'WFC3 F139M', 'zeropoint': 1342.8, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F140W': {'file': 'HST-WFC3_IR_F140W.txt', 'description': 'WFC3 F140W', 'zeropoint': 1324.8, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F153M': {'file': 'HST-WFC3_IR_F153M.txt', 'description': 'WFC3 F153M', 'zeropoint': 1142.0, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F160W': {'file': 'HST-WFC3_IR_F160W.txt', 'description': 'WFC3 F160W', 'zeropoint': 1138.1, 'method': 'vega', 'rsr': False, 'altname': ['wfc3 h']}, \
    'WFC3_F164N': {'file': 'HST-WFC3_IR_F164N.txt', 'description': 'WFC3 F164N', 'zeropoint': 1005.5, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFC3_F167N': {'file': 'HST-WFC3_IR_F167N.txt', 'description': 'WFC3 F167N', 'zeropoint': 1030.0, 'method': 'vega', 'rsr': False, 'altname': []}, \
#    'WFC3_F127M': {'file': 'wfc3_F127M.txt', 'description': 'WFC3 F127M', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
#    'WFC3_F139M': {'file': 'wfc3_F139M.txt', 'description': 'WFC3 F139M', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
#    'WFC3_F164N': {'file': 'wfc3_F164N.txt', 'description': 'WFC3 F164N', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
#    'WFC3_F167N': {'file': 'wfc3_F167N.txt', 'description': 'WFC3 F167N', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_Z': {'file': 'wfcam-z.txt', 'description': 'UKIRT WFCAM Z', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_Y': {'file': 'wfcam-y.txt', 'description': 'UKIRT WFCAM Y', 'zeropoint': 2040.9, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_J': {'file': 'wfcam-j.txt', 'description': 'UKIRT WFCAM J', 'zeropoint': 1548.7, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_H': {'file': 'wfcam-h.txt', 'description': 'UKIRT WFCAM H', 'zeropoint': 1027.1, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_H2': {'file': 'wfcam-h2.txt', 'description': 'UKIRT WFCAM H2', 'zeropoint': 677.1, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WFCAM_BRG': {'file': 'wfcam-brg.txt', 'description': 'UKIRT WFCAM Brackett Gamma', 'zeropoint': 645.5, 'method': 'vega', 'rsr': False, 'altname': ['wfcam brackett gamma']}, \
    'WFCAM_K': {'file': 'wfcam-k.txt', 'description': 'UKIRT WFCAM K', 'zeropoint': 630.0, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_J': {'file': 'wirc_jcont.txt', 'description': 'WIRC J-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_H': {'file': 'wirc_hcont.txt', 'description': 'WIRC H-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_K': {'file': 'wirc_kcont.txt', 'description': 'WIRC K-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_CO': {'file': 'wirc_co.txt', 'description': 'WIRC CO', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_CH4S': {'file': 'wirc_ch4s.txt', 'description': 'WIRC CH4S', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_CH4L': {'file': 'wirc_ch4l.txt', 'description': 'WIRC CH4L', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_FE2': {'file': 'wirc_feii.txt', 'description': 'WIRC Fe II', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRC_BRGAMMA': {'file': 'wirc_brgamma.txt', 'description': 'WIRC H I Brackett Gamma', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': ['wirc brackett gamma']}, \
    'WIRC_PABETA': {'file': 'wirc_pabeta.txt', 'description': 'WIRC H I Paschen Beta', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altname': ['wirc paschen beta']}, \
    'WIRCAM_Y': {'file': 'wircam-cfht-y.txt', 'description': 'CFHT WIRCAM Y', 'zeropoint': 2073.32, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_J': {'file': 'wircam-cfht-j.txt', 'description': 'CFHT WIRCAM J', 'zeropoint': 1551.01, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_H': {'file': 'wircam-cfht-h.txt', 'description': 'CFHT WIRCAM H', 'zeropoint': 1044.35, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_KS': {'file': 'wircam-cfht-ks.txt', 'description': 'CFHT WIRCAM Ks', 'zeropoint': 674.62, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_KCONT': {'file': 'wircam-cfht-kcont.txt', 'description': 'CFHT WIRCAM K-cont', 'zeropoint': 636.17, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_CH4_OFF': {'file': 'wircam-cfht-ch4s.txt', 'description': 'CFHT WIRCAM CH4-off', 'zeropoint': 987.39, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WIRCAM_CH4_ON': {'file': 'wircam-cfht-ch4l.txt', 'description': 'CFHT WIRCAM CH4-on', 'zeropoint': 1076.31, 'method': 'vega', 'rsr': False, 'altname': []}, \
    'WISE_W1': {'file': 'wise_w1.txt', 'description': 'WISE W1 (3.5 micron)', 'zeropoint': 309.54, 'method': 'vega', 'rsr': True, 'altname': ['W1']}, \
    'WISE_W2': {'file': 'wise_w2.txt', 'description': 'WISE W2 (4.6 micron)', 'zeropoint': 171.79, 'method': 'vega', 'rsr': True, 'altname': ['W2']}, \
    'WISE_W3': {'file': 'wise_w3.txt', 'description': 'WISE W3 (13 micron)', 'zeropoint': 31.67, 'method': 'vega', 'rsr': True, 'altname': ['W3']}, \
    'WISE_W4': {'file': 'wise_w4.txt', 'description': 'WISE W4 (22 micron)', 'zeropoint': 8.363, 'method': 'vega', 'rsr': True, 'altname': ['W4']} \
    }
VEGAFILE = 'vega_kurucz.txt'

# some data formats (for future expansion)
INSTRUMENT_DEFAULT_PARAMETERS = {
    'instrument_name': {'altname': ['name','instrument','inst'], 'type': str, 'default': 'UNKNOWN'},
    'altname': {'altname': ['name','instrument','inst'], 'type': str, 'default': []},
    'bibcode': {'altname': ['name','instrument','inst'], 'type': str, 'default': ''},
    'resolution': {'altname': ['name','instrument','inst'], 'type': float, 'default': numpy.nan},
    'wunit': {'altname': ['name','instrument','inst'], 'type': u.quantity.Quantity, 'default': u.micron},
    'funit': {'altname': ['name','instrument','inst'], 'type': u.quantity.Quantity, 'default': u.erg/u.s/u.cm/u.cm/u.micron}
    }

INSTRUMENTS = {
#	'SPEX': {'instrument_name': 'SpeX prism', 'pixelscale': 0.15*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': ['']},
#    'UNKNOWN': {'instrument_name': 'UNKNOWN', 'pixelscale': 1.*u.arcsec, 'slitwidth': 1.*u.arcsec, 'altname': ['UNK'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'RAW': {'instrument_name': 'RAW', 'pixelscale': 1.*u.arcsec, 'slitwidth': 1.*u.arcsec, 'altname': [], 'wunit': u.micron},
    'SED': {'instrument_name': 'SED', 'pixelscale': 1.*u.arcsec, 'wave_range': [0.1,100]*u.micron, 'slitwidth': 1.*u.arcsec, 'altname': ['SPECTRAL_ENERGY_DISTRIBUTION'], 'resolution': 100, 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'APOGEE': {'instrument_name': 'SDSS APOGEE', 'pixelscale': 2./3.*u.arcsec, 'wave_range': [1.51,1.70]*u.micron, 'slitwidth': 2.*u.arcsec, 'resolution': 22500, 'norders': 1, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['APO'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2010SPIE.7735E..1CW'},
    'BOSS': {'instrument_name': 'SDSS BOSS', 'pixelscale': 2./3.*u.arcsec, 'wave_range': [3700,10400]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 2000, 'norders': 1, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['SDSS','BOSS','EBOSS'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'bibcode': '2013AJ....146...32S'},
    'DEIMOS': {'instrument_name': 'Keck DEIMOS', 'pixelscale': 0.1185*u.arcsec, 'wave_range': [5000,9700]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'disperser': '600 l/mm', 'resolution': 2000, 'norders': 1, 'readnoise': 2.5, 'darkcurrent': 0., 'gain': 1.2, 'altname': ['DEIMOS'], 'bibcode': '2003spie.4841.1657F', 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'FIRE': {'instrument_name': 'Magellan FIRE', 'pixelscale': 0.18*u.arcsec, 'wave_range': [0.82,2.51]*u.micron, 'slitwidth': 0.6*u.arcsec, 'resolution': 6000, 'norders': 21, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altname': ['FIRE'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2013PASP..125..270S'},
    'IRS-SL': {'instrument_name': 'Spitzer IRS Short-Low', 'pixelscale': 1.8*u.arcsec, 'wave_range': [5.2,14.5]*u.micron, 'slitwidth': 1.8*u.arcsec, 'resolution': 100, 'norders': 1, 'readnoise': 30., 'darkcurrent': 10., 'gain': 4.6, 'altname': ['IRS','IRS Short Low'], 'wunit': u.micron, 'funit': u.Jy, 'bibcode': '2004ApJS..154...18H'},
    'IRS-LL': {'instrument_name': 'Spitzer IRS Short-Low', 'pixelscale': 5.1*u.arcsec, 'wave_range': [14,38.]*u.micron, 'slitwidth': 5.1*u.arcsec, 'resolution': 90, 'norders': 1, 'readnoise': 30., 'darkcurrent': 40., 'gain': 4.6, 'altname': ['IRS Long Low'], 'wunit': u.micron, 'funit': u.Jy, 'bibcode': '2004ApJS..154...18H'},
    'KAST-RED': {'instrument_name': 'Lick KAST Red Channel', 'pixelscale': 0.43*u.arcsec, 'disperser': '600/7500', 'wave_range': [5700,9200]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 1200, 'norders': 1, 'readnoise': 3.8, 'darkcurrent': 0., 'gain': 1.9, 'altname': ['KAST','KAST-R'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'KAST-BLUE': {'instrument_name': 'Lick KAST Blue Channel', 'pixelscale': 0.43*u.arcsec, 'disperser': '600/4310', 'wave_range': [3300,5520]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 950, 'norders': 1, 'readnoise': 3.7, 'darkcurrent': 0., 'gain': 1.2, 'altname': ['KAST-B'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'LDSS-3': {'instrument_name': 'Magellan LDSS-3', 'pixelscale': 0.189*u.arcsec, 'disperser': 'VPH-RED', 'wave_range': [6000,10000]*u.Angstrom, 'slitwidth': 0.75*u.arcsec, 'resolution': 1810, 'norders': 1, 'readnoise': 4.07, 'darkcurrent': 0., 'gain': 1, 'altname': ['LDSS3'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'bibcode': '1994PASP..106..983A'},
    'LRIS-RED': {'altname': ['LRIS','LRISR'], 'instrument_name': 'Keck LRIS red channel longslit', 'bibcode': '1995PASP..107..375O', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4.6, 'darkcurrent': 0., 'gain': 1.2, 'disperser': '400/8500', 'wave_range': [6300,10100]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1200, 'norders': 1, 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'bibcode': '1995PASP..107..375O'},
    'LRIS-BLUE': {'altname': ['LRISB'], 'instrument_name': 'Keck LRIS blue channel longslit', 'bibcode': '1998SPIE.3355...81M', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4., 'darkcurrent': 0., 'gain': 1.6, 'disperser': '400/3400', 'wave_range': [1270,5740]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 530, 'norders': 1, 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'bibcode': '1995PASP..107..375O'},
    'MAGE': {'instrument_name': 'Magellan MAGE', 'pixelscale': 0.3*u.arcsec, 'wave_range': [3100,10000]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 4100, 'norders': 13, 'readnoise': 2.9, 'darkcurrent': 1.0, 'gain': 1.02, 'altname': ['MagE'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'bibcode': '2008SPIE.7014E..54M'},
    'NIRI-J': {'instrument_name': 'Gemini NIRI J-band (G5202)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [1.05,1.41]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 610, 'norders': 1, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altname': ['NIRI','NIRI G5202'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2003PASP..115.1388H'},
    'NIRI-H': {'instrument_name': 'Gemini NIRI J-band (G5203)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [1.43,1.96]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 825, 'norders': 1, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altname': ['NIRI G5203'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2003PASP..115.1388H'},
    'NIRI-K': {'instrument_name': 'Gemini NIRI J-band (G5204)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [1.90,2.49]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 780, 'norders': 1, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altname': ['NIRI G5204'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2003PASP..115.1388H'},
    'NIRI-L': {'instrument_name': 'Gemini NIRI J-band (G5205)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [2.99,4.15]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 690, 'norders': 1, 'readnoise': 50., 'darkcurrent': 0.25, 'gain': 12.3, 'altname': ['NIRI G5205'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2003PASP..115.1388H'},
    'NIRSPEC': {'instrument_name': 'Keck NIRSPEC', 'pixelscale': 0.43/3.*u.arcsec, 'wave_range': [0.95,5.5]*u.micron, 'slitwidth': 0.43*u.arcsec, 'resolution': 25000, 'norders': 8, 'readnoise': 23., 'darkcurrent': 0.8, 'gain': 5.8, 'altname': ['NIRSPAO'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2000SPIE.4008.1048M'},
    'SPEX-PRISM': {'instrument_name': 'IRTF SpeX prism', 'pixelscale': 0.15*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': ['SPEX','PRISM'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2003PASP..115..362R'},
    'SPEX-SXD': {'instrument_name': 'IRTF SpeX SXD', 'pixelscale': 0.15*u.arcsec, 'wave_range': [0.8,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': ['SXD'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2003PASP..115..362R'},
    'SPEX-LXD1.9': {'instrument_name': 'IRTF SpeX LXD 1.9 micron', 'pixelscale': 0.15*u.arcsec, 'wave_range': [1.95,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': ['SPEX LXD','LXD'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2003PASP..115..362R'},
    'SPEX-LXD2.1': {'instrument_name': 'IRTF SpeX LXD 2.1 micron', 'pixelscale': 0.15*u.arcsec, 'wave_range': [2.15,5.0]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': [], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2003PASP..115..362R'},
    'SPEX-LXD2.3': {'instrument_name': 'IRTF SpeX LXD 2.3 micron', 'pixelscale': 0.15*u.arcsec, 'wave_range': [2.25,5.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altname': [], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'bibcode': '2003PASP..115..362R'},
#	'USPEX': {'instrument_name': 'Updated SpeX prism', 'pixelscale': 0.10*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['']},
#	'USPEX_PRISM': {'instrument_name': 'IRTF Updated SpeX prism', 'pixelscale': 0.10*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['USPEX','UPRISM'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
#    'USPEX_SXD': {'instrument_name': 'IRTF Updated SpeX SXD', 'pixelscale': 0.10*u.arcsec, 'wave_range': [0.7,2.55]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['USXD'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
#    'USPEX_LXD_SHORT': {'instrument_name': 'IRTF Updated SpeX LXD short', 'pixelscale': 0.10*u.arcsec, 'wave_range': [1.67,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500., 'norders': 8, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['ULXD','LXD_SHORT','LXDS'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
#    'USPEX_LXD_LONG': {'instrument_name': 'IRTF Updated SpeX LXD long', 'pixelscale': 0.10*u.arcsec, 'wave_range': [1.98,5.3]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500., 'norders': 7, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altname': ['LXD_LONG','LXDL'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'WFC3-G102': {'instrument_name': 'HST WFC3 IR G102', 'pixelscale': 0.128*u.arcsec, 'wave_range': [0.8,1.15]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 210., 'norders': 7, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altname': ['G102'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'WFC3-G141': {'instrument_name': 'HST WFC3 IR G141', 'pixelscale': 0.128*u.arcsec, 'wave_range': [1.075,1.70]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 130., 'norders': 7, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altname': ['WFC3','HST WFC3','HST WFC3 IR','WFC3 IR','G141'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
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



# spectral model information
#MODEL_PARAMETER_NAMES = ['teff','logg','z','fsed','cld','kzz','slit']
#MODEL_PARAMETERS = {'teff': 1000.0,'logg': 5.0,'z': 0.0,'fsed':'nc','cld':'nc','kzz':'eq','slit':0.5}
#MODEL_PARAMETER_TITLES = {\
#    'teff': '$T_{eff}$',\
#    'logg': '$log\ g$',\
#    'z': '$[M/H]$',\
#    'fsed': '$f_{sed}$',\
#    'cld': '$cld$',\
#    'kzz': '$log\ \kappa_{zz}$',\
#    'slit': '$slit$'}
#MODEL_PARAMETER_UNITS = [u.K,u.cm/u.s/u.s,u.m/u.m,u.m/u.m,u.m/u.m,u.m/u.m,u.arcsec]
#MODEL_PARAMETER_UNITS = {\
#    'teff': u.K, \
#    'logg': u.dex(u.cm/u.s/u.s), \
#    'z': u.dex(), \
#    'fsed': u.m/u.m, \
#    'cld': u.m/u.m, \
#    'kzz': u.m/u.m, \
#    'kzz': u.dex(u.cm*u.cm/u.s), \
#    'slit': u.arcsec}
#DEFINED_MODEL_SET = ['BTSettl2008','burrows06','morley12','morley14','saumon12','drift']
#DEFINED_MODEL_NAMES = {
#    'BTSettl2008': 'BT-Settl (2008)',\
#    'burrows06': 'Burrows et al. (2006)',\
#    'morley12': 'Morley et al. (2012)',\
#    'morley14': 'Morley et al. (2014)',\
#    'saumon12': 'Saumon et al. (2012)',\
#    'drift': 'Witte et al. (2011)'}
SPECTRAL_MODELS = {\
#    'gaia': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/gaia/', 'name': 'AMES GAIA', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altname': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'rawfolder': HOME_FOLDER+'/models/phoenix/nextgen/fullres/', 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
    'btnextgen': {'instruments': {}, 'name': 'BT NextGen', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['nextgen-bt','btnextgen'], 'default': {'teff': 3000., 'logg': 5.0, 'z': 0.0, 'enrich': 0.}}, \
    'btsettl08': {'instruments': {}, 'name': 'BTSettl08', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['allard','allard12','allard2012','btsettl','btsettled','btsettl08','btsettl2008','BTSettl2008'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'enrich': 0.}}, \
    'btsettl15': {'instruments': {}, 'name': 'BTSettl15', 'citation': 'Allard et al. (2015)', 'bibcode': '2015A&A...577A..42B', 'altname': ['allard15','allard2015','btsettl015','btsettl2015','BTSettl2015'],  'default': {'teff': 1500., 'logg': 5.0, 'z': 0.}}, \
    'burrows06': {'instruments': {}, 'name': 'Burrows 2006', 'citation': 'Burrows et al. (2006)', 'bibcode': '2006ApJ...640.1063B', 'altname': ['burrows','burrows2006'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'cld': 'nc'}}, \
    'cond01': {'instruments': {}, 'name': 'AMES Cond', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altname': ['cond','cond-ames','amescond'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'drift': {'instruments': {}, 'name': 'Drift', 'citation': 'Witte et al. (2011)', 'bibcode': '2011A&A...529A..44W', 'altname': ['witte','witte11','witte2011','helling'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
    'dusty01': {'instruments': {}, 'name': 'AMES Dusty', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altname': ['dusty','dusty-ames','amesdusty'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'madhusudhan11': {'instruments': {}, 'name': 'Madhusudhan 2011', 'citation': 'Madhusudhan et al. (2011)', 'bibcode': '2011ApJ...737...34M', 'altname': ['madhu','madhusudhan','madhu11','madhu2011','madhusudhan2011'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.,'cld': 'ae60', 'kzz': 'eq','fsed': 'eq'}}, \
    'morley12': {'instruments': {}, 'name': 'Morley 2012', 'citation': 'Morley et al. (2012)', 'bibcode': '2012ApJ...756..172M', 'altname': ['morley','morley2012'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'fsed': 'f5'}}, \
    'morley14': {'instruments': {}, 'name': 'Morley 2014', 'citation': 'Morley et al. (2014)', 'bibcode': '2014ApJ...787...78M', 'altname': ['morley2014'], 'default': {'teff': 300., 'logg': 5.0, 'z': 0., 'fsed': 'f5', 'cld': 'h50'}}, \
    'nextgen99': {'instruments': {}, 'name': 'Phoenix NextGen', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altname': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
    'saumon12': {'instruments': {}, 'name': 'Saumon 2012', 'citation': 'Saumon et al. (2012)', 'bibcode': '2012ApJ...750...74S', 'altname': ['saumon','saumon2012'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
#    'btcond': {'instruments': {}, 'name': 'BT Cond', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['dusty-cond','bt-cond'], 'rawfolder': '/Volumes/splat/models/btcond/ORIGINAL/', 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.0, 'enrich': 0.0}}, \
#    'btdusty': {'instruments': {}, 'name': 'BT Dusty', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['dusty-bt','bt-dusty'], 'rawfolder': '/Volumes/splat/models/btdusty/ORIGINAL/', 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.0}}, \
}
SPECTRAL_MODEL_PARAMETERS_INORDER = ['teff','logg','z','fsed','cld','kzz','ad','enrich','carbon','oxygen','broad','logpmin','logpmax']
SPECTRAL_MODEL_PARAMETERS = {\
    'teff': {'name': 'temperature', 'prefix': 't', 'unit': u.K, 'default': 1000.0, 'title': '$T_{eff}$ (K)', 'type': 'continuous'}, \
    'logg': {'name': 'gravity', 'prefix': 'g', 'unit': u.dex, 'default': 5.0, 'title': '$\log{g}$ (cgs)', 'type': 'continuous'}, \
    'z': {'name': 'metallicity', 'prefix': 'z', 'unit': u.dex, 'default': 0., 'title': '$[M/H]$', 'type': 'continuous'}, \
    'fsed': {'name': 'rainout', 'prefix': 'f', 'unit': u.m/u.m, 'default': 'nc', 'title': '$f_{sed}$', 'type': 'discrete'}, \
    'cld': {'name': 'cloud', 'prefix': 'c', 'unit': u.m/u.m, 'default': 'nc', 'title': 'Cloud or Condensation Treatment', 'type': 'discrete'}, \
    'kzz': {'name': 'mixing', 'prefix': 'k', 'unit': u.m/u.m, 'default': 'eq', 'title': '$log\ \kappa_{zz}$ (cgs)', 'type': 'discrete'},\
    'ad': {'name': 'adiabat', 'prefix': 'ad', 'unit': u.m/u.m, 'default': 1., 'title': 'Adiabatic Index', 'type': 'continuous'},\
    'enrich': {'name': 'alpha enrichment', 'prefix': 'en', 'unit': u.dex, 'default': 0., 'title': 'Alpha Element Enrichment', 'type': 'continuous'},\
    'carbon': {'name': 'carbon enrichment', 'prefix': 'ca', 'unit': u.dex, 'default': 0., 'title': 'Carbon Enrichment', 'type': 'continuous'},\
    'oxygen': {'name': 'oxygen enrichment', 'prefix': 'ox', 'unit': u.dex, 'default': 0., 'title': 'Oxygen Enrichment', 'type': 'continuous'},\
    'broad': {'name': 'broadening', 'prefix': 'br', 'unit': u.m/u.m, 'default': 'A', 'title': 'Alkali Line Broadening Prescription', 'type': 'discrete'},\
    'logpmin': {'name': 'log pressure top', 'prefix': 'pt', 'unit': u.dex, 'default': -8., 'title': 'log Minimum Pressure (bar)', 'type': 'continuous'},\
    'logpmax': {'name': 'log pressure bottom', 'prefix': 'pb', 'unit': u.dex, 'default': 4., 'title': 'log Maximum Pressure (bar)', 'type': 'continuous'},\
    'radius': {'name': 'radius', 'prefix': 'r', 'unit': u.Rsun, 'default': 0., 'title': 'Radius (R$_{\odot}$)', 'type': 'continuous'},\
    }


# evolutionary model information
EVOLUTIONARY_MODELS = {\
    'baraffe03': {'name': 'Barraffe 2003','citation': 'Baraffe et al. (2003)', 'bibcode': '2003A&A...402..701B', 'altname': ['bar03','baraffe2003']},\
    'baraffe15': {'name': 'Barraffe 2015','citation': 'Baraffe et al. (2015)', 'bibcode': '2015A&A...577A..42B', 'altname': ['baraffe','bar15','baraffe2015']},\
    'burrows01': {'name': 'Burrows 2001','citation': 'Burrows et al. (2001)', 'bibcode': '2015A&A...577A..42B', 'altname': ['burrows','bur01','burrows2001']},\
    'saumon08': {'name': 'Saumon 2008','citation': 'Saumon et al. (2008)', 'bibcode': '2015A&A...577A..42B', 'altname': ['saumon','sau08','saumon2008']}}

EVOLUTIONARY_MODEL_PARAMETERS = {\
    'mass': {'unit': u.solMass, 'default': 0.05, 'title': '$M$'},\
    'age': {'unit': u.Gyr, 'default': 5., 'title': '$\tau$'},\
    'temperature': {'unit': u.K, 'default': 1000.0, 'title': '$T_{eff}$'},\
    'gravity': {'unit': u.dex(u.cm / u.s**2), 'default': 5.0, 'title': '$\log{g}'},\
    'luminosity': {'unit': u.dex(u.solLum), 'default': -5., 'title': '$\log{L_{bol}/L_{\odot}}$'},\
    'radius': {'unit': u.solRad, 'default': 0.1, 'title': '$R_{\odot}$'}}


# Empirical relation contants

INDICES_SETS = {
    'allers': {'bibcode': '2013ApJ...772...79A'},
    'burgasser': {'bibcode': '2006ApJ...637.1067B'},
    'bardalez': {'bibcode': '2014ApJ...794..143B'},
    'geballe': {'bibcode': '2002ApJ...564..466G'},
    'mclean': {'bibcode': '2003ApJ...596..561M'},
    'reid': {'bibcode': '2001AJ....121.1710R'},
    'rojas': {'bibcode': '2012ApJ...748...93R'},
    'slesnick': {'bibcode': '2004ApJ...610.1045S'},
    'testi': {'bibcode': '2001ApJ...552L.147T'},
    'tokunaga': {'bibcode': '1999AJ....117.1010T'},
}    


SPT_TEFF_RELATIONS = {
    'golimowski': {'altname': ['golimowski04','golimowski2004','gol04'], 'reference': 'Golimowski et al. (2004)','bibcode': '2004AJ....127.3516G','sptoffset': 10.,'coeff': [9.5373e-4,-9.8598e-2,4.0323,-8.3099e1,9.0951e2,-5.1287e3,1.4322e4],'range': [16.,38.],'fitunc': 124.},
    'looper': {'altname': ['looper08','looper2008','lop08'], 'reference': 'Looper et al. (2008)','bibcode': '2008ApJ...685.1183L','sptoffset': 20.,'coeff': [9.084e-4,-4.255e-2,6.414e-1,-3.101,1.950,-108.094,2319.92],'range': [20.,38.],'fitunc': 87.},
    'stephens': {'altname': ['stephens09','stephens2009','ste09'], 'reference': 'Stephens et al. (2009)','bibcode': '2009ApJ...702..154S','sptoffset': 10.,'coeff': [-0.0025492,0.17667,-4.4727,54.67,-467.26,4400.],'range': [16.,38.],'fitunc': 100.},
    'stephens-alt': {'altname': ['stephens09-alt','stephens2009-alt','ste09alt'], 'reference': 'Stephens et al. (2009)','bibcode': '2009ApJ...702..154S','sptoffset': 10.,'coeff': [-0.011997,1.2315,-50.472,1031.9,-10560.,44898.],'range': [23.,38.],'fitunc': 100.},
    'marocco': {'altname': ['marocco13','marocco2013','mar13'], 'reference': 'Marocco et al. (2013)','bibcode': '2013AJ....146..161M','sptoffset': 10.,'coeff': [7.4211e-5,-8.43736e-3,3.90319e-1,-9.46896,129.141,-975.953,3561.47,-1613.82],'range': [17.,38.],'fitunc': 140.},
    'filippazzo': {'altname': ['filippazzo15','filippazzo2015','fil15'], 'reference': 'Filippazzo et al. (2015)','bibcode': '2015ApJ...810..158F','sptoffset': 10.,'coeff': [1.546e-4, -1.606e-2, 6.318e-1, -1.191e1, 1.155e2, -7.005e2, 4.747e3], 'range': [16.,39.],'fitunc': 113.},
    'faherty': {'altname': ['faherty16','faherty2016','fah16'], 'reference': 'Faherty et al. (2016)','bibcode': '2016ApJS..225...10F','sptoffset': 10.,'coeff': [1.546e-4,-1.606e-2,6.318e-1,-1.191e1,1.155e2,-7.005e2,4.747e3], 'range': [17.,38.],'fitunc': 113.},
    'faherty-young': {'altname': ['faherty16-young','faherty2016-young','fah16yng'], 'reference': 'Faherty et al. (2016)','bibcode': '2016ApJS..225...10F','sptoffset': 10.,'coeff': [1.330,-6.68637e1,1.23542e3,-1.00688e4,3.27664e4], 'range': [17.,27.],'fitunc': 180.},
    'faherty-young2': {'altname': ['faherty16-young2','faherty2016-young2','fah16yng2'], 'reference': 'Faherty et al. (2016)','bibcode': '2016ApJS..225...10F','sptoffset': 10.,'coeff': [9.106e-4,-1.016e-1,4.578,-1.066e2,1.360e3,-9.183e3,2.795e4], 'range': [17.,27.],'fitunc': 198.},
    'faherty-group': {'altname': ['faherty16-group','faherty2016-group','fah16grp'], 'reference': 'Faherty et al. (2016)','bibcode': '2016ApJS..225...10F','sptoffset': 10.,'coeff': [7.383e0,-3.44522e2,4.87986e3], 'range': [17.,27.],'fitunc': 172.},
    'dupuy-saumon': {'altname': ['dupuy','dupuy17','dupuy2017','dup17','dupuy17-saumon','dupuy2017-saumon','dup17-saumon'], 'reference': 'Dupuy et al. (2017)','bibcode': '','sptoffset': 10.,'coeff': [6.001,-284.52,4544.3], 'range': [21.5,35.],'fitunc': 80.},
    'dupuy-lyon': {'altname': ['dupuy17-lyon','dupuy2017-lyon','dup17-lyon'], 'reference': 'Dupuy et al. (2017)','bibcode': '','sptoffset': 10.,'coeff': [4.582,-238.03,4251.0], 'range': [17.,35.],'fitunc': 90.},
    }


SPT_COLORS_RELATIONS = {
    'skrzypek2015': {'altname': ['skrzypek','skrzypek15'], 'reference': 'Skrzypek et al. (2015)','bibcode': '2015A%26A...574A..78S','range': [15,38],'filters': ['i','z','y','j','h','k','w1','w2'],'scatter': 0.07,
        'values': { 
            'i-z': [0.91,1.45,1.77,1.93,1.99,2.01,2.02,2.04,2.1,2.2,2.33,2.51,2.71,2.93,3.15,3.36,3.55,3.7,3.82,3.9,3.95,3.98,4.01,4.08], \
            'z-y': [0.47,0.6,0.7,0.77,0.82,0.86,0.88,0.9,0.92,0.94,0.97,1.0,1.04,1.09,1.16,1.23,1.33,1.43,1.55,1.68,1.81,1.96,2.11,2.26], \
            'y-j': [0.55,0.67,0.78,0.87,0.96,1.04,1.11,1.18,1.23,1.27,1.31,1.33,1.35,1.21,1.2,1.19,1.19,1.18,1.18,1.17,1.16,1.16,1.15,1.15], \
            'j-h': [0.45,0.53,0.56,0.58,0.6,0.63,0.67,0.73,0.79,0.86,0.91,0.96,0.97,0.96,0.9,0.8,0.65,0.46,0.25,0.02,-0.19,-0.35,-0.43,-0.36], \
            'h-k': [0.32,0.39,0.44,0.47,0.51,0.54,0.58,0.63,0.67,0.71,0.74,0.75,0.75,0.71,0.65,0.56,0.45,0.31,0.16,0.01,-0.11,-0.19,-0.2,-0.09], \
            'k-w1': [0.11,0.22,0.25,0.26,0.27,0.29,0.33,0.4,0.48,0.56,0.65,0.72,0.77,0.79,0.79,0.76,0.71,0.65,0.59,0.55,0.54,0.59,0.7,0.9], \
            'w1-w2': [0.17,0.21,0.24,0.26,0.27,0.27,0.28,0.28,0.29,0.3,0.32,0.36,0.41,0.48,0.57,0.68,0.82,0.99,1.19,1.43,1.7,2.02,2.38,2.79]}
        }
    }


SPT_BC_RELATIONS = {
    'liu2010': {'altname': ['liu10','liu'], 'bibcode': '2010ApJ...722..311L', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.14, 'range' : [16,38.5], 'coeff': [1.462266e-6,-1.558986e-4,6.540717e-3,-1.367605e-1,1.491738e0,-8.053993e0,1.890448e1]},
        'MKO_H': {'fitunc' : 0.07, 'range' : [16,38.5], 'coeff': [1.148133e-6,-1.171595e-4,4.733874e-3,-9.618535e-2,1.027185e0,-5.426683e0,1.366709e1]},
        'MKO_K': {'fitunc' : 0.08, 'range' : [16,38.5], 'coeff': [3.159780e-7,-3.177629e-5,1.282657e-3,-2.647464e-2,2.868014e-1,-1.471358e0,5.795848e0]}}},
    'dupuy2013': {'altname': ['dupuy13','dupuy'], 'bibcode': '2013Sci...341.1492D', 'sptoffset': 10, 'method': 'interpolate', 'filters': {
        'MKO_Y': {'spt': [38,38.5,39,40.,40.5], 'bc': [1.6,1.1,1.0,0.8,-0.7], 'rms': [0.6,0.6,0.6,0.6,0.6]}, 
        'MKO_J': {'spt': [38,38.5,39,39.5,40.,40.5], 'bc': [2.4,2.0,2.0,1.9,0.9,-0.4], 'rms': [0.6,0.6,0.6,0.6,0.6,0.6]}, 
        'MKO_H': {'spt': [38,38.5,39,39.5,40.,40.5], 'bc': [2.1,1.6,1.7,1.5,0.4,-1.5], 'rms': [0.6,0.6,0.6,0.6,0.6,0.6]}}}, 
    'filippazzo2015': {'altname': ['filippazzo','filippazzo-field','filippazzo15','filippazzo15-field','fillippazzo','fillippazzo-field','filipazo','filipazo-field','filippazo','filippazo-field'], 'bibcode': '2015ApJ...810..158F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.163, 'range' : [16,39], 'coeff': [-1.555e-5,1.200e-3,-3.437e-2,4.566e-1,-2.862e0,8.842e0]},
        '2MASS_KS': {'fitunc' : 0.243, 'range' : [16,38], 'coeff': [-5.474e-6,4.549e-4,-1.462e-2,2.204e-1,-1.508e0,6.815e0]}}},
    'filippazzo2015-young': {'altname': ['filippazzo-young','filippazzo15-young','fillippazzo-young','filipazo-young','filippazo-young'], 'bibcode': '2015ApJ...810..158F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.189, 'range' : [17,38], 'coeff': [-1.283e-6,-1.118e-4,1.165e-2,-3.020e-1,2.899e0,-7.352e0]},
        '2MASS_KS': {'fitunc' : 0.126, 'range' : [17,38], 'coeff': [2.742e-6,-2.633e-4,9.711e-3,-1.759e-1,1.565e0,-2.174e0]}}},
        }

SPT_LBOL_RELATIONS = {
    'filippazzo2015': {'altname': ['filippazzo','filippazzo15','fillippazzo','filipazo','filippazo'], 'bibcode': '2015ApJ...810..158F', 'sptoffset': 10., 'method': 'polynomial', 'fitunc' : 0.133, 'range' : [16,39], 'coeff': [2.736e-7,-3.220e-5,1.449e-3,-3.207e-2,3.727e-1,-2.310e0,2.787e0]},
        }

ABSMAG_LBOL_RELATIONS = {
    'dupuy2017': {'altname': ['dupuy17','dupuy'], 'bibcode': '2017ApJS..231...15D',  'method': 'polynomial', 'filters': {
        'MKO_H': {'fitunc' : 0.023, 'range' : [9.6,13.3], 'coeff': [1.06200e-2,-3.51721e-1,3.46876e1,-1.3282e1]},
        'MKO_K': {'fitunc' : 0.05, 'range' : [9.1,17.8], 'coeff': [4.54547e-4,-3.068824e-2,8.162709e-1,-1.0671188e1,6.811147e1,-1.72188e2]},
        '2MASS_H': {'fitunc' : 0.023, 'range' : [9.2,13.3], 'coeff': [8.9222e-3,-2.90797e-1,2.74259e0,-1.0426e1]},
        '2MASS_KS': {'fitunc' : 0.05, 'range' : [8.8,16.6], 'coeff': [9.8821e-5,-7.85837e-3,2.357643e-1,-3.364101e0,2.258776e1,-5.9877e1]}}},
        }

SPT_ABSMAG_RELATIONS = {
    'dahn2002': {'altname': ['dahn','dahn02'], 'bibcode': '2002AJ....124.1170D', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.25, 'range' : [17.,28.], 'coeff': [0.341,8.38]}}},
    'cruz2003': {'altname': ['cruz','cruz03'], 'bibcode': '2003AJ....126.2421C', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.30, 'range' : [16.,28.], 'coeff': [-6.892e-4,3.453e-2,-6.193e-1,5.043,-4.410]}}},
    'burgasser2007': {'altname': ['burgasser','burgasser07'], 'bibcode': '2007ApJ...659..655B', 'sptoffset': 20, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.30, 'range' : [20., 38.], 'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, 
        'MKO_H': {'fitunc' : 0.27, 'range' : [20., 38.], 'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, 
        'MKO_K': {'fitunc' : 0.26, 'range' : [20., 38.], 'coeff': [.0000001051, -.000006985, .0001807, -.002271, .01414, -.04024, .05129, .2322, 10.45]}}},
    'looper2008': {'altname': ['looper','looper08'],'bibcode': '2008ApJ...685.1183L', 'sptoffset': 20, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.29, 'range' : [20., 38.], 'coeff': [-5.462e-6,2.595e-4,-3.915e-3,1.663e-2,3.690e-2,1.255e-1,11.817]},
        '2MASS_H': {'fitunc' : 0.29, 'range' : [20., 38.], 'coeff' : [-4.218e-6,1.987e-4,-2.970e-3,1.261e-2,3.032e-2,1.125e-1,11.010]},
        '2MASS_KS': {'fitunc' : 0.33, 'range' : [20., 38.], 'coeff' : [-4.104e-6,1.911e-4,-2.864e-3,1.299e-2,2.565e-2,7.369e-2,10.521]}}},
    'faherty2012': {'altname': ['faherty12'],'bibcode': '2012ApJ...752...56F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.30, 'range' : [20., 38.], 'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, 
        'MKO_H': {'fitunc' : 0.27, 'range' : [20., 38.], 'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, 
        'MKO_K': {'fitunc' : 0.28, 'range' : [20., 38.], 'coeff' : [.0000816516, -.00469032, .0940816, -.485519, 9.76100]}}},
    'dupuy2012': {'altname': ['dupuy','dupuy12'], 'bibcode': '2012ApJS..201...19D', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
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
        'WISE_W4': {'fitunc': 0.76, 'range':[16., 39.], 'coeff': [-2.16042e-3,1.14630e-1,7.78974e0]}}},
    'dupuy2013': {'altname': ['dupuy13'], 'bibcode': '2013Sci...341.1492D', 'sptoffset': 0, 'method': 'interpolate', 'filters': {
        'MKO_Y': {'spt': [38,38.5,39,40.], 'absmag': [17.4,18.81,19.26,20.24], 'rms': [0.25,0.51,0.88,0.17]}, 
        'MKO_J': {'spt': [38,38.5,39,39.5,40.], 'absmag': [16.43,17.87,18.39,17.68,20.09], 'rms': [0.46,0.44,0.95,0.37,0.25]}, 
        'MKO_H': {'spt': [38,38.5,39,39.5,40.], 'absmag': [16.82,18.2,18.77,18.08,20.6], 'rms': [0.43,0.45,1.08,0.39,0.25]}, 
        'MKO_K': {'spt': [38,38.5,39,40.], 'absmag': [16.93,18.27,18.89,20.7], 'rms': [0.8,0.4,0.57,0.18]}, 
        'IRAC_CH1': {'spt': [38,38.5,39,39.5,40.], 'absmag': [15.11,15.83,16.17,15.58,16.99], 'rms': [0.15,0.22,0.23,0.41,0.21]}, 
        'IRAC_CH2': {'spt': [38,38.5,39,39.5,40.], 'absmag': [13.4,13.79,14.09,13.51,14.66], 'rms': [0.21,0.12,0.2,0.43,0.28]}}}, 
    'tinney2003': {'altname': ['tinney','tinney03'],'bibcode': '2003AJ....126..975T', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'COUSINS_I': {'fitunc' : 0.37, 'range' : [20., 37.5], 'coeff': [-2.49821e-6,1.04398e-3,-6.49719e-2,1.56038,-1.58296e1,7.22089e1]},
        'UKIRT_Z': {'fitunc' : 0.29, 'range' : [20., 37.5], 'coeff': [-9.97226e-7,1.05950e-4,-4.57019e-3,1.02898e-1,-1.29357e0,8.96822e0,-3.08010e1,4.99447e1]},
        'UKIRT_K': {'fitunc' : 0.40, 'range' : [20., 37.5], 'coeff': [1.14139e-5,-8.86885e-4,2.68071e-2,-3.89554e-1,2.95440e0,8.14626e-1]},
        '2MASS_KS': {'fitunc' : 0.38, 'range' : [20., 37.5], 'coeff': [-1.25074e-5,1.63124e-3,-7.42418e-2,1.54509,-1.47407e1,6.27861e1]},
        'UKIRT_J': {'fitunc' : 0.30, 'range' : [20., 37.5], 'coeff': [-9.91110e-7,1.05811e-4,-4.58399e-3,1.03572e-1,-1.30526e0,9.06701e0,-3.13411e1,5.04642e1]},
        '2MASS_J': {'fitunc' : 0.36, 'range' : [20., 37.5], 'coeff': [-2.80824e-6,3.41146e-4,-1.73848e-2,4.82120e-1,-7.86911,7.57222e1,-3.98105e2,8.94012e2]}}},
    'tinney2014': {'altname': ['tinney14'],'bibcode': '2014ApJ...796...39T', 'sptoffset': 0, 'method': 'interpolate', 'filters': {
        'MKO_J': {'spt': [36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,42], 'absmag': [15.22,15.49,16.39,16.66,17.9,18.35,19.08,20.32,22.39,22.18,25.76], 'rms': [0.31,0.37,0.72,0.36,0.46,0.9,0.97,1.25,1.,0.76,3.52]}, 
        'WISE_W2': {'spt': [36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,42], 'absmag': [12.86,13.28,13.39,13.44,13.75,13.92,14.28,14.65,15.2,14.78,15.76], 'rms': [0.17,0.48,0.27,0.23,0.22,0.24,0.46,0.35,1.,0.77,2.15]}}}, 
    'filippazzo2015': {'altname': ['filippazzo','filippazzo15','fillippazzo','filipazo','filippazo'],'bibcode': '2015ApJ...810..158F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [3.478e-5, -2.684e-3, 7.771e-2, -1.058, 7.157, -8.350]}, 
        'WISE_W2': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [8.190e-6, -6.938e-4, 2.283e-2, -3.655e-1, 3.032, -5.043e-1]}}},
    'faherty2016': {'altname': ['faherty','faherty2016','faherty-field'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.402, 'range' : [16., 39.], 'coeff': [3.478e-5,-2.684e-3,7.771e-2,-1.058,7.157,-8.350]},
        '2MASS_H': {'fitunc' : 0.389, 'range' : [16., 39.], 'coeff' : [2.841e-5,-2.217e-3,6.551e-2,-9.174e-1,6.406,-7.496]},
        '2MASS_KS': {'fitunc' : 0.537, 'range' : [16., 39.], 'coeff' : [2.540e-5,-1.997e-3,5.978e-2,-8.481e-1,5.970,-6.704]},
        'WISE_W1': {'fitunc' : 0.365, 'range' : [16., 39.], 'coeff' : [8.337e-6,-6.897e-4,2.258e-2,-3.603e-1,2.991,-1.664e-1]},
        'WISE_W2': {'fitunc' : 0.398, 'range' : [16., 39.], 'coeff' : [8.190e-6,-6.938e-4,2.283e-2,-3.655e-1,3.032,-5.043e-1]},
        'WISE_W3': {'fitunc' : 0.446, 'range' : [16., 39.], 'coeff' : [-1.024e-6,9.477e-5,-2.573e-3,1.520e-2,3.365e-1,6.462]}}},
    'faherty2016-young': {'altname': ['faherty-young'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.647, 'range' : [17., 27.], 'coeff': [4.032e-3,-1.416e-1,2.097,8.478e-1]},
        '2MASS_H': {'fitunc' : 0.634, 'range' : [17., 27.], 'coeff' : [2.642e-3,-1.049e-1,1.753,1.207]},
        '2MASS_KS': {'fitunc' : 0.640, 'range' : [17., 27.], 'coeff' : [-1.585e-2,7.338e-1,4.537]},
        'WISE_W1': {'fitunc' : 0.648, 'range' : [17., 27.], 'coeff' : [-1.397e-2,5.955e-1,5.247]},
        'WISE_W2': {'fitunc' : 0.694, 'range' : [17., 27.], 'coeff' : [-1.507e-2,5.944e-1,5.061]},
        'WISE_W3': {'fitunc' : 0.717, 'range' : [17., 27.], 'coeff' : [-1.003e-4,-1.670e-3,2.023e-1,7.529]}}},
    'faherty2016-group': {'altname': ['faherty-group'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.660, 'range' : [17., 27.], 'coeff': [-3.825e-3,1.370e-1,-9.279e-1,10.141]},
        '2MASS_H': {'fitunc' : 0.603, 'range' : [17., 27.], 'coeff' : [-3.909e-3,1.346e-1,-9.347e-1,9.728]},
        '2MASS_KS': {'fitunc' : 0.556, 'range' : [17., 27.], 'coeff' : [-4.006e-3,1.378e-1,-1.031,9.916]},
        'WISE_W1': {'fitunc' : 0.551, 'range' : [17., 27.], 'coeff' : [-4.483e-3,1.505e-1,-1.208,10.403]},
        'WISE_W2': {'fitunc' : 0.616, 'range' : [17., 27.], 'coeff' : [-6.821e-3,2.322e-1,-2.133,13.322]},
        'WISE_W3': {'fitunc' : 0.427, 'range' : [17., 27.], 'coeff' : [-5.684e-3,1.993e-1,-1.987,13.972]}}},
    'liu2016-ir': {'altname': ['liu','liu-ir','liu2016','liu2016-ir','liu16-ir','liu16-field','liu16-field-ir','liu-field','liu-field-ir'],'bibcode': '2016ApJ...833...96L', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.35, 'range' : [16., 28.], 'coeff': [0.402,8.437]},
        'MKO_J': {'fitunc' : 0.37, 'range' : [16., 28.], 'coeff': [0.368,8.062]},
        'MKO_H': {'fitunc' : 0.33, 'range' : [16., 28.], 'coeff': [0.348,7.516]},
        'MKO_K': {'fitunc' : 0.31, 'range' : [16., 28.], 'coeff': [0.319,7.258]},
        '2MASS_J': {'fitunc' : 0.38, 'range' : [16., 28.], 'coeff': [0.391,7.848]},
        '2MASS_H': {'fitunc' : 0.36, 'range' : [16., 28.], 'coeff': [0.340,7.555]},
        '2MASS_KS': {'fitunc' : 0.33, 'range' : [16., 28.], 'coeff': [0.317,7.292]},
        'WISE_W1': {'fitunc' : 0.26, 'range' : [16., 28.], 'coeff': [0.243,7.698]},
        'WISE_W2': {'fitunc' : 0.24, 'range' : [16., 28.], 'coeff': [0.236,7.466]}}},
    'liu2016-optical': {'altname': ['liu-optical','liu16-optical','liu16-field-optical','liu-field-optical'],'bibcode': '2016ApJ...833...96L', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.29, 'range' : [16., 28.], 'coeff': [0.404,8.677]},
        'MKO_J': {'fitunc' : 0.26, 'range' : [16., 28.], 'coeff': [0.364,8.131]},
        'MKO_H': {'fitunc' : 0.24, 'range' : [16., 28.], 'coeff': [0.331,7.820]},
        'MKO_K': {'fitunc' : 0.25, 'range' : [16., 28.], 'coeff': [0.312,7.433]},
        '2MASS_J': {'fitunc' : 0.30, 'range' : [16., 28.], 'coeff': [0.380,8.021]},
        '2MASS_H': {'fitunc' : 0.30, 'range' : [16., 28.], 'coeff': [0.332,7.723]},
        '2MASS_KS': {'fitunc' : 0.29, 'range' : [16., 28.], 'coeff': [0.308,7.491]},
        'WISE_W1': {'fitunc' : 0.29, 'range' : [16., 28.], 'coeff': [0.255,7.610]},
        'WISE_W2': {'fitunc' : 0.29, 'range' : [16., 28.], 'coeff': [0.231,7.559]}}},
    'liu2016-vlg': {'altname': ['liu-vlg','liu16-vlg'],'bibcode': '2016ApJ...833...96L', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.77, 'range' : [16., 27.], 'coeff': [0.799,3.689]},
        'MKO_J': {'fitunc' : 0.72, 'range' : [16., 27.], 'coeff': [0.731,3.475]},
        'MKO_H': {'fitunc' : 0.64, 'range' : [16., 27.], 'coeff': [0.633,3.778]},
        'MKO_K': {'fitunc' : 0.56, 'range' : [16., 27.], 'coeff': [0.553,3.898]},
        '2MASS_J': {'fitunc' : 0.72, 'range' : [16., 27.], 'coeff': [0.745,3.406]},
        '2MASS_H': {'fitunc' : 0.64, 'range' : [16., 27.], 'coeff': [0.632,3.734]},
        '2MASS_KS': {'fitunc' : 0.55, 'range' : [16., 27.], 'coeff': [0.573,3.699]},
        'WISE_W1': {'fitunc' : 0.46, 'range' : [16., 27.], 'coeff': [0.457,4.392]},
        'WISE_W2': {'fitunc' : 0.48, 'range' : [16., 27.], 'coeff': [0.451,4.039]}}},
    'liu2016-intg': {'altname': ['liu-intg','liu16-intg'],'bibcode': '2016ApJ...833...96L', 'sptoffset': 10., 'method': 'polynomial', 'filters': {
        'MKO_Y': {'fitunc' : 0.36, 'range' : [20., 27.], 'coeff': [0.562,6.778]},
        'MKO_J': {'fitunc' : 0.50, 'range' : [20., 27.], 'coeff': [0.537,5.941]},
        'MKO_H': {'fitunc' : 0.43, 'range' : [20., 27.], 'coeff': [0.415,6.557]},
        'MKO_K': {'fitunc' : 0.39, 'range' : [20., 27.], 'coeff': [0.355,6.555]},
        '2MASS_J': {'fitunc' : 0.50, 'range' : [20., 27.], 'coeff': [0.545,5.914]},
        '2MASS_H': {'fitunc' : 0.47, 'range' : [20., 27.], 'coeff': [0.411,6.552]},
        '2MASS_KS': {'fitunc' : 0.36, 'range' : [20., 27.], 'coeff': [0.356,6.573]},
        'WISE_W1': {'fitunc' : 0.13, 'range' : [20., 27.], 'coeff': [0.237,7.436]},
        'WISE_W2': {'fitunc' : 0.18, 'range' : [20., 27.], 'coeff': [0.193,7.561]}}},
}

# telscope locations
TELESCOPES = {
    'KECK': {'lat': 19.8263*u.deg, 'lon': -155.4783*u.deg, 'height': 4160*u.m, 'altname': ['MAUNA_KEA','SUBARU','IRTF','CHFT','UKIRT','GEMINI_N','GEMINI_NORTH']},
    'VLT': {'lat': -24.6275*u.deg, 'lon': -70.4044*u.deg, 'height': 2636*u.m, 'altname': ['PARANAL','CERRO_PARANAL','VERY_LARGE_TELESCOPE','ESA']},
}

