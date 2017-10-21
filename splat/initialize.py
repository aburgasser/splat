# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
.. note::
         These commands initialize the SPLAT code 
"""

import os
import sys
from astropy import units as u


# things that are constants
VERSION = '0.31'
__version__ = VERSION
SPLAT_URL = 'http://splat.physics.ucsd.edu/splat/'
DOCUMENTATION_URL = 'http://pono.ucsd.edu/~adam/splat/'
GITHUB_URL = 'https://github.com/aburgasser/splat/'
BIBCODE = '2017arXiv170700062B'
EMAIL = 'aburgasser@gmail.com'
DB_SOURCES_FILE = 'source_data.txt'
DB_SPECTRA_FILE = 'spectral_data.txt'
DB_PHOTOMETRY_FILE = 'photometry_data.txt'
BIBFILE = 'splat_bibs.bib'
TMPFILENAME = 'splattmpfile'
HOME_FOLDER = os.path.expanduser('~')
DATA_FOLDER = '/reference/Spectra/'
FILTER_FOLDER = '/reference/Filters/'
SPECTRAL_MODEL_FOLDER = '/reference/SpectralModels/'
EVOLUTIONARY_MODEL_FOLDER = '/reference/EvolutionaryModels/'
DOCS_FOLDER = '/docs/'
DOCS_INDEX_HTML = '/docs/_build/html/index.html'
WEB_HTML_BASE = '/docs/_templates/'
DB_FOLDER = '/db/'
ACCESS_FILE = '.splat_access'
EXTERNAL_SPECTRAL_MODELS_FILE = '.splat_spectral_models'
EXTERNAL_EVOLUTIONARY_MODELS_FILE = '.splat_evolutionary_models'
EXTERNAL_DATA_FILE = '.splat_data'
DEFAULT_WAVE_UNIT = u.micron
DEFAULT_FLUX_UNIT = u.erg/u.s/u.cm/u.cm/u.micron
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

#set user SPLAT model path from environmental variable
SPLAT_USER_DATA = './'
if os.environ.get('SPLAT_DATA') != None:
    SPLAT_USER_DATA = os.environ['SPLAT_DATA']

# Unit standards
BASE_WAVE_UNIT = u.micron
BASE_FLUX_UNIT = u.erg/u.s/u.cm/u.cm/u.micron
BASE_SED_UNIT = u.erg/u.s/u.cm/u.cm


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
    '2MASS_J': {'file': 'j_2mass.txt', 'description': '2MASS J-band', 'zeropoint': 1594.0, 'method': 'vega', 'rsr': True, 'altnames': []}, \
    '2MASS_H': {'file': 'h_2mass.txt', 'description': '2MASS H-band', 'zeropoint': 1024.0, 'method': 'vega', 'rsr': True, 'altnames': []}, \
    '2MASS_KS': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega', 'rsr': True, 'altnames': ['2MASS_K']}, \
#    '2MASS_K': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega'}, \
#    '2MASS_Ks': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega'}, \
    'BESSEL_I': {'file': 'bessel_i.txt', 'description': 'Bessel I-band', 'zeropoint': 2405.3, 'method': 'vega', 'rsr': False, 'altnames': ['I']}, \
    'COUSINS_I': {'file': 'i_cousins.txt', 'description': 'Cousins I-band', 'zeropoint': 2405.3, 'method': 'vega', 'rsr': False, 'altnames': ['IC']}, \
    'FOURSTAR_J': {'file': 'fourstar-j.txt', 'description': 'FOURSTAR J-band', 'zeropoint': 1581.2, 'method': 'vega', 'rsr': False, 'altnames': ['4star j']}, \
    'FOURSTAR_J1': {'file': 'fourstar-j1.txt', 'description': 'FOURSTAR J1-band', 'zeropoint': 1978.7, 'method': 'vega', 'rsr': False, 'altnames': ['4star j1']}, \
    'FOURSTAR_J2': {'file': 'fourstar-j2.txt', 'description': 'FOURSTAR J2-band', 'zeropoint': 1774.5, 'method': 'vega', 'rsr': False, 'altnames': ['4star j2']}, \
    'FOURSTAR_J3': {'file': 'fourstar-j3.txt', 'description': 'FOURSTAR J3-band', 'zeropoint': 1488.8, 'method': 'vega', 'rsr': False, 'altnames': ['4star j3']}, \
    'FOURSTAR_H': {'file': 'fourstar-h.txt', 'description': 'FOURSTAR H-band', 'zeropoint': 1054.9, 'method': 'vega', 'rsr': False, 'altnames': ['4star h']}, \
    'FOURSTAR_H_SHORT': {'file': 'fourstar-hshort.txt', 'description': 'FOURSTAR H short', 'zeropoint': 1119.1, 'method': 'vega', 'rsr': False, 'altnames': ['4star h short','4star h-short','4star hs','fourstar hs']}, \
    'FOURSTAR_H_LONG': {'file': 'fourstar-hlong.txt', 'description': 'FOURSTAR H long', 'zeropoint': 980.7, 'method': 'vega', 'rsr': False, 'altnames': ['4star h long','4star h-long','4star hl','fourstar hl']}, \
    'FOURSTAR_KS': {'file': 'fourstar-ks.txt', 'description': 'FOURSTAR Ks-band', 'zeropoint': 675.7, 'method': 'vega', 'rsr': False, 'altnames': ['4star k','4star ks','fourstar k']}, \
    'HAWK_Y': {'file': 'hawk-y.txt', 'description': 'HAWK Y-band', 'zeropoint': 2092.9, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'HAWK_J': {'file': 'hawk-j.txt', 'description': 'HAWK J-band', 'zeropoint': 1543.5, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'HAWK_H': {'file': 'hawk-h.txt', 'description': 'HAWK H-band', 'zeropoint': 1053.6, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'HAWK_H2': {'file': 'hawk-h2.txt', 'description': 'HAWK H2-band', 'zeropoint': 688.8, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'HAWK_CH4': {'file': 'hawk-ch4.txt', 'description': 'HAWK CH4-band', 'zeropoint': 1093.4, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'HAWK_KS': {'file': 'hawk-ks.txt', 'description': 'HAWK Ks-band', 'zeropoint': 675.3, 'method': 'vega', 'rsr': False, 'altnames': ['hawk k']}, \
    'HAWK_BRG': {'file': 'hawk-brg.txt', 'description': 'HAWK Brackett Gamma', 'zeropoint': 638.9, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'HAWK_NB1060': {'file': 'hawk-nb1060.txt', 'description': 'HAWK Narrow Band 1060', 'zeropoint': 2003.27, 'method': 'vega', 'rsr': False, 'altnames': ['hawk 1060']}, \
    'HAWK_NB1190': {'file': 'hawk-nb1190.txt', 'description': 'HAWK Narrow Band 1190', 'zeropoint': 1697.50, 'method': 'vega', 'rsr': False, 'altnames': ['hawk 1190']}, \
    'HAWK_NB2090': {'file': 'hawk-nb2090.txt', 'description': 'HAWK Narrow Band 2090', 'zeropoint': 706.68, 'method': 'vega', 'rsr': False, 'altnames': ['hawk 2090']}, \
    'IRAC_CH1': {'file': 'irac1.txt', 'description': 'IRAC Channel 1 (3.6 micron)', 'zeropoint': 280.9, 'method': 'vega', 'rsr': True, 'altnames': ['irac 1','irac 3.6']}, \
    'IRAC_CH2': {'file': 'irac2.txt', 'description': 'IRAC Channel 2 (4.5 micron)', 'zeropoint': 179.7, 'method': 'vega', 'rsr': True, 'altnames': ['irac 2','irac 4.5']}, \
    'IRAC_CH3': {'file': 'irac3.txt', 'description': 'IRAC Channel 3 (5.8 micron)', 'zeropoint': 115.0, 'method': 'vega', 'rsr': True, 'altnames': ['irac 3','irac 5.8']}, \
    'IRAC_CH4': {'file': 'irac4.txt', 'description': 'IRAC Channel 4 (8.0 micron)', 'zeropoint': 64.13, 'method': 'vega', 'rsr': True, 'altnames': ['irac 4','irac 8.0']}, \
    'KEPLER': {'file': 'Kepler.txt', 'description': 'Kepler bandpass', 'zeropoint': 3033.1, 'method': 'vega', 'rsr': False, 'altnames': ['kep','kepler k','kp']}, \
    'MKO_J_ATM': {'file': 'j_atm_mko.txt', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'MKO_H_ATM': {'file': 'h_atm_mko.txt', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'MKO_K_ATM': {'file': 'k_atm_mko.txt', 'description': 'MKO K-band + atmosphere', 'zeropoint': 647.7, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'MKO_J': {'file': 'mko_j.txt', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'MKO_H': {'file': 'mko_h.txt', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'MKO_K': {'file': 'mko_ks.txt', 'description': 'MKO K-band', 'zeropoint': 647.7, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'MKO_KP': {'file': 'mko_kp.txt', 'description': 'MKO Kp-band', 'zeropoint': 693.7, 'method': 'vega', 'rsr': False, 'altnames': ['mko k prime']}, \
    'MKO_LP': {'file': 'mko_lp.txt', 'description': 'MKO Lp-band', 'zeropoint': 248.3, 'method': 'vega', 'rsr': False, 'altnames': ['mko l','mko l prime']}, \
    'MKO_MP': {'file': 'mko_mp.txt', 'description': 'MKO Mp-band', 'zeropoint': 164.7, 'method': 'vega', 'rsr': False, 'altnames': ['mko m','mko m prime']}, \
    'NICMOS_F090M': {'file': 'nic1_f090m.txt', 'description': 'NICMOS F090M', 'zeropoint': 2255.0, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F095N': {'file': 'nic1_f095n.txt', 'description': 'NICMOS F095N', 'zeropoint': 2044.6, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F097N': {'file': 'nic1_f097n.txt', 'description': 'NICMOS F097N', 'zeropoint': 2275.4, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F108N': {'file': 'nic1_f108n.txt', 'description': 'NICMOS F108N', 'zeropoint': 1937.3, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F110M': {'file': 'nic1_f110m.txt', 'description': 'NICMOS F110M', 'zeropoint': 1871.8, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F110W': {'file': 'nic1_f110w.txt', 'description': 'NICMOS F110W', 'zeropoint': 1768.5, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F113N': {'file': 'nic1_f113n.txt', 'description': 'NICMOS F113N', 'zeropoint': 1821.0, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F140W': {'file': 'nic1_f140w.txt', 'description': 'NICMOS F140W', 'zeropoint': 1277.1, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F145M': {'file': 'nic1_f145m.txt', 'description': 'NICMOS F145M', 'zeropoint': 1242.0, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F160W': {'file': 'nic1_f160w.txt', 'description': 'NICMOS F160W', 'zeropoint': 1071.7, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F164N': {'file': 'nic1_f164n.txt', 'description': 'NICMOS F164N', 'zeropoint': 1003.0, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F165M': {'file': 'nic1_f165m.txt', 'description': 'NICMOS F165M', 'zeropoint': 1023.6, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F166N': {'file': 'nic1_f166n.txt', 'description': 'NICMOS F166N', 'zeropoint': 1047.7, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F170M': {'file': 'nic1_f170m.txt', 'description': 'NICMOS F170M', 'zeropoint': 979.1, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F187N': {'file': 'nic1_f187n.txt', 'description': 'NICMOS F187N', 'zeropoint': 803.7, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NICMOS_F190N': {'file': 'nic1_f190n.txt', 'description': 'NICMOS F190N', 'zeropoint': 836.5, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NIRC2_J': {'file': 'nirc2-j.txt', 'description': 'NIRC2 J-band', 'zeropoint': 1562.7, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NIRC2_H': {'file': 'nirc2-h.txt', 'description': 'NIRC2 H-band', 'zeropoint': 1075.5, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NIRC2_HCONT': {'file': 'nirc2-hcont.txt', 'description': 'NIRC2 H-continuum band', 'zeropoint': 1044.5, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NIRC2_K': {'file': 'nirc2-k.txt', 'description': 'NIRC2 K-band', 'zeropoint': 648.9, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NIRC2_KP': {'file': 'nirc2-kp.txt', 'description': 'NIRC2 Kp-band', 'zeropoint': 689.3, 'method': 'vega', 'rsr': False, 'altnames': ['nirc2 k prime']}, \
    'NIRC2_KS': {'file': 'nirc2-ks.txt', 'description': 'NIRC2 Ks-band', 'zeropoint': 676.2, 'method': 'vega', 'rsr': False, 'altnames': ['nirc2 k short']}, \
    'NIRC2_KCONT': {'file': 'nirc2-kcont.txt', 'description': 'NIRC2 K continuum-band', 'zeropoint': 605.9, 'method': 'vega', 'rsr': False, 'altnames': ['nirc2 k continuum']}, \
    'NIRC2_FE2': {'file': 'nirc2-fe2.txt', 'description': 'NIRC2 Fe II', 'zeropoint': 1019.7, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NIRC2_LP': {'file': 'nirc2-lp.txt', 'description': 'NIRC2 LP', 'zeropoint': 248.0, 'method': 'vega', 'rsr': False, 'altnames': ['nirc2 l prime','nirc2 l']}, \
    'NIRC2_M': {'file': 'nirc2-ms.txt', 'description': 'NIRC2 M', 'zeropoint': 165.8, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'NIRCAM_F070W': {'file': 'jwst-nircam-F070W.txt', 'description': 'JWST NIRCAM F070W (wide 0.70 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F090W': {'file': 'jwst-nircam-F090W.txt', 'description': 'JWST NIRCAM F090W (wide 0.90 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F115W': {'file': 'jwst-nircam-F115W.txt', 'description': 'JWST NIRCAM F115W (wide 1.15 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F140M': {'file': 'jwst-nircam-F140M.txt', 'description': 'JWST NIRCAM F140M (medium 1.40 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F150W': {'file': 'jwst-nircam-F150W.txt', 'description': 'JWST NIRCAM F150W (wide 1.50 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F150W2': {'file': 'jwst-nircam-F150W2.txt', 'description': 'JWST NIRCAM F150W2 (wide 1.50 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F162M': {'file': 'jwst-nircam-F162M.txt', 'description': 'JWST NIRCAM F162M (medium 1.62 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F164N': {'file': 'jwst-nircam-F164N.txt', 'description': 'JWST NIRCAM F164N (narrow 1.64 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F182M': {'file': 'jwst-nircam-F182M.txt', 'description': 'JWST NIRCAM F182M (medium 1.82 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F187N': {'file': 'jwst-nircam-F187N.txt', 'description': 'JWST NIRCAM F187N (narrow 1.87 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F200W': {'file': 'jwst-nircam-F200W.txt', 'description': 'JWST NIRCAM F200W (wide 2.00 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F210W': {'file': 'jwst-nircam-F210M.txt', 'description': 'JWST NIRCAM F210M (medium 2.10 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F212N': {'file': 'jwst-nircam-F212N.txt', 'description': 'JWST NIRCAM F212N (narrow 2.12 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F250M': {'file': 'jwst-nircam-F250M.txt', 'description': 'JWST NIRCAM F250M (medium 2.50 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F277W': {'file': 'jwst-nircam-F277W.txt', 'description': 'JWST NIRCAM F277W (wide 2.77 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F300M': {'file': 'jwst-nircam-F300M.txt', 'description': 'JWST NIRCAM F300M (medium 3.00 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F322W2': {'file': 'jwst-nircam-F322W2.txt', 'description': 'JWST NIRCAM F322W2 (wide 3.22 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F323N': {'file': 'jwst-nircam-F323N.txt', 'description': 'JWST NIRCAM F323N (narrow 3.23 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F335M': {'file': 'jwst-nircam-F335M.txt', 'description': 'JWST NIRCAM F335M (medium 3.35 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F356W': {'file': 'jwst-nircam-F356W.txt', 'description': 'JWST NIRCAM F356W (wide 3.56 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F360M': {'file': 'jwst-nircam-F360M.txt', 'description': 'JWST NIRCAM F360M (medium 3.60 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F405N': {'file': 'jwst-nircam-F405N.txt', 'description': 'JWST NIRCAM F405N (narrow 4.05 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F410M': {'file': 'jwst-nircam-F410M.txt', 'description': 'JWST NIRCAM F410M (medium 4.10 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F430M': {'file': 'jwst-nircam-F430M.txt', 'description': 'JWST NIRCAM F430M (medium 4.30 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F444W': {'file': 'jwst-nircam-F444W.txt', 'description': 'JWST NIRCAM F444W (wide 4.44 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F460M': {'file': 'jwst-nircam-F460M.txt', 'description': 'JWST NIRCAM F460M (medium 4.60 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F466N': {'file': 'jwst-nircam-F466N.txt', 'description': 'JWST NIRCAM F466N (narrow 4.66 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F470N': {'file': 'jwst-nircam-F470N.txt', 'description': 'JWST NIRCAM F470N (narrow 4.70 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'NIRCAM_F480M': {'file': 'jwst-nircam-F480M.txt', 'description': 'JWST NIRCAM F480M (medium 4.80 micron)', 'zeropoint': 0., 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'PANSTARRS_G': {'file': 'panstarrs-g.txt', 'description': 'PANSTARRS g-band', 'zeropoint': 3909.11, 'method': 'ab', 'rsr': False, 'altnames': []}, \
    'PANSTARRS_R': {'file': 'panstarrs-r.txt', 'description': 'PANSTARRS r-band', 'zeropoint': 3151.44, 'method': 'ab', 'rsr': False, 'altnames': []}, \
    'PANSTARRS_W': {'file': 'panstarrs-w.txt', 'description': 'PANSTARRS w-band', 'zeropoint': 3024.76, 'method': 'ab', 'rsr': False, 'altnames': []}, \
    'PANSTARRS_I': {'file': 'panstarrs-i.txt', 'description': 'PANSTARRS i-band', 'zeropoint': 2584.6, 'method': 'ab', 'rsr': False, 'altnames': []}, \
    'PANSTARRS_Z': {'file': 'panstarrs-z.txt', 'description': 'PANSTARRS z-band', 'zeropoint': 2273.09, 'method': 'ab', 'rsr': False, 'altnames': []}, \
    'PANSTARRS_Y': {'file': 'panstarrs-y.txt', 'description': 'PANSTARRS y-band', 'zeropoint': 2205.95, 'method': 'ab', 'rsr': False, 'altnames': []}, \
    'SDSS_U': {'file': 'sdss-u.txt', 'description': 'SDSS u-band', 'zeropoint': 1568.5, 'method': 'ab', 'rsr': False, 'altnames': ['u']}, \
    'SDSS_G': {'file': 'sdss-g.txt', 'description': 'SDSS g-band', 'zeropoint': 3965.9, 'method': 'ab', 'rsr': False, 'altnames': ['g']}, \
    'SDSS_R': {'file': 'sdss-r.txt', 'description': 'SDSS r-band', 'zeropoint': 3162.0, 'method': 'ab', 'rsr': False, 'altnames': ['r']}, \
    'SDSS_I': {'file': 'sdss-i.txt', 'description': 'SDSS i-band', 'zeropoint': 2602.0, 'method': 'ab', 'rsr': False, 'altnames': ['i']}, \
    'SDSS_Z': {'file': 'sdss-z.txt', 'description': 'SDSS z-band', 'zeropoint': 2244.7, 'method': 'ab', 'rsr': False, 'altnames': ['z']}, \
    'TESS': {'file': 'TESS.txt', 'description': 'TESS bandpass', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'UKIDSS_Z': {'file': 'ukidss-z.txt', 'description': 'UKIDSS Z-band', 'zeropoint': 2261.4, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'UKIDSS_Y': {'file': 'ukidss-y.txt', 'description': 'UKIDSS Y-band', 'zeropoint': 2057.2, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'UKIDSS_J': {'file': 'ukidss-j.txt', 'description': 'UKIDSS J-band', 'zeropoint': 1556.8, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'UKIDSS_H': {'file': 'ukidss-h.txt', 'description': 'UKIDSS H-band', 'zeropoint': 1038.3, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'UKIDSS_K': {'file': 'ukidss-k.txt', 'description': 'UKIDSS K-band', 'zeropoint': 644.1, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'VISTA_Z': {'file': 'vista_z.txt', 'description': 'VISTA Z-band', 'zeropoint': 2263.81, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'VISTA_Y': {'file': 'vista_y.txt', 'description': 'VISTA Y-band', 'zeropoint': 2087.32, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'VISTA_J': {'file': 'vista_j.txt', 'description': 'VISTA J-band', 'zeropoint': 1554.03, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'VISTA_H': {'file': 'vista_h.txt', 'description': 'VISTA H-band', 'zeropoint': 1030.40, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'VISTA_KS': {'file': 'vista_ks.txt', 'description': 'VISTA Ks-band', 'zeropoint': 674.83, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F098M': {'file': 'HST-WFC3_IR_F098M.txt', 'description': 'WFC3 F098M', 'zeropoint': 2154.5, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F105W': {'file': 'HST-WFC3_IR_F105W.txt', 'description': 'WFC3 F105W', 'zeropoint': 1975.2, 'method': 'vega', 'rsr': False, 'altnames': ['wfc3 y']}, \
    'WFC3_F110W': {'file': 'HST-WFC3_IR_F110W.txt', 'description': 'WFC3 F110W', 'zeropoint': 1738.4, 'method': 'vega', 'rsr': False, 'altnames': ['wfc3 yj']}, \
    'WFC3_F125W': {'file': 'HST-WFC3_IR_F125W.txt', 'description': 'WFC3 F125W', 'zeropoint': 1564.3, 'method': 'vega', 'rsr': False, 'altnames': ['wfc3 j']}, \
    'WFC3_F126N': {'file': 'HST-WFC3_IR_F126N.txt', 'description': 'WFC3 F126N', 'zeropoint': 1552.5, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F127M': {'file': 'HST-WFC3_IR_F127M.txt', 'description': 'WFC3 F127M', 'zeropoint': 1496.5, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F128N': {'file': 'HST-WFC3_IR_F128N.txt', 'description': 'WFC3 F128N', 'zeropoint': 1392.6, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F130N': {'file': 'HST-WFC3_IR_F130N.txt', 'description': 'WFC3 F130N', 'zeropoint': 1475.9, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F132N': {'file': 'HST-WFC3_IR_F132N.txt', 'description': 'WFC3 F132N', 'zeropoint': 1466.6, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F139M': {'file': 'HST-WFC3_IR_F139M.txt', 'description': 'WFC3 F139M', 'zeropoint': 1342.8, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F140W': {'file': 'HST-WFC3_IR_F140W.txt', 'description': 'WFC3 F140W', 'zeropoint': 1324.8, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F153M': {'file': 'HST-WFC3_IR_F153M.txt', 'description': 'WFC3 F153M', 'zeropoint': 1142.0, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F160W': {'file': 'HST-WFC3_IR_F160W.txt', 'description': 'WFC3 F160W', 'zeropoint': 1138.1, 'method': 'vega', 'rsr': False, 'altnames': ['wfc3 h']}, \
    'WFC3_F164N': {'file': 'HST-WFC3_IR_F164N.txt', 'description': 'WFC3 F164N', 'zeropoint': 1005.5, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFC3_F167N': {'file': 'HST-WFC3_IR_F167N.txt', 'description': 'WFC3 F167N', 'zeropoint': 1030.0, 'method': 'vega', 'rsr': False, 'altnames': []}, \
#    'WFC3_F127M': {'file': 'wfc3_F127M.txt', 'description': 'WFC3 F127M', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altnames': []}, \
#    'WFC3_F139M': {'file': 'wfc3_F139M.txt', 'description': 'WFC3 F139M', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altnames': []}, \
#    'WFC3_F164N': {'file': 'wfc3_F164N.txt', 'description': 'WFC3 F164N', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altnames': []}, \
#    'WFC3_F167N': {'file': 'wfc3_F167N.txt', 'description': 'WFC3 F167N', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFCAM_Z': {'file': 'wfcam-z.txt', 'description': 'UKIRT WFCAM Z', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFCAM_Y': {'file': 'wfcam-y.txt', 'description': 'UKIRT WFCAM Y', 'zeropoint': 2040.9, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFCAM_J': {'file': 'wfcam-j.txt', 'description': 'UKIRT WFCAM J', 'zeropoint': 1548.7, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFCAM_H': {'file': 'wfcam-h.txt', 'description': 'UKIRT WFCAM H', 'zeropoint': 1027.1, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFCAM_H2': {'file': 'wfcam-h2.txt', 'description': 'UKIRT WFCAM H2', 'zeropoint': 677.1, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WFCAM_BRG': {'file': 'wfcam-brg.txt', 'description': 'UKIRT WFCAM Brackett Gamma', 'zeropoint': 645.5, 'method': 'vega', 'rsr': False, 'altnames': ['wfcam brackett gamma']}, \
    'WFCAM_K': {'file': 'wfcam-k.txt', 'description': 'UKIRT WFCAM K', 'zeropoint': 630.0, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRC_J': {'file': 'wirc_jcont.txt', 'description': 'WIRC J-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRC_H': {'file': 'wirc_hcont.txt', 'description': 'WIRC H-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRC_K': {'file': 'wirc_kcont.txt', 'description': 'WIRC K-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRC_CO': {'file': 'wirc_co.txt', 'description': 'WIRC CO', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRC_CH4S': {'file': 'wirc_ch4s.txt', 'description': 'WIRC CH4S', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRC_CH4L': {'file': 'wirc_ch4l.txt', 'description': 'WIRC CH4L', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRC_FE2': {'file': 'wirc_feii.txt', 'description': 'WIRC Fe II', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRC_BRGAMMA': {'file': 'wirc_brgamma.txt', 'description': 'WIRC H I Brackett Gamma', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': ['wirc brackett gamma']}, \
    'WIRC_PABETA': {'file': 'wirc_pabeta.txt', 'description': 'WIRC H I Paschen Beta', 'zeropoint': 0., 'method': 'vega', 'rsr': False, 'altnames': ['wirc paschen beta']}, \
    'WIRCAM_Y': {'file': 'wircam-cfht-y.txt', 'description': 'CFHT WIRCAM Y', 'zeropoint': 2073.32, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRCAM_J': {'file': 'wircam-cfht-j.txt', 'description': 'CFHT WIRCAM J', 'zeropoint': 1551.01, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRCAM_H': {'file': 'wircam-cfht-h.txt', 'description': 'CFHT WIRCAM H', 'zeropoint': 1044.35, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRCAM_KS': {'file': 'wircam-cfht-ks.txt', 'description': 'CFHT WIRCAM Ks', 'zeropoint': 674.62, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRCAM_KCONT': {'file': 'wircam-cfht-kcont.txt', 'description': 'CFHT WIRCAM K-cont', 'zeropoint': 636.17, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRCAM_CH4_OFF': {'file': 'wircam-cfht-ch4s.txt', 'description': 'CFHT WIRCAM CH4-off', 'zeropoint': 987.39, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WIRCAM_CH4_ON': {'file': 'wircam-cfht-ch4l.txt', 'description': 'CFHT WIRCAM CH4-on', 'zeropoint': 1076.31, 'method': 'vega', 'rsr': False, 'altnames': []}, \
    'WISE_W1': {'file': 'wise_w1.txt', 'description': 'WISE W1 (3.5 micron)', 'zeropoint': 309.54, 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'WISE_W2': {'file': 'wise_w2.txt', 'description': 'WISE W2 (4.6 micron)', 'zeropoint': 171.79, 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'WISE_W3': {'file': 'wise_w3.txt', 'description': 'WISE W3 (13 micron)', 'zeropoint': 31.67, 'method': 'vega', 'rsr': True, 'altnames': []}, \
    'WISE_W4': {'file': 'wise_w4.txt', 'description': 'WISE W4 (22 micron)', 'zeropoint': 8.363, 'method': 'vega', 'rsr': True, 'altnames': []} \
    }
VEGAFILE = 'vega_kurucz.txt'

# some data formats (for future expansion)
INSTRUMENTS = {
#	'SPEX': {'instrument_name': 'SpeX prism', 'pixelscale': 0.15*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['']},
    'UNKNOWN': {'instrument_name': 'UNKNOWN', 'pixelscale': 1.*u.arcsec, 'slitwidth': 1.*u.arcsec, 'altnames': ['UNK'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'RAW': {'instrument_name': 'RAW', 'pixelscale': 1.*u.arcsec, 'slitwidth': 1.*u.arcsec, 'altnames': ['ORIGINAL'], 'wunit': u.micron},
    'SED': {'instrument_name': 'SED', 'pixelscale': 1.*u.arcsec, 'slitwidth': 1.*u.arcsec, 'altnames': ['SPECTRAL_ENERGY_DISTRIBUTION'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'APOGEE': {'instrument_name': 'SDSS APOGEE', 'pixelscale': 2./3.*u.arcsec, 'wave_range': [1.51,1.70]*u.micron, 'slitwidth': 2.*u.arcsec, 'resolution': 22500, 'norders': 1, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altnames': ['APO'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'BOSS': {'instrument_name': 'SDSS BOSS', 'pixelscale': 2./3.*u.arcsec, 'wave_range': [3700,10400]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 2000, 'norders': 1, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altnames': ['SDSS','BOSS','EBOSS'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'DEIMOS': {'instrument_name': 'Keck DEIMOS', 'pixelscale': 0.1185*u.arcsec, 'wave_range': [5000,9700]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'disperser': '600 l/mm', 'resolution': 2000, 'norders': 1, 'readnoise': 2.5, 'darkcurrent': 0., 'gain': 1.2, 'altnames': ['DEIMOS'], 'instrument_reference': '2003spie.4841.1657F', 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'FIRE': {'instrument_name': 'Magellan FIRE', 'pixelscale': 0.18*u.arcsec, 'wave_range': [0.82,2.51]*u.micron, 'slitwidth': 0.6*u.arcsec, 'resolution': 6000, 'norders': 21, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altnames': ['FIRE'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'IRS-SL': {'instrument_name': 'Spitzer IRS Short-Low', 'pixelscale': 1.8*u.arcsec, 'wave_range': [5.2,14.5]*u.micron, 'slitwidth': 1.8*u.arcsec, 'resolution': 100, 'norders': 1, 'readnoise': 30., 'darkcurrent': 10., 'gain': 4.6, 'altnames': ['IRS','IRS Short Low'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'IRS-LL': {'instrument_name': 'Spitzer IRS Short-Low', 'pixelscale': 5.1*u.arcsec, 'wave_range': [14,38.]*u.micron, 'slitwidth': 5.1*u.arcsec, 'resolution': 90, 'norders': 1, 'readnoise': 30., 'darkcurrent': 40., 'gain': 4.6, 'altnames': ['IRS Long Low'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'KAST-RED': {'instrument_name': 'Lick KAST Red Channel', 'pixelscale': 0.43*u.arcsec, 'disperser': '600/7500', 'wave_range': [5700,9200]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 1200, 'norders': 1, 'readnoise': 3.8, 'darkcurrent': 0., 'gain': 1.9, 'altnames': ['KAST','KAST-R'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'KAST-BLUE': {'instrument_name': 'Lick KAST Blue Channel', 'pixelscale': 0.43*u.arcsec, 'disperser': '600/4310', 'wave_range': [3300,5520]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 950, 'norders': 1, 'readnoise': 3.7, 'darkcurrent': 0., 'gain': 1.2, 'altnames': ['KAST-B'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'LDSS-3': {'instrument_name': 'Magellan LDSS-3', 'pixelscale': 0.189*u.arcsec, 'disperser': 'VPH-RED', 'wave_range': [6000,10000]*u.Angstrom, 'slitwidth': 0.75*u.arcsec, 'resolution': 1810, 'norders': 1, 'readnoise': 4.07, 'darkcurrent': 0., 'gain': 1, 'altnames': ['LDSS3'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'LRIS-RED': {'altnames': ['LRIS','LRISR'], 'instrument_name': 'Keck LRIS red channel longslit', 'instrument_reference': '1995PASP..107..375O', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4.6, 'darkcurrent': 0., 'gain': 1.2, 'disperser': '400/8500', 'wave_range': [6300,10100]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1200, 'norders': 1, 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'LRIS-BLUE': {'altnames': ['LRISB'], 'instrument_name': 'Keck LRIS blue channel longslit', 'instrument_reference': '1998SPIE.3355...81M', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4., 'darkcurrent': 0., 'gain': 1.6, 'disperser': '400/3400', 'wave_range': [1270,5740]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 530, 'norders': 1, 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'MAGE': {'instrument_name': 'Magellan MAGE', 'pixelscale': 0.3*u.arcsec, 'wave_range': [3100,10000]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 4100, 'norders': 13, 'readnoise': 2.9, 'darkcurrent': 1.0, 'gain': 1.02, 'altnames': ['MagE'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom},
    'NIRI-J': {'instrument_name': 'Gemini NIRI J-band (G5202)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [1.05,1.41]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 610, 'norders': 1, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altnames': ['NIRI','NIRI G5202'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'NIRI-H': {'instrument_name': 'Gemini NIRI J-band (G5203)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [1.43,1.96]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 825, 'norders': 1, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altnames': ['NIRI G5203'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'NIRI-K': {'instrument_name': 'Gemini NIRI J-band (G5204)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [1.90,2.49]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 780, 'norders': 1, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altnames': ['NIRI G5204'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'NIRI-L': {'instrument_name': 'Gemini NIRI J-band (G5205)', 'pixelscale': 0.47/4.*u.arcsec, 'wave_range': [2.99,4.15]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 690, 'norders': 1, 'readnoise': 50., 'darkcurrent': 0.25, 'gain': 12.3, 'altnames': ['NIRI G5205'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'NIRSPEC': {'instrument_name': 'Keck NIRSPEC', 'pixelscale': 0.43/3.*u.arcsec, 'wave_range': [0.95,5.5]*u.micron, 'slitwidth': 0.43*u.arcsec, 'resolution': 25000, 'norders': 8, 'readnoise': 23., 'darkcurrent': 0.8, 'gain': 5.8, 'altnames': ['NIRSPAO'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'SPEX_PRISM': {'instrument_name': 'IRTF SpeX prism', 'pixelscale': 0.15*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['SPEX','PRISM'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'SPEX_SXD': {'instrument_name': 'IRTF SpeX SXD', 'pixelscale': 0.15*u.arcsec, 'wave_range': [0.8,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['SXD'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'SPEX_LXD1.9': {'instrument_name': 'IRTF SpeX LXD 1.9 micron', 'pixelscale': 0.15*u.arcsec, 'wave_range': [1.95,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['LXD'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'SPEX_LXD2.1': {'instrument_name': 'IRTF SpeX LXD 2.1 micron', 'pixelscale': 0.15*u.arcsec, 'wave_range': [2.15,5.0]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': [], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'SPEX_LXD2.3': {'instrument_name': 'IRTF SpeX LXD 2.3 micron', 'pixelscale': 0.15*u.arcsec, 'wave_range': [2.25,5.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': [], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
#	'USPEX': {'instrument_name': 'Updated SpeX prism', 'pixelscale': 0.10*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['']},
#	'USPEX_PRISM': {'instrument_name': 'IRTF Updated SpeX prism', 'pixelscale': 0.10*u.arcsec, 'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['USPEX','UPRISM'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
#    'USPEX_SXD': {'instrument_name': 'IRTF Updated SpeX SXD', 'pixelscale': 0.10*u.arcsec, 'wave_range': [0.7,2.55]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['USXD'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
#    'USPEX_LXD_SHORT': {'instrument_name': 'IRTF Updated SpeX LXD short', 'pixelscale': 0.10*u.arcsec, 'wave_range': [1.67,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500., 'norders': 8, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['ULXD','LXD_SHORT','LXDS'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
#    'USPEX_LXD_LONG': {'instrument_name': 'IRTF Updated SpeX LXD long', 'pixelscale': 0.10*u.arcsec, 'wave_range': [1.98,5.3]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500., 'norders': 7, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['LXD_LONG','LXDL'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'WFC3_G102': {'instrument_name': 'HST WFC3 IR G102', 'pixelscale': 0.128*u.arcsec, 'wave_range': [0.8,1.15]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 210., 'norders': 7, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altnames': ['G102'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
    'WFC3_G141': {'instrument_name': 'HST WFC3 IR G141', 'pixelscale': 0.128*u.arcsec, 'wave_range': [1.075,1.70]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 130., 'norders': 7, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altnames': ['WFC3','HST WFC3','HST WFC3 IR','WFC3 IR','G141'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron},
}

#experimental
INSTRUMENTS_ALT = { 
    'SED': {'instrument_name': 'SED', 'pixelscale': 0.*u.arcsec, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altnames': ['SPECTRAL_ENERGY_DISTRIBUTION'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'DEFAULT': {'wave_range': [0.1,100.]*u.micron, 'slitwidth': 0.*u.arcsec, 'resolution': 100, 'norders': 1, 'default': True}, 
    }}, 
    'APOGEE': {'instrument_name': 'SDSS APOGEE', 'pixelscale': 2./3.*u.arcsec, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altnames': ['APO'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'DEFAULT': {'wave_range': [1.51,1.70]*u.micron, 'slitwidth': 2.*u.arcsec, 'resolution': 22500, 'norders': 1, 'default': True}, 
    }}, 
    'BOSS': {'instrument_name': 'SDSS BOSS', 'pixelscale': 2./3.*u.arcsec, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altnames': ['BOSS','EBOSS'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': { 
        'DEFAULT': {'wave_range': [3700,10400]*u.Angstrom, 'slitwidth': 2.*u.arcsec, 'resolution': 2000, 'norders': 1, 'default': True}, 
    }}, 
    'FIRE': {'instrument_name': 'Magellan FIRE', 'pixelscale': 0.18*u.arcsec, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altnames': ['FIRE'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': {
        'ECHELLE': {'wave_range': [0.82,2.51]*u.micron, 'slitwidth': 0.6*u.arcsec, 'resolution': 6000, 'norders': 21, 'default': True},
        'PRISM': {'wave_range': [0.82,2.51]*u.micron, 'slitwidth': 0.6*u.arcsec, 'resolution': 450, 'norders': 1, 'default': False},
    }},
    'IRS': {'instrument_name': 'Spitzer IRS', 'readnoise': 30., 'gain': 4.6, 'altnames': ['IRS'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': {
        'SL': {'wave_range': [5.2,14.5]*u.micron, 'slitwidth': 1.8*u.arcsec, 'pixelscale': 1.8*u.arcsec, 'darkcurrent': 10., 'resolution': 100, 'norders': 1, 'default': True},
        'SH': {'wave_range': [9.9,19.6]*u.micron, 'slitwidth': 2.3*u.arcsec, 'pixelscale': 2.3*u.arcsec, 'darkcurrent': 10., 'resolution': 600, 'norders': 1, 'default': False},
        'LL': {'wave_range': [14,38.]*u.micron, 'slitwidth': 5.1*u.arcsec, 'pixelscale': 5.1*u.arcsec, 'darkcurrent': 40., 'resolution': 90, 'norders': 1, 'default': False},
        'LH': {'wave_range': [18.7,37.2]*u.micron, 'slitwidth': 4.5*u.arcsec, 'pixelscale': 4.5*u.arcsec, 'darkcurrent': 40., 'resolution': 90, 'norders': 1, 'default': False},
    }}, 
    'LDSS-3': {'instrument_name': 'Magellan LDSS-3', 'pixelscale': 0.189*u.arcsec, 'readnoise': 4.07, 'darkcurrent': 0., 'gain': 1, 'altnames': ['LDSS3'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': {
        'VPH-RED': {'wave_range': [6000,10000]*u.Angstrom, 'slitwidth': 0.75*u.arcsec, 'resolution': 1810, 'norders': 1, 'default': True}, 
        'VPH-BLUE': {'wave_range': [3800,6200]*u.Angstrom, 'slitwidth': 0.75*u.arcsec, 'resolution': 1900, 'norders': 1, 'default': False}, 
        'VPH-ALL': {'wave_range': [4250,10000]*u.Angstrom, 'slitwidth': 0.75*u.arcsec, 'resolution': 860, 'norders': 1, 'default': False}, 
    }}, 
    'LRIS-BLUE': {'instrument_name': 'Keck LRIS blue channel longslit', 'altnames': ['LRISB'], 'instrument_reference': '1998SPIE.3355...81M', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4., 'darkcurrent': 0., 'gain': 1.6, 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': { 
        '300/5000': {'default': False, 'wave_range': [1570,7420]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 570, 'norders': 1}, 
        '400/3400': {'default': True, 'wave_range': [1270,5740]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 500, 'norders': 1}, 
        '600/4000': {'default': False, 'wave_range': [3010,5600]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1000, 'norders': 1}, 
        '1200/3400': {'default': False, 'wave_range': [2910,3890]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 2200, 'norders': 1}, 
    }}, 
    'LRIS-RED': {'instrument_name': 'Keck LRIS red channel longslit', 'altnames': ['LRIS','LRISR'], 'instrument_reference': '1995PASP..107..375O', 'pixelscale': 0.135*u.arcsec, 'readnoise': 4.6, 'darkcurrent': 0., 'gain': 1.2, 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': { 
        '300/5000': {'default': False, 'wave_range': [3700,8270]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 540, 'norders': 1}, 
        '400/8500': {'default': True, 'wave_range': [6200,10900]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1200, 'norders': 1}, 
        '600/5000': {'default': False, 'wave_range': [3700,6640]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1050, 'norders': 1}, 
        '600/7500': {'default': False, 'wave_range': [6200,9140]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 1600, 'norders': 1}, 
        '600/10000': {'default': False, 'wave_range': [8700,10900]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 2130, 'norders': 1}, 
    }}, 
    'MAGE': {'instrument_name': 'Magellan MAGE', 'pixelscale': 0.3*u.arcsec, 'readnoise': 2.9, 'darkcurrent': 1.0, 'gain': 1.02, 'altnames': ['MagE'], 'wunit': u.Angstrom, 'funit': u.erg/u.s/u.cm/u.cm/u.Angstrom, 'modes': { 
        'DEFAULT': {'wave_range': [3100,10000]*u.Angstrom, 'slitwidth': 1.*u.arcsec, 'resolution': 4100, 'norders': 13, 'default': True}, 
    }}, 
    'NIRI': {'instrument_name': 'Gemini NIRI', 'pixelscale': 0.47/4.*u.arcsec, 'readnoise': 13., 'darkcurrent': 0.25, 'gain': 12.3, 'altnames': ['GEMINI NIRI'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': {
        'J': {'default': True, 'wave_range': [1.05,1.41]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 610, 'norders': 1},
        'H': {'default': False, 'wave_range': [1.43,1.96]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 825, 'norders': 1},
        'K': {'default': False, 'wave_range': [1.90,2.49]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 780, 'norders': 1},
        'L': {'default': False, 'wave_range': [2.99,4.15]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 690, 'norders': 1},
        'M': {'default': False, 'wave_range': [4.45,5.45]*u.micron, 'slitwidth': 0.47*u.arcsec, 'resolution': 770, 'norders': 1},
    }}, 
    'NIRSPEC': {'instrument_name': 'Keck NIRSPEC', 'pixelscale': 0.43/3.*u.arcsec, 'readnoise': 23., 'darkcurrent': 0.8, 'gain': 5.8, 'altnames': ['NIRSPAO'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
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
    'SPEX': {'instrument_name': 'IRTF SpeX', 'pixelscale': 0.15*u.arcsec, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['OLD SPEX'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'PRISM': {'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'default': True}, 
        'SXD': {'wave_range': [0.8,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'default': False}, 
        'LXD1.9': {'wave_range': [1.95,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'default': False}, 
        'LXD2.1': {'wave_range': [2.15,5.0]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'default': False}, 
        'LXD2.3': {'wave_range': [2.25,5.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'default': False}, 
    }}, 
    'USPEX': {'instrument_name': 'IRTF SpeX (upgraded)', 'pixelscale': 0.10*u.arcsec, 'readnoise': 5., 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['SPEX'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'PRISM': {'wave_range': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'default': True}, 
        'SXD': {'wave_range': [0.8,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'default': False}, 
        'LXD_SHORT': {'wave_range': [1.67,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 8, 'default': False}, 
        'LXD_LONG': {'wave_range': [1.98,5.3]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'default': False}, 
    }}, 
    'WFC3-UVIS': {'instrument_name': 'HST WFC3 UVIS', 'pixelscale': 0.0395*u.arcsec, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altnames': ['WFC3-OPT','WFC3-VIS','WFC3-UV'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
        'DEFAULT': {'wave_range': [1900,4500]*u.Angstrom, 'slitwidth': 0.3*u.arcsec, 'resolution': 70, 'norders': 1, 'default': True}, 
    }}, 
    'WFC3-IR': {'instrument_name': 'HST WFC3 IR', 'pixelscale': 0.128*u.arcsec, 'readnoise': 0., 'darkcurrent': 0., 'gain': 1., 'altnames': ['WFC3'], 'wunit': u.micron, 'funit': u.erg/u.s/u.cm/u.cm/u.micron, 'modes': { 
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
#    'atmos-cond': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/atmos-cond/', 'name': 'ATMOS-COND', 'citation': 'Phillips et al. (in prep)', 'bibcode': 'TBD', 'altnames': ['atmos'], 'rawfolder': HOME_FOLDER+'/models/atmos/cond/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
 #   'atmos-rain': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/atmos-rain/', 'name': 'ATMOS-RAINOUT', 'citation': 'Phillips et al. (in prep)', 'bibcode': 'TBD', 'altnames': ['atmos-rainout'], 'rawfolder': HOME_FOLDER+'/models/atmos/rain/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'nextgen99': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/nextgen99/', 'name': 'AMES NextGen', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altnames': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'rawfolder': HOME_FOLDER+'/models/phoenix/nextgen/fullres/', 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
#    'gaia': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/gaia/', 'name': 'AMES GAIA', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altnames': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'rawfolder': HOME_FOLDER+'/models/phoenix/nextgen/fullres/', 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
    'cond01': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/cond01/', 'name': 'AMES Cond', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altnames': ['cond','cond-ames','amescond'], 'rawfolder': HOME_FOLDER+'/models/allard/cond01/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'dusty01': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/dusty01/', 'name': 'AMES Dusty', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altnames': ['dusty','dusty-ames','amesdusty'], 'rawfolder': HOME_FOLDER+'/models/allard/bddusty01/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'drift': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/drift/', 'name': 'Drift', 'citation': 'Witte et al. (2011)', 'bibcode': '2011A&A...529A..44W', 'altnames': ['witte','witte11','witte2011','helling'], 'rawfolder': HOME_FOLDER+'/models/drift/v1.2/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
    'burrows06': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/burrows06/', 'name': 'Burrows 2006', 'citation': 'Burrows et al. (2006)', 'bibcode': '2006ApJ...640.1063B', 'altnames': ['burrows','burrows2006'], 'rawfolder': HOME_FOLDER+'/models/burrows/burrows06/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'cld': 'nc'}}, \
#    'btcond': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/btcond/', 'name': 'BT Cond', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altnames': ['cond-bt','btcond'], 'rawfolder': HOME_FOLDER+'/models/allard/cond01/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
#    'btdusty': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/btdusty/', 'name': 'BT Dusty', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altnames': ['dusty-bt','btdusty'], 'rawfolder': HOME_FOLDER+'/models/allard/bddusty01/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'btnextgen': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/btnextgen/', 'name': 'BT NextGen', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altnames': ['nextgen-bt','btnextgen'], 'rawfolder': HOME_FOLDER+'/models/allard/btnextgen/forsplat/', 'default': {'teff': 3000., 'logg': 5.0, 'z': 0.0, 'enrich': 0.}}, \
    'btsettl08': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/btsettl08/', 'name': 'BTSettl08', 'citation': 'BT-Settl (2008)', 'bibcode': '2012RSPTA.370.2765A', 'altnames': ['allard','allard12','allard2012','btsettl','btsettled','btsettl08','btsettl2008','BTSettl2008'], 'rawfolder': HOME_FOLDER+'/models/allard/cifist2011/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'enrich': 0.}}, \
    'btsettl15': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/btsettl15/', 'name': 'BTSettl15', 'citation': 'BT-Settl (2015)', 'bibcode': '2015A&A...577A..42B', 'altnames': ['allard15','allard2015','btsettl015','btsettl2015','BTSettl2015'], 'rawfolder': HOME_FOLDER+'/models/allard/cifist2015/BT-Settl_M-0.0a+0.0/', 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.}}, \
    'burrows06': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/burrows06/', 'name': 'Burrows 2006', 'citation': 'Burrows et al. (2006)', 'bibcode': '2006ApJ...640.1063B', 'altnames': ['burrows','burrows2006'], 'rawfolder': HOME_FOLDER+'/models/burrows/burrows06/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'cld': 'nc'}}, \
    'madhusudhan11': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/madhusudhan11/', 'name': 'Madhusudhan 2011', 'citation': 'Madhusudhan et al. (2011)', 'bibcode': '2011ApJ...737...34M', 'altnames': ['madhu','madhusudhan','madhu11','madhu2011','madhusudhan2011'], 'rawfolder': HOME_FOLDER+'/models/burrows/burrows11/all/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.,'cld': 'ae60', 'kzz': 'eq','fsed': 'eq'}}, \
    'morley12': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/morley12/', 'name': 'Morley 2012', 'citation': 'Morley et al. (2012)', 'bibcode': '2012ApJ...756..172M', 'altnames': ['morley','morley2012'], 'rawfolder': HOME_FOLDER+'/models/ames/Morley12/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'fsed': 'f5'}}, \
    'morley14': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/morley14/', 'name': 'Morley 2014', 'citation': 'Morley et al. (2014)', 'bibcode': '2014ApJ...787...78M', 'altnames': ['morley2014'], 'rawfolder': HOME_FOLDER+'/models/ames/Morley14/', 'default': {'teff': 300., 'logg': 5.0, 'z': 0., 'fsed': 'f5', 'cld': 'h50'}}, \
    'saumon12': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/saumon12/', 'name': 'Saumon 2012', 'citation': 'Saumon et al. (2012)', 'bibcode': '2012ApJ...750...74S', 'altnames': ['saumon','saumon2012'], 'rawfolder': HOME_FOLDER+'/models/ames/Saumon12/', 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
#    'tremblin16': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/tremblin16/', 'name': 'Tremblin 2016 (test)', 'citation': 'Tremblin et al. (2016)', 'bibcode': '2016ApJ...817L..19T', 'altnames': ['tremblin'], 'rawfolder': HOME_FOLDER+'/models/tremblin/20161120/', 'default': {'teff': 1300., 'logg': 5.0, 'z': 0.,'kzz': '6.0', 'ad': 1.05}}, \
}
SPECTRAL_MODEL_PARAMETERS_INORDER = ['teff','logg','z','fsed','cld','kzz','ad','enrich']
SPECTRAL_MODEL_PARAMETERS = {\
    'teff': {'name': 'temperature', 'prefix': 't', 'unit': u.K, 'default': 1000.0, 'title': '$T_{eff}$ (K)', 'type': 'continuous'}, \
    'logg': {'name': 'gravity', 'prefix': 'g', 'unit': u.dex, 'default': 5.0, 'title': '$\log{g}$ (cgs)', 'type': 'continuous'}, \
    'z': {'name': 'metallicity', 'prefix': 'z', 'unit': u.dex, 'default': 0., 'title': '$[M/H]$', 'type': 'continuous'}, \
    'fsed': {'name': 'rainout', 'prefix': 'f', 'unit': u.m/u.m, 'default': 'nc', 'title': '$f_{sed}$', 'type': 'discrete'}, \
    'cld': {'name': 'cloud', 'prefix': 'c', 'unit': u.m/u.m, 'default': 'nc', 'title': 'Cloud Treatment', 'type': 'discrete'}, \
    'kzz': {'name': 'mixing', 'prefix': 'k', 'unit': u.m/u.m, 'default': 'eq', 'title': '$log\ \kappa_{zz}$ (cgs)', 'type': 'discrete'},\
    'ad': {'name': 'adiabat', 'prefix': 'a', 'unit': u.m/u.m, 'default': 1., 'title': 'Adiabatic Index', 'type': 'continuous'},\
    'enrich': {'name': 'enrichment', 'prefix': 'e', 'unit': u.dex, 'default': 0., 'title': 'Alpha Element Enrichment', 'type': 'continuous'},\
    'radius': {'name': 'radius', 'prefix': 'r', 'unit': u.Rsun, 'default': 0., 'title': 'Radius (R$_{\odot}$)', 'type': 'continuous'}}


# evolutionary model information
EVOLUTIONARY_MODELS = {\
    'baraffe03': {'name': 'Barraffe 2003','citation': 'Baraffe et al. (2003)', 'bibcode': '2003A&A...402..701B', 'altnames': ['bar03','baraffe2003']},\
    'baraffe15': {'name': 'Barraffe 2015','citation': 'Baraffe et al. (2015)', 'bibcode': '2015A&A...577A..42B', 'altnames': ['baraffe','bar15','baraffe2015']},\
    'burrows01': {'name': 'Burrows 2001','citation': 'Burrows et al. (2001)', 'bibcode': '2015A&A...577A..42B', 'altnames': ['burrows','bur01','burrows2001']},\
    'saumon08': {'name': 'Saumon 2008','citation': 'Saumon et al. (2008)', 'bibcode': '2015A&A...577A..42B', 'altnames': ['saumon','sau08','saumon2008']}}

EVOLUTIONARY_MODEL_PARAMETERS = {\
    'mass': {'unit': u.solMass, 'default': 0.05, 'title': '$M$'},\
    'age': {'unit': u.Gyr, 'default': 5., 'title': '$\tau$'},\
    'temperature': {'unit': u.K, 'default': 1000.0, 'title': '$T_{eff}$'},\
    'gravity': {'unit': u.dex(u.cm / u.s**2), 'default': 5.0, 'title': '$\log{g}'},\
    'luminosity': {'unit': u.dex(u.solLum), 'default': -5., 'title': '$\log{L_{bol}/L_{\odot}}$'},\
    'radius': {'unit': u.solRad, 'default': 0.1, 'title': '$R_{\odot}$'}}

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

ABSMAG_SETS = {
    'dahn2002': {'altname': ['dahn','dahn02'], 'bibcode': '2002AJ....124.1170D', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.25, 'range' : [17.,28.], 'coeff': [0.341,8.38]}}},
    'cruz2003': {'altname': ['cruz','cruz03'], 'bibcode': '2003AJ....126.2421C', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.30, 'range' : [16.,28.], 'coeff': [-6.892e-4,3.453e-2,-6.193e-1,5.043,-4.410]}}},
    'burgasser2007': {'altname': ['burgasser','burgasser07'], 'bibcode': '2007ApJ...659..655B', 'sptoffset': 20, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.30, 'range' : [20., 38.], 'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, 
        'MKO_H': {'fitunc' : 0.27, 'range' : [20., 38.], 'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, 
        'MKO_K': {'fitunc' : 0.26, 'range' : [20., 38.], 'coeff': [.0000001051, -.000006985, .0001807, -.002271, .01414, -.04024, .05129, .2322, 10.45]}}},
    'looper2008': {'altname': ['looper','looper08'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 20, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc' : 0.29, 'range' : [20., 38.], 'coeff': [-5.462e-6,2.595e-4,-3.915e-3,1.663e-2,3.690e-2,1.255e-1,11.817]},
        '2MASS_H': {'fitunc' : 0.29, 'range' : [20., 38.], 'coeff' : [-4.218e-6,1.987e-4,-2.970e-3,1.261e-2,3.032e-2,1.125e-1,11.010]},
        '2MASS_KS': {'fitunc' : 0.33, 'range' : [20., 38.], 'coeff' : [-4.104e-6,1.911e-4,-2.864e-3,1.299e-2,2.565e-2,7.369e-2,10.521]}}},
    'faherty2012': {'altname': ['faherty12'],'bibcode': '2012ApJ...752...56F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.30, 'range' : [20., 38.], 'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, 
        'MKO_H': {'fitunc' : 0.27, 'range' : [20., 38.], 'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, 
        'MKO_K': {'fitunc' : 0.28, 'range' : [20., 38.], 'coeff' : [.0000816516, -.00469032, .0940816, -.485519, 9.76100]}}},
    'dupuy2012': {'altname': ['dupuy','dupuy12'], 'bibcode': '2012ApJS..201...19D', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'MKO_J': {'fitunc' : 0.39, 'range' : [16., 39.], 'coeff' : [-.00000194920, .000227641, -.0103332, .232771, -2.74405, 16.3986, -28.3129]}, 
        'MKO_Y': {'fitunc': 0.40, 'range' : [16., 39.], 'coeff': [-.00000252638, .000285027, -.0126151, .279438, -3.26895, 19.5444, -35.1560]}, 
        'MKO_H': {'fitunc': 0.38, 'range' : [16., 39.], 'coeff': [-.00000224083, .000251601, -.0110960, .245209, -2.85705, 16.9138, -29.7306]}, 
        'MKO_K': {'fitunc': 0.40, 'range' : [16., 39.], 'coeff': [-.00000104935, .000125731, -.00584342, .135177, -1.63930, 10.1248, -15.2200]}, 
        'MKO_LP': {'fitunc': 0.28, 'range': [16., 39.], 'coeff': [0.00000, 0.00000, .0000546366, -.00293191, .0530581,  -.196584, 8.89928]}, 
        '2MASS_J': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [-.000000784614, .000100820, -.00482973, .111715, -1.33053, 8.16362, -9.67994]}, 
        '2MASS_H': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [-.00000111499, .000129363, -.00580847, .129202, -1.50370, 9.00279, -11.7526]}, 
        '2MASS_KS': {'fitunc': 0.43, 'range':[16., 39.], 'coeff': [1.06693e-4, -6.42118e-3, 1.34163e-1, -8.67471e-1, 1.10114e1]}, 
        'WISE_W1': {'fitunc': 0.39, 'range':[16., 39.], 'coeff': [1.58040e-5, -3.33944e-4, -4.38105e-3, 3.55395e-1, 7.14765]}, 
        'WISE_W2': {'fitunc': 0.35, 'range':[16., 39.], 'coeff': [1.78555e-5, -8.81973e-4, 1.14325e-2, 1.92354e-1, 7.46564]}}},
    'tinney2013': {'altname': ['tinney','tinney03'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        'COUSINS_I': {'fitunc' : 0.37, 'range' : [20., 37.5], 'coeff': [-2.49821e-6,1.04398e-3,-6.49719e-2,1.56038,-1.58296e1,7.22089e1]},
        '2MASS_J': {'fitunc' : 0.38, 'range' : [20., 37.5], 'coeff': [-1.25074e-5,1.63124e-3,-7.42418e-2,1.54509,-1.47407e1,6.27861e1]},
        '2MASS_KS': {'fitunc' : 0.36, 'range' : [20., 37.5], 'coeff': [-2.80824e-6,3.41146e-4,-1.73848e-2,4.82120e-1,-7.86911,7.57222e1,-3.98105e2,8.94012e2]}}},
    'tinney2014': {'altname': ['tinney15'],'bibcode': '2014ApJ...796...39T', 'sptoffset': 0, 'method': 'interpolate', 'filters': {
        'MKO_J': {'spt': [36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,42], 'absmag': [15.22,15.49,16.39,16.66,17.9,18.35,19.08,20.32,22.39,22.18,25.76], 'rms': [0.31,0.37,0.72,0.36,0.46,0.9,0.97,1.25,1.,0.76,3.52]}, 
        'WISE_W2': {'spt': [36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,42], 'absmag': [12.86,13.28,13.39,13.44,13.75,13.92,14.28,14.65,15.2,14.78,15.76], 'rms': [0.17,0.48,0.27,0.23,0.22,0.24,0.46,0.35,1.,0.77,2.15]}}}, 
    'filippazzo2015': {'altname': ['filippazzo','filippazzo15','fillippazzo','filipazo','filippazo'],'bibcode': '2015ApJ...810..158F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
        '2MASS_J': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [3.478e-5, -2.684e-3, 7.771e-2, -1.058, 7.157, -8.350]}, 
        'WISE_W2': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [8.190e-6, -6.938e-4, 2.283e-2, -3.655e-1, 3.032, -5.043e-1]}}},
    'faherty2016': {'altname': ['faherty','faherty2016,faherty-field'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 10, 'method': 'polynomial', 'filters': {
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
}

# telscope locations
TELESCOPES = {
    'KECK': {'lat': 19.8263*u.deg, 'lon': -155.4783*u.deg, 'height': 4160*u.m, 'altname': ['MAUNA_KEA','SUBARU','IRTF','CHFT','UKIRT','GEMINI_N','GEMINI_NORTH']},
    'VLT': {'lat': -24.6275*u.deg, 'lon': -70.4044*u.deg, 'height': 2636*u.m, 'altname': ['PARANAL','CERRO_PARANAL','VERY_LARGE_TELESCOPE','ESA']},
}

# more faherty2016 relations
#Teff FLD 6.0<SpT <29.0 113.431 4.747e+03 -7.005e+02 1.155e+02 -1.191e+01 6.318e-01 -1.606e-02 1.546e-04
#Teff YNG 7.0<SpT <17.0 180.457 1.330e+00 -66.8637 1235.42 -10068.8 32766.4 · · · · · ·
#Teff YNG2 7.0<SpT <28.0 197.737 2.795e+04 -9.183e+03 1.360e+03 -1.066e+02 4.578e+00 -1.016e-01 9.106e-04
#Teff GRP 7.0<SpT <17.0 172.215 7.383e+00 -344.522 4879.86 · · · · · · · · · · · ·
#Lbol FLD 7.0<SpT <28.0 0.133 2.787e+00 -2.310e+00 3.727e-01 -3.207e-02 1.449e-03 -3.220e-05 2.736e-07
#Lbol YNG 7.0<SpT <17.0 0.335 -6.514e-03 2.448e-01 -3.113e+00 9.492e+00 · · · · · · · · ·
#Lbol YNG2 7.0<SpT <28.0 0.206 2.059e-01 9.585 -3.985 4.923e-01 -3.048e-02 9.134e-04 -1.056e-05
#Lbol GRP 7.0<SpT <17.0 0.221 6.194e-03 -3.757e-01 2.728e-02 · · · · · · · · · · · ·

