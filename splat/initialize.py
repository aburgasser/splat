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
SPLAT_EMAIL = 'aburgasser@gmail.com'
DATA_FOLDER = '/reference/Spectra/'
FILTER_FOLDER = '/reference/Filters/'
SPECTRAL_MODEL_FOLDER = '/reference/SpectralModels/'
EVOLUTIONARY_MODEL_FOLDER = '/reference/EvolutionaryModels/'
DOCS_FOLDER = '/docs/'
DOCS_INDEX_HTML = '/docs/_build/html/index.html'
DB_FOLDER = '/db/'
DB_SOURCES_FILE = 'source_data.txt'
DB_SPECTRA_FILE = 'spectral_data.txt'
DB_PHOTOMETRY_FILE = 'photometry_data.txt'
BIBFILE = 'splat_bibs.bib'
TMPFILENAME = 'splattmpfile'
ACCESS_FILE = '.splat_access'
HOME_FOLDER = os.path.expanduser('~')

MONTHS = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

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
# EMPTY DICTIONARY

# subdwarf spectral standards
STDS_SD_SPEX_KEYS = { \
    'sdM5.5': 11670, #'11670_11134.fits',\
    'sdM6.0': 10265, #'10265_10045.fits',\
    'sdM7.0': 10197, #'10197_11074.fits',\
    'sdM8.0': 10123, #'10123_10145.fits',\
    'sdM9.5': 10188, #'10188_10700.fits',\
    'sdL0.0': 11972, #'11972_10248.fits',\
    'sdL3.5': 10364, #'10364_10946.fits',\
    'sdL4.0': 10203} #'10203_11241.fits'}
# EMPTY DICTIONARY

# extreme subdwarf spectral standards
STDS_ESD_SPEX_KEYS = { \
    'esdM5.0': 10229, #'10229_10163.fits',\
#    'esdM6.5': '_10579.fits',\
    'esdM7.0': 10521, #'10521_10458.fits',\
    'esdM8.5': 10278} #'10278_10400.fits'}
# EMPTY DICTIONARY

# young spectral standards from Allers & Liu (2013) and Cruz et al. (2017)
STDS_VLG_SPEX_KEYS = { \
#    'M5.0': 10228,  # TWA11C - don't have this
#    'M6.0': 10228,  # TWA8B - don't have this
    'M7.0': 11199,  # 0335+2342
    'M8.0': 10799,  # TWA27A
    'M9.0': 10797,  # TWA26
    'L0.0': 10228,  # 0141-4633 - proposed by Cruz et al. 2017
    'L1.0': 11178,  # 0518-2756
    'L2.0': 11696,  # 0536-1920
    'L3.0': 11198,  # 2208+2921
    'L4.0': 11177,  # 0501-0010 - proposed by Cruz et al. 2017
    'L6.0': 10455}  # 2244+2043

STDS_INTG_SPEX_KEYS = { \
    'M8.0': 12572,  # 0019+4614
    'L0.0': 11113,  # 1552+2948
    'L1.0': 10845,  # 0227-1624 - proposed by Cruz et al. 2017
    'L2.0': 11304,  # 0602+3910
    'L3.0': 11070,  # 1726+1538
    'L6.0': 10678}  # 0103+1935


# containters for spectra
#STDS_DWARF_SPEX = {}
#STDS_SD_SPEX = {}
#STDS_ESD_SPEX = {}

# filters
# this information is from the SVO filter profile service: http://svo2.cab.inta-csic.es/svo/theory/fps/
FILTERS = { \
    '2MASS_J': {'file': 'j_2mass.txt', 'description': '2MASS J-band', 'zeropoint': 1594.0, 'method': 'vega', 'rsr': True}, \
    '2MASS_H': {'file': 'h_2mass.txt', 'description': '2MASS H-band', 'zeropoint': 1024.0, 'method': 'vega', 'rsr': True}, \
    '2MASS_KS': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega', 'rsr': True}, \
#    '2MASS_K': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega'}, \
#    '2MASS_Ks': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7, 'method': 'vega'}, \
    'BESSEL_I': {'file': 'bessel_i.txt', 'description': 'Bessel I-band', 'zeropoint': 2405.3, 'method': 'vega', 'rsr': False}, \
    'HAWK_Y': {'file': 'hawk-y.txt', 'description': 'HAWK Y-band', 'zeropoint': 2092.9, 'method': 'vega', 'rsr': False}, \
    'HAWK_J': {'file': 'hawk-j.txt', 'description': 'HAWK J-band', 'zeropoint': 1543.5, 'method': 'vega', 'rsr': False}, \
    'HAWK_H': {'file': 'hawk-h.txt', 'description': 'HAWK H-band', 'zeropoint': 1053.6, 'method': 'vega', 'rsr': False}, \
    'HAWK_H2': {'file': 'hawk-h2.txt', 'description': 'HAWK H2-band', 'zeropoint': 688.8, 'method': 'vega', 'rsr': False}, \
    'HAWK_CH4': {'file': 'hawk-ch4.txt', 'description': 'HAWK CH4-band', 'zeropoint': 1093.4, 'method': 'vega', 'rsr': False}, \
    'HAWK_KS': {'file': 'hawk-ks.txt', 'description': 'HAWK Ks-band', 'zeropoint': 675.3, 'method': 'vega', 'rsr': False}, \
    'HAWK_BRG': {'file': 'hawk-brg.txt', 'description': 'HAWK Brackett Gamma', 'zeropoint': 638.9, 'method': 'vega', 'rsr': False}, \
    'HAWK_NB1060': {'file': 'hawk-nb1060.txt', 'description': 'HAWK Narrow Band 1060', 'zeropoint': 2003.27, 'method': 'vega', 'rsr': False}, \
    'HAWK_NB1190': {'file': 'hawk-nb1190.txt', 'description': 'HAWK Narrow Band 1190', 'zeropoint': 1697.50, 'method': 'vega', 'rsr': False}, \
    'HAWK_NB2090': {'file': 'hawk-nb2090.txt', 'description': 'HAWK Narrow Band 2090', 'zeropoint': 706.68, 'method': 'vega', 'rsr': False}, \
    'FOURSTAR_J': {'file': 'fourstar-j.txt', 'description': 'FOURSTAR J-band', 'zeropoint': 1581.2, 'method': 'vega', 'rsr': False}, \
    'FOURSTAR_J1': {'file': 'fourstar-j1.txt', 'description': 'FOURSTAR J1-band', 'zeropoint': 1978.7, 'method': 'vega', 'rsr': False}, \
    'FOURSTAR_J2': {'file': 'fourstar-j2.txt', 'description': 'FOURSTAR J2-band', 'zeropoint': 1774.5, 'method': 'vega', 'rsr': False}, \
    'FOURSTAR_J3': {'file': 'fourstar-j3.txt', 'description': 'FOURSTAR J3-band', 'zeropoint': 1488.8, 'method': 'vega', 'rsr': False}, \
    'FOURSTAR_H': {'file': 'fourstar-h.txt', 'description': 'FOURSTAR H-band', 'zeropoint': 1054.9, 'method': 'vega', 'rsr': False}, \
    'FOURSTAR_H_SHORT': {'file': 'fourstar-hshort.txt', 'description': 'FOURSTAR H short', 'zeropoint': 1119.1, 'method': 'vega', 'rsr': False}, \
    'FOURSTAR_H_LONG': {'file': 'fourstar-hlong.txt', 'description': 'FOURSTAR H long', 'zeropoint': 980.7, 'method': 'vega', 'rsr': False}, \
    'FOURSTAR_K': {'file': 'fourstar-j.txt', 'description': 'FOURSTAR Ks-band', 'zeropoint': 675.7, 'method': 'vega', 'rsr': False}, \
    'FOURSTAR_KS': {'file': 'fourstar-j.txt', 'description': 'FOURSTAR Ks-band', 'zeropoint': 675.7, 'method': 'vega', 'rsr': False}, \
    'IRAC_CH1': {'file': 'irac1.txt', 'description': 'IRAC Channel 1 (3.6 micron)', 'zeropoint': 280.9, 'method': 'vega', 'rsr': True}, \
    'IRAC_CH2': {'file': 'irac2.txt', 'description': 'IRAC Channel 2 (4.5 micron)', 'zeropoint': 179.7, 'method': 'vega', 'rsr': True}, \
    'IRAC_CH3': {'file': 'irac3.txt', 'description': 'IRAC Channel 3 (5.8 micron)', 'zeropoint': 115.0, 'method': 'vega', 'rsr': True}, \
    'IRAC_CH4': {'file': 'irac4.txt', 'description': 'IRAC Channel 4 (8.0 micron)', 'zeropoint': 64.13, 'method': 'vega', 'rsr': True}, \
    'MKO_J_ATM': {'file': 'j_atm_mko.txt', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3, 'method': 'vega', 'rsr': False}, \
    'MKO_H_ATM': {'file': 'h_atm_mko.txt', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9, 'method': 'vega', 'rsr': False}, \
    'MKO_K_ATM': {'file': 'k_atm_mko.txt', 'description': 'MKO K-band + atmosphere', 'zeropoint': 647.7, 'method': 'vega', 'rsr': False}, \
    'MKO_J': {'file': 'mko_j.txt', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3, 'method': 'vega', 'rsr': False}, \
    'MKO_H': {'file': 'mko_h.txt', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9, 'method': 'vega', 'rsr': False}, \
    'MKO_K': {'file': 'mko_ks.txt', 'description': 'MKO K-band', 'zeropoint': 647.7, 'method': 'vega', 'rsr': False}, \
    'MKO_KP': {'file': 'mko_kp.txt', 'description': 'MKO Kp-band', 'zeropoint': 693.7, 'method': 'vega', 'rsr': False}, \
    'MKO_LP': {'file': 'mko_lp.txt', 'description': 'MKO Lp-band', 'zeropoint': 248.3, 'method': 'vega', 'rsr': False}, \
    'MKO_MP': {'file': 'mko_mp.txt', 'description': 'MKO Mp-band', 'zeropoint': 164.7, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F090M': {'file': 'nic1_f090m.txt', 'description': 'NICMOS F090M', 'zeropoint': 2255.0, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F095N': {'file': 'nic1_f095n.txt', 'description': 'NICMOS F095N', 'zeropoint': 2044.6, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F097N': {'file': 'nic1_f097n.txt', 'description': 'NICMOS F097N', 'zeropoint': 2275.4, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F108N': {'file': 'nic1_f108n.txt', 'description': 'NICMOS F108N', 'zeropoint': 1937.3, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F110M': {'file': 'nic1_f110m.txt', 'description': 'NICMOS F110M', 'zeropoint': 1871.8, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F110W': {'file': 'nic1_f110w.txt', 'description': 'NICMOS F110W', 'zeropoint': 1768.5, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F113N': {'file': 'nic1_f113n.txt', 'description': 'NICMOS F113N', 'zeropoint': 1821.0, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F140W': {'file': 'nic1_f140w.txt', 'description': 'NICMOS F140W', 'zeropoint': 1277.1, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F145M': {'file': 'nic1_f145m.txt', 'description': 'NICMOS F145M', 'zeropoint': 1242.0, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F160W': {'file': 'nic1_f160w.txt', 'description': 'NICMOS F160W', 'zeropoint': 1071.7, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F164N': {'file': 'nic1_f164n.txt', 'description': 'NICMOS F164N', 'zeropoint': 1003.0, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F165M': {'file': 'nic1_f165m.txt', 'description': 'NICMOS F165M', 'zeropoint': 1023.6, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F166N': {'file': 'nic1_f166n.txt', 'description': 'NICMOS F166N', 'zeropoint': 1047.7, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F170M': {'file': 'nic1_f170m.txt', 'description': 'NICMOS F170M', 'zeropoint': 979.1, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F187N': {'file': 'nic1_f187n.txt', 'description': 'NICMOS F187N', 'zeropoint': 803.7, 'method': 'vega', 'rsr': False}, \
    'NICMOS_F190N': {'file': 'nic1_f190n.txt', 'description': 'NICMOS F190N', 'zeropoint': 836.5, 'method': 'vega', 'rsr': False}, \
    'NIRC2_J': {'file': 'nirc2-j.txt', 'description': 'NIRC2 J-band', 'zeropoint': 1562.7, 'method': 'vega', 'rsr': False}, \
    'NIRC2_H': {'file': 'nirc2-h.txt', 'description': 'NIRC2 H-band', 'zeropoint': 1075.5, 'method': 'vega', 'rsr': False}, \
    'NIRC2_HCONT': {'file': 'nirc2-hcont.txt', 'description': 'NIRC2 H-continuum band', 'zeropoint': 1044.5, 'method': 'vega', 'rsr': False}, \
    'NIRC2_K': {'file': 'nirc2-k.txt', 'description': 'NIRC2 K-band', 'zeropoint': 648.9, 'method': 'vega', 'rsr': False}, \
    'NIRC2_KP': {'file': 'nirc2-kp.txt', 'description': 'NIRC2 Kp-band', 'zeropoint': 689.3, 'method': 'vega', 'rsr': False}, \
    'NIRC2_KS': {'file': 'nirc2-ks.txt', 'description': 'NIRC2 Ks-band', 'zeropoint': 676.2, 'method': 'vega', 'rsr': False}, \
    'NIRC2_KCONT': {'file': 'nirc2-kcont.txt', 'description': 'NIRC2 K continuum-band', 'zeropoint': 605.9, 'method': 'vega', 'rsr': False}, \
    'NIRC2_FE2': {'file': 'nirc2-fe2.txt', 'description': 'NIRC2 Fe II', 'zeropoint': 1019.7, 'method': 'vega', 'rsr': False}, \
    'NIRC2_LP': {'file': 'nirc2-lp.txt', 'description': 'NIRC2 LP', 'zeropoint': 248.0, 'method': 'vega', 'rsr': False}, \
    'NIRC2_M': {'file': 'nirc2-ms.txt', 'description': 'NIRC2 M', 'zeropoint': 165.8, 'method': 'vega', 'rsr': False}, \
    'PANSTARRS_G': {'file': 'panstarrs-g.txt', 'description': 'PANSTARRS g-band', 'zeropoint': 3909.11, 'method': 'ab', 'rsr': False}, \
    'PANSTARRS_R': {'file': 'panstarrs-r.txt', 'description': 'PANSTARRS r-band', 'zeropoint': 3151.44, 'method': 'ab', 'rsr': False}, \
    'PANSTARRS_W': {'file': 'panstarrs-w.txt', 'description': 'PANSTARRS w-band', 'zeropoint': 3024.76, 'method': 'ab', 'rsr': False}, \
    'PANSTARRS_I': {'file': 'panstarrs-i.txt', 'description': 'PANSTARRS i-band', 'zeropoint': 2584.6, 'method': 'ab', 'rsr': False}, \
    'PANSTARRS_Z': {'file': 'panstarrs-z.txt', 'description': 'PANSTARRS z-band', 'zeropoint': 2273.09, 'method': 'ab', 'rsr': False}, \
    'PANSTARRS_Y': {'file': 'panstarrs-y.txt', 'description': 'PANSTARRS y-band', 'zeropoint': 2205.95, 'method': 'ab', 'rsr': False}, \
    'SDSS_U': {'file': 'sdss-u.txt', 'description': 'SDSS u-band', 'zeropoint': 1568.5, 'method': 'ab', 'rsr': False}, \
    'SDSS_G': {'file': 'sdss-g.txt', 'description': 'SDSS g-band', 'zeropoint': 3965.9, 'method': 'ab', 'rsr': False}, \
    'SDSS_R': {'file': 'sdss-r.txt', 'description': 'SDSS r-band', 'zeropoint': 3162.0, 'method': 'ab', 'rsr': False}, \
    'SDSS_I': {'file': 'sdss-i.txt', 'description': 'SDSS i-band', 'zeropoint': 2602.0, 'method': 'ab', 'rsr': False}, \
    'SDSS_Z': {'file': 'sdss-z.txt', 'description': 'SDSS z-band', 'zeropoint': 2244.7, 'method': 'ab', 'rsr': False}, \
    'UKIDSS_Z': {'file': 'ukidss-z.txt', 'description': 'UKIDSS Z-band', 'zeropoint': 2261.4, 'method': 'vega', 'rsr': False}, \
    'UKIDSS_Y': {'file': 'ukidss-y.txt', 'description': 'UKIDSS Y-band', 'zeropoint': 2057.2, 'method': 'vega', 'rsr': False}, \
    'UKIDSS_J': {'file': 'ukidss-j.txt', 'description': 'UKIDSS J-band', 'zeropoint': 1556.8, 'method': 'vega', 'rsr': False}, \
    'UKIDSS_H': {'file': 'ukidss-h.txt', 'description': 'UKIDSS H-band', 'zeropoint': 1038.3, 'method': 'vega', 'rsr': False}, \
    'UKIDSS_K': {'file': 'ukidss-k.txt', 'description': 'UKIDSS K-band', 'zeropoint': 644.1, 'method': 'vega', 'rsr': False}, \
    'VISTA_Z': {'file': 'vista_z.txt', 'description': 'VISTA Z-band', 'zeropoint': 2263.81, 'method': 'vega', 'rsr': False}, \
    'VISTA_Y': {'file': 'vista_y.txt', 'description': 'VISTA Y-band', 'zeropoint': 2087.32, 'method': 'vega', 'rsr': False}, \
    'VISTA_J': {'file': 'vista_j.txt', 'description': 'VISTA J-band', 'zeropoint': 1554.03, 'method': 'vega', 'rsr': False}, \
    'VISTA_H': {'file': 'vista_h.txt', 'description': 'VISTA H-band', 'zeropoint': 1030.40, 'method': 'vega', 'rsr': False}, \
    'VISTA_KS': {'file': 'vista_ks.txt', 'description': 'VISTA Ks-band', 'zeropoint': 674.83, 'method': 'vega', 'rsr': False}, \
    'WFC3_F127M': {'file': 'wfc3_F127M.txt', 'description': 'WFC3 F127M', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False}, \
    'WFC3_F139M': {'file': 'wfc3_F139M.txt', 'description': 'WFC3 F139M', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False}, \
    'WFC3_F164N': {'file': 'wfc3_F164N.txt', 'description': 'WFC3 F164N', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False}, \
    'WFC3_F167N': {'file': 'wfc3_F167N.txt', 'description': 'WFC3 F167N', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False}, \
    'WFCAM_Z': {'file': 'wfcam-z.txt', 'description': 'UKIRT WFCAM Z', 'zeropoint': 2261.3, 'method': 'vega', 'rsr': False}, \
    'WFCAM_Y': {'file': 'wfcam-y.txt', 'description': 'UKIRT WFCAM Y', 'zeropoint': 2040.9, 'method': 'vega', 'rsr': False}, \
    'WFCAM_J': {'file': 'wfcam-j.txt', 'description': 'UKIRT WFCAM J', 'zeropoint': 1548.7, 'method': 'vega', 'rsr': False}, \
    'WFCAM_H': {'file': 'wfcam-h.txt', 'description': 'UKIRT WFCAM H', 'zeropoint': 1027.1, 'method': 'vega', 'rsr': False}, \
    'WFCAM_H2': {'file': 'wfcam-h2.txt', 'description': 'UKIRT WFCAM H2', 'zeropoint': 677.1, 'method': 'vega', 'rsr': False}, \
    'WFCAM_BRG': {'file': 'wfcam-brg.txt', 'description': 'UKIRT WFCAM Brackett Gamma', 'zeropoint': 645.5, 'method': 'vega', 'rsr': False}, \
    'WFCAM_K': {'file': 'wfcam-k.txt', 'description': 'UKIRT WFCAM K', 'zeropoint': 630.0, 'method': 'vega', 'rsr': False}, \
    'WIRCAM_Y': {'file': 'wircam-cfht-y.txt', 'description': 'CFHT WIRCAM Y', 'zeropoint': 2073.32, 'method': 'vega', 'rsr': False}, \
    'WIRCAM_J': {'file': 'wircam-cfht-j.txt', 'description': 'CFHT WIRCAM J', 'zeropoint': 1551.01, 'method': 'vega', 'rsr': False}, \
    'WIRCAM_H': {'file': 'wircam-cfht-h.txt', 'description': 'CFHT WIRCAM H', 'zeropoint': 1044.35, 'method': 'vega', 'rsr': False}, \
    'WIRCAM_KS': {'file': 'wircam-cfht-ks.txt', 'description': 'CFHT WIRCAM Ks', 'zeropoint': 674.62, 'method': 'vega', 'rsr': False}, \
    'WIRCAM_KCONT': {'file': 'wircam-cfht-kcont.txt', 'description': 'CFHT WIRCAM K-cont', 'zeropoint': 636.17, 'method': 'vega', 'rsr': False}, \
    'WIRCAM_CH4_OFF': {'file': 'wircam-cfht-ch4s.txt', 'description': 'CFHT WIRCAM CH4-off', 'zeropoint': 987.39, 'method': 'vega', 'rsr': False}, \
    'WIRCAM_CH4_ON': {'file': 'wircam-cfht-ch4l.txt', 'description': 'CFHT WIRCAM CH4-on', 'zeropoint': 1076.31, 'method': 'vega', 'rsr': False}, \
    'WIRC_J': {'file': 'wirc_jcont.txt', 'description': 'WIRC J-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False}, \
    'WIRC_H': {'file': 'wirc_hcont.txt', 'description': 'WIRC H-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False}, \
    'WIRC_K': {'file': 'wirc_kcont.txt', 'description': 'WIRC K-cont', 'zeropoint': 0., 'method': 'vega', 'rsr': False}, \
    'WIRC_CO': {'file': 'wirc_co.txt', 'description': 'WIRC CO', 'zeropoint': 0., 'method': 'vega', 'rsr': False}, \
    'WIRC_CH4S': {'file': 'wirc_ch4s.txt', 'description': 'WIRC CH4S', 'zeropoint': 0., 'method': 'vega', 'rsr': False}, \
    'WIRC_CH4L': {'file': 'wirc_ch4l.txt', 'description': 'WIRC CH4L', 'zeropoint': 0., 'method': 'vega', 'rsr': False}, \
    'WIRC_FE2': {'file': 'wirc_feii.txt', 'description': 'WIRC Fe II', 'zeropoint': 0., 'method': 'vega', 'rsr': False}, \
    'WIRC_BRGAMMA': {'file': 'wirc_brgamma.txt', 'description': 'WIRC H I Brackett Gamma', 'zeropoint': 0., 'method': 'vega', 'rsr': False}, \
    'WIRC_PABETA': {'file': 'wirc_pabeta.txt', 'description': 'WIRC H I Paschen Beta', 'zeropoint': 0., 'method': 'vega', 'rsr': False}, \
    'WISE_W1': {'file': 'wise_w1.txt', 'description': 'WISE W1 (3.5 micron)', 'zeropoint': 309.54, 'method': 'vega', 'rsr': True}, \
    'WISE_W2': {'file': 'wise_w2.txt', 'description': 'WISE W2 (4.6 micron)', 'zeropoint': 171.79, 'method': 'vega', 'rsr': True}, \
    'WISE_W3': {'file': 'wise_w3.txt', 'description': 'WISE W3 (13 micron)', 'zeropoint': 31.67, 'method': 'vega', 'rsr': True}, \
    'WISE_W4': {'file': 'wise_w4.txt', 'description': 'WISE W4 (22 micron)', 'zeropoint': 8.363, 'method': 'vega', 'rsr': True} \
    }
VEGAFILE = 'vega_kurucz.txt'

# some data formats (for future expansion)
INSTRUMENTS = {
#	'SPEX': {'instrument_name': 'SpeX prism', 'pixelscale': 0.15*u.arcsec, 'waverange': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['']},
    'SED': {'instrument_name': 'SED', 'pixelscale': 0.*u.arcsec, 'waverange': [0.1,100.]*u.micron, 'slitwidth': 0.*u.arcsec, 'resolution': 100, 'norders': 1, 'readnoise': 0, 'darkcurrent': 0., 'gain': 1, 'altnames': ['SPECTRAL_ENERGY_DISTRIBUTION']},
    'SPEX_PRISM': {'instrument_name': 'SpeX prism', 'pixelscale': 0.15*u.arcsec, 'waverange': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['SPEX']},
    'SPEX_SXD': {'instrument_name': 'SpeX SXD', 'pixelscale': 0.15*u.arcsec, 'waverange': [0.8,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['']},
    'SPEX_LXD1.9': {'instrument_name': 'SpeX LXD 1.9 micron', 'pixelscale': 0.15*u.arcsec, 'waverange': [1.95,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['']},
    'SPEX_LXD2.1': {'instrument_name': 'SpeX LXD 2.1 micron', 'pixelscale': 0.15*u.arcsec, 'waverange': [2.15,5.0]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['']},
    'SPEX_LXD2.3': {'instrument_name': 'SpeX LXD 2.3 micron', 'pixelscale': 0.15*u.arcsec, 'waverange': [2.25,5.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 12, 'darkcurrent': 0.2, 'gain': 12, 'altnames': ['']},
#	'USPEX': {'instrument_name': 'Updated SpeX prism', 'pixelscale': 0.10*u.arcsec, 'waverange': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['']},
	'USPEX_PRISM': {'instrument_name': 'Updated SpeX prism', 'pixelscale': 0.10*u.arcsec, 'waverange': [0.7,2.5]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 200, 'norders': 1, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['USPEX']},
    'USPEX_SXD': {'instrument_name': 'Updated SpeX SXD', 'pixelscale': 0.10*u.arcsec, 'waverange': [0.7,2.55]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2000, 'norders': 7, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['']},
    'USPEX_LXD_SHORT': {'instrument_name': 'Updated SpeX LXD short', 'pixelscale': 0.10*u.arcsec, 'waverange': [1.67,4.2]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 8, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['']},
    'USPEX_LXD_LONG': {'instrument_name': 'Updated SpeX LXD long', 'pixelscale': 0.10*u.arcsec, 'waverange': [1.98,5.3]*u.micron, 'slitwidth': 0.3*u.arcsec, 'resolution': 2500, 'norders': 7, 'readnoise': 5, 'darkcurrent': 0.05, 'gain': 1.5, 'altnames': ['']},
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
    'burrows06': {'name': 'Burrows et al. (2006)', 'bibcode': '2006ApJ...640.1063B', 'altnames': ['burrows','burrows2006'], 'rawfolder': HOME_FOLDER+'/models/burrows/burrows06/'}, \
    'btsettl08': {'name': 'BT-Settl (2008)', 'bibcode': '2012RSPTA.370.2765A', 'altnames': ['allard','allard12','allard2012','btsettl','btsettled','btsettl08','btsettl2008','BTSettl2008'], 'rawfolder': HOME_FOLDER+'/models/allard/cifist2011/'}, \
    'btsettl15': {'name': 'BT-Settl (2015)', 'bibcode': '2015A&A...577A..42B', 'altnames': ['allard15','allard2015','btsettl015','btsettl2015','BTSettl2015'], 'rawfolder': HOME_FOLDER+'/models/allard/cifist2015/BT-Settl_M-0.0a+0.0/'}, \
    'madhusudhan11': {'name': 'Madhusudhan et al. (2011)', 'bibcode': '2011ApJ...737...34M', 'altnames': ['madhu','madhusudhan','madhu11','madhu2011','madhusudhan2011'], 'rawfolder': HOME_FOLDER+'/models/burrows/burrows11/all/'}, \
    'morley12': {'name': 'Morley et al. (2012)', 'bibcode': '2012ApJ...756..172M', 'altnames': ['morley','morley2012'], 'rawfolder': HOME_FOLDER+'/models/ames/Morley12/'}, \
    'morley14': {'name': 'Morley et al. (2014)', 'bibcode': '2014ApJ...787...78M', 'altnames': ['morley2014'], 'rawfolder': HOME_FOLDER+'/models/ames/Morley14/'}, \
    'saumon12': {'name': 'Saumon et al. (2012)', 'bibcode': '2012ApJ...750...74S', 'altnames': ['saumon','saumon2012'], 'rawfolder': HOME_FOLDER+'/models/ames/Saumon12/'}, \
    'drift': {'name': 'Witte et al. (2011)', 'bibcode': '2011A&A...529A..44W', 'altnames': ['witte','witte11','witte2011','helling'], 'rawfolder': HOME_FOLDER+'/models/drift/v1.2/'}}
SPECTRAL_MODEL_PARAMETERS_INORDER = ['teff','logg','z','fsed','cld','kzz']
SPECTRAL_MODEL_PARAMETERS = {\
    'teff': {'prefix': 't', 'unit': u.K, 'default': 1000.0, 'title': '$T_{eff}$', 'type': 'continuous'}, \
    'logg': {'prefix': 'g', 'unit': u.dex(), 'default': 5.0, 'title': '$\log{g}$', 'type': 'continuous'}, \
    'z': {'prefix': 'z', 'unit': u.dex(), 'default': 0., 'title': '$[M/H]$', 'type': 'continuous'}, \
    'fsed': {'prefix': 'f', 'unit': u.m/u.m, 'default': 'nc', 'title': '$f_{sed}$', 'type': 'discrete'}, \
    'cld': {'prefix': 'c', 'unit': u.m/u.m, 'default': 'nc', 'title': '$cld$', 'type': 'discrete'}, \
    'kzz': {'prefix': 'k', 'unit': u.m/u.m, 'default': 'eq', 'title': '$log\ \kappa_{zz}$', 'type': 'discrete'}}
#    'kzz': u.dex(u.cm*u.cm/u.s), \
#    'slit': {'unit': u.arcsec, 'default': 0.5, 'title': 'slit'}}
SPECTRAL_MODEL_FLUX_UNIT = u.erg/(u.s*u.micron*(u.cm**2))
SPECTRAL_MODEL_SED_UNIT = u.erg/(u.s*(u.cm**2))
SPECTRAL_MODEL_WAVE_UNIT = u.micron
#DEFINED_MODEL_BIBCODES = {\
#    'BTSettl2008': '', \
#    'burrows06': '2006ApJ...640.1063B',\
#    'morley12': '2012ApJ...756..172M',\
#    'morley14': '2014ApJ...787...78M',\
#    'saumon12': '2012ApJ...750...74S',\
#    'drift': '2011A&A...529A..44W'}

# evolutionary model information
EVOLUTIONARY_MODELS = {\
    'baraffe03': {'name': 'Baraffe et al. (2003)', 'bibcode': '2003A&A...402..701B', 'altnames': ['bar03','baraffe2003']},\
    'baraffe15': {'name': 'Baraffe et al. (2015)', 'bibcode': '2015A&A...577A..42B', 'altnames': ['baraffe','bar15','baraffe2015']},\
    'burrows01': {'name': 'Burrows et al. (2001)', 'bibcode': '2015A&A...577A..42B', 'altnames': ['burrows','bur01','burrows2001']},\
    'saumon08': {'name': 'Saumon et al. (2008)', 'bibcode': '2015A&A...577A..42B', 'altnames': ['saumon','sau08','saumon2008']}}

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
    'faherty2012': {'altname': ['faherty12'],'bibcode': '2012ApJ...752...56F', 'sptoffset': 10, 'filters': {
        'MKO J': {'fitunc' : 0.30, 'range' : [20., 38.], 'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, 
        'MKO H': {'fitunc' : 0.27, 'range' : [20., 38.], 'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, 
        'MKO K': {'fitunc' : 0.28, 'range' : [20., 38.], 'coeff' : [.0000816516, -.00469032, .0940816, -.485519, 9.76100]}}},
    'burgasser2007': {'altname': ['burgasser','burgasser07'], 'bibcode': '2007ApJ...659..655B', 'sptoffset': 20, 'filters': {
        'MKO J': {'fitunc' : 0.30, 'range' : [20., 38.], 'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, 
        'MKO H': {'fitunc' : 0.27, 'range' : [20., 38.], 'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, 
        'MKO K': {'fitunc' : 0.26, 'range' : [20., 38.], 'coeff': [.0000001051, -.000006985, .0001807, -.002271, .01414, -.04024, .05129, .2322, 10.45]}}},
    'dupuy2012': {'altname': ['dupuy','dupuy12'], 'bibcode': '2012ApJS..201...19D', 'sptoffset': 10, 'filters': {
        'MKO J': {'fitunc' : 0.39, 'range' : [16., 39.], 'coeff' : [-.00000194920, .000227641, -.0103332, .232771, -2.74405, 16.3986, -28.3129]}, 
        'MKO Y': {'fitunc': 0.40, 'range' : [16., 39.], 'coeff': [-.00000252638, .000285027, -.0126151, .279438, -3.26895, 19.5444, -35.1560]}, 
        'MKO H': {'fitunc': 0.38, 'range' : [16., 39.], 'coeff': [-.00000224083, .000251601, -.0110960, .245209, -2.85705, 16.9138, -29.7306]}, 
        'MKO K': {'fitunc': 0.40, 'range' : [16., 39.], 'coeff': [-.00000104935, .000125731, -.00584342, .135177, -1.63930, 10.1248, -15.2200]}, 
        'MKO LP': {'fitunc': 0.28, 'range': [16., 39.], 'coeff': [0.00000, 0.00000, .0000546366, -.00293191, .0530581,  -.196584, 8.89928]}, 
        '2MASS J': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [-.000000784614, .000100820, -.00482973, .111715, -1.33053, 8.16362, -9.67994]}, 
        '2MASS H': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [-.00000111499, .000129363, -.00580847, .129202, -1.50370, 9.00279, -11.7526]}, 
        '2MASS KS': {'fitunc': 0.43, 'range':[16., 39.], 'coeff': [1.06693e-4, -6.42118e-3, 1.34163e-1, -8.67471e-1, 1.10114e1]}, 
        'WISE W1': {'fitunc': 0.39, 'range':[16., 39.], 'coeff': [1.58040e-5, -3.33944e-4, -4.38105e-3, 3.55395e-1, 7.14765]}, 
        'WISE W2': {'fitunc': 0.35, 'range':[16., 39.], 'coeff': [1.78555e-5, -8.81973e-4, 1.14325e-2, 1.92354e-1, 7.46564]}}},
    'filippazzo2015': {'altname': ['filippazzo','filippazzo15'],'bibcode': '2015ApJ...810..158F', 'sptoffset': 10, 'filters': {
        '2MASS J': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [3.478e-5, -2.684e-3, 7.771e-2, -1.058e0, 7.157e0, -8.350e0]}, 
        'WISE W2': {'fitunc': 0.40, 'range': [16., 39.], 'coeff': [8.190e-6, -6.938e-4, 2.283e-2, -3.655e-1, 3.032e0, -5.043e-1]}}},
    'faherty2016': {'altname': ['faherty','faherty2016,faherty-field'],'bibcode': '2016ApJS..225...10F', 'sptoffset': 10, 'filters': {
        '2MASS J': {'fitunc' : 0.402, 'range' : [16., 39.], 'coeff': [3.478e-5,-2.684e-3,7.771e-2,-1.058e0,7.157e0,-8.350e0]},
        '2MASS H': {'fitunc' : 0.389, 'range' : [16., 39.], 'coeff' : [2.841e-5,-2.217e-3,6.551e-2,-9.174e-1,6.406e0,-7.496e0]},
        '2MASS KS': {'fitunc' : 0.537, 'range' : [16., 39.], 'coeff' : [2.540e-05,-1.997e-03,5.978e-02,-8.481e-01,5.970e+00,-6.704e+00]}}}
}


#    7.0 < SpT < 17.0    0.647   4.032e-03   -1.416e-01  2.097e+00   8.478e-01   cdots   cdots   cdots   
#    7.0 < SpT < 17.0    0.660   -3.825e-03  1.370e-01   -9.279e-01  10.141e+00  cdots   cdots   cdots   
#   7.0 < SpT < 17.0    0.634   2.642e-03   -1.049e-01  1.753e+00   1.207e+00   cdots   cdots   cdots   
#   7.0 < SpT < 17.0    0.603   -3.909e-03  1.346e-01   -9.347e-01  9.728e+00   cdots   cdots   cdots   
#   7.0 < SpT < 17.0    0.640   -1.585e-02  7.338e-01   4.537e+00   cdots   cdots   cdots   cdots   
#   7.0 < SpT < 17.0    0.556   -4.006e-03  1.378e-01   -1.031e+00  9.916e+00   cdots   cdots   cdots   
#   6.0 < SpT < 29.0    0.365   -1.664e-01  2.991e+00   -3.603e-01  2.258e-02   -6.897e-04  8.337e-06   cdots   
#   7.0 < SpT < 17.0    0.648   -1.397e-02  5.955e-01   5.247e+00   cdots   cdots   cdots   cdots   
#   7.0 < SpT < 17.0    0.551   -4.483e-03  1.505e-01   -1.208e+00  10.403e+00  cdots   cdots   cdots   
#   6.0 < SpT < 29.0    0.398   -5.043e-01  3.032e+00   -3.655e-01  2.283e-02   -6.938e-04  8.190e-06   cdots   
#   7.0 < SpT < 17.0    0.694   -1.507e-02  5.944e-01   5.061e+00   cdots   cdots   cdots   cdots   
#   7.0 < SpT < 17.0    0.616   -6.821e-03  2.322e-01   -2.133e+00  13.322e+00  cdots   cdots   cdots   
#   6.0 < SpT < 29.0    0.446   6.462e+00   3.365e-01   1.520e-02   -2.573e-03  9.477e-05   -1.024e-06  cdots   
#   7.0 < SpT < 17.0    0.717   -1.003e-04  -1.670e-03  2.023e-01   7.529e+00   cdots   cdots   cdots   
#   7.0 < SpT < 17.0    0.427   -5.684e-03  1.993e-01   -1.987e+00  13.972e+00  cdots   cdots   cdots   
#   6.0 < SpT < 29.0    113.431 4.747e+03   -7.005e+02  1.155e+02   -1.191e+01  6.318e-01   -1.606e-02  1.546e-04   
#   7.0 < SpT < 17.0    180.457 1.330e+00   -66.8637    1235.42 -10068.8    32766.4 cdots   cdots   
#   7.0 < SpT < 28.0    197.737 2.795e+04   -9.183e+03  1.360e+03   -1.066e+02  4.578e+00   -1.016e-01  9.106e-04   
#   7.0 < SpT < 17.0    172.215 7.383e+00   -344.522    4879.86 cdots   cdots   cdots   cdots   
#   7.0 < SpT < 28.0    0.133   2.787e+00   -2.310e+00  3.727e-01   -3.207e-02  1.449e-03   -3.220e-05  2.736e-07   
#   7.0 < SpT < 17.0    0.335   -6.514e-03  2.448e-01   -3.113e+00  9.492e+00   cdots   cdots   cdots   
#   7.0 < SpT < 28.0    0.206   2.059e-01   9.585   -3.985  4.923e-01   -3.048e-02  9.134e-04   -1.056e-05  
#   7.0 < SpT < 17.0    0.221   6.194e-03   -3.757e-01  2.728e-02   cdots   cdots   cdots   cdots   
#   7.0 < SpT < 28.0        -3.46623e-01    3.40366e-02 -3.072e-03  

