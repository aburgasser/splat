# -*- coding: utf-8 -*-
from __future__ import print_function, division

# WORKING COPY OF SPLAT CODE LIBRARY
# based on routines developed by:
#    Christian Aganze
#    Daniella Bardalez Gagliuffi
#    Adam Burgasser
#    Caleb Choban
#    Ivanna Escala
#    Aishwarya Iyer
#    Yuhui Jin
#    Michael Lopez
#    Alex Mendez
#    Gretel Mercado
#    Jonathan Parra
#    Maitrayee Sahi
#    Adrian Suarez
#    Melisa Tallis
#    Tomoki Tamiya

#
# CURRENT STATUS (3/12/2015)
# URGENT
# fails when reading in unpublished (online) data
#
# LESS URGENT
# reformat fits files so headers have valid information, and spectra are flux calibrated
# update information for sources in source database
# add help sections to all programs
# plotspectrum => multiplot, multipage files
# proper SQL search for database
# have classifybyindex return individual index classifications (e.g., 'full' keyword)

# verify the version is correct
import sys
if sys.version_info.major != 2 and sys.version_info.major != 3:
    raise NameError('\nSPLAT only works on Python 2.7 and 3.X\n')

# imports
import astropy
import base64
import copy
import matplotlib.pyplot as plt
import numpy
import os
import random
import re
import requests
import scipy
if sys.version_info.major == 2:     # switch for those using python 3
    import string
import warnings

from astropy.io import fits            # for reading in spreadsheet
from astropy.table import Table, join            # for reading in table files
from astropy.coordinates import SkyCoord      # coordinate conversion
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units
from scipy import stats, signal
from scipy.integrate import trapz        # for numerical integration
from scipy.interpolate import interp1d

# suppress warnings - probably not an entirely safe approach!
numpy.seterr(all='ignore')
warnings.simplefilter("ignore")
#from splat._version import __version__

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

# local application/library specific import
#import bdevopar as splevol
from splat_db import *
from splat_model import *
from splat_plot import *
#import splat_db

#################### CONSTANTS ####################
SPLAT_URL = 'http://pono.ucsd.edu/~adam/splat/'
DATA_FOLDER = '/reference/Spectra/'

# explicitly read in source and spectral databases
DB_SOURCES = fetchDatabase(splat.DB_SOURCES_FILE)
DB_SPECTRA = fetchDatabase(splat.DB_SPECTRA_FILE)

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
spex_pixel_scale = 0.15            # spatial scale in arcseconds per pixel
uspex_pixel_scale = 0.10            # spatial scale in arcseconds per pixel
spex_wave_range = [0.65,2.45]*u.micron    # default wavelength range
max_snr = 1000.0                # maximum S/N ratio permitted
TMPFILENAME = 'splattmpfile'

SPEX_STDFILES = { \
    'M0.0': '11335_10505.fits',\
    'M1.0': '11364_10806.fits',\
    'M2.0': '11181_10187.fits',\
    'M3.0': '10823_11422.fits',\
    'M4.0': '12004_10444.fits',\
    'M5.0': '10829_10104.fits',\
    'M6.0': '11182_10188.fits',\
    'M7.0': '10822_11283.fits',\
    'M8.0': '10824_11423.fits',\
    'M9.0': '10821_11058.fits',\
    'L0.0': '10107_10315.fits',\
    'L1.0': '11072_11527.fits',\
    'L2.0': '10600_10957.fits',\
    'L3.0': '10592_11111.fits',\
    'L4.0': '10675_11572.fits',\
    'L5.0': '10351_10583.fits',\
    'L6.0': '10375_10696.fits',\
    'L7.0': '10678_10105.fits',\
    'L8.0': '10115_11254.fits',\
    'L9.0': '10268_10237.fits',\
    'T0.0': '10771_10871.fits',\
    'T1.0': '10767_10591.fits',\
    'T2.0': '10017_10945.fits',\
    'T3.0': '10034_10874.fits',\
    'T4.0': '10143_11632.fits',\
    'T5.0': '10021_11106.fits',\
    'T6.0': '10200_11236.fits',\
    'T7.0': '10159_10513.fits',\
    'T8.0': '10126_10349.fits',\
    'T9.0': '11536_10509.fits'}
# EMPTY DICTIONARY
SPEX_STDS = {}

SPEX_SD_STDFILES = { \
    'sdM5.5': '11670_11134.fits',\
    'sdM6.0': '10265_10045.fits',\
    'sdM7.0': '10197_11074.fits',\
    'sdM8.0': '10123_10145.fits',\
    'sdM9.5': '10188_10700.fits',\
    'sdL0.0': '11972_10248.fits',\
    'sdL3.5': '10364_10946.fits',\
    'sdL4.0': '10203_11241.fits'}
# EMPTY DICTIONARY
SPEX_SD_STDS = {}

SPEX_ESD_STDFILES = { \
    'esdM5.0': '10229_10163.fits',\
#    'esdM6.5': '_10579.fits',\
    'esdM7.0': '10521_10458.fits',\
    'esdM8.5': '10278_10400.fits'}
# EMPTY DICTIONARY
SPEX_ESD_STDS = {}


# filters
FILTER_FOLDER = '/reference/Filters/'
FILTERS = { \
    '2MASS_J': {'file': 'j_2mass.txt', 'description': '2MASS J-band', 'zeropoint': 1594.0}, \
    '2MASS_H': {'file': 'h_2mass.txt', 'description': '2MASS H-band', 'zeropoint': 1024.0}, \
    '2MASS_KS': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band', 'zeropoint': 666.7}, \
    'BESSEL_I': {'file': 'bessel_i.txt', 'description': 'Bessel I-band', 'zeropoint': 2405.3}, \
    'HAWK_Y': {'file': 'hawk-y.txt', 'description': 'HAWK Y-band', 'zeropoint': 2092.9}, \
    'HAWK_J': {'file': 'hawk-j.txt', 'description': 'HAWK J-band', 'zeropoint': 1543.5}, \
    'HAWK_H': {'file': 'hawk-h.txt', 'description': 'HAWK H-band', 'zeropoint': 1053.6}, \
    'HAWK_H2': {'file': 'hawk-h2.txt', 'description': 'HAWK H2-band', 'zeropoint': 688.8}, \
    'HAWK_CH4': {'file': 'hawk-ch4.txt', 'description': 'HAWK CH4-band', 'zeropoint': 1093.4}, \
    'HAWK_KS': {'file': 'hawk-ks.txt', 'description': 'HAWK Ks-band', 'zeropoint': 675.3}, \
    'HAWK_BRG': {'file': 'hawk-brg.txt', 'description': 'HAWK Brackett Gamma', 'zeropoint': 638.9}, \
    'HAWK_NB1060': {'file': 'hawk-nb1060.txt', 'description': 'HAWK Narrow Band 1060', 'zeropoint': 2003.27}, \
    'HAWK_NB1190': {'file': 'hawk-nb1190.txt', 'description': 'HAWK Narrow Band 1190', 'zeropoint': 1697.50}, \
    'HAWK_NB2090': {'file': 'hawk-nb2090.txt', 'description': 'HAWK Narrow Band 2090', 'zeropoint': 706.68}, \
    'FOURSTAR_J': {'file': 'fourstar-j.txt', 'description': 'FOURSTAR J-band', 'zeropoint': 1581.2}, \
    'FOURSTAR_J1': {'file': 'fourstar-j1.txt', 'description': 'FOURSTAR J1-band', 'zeropoint': 1978.7}, \
    'FOURSTAR_J2': {'file': 'fourstar-j2.txt', 'description': 'FOURSTAR J2-band', 'zeropoint': 1774.5}, \
    'FOURSTAR_J3': {'file': 'fourstar-j3.txt', 'description': 'FOURSTAR J3-band', 'zeropoint': 1488.8}, \
    'FOURSTAR_H': {'file': 'fourstar-h.txt', 'description': 'FOURSTAR H-band', 'zeropoint': 1054.9}, \
    'FOURSTAR_H_SHORT': {'file': 'fourstar-hshort.txt', 'description': 'FOURSTAR H short', 'zeropoint': 1119.1}, \
    'FOURSTAR_H_LONG': {'file': 'fourstar-hlong.txt', 'description': 'FOURSTAR H long', 'zeropoint': 980.7}, \
    'FOURSTAR_KS': {'file': 'fourstar-j.txt', 'description': 'FOURSTAR Ks-band', 'zeropoint': 675.7}, \
    'IRAC_CH1': {'file': 'irac1.txt', 'description': 'IRAC Channel 1 (3.6 micron)', 'zeropoint': 280.9}, \
    'IRAC_CH2': {'file': 'irac2.txt', 'description': 'IRAC Channel 2 (4.5 micron)', 'zeropoint': 179.7}, \
    'IRAC_CH3': {'file': 'irac3.txt', 'description': 'IRAC Channel 3 (5.8 micron)', 'zeropoint': 115.0}, \
    'IRAC_CH4': {'file': 'irac4.txt', 'description': 'IRAC Channel 4 (8.0 micron)', 'zeropoint': 64.13}, \
    'MKO_J_ATM': {'file': 'j_atm_mko.txt', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3}, \
    'MKO_H_ATM': {'file': 'h_atm_mko.txt', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9}, \
    'MKO_K_ATM': {'file': 'k_atm_mko.txt', 'description': 'MKO K-band + atmosphere', 'zeropoint': 647.7}, \
    'MKO_J': {'file': 'mko_j.txt', 'description': 'MKO J-band + atmosphere', 'zeropoint': 1562.3}, \
    'MKO_H': {'file': 'mko_h.txt', 'description': 'MKO H-band + atmosphere', 'zeropoint': 1045.9}, \
    'MKO_K': {'file': 'mko_ks.txt', 'description': 'MKO K-band', 'zeropoint': 647.7}, \
    'MKO_KP': {'file': 'mko_kp.txt', 'description': 'MKO Kp-band', 'zeropoint': 693.7}, \
    'MKO_LP': {'file': 'mko_lp.txt', 'description': 'MKO Lp-band', 'zeropoint': 248.3}, \
    'MKO_MP': {'file': 'mko_mp.txt', 'description': 'MKO Mp-band', 'zeropoint': 164.7}, \
    'NICMOS_F090M': {'file': 'nic1_f090m.txt', 'description': 'NICMOS F090M', 'zeropoint': 2255.0}, \
    'NICMOS_F095N': {'file': 'nic1_f095n.txt', 'description': 'NICMOS F095N', 'zeropoint': 2044.6}, \
    'NICMOS_F097N': {'file': 'nic1_f097n.txt', 'description': 'NICMOS F097N', 'zeropoint': 2275.4}, \
    'NICMOS_F108N': {'file': 'nic1_f108n.txt', 'description': 'NICMOS F108N', 'zeropoint': 1937.3}, \
    'NICMOS_F110M': {'file': 'nic1_f110m.txt', 'description': 'NICMOS F110M', 'zeropoint': 1871.8}, \
    'NICMOS_F110W': {'file': 'nic1_f110w.txt', 'description': 'NICMOS F110W', 'zeropoint': 1768.5}, \
    'NICMOS_F113N': {'file': 'nic1_f113n.txt', 'description': 'NICMOS F113N', 'zeropoint': 1821.0}, \
    'NICMOS_F140W': {'file': 'nic1_f140w.txt', 'description': 'NICMOS F140W', 'zeropoint': 1277.1}, \
    'NICMOS_F145M': {'file': 'nic1_f145m.txt', 'description': 'NICMOS F145M', 'zeropoint': 1242.0}, \
    'NICMOS_F160W': {'file': 'nic1_f160w.txt', 'description': 'NICMOS F160W', 'zeropoint': 1071.7}, \
    'NICMOS_F164N': {'file': 'nic1_f164n.txt', 'description': 'NICMOS F164N', 'zeropoint': 1003.0}, \
    'NICMOS_F165M': {'file': 'nic1_f165m.txt', 'description': 'NICMOS F165M', 'zeropoint': 1023.6}, \
    'NICMOS_F166N': {'file': 'nic1_f166n.txt', 'description': 'NICMOS F166N', 'zeropoint': 1047.7}, \
    'NICMOS_F170M': {'file': 'nic1_f170m.txt', 'description': 'NICMOS F170M', 'zeropoint': 979.1}, \
    'NICMOS_F187N': {'file': 'nic1_f187n.txt', 'description': 'NICMOS F187N', 'zeropoint': 803.7}, \
    'NICMOS_F190N': {'file': 'nic1_f190n.txt', 'description': 'NICMOS F190N', 'zeropoint': 836.5}, \
    'NIRC2_J': {'file': 'nirc2-j.txt', 'description': 'NIRC2 J-band', 'zeropoint': 1562.7}, \
    'NIRC2_H': {'file': 'nirc2-h.txt', 'description': 'NIRC2 H-band', 'zeropoint': 1075.5}, \
    'NIRC2_HCONT': {'file': 'nirc2-hcont.txt', 'description': 'NIRC2 H-continuum band', 'zeropoint': 1044.5}, \
    'NIRC2_K': {'file': 'nirc2-k.txt', 'description': 'NIRC2 K-band', 'zeropoint': 648.9}, \
    'NIRC2_KP': {'file': 'nirc2-kp.txt', 'description': 'NIRC2 Kp-band', 'zeropoint': 689.3}, \
    'NIRC2_KS': {'file': 'nirc2-ks.txt', 'description': 'NIRC2 Ks-band', 'zeropoint': 676.2}, \
    'NIRC2_KCONT': {'file': 'nirc2-kcont.txt', 'description': 'NIRC2 K continuum-band', 'zeropoint': 605.9}, \
    'NIRC2_FE2': {'file': 'nirc2-fe2.txt', 'description': 'WIRC Fe II', 'zeropoint': 1019.7}, \
    'NIRC2_LP': {'file': 'nirc2-lp.txt', 'description': 'WIRC Fe II', 'zeropoint': 248.0}, \
    'NIRC2_M': {'file': 'nirc2-ms.txt', 'description': 'WIRC Fe II', 'zeropoint': 165.8}, \
    'PANSTARRS_I': {'file': 'panstarrs-i.txt', 'description': 'PANSTARRS i-band', 'zeropoint': 2584.6}, \
    'PANSTARRS_Z': {'file': 'panstarrs-z.txt', 'description': 'PANSTARRS z-band', 'zeropoint': 2584.6}, \
    'PANSTARRS_Y': {'file': 'panstarrs-y.txt', 'description': 'PANSTARRS y-band', 'zeropoint': 2584.6}, \
    'UKIDSS_Z': {'file': 'ukidss-z.txt', 'description': 'UKIDSS Z-band', 'zeropoint': 2261.4}, \
    'UKIDSS_Y': {'file': 'ukidss-y.txt', 'description': 'UKIDSS Y-band', 'zeropoint': 2057.2}, \
    'UKIDSS_J': {'file': 'ukidss-j.txt', 'description': 'UKIDSS J-band', 'zeropoint': 1556.8}, \
    'UKIDSS_H': {'file': 'ukidss-h.txt', 'description': 'UKIDSS H-band', 'zeropoint': 1038.3}, \
    'UKIDSS_K': {'file': 'ukidss-k.txt', 'description': 'UKIDSS K-band', 'zeropoint': 644.1}, \
    'VISTA_Z': {'file': 'vista_z.txt', 'description': 'VISTA Z-band', 'zeropoint': 2263.81}, \
    'VISTA_Y': {'file': 'vista_y.txt', 'description': 'VISTA Y-band', 'zeropoint': 2087.32}, \
    'VISTA_J': {'file': 'vista_j.txt', 'description': 'VISTA J-band', 'zeropoint': 1554.03}, \
    'VISTA_H': {'file': 'vista_h.txt', 'description': 'VISTA H-band', 'zeropoint': 1030.40}, \
    'VISTA_KS': {'file': 'vista_ks.txt', 'description': 'VISTA Ks-band', 'zeropoint': 674.83}, \
    'WFC3_F127M': {'file': 'wfc3_F127M.txt', 'description': 'WFC3 F127M', 'zeropoint': 2261.3}, \
    'WFC3_F139M': {'file': 'wfc3_F139M.txt', 'description': 'WFC3 F139M', 'zeropoint': 2261.3}, \
    'WFC3_F164N': {'file': 'wfc3_F164M.txt', 'description': 'WFC3 F164N', 'zeropoint': 2261.3}, \
    'WFC3_F167N': {'file': 'wfc3_F160W.txt', 'description': 'WFC3 F160W', 'zeropoint': 2261.3}, \
    'WFCAM_Z': {'file': 'wfcam-z.txt', 'description': 'UKIRT WFCAM Z', 'zeropoint': 2261.3}, \
    'WFCAM_Y': {'file': 'wfcam-y.txt', 'description': 'UKIRT WFCAM Y', 'zeropoint': 2040.9}, \
    'WFCAM_J': {'file': 'wfcam-j.txt', 'description': 'UKIRT WFCAM J', 'zeropoint': 1548.7}, \
    'WFCAM_H': {'file': 'wfcam-h.txt', 'description': 'UKIRT WFCAM H', 'zeropoint': 1027.1}, \
    'WFCAM_H2': {'file': 'wfcam-h2.txt', 'description': 'UKIRT WFCAM H2', 'zeropoint': 677.1}, \
    'WFCAM_BRG': {'file': 'wfcam-brg.txt', 'description': 'UKIRT WFCAM Brackett Gamma', 'zeropoint': 645.5}, \
    'WFCAM_K': {'file': 'wfcam-k.txt', 'description': 'UKIRT WFCAM K', 'zeropoint': 630.0}, \
    'WIRCAM_Y': {'file': 'wircam-cfht-y.txt', 'description': 'CFHT WIRCAM Y', 'zeropoint': 2073.32}, \
    'WIRCAM_J': {'file': 'wircam-cfht-j.txt', 'description': 'CFHT WIRCAM J', 'zeropoint': 1551.01}, \
    'WIRCAM_H': {'file': 'wircam-cfht-h.txt', 'description': 'CFHT WIRCAM H', 'zeropoint': 1044.35}, \
    'WIRCAM_KS': {'file': 'wircam-cfht-ks.txt', 'description': 'CFHT WIRCAM Ks', 'zeropoint': 674.62}, \
    'WIRCAM_KCONT': {'file': 'wircam-cfht-kcont.txt', 'description': 'CFHT WIRCAM K-cont', 'zeropoint': 636.17}, \
    'WIRCAM_CH4_OFF': {'file': 'wircam-cfht-ch4s.txt', 'description': 'CFHT WIRCAM CH4-on', 'zeropoint': 987.39}, \
    'WIRCAM_CH4_ON': {'file': 'wircam-cfht-ch4l.txt', 'description': 'CFHT WIRCAM CH4-off', 'zeropoint': 1076.31}, \
    'WIRC_J': {'file': 'wirc_jcont.txt', 'description': 'WIRC J-cont', 'zeropoint': 0.}, \
    'WIRC_H': {'file': 'wirc_hcont.txt', 'description': 'WIRC H-cont', 'zeropoint': 0.}, \
    'WIRC_K': {'file': 'wirc_kcont.txt', 'description': 'WIRC K-cont', 'zeropoint': 0.}, \
    'WIRC_CO': {'file': 'wirc_co.txt', 'description': 'WIRC CO', 'zeropoint': 0.}, \
    'WIRC_CH4S': {'file': 'wirc_ch4s.txt', 'description': 'WIRC CH4S', 'zeropoint': 0.}, \
    'WIRC_CH4L': {'file': 'wirc_ch4l.txt', 'description': 'WIRC CH4L', 'zeropoint': 0.}, \
    'WIRC_FE2': {'file': 'wirc_feii.txt', 'description': 'WIRC Fe II', 'zeropoint': 0.}, \
    'WIRC_BRGAMMA': {'file': 'wirc_brgamma.txt', 'description': 'WIRC H I Brackett Gamma', 'zeropoint': 0.}, \
    'WIRC_PABETA': {'file': 'wirc_pabeta.txt', 'description': 'WIRC H I Paschen Beta', 'zeropoint': 0.}, \
    'WISE_W1': {'file': 'wise_w1.txt', 'description': 'WISE W1 (3.5 micron)', 'zeropoint': 309.54}, \
    'WISE_W2': {'file': 'wise_w2.txt', 'description': 'WISE W2 (4.6 micron)', 'zeropoint': 171.79}, \
    'WISE_W3': {'file': 'wise_w3.txt', 'description': 'WISE W3 (13 micron)', 'zeropoint': 31.67}, \
    'WISE_W4': {'file': 'wise_w4.txt', 'description': 'WISE W4 (22 micron)', 'zeropoint': 8.363} \
    }

# Index sets
index_sets = ['burgasser','bardalez','tokunaga','reid','geballe','allers','testi','slesnick','mclean','rojas']

# change the command prompt
sys.ps1 = 'splat> '

#####################################################


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
    '''
    :Description: Primary class for containing spectral and source data for SpeX Prism Library.

    :param model:
    :type model: optional, default = False
    :param wlabel: label of wavelength
    :type wlabel: optional, default = 'Wavelength'
    :param wunit: unit in which wavelength is measured
    :type wunit: optional, default = ``u.micron``
    :param wunit_label: label of the unit of wavelength
    :type wunit_label: optional, default = :math:`\\mu m`
    :param flabel: label of flux density
    :type flabel: optional, default = :math:`F_\\lambda`
    :param fscale: string describing how flux density is scaled
    :type fscale: optional, default = ''
    :param funit: unit in which flux density is measured
    :type funit: optional, default = ``u.erg/(u.cm**2 * u.s * u.micron)``
    :param funit_label: label of the unit of flux density
    :type funit_label: optional, default = :math:`erg\;cm^{-2} s^{-1} \\mu m^{-1}`
    :param resolution:
    :type resolution: optional, default = 150
    :param slitpixelwidth: Width of the slit measured in subpixel values.
    :type slitpixelwidth: optional, default = 3.33
    :param slitwidth: Actual width of the slit, measured in arc seconds. Default value is the ``slitpixelwidth`` multiplied by the spectrograph pixel scale of 0.15 arcseconds.
    :type slitwidth: optional, default = ``slitpixelwidth`` * 0.15
    :param header: header info of the spectrum
    :type header: optional, default = Table()
    :param filename: a string containing the spectrum's filename.
    :type filename: optional, default = ''
    :param file: same as filename
    :type file: optional, default = ''
    :param idkey: spectrum key of the desired spectrum
    :type idkey: optional, default = False

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
        self.model = kwargs.get('model',False)
        self.wlabel = kwargs.get('wlabel',r'Wavelength')
        self.wunit = kwargs.get('wunit',u.micron)
        self.wunit_label = kwargs.get('wunit_label',r'$\mu$m')
        self.flabel = kwargs.get('flabel',r'F$_{\lambda}$')
        self.fscale = kwargs.get('fscale','')
        self.funit = kwargs.get('funit',u.erg/(u.cm**2 * u.s * u.micron))
        self.funit_label = kwargs.get('funit_label',r'erg~cm$^{-2}$~s$^{-1}$~$\mu$m$^{-1}$')
        self.resolution = kwargs.get('resolution',150)    # default placeholder
        self.slitpixelwidth = kwargs.get('slitwidth',3.33)        # default placeholder
        self.slitwidth = self.slitpixelwidth*spex_pixel_scale
#        self.header = kwargs.get('header',fits.PrimaryHDU())
        self.header = kwargs.get('header',{})
        self.filename = kwargs.get('file','')
        self.filename = kwargs.get('filename',self.filename)
        self.idkey = kwargs.get('idkey',False)

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
        elif self.model == False and self.filename != '':
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
            if ~isinstance(self.wave,astropy.units.quantity.Quantity):
                self.wave = numpy.array(self.wave)*self.wunit
            if ~isinstance(self.flux,astropy.units.quantity.Quantity):
                self.flux = numpy.array(self.flux)*self.funit
            if ~isinstance(self.wave,astropy.units.quantity.Quantity):
                self.noise = numpy.array(self.noise)*self.funit
# some conversions
            self.flam = self.flux
            self.nu = self.wave.to('Hz',equivalencies=u.spectral())
            self.fnu = self.flux.to('Jy',equivalencies=u.spectral_density(self.wave))
            self.fnu_unit = u.Jansky
# calculate variance
            self.variance = self.noise**2
            self.dof = numpy.round(len(self.wave)/self.slitpixelwidth)
# signal to noise
            w = numpy.where(self.flux.value > numpy.median(self.flux.value))
            self.snr = numpy.nanmedian(self.flux.value[w]/self.noise.value[w])

# preserve original values
            self.wave_original = copy.deepcopy(self.wave)
            self.flux_original = copy.deepcopy(self.flux)
            self.noise_original = copy.deepcopy(self.noise)
            self.variance_original = copy.deepcopy(self.variance)
#        self.resolution = copy.deepcopy(self.resolution)
#        self.slitpixelwidth = copy.deepcopy(self.slitpixelwidth)
        else:
            print ('Warning: not information provided, creating an empty Spectrum object')

# populate information on source and spectrum from database
        if sdb != False:
            for k in sdb.keys():
                setattr(self,k.lower(),sdb[k][0])
            self.shortname = designationToShortName(self.designation)
            self.date = self.observation_date
# convert some data into numbers
            kconv = ['ra','dec','julian_date','median_snr','resolution','airmass',\
            'jmag','jmag_error','hmag','hmag_error','kmag','kmag_error','source_key']
            for k in kconv:
                try:
                    setattr(self,k,float(getattr(self,k)))
                except:
                    setattr(self,k,numpy.nan)
#                print(getattr(self,k))
        

# information on model
        if self.model == True:
            self.teff = kwargs.get('teff',numpy.nan)
            self.logg = kwargs.get('logg',numpy.nan)
            self.z = kwargs.get('z',numpy.nan)
            self.fsed = kwargs.get('fsed',numpy.nan)
            self.cld = kwargs.get('cld',numpy.nan)
            self.kzz = kwargs.get('kzz',numpy.nan)
            self.slit = kwargs.get('slit',numpy.nan)
            self.modelset = kwargs.get('set','')
            self.name = self.modelset+' Teff='+str(self.teff)+' logg='+str(self.logg)+' [M/H]='+str(self.z)
            self.fscale = 'Surface'

# populate header            
        kconv = ['designation','name','shortname','ra','dec','slitwidth','source_key','data_key','observer',\
            'data_reference','discovery_reference','program_pi','program_number','airmass','reduction_spextool_version',\
            'reduction_person','reduction_date','observation_date','julian_date','median_snr','resolution','airmass']
        for k in kconv:
            try:
                self.header[k] = getattr(self,k)
            except:
                self.header[k] = ''

        self.history = ['Loaded']


    def __copy__(self):
            s = type(self)()
            s.__dict__.update(self.__dict__)
            return s

# backup version
    def copy(self):
            s = type(self)()
            s.__dict__.update(self.__dict__)
            return s

    def __repr__(self):
        '''
        :Purpose: A simple representation of an object is to just give it a name
        '''
        return 'Spectrum of {}'.format(self.name)

    def __add__(self,other):
        '''
        :Purpose: Adds two spectra
        :param other: the other spectrum to add to the object
        '''
        sp = self.copy()
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.add(self.flux,f(self.wave)*other.funit)
        sp.variance = sp.variance+n(self.wave)*(other.funit**2)
        sp.noise = sp.variance**0.5
        sp.flux_original=sp.flux
        sp.noise_original=sp.noise
        sp.variance_original=sp.variance
        sp.name = self.name+' + '+other.name
        return sp

    def __sub__(self,other):
        '''
        :Purpose: Subtracts two spectra
        :param other: the other spectrum to subtract from the object
        '''
        sp = self.copy()
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.subtract(self.flux,f(self.wave)*other.funit)
        sp.variance = sp.variance+n(self.wave)*(other.funit**2)
        sp.noise = sp.variance**0.5
        sp.flux_original=sp.flux
        sp.noise_original=sp.noise
        sp.variance_original=sp.variance
        sp.name = self.name+' - '+other.name
        return sp

    def __mul__(self,other):
        '''
        :Purpose: Multiplies two spectra
        :param other: the other spectrum to multiply by the object
        '''
        sp = self.copy()
        f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
        n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=numpy.nan)
        sp.flux = numpy.multiply(self.flux,f(self.wave)*other.funit)
        sp.variance = numpy.multiply(sp.flux**2,(\
            numpy.divide(self.variance.value,sp.flux.value**2)+numpy.divide(n(self.wave),f(self.wave)**2)))
        sp.noise = sp.variance**0.5
        sp.flux_original=sp.flux
        sp.noise_original=sp.noise
        sp.variance_original=sp.variance
        sp.name = self.name+' x '+other.name
        sp.funit = sp.flux.unit
        return sp

    def __div__(self,other):
        '''
        :Purpose: Divides two spectra
        :param other: the other spectrum to divide by the object
        '''
        sp = self.copy()
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
        sp.flux_original=sp.flux
        sp.noise_original=sp.noise
        sp.variance_original=sp.variance
        sp.name = self.name+' / '+other.name
        sp.funit = sp.flux.unit
        return sp

    def info(self):
        '''
        :Purpose: Reports some information about this spectrum.
        '''
        if (self.model):
            print('\n{} model with the following parmeters:'.format(self.modelset))
            print('Teff = {} K'.format(self.teff))
            print('logg = {} cm/s2'.format(self.logg))
            print('z = {}'.format(self.z))
            print('fsed = {}'.format(self.fsed))
            print('cld = {}'.format(self.cld))
            print('kzz = {}'.format(self.kzz))
            print('Smoothed to slit width {} arcseconds\n'.format(self.slit))
        else:
#            print('\nSpectrum of {0} observed on {1}'''.format(self.name, self.date))
            text = ['Spectrum of','Observed on','Observed by','at an airmass of','Full source designation is', 'Median S/N =','SPLAT source key is','SPLAT spectrum key is']
            ref = ['name','date','observer','airmass','designation','median_snr','source_key','data_key']
            for i,k in enumerate(ref):
                try:
                    if getattr(self,k) != '':
                        print('\t{} {}'.format(text[i],getattr(self,k)))
                except:
                    pass
            try:
                bib = getBibTex(self.data_reference)
                print('\tData published in {}'.format(shortRef(bib)))
            except:
                pass
#        print('\nPlot spectrum using .plot()')
        return

    def flamToFnu(self):
         '''
         :Purpose: Converts flux density from :math:`F_\\lambda` to :math:`F_\\nu`, the later in Jy
         '''
         self.funit = u.Jy
         self.flabel = 'F_nu'
         self.flux.to(self.funit,equivalencies=u.spectral_density(self.wave))
         self.noise.to(self.funit,equivalencies=u.spectral_density(self.wave))
         return

    def fluxCalibrate(self,f,mag,**kwargs):
        '''
        :Purpose: Calibrates spectrum to input magnitude
        '''
        absolute = kwargs.get('absolute',False)
        apparent = kwargs.get('apparent',False)
        self.normalize()
        apmag,apmag_e = filterMag(self,f,**kwargs)
# NOTE: NEED TO INCORPORATE UNCERTAINTY INTO SPECTRAL UNCERTAINTY
        if (~numpy.isnan(apmag)):
            self.scale(10.**(0.4*(apmag-mag)))
            if (absolute):
                self.fscale = 'Absolute'
            if (apparent):
                self.fscale = 'Apparent'
        return

# determine maximum flux, by default in non telluric regions
    def fluxMax(self,**kwargs):
        '''
        :Purpose: Determines maximum flux of spectrum
        :param maskTelluric: masks telluric regions
        :type maskTelluric: optional, default = True
        '''
        if kwargs.get('maskTelluric',True):
            return numpy.nanmax(self.flux.value[numpy.where(\
                numpy.logical_or(\
                    numpy.logical_and(self.wave > 0.9*u.micron,self.wave < 1.35*u.micron),
                    numpy.logical_and(self.wave > 1.42*u.micron,self.wave < 1.8*u.micron),
                    numpy.logical_and(self.wave > 1.92*u.micron,self.wave < 2.3*u.micron)))])*self.funit
        else:
            return numpy.nanmax(self.flux.value[numpy.where(\
                numpy.logical_and(self.wave > 0.9*u.micron,self.wave < 2.3*u.micron))])*self.funit

    def fnuToFlam(self):
         '''
         :Purpose: Converts flux density from :math:`F_\\nu` to :math:`F_\\lambda`, the later in erg/s/cm2/Hz
         '''
         self.funit = u.erg/(u.cm**2 * u.s * u.micron)
         self.flabel = 'F_lam'
         self.flux.to(self.funit,equivalencies=u.spectral_density(self.wave))
         self.noise.to(self.funit,equivalencies=u.spectral_density(self.wave))
         self.variance = self.noise**2
         return

    def normalize(self,**kwargs):
        '''
        :Purpose: Normalizes spectrum to its maximum flux.
        :param maskTelluric: masks telluric regions
        :type maskTelluric: optional, default = True
        '''
        if kwargs.get('waveRange',False) != False:
            rng = kwargs.get('waveRange')
            if not isinstance(rng,list):
                rng = [rng]
            if len(rng) < 2:
                rng = [rng[0]-0.02,rng[0]+0.02]
            self.scale(1./numpy.nanmax(self.flux.value[numpy.where(numpy.logical_and(self.wave > rng[0]*u.micron,self.wave < rng[1]*u.micron))]))
        else:
            self.scale(1./self.fluxMax(**kwargs).value)
        self.fscale = 'Normalized'
        return

    def plot(self,**kwargs):
        '''
        :Purpose: Plots spectrum. See SPLAT Plotting Routines page for more details.
        '''
        plotSpectrum(self,showNoise=True,showZero=True,**kwargs)

    def reset(self):
        '''
        :Purpose: Resets to original spectrum
        '''
        self.wave = copy.deepcopy(self.wave_original)
        self.flux = copy.deepcopy(self.flux_original)
        self.noise = copy.deepcopy(self.noise_original)
        self.variance = copy.deepcopy(self.variance_original)
        self.resolution = copy.deepcopy(self.resolution_original)
        self.slitpixelwidth = copy.deepcopy(self.slitpixelwidth_original)
        self.slitwidth = self.slitpixelwidth*spex_pixel_scale
        self.fscale = ''
        return

    def export(self,*args,**kwargs):
        '''
        :Purpose: export spectrum object to a file, either fits or ascii depending on file extension
        '''
        filename = self.simplefilename
        if len(args) > 0:
            filename = args[0]
        filename = kwargs.get('filename',filename)
        filename = kwargs.get('file',filename)

# determine which type of file
        ftype = filename.split('.')[-1]

# fits file
        if (ftype == 'fit' or ftype == 'fits'):
            try:
                data = numpy.vstack((self.wave.value,self.flux.value,self.noise.value)).T
                hdu = fits.PrimaryHDU(data)
                for k in self.header.keys():
                    hdu.header[k] = self.header[k]
                hdu.writeto(filename,clobber=True)
            except:
                raise NameError('Problem saving spectrum object to file {}'.format(filename))

# ascii file - by default space delimited (could make this more flexible)
        else:
            try:
                t = Table([self.wave.value,self.flux.value,self.noise.value],names=['wavelength','flux','uncertainty'])
                if kwargs.get('header',True):
                    hd = ['{} = {}'.format(k,self.header[k]) for k in self.header.keys()]
                    hd = list(set(hd))
                    hd.sort()
                    t.meta['comments']=hd
                ascii.write(t,output=filename,format=kwargs.get('format','commented_header'))
            except:
                raise NameError('Problem saving spectrum object to file {}'.format(filename))
        return

    def scale(self,factor):
         '''
        :Purpose: Smooths spectrum to a constant slit width (smooths by pixels)
        :param method: the type of window to create. See http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details.
        :type method: optional, default = 'hanning'
        :param slitpixelwidth: Width of the slit measured in subpixel values.
        :type slitpixelwidth: optional, default = 3.33
        :param slitwidth: Actual width of the slit, measured in arc seconds. Default value is the ``slitpixelwidth`` multiplied by the spectrograph pixel scale of 0.15 arcseconds.
        :param resolution:
        :type resolution: optional, default = 150
        :type slitwidth: optional, default = slitpixelwidth * 0.15
        '''
         self.flux = self.flux*factor
         self.noise = self.noise*factor
         self.variance = self.noise**2
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
        '''
        :Purpose: Smooths spectrum to a constant resolution
        :param resolution: resolution for smoothing the spectrum
        :param overscale: used for computing number of samples in the window
        :type overscale: optional, default = 10.
        :param method: the type of window to create. See http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details.
        :type method: optional, default = 'hanning'
        '''
        overscale = kwargs.get('overscale',10.)
        method = kwargs.get('method','hamming')
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
             self.flux = f(self.wave)*self.funit
             self.variance = v(self.wave)*self.funit**2
             self.noise = numpy.array([n**0.5 for n in self.variance])
             self.slitpixelwidth = self.slitpixelwidth*self.resolution/resolution
             self.resolution = resolution
             self.slitwidth = self.slitpixelwidth*spex_pixel_scale
             self.history = ['Smoothed to constant resolution {}'.format(self.resolution)]
        return

    def smoothToSlitPixelWidth(self,width,**kwargs):
        '''
        :Purpose: Smooths spectrum to a constant slit width
        :param width: width for smoothing the spectrum
        :param method: the type of window to create. See http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details.
        :type method: optional, default = 'hanning'
        '''
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
        '''
        :Purpose: Smooths spectrum to a constant slit width (smooths by pixels)
        :param width: width for smoothing the spectrum
        :param method: the type of window to create. See http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.get_window.html for more details.
        :type method: optional, default = 'hanning'
        '''
        method = kwargs.get('method','hanning')
        kwargs['method'] = method
        pwidth = width/spex_pixel_scale
        self.smoothToSlitPixelWidth(pwidth,**kwargs)
        return

    def snr(self):
        '''
        :Purpose: Compute a representative S/N value
        .. note:: Unfinished
        '''
        pass
        return

    def surface(self,radius):
         '''
         :Purpose: Convert to surface fluxes given a radius, assuming at absolute fluxes
         .. note:: Unfinished
         '''
         pass
         return

    def waveRange(self):
        ii = numpy.where(self.flux.value > 0)
        return [numpy.nanmin(self.wave[ii]), numpy.nanmax(self.wave[ii])]




# FUNCTIONS FOR SPLAT
def caldateToDate(d):
    '''
    :Purpose: Convert from numeric date to calendar date, and vice-versa.
    :param d: A numeric date of the format '20050412', or a date in the
                calendar format '2005 Jun 12'
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
    :Purpose: Checks if a spectrum file exists in the SPLAT's library.
    :param filename: A string containing the spectrum's filename.
    :Example:
       >>> import splat
       >>> spectrum1 = 'spex_prism_1315+2334_110404.fits'
       >>> print splat.checkFile(spectrum1)
       True
       >>> spectrum2 = 'fake_name.fits'
       >>> print splat.checkFile(spectrum2)
       False
    '''
    url = kwargs.get('url',SPLAT_URL)+DATA_FOLDER
    return requests.get(url+filename).status_code == requests.codes.ok
#    flag = checkOnline()
#    if (flag):
#        try:
#            r = requests.get(url+filename)
#            open(os.path.basename(filename), 'wb').write(r.content)
#            open(os.path.basename(filename), 'wb').write(urllib2.urlopen(url+filename).read())
#        except:
#            flag = False
#    return flag


def checkAccess(**kwargs):
    '''
    :Purpose: Checks if user has access to unpublished spectra in SPLAT library.
    :Example:
       >>> import splat
       >>> print splat.checkAccess()
       True
    :Note: Must have the file .splat_access in your home directory with the correct passcode to use.
    '''
    access_file = '.splat_access'
    result = False

    try:
        home = os.environ.get('HOME')
        if home == None:
            home = './'
#        print(home)
#        print(SPLAT_URL+access_file)
        bcode = requests.get(SPLAT_URL+access_file).content
#        print(bcode)
        lcode = base64.b64encode(open(home+'/'+access_file,'r').read())
#        print(lcode)
        if (bcode in lcode):        # changed to partial because of EOL variations
            result = True
    except:
        result = False

    if (kwargs.get('report','') != ''):
        if result == True:
            print('You have full access to all SPLAT data')
        else:
            print('You have access only to published data')
    return result


def checkLocal(inputfile):
    '''
    :Purpose: Checks if a file is present locally or within the SPLAT
                code directory
    :Example:
       >>> import splat
       >>> splat.checkLocal('splat.py')
       True  # found the code
       >>> splat.checkLocal('parameters.txt')
       False  # can't find this file
       >>> splat.checkLocal('SpectralModels/BTSettl08/parameters.txt')
       True  # found it
    '''
    if not os.path.exists(inputfile):
        if not os.path.exists(SPLAT_PATH+inputfile):
            return ''
        else:
            return SPLAT_PATH+inputfile
    else:
        return inputfile


def checkOnline(*args):
    '''
    :Purpose: Checks if SPLAT's URL is accessible from your machine--
                that is, checks if you and the host are online. Alternately
                checks if a given filename is present locally or online
    :Example:
       >>> import splat
       >>> splat.checkOnline()
       True  # SPLAT's URL was detected.
       >>> splat.checkOnline()
       False # SPLAT's URL was not detected.
       >>> splat.checkOnline('SpectralModels/BTSettl08/parameters.txt')
       '' # Could not find this online file.
    '''
    try:
        if (len(args) != 0):
            if 'http://' in args[0]:
                if requests.get(args[0]).status_code == requests.codes.ok:
                	return args[0]
                return ''
            else:
                if requests.get(SPLAT_URL+args[0]).status_code == requests.codes.ok:
                    return SPLAT_URL+args[0]
                return ''
        else:
        	return requests.get(SPLAT_URL).status_code == requests.codes.ok
    except:
        return ''




def classifyByIndex(sp, *args, **kwargs):
    '''
    :Purpose: Determine the spectral type and uncertainty for a spectrum
                based on indices. Makes use of published index-SpT relations
                from `Reid et al. (2001) <http://adsabs.harvard.edu/abs/2001AJ....121.1710R>`_;
                `Testi et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...552L.147T>`_;
                `Allers et al. (2007) <http://adsabs.harvard.edu/abs/2007ApJ...657..511A>`_;
                and `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_. Returns 2-element tuple
                containing spectral type (numeric or string) and
                uncertainty.

    :param sp: Spectrum class object, which should contain wave, flux and
               noise array elements.

    :param set: named set of indices to measure and compute spectral type

        - *'allers'*: H2O from `Allers et al. (2007) <http://adsabs.harvard.edu/abs/2007ApJ...657..511A>`_
        - *'burgasser'*: H2O-J, CH4-J, H2O-H, CH4-H, CH4-K from `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_
        - *'reid'*:H2O-A and H2O-B from `Reid et al. (2001) <http://adsabs.harvard.edu/abs/2001AJ....121.1710R>`_
        - *'testi'*: sHJ, sKJ, sH2O_J, sH2O_H1, sH2O_H2, sH2O_K from `Testi et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...552L.147T>`_

    :type set: optional, default = 'burgasser'
    :param string: return spectral type as a string (uses typeToNum)
    :type string: optional, default = False
    :param round: rounds off to nearest 0.5 subtypes
    :type round: optional, default = False
    :param remeasure: force remeasurement of indices
    :type remeasure: optional, default = True
    :param nsamples: number of Monte Carlo samples for error computation
    :type nsamples: optional, default = 100
    :param nloop: number of testing loops to see if spectral type is within a certain range
    :type nloop: optional, default = 5

    :Example:
    >>> import splat
    >>> spc = splat.getSpectrum(shortname='0559-1404')[0]
    >>> print splat.classifyByIndex(spc, string=True, set='burgasser', round=True)
        ('T4.5', 0.2562934083414341)

    .. note::
        * Need to allow output of individual spectral types from individual indices
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
    allowed_sets = ['aganze','burgasser','reid','testi','allers']
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
        vals = numpy.polyval(coeffs[index]['coeff'],numpy.random.normal(indices[index][0],indices[index][1],nsamples))*sptfact
        coeffs[index]['spt'] = numpy.nanmean(vals)+sptoffset
        coeffs[index]['sptunc'] = (numpy.nanstd(vals)**2+coeffs[index]['fitunc']**2)**0.5
#        print(index, coeffs[index]['spt'], coeffs[index]['range'], coeffs[index]['spt'] < coeffs[index]['range'][0], coeffs[index]['spt'] > coeffs[index]['range'][1])
        if (coeffs[index]['spt'] < coeffs[index]['range'][0] or coeffs[index]['spt'] > coeffs[index]['range'][1]):
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

    return spt, sptn_e



def classifyByStandard(sp, *args, **kwargs):
    '''
    :Purpose: Determine the spectral type and uncertainty for a
                spectrum by direct comparison to spectral standards.
                Standards span M0-T9 and include the standards listed in
                `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_
                with addition of UGPS 0722-0540  as the T9 standard.  Returns the best
                match or an F-test weighted mean and uncertainty. There is an option
                to follow the procedure of `Kirkpatrick et al. (2010)
                <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_, fitting only in
                the 0.9-1.4 micron region.

    .. :Usage: spt,unc = splat.classifyByStandard(sp, \**kwargs)

    :param sp: spectrum class object, which should contain wave, flux and
               noise array elements.
    :param best: return the best fit standard type only
    :type best: optional, default = True
    :param average: return an chi-square weighted type only
    :type average: optional, default = True
    :param compareto: compare to a single standard (string or number)
    :type compareto: optional, default = False
    :param plot: generate a plot comparing best fit standard to source, can be saved to a file using the ``file`` keyword
    :type plot: optional, default = False
    :param file: output spectrum plot to a file
    :type file: optional, default = ''
    :param method: set to ``'kirkpatrick'`` to follow the `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_  method, fitting only to the 0.9-1.4 micron band
    :type method: optional, default = ''
    :param sptrange: constraint spectral type range to fit, can be strings or numbers
    :type sptrange: optional, default = ['M0','T9']
    :param string: return spectral type as a string
    :type string: optional, default = True
    :param verbose: give lots of feedback
    :type verbose: optional, default = False

    :Example:
    >>> import splat
    >>> spc = splat.getSpectrum(shortname='1507-1627')[0]
    >>> print splat.classifyByStandard(spc,string=True,method='kirkpatrick',plot=True)
        ('L4.5', 0.7138959194725174)
    '''

    verbose = kwargs.get('verbose',False)
    method = kwargs.get('method','')
    best_flag = kwargs.get('best',True)
    average_flag = kwargs.get('average',not best_flag)
    sptrange = kwargs.get('sptrange',[10,39])
    sptrange = kwargs.get('range',sptrange)
    sptrange = kwargs.get('spt',sptrange)
    if not isinstance(sptrange,list):
        sptrange = [sptrange,sptrange]
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

# read in standards is necessary
    initiateStandards(**kwargs)

# get standards
#    stds = getStandard(sptrange,**kwargs)

# assign subclasses
    if kwargs.get('sd',False):
        stds = splat.SPEX_SD_STDS
        subclass = 'sd'
    elif kwargs.get('esd',False):
        stds = splat.SPEX_ESD_STDS
        subclass = 'esd'
    else:
        stds = splat.SPEX_STDS
        subclass = ''

# select desired spectral range
    spt_allowed = numpy.array([typeToNum(s) for s in stds.keys()])
    spt_sample = spt_allowed[numpy.where(spt_allowed >= sptrange[0])]
    spt_sample = spt_sample[numpy.where(spt_sample <= sptrange[1])]

# determine comparison range based on method
    if (method == 'kirkpatrick'):
        comprng = [0.9,1.4]*u.micron         # as prescribed in Kirkpatrick et al. 2010, ApJS,
    else:
        comprng = [0.7,2.45]*u.micron       # by default, compare whole spectrum

# which comparison statistic to use
    if numpy.isnan(numpy.median(sp.noise)):
        compstat = 'stddev_norm'
    else:
        compstat = 'chisqr'

# compute fitting statistics
    stat = []
    sspt = []

    for t in spt_sample:
        chisq,scale = compareSpectra(sp,stds[typeToNum(t,subclass=subclass)],fit_ranges=[comprng],stat=compstat,novar2=True)
        stat.append(chisq)
        sspt.append(t)
#        if (verbose):
#            print(t, chisq, scale)

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
        spt = typeToNum(sptn,uncertainty=sptn_e,subclass=subclass)
    else:
        spt = sptn

# plot spectrum compared to best spectrum
    if (kwargs.get('plot',False) != False):
#        spstd = Spectrum(file=sorted_stdfiles[0])
#        print(typeToNum(sorted_stdsptnum[0],subclass=subclass))
        spstd = getStandard(typeToNum(sorted_stdsptnum[0],subclass=subclass))[0]
        chisq,scale = compareSpectra(sp,spstd,fit_ranges=[comprng],stat=compstat)
        spstd.scale(scale)
        if kwargs.get('colors',False) == False:
            kwargs['colors'] = ['k','r']
        if kwargs.get('labels',False) == False:
            kwargs['labels'] = [sp.name,typeToNum(sorted_stdsptnum[0],subclass=subclass)+' Standard']
        plotSpectrum(sp,spstd,**kwargs)

    return spt, sptn_e




def classifyByTemplate(sp, *args, **kwargs):
    '''
    :Purpose: Determine the spectral type and uncertainty for a
                spectrum by direct comparison to a large set of spectra in
                the library. One can select down the spectra by using the set
                command. Returns the best match or an F-test weighted mean and
                uncertainty. There is an option to follow  the procedure of
                `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_,
                fitting only in the 0.9-1.4 micron region.

    .. :Usage: result = splat.classifyByTemplate(sp, \*args, \**kwargs)

    :Output: Returns a dictionary containing the following keys:

                    - **result**: a tuple containing the spectral type and its uncertainty)
                    - **chisquare**: array of nbest chi-square values
                    - **name**: array of nbest source names
                    - **scale**: array of nbest optimal scale factors
                    - **spectra**: array of nbest Spectrum objects
                    - **spt**: array of nbest spectral types

    :param sp: Spectrum class object, which should contain wave, flux and
               noise array elements.
    :param spt or spt_range: restrict the spectral type range over which templates are chosen
    :type best: optional, default = None
    :param best: return only the best fit template type
    :type best: optional, default = False
    :param nbest: number of best fitting spectra to return
    :type nbest: optional, default = 1
    :param plot: generate a plot comparing best fit standard to source, can be save to a file using the file keyword
    :type plot: optional, default = False
    :param file: output spectrum plot to a file
    :type file: optional, default = ''
    :param method: set to ``'kirkpatrick'`` to follow the `Kirkpatrick et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..190..100K>`_ method, fitting only to the 0.9-1.4 micron band
    :type method: optional, default = ''
    :param set: string defining which spectral template set you want to compare to; several options which can be combined:

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
            - *standard*: only spectral standards (use classifyByStandard instead)

    :type set: optional, default = ''

    :param string: return spectral type as a string
    :type string: optional, default = True
    :param spt_type: specify which spectral classification type to return; can be 'spex', 'opt', 'nir', or 'lit'
    :type spt_type: optional, default = 'literature'
    :param verbose: give lots of feedback
    :type verbose: optional, default = False

    :Example:
    >>> import splat
    >>> spc = splat.getSpectrum(shortname='1507-1627')[0]
    >>> print splat.classifyByTemplate(spc,string=True,set='l dwarf, high sn', spt_type='spex', plot=True)
        ('L4.5', 0.7138959194725174)
    '''

#
    spt_type = kwargs.get('spt_type','literature')
    spt_range = kwargs.get('spt_range',[10.,39.9])
    spt_range = kwargs.get('spt',spt_range)
    nbest = kwargs.get('nbest',1)
    verbose = kwargs.get('verbose',False)
    published = kwargs.get('published','')
    published = kwargs.get('public',published)
    set = kwargs.get('select','')
#   placeholder for a systematic uncertainty term
    unc_sys = 0.
    if (kwargs.get('method','') == 'kirkpatrick'):
        comprng = [0.9,1.4]*u.micron         # as prescribed in Kirkpatrick et al. 2010, ApJS,
    else:
        comprng = [0.7,2.45]*u.micron       # by default, compare whole spectrum

#  canned searches
#  constrain spectral types
    if ('lit' in spt_type.lower()):
        spt_type = 'LIT_TYPE'
    elif ('opt' in spt_type.lower() or 'optical' in set):
        spt_type = 'OPT_TYPE'
    elif ('nir' in spt_type.lower()):
        spt_type = 'NIR_TYPE'
    else:
        spt_type = 'LIT_TYPE'

    if ('m dwarf' in set.lower()):
        spt_range = [10,19.9]
    if ('l dwarf' in set.lower()):
        spt_range = [20,29.9]
    if ('t dwarf' in set.lower()):
        spt_range = [30,39.9]
    if ('vlm' in set.lower()):
        spt_range = [numpy.max([17,spt_range[0]]),numpy.min([39.9,spt_range[-1]])]

#  constrain S/N
    snr = 0.
    if ('high sn' in set.lower()):
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
    if 'giant' in set.lower():
        giant = True
    if 'not giant' in set.lower():
        giant = False
    companion = ''
    if 'companion' in set.lower():
        companion = True
    if 'not companion' in set.lower():
        companion = False
    young = ''
    if 'young' in set.lower():
        young = True
    if 'not young' in set.lower():
        young = False
    binary = ''
    if 'binary' in set.lower():
        binary = True
    if 'not binary' in set.lower():
        binary = False
    spbinary = ''
    if 'spectral binary' in set.lower():
        spbinary = True
    if 'not spectral binary' in set.lower():
        spbinary = False

# REARRANGE THIS - SEND IN KWARGS WITH OUTPUT, LOGIC SET, AND THE REST ARE UP TO USER?

    lib = searchLibrary(excludefile=excludefile,excludekey=excludekey,excludeshortname=excludeshortname, \
        snr=snr,spt_type=spt_type,spt_range=spt_range,published=published, \
        giant=giant,companion=companion,young=young,binary=binary,spbinary=spbinary,output='all',logic='and')

#    print([x for x in compsp])
#    elif ('companion' in set):
#        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,companion=True,published=published,giant=False,logic='and')
#    elif ('young' in set):
#        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,young=True,published=published,logic='and')
#    elif ('subdwarf' in set):
#        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,subdwarf=True,published=published,logic='and')
#    elif ('single' in set):
#        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,spbinary=False,published=published,binary=False,giant=False,logic='and')
#    elif ('spectral binaries' in set):
#        lib = searchLibrary(output='all',excludefile=excludefile,snr=snr,spt=spt,spbinary=True,published=published,logic='and')
#    elif (set != ''):
#        lib = searchLibrary(output='all',excludefile=excludefile,published=published,snr=snr,spt=spt)
#    else:
#        lib = searchLibrary(output='all',excludefile=excludefile,published=published,**kwargs)

#    lib = lib[:][numpy.where(lib[sptref] != '')]

# first search for the spectra desired - parameters are set by user
    if len(lib) == 0:
        print('\nNo templates available for comparison\n\n')
        return numpy.nan, numpy.nan

    files = lib['DATA_FILE']
    dkey = lib['DATA_KEY']
    sspt = [typeToNum(s) for s in lib[spt_type]]

    if (verbose):
        print('\nComparing to {} templates\n\n'.format(len(files)))
    if len(files) > 100:
        print('\nComparing to {} templates\n\n'.format(len(files)))
        print('This may take some time!\n\n'.format(len(files)))

# do comparison
    stat = []
    scl = []
    for i,d in enumerate(dkey):

# INSERT TRY STATEMNT HERE?

        s = Spectrum(idkey=d)
        chisq,scale = compareSpectra(sp,s,fit_ranges=[comprng],stat='chisqr',novar2=True,*kwargs)
        stat.append(chisq)
        scl.append(scale)
        if (verbose):
            print(keySpectrum(d)['NAME'][0], typeToNum(sspt[i]), chisq, scale)

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
        plotSpectrum(sp,s,colors=['k','r'],title=sp.name+' vs '+s.name,**kwargs)

# string or not?
    if (kwargs.get('string', True) == True):
        spt = typeToNum(sptn,uncertainty=sptn_e)
    else:
        spt = sptn

# return dictionary of results
    return {'result': (spt,sptn_e), \
        'chisquare': sorted(stat)[0:nbest], 'spt': sorted_spt[0:nbest], 'scale': sorted_scale[0:nbest], \
        'name': [keySpectrum(d)['NAME'][0] for d in sorted_dkey[0:nbest]], \
        'spectra': [Spectrum(idkey=d) for d in sorted_dkey[0:nbest]]}



def classifyGravity(sp, *args, **kwargs):
    '''
    :Purpose: Determine the gravity classification of a brown dwarf using the method of `Allers & Liu (2013) <http://adsabs.harvard.edu/abs/2013ApJ...772...79A>`_

    :param sp: Spectrum class object, which should contain wave, flux and
               noise array elements. Must be between M6.0 and L7.0.
    :param indices: specify indices set using ``measureIndexSet``.
    :type indices: optional, default = False
    :param spt: spectral type of ``sp``. Must be between M6.0 and L7.0
    :type spt: optional, default = False
    :param plot: plot against spectral standard
    :type plot: optional, default = False
    :param allscores: returns the full result, including the gravity scores from different indices
    :type allscores: optional, default = False
    :param verbose: give lots of feedback
    :type verbose: optional, default = False

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
    >>> print splat.classifyGravity(sp)
        FLD-G
    >>> print splat.classifyGravity(sp, allscores = True)
        {'VO-z': 0.0, 'FeH-z': 1.0, 'gravity_class': 'FLD-G', 'H-cont': 0.0, 'KI-J': 1.0, 'score': 0.5}
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
        gravscore['score'] = scipy.stats.nanmean(medgrav)
    else:
        gravscore['score'] = scipy.stats.nanmedian(medgrav)

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
    :Purpose: Compare two spectra against each other. Returns the stat value as well as
                the scale factor. Minimum possible value is 1.e-9.

    :param sp1: Spectrum class object
    :param sp2: Spectrum class object to compare with ``sp1``
    :param weights: set weights for flux values of ``sp1``
    :type weights: optional, default = [1, ..., 1] for len(sp1.wave)
    :param mask: mask any flux value of ``sp1``; has to be an array with length equal as ``sp1`` with only 0 (unmask) or 1 (mask).
    :type mask: optional, default = [0, ..., 0] for len(sp1.wave)
    :param fit_ranges: wavelength range, measured in microns
    :type fit_ranges: optional, default = [0.65,2.45]
    :param mask_ranges: mask any flux value of ``sp1`` by specifying the wavelength range. Must be in microns
    :type mask_ranges: optional, default = []
    :param mask_telluric: masks certain wavelengths to avoid effects from telluric absorption
    :type mask_telluric: optional, default = False
    :param mask_standard: same as ``mask_telluric``
    :type mask_standard: optional, default = False
    :param novar2: compute without using variance of ``sp2``
    :type novar2: optional, default = False
    :param stat: string defining which statistical method to use; available options are:

            - *'chisqr'*: compare by computing chi squared value
            - *'stddev'*: compare by computing standard deviation
            - *'stddev_norm'*: compare by computing normalized standard deviation
            - *'absdev'*: compare by computing absolute deviation

    :type stat: optional, default = 'chisqr'
    :param plot: plot ``sp1`` with scaled ``sp2``
    :type plot: optional, default = False

    :Example:
    >>> import splat
    >>> import numpy
    >>> sp1 = splat.getSpectrum(shortname = '2346-3153')[0]
    >>> sp2 = splat.getSpectrum(shortname = '1421+1827')[0]
    >>> print splat.compareSpectra(sp1, sp2)
        (<Quantity 17633.117108524682>, 1.1354366238893992e-15)
    >>> m = numpy.zeros(len(sp1.wave))
    >>> m[200:300] = 1
    >>> m_range = [[1.3, 1.5]]
    >>> splat.compareSpectra(sp1, sp2, mask = m, mask_ranges = m_range, novar2 = True)
        (<Quantity 18980.062433571456>, 1.1694685419610533e-15)
    '''
    weights = kwargs.get('weights',numpy.zeros(len(sp1.wave)))
    mask = kwargs.get('mask',numpy.zeros(len(sp1.wave)))    # mask = 1 -> ignore
    fit_ranges = kwargs.get('fit_ranges',[spex_wave_range])
    mask_ranges = kwargs.get('mask_ranges',[])
    mask_standard = kwargs.get('mask_standard',False)
    mask_telluric = kwargs.get('mask_telluric',mask_standard)
    var_flag = kwargs.get('novar2',True)
    stat = kwargs.get('stat','chisqr')
    minreturn = 1.e-9
    if ~isinstance(fit_ranges[0],astropy.units.quantity.Quantity):
        fit_ranges*=u.micron

    if (mask_standard):
        mask_telluric == True

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
# telluric absorption
    if (mask_telluric):
        mask_ranges.append([0.,0.65])        # meant to clear out short wavelengths
        mask_ranges.append([1.35,1.42])
        mask_ranges.append([1.8,1.92])
        mask_ranges.append([2.45,99.])        # meant to clear out long wavelengths

    if (mask_standard):
        mask_ranges.append([0.,0.8])        # standard short cut
        mask_ranges.append([2.35,99.])        # standard long cut

    if ~isinstance(fit_ranges[0],astropy.units.quantity.Quantity):
        mask_ranges*=u.micron

    for ranges in mask_ranges:
        mask[numpy.where(((sp1.wave >= ranges[0]) & (sp1.wave <= ranges[1])))] = 1

# set the weights
    for ranges in fit_ranges:
        weights[numpy.where(((sp1.wave >= ranges[0]) & (sp1.wave <= ranges[1])))] = 1

# mask flux < 0
    mask[numpy.where(numpy.logical_or(sp1.flux < 0,f(sp1.wave) < 0))] = 1

    weights = weights*(1.-mask)

# comparison statistics

# switch to standard deviation if no uncertainty
    if numpy.isnan(numpy.nanmax(vtot)):
        stat = 'stddev'

# chi^2
    if (stat == 'chisqr'):
# compute scale factor
        scale = numpy.nansum(weights*sp1.flux.value*f(sp1.wave)/vtot)/ \
            numpy.nansum(weights*f(sp1.wave)*f(sp1.wave)/vtot)
# correct variance
        vtot = numpy.nanmax([sp1.variance.value,sp1.variance.value+v(sp1.wave)*scale**2],axis=0)
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave)*scale)**2/vtot)
        unit = sp1.funit/sp1.funit

# normalized standard deviation
    elif (stat == 'stddev_norm'):
# compute scale factor
        scale = numpy.nansum(weights*sp1.flux.value)/ \
            numpy.nansum(weights*f(sp1.wave))
# correct variance
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave)*scale)**2)/ \
            numpy.median(sp1.flux.value)**2
        unit = sp1.funit/sp1.funit

# standard deviation
    elif (stat == 'stddev'):
# compute scale factor
        scale = numpy.nansum(weights*sp1.flux.value*f(sp1.wave))/ \
            numpy.nansum(weights*f(sp1.wave)*f(sp1.wave))
# correct variance
        stat = numpy.nansum(weights*(sp1.flux.value-f(sp1.wave)*scale)**2)
        unit = sp1.funit**2

# absolute deviation
    elif (stat == 'absdev'):
# compute scale factor
        scale = numpy.nansum(weights*sp1.flux.value)/ \
            numpy.nansum(weights*f(sp1.wave))
# correct variance
        stat = numpy.nansum(weights*abs(sp1.flux.value-f(sp1.wave)*scale))
        unit = sp1.funit

# error
    else:
        print('Error: statistic {} for compareSpectra not available'.format(stat))
        return numpy.nan, numpy.nan

# plot spectrum compared to best spectrum
    if (kwargs.get('plot',False) != False):
        spcomp = sp2.copy()
        spcomp.scale(scale)
        kwargs['colors'] = kwargs.get('colors',['k','r','b'])
        kwargs['title'] = kwargs.get('title',sp1.name+' vs '+sp2.name)
        plotSpectrum(sp1,spcomp,sp1-spcomp,**kwargs)

    return numpy.nanmax([stat,minreturn])*unit, scale


def coordinateToDesignation(c):
    '''
    :Purpose: Converts right ascension and declination into a designation string
    :param c: RA and Dec to be converted; can be a SkyCoord object with units of degrees,
              a list with RA and Dec in degrees, or a string with RA measured in hour
              angles and Dec in degrees
    :Example:
    >>> import splat
    >>> from astropy.coordinates import SkyCoord
    >>> c = SkyCoord(238.86, 9.90, unit="deg")
    >>> print splat.coordinateToDesignation(c)
        J15552640+0954000
    >>> print splat.coordinateToDesignation([238.86, 9.90])
        J15552640+0954000
    >>> print splat.coordinateToDesignation('15:55:26.4 +09:54:00.0')
        J15552640+0954000
    '''
# input is ICRS
    if isinstance(c,SkyCoord):
        cc = c
    else:
        cc = properCoordinates(c)
# input is [RA,Dec] pair in degrees
    if sys.version_info.major == 2:
        return string.replace('J{0}{1}'.format(cc.ra.to_string(unit=u.hour, sep='', precision=2, pad=True), \
        cc.dec.to_string(unit=u.degree, sep='', precision=1, alwayssign=True, pad=True)),'.','')
    else:
        return replace('J{0}{1}'.format(cc.ra.to_string(unit=u.hour, sep='', precision=2, pad=True), \
        cc.dec.to_string(unit=u.degree, sep='', precision=1, alwayssign=True, pad=True)),'.','')


def dateToCaldate(d):
    '''
    :Purpose: Converts numeric date to calendar date
    :param d: string in the form 'YYYYMMDD'
    :Example:
    >>> import splat
    >>> print splat.dateToCaldate('19940523')
        1994 May 23
    '''
    return d[:4]+' '+months[int(d[5:6])-1]+' '+d[-2:]


def designationToCoordinate(value, **kwargs):
    '''
    :Purpose: Convert a designation into a RA, Dec tuple or ICRS
    :param value: string with RA measured in hour angles and Dec in degrees
    :param icrs: returns coordinate in ICRS frame if ``True``
    :type icrs: optional, defualt = True
    :Example:
    >>> import splat
    >>> print splat.designationToCoordinate('1555.2640+0954.000')
        <SkyCoord (ICRS): ra=238.86 deg, dec=9.9 deg>
    '''
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
    '''
    :Purpose: Produce a shortened version of designation
    :param value: string with RA measured in hour angles and Dec in degrees
    :Example:
    >>> import splat
    >>> print splat.designationToShortName('1555.2640+0954.000')
        J1555+0954
    '''
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





# DEPRECATED
#def filenameToNameDate(filename):
#    '''Extract from a SPLAT filename the source name and observation date'''
#    ind = filename.rfind('.')
#    base = filename[:ind]
#    spl = base.split('_')
#    if (len(spl) < 2):
#        return '', ''
#    else:
#        name = spl[-2]
#        d = spl[-1]
#        try:
#            float(d)
#            date = '20'+d
#        except ValueError:
#            print(filename+' does not contain a date')
#            date = ''
#
#        return name, date
#


def filterMag(sp,f,*args,**kwargs):
    '''
    :Purpose: Determine the photometric magnitude of a source based on the
                spectrum. Spectra are convolved with the filter specified by
                the ``filter`` input.  By default this filter is also
                convolved with a model of Vega to extract Vega magnitudes,
                but the user can also specify AB magnitudes, photon flux or
                energy flux.

    .. :Usage: mag, unc = splat.filterMag(sp, 'filter name',\*args, \**kwargs) (Unneeded)

    :param sp: Spectrum class object, which should contain wave, flux and
                 noise array elements.
    :param filter: name of filter
    
    :param info: give the filter names available
    :type info: optional, default = False
    :param custom: specify to a 2 x N vector array specifying the wavelengths and transmissions for a custom filter
    :type custom: optional, default = False
    :param ab: compute AB magnitudes
    :type ab: optional, default = False
    :param vega: compute Vega magnitudes
    :type vega: optional, default = True
    :param energy: compute energy flux
    :type energy: optional, default = False
    :param photon: compute photon flux
    :type photon: optional, default = False
    :param filterFolder: folder containing the filter transmission files
    :type filterFolder: optional, default = ''
    :param vegaFile: name of file containing Vega flux file, must be within ``filterFolder``
    :type vegaFile: optional, default = ''
    :param nsamples: number of samples to use in Monte Carlo error estimation
    :type nsamples: optional, default = 100

    :Example:
    >>> import splat
    >>> spc = splat.getSpectrum(shortname='1507-1627')[0]
    >>> spc.fluxCalibrate('2MASS J',14.5)
    >>> print splat.filterMag(spc,'MKO J')
        (14.345894376898123, 0.027596454828421831)
    '''
# keyword parameters
    filterFolder = kwargs.get('filterFolder',SPLAT_PATH+FILTER_FOLDER)
    if not os.path.exists(filterFolder):
        filterFolder = SPLAT_URL+FILTERFOLDER
    vegaFile = kwargs.get('vegaFile','vega_kurucz.txt')
    info = kwargs.get('info',False)
    custom = kwargs.get('custom',False)
    vega = kwargs.get('vega',True)
    ab = kwargs.get('ab',False)
    photons = kwargs.get('photons',False)
    energy = kwargs.get('energy',False)
    nsamples = kwargs.get('nsamples',100)


# check that requested filter is in list
    f = f.replace(' ','_')
    f.upper()
    if (f not in FILTERS.keys() and isinstance(custom,bool)):
        print('Filter '+f+' not included in filterMag')
        info = True

# print out what's available
    if (info):
        print('Filter names:')
        for x in FILTERS.keys():
            print(x+': '+FILTERS[x]['description'])
        return numpy.nan, numpy.nan

# Read in filter
    if isinstance(custom,bool):
        fwave,ftrans = numpy.genfromtxt(filterFolder+FILTERS[f]['file'], comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
    else:
        fwave,ftrans = custom[0],custom[1]
    fwave = fwave[~numpy.isnan(ftrans)]*u.micron   # temporary fix
    ftrans = ftrans[~numpy.isnan(ftrans)]

# interpolate spectrum onto filter wavelength function
    wgood = numpy.where(~numpy.isnan(sp.noise))
    if len(sp.wave[wgood]) > 0:
        d = interp1d(sp.wave[wgood],sp.flux[wgood],bounds_error=False,fill_value=0.)
        n = interp1d(sp.wave[wgood],sp.noise[wgood],bounds_error=False,fill_value=numpy.nan)
# catch for models
    else:
        print(f+': no good points')
        d = interp1d(sp.wave,sp.flux,bounds_error=False,fill_value=0.)
        n = interp1d(sp.wave,sp.flux*1.e-9,bounds_error=False,fill_value=numpy.nan)

    result = []
    if (vega):
# Read in Vega spectrum
        vwave,vflux = numpy.genfromtxt(filterFolder+vegaFile, comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
        vwave = vwave[~numpy.isnan(vflux)]*u.micron
        vflux = vflux[~numpy.isnan(vflux)]*(u.erg/(u.cm**2 * u.s * u.micron))
        vflux.to(sp.funit,equivalencies=u.spectral_density(vwave))
# interpolate Vega onto filter wavelength function
        v = interp1d(vwave,vflux,bounds_error=False,fill_value=0.)
        for i in numpy.arange(nsamples):
#            result.append(-2.5*numpy.log10(trapz(ftrans*numpy.random.normal(d(fwave),n(fwave))*sp.funit,fwave)/trapz(ftrans*v(fwave)*sp.funit,fwave)))
            result.append(-2.5*numpy.log10(trapz(ftrans*(d(fwave)+numpy.random.normal(0,1.)*n(fwave))*sp.funit,fwave)/trapz(ftrans*v(fwave)*sp.funit,fwave)))

    if (energy or photons):
        for i in numpy.arange(nsamples):
#            result.append(trapz(ftrans*numpy.random.normal(d(fwave),n(fwave))*sp.funit,fwave))
            result.append(trapz(ftrans*(d(fwave)+numpy.random.normal(0,1.)*n(fwave))*sp.funit,fwave))
        if (photons):
            convert = const.h.to('erg s')*const.c.to('micron/s')
            result = result/convert.value
    if (ab):
        fnu = fwave.to('Hz',equivalencies=u.spectral())
        fconst = 3631*u.jansky
        fconst = fconst.to(sp.fnu_unit,equivalencies=u.spectral())
        d = interp1d(sp.nu,sp.fnu,bounds_error=False,fill_value=0.)
        n = interp1d(sp.nu,sp.noise,bounds_error=False,fill_value=0.)
        b = trapz((ftrans/fnu)/fconst,fnu)
        for i in numpy.arange(nsamples):
            a = trapz(ftrans*(d(fnu)+numpy.random.normal(0,1)*n(fnu))*sp.fnu_unit/fnu,fnu)
            result.append(-2.5*numpy.log10(a/b))

    val = numpy.nanmean(result)
    err = numpy.nanstd(result)
    if len(sp.wave[wgood]) == 0:
        err = 0.
    return val,err


def filterProperties(f,**kwargs):
    '''
    :Purpose: Returns a dictionary containing key factors about a filter profile.
    .. :Usage: info = splat.filterProperties('filter name',**kwargs)

    :param 'filter name': name of filter, must be one of the specifed filters given by splat.FILTERS
    :param verbose: print out information about filter to screen
    :type custom: optional, default = False
    '''

    filt = f.replace(' ','_')
    filterFolder = kwargs.get('filterFolder',SPLAT_PATH+FILTER_FOLDER)
    if not os.path.exists(filterFolder):
        filterFolder = SPLAT_URL+FILTERFOLDER

    if (filt not in FILTERS.keys()):
        print('Filter '+f+' not among the available filters:')
        for k in FILTERS.keys():
            print('  '+k.replace('_',' ')+': '+FILTERS[k]['description'])
        return None
    else:
        report = {}
        report['name'] = f
        report['description'] = FILTERS[filt]['description']
        report['zeropoint'] = FILTERS[filt]['zeropoint']
        fwave,ftrans = numpy.genfromtxt(filterFolder+FILTERS[filt]['file'], comments='#', unpack=True, \
            missing_values = ('NaN','nan'), filling_values = (numpy.nan))
        report['lambda_mean'] = trapz(fwave*ftrans,fwave)/trapz(ftrans,fwave)
        report['lambda_pivot'] = numpy.sqrt(trapz(fwave*ftrans,fwave)/trapz(ftrans/fwave,fwave))
        fw = fwave[numpy.where(ftrans > 0.5*ftrans)]
        report['lambda_central'] = 0.5*(numpy.max(fw)+numpy.min(fw))
        report['lambda_fwhm'] = numpy.max(fw)-numpy.min(fw)
        fw = fwave[numpy.where(ftrans > 0.01*ftrans)]
        report['lambda_min'] = numpy.min(fw)
        report['lambda_max'] = numpy.max(fw)
        if kwargs.get('verbose',False):
            print('\nFilter '+f+': '+report['description'])
            print('Filter zeropoint = {} Jy'.format(report['zeropoint']))
#            print('Filter mid-points: mean = {:.3f} micron, pivot = {:.3f} micron, central = {:.3f} micron'.format(report['lambda_mean'],report['lambda_pivot'],report['lambda_central'],))
            print('Filter mid-point: = {:.3f} micron'.format(report['lambda_mean'],report['lambda_pivot'],report['lambda_central'],))
            print('Filter FWHM = {:.3f} micron'.format(report['lambda_fwhm']))
            print('Filter range = {:.3f} to {:.3f} micron\n'.format(report['lambda_min'],report['lambda_max']))
        return report


def getSpectrum(*args, **kwargs):
    '''
    :Purpose: Gets a spectrum from the SPLAT library using various selection criteria.
    .. :Usage: [sp] = splat.getSpectrum({search commands},**kwargs)

    :param sp: array of Spectrum class objects, each of which should contain wave, flux and
                 noise array elements.
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
    :param sort: sort results based on Right Ascension
    :type sort: optional, default = True
    :param list: if True, return just a list of the data files (can be done with searchLibrary as well)
    :type list: optional, default = False
    :param lucky: if True, return one randomly selected spectrum from the selected sample
    :type lucky: optional, default = False

    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1507-1627')[0]
        Retrieving 1 file
    >>> sparr = splat.getSpectrum(spt='M7')
        Retrieving 120 files
    >>> sparr = splat.getSpectrum(spt='T5',young=True)
        No files match search criteria
    '''

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
                ind = random.choice(range(len(files)))
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
    NOTE: THIS PROGRAM MAY BE OBSOLETE

    .. :Usage: [sp] = splat.getStandard(spt,**kwargs)

    :param sp: array of Spectrum class objects
    :param spt: spectral type of standard desired ('M7')

    :param optional sd: get subdwarf standard
    :type sd: optional, default = False
    :param optional esd: get extreme subdwarf standard
    :type esd: optional, default = False

    :Example:
    >>> import splat
    >>> sp = splat.getStandard('M7')[0]
        Spectrum of VB 8
    >>> sparr = splat.getStandard('T5',esd=True)
        Type esdT5.0 is not in esd standards: try one of the following:
        ['esdM5.0', 'esdM7.0', 'esdM8.5']
    '''

#    if not isinstance(sptrange,list):
#        sptrange = [sptrange,sptrange]

# make sure standards are read in
    initiateStandards(**kwargs)

# get standards
    if kwargs.get('sd',False):
        stds = splat.SPEX_SD_STDS
        subclass = 'sd'
    elif kwargs.get('esd',False):
        stds = splat.SPEX_ESD_STDS
        subclass = 'esd'
    else:
        stds = splat.SPEX_STDS
        subclass = ''

# set up subtype to use, convert to number then back to string
    if (isinstance(spt,str) != False):
        spt = typeToNum(spt)
    spt = typeToNum(spt,subclass=subclass)

# select among defined spectra - just getting closest match
#    spt_allowed = numpy.array([typeToNum(s) for s in stds.keys()])
#    spt_sample = spt_allowed[numpy.where(spt_allowed >= sptrange[0])]
#    spt_sample = spt_sample[numpy.where(spt_sample <= sptrange[1])]

# nothing there, return
    if spt not in stds.keys():
        print('Type {} is not in {} standards: try one of the following:'.format(spt,subclass))
        print(sorted(stds.keys()))
        return []
    else:
        return [stds[spt]]

# build up file or Spectrum list
#    result = {}
#    for t in spt_sample:
#        result[typeToNum(t,subclass=subclass)] = stds[typeToNum(t,subclass=subclass)]
#
#    return result


def initiateStandards(**kwargs):
    '''
    :Purpose: Initiates the spectral standards (dwarf, subdwarf or extreme subdwarf) in the SpeX library.
    .. :Usage: splat.initateStandards(**kwargs)

    :param optional sd: get subdwarf standard
    :type sd: optional, default = False
    :param optional esd: get extreme subdwarf standard
    :type esd: optional, default = False

    :Example:
    >>> import splat
    >>> splat.initiateStandards()
    >>> splat.SPEX_SDS['M5']
    '''

# choose what kind of standards desired - d, sd, esd
# and read in standards into dictionary if they haven't been already
    if kwargs.get('sd',False):
        if len(splat.SPEX_SD_STDS) == 0:
            for t in splat.SPEX_SD_STDFILES.keys():
                splat.SPEX_SD_STDS[t] = Spectrum(file=splat.SPEX_SD_STDFILES[t])
    elif kwargs.get('esd',False):
        if len(splat.SPEX_ESD_STDS) == 0:
            for t in splat.SPEX_ESD_STDFILES.keys():
                splat.SPEX_ESD_STDS[t] = Spectrum(file=splat.SPEX_ESD_STDFILES[t])
    else:
        if len(splat.SPEX_STDS) == 0:
            for t in splat.SPEX_STDFILES.keys():
                splat.SPEX_STDS[t] = Spectrum(file=splat.SPEX_STDFILES[t])

    return


# simple number checker
def isNumber(s):
    '''
    :Purpose: Checks if something is a number.
    :param s: object to be checked
    :Example:
    >>> import splat
    >>> print splat.isNumber(3)
        True
    >>> print splat.isNumber('hello')
        False
    '''
    try:
        t = float(s)
        return (True and ~numpy.isnan(t))
    except ValueError:
        return False



def estimateDistance(sp, **kwargs):
    '''
    :Purpose: Takes the apparent magnitude and either takes or determines the absolute
                magnitude, then uses the magnitude/distance relation to estimate the
                distance to the object in parsecs. Returns estimated distance and
                uncertainty in parsecs

    :param sp: Spectrum class object, which should be flux calibrated to its empirical apparent magnitude
    :param mag: apparent magnitude of ``sp``
    :type mag: optional, default = False
    :param mag_unc: uncertainty of the apparent magnitude
    :type mag_unc: optional, default = 0
    :param absmag: absolute magnitude of ``sp``
    :type absmag: optional, default = False
    :param absmag_unc: uncertainty of the absolute magnitude
    :type absmag_unc: optional, default = 0
    :param spt: spectral type of ``sp``
    :type spt: optional, default = False
    :param spt_e: uncertainty of the spectral type
    :type spt_e: optional, default = 0
    :param nsamples: number of samples to use in Monte Carlo error estimation
    :type nsamples: optional, default = 100
    :param filter: Name of filter, must be one of the following:

                    - '2MASS J', '2MASS H', '2MASS Ks'
                    - 'MKO J', 'MKO H', 'MKO K', MKO Kp', 'MKO Ks'
                    - 'NICMOS F090M', 'NICMOS F095N', 'NICMOS F097N', 'NICMOS F108N'
                    - 'NICMOS F110M', 'NICMOS F110W', 'NICMOS F113N', 'NICMOS F140W'
                    - 'NICMOS F145M', 'NICMOS F160W', 'NICMOS F164N', 'NICMOS F165M'
                    - 'NICMOS F166N', 'NICMOS F170M', 'NICMOS F187N', 'NICMOS F190N'
                    - 'NIRC2 J', 'NIRC2 H', 'NIRC2 Kp', 'NIRC2 Ks'
                    - 'WIRC J', 'WIRC H', 'WIRC K', 'WIRC CH4S', 'WIRC CH4L'
                    - 'WIRC CO', 'WIRC PaBeta', 'WIRC BrGamma', 'WIRC Fe2'
                    - 'WISE W1', 'WISE W2'

    :type filter: optional, default = False
    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1555+0954')[0]
    >>> print splat.estimateDistance(sp)
        Please specify the filter used to determine the apparent magnitude
        (nan, nan)
    >>> print splat.estimateDistance(sp, mag = 12.521, mag_unc = 0.022, absmag = 7.24, absmag_unc = 0.50, spt = 'M3')
        (116.36999172188771, 33.124820555524224)
    '''

    mag = kwargs.get('mag', False)
    mag_unc = kwargs.get('mag_unc', 0.)
    absmag = kwargs.get('absmag', False)
    absmag_unc = kwargs.get('absmag_unc', 0.)
    spt = kwargs.get('spt', False)
    spt_unc = kwargs.get('spt_e', 0.)
    nsamples = kwargs.get('nsamples', 100)
    filt = kwargs.get('filter', False)

# if no apparent magnitude then calculate from spectrum
    if (mag == False):
        if (filt == False):
            sys.stderr.write('\nPlease specify the filter used to determine the apparent magnitude\n')
            return numpy.nan, numpy.nan
        mag, mag_unc = filterMag(sp,filt)

# if no spt then calculate from spectrum
    if spt == False:
        spt, spt_unc = classifyByIndex(sp)


# if no absolute magnitude then estimate from spectral type
    if absmag == False:
        if filt == False:
            sys.stderr.write('\nPlease specify the filter used to determine the absolute magnitude\n')
            return numpy.nan, numpy.nan
        absmag, absmag_unc = typeToMag(spt,filt,unc=spt_unc)
        print(absmag, absmag_unc)

# create Monte Carlo sets
    if mag_unc > 0.:
        mags = numpy.random.normal(mag, mag_unc, nsamples)
    else:
        mags = nsamples*[mag]

    if absmag_unc > 0.:
        absmags = numpy.random.normal(absmag, absmag_unc, nsamples)
    else:
        absmags = nsamples*[absmag]

# calculate
    distances = 10.**(numpy.subtract(mags,absmags)/5. + 1.)
    d = numpy.mean(distances)
    unc = numpy.std(distances)

    return d, unc


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

    if noiseFlag == True:
#linear fit to continuum
        pCont = numpy.poly1d(numpy.polyfit(wave_cont,fCont,1))
        fContFit = pCont(wLine)
        ew = trapz((numpy.ones(len(fLine))-(fLine/fContFit)), wLine)*1e4
        return ew*u.angstrom, numpy.nan
#monte carlo
    else:
        ew=[]
        for i in range(nsamples):
#generate simulated fluxes
#            fContVar = fCont+numpy.random.normal(0.,1.)*nCont
#            fLineVar = fLine+numpy.random.normal(0.,1.)*nLine
            fContVar = numpy.random.normal(fCont,nCont)
            fLineVar = numpy.random.normal(fLine,nLine)

#linear fit to continuum
            pCont = numpy.poly1d(numpy.polyfit(wave_cont,fContVar,1))
            fContFit = pCont(wLine)
            ew.append(trapz((numpy.ones(len(fLineVar))-(fLineVar/fContFit)), wLine)*1e4)

# some error checking
#            plt.plot(wLine,fContFit,color='r')
#            plt.plot(wLine,fLine,color='k')
#            plt.show()

# following line is more correct but having problem with output
#       return numpy.nanmean(ew)*u.angstrom, numpy.nanstd(ew)*u.angstrom
        return numpy.nanmean(ew), numpy.nanstd(ew)


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
    nsamples = kwargs.get('nsamples',100)
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
    values = numpy.zeros((len(args),nsamples))

# loop over all sampling regions
    for i,waveRng in enumerate(args):
        xNum = (numpy.arange(0,nsamples+1.0)/nsamples)* \
            (numpy.nanmax(waveRng)-numpy.nanmin(waveRng))+numpy.nanmin(waveRng)
        yNum = f(xNum)
        yNum_e = s(xNum)

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
        vals = 0.5*(values[0,:]+values[1,:])/values[2,:]
    elif (method == 'inverse_line'):
        vals = 2.*values[0,:]/(values[1,:]+values[2,:])
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
        return numpy.nanmean(vals), numpy.nan
    else:
        return numpy.nanmean(vals), numpy.nanstd(vals)


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
        refcode = 'ALL13'
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
        refcode = 'BAR14'
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
        refcode = 'BUR06'
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
        refcode = 'GEB02'
        names = ['H2O-1.2','H2O-1.5','CH4-2.2']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.26,1.29],[1.13,1.16],method='ratio',sample='integrate',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[1.57,1.59],[1.46,1.48],method='ratio',sample='integrate',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[2.08,2.12],[2.215,2.255],method='ratio',sample='integrate',**kwargs)
    elif ('mclean' in set.lower()):
        reference = 'Indices from McLean et al. (2003)'
        refcode = 'MCL03'
        names = ['H2OD']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.951,1.977],[2.062,2.088],method='ratio',sample='average',**kwargs)
    elif ('reid' in set.lower()):
        reference = 'Indices from Reid et al. (2001)'
        refcode = 'REI01'
        names = ['H2O-A','H2O-B']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.33,1.35],[1.28,1.30],method='ratio',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[1.47,1.49],[1.59,1.61],method='ratio',sample='average',**kwargs)
    elif ('rojas' in set.lower()):
        reference = 'Indices from Rojas-Ayala et al.(2012)'
        refcode = 'ROJ12'
        names = ['H2O-K2']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        num, er1= measureIndex(sp,[2.070,2.090],[2.235,2.255],method='ratio',sample='average',**kwargs)
        den, er2= measureIndex(sp,[2.235,2.255],[2.360,2.380],method='ratio',sample='average',**kwargs)
        inds[0]= num/den
        errs[0]= inds[0]*numpy.sqrt((er1/num)**2+(er2/den)**2)
    elif ('slesnick' in set.lower()):
        reference = 'Indices from Slesnick et al. (2004)'
        refcode = 'SEL04'
        names = ['H2O-1','H2O-2','FeH']
        inds = numpy.zeros(len(names))
        errs = numpy.zeros(len(names))
        inds[0],errs[0] = measureIndex(sp,[1.335,1.345],[1.295,1.305],method='ratio',sample='average',**kwargs)
        inds[1],errs[1] = measureIndex(sp,[2.035,2.045],[2.145,2.155],method='ratio',sample='average',**kwargs)
        inds[2],errs[2] = measureIndex(sp,[1.1935,1.2065],[1.2235,1.2365],method='ratio',sample='average',**kwargs)
    elif ('testi' in set.lower()):
        reference = 'Indices from Testi et al. (2001)'
        refcode = 'TES01'
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
        refcode = 'TOK99'
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
    nsamples = kwargs.get('nsamples',100)

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


    mh = numpy.ones(nsamples)*coeff_mh[0]+\
        (numpy.random.normal(nai,nai_e,nsamples)/numpy.random.normal(h2ok2,h2ok2_e,nsamples))*coeff_mh[1]+\
        (numpy.random.normal(cai,cai_e,nsamples)/numpy.random.normal(h2ok2,h2ok2_e,nsamples))*coeff_mh[2]

    return numpy.nanmean(mh), numpy.sqrt(numpy.nanstd(mh)**2+mh_unc**2)


def properCoordinates(c):
    '''
    :Purpose: Converts various coordinate forms to the proper SkyCoord format. Convertible forms include lists and strings.

    :param c: coordinate to be converted. Can be a list (ra, dec) or a string.

    :Example:
    >>> import splat
    >>> print splat.properCoordinates([104.79, 25.06])
        <SkyCoord (ICRS): ra=104.79 deg, dec=25.06 deg>
    >>> print splat.properCoordinates('06:59:09.60 +25:03:36.0')
        <SkyCoord (ICRS): ra=104.79 deg, dec=25.06 deg>
    >>> print splat.properCoordinates('J06590960+2503360')
        <SkyCoord (ICRS): ra=104.79 deg, dec=25.06 deg>
    '''
    if isinstance(c,SkyCoord):
        return c
    elif isinstance(c,list):
        return SkyCoord(c[0]*u.deg,c[1]*u.deg,frame='icrs')
# input is sexigessimal string
    elif isinstance(c,str):
        if c[0] == 'J':
            return designationToCoordinate(c)
        else:
            return SkyCoord(c,'icrs', unit=(u.hourangle, u.deg))
    else:
        raise ValueError('\nCould not parse input format\n\n')


def readSpectrum(*args,**kwargs):
    '''
    .. will come back to this one
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
        file = checkLocal(kwargs['filename'])
        if file=='':
            file = checkLocal(kwargs['folder']+os.path.basename(kwargs['filename']))
            if file=='':
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
        with fits.open(file) as data:
            if 'NAXIS3' in data[0].header.keys():
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
#    print(flux, numpy.median(flux))
    if (catchSN):
          w = numpy.where(flux > stats.nanmedian(flux))
          if (stats.nanmedian(flux[w]/noise[w]) < 1.):
              noise = flux/noise
              w = numpy.where(numpy.isnan(noise))
              noise[w] = stats.nanmedian(noise)

# clean up
#    if url != '' and not local:
#        os.remove(os.path.basename(TMPFILENAME))

    return {'wave':wave,'flux':flux,'noise':noise,'header':header}



def redden(sp, **kwargs):
    '''
    Description:
      Redden a spectrum based on an either Mie theory or a standard interstellar profile
      using Cardelli, Clayton, and Mathis (1989 ApJ. 345, 245)

    **Usage**

       >>> import splat
       >>> sp = splat.Spectrum(10001)                   # read in a source
       >>> spr = splat.redden(sp,av=5.,rv=3.2)          # redden to equivalent of AV=5

    **Note**
      This routine is still in beta form; only the CCM89 currently works

    '''
    w = sp.wave.value                           # assuming in microns!
    av = kwargs.get('av',0.0)


    if kwargs.get('mie',False):                 # NOT CURRENTLY FUNCTIONING
        a = kwargs.get('a',10.)                 # grain size
        n = kwargs.get('n',1.33)                # complex index of refraction
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
    else:
        x = 1./w
        a = 0.574*(x**1.61)
        b = -0.527*(x**1.61)
        rv = kwargs.get('rv',3.1)
        absfrac = 10.**(-0.4*av*(a+b/rv))

    if kwargs.get('normalize',False):
        absfrac = absfrac/numpy.median(absfrac)

    print(tau0, min(tau), max(tau), max(absfrac), min(absfrac))
    spabs = splat.Spectrum(wave=w,flux=absfrac)
    return sp*spabs




def test():
    '''
    :Purpose: Tests the SPLAT Code
    :Checks the following:

        - If you are online and can see the SPLAT website
        - If you have access to unpublished spectra
        - If you can search for and load a spectrum
        - If ``searchLibrary`` functions properly
        - If index measurement routines functions properly
        - If classification routines function properly
        - If ``typeToTeff`` functions properly
        - If flux calibration and normalization function properly
        - If ``loadModel`` functions properly
        - If ``compareSpectra`` functions properly
        - If ``plotSpectrum`` functions properly
    '''

    test_src = 'Random'

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
#    sp = getSpectrum(shortname=test_src)[0]
    sp = getSpectrum(spt=['L5','T5'],lucky=True)[0]
    sp.info()
    sys.stderr.write('\n...getSpectrum successful\n')

# check searchLibrary
    list = searchLibrary(young=True,output='DATA_FILE')
    sys.stderr.write('\n{} young spectra in the SPL  ...searchLibrary successful\n'.format(len(list)))

# check index measurement
    ind = measureIndexSet(sp,set='burgasser')
    sys.stderr.write('\nSpectral indices for '+test_src+':\n')
    for k in ind.keys():
        print('\t{:s}: {:4.3f}+/-{:4.3f}'.format(k, ind[k][0], ind[k][1]))
    sys.stderr.write('...index measurement successful\n')

# check classification
    spt, spt_e = classifyByIndex(sp,set='burgasser')
    sys.stderr.write('\n...index classification of '+test_src+' = {:s}+/-{:2.1f}; successful\n'.format(spt,spt_e))

    spt, spt_e = classifyByStandard(sp,method='kirkpatrick')
    sys.stderr.write('\n...standard classification of '+test_src+' = {:s}+/-{:2.1f}; successful\n'.format(spt,spt_e))

    grav = classifyGravity(sp)
    sys.stderr.write('\n...gravity class of '+test_src+' = {}; successful\n'.format(grav))

# check SpT -> Teff
    teff, teff_e = typeToTeff(spt,unc=spt_e)
    sys.stderr.write('\n...Teff of '+test_src+' = {:.1f}+/-{:.1f} K; successful\n'.format(teff,teff_e))

# check flux calibration
    sp.normalize()
    sp.fluxCalibrate('2MASS J',15.0,apparent=True)
    mag,mag_e = filterMag(sp,'MKO J')
    sys.stderr.write('\n...apparent magnitude MKO J = {:3.2f}+/-{:3.2f} from 2MASS J = 15.0; filter calibration successful\n'.format(mag,mag_e))

# check models
#   mdl = loadModel(teff=1000,logg=5.0,set='BTSettl2008')
    logg = 5.2
    if grav == 'VL-G':
        logg = 4.2
    elif grav == 'INT-G':
        logg = 4.6
    mdl = loadModel(teff=teff,logg=logg,set='BTSettl2008')
    sys.stderr.write('\n...interpolated model generation successful\n')

# check normalization
    sys.stderr.write('\n...normalization successful\n')

# check compareSpectrum
    chi, scale = compareSpectra(sp,mdl,mask_standard=True,stat='chisqr')
    sys.stderr.write('\nScaling model: chi^2 = {}, scale = {}'.format(chi,scale))
    sys.stderr.write('\n...compareSpectra successful\n'.format(chi,scale))

# check plotSpectrum
    mdl.scale(scale)
    plotSpectrum(sp,mdl,colors=['k','r'],title='If this appears everything is OK: close window')
    sys.stderr.write('\n...plotSpectrum successful\n')
    sys.stderr.write('\n>>>>>>>>>>>> SPLAT TEST SUCCESSFUL; HAVE FUN! <<<<<<<<<<<<\n\n')


def typeToColor(spt,color, **kwargs):
    """
    :Purpose: Takes a spectral type and optionally a color (string) and returns the typical color of the source. 
    :param spt: string or integer of the spectral type
    :param color: string indicating color; e.g., color='i-z' (note that case does not matter)
    :type color: optional, default = 'J-K'
    :param ref: Abs Mag/SpT relation used to compute the absolute magnitude. Options are:

        - *skrzypek* (default): Color trends from `Skryzpek et al. (2015) <http://adsabs.harvard.edu/abs/2015A%26A...574A..78S>`_.
          Spectral type range is M5 to T8
          Colors include i-z, z-Y, Y-J, J-H, H-K, K-W1, W1-W2, and combinations therein.


    :type ref: optional, default = 'dupuy'
    :param nsamples: number of Monte Carlo samples for error computation
    :type nsamples: optional, default = 100
    :param unc: uncertainty of ``spt``; if included, returns a tuple with color and uncertainty
    :type unc: optional, default = 0.
    :param verbose: Give feedback while in operation
    :type verbose: optional, default = False
    :Example:
        >>> import splat
        >>> print splat.typeToColor('L3', 'J-K')
            XXXX
        >>> print splat.typeToMag('M5', 'i-K', ref = 'skrzypek', unc=0.5)
            (XXXX, XXXX)
        >>> print splat.typeToMag('M0', 'i-z', ref = 'skrzypek')
            Spectral Type is out of range for Color/SpT trends from Skryzpek et al. (2015)
            nan
    """

#Keywords
    nsamples = kwargs.get('nsamples', 100)
    ref = kwargs.get('ref', 'skrzypek')
    unc = kwargs.get('unc', 0.)

#Convert spectral type string to number
    if (type(spt) == str):
        spt = typeToNum(spt, uncertainty=unc)
    else:
        spt = copy.deepcopy(spt)

#Faherty
    if (ref.lower() == 'skrzypek'):
        reference = 'Skrzypek et al. (2015)'
        rng = [15,38]
        filters = ['i','z','y','j','h','k','w1','w2']
        values = { \
            'i-z': [0.91,1.45,1.77,1.93,1.99,2.01,2.02,2.04,2.1,2.2,2.33,2.51,2.71,2.93,3.15,3.36,3.55,3.7,3.82,3.9,3.95,3.98,4.01,4.08], \
            'z-y': [0.47,0.6,0.7,0.77,0.82,0.86,0.88,0.9,0.92,0.94,0.97,1.0,1.04,1.09,1.16,1.23,1.33,1.43,1.55,1.68,1.81,1.96,2.11,2.26], \
            'y-j': [0.55,0.67,0.78,0.87,0.96,1.04,1.11,1.18,1.23,1.27,1.31,1.33,1.35,1.21,1.2,1.19,1.19,1.18,1.18,1.17,1.16,1.16,1.15,1.15], \
            'j-h': [0.45,0.53,0.56,0.58,0.6,0.63,0.67,0.73,0.79,0.86,0.91,0.96,0.97,0.96,0.9,0.8,0.65,0.46,0.25,0.02,-0.19,-0.35,-0.43,-0.36], \
            'h-k': [0.32,0.39,0.44,0.47,0.51,0.54,0.58,0.63,0.67,0.71,0.74,0.75,0.75,0.71,0.65,0.56,0.45,0.31,0.16,0.01,-0.11,-0.19,-0.2,-0.09], \
            'k-w1': [0.11,0.22,0.25,0.26,0.27,0.29,0.33,0.4,0.48,0.56,0.65,0.72,0.77,0.79,0.79,0.76,0.71,0.65,0.59,0.55,0.54,0.59,0.7,0.9], \
            'w1-w2': [0.17,0.21,0.24,0.26,0.27,0.27,0.28,0.28,0.29,0.3,0.32,0.36,0.41,0.48,0.57,0.68,0.82,0.99,1.19,1.43,1.7,2.02,2.38,2.79]}

    else:
        sys.stderr.write('\n Do not have color set {}\n\n'.format(ref))
        return numpy.nan

    if kwargs.get('verbose',False):
        print('\nUsing the SpT/color trends from {}\n'.format(reference))

# spectral type array
    if (rng[0] <= spt <= rng[1]):

# fill in extra colors - a little inefficient right now  
        if color.lower() not in values.keys():
            for i in numpy.arange(len(values.keys())):
                for a in values.keys():
                    for b in values.keys():
                        f1 = a.split('-')
                        f2 = b.split('-')
                        if f1[-1] == f2[0]:
                            k = '{}-{}'.format(f1[0],f2[-1])
                            if k not in values.keys():
                                values[k] = [sum(x) for x in zip(values[a], values[b])]

        if color.lower() in values.keys():
            f = interp1d(numpy.arange(rng[0],rng[1]+1),values[color.lower()],bounds_error=False,fill_value=0.)
            if (unc > 0.):
                vals = [f(x) for x in numpy.random.normal(spt, unc, nsamples)]
                return numpy.nanmean(vals), numpy.nanstd(vals)
            else:
                return float(f(spt))
        else:
            sys.stderr.write('\n Color {} is not in reference set for {}\n\n'.format(color,reference))
            return numpy.nan

    else:
        sys.stderr.write('\n Spectral type {} is outside the range for reference set {}\n\n'.format(typeToNum(spt),reference))
        return numpy.nan




def typeToMag(spt, filt, **kwargs):
    """
    :Purpose: Takes a spectral type and a filter, and returns absolute magnitude
    :param spt: string or integer of the spectral type
    :param filter: filter of the absolute magnitude. Options are MKO K, MKO H, MKO J, MKO Y, MKO LP, 2MASS J, 2MASS K, or 2MASS H
    :param nsamples: number of Monte Carlo samples for error computation
    :type nsamples: optional, default = 100
    :param unc: uncertainty of ``spt``
    :type unc: optional, default = 0.
    :param ref: Abs Mag/SpT relation used to compute the absolute magnitude. Options are:

        - *burgasser*: Abs Mag/SpT relation from `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_.
          Allowed spectral type range is L0 to T8, and allowed filters are MKO K.
        - *faherty*: Abs Mag/SpT relation from `Faherty et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...752...56F>`_.
          Allowed spectral type range is L0 to T8, and allowed filters are MKO J, MKO H and MKO K.
        - *dupuy*: Abs Mag/SpT relation from `Dupuy & Liu (2012) <http://adsabs.harvard.edu/abs/2012ApJS..201...19D>`_.
          Allowed spectral type range is M6 to T9, and allowed filters are MKO J, MKO Y, MKO H, MKO K, MKO LP, 2MASS J, 2MASS H, and 2MASS K.
        - *filippazzo*: Abs Mag/SpT relation from Filippazzo et al. (2015). Allowed spectral type range is M6 to T9, and allowed filters are 2MASS J and WISE W2.


    :type ref: optional, default = 'dupuy'
    :Example:
        >>> import splat
        >>> print splat.typeToMag('L3', '2MASS J')
            (12.730064813273996, 0.4)
        >>> print splat.typeToMag(21, 'MKO K', ref = 'burgasser')
            (10.705292820099999, 0.26)
        >>> print splat.typeToMag(24, '2MASS J', ref = 'faherty')
            Invalid filter given for Abs Mag/SpT relation from Faherty et al. (2012)
            (nan, nan)
        >>> print splat.typeToMag('M0', '2MASS H', ref = 'dupuy')
            Spectral Type is out of range for Abs Mag/SpT relation from Dupuy & Liu (2012) Abs Mag/SpT relation
            (nan, nan)
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
            '2MASS H': {'fitunc': 0.40, 'range': [16., 39.], \
                'coeff': [-.00000111499, .000129363, -.00580847, .129202, -1.50370, 9.00279, -11.7526]}, \
            '2MASS KS': {'fitunc': 0.43, 'range':[16., 39.], \
                'coeff': [1.06693e-4, -6.42118e-3, 1.34163e-1, -8.67471e-1, 1.10114e1]}, \
            'WISE W1': {'fitunc': 0.39, 'range':[16., 39.], \
                'coeff': [1.58040e-5, -3.33944e-4, -4.38105e-3, 3.55395e-1, 7.14765]}, \
            'WISE W2': {'fitunc': 0.35, 'range':[16., 39.], \
                'coeff': [1.78555e-5, -8.81973e-4, 1.14325e-2, 1.92354e-1, 7.46564]}}

    elif (ref.lower() == 'filippazzo'):
        reference = 'Abs Mag/SpT relation from Filippazzo et al. (2015)'
        sptoffset = 10.
        coeffs = { \
            '2MASS J': {'fitunc': 0.40, 'range': [16., 39.], \
                'coeff': [3.478e-5, -2.684e-3, 7.771e-2, -1.058e0, 7.157e0, -8.350e0]}, \
            'WISE W2': {'fitunc': 0.40, 'range': [16., 39.], \
                'coeff': [8.190e-6, -6.938e-4, 2.283e-2, -3.655e-1, 3.032e0, -5.043e-1]}}

    else:
        sys.stderr.write('\nInvalid Abs Mag/SpT relation given: %s\n' % ref)
        return numpy.nan, numpy.nan

    if (filt.upper() in coeffs.keys()) == 1:
        for f in coeffs.keys():
            if filt.upper() == f:
                coeff = coeffs[f]['coeff']
                fitunc = coeffs[f]['fitunc']
                rng = coeffs[f]['range']
    else:
        sys.stderr.write('\n Invalid filter {} given for {}\n'.format(filt,reference))
        return numpy.nan, numpy.nan

# compute magnitude if its in the right spectral type range
    if (rng[0] <= spt <= rng[1]):
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
    '''
    :Purpose: Converts between string and numeric spectral types, and vise versa.
    :param input: Spectral type to convert. Can convert a number or a string from 0 (K0) and 49.0 (Y9).
    :param error: magnitude of uncertainty. ':' for uncertainty > 1 and '::' for uncertainty > 2.
    :type error: optional, default = ''
    :param uncertainty: uncertainty of spectral type
    :type uncertainty: optional, default = 0
    :param subclass: subclass of object. Options include:

        - *sd*: object is a subdwarf
        - *esd*: object is an extreme subdwarf
        - *usd*: object is an ultra subdwarf

    :type subclass: optional, default = ''
    :param lumclass: luminosity class of object represented by roman numerals
    :type lumclass: optional, default = ''
    :param ageclass: age class of object
    :type ageclass: optional, default = ''
    :param colorclass: color class of object
    :type colorclass: optional, default = ''
    :param peculiar: if object is peculiar or not
    :type peculiar: optional, default = False

    .. not too sure how colorclass and ageclass work

    :Example:
        >>> import splat
        >>> print splat.typeToNum(30)
            T0.0
        >>> print splat.typeToNum('T0.0')
            30.0
        >>> print splat.typeToNum(27, peculiar = True, uncertainty = 1.2, lumclass = 'II')
            L7.0IIp:
        >>> print splat.typeToNum(50)
            Spectral type number must be between 0 (K0) and 49.0 (Y9)
            nan
    '''
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
            print('Spectral type number must be between 0 ({}0) and {} ({}9)'.format(spletter[0],len(spletter)*10.-1.,spletter[-1]))
            output = numpy.nan

# spectral type -> number
    elif isinstance(input,str):
        if (sys.version_info.major == 2):
            input = string.split(input,sep='+')[0]    # remove +/- sides
            input = string.split(input,sep='-')[0]    # remove +/- sides
        else:
            input = input.split('+')[0]    # remove +/- sides
            input = input.split('-')[0]    # remove +/- sides

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
#            print('Only spectral classes {} are handled with this routine'.format(spletter))
            output = numpy.nan
# none of the above - return the input
    else:
        output = input
    return output


def typeToTeff(input, **kwargs):
    '''
    :Purpose: Returns an effective temperature (Teff) and its uncertainty for a given spectral type
    :param input: Spectral type; can be a number or a string from 0 (K0) and 49.0 (Y9).
    :param uncertainty: uncertainty of spectral type
    :type uncertainty: optional, default = 0.001
    :param unc: same as ``uncertainty``
    :type unc: optional, default = 0.001
    :param spt_e: same as ``uncertainty``
    :type spt_e: optional, default = 0.001
    :param ref: Teff/SpT relation used to compute the effective temperature. Options are:

        - *golimowski*: Teff/SpT relation from `Golimowski et al. (2004) <http://adsabs.harvard.edu/abs/2004AJ....127.3516G>`_.
          Allowed spectral type range is M6 to T8.
        - *looper*: Teff/SpT relation from `Looper et al. (2008) <http://adsabs.harvard.edu/abs/2008ApJ...685.1183L>`_.
          Allowed spectral type range is L0 to T8.
        - *stephens*: Teff/SpT relation from `Stephens et al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...702..154S>`_.
          Allowed spectral type range is M6 to T8 and uses alternate coefficients for L3 to T8.
        - *marocco*: Teff/SpT relation from `Marocco et al. (2013) <http://adsabs.harvard.edu/abs/2013AJ....146..161M>`_.
          Allowed spectral type range is M7 to T8.
        - *filippazzo*: Teff/SpT relation from Filippazzo et al. (2015). Allowed spectral type range is M6 to T9.

    :type ref: optional, default = 'stephens2009'
    :param set: same as ``ref``
    :type set: optional, default = 'stephens2009'
    :param method: same as ``ref``
    :type method: optional, default = 'stephens2009'
    :param nsamples: number of samples to use in Monte Carlo error estimation
    :type nsamples: optional, default = 100

    :Example:
        >>> import splat
        >>> print splat.typeToTeff(20)
            (2233.4796740905499, 100.00007874571999)
        >>> print splat.typeToTeff(20, unc = 0.3, ref = 'golimowski')
            (2305.7500497902788, 127.62548366132124)
    '''
# keywords
    nsamples = kwargs.get('nsamples',100)
    unc = kwargs.get('uncertainty',0.001)
    unc = kwargs.get('unc',unc)
    unc = kwargs.get('spt_e',unc)
    ref = kwargs.get('ref','stephens2009')
    ref = kwargs.get('set',ref)
    ref = kwargs.get('method',ref)

# convert spectral type string to number
    if (type(input) == str):
        spt = typeToNum(input,uncertainty=unc)
    else:
        spt = copy.deepcopy(input)

#    if spt < 20. and 'marocco' not in ref.lower():
#        ref='stephens2009'

# choose among possible options

# Golimowski et al. (2004, AJ, 127, 3516)
    if ('golimowski' in ref.lower()):
        reference = 'Teff/SpT relation from Golimowski et al. (2004)'
        sptoffset = 10.
        coeff = [9.5373e-4,-9.8598e-2,4.0323,-8.3099e1,9.0951e2,-5.1287e3,1.4322e4]
        range = [16.,38.]
        fitunc = 124.

# Looper et al. (2008, ApJ, 685, 1183)
    elif ('looper' in ref.lower()):
        reference = 'Teff/SpT relation from Looper et al. (2008)'
        sptoffset = 20.
        coeff = [9.084e-4,-4.255e-2,6.414e-1,-3.101,1.950,-108.094,2319.92]
        range = [20.,38.]
        fitunc = 87.

# Stephens et al. (2009, ApJ, 702, 1545); using OPT/IR relation for M6-T8
# plus alternate coefficients for L3-T8
    elif ('stephens' in ref.lower()):
        reference = 'Teff/SpT relation from Stephens et al. (2009)'
        sptoffset = 10.
        coeff = [-0.0025492,0.17667,-4.4727,54.67,-467.26,4400.]
        range = [16.,38.]
        fitunc = 100.
        coeff_alt = [-0.011997,1.2315,-50.472,1031.9,-10560.,44898.]
        range_alt = [23.,38.]

# Marocco et al. (2013, AJ, 146, 161)
    elif ('marocco' in ref.lower()):
        reference = 'Teff/SpT relation from Marocco et al. (2013)'
        sptoffset = 10.
        coeff = [7.4211e-5,-8.43736e-3,3.90319e-1,-9.46896,129.141,-975.953,3561.47,-1613.82]
        range = [17.,38.]
        fitunc = 140.

    elif ('filippazzo' in ref.lower()):
        reference = 'Teff/SpT relation from Filippazzo et al. (2015)'
        sptoffset = 10.
        coeff = [1.546e-4, -1.606e-2, 6.318e-1, -1.191e1, 1.155e2, -7.005e2, 4.747e3]
        range = [16., 39.]
        fitunc = 113.

    else:
        sys.stderr.write('\nInvalid Teff/SpT relation given ({})\n'.format(ref))
        return numpy.nan, numpy.nan

    if (range[0] <= spt <= range[1]):
        vals = numpy.polyval(coeff,numpy.random.normal(spt-sptoffset,unc,nsamples))
        if ('stephens' in ref.lower()):
            if (range_alt[0] <= spt <= range_alt[1]):
                vals = numpy.polyval(coeff_alt,numpy.random.normal(spt-sptoffset,unc,nsamples))
        teff = numpy.nanmean(vals)
        teff_e = (numpy.nanstd(vals)**2+fitunc**2)**0.5
        return teff, teff_e
    else:
        sys.stderr.write('\nSpectral Type is out of range for {:s} Teff/SpT relation\n'.format(reference))
        return numpy.nan, numpy.nan


def weightedMeanVar(vals, winp, *args, **kwargs):
    '''
    :Purpose: Computes weighted mean of an array of values through various methods. Returns weighted mean and weighted uncertainty.
    :param vals: array of values
    :param winp: array of weights associated with ``vals``
    :param method: input type of weights. Default is where ``winp`` is the actual weights of ``vals``. Options include:

        - *uncertainty*: uncertainty weighting, where ``winp`` is the uncertainties of ``vals``
        - *ftest*: ftest weighting, where ``winp`` is the chi squared values of ``vals``

    :type method: optional, default = ''
    :param weight_minimum: minimum possible weight value
    :type weight_minimum: optional, default = 0.
    :param dof: effective degrees of freedom
    :type dof: optional, default = len(vals) - 1

    .. note:: When using ``ftest`` method, extra ``dof`` value is required

    :Example:
        >>> import splat
        >>> print splat.weightedMeanVar([3.52, 5.88, 9.03], [0.65, 0.23, 0.19])
            (5.0057009345794379, 4.3809422657000594)
        >>> print splat.weightedMeanVar([3.52, 5.88, 9.03], [1.24, 2.09, 2.29], method = 'uncertainty')
            (5.0069199363443841, 4.3914329968409946)
    '''

    method = kwargs.get('method','')
    minwt = kwargs.get('weight_minimum',0.)
    dof = kwargs.get('dof',len(vals)-1)
    if (numpy.nansum(winp) <= 0.):
        weights = numpy.ones(len(vals))
    if isinstance(winp,astropy.units.quantity.Quantity):
        winput = winp.value
    else:
        winput = copy.deepcopy(winp)

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
    if (method == 'uncertainty'):
        var+=numpy.nansum([w**2 for w in winput])/(len(winput)**2)

    return mn,numpy.sqrt(var)


# run test program if calling from command line
if __name__ == '__main__':
    splat.test()
