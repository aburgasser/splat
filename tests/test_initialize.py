# -*- coding: utf-8 -*-
from __future__ import print_function

# this is the test splat initialization

# imports - internal
import copy
import glob
import os
# # imports - external
import numpy
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units
from astropy import coordinates as coord      # coordinate conversion
from astropy.io import fits
from numpy.testing import assert_allclose

# splat functions and constants
import splat


#####################
# TESTING FUNCTIONS #
#####################


def test_folders_files():
    assert os.access(splat.SPLAT_PATH,os.R_OK)
    print(splat.SPLAT_PATH)
    print(os.path.join(splat.SPLAT_PATH,splat.DATA_FOLDER))
    assert os.access(os.path.join(splat.SPLAT_PATH,splat.DATA_FOLDER),os.R_OK)
    assert os.access(os.path.join(splat.SPLAT_PATH,splat.FILTER_FOLDER),os.R_OK)
    assert os.access(os.path.join(splat.SPLAT_PATH,splat.SPECTRAL_MODEL_FOLDER),os.R_OK)
    assert os.access(os.path.join(splat.SPLAT_PATH,splat.EVOLUTIONARY_MODEL_FOLDER),os.R_OK)
    assert os.access(os.path.join(splat.SPLAT_PATH,splat.DB_FOLDER),os.R_OK)
    assert os.access(os.path.join(splat.SPLAT_PATH,splat.DB_FOLDER,splat.DB_SOURCES_FILE),os.R_OK)
    assert os.access(os.path.join(splat.SPLAT_PATH,splat.DB_FOLDER,splat.DB_SPECTRA_FILE),os.R_OK)
#    assert os.access(SPLAT_PATH+DB_FOLDER+DB_PHOTOMETRY_FILE,os.R_OK)
    assert os.access(os.path.join(splat.SPLAT_PATH,splat.DB_FOLDER,splat.BIBFILE),os.R_OK)
    # if os.access(os.path.join(splat.SPLAT_PATH,splat.ACCESS_FILE),os.R_OK) != True:
    # 	print('Warning: {} is not found'.format(splat.ACCESS_FILE))
    return

def test_filter_files():
# getBibTex
# getBibTexOnline
    for f in list(splat.FILTERS.keys()):
    	assert os.access(os.path.join(splat.SPLAT_PATH,splat.FILTER_FOLDER,splat.FILTERS[f]['file']),os.R_OK)
    return

def test_standards_present():
    for s in list(splat.STDS_DWARF_SPEX_KEYS.keys()):
    	sp = splat.Spectrum(splat.STDS_DWARF_SPEX_KEYS[s])
    	assert type(sp) == splat.Spectrum
    for s in list(splat.STDS_SD_SPEX_KEYS.keys()):
    	sp = splat.Spectrum(splat.STDS_SD_SPEX_KEYS[s])
    	assert type(sp) == splat.Spectrum
    for s in list(splat.STDS_ESD_SPEX_KEYS.keys()):
    	sp = splat.Spectrum(splat.STDS_ESD_SPEX_KEYS[s])
    	assert type(sp) == splat.Spectrum
    return


def test_parse_bibtex():
# getBibTex
# getBibTexOnline
# bibTexParser
    pass

