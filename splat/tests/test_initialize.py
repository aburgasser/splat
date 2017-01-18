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
from numpy.testing.utils import assert_allclose

# splat functions and constants
from splat.initialize import *
import splat


#####################
# TESTING FUNCTIONS #
#####################


def test_folders_files():
    assert os.access(SPLAT_PATH,os.R_OK)
    assert os.access(SPLAT_PATH+DATA_FOLDER,os.R_OK)
    assert os.access(SPLAT_PATH+FILTER_FOLDER,os.R_OK)
    assert os.access(SPLAT_PATH+SPECTRAL_MODEL_FOLDER,os.R_OK)
    assert os.access(SPLAT_PATH+EVOLUTIONARY_MODEL_FOLDER,os.R_OK)
    assert os.access(SPLAT_PATH+DB_FOLDER,os.R_OK)
    assert os.access(SPLAT_PATH+DB_FOLDER+DB_SOURCES_FILE,os.R_OK)
    assert os.access(SPLAT_PATH+DB_FOLDER+DB_SPECTRA_FILE,os.R_OK)
    assert os.access(SPLAT_PATH+DB_FOLDER+DB_PHOTOMETRY_FILE,os.R_OK)
    assert os.access(SPLAT_PATH+DB_FOLDER+BIBFILE,os.R_OK)
    if os.access(SPLAT_PATH+ACCESS_FILE,os.R_OK) != True:
    	print('Warning: {} is not found'.format(ACCESS_FILE))

def test_filter_files():
# getBibTex
# getBibTexOnline
    for f in list(FILTERS.keys()):
    	assert os.access(SPLAT_PATH+FILTER_FOLDER+FILTERS[f]['file'],os.R_OK)

def test_standards_present():
# getBibTex
# getBibTexOnline
    for s in list(STDS_DWARF_SPEX_KEYS.keys()):
    	sp = splat.Spectrum(STDS_DWARF_SPEX_KEYS[s])
    	assert type(sp) == splat.Spectrum
    for s in list(STDS_SD_SPEX_KEYS.keys()):
    	sp = splat.Spectrum(STDS_SD_SPEX_KEYS[s])
    	assert type(sp) == splat.Spectrum
    for s in list(STDS_ESD_SPEX_KEYS.keys()):
    	sp = splat.Spectrum(STDS_ESD_SPEX_KEYS[s])
    	assert type(sp) == splat.Spectrum


def test_parse_bibtex():
# bibTexParser
    pass

