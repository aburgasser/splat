# -*- coding: utf-8 -*-
from __future__ import print_function

# this is the test function set for splat plot functions

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
from splat.utilities import *
import splat

output_folder = splat.SPLAT_PATH+splat.DOCS_FOLDER+'/_images/'


#####################
# TESTING FUNCTIONS #
#####################


def test_plotSpectrum():
# plotSpectrum
    pass

def test_plotBatch():
# plotBatch
    pass

def test_plotSequence():
# plotSequence
	sp = splat.getSpectrum(lucky=True)[0]
	splot.plotSequence(sp,output=output_folder+'example_plotsequence_dwarf.png')
	sp = splat.getSpectrum(subdwarf=True,lucky=True)[0]
	splot.plotSequence(sp,type='sd',output=output_folder+'example_plotsequence_subdwarf.png')
	sp = splat.getSpectrum(young=True,lucky=True)[0]
	splot.plotSequence(sp,type='lowg',output=output_folder+'example_plotsequence_lowg.png')
	return True

def test_plotIndices():
# plotIndices
    pass

def test_plotSED():
# plotSED
    pass

def test_plotMap():
# plotSED
    s = splat.searchLibrary(young=True)
	c = [splat.properCoordinates(x) for x in s['DESIGNATION']]
	splot.plotMap(c,output=output_folder+'example_plotmap_young.png')

	sm = splat.searchLibrary(spt=['M0','M9.5'],jmag=13)
	sl = splat.searchLibrary(spt=['L0','L9.5'],jmag=13)
	cm = [splat.properCoordinates(x) for x in sm['DESIGNATION']]
	cl = [splat.properCoordinates(x) for x in sl['DESIGNATION']]
	splot.plotMap(cm,cl,colors=['k','r'],markers=['o','o'],alpha=[0.3,0.9],output=output_folder+'example_plotmap_brightml.png')

