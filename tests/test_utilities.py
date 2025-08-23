# -*- coding: utf-8 -*-
from __future__ import print_function

# this is the test function set for splat utility functions

# imports - internal
import copy
# # imports - external
import numpy
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units
from astropy import coordinates as coord      # coordinate conversion
from astropy.io import fits
from numpy.testing import assert_allclose
import os
# splat functions and constants
import splat
from splat.initialize import *
#import splat as splat

# things to test
# isNumber 


#####################
# TESTING FUNCTIONS #
#####################


# def test_files():
# # checkFile
# # checkAccess
# # checkLocal
# # checkOnline
# # checkOnlineFile
#     # assert splat.checkOnline()
#     # assert splat.checkOnlineFile()
#     # assert len(splat.checkOnlineFile('/Reference/SpectralModels/'))>0
#     # assert len(splat.checkOnlineFile('/Reference/Spectra/11221_10235.fits'))>0
#     # assert len(splat.checkOnlineFile('/Reference/Spectra/11220_10166.fits'))>0
#     assert len(splat.checkLocal(FILTER_FOLDER))>0
#     assert len(splat.checkLocal(splat.SPECTRAL_MODEL_FOLDER))>0
#     assert len(splat.checkLocal(splat.EVOLUTIONARY_MODEL_FOLDER))>0
#     assert len(splat.checkLocal(splat.TELLURIC_MODEL_FOLDER))>0
#     assert len(splat.checkLocal(splat.CITATION_RESOURCES_FOLDER))>0
#     assert len(splat.checkLocal(splat.DATA_FOLDERS))>0
#     assert len(splat.checkLocal(os.path.join(splat.DATA_FOLDERS,splat.DB_SPECTRA_FILE)))>0
# #    assert len(splat.checkLocal(os.path.join(splat.DATA_FOLDER,'SPEX-PRISM',splat.DB_SOURCES_FILE)))>0
#     assert len(splat.checkLocal(os.path.join(splat.DATA_FOLDER,'11221_10235.fits')))>0
# #    assert len(splat.checkLocal(os.path.join(splat.DATA_FOLDER,'11220_10166.fits')))==0
#     assert splat.checkFile('11221_10235.fits')
# #    assert splat.checkFile('11220_10166.fits')
# #    assert splat.checkAccess()

def test_isnumber():
# checkKeys
    assert splat.isNumber(3)
    assert splat.isNumber(3.0)
    assert splat.isNumber(3.0*u.arcsec)
    assert splat.isNumber('3')
    assert not splat.isNumber('hello')
    assert not splat.isNumber(numpy.nan)
    assert not splat.isNumber(True)

def test_coordinate_conversion():
# coordinateToDesignation 
# designationToCoordinate 
# designationToShortName 
# properCoordinates 
    d = 'J05591914-1404488'
    assert splat.designationToShortName(d) == 'J0559-1404'
    c = splat.designationToCoordinate(d)
    assert isinstance(c,coord.sky_coordinate.SkyCoord)
    assert_allclose(c.ra.value,89.82975,rtol=1.e-5)
    assert_allclose(c.dec.value,-14.08002,rtol=1.e-5)
    c1 = splat.properCoordinates(d)
    assert isinstance(c1,coord.sky_coordinate.SkyCoord)
    assert_allclose(c1.ra.value,c.ra.value,rtol=1.e-5)
    assert_allclose(c1.dec.value,c.dec.value,rtol=1.e-5)
    d1 = splat.coordinateToDesignation(c)
    assert splat.designationToShortName(d1) == 'J0559-1404'
    c2 = splat.designationToCoordinate(d1)
    assert isinstance(c2,coord.sky_coordinate.SkyCoord)
    assert_allclose(c2.ra.value,c.ra.value,rtol=1.e-5)
    assert_allclose(c2.dec.value,c.dec.value,rtol=1.e-5)  
    c3 = splat.properCoordinates('05 59 19.14 -14 04 48.8') 
    assert_allclose(c3.ra.value,c.ra.value,rtol=1.e-4)
    assert_allclose(c3.dec.value,c.dec.value,rtol=1.e-4)  
    c4 = splat.properCoordinates('05:59:19.14 -14:04:48.8') 
    assert_allclose(c4.ra.value,c.ra.value,rtol=1.e-4)
    assert_allclose(c4.dec.value,c.dec.value,rtol=1.e-4)  
    c5 = splat.properCoordinates('05h59m19.14s -14d04m48.8s') 
    assert_allclose(c5.ra.value,c.ra.value,rtol=1.e-4)
    assert_allclose(c5.dec.value,c.dec.value,rtol=1.e-4)  
    c6 = splat.properCoordinates([89.82975,-14.08002]) 
    assert_allclose(c6.ra.value,c.ra.value,rtol=1.e-4)
    assert_allclose(c6.dec.value,c.dec.value,rtol=1.e-4)  


def test_type_conversion():
# typeToNum 
    assert splat.typeToNum('L5') == 25.
    assert splat.typeToNum(25.0) == 'L5.0'
    splat.typeToNum(25.,subclass='sd') == 'sdL5.0'
    assert splat.typeToNum('K0') == 0.
    assert splat.typeToNum('m0') == 10.
    assert splat.typeToNum('L0') == 20.
    assert splat.typeToNum('t0') == 30.
    assert splat.typeToNum('Y0') == 40.
#    assert splat.typeToNum('A0') == numpy.nan
#    assert splat.typeToNum(25.,error=':') == 'L5.0:'
    assert splat.typeToNum(25.,uncertainty=1.5) == 'L5.0:'
    assert splat.typeToNum(25.,lumclass='I') == 'L5.0I'
#    assert splat.typeToNum(25.,peculiar=True) == 'L5.0p'
    assert splat.typeToNum(25.,ageclass='beta',uncertainty=2.5) == 'L5.0beta::'
#    t = splat.typeToNum([12,18,26])
    # assert isinstance(t,list)
    # assert 'M8.0' in t

def test_date_conversion():
# properdate
    assert splat.properDate('20050523') == '2005-05-23'
    assert splat.properDate(20050523) == '2005-05-23'
    assert splat.properDate('050523') == '2005-05-23'
    assert splat.properDate(990523) == '1999-05-23'
    assert splat.properDate('5/23/2005') == '2005-05-23'
    assert splat.properDate('5/23/05') == '2005-05-23'
    assert splat.properDate('5/23/85') == '1985-05-23'
    assert splat.properDate('23 May 2005') == '2005-05-23'
    assert splat.properDate('2005 May 23') == '2005-05-23'
    assert splat.properDate('20050523',output='YYYYMMDD') == '20050523'
    assert splat.properDate('20050523',output='YYMMDD') == '050523'
    assert splat.properDate('20050523',output='MM/DD/YY') == '05/23/05'
    assert splat.properDate('20050523',output='MM/DD/YYYY') == '05/23/2005'
    assert splat.properDate('20050523',output='YYYY/MM/DD') == '2005/05/23'
    assert splat.properDate('20050523',output='DD/MM/YYYY') == '23/05/2005'
    assert splat.properDate('20050523',output='DD MMM YYYY') == '23 May 2005'
    assert splat.properDate('20050523',output='YYYY MMM DD') == '2005 May 23'


def test_stats():
# distributionStats 
# weightedMeanVar 
    values = [0.3,0.5,2.4,6.1,5.6,3.2,1.8,5.0,2.2,-4.5,10.2,-0.5,numpy.nan]
    weights = [1,1,2,3,3,2,2,3,2,0,4,0,0]
    x = splat.distributionStats(values)
    assert isinstance(x,numpy.ndarray)
    assert_allclose(x[1],numpy.nanmedian(values))
    x1 = splat.distributionStats(values,sigma=1.5)
    assert x[0] >= x1[0]
    assert x[-1] <= x1[-1]
    assert x[1] == x1[1]
    x2 = splat.distributionStats(values,q=[0.25,0.5,0.75])
    assert x[0] <= x2[0]
    assert x[-1] >= x2[-1]
    assert x[1] == x2[1]
    x3 = splat.distributionStats(values,weights=weights)
    assert x[0] <= x3[0]
    assert x[-1] <= x3[-1]
    assert x[1] <= x3[1]

    weights0 = numpy.ones(len(weights))
    x = splat.weightedMeanVar(values,weights0)
    assert isinstance(x,tuple)
    x2 = splat.weightedMeanVar(values,weights)
    assert x2[0] > x[0]
    assert x2[-1] < x[-1]
    

# incomplete tests
def test_check_keys():
# checkKeys
    pass

