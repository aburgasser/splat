# -*- coding: utf-8 -*-
from __future__ import print_function

# this is the test function set for splat core functions

# imports - internal
import copy
import glob
import os
# # imports - external
import numpy
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units
from astropy.io import fits
from numpy.testing import assert_allclose

# splat functions and constants
import splat
import splat.model as spmdl
from splat.initialize import *
from splat.utilities import directoryTree


#####################
# TESTING FUNCTIONS #
#####################

def test_read():
# readSpectrum
    file = os.path.join(directoryTree(SPECTRAL_DATA_FOLDER)[0],'10001_10443.fits')
    dstr = splat.readSpectrum(file)
    assert 'wave' in list(dstr.keys())
    assert 'flux' in list(dstr.keys())
    assert 'noise' in list(dstr.keys())
    assert 'header' in list(dstr.keys())
    assert len(dstr['wave']) > 0
    assert len(dstr['flux']) > 0
    assert len(dstr['noise']) > 0
    assert type(dstr['header']) == fits.header.Header

def test_spectrum_class():
# Spectrum class
# __init__
# __copy__
# copy
# info
# export
# save
    ofile = TMPFILENAME+'.fits'
    s = splat.Spectrum(10001)
    s2 = s.copy()
    assert_allclose(s2.flux[0].value,s.flux[0].value)
    assert isinstance(s.info(printout=False),str)
    if len(glob.glob(ofile)) > 0:
        os.remove(ofile)
    s.export(ofile)
    s2 = splat.Spectrum(filename=ofile)
    assert_allclose(s2.flux[0].value,s.flux[0].value)
    os.remove(ofile)
    s.save(ofile)
    s2 = splat.Spectrum(filename=ofile)
    assert_allclose(s2.flux[0].value,s.flux[0].value)
    os.remove(ofile)


def test_spectrum_trim():
# waverange
# trim
    s = splat.Spectrum(10001)
    assert_allclose(s.waveRange()[0].value,0.6450077295303345)
    un = s.waveRange()[0].unit
    s.trim([1.0,2.0])
    assert_allclose(s.waveRange()[0].value,1.0,rtol=1.e-2)
    assert_allclose(s.waveRange()[1].value,2.0,rtol=1.e-2)
    assert s.waveRange()[0].unit == un
    s.reset()
    assert_allclose(s.waveRange()[0].value,0.6450077295303345)
    assert s.waveRange()[0].unit == un

    
def test_spectrum_math():
# Spectrum class
# __add__
# __sub__
# __mul__
# __div__
    s1 = splat.Spectrum(10001)
    s1.normalize()
    assert_allclose(s1.fluxMax().value,1.)
    s2 = splat.Spectrum(10001)
    s2.normalize()
    assert_allclose(s2.fluxMax().value,1.)
    assert s1.flux.unit == s2.flux.unit

    s3 = s1+s2
    assert_allclose(s3.fluxMax().value,2.)
    assert s3.flux.unit == s1.flux.unit
#    assert_allclose(s3.snr,s1.snr*(2.**0.5))
    s4 = s3-s1
    assert_allclose(s4.fluxMax().value,1.)
    assert s4.flux.unit == s1.flux.unit
#    assert_allclose(s4.snr,s1.snr/(3.**0.5))
    s5 = s1*s2
    assert_allclose(s5.fluxMax().value,1.)
    assert s5.flux.unit == s1.flux.unit**2
#    assert s5.snr < s1.snr

#
# NOTE: THIS IS NOT WORKING AT THE MOMENT
#
    s6 = s3/s1
    assert_allclose(s6.fluxMax().value,2.)
    assert s6.flux.unit == s1.flux.unit/s1.flux.unit


    

def test_spectrum_scaling():
# normalize
# scale
# computeSN
# fluxCalibrate
# fluxMax
# reset
    s = splat.Spectrum(10001)
    omx = s.fluxMax().value
    snr0 = s.computeSN()

    s.normalize()
    assert_allclose(s.fluxMax().value,1.)
    assert_allclose(snr0,s.computeSN())
#    assert s.fscale == 'Normalized'
    s.scale(2.0)
    assert_allclose(s.fluxMax().value,2.)
    assert_allclose(snr0,s.computeSN())
#    assert s.fscale == 'Scaled'

    s.fluxCalibrate('2MASS J',15.0)
#    assert s.fscale == 'Apparent'
    s.fluxCalibrate('2MASS J',15.0,apparent=True)
#    assert s.fscale == 'Apparent'
    s.fluxCalibrate('2MASS J',15.0,absolute=True)
#    assert s.fscale == 'Absolute'
    mag,mag_e = splat.filterMag(s,'2MASS J')
    assert_allclose(mag,15.0)

    s.reset()
    assert_allclose(s.fluxMax().value,omx)
    assert_allclose(snr0,s.computeSN())

    
def test_spectrum_conversions():
# toFnu
# toFlam
# surface
# brightnessTemperature
# temperature
    s = splat.Spectrum(10001)
    bunit = s.flux.unit
    omx = s.fluxMax()
    assert s.flux.unit.is_equivalent(u.watt/u.m**3)
    s.flamToFnu()
    assert s.flux.unit.is_equivalent(u.Jy)
    s.fnuToFlam()
    assert s.flux.unit.is_equivalent(u.watt/u.m**3)

    s.fluxCalibrate('2MASS J',15.0,absolute=True)
#    assert s.fscale == 'Absolute'
    mx = s.fluxMax()
    s.toSurface(0.1)
#    assert s.fscale == 'Surface'
    assert_allclose(s.fluxMax()/mx,((100.*u.pc/const.R_sun)**2).to(u.m/u.m))

    s.toBrightnessTemperature()
#    assert s.fscale == 'Temperature'
    assert s.flux.unit == u.K

    s.reset()
    assert s.flux.unit == bunit
    assert s.fluxMax() == omx

    
def test_spectrum_plot():
# plot
    file = TMPFILENAME+'.png'
    s = splat.Spectrum(10001)
    s.plot(file=file)
    assert os.access(file,os.R_OK)
    os.remove(file)


# THIS TEST NEEDS TO BE RE-WORKED
# def test_spectrum_smooth():
# # smooth
# # _smoothToResolution
# # _smoothToSlitPixelWidth
# # _smoothToSlitWidth
#     s = splat.Spectrum(10001)
#     res0 = copy.deepcopy(s.resolution)
#     slit0 = copy.deepcopy(s.slitwidth)
#     slitp0 = copy.deepcopy(s.slitpixelwidth)
#     assert slit0.unit == u.arcsec

#     s.smooth(resolution = 50)
#     assert_allclose(s.resolution,50)
#     assert_allclose(s.slitwidth,slit0*res0/50)
#     assert_allclose(s.slitpixelwidth,slitp0*res0/50)
#     s.reset()

#     s.smooth(slitwidth = 3*u.arcsec)
#     assert_allclose(s.slitwidth,3.*u.arcsec)
#     assert_allclose(s.resolution,res0*(slit0/s.slitwidth).value)
#     assert_allclose(s.slitpixelwidth,slitp0*(s.slitwidth/slit0).value)
#     s.reset()

#     s.smooth(slitwidth = 3.)
#     assert_allclose(s.slitwidth,3.*u.arcsec)
#     assert_allclose(s.resolution,res0*(slit0/s.slitwidth).value)
#     assert_allclose(s.slitpixelwidth,slitp0*(s.slitwidth/slit0).value)
#     s.reset()

#     s.smooth(slitpixelwidth = 5)
#     assert_allclose(s.slitpixelwidth,5)
#     assert_allclose(s.resolution,res0*slitp0/s.slitpixelwidth)
#     assert_allclose(s.slitwidth,slit0*s.slitpixelwidth/slitp0)
#     s.reset()

#     assert_allclose(s.slitpixelwidth,slitp0)
#     assert_allclose(s.resolution,res0)
#     assert_allclose(s.slitwidth,slit0)

    
def test_classify_gravity():
# classifyByTemplate
# classifyGravity
    sp = splat.Spectrum(10001)  
    grav = splat.classifyGravity(sp)
    assert grav == 'FLD-G'
    grav = splat.classifyGravity(sp,output='allmeasures')
    assert isinstance(grav,dict)
    assert grav['gravity_class'] == 'FLD-G'
    #assert grav['spt'] == 'L4.0'
    grav = splat.classifyGravity(sp,output='allmeasures',spt='L5')
    assert grav['gravity_class'] == 'FLD-G'
    assert grav['spt'] == 'L5.0'

    sp = splat.getSpectrum(name='G 196-3B')[0]
    grav = splat.classifyGravity(sp,output='allmeasures')
    assert grav['gravity_class'] == 'VL-G'
    assert grav['spt'] == 'L3.0'
    assert grav['score'] >= 2.



def test_compare():
# compareSpectra
    sp = splat.Spectrum(10001)
    mdl = spmdl.getModel(teff=2000,logg=5.0)
    chi, scale = splat.compareSpectra(sp,mdl,mask_standard=True,stat='chisqr')
#    sys.stderr.write('\nScaling model: chi^2 = {}, scale = {}'.format(chi,scale))
#    sys.stderr.write('\n...compareSpectra successful\n'.format(chi,scale))
    
    
def test_indices():
# measureIndex
# measureIndexSet
# classifyByIndex
    sp = splat.Spectrum(10001)
    ind = splat.measureIndexSet(sp,set='burgasser')
    assert 'CH4-H' in list(ind.keys())
    assert_allclose(ind['CH4-H'][0],1.04,rtol=ind['CH4-H'][1])

    for st in list(INDEX_SETS.keys()):
        ind = splat.measureIndexSet(sp,set=st)
        assert len(ind) > 0
    
    spt, spt_e = splat.classifyByIndex(sp,'allers',string_flag=False)
    assert_allclose(spt,24.5,rtol=spt_e)
    spt, spt_e = splat.classifyByIndex(sp,'burgasser',string_flag=False)
    assert_allclose(spt,26.0,rtol=spt_e)
    spt, spt_e = splat.classifyByIndex(sp,'reid',string_flag=False)
    assert_allclose(spt,24.6,rtol=spt_e)



    
def test_standards():
# initiateStandards
# getStandard
# classifyByStandard
    std = splat.getStandard('M1.0')
    assert 'M1.0' in list(STDS_DWARF_SPEX.keys())
    assert std.name == 'Gl424'
    assert std.spex_type == 'M1.0'

    std = splat.getStandard('esdM5.0')
    assert 'esdM5.0' in list(STDS_ESD_SPEX.keys())
    assert std.name == 'LP 589-7'
    assert std.opt_type == 'esdM5'

    std = splat.getStandard(16,sd=True)
    assert 'sdM6.0' in list(STDS_SD_SPEX.keys())
    assert std.name == 'LHS 1074'
    assert std.opt_type == 'sdM6'

    splat.initiateStandards()
    for typ in list(STDS_DWARF_SPEX_KEYS.keys()):
        assert type(STDS_DWARF_SPEX[typ]) == splat.Spectrum
    splat.initiateStandards(sd=True)
    for typ in list(STDS_SD_SPEX_KEYS.keys()):
        assert type(STDS_SD_SPEX[typ]) == splat.Spectrum
    splat.initiateStandards(esd=True)
    for typ in list(STDS_ESD_SPEX_KEYS.keys()):
        assert type(STDS_ESD_SPEX[typ]) == splat.Spectrum

    sp = splat.Spectrum(10001)
    spt, spt_e = splat.classifyByStandard(sp,method='kirkpatrick',verbose=True)
    assert spt == 'L5.0'
    spt, spt_e = splat.classifyByStandard(sp,method='kirkpatrick',string=False)
    assert spt == 25.0
    spt, spt_e = splat.classifyByStandard(sp,sptrange=['M1.0','M8.0'],verbose=True)
    assert spt == 'M8.0'
    spt, spt_e = splat.classifyByStandard(sp,sd=True)
    assert spt == 'sdL5.0'
    spt, spt_e = splat.classifyByStandard(sp,esd=True)
    assert spt == 'esdM5.0'

    file = TMPFILENAME+'.png'
    spt, spt_e = splat.classifyByStandard(sp,method='kirkpatrick',plot=True,file=file)
    assert os.access(file,os.R_OK)
    os.remove(file)

    
def test_getspectrum():
# getSpecturm
# get random spectrum
    sp = splat.getSpectrum(lucky=True,published=True)
    assert type(sp) == list
    assert type(sp[0]) == splat.Spectrum
#    print(type(sp[0]))

# get spectrum by shortname
    sp = splat.getSpectrum(shortname='J0559-1404')
    assert sp[0].name == '2MASS J05591914-1404488'

# get spectrum by source key
    sp = splat.getSpectrum(source_key=10001)    
    assert int(sp[0].source_key) == 10001
    print(sp[0].name)
    assert sp[0].name == 'SDSS J000013.54+255418.6'

# get spectrum by data key
    sp = splat.getSpectrum(data_key=10001)
    assert int(sp[0].data_key) == 10001
    assert sp[0].name == 'SDSSp J053951.99-005902.0'

# get spectrum in a range of subtypes
    sp = splat.getSpectrum(spt=['L5','T5'],spt_type='spex',lucky=True)[0]
    n = splat.typeToNum(sp.spex_type)
    assert (n>=25 and n<=35)

# make sure getSpectrum fails properly
    sp = splat.getSpectrum(shortname='9999')
    assert len(sp) == 0


# incomplete tests
def test_database():
# fetchDatabase
# keySource
# keySpectrum
# searchLibrary
    s = splat.searchLibrary(published=True)
    assert len(s) > 1500
    s = splat.searchLibrary(young=True,output='DATA_FILE')
    assert '.fits' in s[0]


def test_masking():
# _generateMask
    pass

def test_classify_template():
# classifyByTemplate
    pass

def test_ews():
# measureEW
# measureEWSet
    pass

def test_metallicity():
# metallicity
    pass

