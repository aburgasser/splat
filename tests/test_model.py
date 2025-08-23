# -*- coding: utf-8 -*-
from __future__ import print_function

# this is the test function set for splat model functions

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
import splat.model as spmdl
import splat.empirical as spemp
from splat.initialize import *


#####################
# TESTING FUNCTIONS #
#####################


def test_loadmodel():
# grid model
    mdl = spmdl.loadModel(modelset='burrows',teff=1000,logg=5.0)
    mdl.info()
    mdl.normalize()
    assert numpy.nanmax(mdl.flux.value)==1.0
# interpolated model
    mdl = spmdl.loadModel(modelset='burrows',teff=1025,logg=4.7)
    mdl.normalize()
    assert numpy.nanmax(mdl.flux.value)==1.0
    return 

def test_modelfitgrid():
    sp = splat.getSpectrum(shortname='J1507-1627')[0]
    sp.fluxCalibrate('2MASS J',12.32,absolute=True)
    res = spmdl.modelFitGrid(sp,teff_range=[1000,2200],logg_range=[4.5,5.5],modelset='btsettl',plot=False)
    assert res['logg']==5.5
    assert res['teff']==1800
    assert res['z']==0.0
    return

# THESE TESTS ARE FAILING, NEED TO REBUILD 

# def OLD_test_modelfitEMCEE(folder):
#     sp = splat.getSpectrum(shortname='1507-1627')[0]
#     sp.fluxCalibrate('2MASS J',12.32,absolute=True)
#     spt,spt_e = splat.classifyByStandard(sp,method='kirkpatrick')
#     teff,teff_e = spemp.typeToTeff('L5')
#     print('\nPerforming emcee model fit of {} with SpT = {} and initial Teff = {}\n'.format(sp.name,spt,teff))
# # this takes about 1 hour
# #    spmdl.modelFitEMCEE(sp,t0=teff,g0=5.0,z0=0.,noprompt=True,use_weights=True,fit_metallicity=False,nwalkers=10,nsamples=100,output=folder+'test_modelfitEMCEE',verbose=True)
#     return

# def OLD_test_modelfitMCMC(folder):
#     sp = splat.getSpectrum(shortname='1047+2124')[0]        # T6.5 radio emitter
#     spt,spt_e = splat.classifyByStandard(sp,spt=['T2','T8'])
#     teff,teff_e = spemp.typeToTeff(spt)
#     sp.fluxCalibrate('MKO J',splat.typeToMag(spt,'MKO J')[0],absolute=True)
#     spmdl.modelFitMCMC(sp, mask_standard=True, initial_guess=[teff, 5.3, 0.], zstep=0.1, nsamples=100,savestep=0,filebase=basefolder+'fit1047',verbose=True)
#     return

