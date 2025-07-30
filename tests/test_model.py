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



#####################
# TESTING FUNCTIONS #
#####################


def OLD_test_loadmodel(model='burrows',teff=1000,logg=5.0,**kwargs):
    mdl = loadModel(model=model,teff=teff,logg=logg,**kwargs)
    if len(mdl.wave) > 0:
        mdl.info()
        mdl.scale(1.e-24)
        print(mdl.fluxMax())
#        mdl.plot()
    return True

def OLD_test_loadinterpolatedmodel(model='burrows',teff=1025,logg=4.75,**kwargs):
    mdl = loadInterpolatedModel(model=model,teff=teff,logg=logg,**kwargs)
    if len(mdl.wave) > 0:
        mdl.info()
        mdl.scale(1.e-24)
        print(mdl.fluxMax())
#        mdl.plot()
    return True

def OLD_test_modelfitgrid(shname='1507-1627',model='BTSettl2008',teff_range=[1000,2200],logg_range=[4.5,5.5],**kwargs):
#    tbl = splat.searchLibrary(spt=['M7','T8'])
#    sp = splat.Spectrum(numpy.random.choice(tbl['DATA_KEY']))
    from .splat import getSpectrum
    sp = getSpectrum(shortname=shname)[0]
    sp.fluxCalibrate('2MASS J',12.32,absolute=True)
    bp = modelFitGrid(sp,teff_range=teff_range,logg_range=logg_range,model=model,file=kwargs.get('folder','')+'test_modelfitgrid.pdf',**kwargs)
    print(bp)

    return

def OLD_test_modelfitEMCEE(folder):
#    tbl = splat.searchLibrary(spt=['M7','T8'])
#    sp = splat.Spectrum(numpy.random.choice(tbl['DATA_KEY']))
    folder='/Users/adam/projects/splat/code/testing/'
    from .splat import getSpectrum, classifyByStandard
    sp = getSpectrum(shortname='1507-1627')[0]
    sp.fluxCalibrate('2MASS J',12.32,absolute=True)
    spt,spt_e = classifyByStandard(sp,method='kirkpatrick')
    teff,teff_e = spemp.typeToTeff('L5')
    print('\nPerforming emcee model fit of {} with SpT = {} and initial Teff = {}\n'.format(sp.name,spt,teff))
# this takes about 1 hour
    return modelFitEMCEE(sp,t0=teff,g0=5.0,z0=0.,noprompt=True,use_weights=True,fit_metallicity=False,nwalkers=10,nsamples=100,output=folder+'test_modelfitEMCEE',verbose=True)

def OLD_test_modelfitMCMC(folder):
    from .splat import getSpectrum, classifyByStandard
    sp = getSpectrum(shortname='1047+2124')[0]        # T6.5 radio emitter
    spt,spt_e = classifyByStandard(sp,spt=['T2','T8'])
    teff,teff_e = spemp.typeToTeff(spt)
    sp.fluxCalibrate('MKO J',splat.typeToMag(spt,'MKO J')[0],absolute=True)
    return modelFitMCMC(sp, mask_standard=True, initial_guess=[teff, 5.3, 0.], zstep=0.1, nsamples=100,savestep=0,filebase=basefolder+'fit1047',verbose=True)

