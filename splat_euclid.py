from __future__ import print_function, division

"""
.. note::
         These are the database functions for SPLAT 
"""

import astropy
#import copy
#from datetime import datetime
import os
import re
import requests
#from splat import SPLAT_PATH, SPLAT_URL
import sys
#from scipy import stats
from astropy.io import ascii, fits            # for reading in spreadsheet
from astropy.table import Table, join            # for reading in table files
from astropy.coordinates import SkyCoord

import copy
import splat
import numpy
from astropy import units as u            # standard units
from scipy.interpolate import interp1d

# some euclid parameters
EUCLID_WAVERANGE = [1.25,1.85]
EUCLID_RESOLUTION = 250
EUCLID_NOISE = 3e-15*u.erg/u.s/u.cm**2/u.micron	
EUCLID_SAMPLING = 0.0013 # micron/pixel
EUCLID_SLITWIDTH = numpy.mean(EUCLID_WAVERANGE)/EUCLID_RESOLUTION/EUCLID_SAMPLING # micron/pixel

sys.ps1 = 'splat euclid> '

# program to convert spex spectra in to euclid spectra

def spexToEuclid(sp):
	sp.resolution = EUCLID_RESOLUTION
	sp.slitpixelwidth = EUCLID_SLITWIDTH
	f = interp1d(sp.wave,sp.flux)
	n = interp1d(sp.wave,sp.noise)
	sp.wave = numpy.arange(EUCLID_WAVERANGE[0],EUCLID_WAVERANGE[1],EUCLID_SAMPLING)*sp.wunit
	sp.flux = f(sp.wave/sp.wunit)*sp.funit
	sp.noise = n(sp.wave/sp.wunit)*sp.funit
# update other spectrum elements
	sp.snr = sp.computeSN()
	sp.flam = sp.flux
	sp.nu = sp.wave.to('Hz',equivalencies=u.spectral())
	sp.fnu = sp.flux.to('Jy',equivalencies=u.spectral_density(sp.wave))
	sp.variance = sp.noise**2
	sp.dof = numpy.round(len(sp.wave)/sp.slitpixelwidth)
	sp.history.append('Converted to EUCLID format')

# adds noise to spectrum
def addEuclidNoise(sp):
	sp2 = copy.deepcopy(sp)
	bnoise = numpy.zeros(len(sp2.noise))+EUCLID_NOISE
	bnoise.to(sp2.funit,equivalencies=u.spectral())
	anoise = numpy.random.normal(0,EUCLID_NOISE/sp2.funit,len(bnoise))*sp2.funit
	sp2.flux = sp2.flux+anoise
	sp2.variance = sp2.variance+bnoise**2
	sp2.noise = sp2.variance**0.5
	sp2.snr = sp2.computeSN()
	sp2.history.append('Added spectral noise based on Euclid sensitivity')
	return sp2


if __name__ == '__main__':
    ofold = '/Users/adam/projects/splat/euclid/'
    sp = splat.getSpectrum(shortname='J0559-1404')[0]
    spt = 'T4.5'
    filter = 'MKO H'
    m1 = splat.typeToMag(spt,filter)[0]
    m2 = 21
    d2 = splat.estimateDistance(sp,spt=spt,mag=m2, absmag=m1)[0]
    m3 = 19
    d3 = splat.estimateDistance(sp,spt=spt,mag=m3, absmag=m1)[0]

    spexToEuclid(sp)
    sp.normalize()
    sp.fluxCalibrate(filter,m2)
    print(sp.snr)
    sp2 = addEuclidNoise(sp)
    print(sp2.snr)
    sp.fluxCalibrate(filter,m3)
    sp3 = addEuclidNoise(sp)
    print(sp3.snr)
    sp.normalize()
    sp2.normalize()
    sp3.normalize()
    cls1 = splat.classifyByIndex(sp)
    cls2 = splat.classifyByIndex(sp2)
    cls3 = splat.classifyByIndex(sp3)
    splat.plotSpectrum(sp,sp3,sp2,colors=['k','b','r'],stack=0.5,xrange=EUCLID_WAVERANGE,\
    	yrange=[-0.3,2.2],output=ofold+'spectral_degradation.eps',legend=[\
    	'J0559-1404 MH = {:.1f}, SpT = {}+/-{:.1f}'.format(m1,cls1[0],cls1[1]),\
    	'H = {:.1f}, d = {:.1f} pc, SpT = {}+/-{:.1f}'.format(m3,d3,cls3[0],cls3[1]),\
    	'H = {:.1f}, d = {:.1f} pc, SpT = {}+/-{:.1f}'.format(m2,d2,cls2[0],cls2[1])])
    print('For H = {}, SpT = {}+/-{}'.format(m1,spt[0],spt[1]))
    print('For H = {}, SpT = {}+/-{}'.format(m3,spt[0],spt[1]))
    print('For H = {}, SpT = {}+/-{}'.format(m2,spt[0],spt[1]))

