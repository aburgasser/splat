from __future__ import print_function, division

"""
.. note::
         These are functions related to the EUCLID analysis based on SPLAT tools 
"""

#import astropy
#import copy
#from datetime import datetime
#import os
#import re
#import requests
#from splat import SPLAT_PATH, SPLAT_URL
#from scipy import stats
#from astropy.io import ascii, fits            # for reading in spreadsheet
#from astropy.table import Table, join            # for reading in table files
#from astropy.coordinates import SkyCoord

# imports: internal
import copy
import os
import time

# imports: external
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import pandas
from astropy import units as u            # standard units
from scipy.interpolate import interp1d
from scipy.integrate import trapz        # for numerical integration

# imports: splat
from .core import getSpectrum, classifyByIndex
import splat.empirical as spemp
import splat.evolve as spev
import splat.plot as splot

# some euclid parameters
EUCLID_WAVERANGE = [1.25,1.85]
EUCLID_RESOLUTION = 250
EUCLID_NOISE = 3e-15*u.erg/u.s/u.cm**2/u.micron	
EUCLID_SAMPLING = 0.0013 # micron/pixel
EUCLID_SLITWIDTH = numpy.mean(EUCLID_WAVERANGE)/EUCLID_RESOLUTION/EUCLID_SAMPLING # micron/pixel

# program to convert spex spectra in to euclid spectra

def spexToEuclid(sp):
	'''
    :Purpose: Convert a SpeX file into EUCLID form, using the resolution and wavelength coverage
    				defined from the Euclid Red Book (`Laurijs et al. 2011 <http://sci.esa.int/euclid/48983-euclid-definition-study-report-esa-sre-2011-12/>`_). This function changes the input Spectrum
    				objects, which can be restored by the Spectrum.reset() method.

    :param sp: Spectrum class object, which should contain wave, flux and noise array elements

    :Example:
       >>> import splat
       >>> sp = splat.getSpectrum(lucky=True)[0]   				# grab a random file
       >>> splat.spexToEuclid(sp)
       >>> min(sp.wave), max(sp.wave)
		(<Quantity 1.25 micron>, <Quantity 1.8493000000000364 micron>)
       >>> sp.history
		[``'Spectrum successfully loaded``',
		 ``'Converted to EUCLID format``']
       >>> sp.reset()
       >>> min(sp.wave), max(sp.wave)
		(<Quantity 0.6454827785491943 micron>, <Quantity 2.555659770965576 micron>)
	'''

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
	'''
    :Purpose: Adds Gaussian noise to a EUCLID-formatted spectrum assuming a constant noise
    			model of 3e-15 erg/s/cm2/micron (as extrapolated from the Euclid Red Book; 
    			Laurijs et al. 2011 <http://sci.esa.int/euclid/48983-euclid-definition-study-report-esa-sre-2011-12/>`_).
    			Note that noise is added to both flux and (in quadrature) variance. This function creates a
    			new Spectrum object so as not to corrupt the original data. 

    :param sp: Spectrum class object, which should contain wave, flux and noise array elements

    :Output: Spectrum object with Euclid noise added in

    :Example:
       >>> import splat
       >>> sp = splat.getSpectrum(lucky=True)[0]   				# grab a random file
       >>> splat.spexToEuclid(sp)
       >>> sp.normalize()
       >>> sp.scale(1.e-14)
       >>> sp.computeSN()
       		115.96374031163553
       >>> sp_noisy = splat.addEculidNoise(sp)
       >>> sp_noisy.computeSN()
       		3.0847209519763172
	'''

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
	'''
	Test function for splat_euclid functions, taking an 0559-1404 spectrum and plotting a 3 
	apparent magnitudes with corresponding distances
	'''

	ofold = '/Users/adam/projects/splat/euclid/'
	sp = getSpectrum(shortname='J0559-1404')[0]
	spt = 'T4.5'
	filter = 'MKO H'
	m1 = spemp.typeToMag(spt,filter)[0]
	m2 = 21
	d2 = spemp.estimateDistance(sp,spt=spt,mag=m2, absmag=m1)[0]
	m3 = 19
	d3 = spemp.estimateDistance(sp,spt=spt,mag=m3, absmag=m1)[0]

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
	cls1 = classifyByIndex(sp)
	cls2 = classifyByIndex(sp2)
	cls3 = classifyByIndex(sp3)
	splot.plotSpectrum(sp,sp3,sp2,colors=['k','b','r'],stack=0.5,xrange=EUCLID_WAVERANGE,\
		yrange=[-0.3,2.2],output=ofold+'spectral_degradation.eps',legend=[\
		'J0559-1404 MH = {:.1f}, SpT = {}+/-{:.1f}'.format(m1,cls1[0],cls1[1]),\
		'H = {:.1f}, d = {:.1f} pc, SpT = {}+/-{:.1f}'.format(m3,d3,cls3[0],cls3[1]),\
		'H = {:.1f}, d = {:.1f} pc, SpT = {}+/-{:.1f}'.format(m2,d2,cls2[0],cls2[1])])
	print('For H = {}, SpT = {}+/-{}'.format(m1,spt[0],spt[1]))
	print('For H = {}, SpT = {}+/-{}'.format(m3,spt[0],spt[1]))
	print('For H = {}, SpT = {}+/-{}'.format(m2,spt[0],spt[1]))



