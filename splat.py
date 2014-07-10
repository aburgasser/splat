# WORKING COPY OF SPLAT CODE LIBRARY
# based on routines developed by:
#	Christian Aganze
#	Daniella Bardalez Gagliuffi
# 	Adam Burgasser
#	Caleb Choban
#	Ivanna Escala
#	Aishwarya Iyer
# 	Yuhui Jin
#	Michael Lopez
#	Alex Mendez
#	Julian Pilate-Hutcherson
#	Maitrayee Sahi
#	Melisa Tallis

#
# CURRENT STATUS (5/6/2014)
# can now load up spectra and models from online sources
# source spectra can be selected by name, designation, young/subdwarf/red/blue/binary/spbin
# returned as array of spectra
#
# There is an odd error that is coming in when reading in multiple spectra:
#	/Users/adam/projects/splat/exercises/ex4/splat.py:746: RuntimeWarning: invalid value encountered in greater
#	w = numpy.where(flux > numpy.median(flux))
# only occurs on first go
#

# imports
import sys
import os
import copy
import numpy
import scipy
import astropy
import matplotlib.pyplot as plt
import re
import urllib2
import string
from scipy import stats, signal
from scipy.integrate import trapz		# for numerical integration
from scipy.interpolate import interp1d, griddata
from astropy.io import ascii, fits			# for reading in spreadsheet
from astropy.table import Table, join			# for reading in table files
from astropy.coordinates import ICRS, Galactic		# coordinate conversion
from astropy import units as u			# standard units
from astropy import constants as const		# physical constants in SI units

numpy.seterr(divide='warn')		# this is probably not an entirely safe approach

############ PARAMETERS ############
SPLAT_URL = 'http://pono.ucsd.edu/~adam/splat/'
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
parameter_names = ['teff','logg','z','fsed','kzz']
spex_pixel_scale = 0.15		# spatial scale in arcseconds per pixel
####################################

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
# 		if not quiet:
#				self.info()
#		  if not noplot:
#				self.plot(**kwargs)
		  return tmp
	 return _show

def Copy(fn):
	 def _copy(self, *args, **kwargs):
		  out = copy.copy(self)
		  return fn(out, *args, **kwargs)
	 return _copy


# define the Spectrum class which contains the relevant information
class Spectrum(object):
	@Show
	def __init__(self, **kwargs):
		'''Load the file'''
		self.model = False
		self.filename = kwargs.get('filename','')
		self.wlabel = kwargs.get('wlabel','Wavelength')
		self.wunit = kwargs.get('wunit',u.micron)
		self.flabel = kwargs.get('flabel','F_lambda')
		self.fscale = kwargs.get('fscale','')
		self.funit = kwargs.get('funit',u.erg/(u.cm**2 * u.s * u.micron))
		self.resolution = kwargs.get('resolution',250)	# default placeholder
		self.slitpixelwidth = kwargs.get('slitwidth',2)		# default placeholder
		self.slitwidth = self.slitpixelwidth*spex_pixel_scale
		self.simplefilename = os.path.basename(self.filename)
# wave and flux given
		if len(kwargs.get('wave','')) > 0 and len(kwargs.get('flux','')) > 0:
			self.wave = kwargs['wave']
			self.flux = kwargs['flux']
			if len(kwargs.get('noise','')) > 0:
				self.noise = kwargs['noise']
			else:
				self.noise = [numpy.nan for i in self.wave]
		else:
# filename given
			try:
				self.wave, self.flux, self.noise = readSpectrum(**kwargs)
			except:
				raise NameError('\nCould not load up spectral file')
		self.nu = const.c.to('micron/s').value/self.wave
# calculate variance
		self.variance = self.noise
		if (self.variance[0] != numpy.nan):
			self.variance = [n**2 for n in self.noise]

# preserve original values
		self.wave_original = copy.deepcopy(self.wave)
		self.flux_original = copy.deepcopy(self.flux)
		self.noise_original = copy.deepcopy(self.noise)
		self.variance_original = copy.deepcopy(self.variance)
		self.resolution = copy.deepcopy(self.resolution)
		self.slitpixelwidth = copy.deepcopy(self.slitpixelwidth)
		
# information on source spectrum
		if not (kwargs.get('model',False)):
			x,y = filenameToNameDate(self.filename)
			self.name = kwargs.get('name',x)
			self.date = kwargs.get('date',y)
			try:
				self.caldate = dateToCaldate(self.date)
			except:
				self.caldate = ''
		else:
# information on model
			self.model = True
			self.teff = kwargs.get('teff',numpy.nan)
			self.logg = kwargs.get('logg',numpy.nan)
			self.z = kwargs.get('z',numpy.nan)
			self.cloud = kwargs.get('cloud',numpy.nan)
			self.modelset = kwargs.get('set','')
			self.name = self.modelset+' Teff='+str(self.teff)+' logg='+str(self.logg)
			self.fscale = 'Surface'
		self.history = ['Loaded']

				
	def __repr__(self):
		'''A simple representation of an object is to just give it a name'''
		return 'Spectra Object for {}'.format(self.name)

	def __add__(self,other):
		'''Adding two spectra '''
		sp = copy.deepcopy(self)
		f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
		n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=0.)
		sp.flux = numpy.add(self.flux,f(self.wave))
		sp.variance = sp.variance+n(self.wave)
		sp.noise = [n**0.5 for n in sp.variance]
		sp.flux_original=sp.flux
		sp.noise_original=sp.noise
		sp.variance_original=sp.variance
		return sp

	def __sub__(self,other):
		'''Subtracting two spectra '''
		sp = copy.deepcopy(self)
		f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
		n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=0.)
		sp.flux = numby.subtract(self.flux,f(self.wave))
		sp.variance = sp.variance+n(self.wave)
		sp.noise = [n**0.5 for n in sp.variance]
		sp.flux_original=sp.flux
		sp.noise_original=sp.noise
		sp.variance_original=sp.variance
		return sp

	def __mul__(self,other):
		'''Multiplying two spectra'''
		sp = copy.deepcopy(self)
		f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
		n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=0.)
		sp.flux = numpy.multiply(self.flux,f(self.wave))
		sp.variance = numpy.multiply(numpy.power(sp.flux,2),(\
			numpy.divide(self.variance,numpy.power(sp.flux,2))+\
			numpy.divide(n,numpy.power(f(self.wave),2))))
		sp.noise = [n**0.5 for n in sp.variance]
		sp.flux_original=sp.flux
		sp.noise_original=sp.noise
		sp.variance_original=sp.variance
		return sp

	def __div__(self,other):
		'''Dividing two spectra'''
		sp = copy.deepcopy(self)
		f = interp1d(other.wave,other.flux,bounds_error=False,fill_value=0.)
		n = interp1d(other.wave,other.variance,bounds_error=False,fill_value=0.)
		sp.flux = numpy.divide(self.flux,f(self.wave))
		sp.variance = numpy.multiply(numpy.power(sp.flux,2),(\
			numpy.divide(self.variance,numpy.power(sp.flux,2))+\
			numpy.divide(n,numpy.power(f(self.wave),2))))
		sp.noise = [n**0.5 for n in sp.variance]
		sp.flux_original=sp.flux
		sp.noise_original=sp.noise
		sp.variance_original=sp.variance
		return sp

	def copy(self):
		  '''Make a copy of the current spectrum'''
		  other = copy.deepcopy(self)
		  return other

	def info(self):
		  '''Report some information about this spectrum'''
		  if (self.model):
		  	print '''{0} model with Teff = {1} and log g = {2}'''.format(self.modelset, self.teff, self.logg)
		  else:
		  	print '''Spectrum of {0} taken on {1}'''.format(self.name, self.date)
		  return

	def flamToFnu(self):
	 	'''Convert flux density from F_lam to F_nu, the later in Jy'''
	 	self.funit = u.Jy
	 	self.flabel = 'F_nu'
	 	self.flux.to(funit,equivalencies=u.spectral_density(self.wave))
	 	self.noise.to(funit,equivalencies=u.spectral_density(self.wave))
	 	return

	def fluxCalibrate(self,filter,mag,**kwargs):
	 	'''Calibrate spectrum to input magnitude'''
		absolute = kwargs.get('absolute',False)
		apparent = kwargs.get('apparent',False)
	 	self.normalize()
	 	apmag = filterMag(self,filter,**kwargs)
	 	if (~numpy.isnan(apmag)):
	 		self.scale(10.**(0.4*(apmag-mag)))
			if (absolute):
				self.fscale = 'Absolute'
			if (apparent):
				self.fscale = 'Apparent'
		return

	def fluxMax(self):
		return numpy.nanmax(self.flux[numpy.where(\
			numpy.logical_and(self.wave > 0.8,self.wave < 2.3))])

	def fnuToFlam(self):
	 	'''Convert flux density from F_nu to F_lam, the later in erg/s/cm2/Hz'''
	 	self.funit = u.erg/(u.cm**2 * u.s * u.micron)
	 	self.flabel = 'F_lam'
	 	self.flux.to(funit,equivalencies=u.spectral_density(self.wave))
	 	self.noise.to(funit,equivalencies=u.spectral_density(self.wave))
	 	self.variance = [n**2 for n in self.noise]
	 	return

	def normalize(self):
	 	'''Normalize spectrum'''
		self.scale(1./self.fluxMax())
		self.fscale = 'Normalized'
		return

	def reset(self):
	 	'''Reset to original spectrum'''
		self.wave = copy.deepcopy(self.wave_original)
		self.flux = copy.deepcopy(self.flux_original)
		self.noise = copy.deepcopy(self.noise_original)
		self.variance = copy.deepcopy(self.variance_original)
		self.resolution = copy.deepcopy(self.resolution_original)
		self.slitpixelwidth = copy.deepcopy(self.slitpixelwidth_original)
		self.slitwidth = self.slitpixelwidth*spex_pixel_scale
		self.fscale = ''
		return

	def scale(self,factor):
	 	'''Scale spectrum and noise by a constant factor'''
	 	self.flux = self.flux*factor
	 	self.noise = [n*factor for n in self.noise]
	 	self.variance = [n**2 for n in self.noise]
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
	 	'''Smooth spectrum to a constant resolution'''
	 	overscale = kwargs.get('overscale',10.)
		method = kwargs.get('method','hanning')
		kwargs['method'] = method

# do nothing if requested resolution is higher than current resolution
	 	if (resolution < self.resolution):
# sample onto a constant resolution grid at 5x current resolution
		 	r = resolution*overscale
		 	waveRng = self.waveRange()
		 	npix = numpy.floor(numpy.log(waveRng[1]/waveRng[0])/numpy.log(1.+1./r))
		 	wave_sample = [waveRng[0]*(1.+1./r)**i for i in numpy.arange(npix)]
		 	f = interp1d(self.wave,self.flux,bounds_error=False,fill_value=0.)
		 	v = interp1d(self.wave,self.variance,bounds_error=False,fill_value=0.)
		 	flx_sample = f(wave_sample)
		 	var_sample = v(wave_sample)
# now convolve a function to smooth resampled spectrum
			window = signal.get_window(method,numpy.round(overscale))
			neff = numpy.sum(window)/numpy.nanmax(window)		# effective number of pixels
			flx_smooth = signal.convolve(flx_sample, window/numpy.sum(window), mode='same')
			var_smooth = signal.convolve(var_sample, window/numpy.sum(window), mode='same')/neff
# resample back to original wavelength grid
			f = interp1d(wave_sample,flx_smooth,bounds_error=False,fill_value=0.)
			v = interp1d(wave_sample,var_smooth,bounds_error=False,fill_value=0.)
			self.flux = f(self.wave)
			self.variance = v(self.wave)
			self.noise = [n**0.5 for n in self.variance]
			self.slitpixelwidth = self.slitpixelwidth*self.resolution/resolution
			self.resolution = resolution
			self.slitwidth = self.slitpixelwidth*spex_pixel_scale
			self.history = ['Smoothed to constant resolution {}'.format(self.resolution)]		
	 	return

	def smoothToSlitPixelWidth(self,width,**kwargs):
	 	'''Smooth spectrum to a constant slit width (smooth by pixels)'''
		method = kwargs.get('method','hanning')
		kwargs['method'] = method
# do nothing if requested resolution is higher than current resolution
	 	if (width > self.slitpixelwidth):
# convolve a function to smooth spectrum
			window = signal.get_window(method,numpy.round(width))
			neff = numpy.sum(window)/numpy.nanmax(window)		# effective number of pixels
			self.flux = signal.convolve(self.flux, window/numpy.sum(window), mode='same')
			self.variance = signal.convolve(self.variance, window/numpy.sum(window), mode='same')/neff
			self.noise = [n**0.5 for n in self.variance]
			self.resolution = self.resolution*self.slitpixelwidth/width
			self.slitpixelwidth = width
			self.slitwidth = self.slitpixelwidth*spex_pixel_scale
			self.history = ['Smoothed to slit width of {}'.format(self.slitwidth)]		
	 	return

	def smoothToSlitWidth(self,width,**kwargs):
		method = kwargs.get('method','hanning')
		kwargs['method'] = method
	 	'''Smooth spectrum to a constant slit width (smooth by pixels)'''
		pwidth = width/spex_pixel_scale
		self.smoothToSlitPixelWidth(pwidth,**kwargs)
	 	return

	def snr(self):
	 	'''Compute a representative S/N value'''
	 	pass
	 	return

	def surface(self,radius):
	 	'''Convert to surface fluxes given a radius, assuming at absolute fluxes'''
	 	pass
	 	return

	def waveRange(self):
		ii = numpy.where(self.flux > 0)
		return [numpy.nanmin(self.wave[ii]), numpy.nanmax(self.wave[ii])]
	 
		
							 

# FUNCTIONS FOR SPLAT
def caldateToDate(d):
	'''Convert from numeric date to calendar date'''
	return d[:4]+str((months.index(d[5:8])+1)/100.)[2:4]+d[-2:]


def checkFile(filename,**kwargs):
	'''Check if a particular file is present in the online database'''
	url = kwargs.get('url',SPLAT_URL)+'/Spectra/'
	flag = checkOnline()
	if (flag):
		try:
			open(os.path.basename(filename), 'wb').write(urllib2.urlopen(url+filename).read())
		except urllib2.URLError, ex:
			flag = False
	return flag


def checkOnline():
	'''Check if you are online'''
	try:
		urllib2.urlopen(SPLAT_URL)
		return True
	except urllib2.URLError, ex:
		return False



def classifyByIndex(sp, *args, **kwargs):
	'''Classify a spectrum based on its spectral indices'''
	
	str_flag = kwargs.get('string', False)
	rnd_flag = kwargs.get('round', False)
	rem_flag = kwargs.get('remeasure', True)
	nsamples = kwargs.get('nsamples', 100)
	nloop = kwargs.get('nloop', 5)
	set = kwargs.get('set','burgasser')
	allowed_sets = ['burgasser','reid','testi','allers']

# measure indices if necessary
	if (len(args) != 0):
		indices = args[0]

# Burgasser (2007, ApJ, 659, 655) calibration
	if (set.lower() == 'burgasser'):
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

# Allers et al. (2013, ApJ, 657, 511)
	elif (set.lower() == 'allers'):
		if (rem_flag or len(args) == 0):
			kwargs['set'] = 'mclean'
			i1 = measureIndexSet(sp, **kwargs)
			kwargs['set'] = 'slesnick'
			i2 = measureIndexSet(sp, **kwargs)
			kwargs['set'] = 'allers'
			i3 = measureIndexSet(sp, **kwargs)
			indices = dict(i1.items() + i2.items() + i3.items())
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

	else:
		sys.stderr.write('\nWarning: '+set.lower()+' SpT-index relation not in measureSpT code\n\n')
		return numpy.nan, numpy.nan


	for index in coeffs.keys():
		vals = numpy.polyval(coeffs[index]['coeff'],numpy.random.normal(indices[index][0],indices[index][1],nsamples))*sptfact
		coeffs[index]['spt'] = numpy.nanmean(vals)+sptoffset
		coeffs[index]['sptunc'] = (numpy.nanstd(vals)**2+coeffs[index]['fitunc']**2)**0.5
		print index, indices[index], numpy.nanmean(vals), numpy.nanstd(vals), coeffs[index]['spt']
	
#	print indices[index][0], numpy.polyval(coeffs[index]['coeff'],indices[index][0]), coeffs[index]
	mask = numpy.ones(len(coeffs.keys()))
	result = numpy.zeros(2)
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
			coeffs[index]['mask'] = numpy.where( \
				coeffs[index]['range'][0] <= sptn <= coeffs[index]['range'][1],1,0)

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
	'''Classify a spectrum by comparing to spectral standards'''
	return numpy.nan, numpy.nan
	

def classifyByTemplate(sp, *args, **kwargs):
	'''Classify a spectrum by comparing to spectral templates'''
	return numpy.nan, numpy.nan
	
	
def compareSpectra(sp1, sp2, *args, **kwargs):
	'''Compare two spectra against each other'''
	return numpy.nan, numpy.nan
	
	

def coordinateToDesignation(c):
	'''Convert RA, Dec into designation string'''
# input is ICRS
	if isinstance(c,ICRS):
		return string.replace('J{0}{1}'.format(c.ra.to_string(u.hour, sep='', precision=2, pad=True), \
			c.dec.to_string(u.degree, sep='', precision=1, alwayssign=True, pad=True)),'.','')
# input is [RA,Decl] pair in degrees
	elif isinstance(c,list):
		cc = ICRS(ra=c[0],dec=c[1],unit=(u.deg,u.deg))
		return string.replace('J{0}{1}'.format(cc.ra.to_string(u.hour, sep='', precision=2, pad=True), \
			cc.dec.to_string(u.degree, sep='', precision=1, alwayssign=True, pad=True)),'.','')
# input is sexigessimal string
	elif isinstance(c,str):
		cc = ICRS(c)
		return string.replace('J{0}{1}'.format(cc.ra.to_string(u.hour, sep='', precision=2, pad=True), \
			cc.dec.to_string(u.degree, sep='', precision=1, alwayssign=True, pad=True)),'.','')
	else:
		raise ValueError('\nCould not parse input format\n\n')


def dateToCaldate(d):
	'''Convert from numeric date to calendar date'''
	return d[:4]+' '+months[int(d[5:6])-1]+' '+d[-2:]


def designationToCoordinate(value, **kwargs):
	'''Convert a designation into a RA, Dec tuple or ICRS'''
	icrsflag = kwargs.get('ICRS',True)

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
		return ICRS(ra=ra, dec=dec, unit=(u.degree, u.degree))
	else:
		return [ra,dec]


def designationToShortName(value):
	'''Produce a shortened version of designation'''
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

		

def fetchDatabase(*args, **kwargs):	
	'''Get the SpeX Database from either online repository or local drive'''
	dataFile = kwargs.get('dataFile','db_spexprism.txt')
	folder = kwargs.get('folder','~/')
	url = kwargs.get('url',SPLAT_URL)
	local = kwargs.get('local',False)

# check if online
	local = local or (not checkOnline())
		
# first try online
	if not local:
		try:
			open(os.path.basename(dataFile), 'wb').write(urllib2.urlopen(url+dataFile).read())
			data = ascii.read(dataFile, delimiter='	',fill_values='-99.')
			os.remove(os.path.basename(dataFile))
#			return data
		except urllib2.URLError, ex:
			sys.stderr.write('\nReading local '+dataFile+'\n\n')
			local = True
	
# now try local drive - NOT WORKING
	file = dataFile
	if local:
		if (os.path.exists(file) == False):
			file = folder+os.path.basename(file)
		if (os.path.exists(file) == False):
			raise NameError('\nCould not find '+dataFile+' or '+file+' locally\n\n')
		else:
			data = ascii.read(file, delimiter='	')

# TEMPORARY ADD-ONS UNTIL DATABASE IS COMPLETED
# add in RA/Dec (TEMPORARY)
	ra = []
	dec = []
	for x in data['designation']:
		c = designationToCoordinate(x,ICRS=False)
		ra.append(c[0])
		dec.append(c[1])
	data['ra'] = ra
	data['dec'] = dec

# add in young, subdwarf, binary, sbinary categories (TEMPORARY)
	data['young'] = ['young' in x for x in data['library']]
	data['subdwarf'] = ['subdwarf' in x for x in data['library']]
	data['binary'] = ['binary' in x for x in data['library']]
	data['spbin'] = ['spbin' in x for x in data['library']]
	data['blue'] = ['blue' in x for x in data['library']]
	data['red'] = ['red' in x for x in data['library']]

# add in shortnames (TEMPORARY)
	data['shortname'] = [designationToShortName(x) for x in data['designation']]

	return data


def filenameToNameDate(filename):
	'''Extract from a SPLAT filename the source name and observation date'''
	ind = filename.rfind('.')
	base = filename[:ind]
	spl = base.split('_')
	if (len(spl) < 2):
		return '', ''
	else:
		name = spl[-2]
		d = spl[-1]
		try:
			float(d)
			date = '20'+d
		except ValueError:
			print filename+' does not contain a date'
			date = d
		
		return name, date



def filterMag(sp,filter,*args,**kwargs):
	'''Compute the magnitudes of a spectrum given a filter name'''
# keyword parameters
	filterFolder = kwargs.get('filterFolder',SPLAT_URL+'Filters/')
	vegaFile = kwargs.get('vegaFile','vega_kurucz.txt')
	info = kwargs.get('info',False)
	units = kwargs.get('units',sp.funit)
	custom = kwargs.get('custom',False)
	vega = kwargs.get('vega',True)
	ab = kwargs.get('ab',False)
	photons = kwargs.get('photons',False)
	energy = kwargs.get('energy',False)

# filter file assignments
	filters = { \
		'2MASS J': {'file': 'j_2mass.txt', 'description': '2MASS J-band'}, \
		'2MASS H': {'file': 'h_2mass.txt', 'description': '2MASS H-band'}, \
		'2MASS Ks': {'file': 'ks_2mass.txt', 'description': '2MASS Ks-band'}, \
		'MKO J': {'file': 'j_atm_mko.txt', 'description': 'MKO J-band + atmosphere'}, \
		'MKO H': {'file': 'h_atm_mko.txt', 'description': 'MKO H-band + atmosphere'}, \
		'MKO K': {'file': 'k_atm_mko.txt', 'description': 'MKO K-band + atmosphere'}, \
		'MKO Ks': {'file': 'k_atm_mko.txt', 'description': 'MKO Ks-band'} \
		}

# check that requested filter is in list
	if (filter not in filters.keys()):
		print 'Filter '+filter+' not included in filterMag'
		info = True
		
# print out what's available
	if (info):
		print 'Filter names:'
		for x in filters.keys():
			print x+': '+filters[x]['description']	
		return numpy.nan

# convert units to erg/cm2/s/um if needed - TO BE DONE

# Read in filter
	if (custom == False):
		fwave,ftrans = numpy.genfromtxt(filterFolder+filters[filter]['file'], comments='#', unpack=True, \
			missing_values = ('NaN','nan'), filling_values = (numpy.nan))
	else:
		fwave,ftrans = custom[0],custom[1]
	fnu = const.c.to('micron/s').value/fwave
	fwave = fwave[~numpy.isnan(ftrans)]			# temporary fix
	ftrans = ftrans[~numpy.isnan(ftrans)]
				
# interpolate spectrum and vega spectrum onto filter wavelength function
	d = interp1d(sp.wave,sp.flux)
	result = numpy.nan
	
# compute relevant quantity
	if (vega):
# Read in Vega spectrum
		vwave,vtrans = numpy.genfromtxt(filterFolder+vegaFile, comments='#', unpack=True, \
			missing_values = ('NaN','nan'), filling_values = (numpy.nan))
		vwave = vwave[~numpy.isnan(vtrans)]			# temporary fix
		vtrans = vtrans[~numpy.isnan(vtrans)]
		v = interp1d(vwave,vtrans)
		result = -2.5*numpy.log10(trapz(ftrans*d(fwave),fwave)/trapz(ftrans*v(fwave),fwave))
	if (energy):
		result = trapz(ftrans*d(fwave),fwave)
	if (photons):
		convert = const.h.to('erg s')*const.c.to('micron/s')
		result = trapz(ftrans*d(fwave)*fwave,fwave) / convert.value
	if (ab):
# NOTE: THIS IS CURRENTLY NOT WORKING
		sp.flamToFnu()
		dd = interp1d(sp.nu,sp.flux)
		a = trapz(ftrans*dd(fnu),numpy.log10(fnu))
		b = trapz(ftrans,numpy.log10(fnu))
		result = -2.5*numpy.log10(a/b)-48.6
		sp.fnuToFlam()
	
	return result



def getSourceKey(*args, **kwargs):
	'''Search the SpeX database to extract the key reference for that Source'''

	pass
	return


def getSpectrum(*args, **kwargs):
	'''Get specific spectrs from online library'''

	result = []
	kwargs['output'] = 'data_file'
	files = searchLibrary(*args, **kwargs)
	if len(files) > 0:
		for x in files:
			result.append(loadSpectrum(x))
	else:
		sys.stderr.write('\nNo files match search criteria\n\n')
	return result
		


# simple number checker
def isNumber(s):
	'''check if something is a number'''
	try:
		float(s)
		return True
	except ValueError:
		return False


def loadInterpolatedModel(*args,**kwargs):
# interpolates within a 2x2 grid of teff and logg ONLY
	set = kwargs.get('set','btsettl08')
	kwargs['set'] = set
	teff = kwargs.get('teff',1000)
	kwargs['teff'] = teff
	logg = kwargs.get('logg',5.0)
	kwargs['logg'] = logg
	url = kwargs.get('folder',SPLAT_URL+'/Models/')
#	local = kwargs.get('local',False)
	kwargs['model'] = True

# first get model parameters
	param = loadModelParameters(set)
#	print param
	
# check if input parameters are within range
#	for t in param.colnames:
#		flg = kwargs.get(t,False)
#		if flg != False:
#			if (flg < param[t][0] or flg > param[t][1]):
#				raise ValueError('\n\nValue {}={} outside model range'.format(t,str(flg)))

# try simply loading model (on grid)
#	try:
#		return loadModel(**kwargs)
#	except:
#		sys.stderr.write('\nCalculating interpolated model\n')

# identify grid points around input parameters
# NOTE: THIS IS JUST FOR TEFF AND LOGG GRIDDING FOR NOW

	t = teff - teff%param['teff'][2]
	trng = [t,t+param['teff'][2]]
	g = logg - logg%param['logg'][2]
	grng = [g,g+param['logg'][2]]
	x,y = numpy.meshgrid(trng,grng)

	mkwargs = kwargs.copy()
	try:
		mkwargs['teff'] = trng[0]
		mkwargs['logg'] = grng[0]
		md11 = loadModel(**mkwargs)
		mkwargs['teff'] = trng[1]
		md21 = loadModel(**mkwargs)
		mkwargs['teff'] = trng[0]
		mkwargs['logg'] = grng[1]
		md12 = loadModel(**mkwargs)
		mkwargs['teff'] = trng[1]
		md22 = loadModel(**mkwargs)
	except:
		raise NameError('\nProblem loading grid of models in loadInterpolatedModel\n\n')

	mflx = numpy.zeros(len(md11.wave))
	val = numpy.array([[None]*2]*2)
	for i,w in enumerate(md11.wave):
		val = numpy.array([ \
			[numpy.log10(md11.flux[i]),numpy.log10(md21.flux[i])], \
			[numpy.log10(md12.flux[i]),numpy.log10(md22.flux[i])]])
		mflx[i] = 10.**(griddata((x.flatten(),y.flatten()),val.flatten(),(teff,logg),'linear'))
	
	return Spectrum(wave=md11.wave,flux=mflx,**kwargs)


def loadModel(*args, **kwargs):
	'''load up a model spectrum based on parameters'''
# keyword parameters
	set = kwargs.get('set','btsettl08')
	teff = kwargs.get('teff',1000)
	logg = kwargs.get('logg',5.0)
	z = kwargs.get('z',0.0)
	kzz = kwargs.get('kzz',0)
	fsed = kwargs.get('fsed',0)
	folder = kwargs.get('folder','')
	url = kwargs.get('folder',SPLAT_URL+'/Models/')
	local = kwargs.get('local',False)
	kwargs['model'] = True
	fileFlag = False
	
# check if online
	local = local or (not checkOnline())

# a filename has been passed - simply read this file
	if (len(args) > 0):
		kwargs['filename'] = args[0]
		fileFlag = True

# determine model set
	else:

# get model parameters
		param = loadModelParameters(set)

		if (set == 'btsettl08'):
			url = url+'/'+set+'/'
			kwargs['filename'] = set+'_'+str(teff)+'_'+str(logg)+'_-0.0_nc_nc_eq_0.5.txt'
#			kwargs['filename'] = 'lte'+'{:5.3f}'.format(teff/100000.)[2:]+'-'+str(logg)[0:3]+'-0.0.BT-Settl.7_r120.txt'		
		else: 
			raise NameError('\nCurrently only have BTSettl models\n\n')

		if (teff < param['teff'][0] or teff > param['teff'][1]):
			raise NameError('\nInput Teff out of range\n\n')

		if (logg < param['logg'][0] or logg > param['logg'][1]):
			raise NameError('\nInput logg out of range\n\n')

	
# if model is on grid read in that model
	if (teff in numpy.arange(param['teff'][0],param['teff'][1],param['teff'][2]) and \
		logg in numpy.arange(param['logg'][0],param['logg'][1],param['logg'][2])):
		fileFlag = True

# simple read in a file
	if fileFlag:
	
# first try online
		if not local:
			try:
				open(os.path.basename(kwargs['filename']), 'wb').write(urllib2.urlopen(url+kwargs['filename']).read())
				kwargs['filename'] = os.path.basename(kwargs['filename'])
				sp = Spectrum(**kwargs)
				os.remove(kwargs['filename'])
				return sp
			except urllib2.URLError, ex:
				sys.stderr.write('\nCould not find model file '+kwargs['filename']+' on SPLAT website\n\n')
				os.remove(os.path.basename(kwargs['filename']))
				local = True
	
	# now try local drive
		if (os.path.exists(kwargs['filename']) == False):
			kwargs['filename'] = folder+os.path.basename(kwargs['filename'])
			if (os.path.exists(kwargs['filename']) == False):
				raise NameError('\nCould not find '+kwargs['filename']+' locally\n\n')
			else:
				return Spectrum(**kwargs)
		else:
			return Spectrum(**kwargs)

# or do an interpolated Model
	else:
		return loadInterpolatedModel(**kwargs)
		


def loadModelParameters(*args, **kwargs):
	'''Load up Model Parameters based on help file'''
# keyword parameters
	set = args[0]
	pfile = kwargs.get('parameterFile','parameters.txt')
	
# check if online
	if not checkOnline():
		raise urllib2.URLError('\n\nCannot access online models\n\n')

	url = kwargs.get('folder',SPLAT_URL+'/Models/'+set+'/')
	try:
		open(os.path.basename(pfile), 'wb').write(urllib2.urlopen(url+pfile).read())
		parameters = ascii.read(pfile)
		os.remove(os.path.basename(pfile))
	except urllib2.URLError, ex:
		raise urllib2.URLError('\n\nCannot access models from SPLAT website\n\n')
	
	return parameters



def loadSpectrum(*args, **kwargs):
	'''load up a SpeX spectrum based name, shortname and/or date'''
# keyword parameters
#	name = kwargs.get('name','')
#	shname = kwargs.get('shname','')
#	date = kwargs.get('date','')
	local = kwargs.get('local',False)
	tempfilename = 'temp_model.txt'
	folder = kwargs.get('folder','')
	url = kwargs.get('folder',SPLAT_URL+'/Spectra/')
	kwargs['model'] = False

# check if online
	local = local or (not checkOnline())
	
# CODE NEEDED HERE TO SET UP FILE NAME; FOR NOW JUST ERROR
	
# a filename has been passed
	if (len(args) > 0):
		kwargs['filename'] = args[0]
	else:
		raise NameError('\nNeed to pass in filename for spectral data')

# first try online
	if not local:
		try:
			open(os.path.basename(kwargs['filename']), 'wb').write(urllib2.urlopen(url+kwargs['filename']).read())
			kwargs['filename'] = os.path.basename(kwargs['filename'])
			sp = Spectrum(**kwargs)
			os.remove(os.path.basename(kwargs['filename']))
			return sp
		except urllib2.URLError, ex:
			sys.stderr.write('\nCould not find data file '+kwargs['filename']+' at '+url+'\n\n')
			local = True
	
# now try local drive
	if (os.path.exists(kwargs['filename']) == False):
		kwargs['filename'] = folder+os.path.basename(kwargs['filename'])
		if (os.path.exists(kwargs['filename']) == False):
			raise NameError('\nCould not find '+kwargs['filename']+' locally\n\n')
		else:
			return Spectrum(**kwargs)
	else:
		return Spectrum(**kwargs)


# code to measure a defined index from a spectrum using Monte Carlo noise estimate
# measure method can be mean, median, integrate
# index method can be ratio = 1/2, valley = 1-2/3, OTHERS
# output is index value and uncertainty
def measureIndex(sp,*args,**kwargs):
	'''Measure an index on a spectrum based on defined methodology'''

# keyword parameters
	method = kwargs.get('method','ratio')
	sample = kwargs.get('sample','integrate')
	nsamples = kwargs.get('nsamples',100)
			
# create interpolation functions
	w = numpy.where(sp.flux*sp.noise != numpy.nan)
	f = interp1d(sp.wave[w],sp.flux[w])
	s = interp1d(sp.wave[w],sp.noise[w])
				
# error checking on number of arguments provided
	if (len(args) < 2):
		print 'measureIndex needs at least two samples to function'
		return numpy.nan, numpy.nan
	elif (len(args) < 3 and (method == 'line' or method == 'allers')):
		print method+' requires at least 3 sample regions'
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

# choose function for measuring indices
			if (sample == 'integrate'):
				values[i,j] = trapz(numpy.random.normal(yNum,yNum_e),xNum)
			elif (sample == 'average'):
				values[i,j] = numpy.nanmean(numpy.random.normal(yNum,yNum_e))
			elif (sample == 'median'):
				values[i,j] = numpy.median(numpy.random.normal(yNum,yNum_e))
			elif (sample == 'maximum'):
				values[i,j] = numpy.nanmax(numpy.random.normal(yNum,yNum_e))
			elif (sample == 'minimum'):
				values[i,j] = numpy.nanmin(numpy.random.normal(yNum,yNum_e))
			else:
				values[i,j] = numpy.nanmean(numpy.random.normal(yNum,yNum_e))

# compute index based on defined method
# default is a simple ratio
	if (method == 'ratio'):
		vals = values[0,:]/values[1,:]
	elif (method == 'line'):
		vals = (values[0,:]+values[1,:])/values[2,:]
	elif (method == 'change'):
		vals = 2.*(values[0,:]-values[1,:])/(values[0,:]+values[1,:])
	elif (method == 'allers'):
		vals = (((numpy.mean(args[0])-numpy.mean(args[1]))/(numpy.mean(args[2])-numpy.mean(args[1])))*values[2,:] \
			+ ((numpy.mean(args[2])-numpy.mean(args[0]))/(numpy.mean(args[2])-numpy.mean(args[1])))*values[1,:]) \
			/values[0,:]
	else:
		vals = values[0,:]/values[1,:]
			
# output mean, standard deviation
	return numpy.mean(vals), numpy.std(vals)


# wrapper function for measuring specific sets of indices

def measureIndexSet(sp,**kwargs):

# keyword parameters
	set = kwargs.get('set','burgasser')

# determine combine method
	if (set.lower() == 'burgasser'):
		reference = 'Indices from Burgasser et al. (2006)'
		names = ['H2O-J','CH4-J','H2O-H','CH4-H','H2O-K','CH4-K','K/J']
		inds = numpy.zeros(len(names))
		errs = numpy.zeros(len(names))
		inds[0],errs[0] = measureIndex(sp,[1.14,1.165],[1.26,1.285],method='ratio',sample='integrate',**kwargs)
		inds[1],errs[1] = measureIndex(sp,[1.315,1.335],[1.26,1.285],method='ratio',sample='integrate',**kwargs)
		inds[2],errs[2] = measureIndex(sp,[1.48,1.52],[1.56,1.60],method='ratio',sample='integrate',**kwargs)
		inds[3],errs[3] = measureIndex(sp,[1.635,1.675],[1.56,1.60],method='ratio',sample='integrate',**kwargs)
		inds[4],errs[4] = measureIndex(sp,[1.975,1.995],[2.08,2.12],method='ratio',sample='integrate',**kwargs)
		inds[5],errs[5] = measureIndex(sp,[2.215,2.255],[2.08,2.12],method='ratio',sample='integrate',**kwargs)
		inds[6],errs[6] = measureIndex(sp,[2.06,2.10],[1.25,1.29],method='ratio',sample='integrate',**kwargs)
	elif (set.lower() == 'tokunaga'):
		reference = 'Indices from Tokunaga & Kobayashi (1999)'
		names = ['K1','K2']
		inds = numpy.zeros(len(names))
		errs = numpy.zeros(len(names))
		inds[0],errs[0] = measureIndex(sp,[2.1,2.18],[1.96,2.04],method='change',sample='average',**kwargs)
		inds[1],errs[1] = measureIndex(sp,[2.2,2.28],[2.1,2.18],method='change',sample='average',**kwargs)
	elif (set.lower() == 'reid'):
		reference = 'Indices from Reid et al. (2001)'
		names = ['H2O-A','H2O-B']
		inds = numpy.zeros(len(names))
		errs = numpy.zeros(len(names))
		inds[0],errs[0] = measureIndex(sp,[1.33,1.35],[1.28,1.30],method='ratio',sample='average',**kwargs)
		inds[1],errs[1] = measureIndex(sp,[1.47,1.49],[1.59,1.61],method='ratio',sample='average',**kwargs)
	elif (set.lower() == 'geballe'):
		reference = 'Indices from Geballe et al. (2002)'
		names = ['H2O-1.2','H2O-1.5','CH4-2.2']
		inds = numpy.zeros(len(names))
		errs = numpy.zeros(len(names))
		inds[0],errs[0] = measureIndex(sp,[1.26,1.29],[1.13,1.16],method='ratio',sample='integrate',**kwargs)
		inds[1],errs[1] = measureIndex(sp,[1.57,1.59],[1.46,1.48],method='ratio',sample='integrate',**kwargs)
		inds[2],errs[2] = measureIndex(sp,[2.08,2.12],[2.215,2.255],method='ratio',sample='integrate',**kwargs)
	elif (set.lower() == 'allers'):
		reference = 'Indices from Allers et al. (2007), Allers & Liu (2013)'
		names = ['H2O','FeH-z','VO-z','FeH-J','KI-J','H-cont']
		inds = numpy.zeros(len(names))
		errs = numpy.zeros(len(names))
		inds[0],errs[0] = measureIndex(sp,[1.55,1.56],[1.492,1.502],method='ratio',sample='average',**kwargs)
		inds[1],errs[1] = measureIndex(sp,[0.99135,1.00465],[0.97335,0.98665],[1.01535,1.02865],method='allers',sample='average',**kwargs)
		inds[2],errs[2] = measureIndex(sp,[1.05095,1.06505],[1.02795,1.04205],[1.07995,1.09405],method='allers',sample='average',**kwargs)
		inds[3],errs[3] = measureIndex(sp,[1.19880,1.20120],[1.19080,1.19320],[1.20680,1.20920],method='allers',sample='average',**kwargs)
		inds[4],errs[4] = measureIndex(sp,[1.23570,1.25230],[1.21170,1.22830],[1.26170,1.27830],method='allers',sample='average',**kwargs)
		inds[5],errs[5] = measureIndex(sp,[1.54960,1.57040],[1.45960,1.48040],[1.65960,1.68040],method='allers',sample='average',**kwargs)
	elif (set.lower() == 'slesnick'):
		reference = 'Indices from Slesnick et al. (2004)'
		names = ['H2O-1','H2O-2','FeH']
		inds = numpy.zeros(len(names))
		errs = numpy.zeros(len(names))
		inds[0],errs[0] = measureIndex(sp,[1.335,1.345],[1.295,1.304],method='ratio',sample='average',**kwargs)
		inds[1],errs[1] = measureIndex(sp,[2.035,2.045],[2.145,2.155],method='ratio',sample='average',**kwargs)
		inds[2],errs[2] = measureIndex(sp,[1.1935,1.2065],[1.2235,1.2365],method='ratio',sample='average',**kwargs)
	elif (set.lower() == 'mclean'):
		reference = 'Indices from McLean et al. (2003)'
		names = ['H2OD']
		inds = numpy.zeros(len(names))
		errs = numpy.zeros(len(names))
		inds[0],errs[0] = measureIndex(sp,[1.951,1.977],[2.062,2.088],method='ratio',sample='average',**kwargs)

# output dictionary of indices
	result = {names[i]: (inds[i],errs[i]) for i in numpy.arange(len(names))}
#	result['reference'] = reference
#	return inds,errs,names

	return result



# To do:
# 	masking telluric regions
#	labeling features

def plotSpectrum(*args, **kwargs):

# error check - make sure you're plotting something
	if (len(args) < 1):
		print 'plotSpectrum needs at least on Spectrum object to plot'
		return

# this loop solves issues if a list of spectrum objects is provided accidentally
	sp = []
	for i,a in enumerate(args):
		if (type(a) == list):
			sp.append(a[0])
		else:
			sp.append(a)

# keyword parameters
	title = kwargs.get('title','')
	zeropoint = kwargs.get('zeropoint',[0. for x in range(len(sp))])
	xlabel = kwargs.get('xlabel','{} ({})'.format(sp[0].wlabel,sp[0].wunit))
	ylabel = kwargs.get('ylabel','{} {} ({})'.format(sp[0].fscale,sp[0].flabel,sp[0].funit))
	xrange = kwargs.get('xrange',sp[0].waveRange())
	bound = xrange
	ymax = [s.fluxMax() for s in sp]
	yrange = kwargs.get('yrange',[0,numpy.nanmax(ymax)+numpy.nanmax(zeropoint)])
	bound.extend(yrange)
	grid = kwargs.get('grid',False)
	colors = kwargs.get('colors',['k' for x in range(len(sp))])
	if (len(colors) < len(args)):
		colors.extend(['k' for x in range(len(sp)-len(colors))])
	colorsUnc = kwargs.get('colors',['k' for x in range(len(sp))])
	if (len(colorsUnc) < len(sp)):
		colorsUnc.extend(['k' for x in range(len(sp)-len(colorsUnc))])
	linestyle = kwargs.get('linestyle',['steps' for x in range(len(sp))])
	if (len(linestyle) < len(sp)):
		linestyle.extend(['steps' for x in range(len(sp)-len(linestyle))])
	filename = kwargs.get('filename','')
	format = kwargs.get('format',filename.split('.')[-1])
	showNoise = kwargs.get('showNoise',[False for x in range(len(sp))])
	if not isinstance(showNoise, list):
		showNoise = [showNoise]
	if (len(showNoise) < len(sp)):
		showNoise.extend([True for x in range(len(sp)-len(showNoise))])
	showZero = kwargs.get('showZero',[False for x in numpy.arange(len(sp))])
	if not isinstance(showZero, list):
		showZero = [showZero]
	if (len(showZero) < len(sp)):
		showZero.extend([True for x in range(len(sp)-len(showZero))])
#	mask = kwargs.get('mask',False)				# not yet implemented
#	labels = kwargs.get('labels','')			# not yet implemented
#	features = kwargs.get('features','')		# not yet implemented


#	plt.clf()
# loop through sources
	plt.subplots(1)
	for ii,a in enumerate(sp):
		flx = [i+zeropoint[ii] for i in a.flux]
		plt.plot(a.wave,flx,color=colors[ii],linestyle=linestyle[ii])
# show noise
		if (showNoise[ii]):
			ns = [i+zeropoint[ii] for i in a.noise]
			plt.plot(a.wave,ns,color=colorsUnc[ii],linestyle=linestyle[ii],alpha=0.3)
# zeropoint
		if (showZero[ii]):
			ze = numpy.ones(len(a.flux))*zeropoint[ii]
			plt.plot(a.wave,ze,color=colors[ii],linestyle=':',alpha=0.3)

# grid
	if (grid):
		plt.grid()
# labels
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.axis(bound)
	plt.title(title)
	
# save to filen or display
	if (len(filename) > 0): 
		plt.savefig(filename, format=format)
	
	else:
		plt.show()
		plt.ion()		# make window interactive by default
	
	return


def readSpectrum(**kwargs):

# TO BE DONE:
# FIX IF THERE IS NO NOISE CHANNEL
# PRODUCE AND RETURN HEADER => CHANGE OUTPUT TO DICTIONARY?
	
# keyword parameters
	folder = kwargs.get('folder','./')
	catchSN = kwargs.get('catchSN',True)
	model = kwargs.get('model',False)
	uncertainty = kwargs.get('uncertainty',not model)
	file = kwargs.get('filename','')
	if (os.path.exists(file) == False):
		file = folder+os.path.basename(filename)
	if (os.path.exists(file) == False):
		raise NameError('\nCould not find ' + filename+'\n\n')
	
# determine which type of file
	ftype = file.split('.')[-1]

# fits file	
	if (ftype == 'fit' or ftype == 'fits'):
		data = fits.open(file)
		wave = data[0].data[0,:]
		flux = data[0].data[1,:]
		if (len(data[0].data[:,0]) > 2):
			noise = data[0].data[2,:]
		data.close()

# ascii file	
	else:
		if (uncertainty == True):
			try:
				wave,flux,noise = numpy.genfromtxt(file, comments='#', unpack=True, \
					missing_values = ('NaN','nan'), filling_values = (numpy.nan))
			except ValueError:
				wave,flux,noise = numpy.genfromtxt(file, comments=';', unpack=True, \
	 				missing_values = ('NaN','nan'), filling_values = (numpy.nan))
	 	if (uncertainty == False):
	 		try:
	 			wave,flux = numpy.genfromtxt(file, comments='#', unpack=True, \
	 				missing_values = ('NaN','nan'), filling_values = (numpy.nan))
	 		except ValueError:
	 			wave,flux = numpy.genfromtxt(file, comments=';', unpack=True, \
	 				missing_values = ('NaN','nan'), filling_values = (numpy.nan))

# add in fake uncertainty vector if needed
	if (not uncertainty):
		noise = numpy.zeros(len(flux))
		noise[:] = numpy.nan
  			

# fix places where noise is claimed to be 0
	w = numpy.where(noise == 0.)
	noise[w] = numpy.nan

# fix to catch badly formatted files where noise column is S/N	 			
#	print flux, numpy.median(flux)
	if (catchSN):
  		w = numpy.where(flux > stats.nanmedian(flux))
  		if (stats.nanmedian(flux[w]/noise[w]) < 1.):
  			noise = flux/noise
  			w = numpy.where(numpy.isnan(noise))
  			noise[w] = stats.nanmedian(noise)

	return wave, flux, noise



def searchLibrary(*args, **kwargs):
	'''Search the SpeX database to extract the key reference for that Spectrum
		Note that this is currently only and AND search - need to figure out
		how to a full SQL style search'''

# get database
	data = fetchDatabase(**kwargs)
	sql = Table()
	ref = kwargs.get('output','data_file')

# search parameters
	if kwargs.get('name',False) != False:
		nm = kwargs['name']
		if isinstance(nm,str):
			nm = [nm]
		sql['name'] = nm
	if kwargs.get('designation',False) != False:
		desig = kwargs['designation']
		if isinstance(desig,str):
			desig = [desig]
		sql['designation'] = desig
	if kwargs.get('shortname',False) != False:
		sname = kwargs['shortname']
		if isinstance(sname,str):
			sname = [sname]
		for i,sn in enumerate(sname):
			if sn[0].lower() != 'j':
				sname[i] = 'J'+sname[i]
		sql['shortname'] = sname
	if kwargs.get('date',False) != False:
		sql['observation_date'] = [kwargs['date']]
	if kwargs.get('young',False) != False:
		sql['young'] = [True]
	if kwargs.get('subdwarf',False) != False:
		sql['subdwarf'] = [True]
	if kwargs.get('binary',False) != False:
		sql['binary'] = [True]
	if kwargs.get('spbin',False) != False:
		sql['spbin'] = [True]
	if kwargs.get('red',False) != False:
		sql['red'] = [True]
	if kwargs.get('blue',False) != False:
		sql['blue'] = [True]

# NEED TO ADD IN SEARCH IN AREA, MAGNITUDE RANGE, OBS DATE RANGE, REFERENCES

	if len(sql) > 0:
		result = join(data,sql)
	else:
		result = data
		
	return result[ref]

	


def typeToNum(input, **kwargs):
	'''convert between string and numeric spectral types'''
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
			print 'Spectral type number must be between 0 ({}0) and {} ({}9)'.format(spletter[0],len(spletter)*10.-1.,spletter[-1])
			output = numpy.nan

# spectral type -> number
	else:
		sptype = re.findall('[{}]'.format(spletter),input)
		if (len(sptype) == 1):
			output = spletter.find(sptype[0])*10.
			spind = input.find(sptype[0])+1
			if (input.find('.') < 0):
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
			print 'Only spectral classes {} are handled with this routine'.format(spletter)
			output = numpy.nan
	return output
