# -*- coding: utf-8 -*-
from __future__ import print_function

# this is the test function set for splat evolutionary functions

# imports - internal
import copy
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
import splat.evolve as spel
#import splat as splat

# things to test
# isNumber 


#####################
# TESTING FUNCTIONS #
#####################


def test_loadEvolModel():
# loadEvolModel

# check that models load properly
    for e in EMODELS:
    	m = spel.loadEvolModel(e)
    	for k in list(m.keys()): assert k in EPARAMETERS
    	assert isinstance(m['mass'],list)
    	assert isinstance(m['mass'][0],numpy.ndarray)
    cld = ['hybrid','cloud-free','f2']
    for c in cld:
    	m = spel.loadEvolModel('saumon',z=0.,cloud=c)
    	for k in list(m.keys()): assert k in EPARAMETERS
    	assert isinstance(m['mass'],list)
    	assert isinstance(m['mass'][0],numpy.ndarray)
    metal = [-0.3,0.,0.3]
    for z in metal:
    	m = spel.loadEvolModel('saumon',z=z,cloud='cloud-free')
    	for k in list(m.keys()): assert k in EPARAMETERS
    	assert isinstance(m['mass'],list)
    	assert isinstance(m['mass'][0],numpy.ndarray)

def test_modelParameters():
# _modelParametersSingle (implicit)
# modelParameters

	mkwargs = {'model': 'baraffe'}
	mdl = spel.loadEvolModel(**mkwargs)
	masses = numpy.linspace(numpy.min(mdl['mass'][-1]),numpy.max(mdl['mass'][0]),len(mdl['age']))
	p = spel.modelParameters(mass=masses,age=mdl['age'],**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,numpy.ndarray)
	p = spel.modelParameters(teff=1220,logg=5.1,**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,float)
	l = list(p.keys())
	for k in range(10):
		trl = numpy.random.choice(l,2,replace=False)
		kwargs = {trl[0]: p[trl[0]].value, trl[1]: p[trl[1]].value}
		for mk in list(mkwargs.keys()): kwargs[mk] = mkwargs[mk]
		p2 = spel.modelParameters(**kwargs)
		assert isinstance(p2,dict)
		for e in EPARAMETERS: 
			assert isinstance(p2[e],u.quantity.Quantity)
			assert isinstance(p2[e].value,float)
			assert_allclose(p[e].value,p2[e].value,rtol=1.e-1)

	mkwargs = {'model': 'burrows'}
	mdl = spel.loadEvolModel(**mkwargs)
	masses = numpy.linspace(numpy.min(mdl['mass'][-1]),numpy.max(mdl['mass'][0]),len(mdl['age']))
	p = spel.modelParameters(mass=masses,age=mdl['age'],**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,numpy.ndarray)
	p = spel.modelParameters(teff=1040,logg=5.1,**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,float)
	l = list(p.keys())
	for k in range(10):
		trl = numpy.random.choice(l,2,replace=False)
		kwargs = {trl[0]: p[trl[0]].value, trl[1]: p[trl[1]].value}.update(mkwargs)
		for mk in list(mkwargs.keys()): kwargs[mk] = mkwargs[mk]
		p2 = spel.modelParameters(**kwargs)
		assert isinstance(p2,dict)
		for e in EPARAMETERS: 
			assert isinstance(p2[e],u.quantity.Quantity)
			assert isinstance(p2[e].value,float)
			assert_allclose(p[e].value,p2[e].value,rtol=1.e-1)

	mkwargs = {'model': 'saumon', 'cloud': 'cloud-free'}
	mdl = spel.loadEvolModel(**mkwargs)
	masses = numpy.linspace(numpy.min(mdl['mass'][-1]),numpy.max(mdl['mass'][0]),len(mdl['age']))
	p = spel.modelParameters(mass=masses,age=mdl['age'],**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,numpy.ndarray)
	p = spel.modelParameters(teff=1310,logg=5.3,**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,float)
	l = list(p.keys())
	for k in range(10):
		trl = numpy.random.choice(l,2,replace=False)
		kwargs = {trl[0]: p[trl[0]].value, trl[1]: p[trl[1]].value}.update(mkwargs)
		for mk in list(mkwargs.keys()): kwargs[mk] = mkwargs[mk]
		p2 = spel.modelParameters(**kwargs)
		assert isinstance(p2,dict)
		for e in EPARAMETERS: 
			assert isinstance(p2[e],u.quantity.Quantity)
			assert isinstance(p2[e].value,float)
			assert_allclose(p[e].value,p2[e].value,rtol=1.e-1)

	mkwargs = {'model': 'saumon', 'cloud': 'f2'}
	mdl = spel.loadEvolModel(**mkwargs)
	masses = numpy.linspace(numpy.min(mdl['mass'][-1]),numpy.max(mdl['mass'][0]),len(mdl['age']))
	p = spel.modelParameters(mass=masses,age=mdl['age'],**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,numpy.ndarray)
	p = spel.modelParameters(teff=1310,logg=5.3,**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,float)
	l = list(p.keys())
	for k in range(10):
		trl = numpy.random.choice(l,2,replace=False)
		kwargs = {trl[0]: p[trl[0]].value, trl[1]: p[trl[1]].value}.update(mkwargs)
		for mk in list(mkwargs.keys()): kwargs[mk] = mkwargs[mk]
		p2 = spel.modelParameters(**kwargs)
		assert isinstance(p2,dict)
		for e in EPARAMETERS: 
			assert isinstance(p2[e],u.quantity.Quantity)
			assert isinstance(p2[e].value,float)
			assert_allclose(p[e].value,p2[e].value,rtol=1.e-1)

	mkwargs = {'model': 'saumon', 'cloud': 'hybrid'}
	mdl = spel.loadEvolModel(**mkwargs)
	masses = numpy.linspace(numpy.min(mdl['mass'][-1]),numpy.max(mdl['mass'][0]),len(mdl['age']))
	p = spel.modelParameters(mass=masses,age=mdl['age'],**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,numpy.ndarray)
	p = spel.modelParameters(teff=1310,logg=5.3,**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,float)
	l = list(p.keys())
	for k in range(10):
		trl = numpy.random.choice(l,2,replace=False)
		kwargs = {trl[0]: p[trl[0]].value, trl[1]: p[trl[1]].value}.update(mkwargs)
		for mk in list(mkwargs.keys()): kwargs[mk] = mkwargs[mk]
		p2 = spel.modelParameters(**kwargs)
		assert isinstance(p2,dict)
		for e in EPARAMETERS: 
			assert isinstance(p2[e],u.quantity.Quantity)
			assert isinstance(p2[e].value,float)
			assert_allclose(p[e].value,p2[e].value,rtol=1.e-1)

	mkwargs = {'model': 'saumon', 'z': -0.3}
	mdl = spel.loadEvolModel(**mkwargs)
	masses = numpy.linspace(numpy.min(mdl['mass'][-1]),numpy.max(mdl['mass'][0]),len(mdl['age']))
	p = spel.modelParameters(mass=masses,age=mdl['age'],**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,numpy.ndarray)
	p = spel.modelParameters(teff=1310,logg=5.3,**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,float)
	l = list(p.keys())
	for k in range(10):
		trl = numpy.random.choice(l,2,replace=False)
		kwargs = {trl[0]: p[trl[0]].value, trl[1]: p[trl[1]].value}.update(mkwargs)
		for mk in list(mkwargs.keys()): kwargs[mk] = mkwargs[mk]
		p2 = spel.modelParameters(**kwargs)
		assert isinstance(p2,dict)
		for e in EPARAMETERS: 
			assert isinstance(p2[e],u.quantity.Quantity)
			assert isinstance(p2[e].value,float)
			assert_allclose(p[e].value,p2[e].value,rtol=1.e-1)

	mkwargs = {'model': 'saumon', 'z': 0.3}
	mdl = spel.loadEvolModel(**mkwargs)
	masses = numpy.linspace(numpy.min(mdl['mass'][-1]),numpy.max(mdl['mass'][0]),len(mdl['age']))
	p = spel.modelParameters(mass=masses,age=mdl['age'],**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,numpy.ndarray)
	p = spel.modelParameters(teff=1310,logg=5.3,**mkwargs)
	assert isinstance(p,dict)
	for e in EPARAMETERS: 
		assert isinstance(p[e],u.quantity.Quantity)
		assert isinstance(p[e].value,float)
	l = list(p.keys())
	for k in range(10):
		trl = numpy.random.choice(l,2,replace=False)
		kwargs = {trl[0]: p[trl[0]].value, trl[1]: p[trl[1]].value}.update(mkwargs)
		for mk in list(mkwargs.keys()): kwargs[mk] = mkwargs[mk]
		p2 = spel.modelParameters(**kwargs)
		assert isinstance(p2,dict)
		for e in EPARAMETERS: 
			assert isinstance(p2[e],u.quantity.Quantity)
			assert isinstance(p2[e].value,float)
			assert_allclose(p[e].value,p2[e].value,rtol=1.e-1)



def test_plotModelParameters():
# plotModelParameters
    file = TMPFILENAME+'.png'
    age_samp = 10.**numpy.random.normal(numpy.log10(1.),0.3,50)
    mass_samp = numpy.random.uniform(0.001,0.1,50)
    p = spel.modelParameters(model='baraffe',age=age_samp,mass=mass_samp)
    l = list(p.keys())
    for k in numpy.arange(10):
    	trl = numpy.random.choice(l,2,replace=False)
    	plot = spel.plotModelParameters(p,*trl,model='baraffe',file=file)
    	assert os.access(file,os.R_OK)
    	os.remove(file)


def test_simulate():
# simulateAges
# simulateMasses
# simulateMassRatios
# simulateSpatialDistribution
# simulateKinematics
# simulateBinaryOrbits
# simulateGalacticOrbits
# simulatePhotometry
# simulatePopulation
    num = 10000
    ages = spel.simulateAges(num,distribution='uniform')
    
    masses = spel.simulateMasses(num,distribution='uniform')
    masses = spel.simulateMasses(num,distribution='flat')
    masses = spel.simulateMasses(num,distribution='broken-power-law')
    masses = spel.simulateMasses(num,distribution='chabrier')
    masses = spel.simulateMasses(num,distribution='lognormal')
    for a in numpy.linspace(-2.,2.,10):
    	masses = spel.simulateMasses(num,distribution='power-law',alpha=a)



