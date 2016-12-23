# -*- coding: utf-8 -*-
from __future__ import print_function

# this is the test function set for splat evolutionary functions

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
from ..initialize import *
from ..utilities import *
from .. import core as splat
#import splat as splat

# things to test
# isNumber 


#####################
# TESTING FUNCTIONS #
#####################


def test_evolModels():
# loadEvolModel
# general model properties
# modelParametersSingle
# modelParameters
    pass

def test_plotModelParameters():
# plotModelParameters
    pass

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
    pass



def OLD_test_readmodel():
    m = loadEvolModel('baraffe')
    print('\nBaraffe')
    for k in list(m.keys()):
        print('\n{}: {}'.format(k,m[k]))
    m = loadEvolModel('burrows')
    print('\nBurrows')
    for k in list(m.keys()):
        print('\n{}: {}'.format(k,m[k]))
    m = loadEvolModel('saumon',z=0.,cloud='hybrid')
    print('\nSaumon z=0 Cloud=hybrid')
    for k in list(m.keys()):
        print('\n{}: {}'.format(k,m[k]))
    m = loadEvolModel('saumon',z=0.,cloud='cloud-free')
    print('\nSaumon z=0 Cloud-free')
    for k in list(m.keys()):
        print('\n{}: {}'.format(k,m[k]))
    m = loadEvolModel('saumon',z=0.,cloud='f2')
    print('\nSaumon z=0 Cloud=f2')
    for k in list(m.keys()):
        print('\n{}: {}'.format(k,m[k]))
    m = loadEvolModel('saumon',z=-0.3,cloud='cloud-free')
    print('\nSaumon z=-0.3 Cloud-free')
    for k in list(m.keys()):
        print('\n{}: {}'.format(k,m[k]))
    m = loadEvolModel('saumon',z=0.3,cloud='cloud-free')
    print('\nSaumon z=+0.3 Cloud-free')
    for k in list(m.keys()):
        print('\n{}: {}'.format(k,m[k]))

def OLD_test_evolve_basic():
    print('\nTesting known Teff and known logg with Baraffe models\n')
    p = modelParameters(temperature=1200, gravity=4.5,model='baraffe')
    for k in p.keys():
        print('{} = {}'.format(k,p[k]))

    print('\nTesting known mass and known logg with Burrows models\n')
    p = modelParameters(mass=0.05, gravity=4.5,model='burrows')
    for k in p.keys():
        print('{} = {}'.format(k,p[k]))

    print('\nTesting known temperature and known age with Saumon solar hybrid models\n')
    p = modelParameters(temperature=1200, age=0.8,model='saumon',metallicity=0.,cloud='hybrid')
    for k in p.keys():
        print('{} = {}'.format(k,p[k]))

    print('\nTesting known temperature and known luminosity with Saumon metal-poor cloud-free models\n')
    p = modelParameters(temperature=1500, lbol=-4.0,model='saumon',z=-0.3,cloud='cloud-free')
    for k in p.keys():
        print('{} = {}'.format(k,p[k]))

    print('\nTesting known temperature and known logg with Saumon metal-rich cloud-free models\n')
    p = modelParameters(temperature=1000, logg=4.4,model='saumon',metallicity=0.3,cloud='cloud-free')
    for k in p.keys():
        print('{} = {}'.format(k,p[k]))

    print('\n')
    return 

def OLD_test_evolve_accuracy(modelname='baraffe', parameter='temperature',metallicity=0.,clouds='nocloud'):
    model = loadEvolModel(modelname,metallicity=metallicity,clouds=clouds)
    masses = model['mass'][-2]
    vals = []
    for i,age in enumerate(model['age']):
        t = []
        for j,m in enumerate(masses):
            if m in model['mass'][i]:
                t.append(numpy.array(model[parameter][i])[numpy.where(model['mass'][i]==m)])
            else:
                t.append(numpy.nan)
        vals.append(t)
    vl = numpy.array(vals).T
    for j,m in enumerate(masses):
        if parameter=='temperature' or parameter=='radius':
            plt.loglog(model['age'],vl[j],color='grey')
        else:
            plt.semilogx(model['age'],vl[j],color='grey')

    age_samp = 10.**(numpy.arange(1,100)/100.*4-3.)
    mass_samp = [0.003,0.005,0.01,0.02,0.05,0.08,0.1,0.2]
    for m in mass_samp:
        p = modelParameters(model,age=age_samp,mass=[m for i in range(len(age_samp))])
        if parameter=='temperature' or parameter=='radius':
            plt.loglog(age_samp,p[parameter],color='r')
        else:
            plt.semilogx(age_samp,p[parameter],color='r')
        plt.ylabel(parameter)
        plt.xlabel('Age')

    plt.show()
    return True


def OLD_test_evolve_accuracy_plotting(xparam,yparam,modelname='baraffe',metallicity=0.,clouds='nocloud',**kwargs):
    model = loadEvolModel(modelname,metallicity=metallicity,clouds=clouds)
#    age_samp = (10.**numpy.random.normal(numpy.log10(1.),0.3,50)).tolist() 
#    age_samp.extend([3]*50)
#    mass_samp = [0.05]*50
#    mass_samp.extend(numpy.random.uniform(0.001,0.1,50).tolist())
    age_samp = 10.**numpy.random.normal(numpy.log10(1.),0.3,50)
    mass_samp = numpy.random.uniform(0.001,0.1,50)
    p = modelParameters(model,age=age_samp,mass=mass_samp)
    plot = plotModelParameters(p,xparam,yparam,model=model,**kwargs)
    return True

def OLD_test_ages(num,**kwargs):
    ages = generateAges(num,**kwargs)
    plt.hist(ages)
    plt.ylabel('Number')
    plt.xlabel('Age')
    plt.show()
    return True

def OLD_test_masses1():
    num = 100000
    mflat1 = generateMasses(num,distribution='uniform')
    mflat2 = generateMasses(num,distribution='flat')
    mflat3 = generateMasses(num,distribution='power-law',alpha=0.)
    plt.hist(mflat1,color='yellow',alpha=0.7)
    plt.hist(mflat2,color='blue',alpha=0.7)
    plt.hist(mflat3,color='red',alpha=0.7)
    plt.show()
    return True

def OLD_test_masses2():
    num = 100000
    mflat1 = generateMasses(num,distribution='power-law',alpha=-2)
    mflat2 = generateMasses(num,distribution='power-law',alpha=-1)
    mflat3 = generateMasses(num,distribution='power-law',alpha=1.)
    mflat4 = generateMasses(num,distribution='power-law',alpha=2.)
    x = numpy.linspace(0.1,1.,1000)
    n,b,p = plt.hist(mflat1,color='yellow',alpha=0.9)
    plt.plot(x,n[-1]*x**2.,color='yellow')
    n,b,p = plt.hist(mflat2,color='blue',alpha=0.9)
    plt.plot(x,n[-1]*x,color='blue')
    n,b,p = plt.hist(mflat3,color='red',alpha=0.9)
    plt.plot(x,n[-1]*x**(-1.),color='red')
    n,b,p = plt.hist(mflat4,color='green',alpha=0.9)
    plt.plot(x,n[-1]*x**(-2.),color='green')
    plt.ylim([0,num/3.])
    plt.show()
    return True

def OLD_test_masses3():
    num = 100000
    mf1 = generateMasses(num,distribution='broken-power-law')
    mf2 = generateMasses(num,distribution='chabrier')
    mf3 = generateMasses(num,distribution='lognormal')
    x = numpy.linspace(0.1,1.,1000)
    n,b,p = plt.hist(mf1,color='green',alpha=0.9)
    n,b,p = plt.hist(mf2,color='blue',alpha=0.9)
    n,b,p = plt.hist(mf3,color='red',alpha=0.9)
    plt.show()
    return True

