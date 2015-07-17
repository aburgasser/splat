# -*- coding: utf-8 -*-
"""
Created on Thu Jul 09 16:53:07 2015

@author: Caleb
"""


import numpy as numpy
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import copy
import sys
from astropy.io import ascii
from matplotlib import cm
from scipy import stats
import os
import pickle
import urllib2
import splat
   
TEN_PARSEC = 443344480  
   
def modelFitMCMC(spec, **kwargs):
    return
   
   
def markovModelDetermination(spec, **kwargs):

    nsample = kwargs.get('nsamples', 1000)
    cutout = kwargs.get('initial_cut', 0.1)  # what fraction of the initial steps are to be discarded
    m_set = kwargs.get('set', 'BTSettl2008')
    step_size = kwargs.get('step_size', [100,0.5])  # the std of the jump size
    plot = kwargs.get('plot', False)
    contour = kwargs.get('contour', False)
    landscape = kwargs.get('landscape', False)
    mask_ranges = kwargs.get('mask_ranges',[])
    mask_telluric = kwargs.get('mask_telluric',False)
    mask_standard = kwargs.get('mask_standard',False)
    xstep = kwargs.get('xstep', 20)
    ystep = kwargs.get('ystep', 20)
    mask = kwargs.get('mask',numpy.zeros(len(spec.wave)))
    calcRadius = kwargs.get('radius', spec.fscale == 'Absolute')
    filename = kwargs.get('filename', spec.filename[:-3] + m_set + '.dat')
    
    if (mask_standard == True):
        mask_telluric == True
   
    if mask_telluric:
        mask_ranges.append([0.,0.65])        # meant to clear out short wavelengths
        mask_ranges.append([1.35,1.42])
        mask_ranges.append([1.8,1.92])
        mask_ranges.append([2.45,99.]) 
        
    if (mask_standard):
        mask_ranges.append([0.,0.8])        # standard short cut
        mask_ranges.append([2.35,99.])      # standard long cut
        
# Mask certain wavelengths and find effective degrees of freedom
    for ranges in mask_ranges:
        mask[numpy.where(((spec.wave.value >= ranges[0]) & (spec.wave.value <= ranges[1])))] = 1
    
    slit_weight = 3
    eff_dof = numpy.round((numpy.nansum(mask) / slit_weight) - 3)
    
    
    # NEED TO MAKE A GUESS
    rang = splat.loadModelParameters(set = m_set) # Range parameters can fall in
    temp_range = rang['teff'][0:2]
    #temp_range[1] = 2200
    grav_range = rang['logg'][0:2]
    z_range = rang['z'][0:2]
    
    tg0 = kwargs.get("initial_guess", [numpy.random.uniform(temp_range[0], \
          temp_range[1]), numpy.random.uniform(grav_range[0], grav_range[1]),\
          numpy.random.uniform(z_range[0], z_range[1])])
    print "initial guess", tg0
    # Checks if initial guess is within range
    if not (temp_range[0] <= tg0[0] <= temp_range[1] and grav_range[0] <= \
            tg0[1] <= grav_range[1] and z_range[0] <= tg0[2] <= z_range[1]):
        sys.stderr.write("Initial guess is out of model range and so it will" + \
                            "default to a random guess.")
        tg0 = [numpy.random.random_integers(temp_range[0], high = temp_range[1]), 
          numpy.random.random_integers(grav_range[0], high = grav_range[1]), 
          numpy.random.random_integers(z_range[0], high = z_range[1])]
 
    model = splat.loadModel(teff = tg0[0], logg = tg0[1], z = tg0[2], 
            mask_telluric = mask_telluric, mask_ranges = mask_ranges, 
            mask = mask, mask_standard = mask_standard, set = m_set)
    
    chisqr0,alpha0 = splat.compareSpectra(spec, model)
    temps = [tg0[0]]
    gravs = [tg0[1]]
    z = [tg0[2]]
    radii = []
    chisqrs = [chisqr0]    
    for i in range(nsample):
        for j in range(len(tg0)):
            # Needed in order to catch some models that do not exist but are 
            # still in the allowed ranges
            try:
            
                tg1 = copy.deepcopy(tg0)
                tg1[j] = numpy.random.normal(tg1[j],step_size[j])
                
                # Makes sure log g and teff stay within model ranges
                while not (temp_range[0] <= tg1[0] <= temp_range[1] and
                grav_range[0] <= tg1[1] <= grav_range[1] and 
                z_range[0] <= tg1[2] <= z_range[1]):
                    tg1[j] = numpy.random.normal(tg1[j],step_size[j])  
            
            
                model = splat.loadModel(teff = tg1[0], logg = tg1[1], 
                        z = tg1[2],mask_telluric = mask_telluric,
                        mask_ranges = mask_ranges, mask = mask, 
                        mask_standard = mask_standard, set = m_set)
                chisqr1,alpha1 = splat.compareSpectra(spec,model, **kwargs)  
                # Probability that it will jump to this new point
                h = 1 - stats.f.cdf(chisqr1/chisqr0, eff_dof, eff_dof)
                # Determines if step will be taken
                if numpy.random.uniform(0,1) < h:
                    tg0[j] = tg1[j]
                    chisqr0 = chisqr1
                    alpha0 = alpha1
            
                # Adds new temp, log g, and chisqr to lists even if they did not change
                temps.append(tg0[0])
                gravs.append(tg0[1])
                z.append(tg0[2])
                if calcRadius:
                    radii.append(TEN_PARSEC*numpy.sqrt(alpha0))
                chisqrs.append(chisqr0)
                
            except:
                continue

    # report results
    cut = int(cutout*len(temps)) # Cuts out intial cutout percent of steps 
    print "Effective Temp", numpy.mean(temps[cut:]),numpy.std(temps[cut:])
    print "Log G", numpy.mean(gravs[cut:]),numpy.std(gravs[cut:])
    print "Metallicity", numpy.mean(z[cut:]),numpy.std(z[cut:])
    if calcRadius:
        print "Radius", numpy.mean(radii[cut:]),numpy.std(radii[cut:])

    if calcRadius:
        data = {'temps': temps[cut:], 'gravs': gravs[cut:], 'zs': z[cut:], 
                'radii': radii[cut:], 'chis': chisqrs[cut:], 
                'temp': numpy.mean(temps[cut:]), 'grav': numpy.mean(gravs[cut:]), 
                'z': numpy.mean(z[cut:]), 'radius': numpy.mean(radii[cut:])}
    else:
        data = {'temps': temps[cut:], 'gravs': gravs[cut:], 'zs': z[cut:], 
                'chis': chisqrs[cut:], 'temp': numpy.mean(temps[cut:]), 
                'grav': numpy.mean(gravs[cut:]), 'z': numpy.mean(z[cut:])}
    
    with open(filename, 'wb') as f:
      pickle.dump(data, f)
      
    if calcRadius:
        return [numpy.mean(temps[cut:]),numpy.std(temps[cut:])], \
               [numpy.mean(gravs[cut:]),numpy.std(gravs[cut:])], \
               [numpy.mean(z[cut:]),numpy.std(z[cut:])], \
               [numpy.mean(radii[cut:]),numpy.std(radii[cut:])]
    else:
        return [numpy.mean(temps[cut:]),numpy.std(temps[cut:])], \
               [numpy.mean(gravs[cut:]),numpy.std(gravs[cut:])], \
               [numpy.mean(z[cut:]),numpy.std(z[cut:])]


def main():

  sp = splat.Spectrum('10548_10753.fits')
  markovModelDetermination(sp, mask_standard=True, initial_guess=[850, 5.3, 1.0], nsamples=10, contour=True)
    
if __name__ == '__main__':
    main() 
