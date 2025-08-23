# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
.. note::
         These are the empirical relations functions for SPLAT 
"""

# imports: internal
import copy
import sys

# imports - external
import numpy

from astropy import units as u            # standard units
#from astropy import constants as const        # physical constants in SI units
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# splat functions
from splat.initialize import *
from splat.utilities import *
from splat.photometry import filterMag
from splat.core import classifyByIndex
from splat.citations import shortRef

# Python 2->3 fix for input
try: input=raw_input
except NameError: pass


#Constants.SPLAT_URL = 'http://pono.ucsd.edu/~adam/splat/'
#DATA_FOLDER = '/reference/Spectra/'

#Constants.DB_SOURCES = fetchDatabase(Constants.DB_SOURCES_FILE)
#print(Constants.DB_SOURCES)
#Constants.DB_SPECTRA = fetchDatabase(Constants.DB_SPECTRA_FILE)

# change the command prompt
#sys.ps1 = 'splat empirical> '



#####################################################
##############   DISTANCE ESTIMATION   ##############
#####################################################

def estimateDistance(*args, **kwargs):
    '''
    :Purpose: Takes the apparent magnitude and either takes or determines the absolute
                magnitude, then uses the magnitude/distance relation to estimate the
                distance to the object in parsecs. Returns estimated distance and
                uncertainty in parsecs

    :param sp: Spectrum class object, which should be flux calibrated to its empirical apparent magnitude
    :param mag: apparent magnitude of ``sp``
    :type mag: optional, default = False
    :param mag_unc: uncertainty of the apparent magnitude
    :type mag_unc: optional, default = 0
    :param absmag: absolute magnitude of ``sp``
    :type absmag: optional, default = False
    :param absmag_unc: uncertainty of the absolute magnitude
    :type absmag_unc: optional, default = 0
    :param spt: spectral type of ``sp``
    :type spt: optional, default = False
    :param spt_e: uncertainty of the spectral type
    :type spt_e: optional, default = 0
    :param nsamples: number of samples to use in Monte Carlo error estimation
    :type nsamples: optional, default = 100
    :param filter: Name of filter, must be one of the following:

                    - '2MASS J', '2MASS H', '2MASS Ks'
                    - 'MKO J', 'MKO H', 'MKO K', MKO Kp', 'MKO Ks'
                    - 'NICMOS F090M', 'NICMOS F095N', 'NICMOS F097N', 'NICMOS F108N'
                    - 'NICMOS F110M', 'NICMOS F110W', 'NICMOS F113N', 'NICMOS F140W'
                    - 'NICMOS F145M', 'NICMOS F160W', 'NICMOS F164N', 'NICMOS F165M'
                    - 'NICMOS F166N', 'NICMOS F170M', 'NICMOS F187N', 'NICMOS F190N'
                    - 'NIRC2 J', 'NIRC2 H', 'NIRC2 Kp', 'NIRC2 Ks'
                    - 'WIRC J', 'WIRC H', 'WIRC K', 'WIRC CH4S', 'WIRC CH4L'
                    - 'WIRC CO', 'WIRC PaBeta', 'WIRC BrGamma', 'WIRC Fe2'
                    - 'WISE W1', 'WISE W2'

    :type filter: optional, default = False
    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1555+0954')[0]
    >>> print splat.estimateDistance(sp)
        Please specify the filter used to determine the apparent magnitude
        (nan, nan)
    >>> print splat.estimateDistance(sp, mag = 12.521, mag_unc = 0.022, absmag = 7.24, absmag_unc = 0.50, spt = 'M3')
        (116.36999172188771, 33.124820555524224)
    '''

    mag = kwargs.get('mag', False)
    mag_unc = kwargs.get('mag_unc', 0.)
    mag_unc = kwargs.get('mag_e', mag_unc)
    absmag = kwargs.get('absmag', False)
    absmag_unc = kwargs.get('absmag_unc', 0.)
    absmag_unc = kwargs.get('absmag_e', absmag_unc)
    spt = kwargs.get('spt', False)
    spt_unc = kwargs.get('spt_unc', 0.)
    spt_unc = kwargs.get('spt_e', spt_unc)
    nsamples = kwargs.get('nsamples', 100)
    filt = kwargs.get('filter', False)

# require spectum object if filter, magnitude and spt not all provided
    if mag == False or filt == False or spt == False:
        if len(args) == 0:
            sys.stderr.write('\nYou must include the Spectrum object if you do not specify filter, magnitude and spt\n')
            return numpy.nan, numpy.nan
        else:
            sp = args[0]

# if no apparent magnitude then calculate from spectrum
    if (mag == False):
        if (filt == False):
            sys.stderr.write('\nPlease specify the filter used to determine the apparent magnitude\n')
            return numpy.nan, numpy.nan
        mag, mag_unc = filterMag(sp,filt)

# if no spt then calculate from spectrum
    if spt == False:
        spt, spt_unc = classifyByIndex(sp)


# if no absolute magnitude then estimate from spectral type
    if absmag == False:
        if filt == False:
            sys.stderr.write('\nPlease specify the filter used to determine the absolute magnitude\n')
            return numpy.nan, numpy.nan
        absmag, absmag_unc = typeToMag(spt,filt,unc=spt_unc)
#        print(absmag, absmag_unc)

# create Monte Carlo sets
    if mag_unc > 0.:
        mags = numpy.random.normal(mag, mag_unc, nsamples)
    else:
        mags = nsamples*[mag]

    if absmag_unc > 0.:
        absmags = numpy.random.normal(absmag, absmag_unc, nsamples)
    else:
        absmags = nsamples*[absmag]

# calculate
    distances = 10.**(numpy.subtract(mags,absmags)/5. + 1.)
    d = numpy.mean(distances)
    unc = numpy.std(distances)

    return d, unc



#####################################################
###############   SPT -> PARAMETER   ################
#####################################################



def typeToColor(spt,color,reference='skrzypek2015',uncertainty=0.,nsamples=100,verbose=False,forgiving=True,**kwargs):
    """
    :Purpose: 

        Takes a spectral type and optionally a color (string) and returns the typical color of the source. 

    :Required Inputs:

        :param spt: string, integer, float, list or numpy array specify the spectral type(s) to be evaluated
        :param color: string indicating color; e.g., color='SDSS_I-SDSS_Z'; filter names are checked against splat.FILTERS structure

    :Optional Inputs:

        :param reference: Abs Mag/SpT relation used to compute the absolute magnitude. Options are:

            - *skrzypek* (default): Color trends from `Skryzpek et al. (2015) <http://adsabs.harvard.edu/abs/2015A%26A...574A..78S>`_.
              Spectral type range is M5 to T8
              Colors include SDSS_I-SDSS_Z, SDSS_Z-UKIDSS_Y, UKIDSS_Y-MKO_J, MKO_J-MKO_H, MKO_H-MKO_K, MKO_K-WISE_W1, WISE_W1-WISE_W2, and combinations therein.
            - *leggett*: Color trends from `Leggett et al. (2017) <http://adsabs.harvard.edu/abs/2017ApJ...842..118L>`_.
              Spectral type range is T7.5 to Y2
              Colors include MKO_Y-WFC3_F105W, MKO_J-WFC3_F125W, MKO_J-WFC3_F127M, MKO_H-WFC3_F160W, MKO_H-WIRC_CH4S

        :param uncertainty: uncertainty on input spectral types, should be float, list of floats or numpy array (default = 0)
        :param nsamples: number of Monte Carlo samples for error computation (default = 100)
        :param verbose: Give feedback while in operation (default = False)

    :Output:
        A tuple containing the value and uncertainty (including systematic) of the resulting color.
        If more than one spectral type is given, this tuple contains two numpy arrays

    :Example:
        >>> import splat
        >>> import splat.empirical as spem
        >>> spem.typeToColor('L3', 'MKO_J-MKO_K')
            (1.4266666666666665, 0.07)
        >>> spem.typeToColor(['M5','M6','M7'], 'SDSS_I-SDSS_Z', reference = 'skrzypek', uncertainty=0.5)
            (array([0.91      , 1.4275    , 1.74333333]),
             array([0.5624636 , 0.27662768, 0.13987478]))
        >>> spem.typeToColor('M0', 'SDSS_I-SDSS_Z', ref = 'skrzypek', verbose=True)
            Using the SpT/color trends from skrzypek2015
            Rejected 1 spectral type(s) for being outside relation range
            (nan, nan)
    """

# Keywords alternatives
    for f in ['unc','spt_e','error']:
        if f in list(kwargs.keys()):
            uncertainty = kwargs.get(f,uncertainty)
    for f in ['ref','set','method','model','relation']:
        if f in list(kwargs.keys()):
            reference = kwargs.get(f,reference)
    ref = checkEmpiricalRelation(reference.lower().replace(' ',''),splat.SPT_COLORS_RELATIONS)
    if ref == False:
        print('\nColor set from {} has not be integrated into SPLAT\n\n'.format(reference))
        return numpy.nan, numpy.nan
    if verbose==True: print('\nUsing the SpT/color trends from {}\n'.format(ref))


# Check and convert spectral type variable to an array
#    spt_type = type(spt)
    sptn = copy.deepcopy(spt)
    if isinstance(sptn,str): sptn = [sptn]
    try:
        sptn = list(sptn)
    except:
        sptn = [sptn]
    if isinstance(sptn[0],str):
        sptn = [typeToNum(s) for s in sptn]
    try:
        sptn = numpy.array(sptn)
    except:
        raise ValueError('\nInput spectral type {} must be a string, float, int, list or numpy array'.format(spt))

# Check uncertainties
    uncn = copy.deepcopy(uncertainty)
    if not isinstance(uncn,list) and not isinstance(uncn,numpy.ndarray): uncn = [uncn]
    if len(uncn) == 1:
        uncn = numpy.zeros(len(sptn))+float(uncn[0])

    try:
        uncn = numpy.array(uncn)
    except:
        raise ValueError('\nInput spectral type uncertainty {} must be a float, int, list or numpy array'.format(unc))
    uncn = numpy.abs(uncn)

# Process color requested
    cfilts = (color.upper()).split('-')
    cfilts = [c.strip().replace(' ','_') for c in cfilts]
    cfilts_use = []
    for i,c in enumerate(cfilts):
        filtcheck = checkFilterName(c,verbose=verbose)
        if filtcheck == False: 
            print('\nDid not recognize filter {} in color {} requested'.format(c,color))
            return numpy.nan,numpy.nan
        else: cfilts_use.append(filtcheck)
    col = '{}-{}'.format(cfilts_use[0],cfilts_use[1])

#    cols = .replace('2mass','').replace('sdss','').replace('wise','').replace('denis','')

# NEED TO PUT IN CODE FOR "FORGIVING" FILTER NAMES: MKO_J,2MASS_J -> J ETC.
# AND ALSO INTELLIGENCE TO REVERSE COLOR

# fill in extra colors if request color is not present - a little inefficient right now  
    if col not in list(splat.SPT_COLORS_RELATIONS[ref]['colors'].keys()) and splat.SPT_COLORS_RELATIONS[ref]['method'] == 'interpolate':
# base color 
        c1 = (col.split('-'))[0]
        refcol = ''
        for x in list(splat.SPT_COLORS_RELATIONS[ref]['colors'].keys()):
            if x.split('-')[0] == c1: refcol = x
        if refcol == '': 
            print('\nUnable to constuct color {} for reference set {} which has colors {}\n'.format(color,ref,list(splat.SPT_COLORS_RELATIONS[ref]['colors'].keys())))
            return numpy.nan, numpy.nan
        refspt = numpy.array(splat.SPT_COLORS_RELATIONS[ref]['colors'][refcol]['spt'])
        refcolors = numpy.array(splat.SPT_COLORS_RELATIONS[ref]['colors'][refcol]['values'])
# now run through colors until you create the correct match
        cntr = 0
        maxcntr = len(list(splat.SPT_COLORS_RELATIONS[ref]['colors'].keys()))
        while refcol != col and cntr < maxcntr:
            refadd = ''
            for x in list(splat.SPT_COLORS_RELATIONS[ref]['colors'].keys()):
                if x.split('-')[0] == (refcol.split('-'))[-1]: refadd = x
            if refadd == '': 
                print('\nUnable to constuct color {} for reference set {} which has colors {}\n'.format(color,ref,list(splat.SPT_COLORS_RELATIONS[ref]['colors'].keys())))
                return numpy.nan, numpy.nan
            refcol='{}-{}'.format((refcol.split('-'))[0],(refadd.split('-'))[-1])
# select subset of spts that overlap
            refcolors2 = numpy.array([numpy.nan]*len(refspt))
            for i,s in enumerate(refspt):
                if s in splat.SPT_COLORS_RELATIONS[ref]['colors'][refadd]['spt']:
                    refcolors2[i] = splat.SPT_COLORS_RELATIONS[ref]['colors'][refadd]['values'][splat.SPT_COLORS_RELATIONS[ref]['colors'][refadd]['spt'].index(s)]


            refcolors = refcolors+refcolors2
            splat.SPT_COLORS_RELATIONS[ref]['colors'][refcol] = {'spt': refspt, 'values': refcolors}
            # for k in list(splat.SPT_COLORS_RELATIONS[ref]['colors'][refadd].keys()): splat.SPT_COLORS_RELATIONS[ref]['colors'][refcol][k] = splat.SPT_COLORS_RELATIONS[ref]['colors'][refadd][k]
            # splat.SPT_COLORS_RELATIONS[ref]['colors'][refcol]['values'] = refcolors
            cntr=cntr+1
        if cntr >= maxcntr:
            print('\nUnable to constuct color {} for reference set {} which has colors {}\n'.format(color,ref,list(splat.SPT_COLORS_RELATIONS[ref]['colors'].keys())))
            return numpy.nan, numpy.nan

# still not there? quit
    if col not in list(splat.SPT_COLORS_RELATIONS[ref]['colors'].keys()):
        print('\nUnable to constuct color {} for reference set {} which has colors {}\n'.format(color,ref,list(splat.SPT_COLORS_RELATIONS[ref]['colors'].keys())))
        return numpy.nan, numpy.nan

# conduct calculation for interpolating

    if splat.SPT_COLORS_RELATIONS[ref]['method'] == 'interpolate':
#        f = interp1d(numpy.linspace(splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['range'][0],splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['range'][1]+1,num=len(splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['values'])),splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['values'],bounds_error=False,fill_value=0.)
        f = interp1d(splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['spt'],splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['values'],bounds_error=False,fill_value=numpy.nan)
        output = f(sptn)
# *******************
# THIS NEEDS TO BE UPDATED SO THAT EMPIRICAL UNCERTAINTIES MAPPED TO SPECTRAL TYPES ARE INLCUDED
# *******************
        if 'error' not in list(splat.SPT_COLORS_RELATIONS[ref]['colors'][col].keys()): 
            errors = numpy.zeros(len(output))+splat.SPT_COLORS_RELATIONS[ref]['scatter']
            errors[numpy.isfinite(output)==False] = numpy.nan
        else: 
            g = interp1d(splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['spt'],splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['values'],bounds_error=False,fill_value=numpy.nan)
            errors = g(sptn)

        if uncn[0] > 0.:
            for i,s in enumerate(sptn):
                vals = f(numpy.random.normal(s, uncn[i], nsamples))
                try:
                    std = numpy.nanstd(vals)
                    errors[i] = (errors[i]**2+std**2)**0.5
                except:
                    pass
                    
    elif splat.SPT_COLORS_RELATIONS[ref]['method'] == 'polynomial':
        s = [x-splat.SPT_COLORS_RELATIONS[ref]['sptoffset'] for x in sptn]
        output = numpy.polyval(splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['coeff'],s)
        errors = [splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['fitunc'] for x in sptn]
        if uncn[0] > 0.:
            errors = []
            for i,s in enumerate(sptn):
                vals = numpy.polyval(splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['coeff'],numpy.random.normal(s, uncn[i], nsamples))
                try:
                    std = numpy.nanstd(vals)
                    errors.append((std**2+splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['fitunc']**2)**0.5)
                except:
                    pass
        output[sptn < splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['range'][0]] = numpy.nan
        output[sptn > splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['range'][1]] = numpy.nan
        errors = numpy.array(errors)
        errors[sptn < splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['range'][0]] = numpy.nan
        errors[sptn > splat.SPT_COLORS_RELATIONS[ref]['colors'][col]['range'][1]] = numpy.nan

    else:
        print('\nInterpolation method {} not allowed for typeToColor()\n'.format(color,ref,list(splat.SPT_COLORS_RELATIONS[ref]['method'])))
        return numpy.nan, numpy.nan

# cover out of range values
    if verbose == True:
        rej = [int(numpy.isnan(x)) for x in output]
        if numpy.nansum(rej) > 0:
            print('Rejected {} spectral type(s) for being outside relation range'.format(numpy.nansum(rej)))

    if len(sptn) == 1:
        return output[0],errors[0]
    else:
        return output, errors


def typeToColor_old(spt,color,reference='skrzypek2015',uncertainty=0.,nsamples=100,verbose=False,**kwargs):
    """
    :Purpose: Takes a spectral type and optionally a color (string) and returns the typical color of the source. 
    :param spt: string or integer of the spectral type
    :param color: string indicating color; e.g., color='i-z' (note that case does not matter)
    :type color: optional, default = 'J-K'
    :param ref: Abs Mag/SpT relation used to compute the absolute magnitude. Options are:

        - *skrzypek* (default): Color trends from `Skryzpek et al. (2015) <http://adsabs.harvard.edu/abs/2015A%26A...574A..78S>`_.
          Spectral type range is M5 to T8
          Colors include i-z, z-Y, Y-J, J-H, H-K, K-W1, W1-W2, and combinations therein.


    :type ref: optional, default = 'dupuy'
    :param nsamples: number of Monte Carlo samples for error computation
    :type nsamples: optional, default = 100
    :param unc: uncertainty of ``spt``; if included, returns a tuple with color and uncertainty
    :type unc: optional, default = 0.
    :param verbose: Give feedback while in operation
    :type verbose: optional, default = False
    :Example:
        >>> import splat
        >>> print splat.typeToColor('L3', 'J-K')
            (1.46, nan)
        >>> print splat.typeToColor('M5', 'i-z', ref = 'skrzypek', unc=0.5)
            (0.91, 0.57797809947624645)
        >>> print splat.typeToColor('M0', 'i-z', ref = 'skrzypek')
            Spectral type M0.0 is outside the range for reference set Skrzypek et al. (2015)
            (nan, nan)
    """

# Keywords alternatives
    for f in ['unc','spt_e','error']:
        if f in list(kwargs.keys()):
            uncertainty = kwargs.get(f,uncertainty)
    for f in ['ref','set','method','model','relation']:
        if f in list(kwargs.keys()):
            reference = kwargs.get(f,reference)
    ref = checkEmpiricalRelation(reference.lower().replace(' ',''),splat.SPT_COLORS_RELATIONS)
    if ref == False:
        print('\nColor set from {} has not be integrated into SPLAT\n\n'.format(reference))
        return numpy.nan, numpy.nan
    if verbose==True: print('\nUsing the SpT/color trends from {}\n'.format(ref))

# Convert spectral type string to number
    if isinstance(spt,str): sptn = typeToNum(spt)
    else: sptn = copy.deepcopy(spt)

    col = color.lower().replace(' ','').replace('2mass','').replace('sdss','').replace('wise','').replace('denis','')

# check spt is in range
    if not (splat.SPT_COLORS_RELATIONS[ref]['range'][0] <= sptn <= splat.SPT_COLORS_RELATIONS[ref]['range'][1]):
        print('\n Spectral type {} is outside the range for reference set {}: {} to {}\n\n'.format(spt,ref,splat.SPT_COLORS_RELATIONS[ref]['range'][0],splat.SPT_COLORS_RELATIONS[ref]['range'][1]))
        return numpy.nan, numpy.nan

# fill in extra colors - a little inefficient right now  
    if col not in list(splat.SPT_COLORS_RELATIONS[ref]['values'].keys()):
# base color 
        c1 = (col.split('-'))[0]
        refcol = ''
        for x in list(splat.SPT_COLORS_RELATIONS[ref]['values'].keys()):
            if x.split('-')[0] == c1: refcol = x
        if refcol == '': 
            print('\nUnable to constuct color {} for reference set {} which has colors {}\n'.format(color,ref,list(splat.SPT_COLORS_RELATIONS[ref]['values'].keys())))
            return numpy.nan, numpy.nan
        refcolors = numpy.array(splat.SPT_COLORS_RELATIONS[ref]['values'][refcol])
# now run through colors until you create the correct match
        cntr = 0
        maxcntr = len(list(splat.SPT_COLORS_RELATIONS[ref]['values'].keys()))
        while refcol != col and cntr < maxcntr:
            refadd = ''
            for x in list(splat.SPT_COLORS_RELATIONS[ref]['values'].keys()):
                if x.split('-')[0] == (refcol.split('-'))[-1]: refadd = x
            if refadd == '': 
                print('\nUnable to constuct color {} for reference set {} which has colors {}\n'.format(color,ref,list(splat.SPT_COLORS_RELATIONS[ref]['values'].keys())))
                return numpy.nan, numpy.nan
            refcol='{}-{}'.format((refcol.split('-'))[0],(refadd.split('-'))[-1])
            refcolors = refcolors+numpy.array(splat.SPT_COLORS_RELATIONS[ref]['values'][refadd])
            splat.SPT_COLORS_RELATIONS[ref]['values'][refcol] = refcolors
            cntr=cntr+1
        if cntr >= maxcntr:
            print('\nUnable to constuct color {} for reference set {} which has colors {}\n'.format(color,ref,list(splat.SPT_COLORS_RELATIONS[ref]['values'].keys())))
            return numpy.nan, numpy.nan

    if col in list(splat.SPT_COLORS_RELATIONS[ref]['values'].keys()):
        f = interp1d(numpy.arange(splat.SPT_COLORS_RELATIONS[ref]['range'][0],splat.SPT_COLORS_RELATIONS[ref]['range'][1]+1),splat.SPT_COLORS_RELATIONS[ref]['values'][col],bounds_error=False,fill_value=0.)
        if uncertainty > 0.:
            vals = f(numpy.random.normal(sptn, uncertainty, nsamples))
            return float(f(sptn)), (numpy.nanstd(vals)**2+splat.SPT_COLORS_RELATIONS[ref]['scatter']**2)**0.5
        else:
            return float(f(sptn)), splat.SPT_COLORS_RELATIONS[ref]['scatter']
    else:
        print('\nUnable to constuct color {} for reference set {} which has colors {}\n'.format(color,ref,list(splat.SPT_COLORS_RELATIONS[ref]['values'].keys())))
        return numpy.nan, numpy.nan



def typeToMag_old(spt, filt, **kwargs):
    """
    :Purpose: Takes a spectral type and a filter, and returns absolute magnitude
    :param spt: string or integer of the spectral type
    :param filter: filter of the absolute magnitude. Options are MKO K, MKO H, MKO J, MKO Y, MKO LP, 2MASS J, 2MASS K, or 2MASS H
    :param nsamples: number of Monte Carlo samples for error computation
    :type nsamples: optional, default = 100
    :param unc: uncertainty of ``spt``
    :type unc: optional, default = 0.
    :param ref: Abs Mag/SpT relation used to compute the absolute magnitude. Options are:

        - *burgasser*: Abs Mag/SpT relation from `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_.
          Allowed spectral type range is L0 to T8, and allowed filters are MKO K.
        - *faherty*: Abs Mag/SpT relation from `Faherty et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...752...56F>`_.
          Allowed spectral type range is L0 to T8, and allowed filters are MKO J, MKO H and MKO K.
        - *dupuy*: Abs Mag/SpT relation from `Dupuy & Liu (2012) <http://adsabs.harvard.edu/abs/2012ApJS..201...19D>`_.
          Allowed spectral type range is M6 to T9, and allowed filters are MKO J, MKO Y, MKO H, MKO K, MKO LP, 2MASS J, 2MASS H, and 2MASS K.
        - *filippazzo*: Abs Mag/SpT relation from Filippazzo et al. (2015). Allowed spectral type range is M6 to T9, and allowed filters are 2MASS J and WISE W2.


    :type ref: optional, default = 'dupuy'
    :Example:
        >>> import splat
        >>> print splat.typeToMag('L3', '2MASS J')
            (12.730064813273996, 0.4)
        >>> print splat.typeToMag(21, 'MKO K', ref = 'burgasser')
            (10.705292820099999, 0.26)
        >>> print splat.typeToMag(24, '2MASS J', ref = 'faherty')
            Invalid filter given for Abs Mag/SpT relation from Faherty et al. (2012)
            (nan, nan)
        >>> print splat.typeToMag('M0', '2MASS H', ref = 'dupuy')
            Spectral Type is out of range for Abs Mag/SpT relation from Dupuy & Liu (2012) Abs Mag/SpT relation
            (nan, nan)
    """

#Keywords
    verbose = kwargs.get('verbose',False)
    nsamples = kwargs.get('nsamples', 100)
    ref = kwargs.get('ref', 'dupuy')
    ref = kwargs.get('set', ref)
    unc = kwargs.get('unc', 0.)

#Convert spectral type string to number
    if isinstance(spt,str):
        spt = typeToNum(spt, uncertainty=unc)
    else:
        spt = copy.deepcopy(spt)

#Faherty
    if (ref.lower() == 'faherty'):
        sptoffset = 10.
        reference = 'Abs Mag/SpT relation from Faherty et al. (2012)'
        coeffs = { \
            'MKO J': {'fitunc' : 0.30, 'range' : [20., 38.],  \
                'coeff': [.000203252, -.0129143, .275734, -1.99967, 14.8948]}, \
            'MKO H': {'fitunc' : 0.27, 'range' : [20., 38.], \
                'coeff' : [.000175368, -.0108205, .227363, -1.60036, 13.2372]}, \
            'MKO K': {'fitunc' : 0.28, 'range' : [20., 38.], \
                'coeff' : [.0000816516, -.00469032, .0940816, -.485519, 9.76100]}}

# Burgasser
    elif (ref.lower() == 'burgasser'):
        sptoffset = 20.
        reference = 'Abs Mag/SpT relation from Burgasser (2007)'
        coeffs = { \
            'MKO K': {'fitunc' : 0.26, 'range' : [20., 38.], \
                'coeff': [.0000001051, -.000006985, .0001807, -.002271, .01414, -.04024, .05129, .2322, 10.45]}}

# Dupuy & Liu, default reference
    elif (ref.lower() == 'dupuy'):
        reference = 'Abs Mag/SpT relation from Dupuy & Liu (2012)'
        sptoffset = 10.
        coeffs = { \
            'MKO J': {'fitunc' : 0.39, 'range' : [16., 39.], \
                'coeff' : [-.00000194920, .000227641, -.0103332, .232771, -2.74405, 16.3986, -28.3129]}, \
            'MKO Y': {'fitunc': 0.40, 'range' : [16., 39.], \
                'coeff': [-.00000252638, .000285027, -.0126151, .279438, -3.26895, 19.5444, -35.1560]}, \
            'MKO H': {'fitunc': 0.38, 'range' : [16., 39.], \
                'coeff': [-.00000224083, .000251601, -.0110960, .245209, -2.85705, 16.9138, -29.7306]}, \
            'MKO K': {'fitunc': 0.40, 'range' : [16., 39.], \
                'coeff': [-.00000104935, .000125731, -.00584342, .135177, -1.63930, 10.1248, -15.2200]}, \
            'MKO LP': {'fitunc': 0.28, 'range': [16., 39.], \
                'coeff': [0.00000, 0.00000, .0000546366, -.00293191, .0530581,  -.196584, 8.89928]}, \
            '2MASS J': {'fitunc': 0.40, 'range': [16., 39.], \
                'coeff': [-.000000784614, .000100820, -.00482973, .111715, -1.33053, 8.16362, -9.67994]}, \
            '2MASS H': {'fitunc': 0.40, 'range': [16., 39.], \
                'coeff': [-.00000111499, .000129363, -.00580847, .129202, -1.50370, 9.00279, -11.7526]}, \
            '2MASS KS': {'fitunc': 0.43, 'range':[16., 39.], \
                'coeff': [1.06693e-4, -6.42118e-3, 1.34163e-1, -8.67471e-1, 1.10114e1]}, \
            'WISE W1': {'fitunc': 0.39, 'range':[16., 39.], \
                'coeff': [1.58040e-5, -3.33944e-4, -4.38105e-3, 3.55395e-1, 7.14765]}, \
            'WISE W2': {'fitunc': 0.35, 'range':[16., 39.], \
                'coeff': [1.78555e-5, -8.81973e-4, 1.14325e-2, 1.92354e-1, 7.46564]}}

    elif (ref.lower() == 'filippazzo'):
        reference = 'Abs Mag/SpT relation from Filippazzo et al. (2015)'
        sptoffset = 10.
        coeffs = { \
            '2MASS J': {'fitunc': 0.40, 'range': [16., 39.], \
                'coeff': [3.478e-5, -2.684e-3, 7.771e-2, -1.058e0, 7.157e0, -8.350e0]}, \
            'WISE W2': {'fitunc': 0.40, 'range': [16., 39.], \
                'coeff': [8.190e-6, -6.938e-4, 2.283e-2, -3.655e-1, 3.032e0, -5.043e-1]}}

    else:
        sys.stderr.write('\nInvalid Abs Mag/SpT relation given: %s\n' % ref)
        return numpy.nan, numpy.nan

    if (filt.upper() in coeffs.keys()) == 1:
        for f in coeffs.keys():
            if filt.upper() == f:
                coeff = coeffs[f]['coeff']
                fitunc = coeffs[f]['fitunc']
                rng = coeffs[f]['range']
    else:
        sys.stderr.write('\n Invalid filter {} given for {}\n'.format(filt,reference))
        return numpy.nan, numpy.nan

# compute magnitude if its in the right spectral type range
    if (rng[0] <= spt <= rng[1]):
        abs_mag = numpy.polyval(coeff, spt-sptoffset)
        abs_mag_error = fitunc
        if (unc > 0.):
            vals = numpy.polyval(coeff, numpy.random.normal(spt - sptoffset, unc, nsamples))
#            abs_mag = numpy.nanmean(vals)
            abs_mag_error = (numpy.nanstd(vals)**2+fitunc**2)**0.5
        return abs_mag, abs_mag_error
    else:
        if verbose: sys.stderr.write('\nSpectral Type {} is out of range for {}\n'.format(typeToNum(spt),reference))
        return numpy.nan, numpy.nan


def typeToMag_old2(spt, filt, unc=0.,ref='filippazzo2015',verbose=False,nsamples=100,**kwargs):
    """
    :Purpose: 

    Takes a spectral type and a filter, and returns the expected absolute magnitude based on empirical relations

    :Required Inputs: 

        :param spt: string or integer of the spectral type
        :param filter: filter for which to retrieve absolute magnitude, which must be defined for the given reference set.
        You can check what filters are available by printing splat.SPT_ABSMAG_RELATIONS[reference]['filters'].keys(), where reference is, e.g., 'filippazzo2015'

    :Optional Inputs: 

        :param reference: Abs Mag/SpT relation used to compute the absolute magnitude (also 'ref' and 'set'). These are defined in splat.SPT_ABSMAG_RELATIONS and are currently as follows:

            - *dahn2002*: Abs Mag/SpT relation from `Dahn et al. (2002) <http://adsabs.harvard.edu/abs/2002AJ....124.1170D>`_
              Allowed spectral type range is M7 to L8, and allowed filters are 2MASS J
            - *cruz2003*: Abs Mag/SpT relation from `Cruz et al. (2003) <http://adsabs.harvard.edu/abs/2003AJ....126.2421C>`_
              Allowed spectral type range is M6 to L8, and allowed filters are 2MASS J
            - *tinney2003*: Abs Mag/SpT relation from `Tinney et al. (2003) <http://adsabs.harvard.edu/abs/2003AJ....126..975T>`_.
              Allowed spectral type range is L0 to T7.5, and allowed filters are Cousins I, UKIRT Z, J, K and 2MASS J, Ks
            - *burgasser2007*: Abs Mag/SpT relation from `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_.
              Allowed spectral type range is L0 to T8, and allowed filters are MKO K.
            - *looper2008*: Abs Mag/SpT relation from `Looper et al. (2008) <http://adsabs.harvard.edu/abs/2008ApJ...685.1183L>`_.
              Allowed spectral type range is L0 to T8, and allowed filters are 2MASS J, H, Ks.
            - *dupuy2012*: Abs Mag/SpT relation from `Dupuy & Liu (2012) <http://adsabs.harvard.edu/abs/2012ApJS..201...19D>`_.
              Allowed spectral type range is M6 to T9, and allowed filters are MKO Y, J, H,K, LP, 2MASS J, H, Ks, and WISE W1, W2.
            - *faherty2012*: Abs Mag/SpT relation from `Faherty et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...752...56F>`_.
              Allowed spectral type range is L0 to T8, and allowed filters are MKO J, H, K.
            - *tinney2014*: Abs Mag/SpT relation from `Tinney et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...796...39T>`_.
              Allowed spectral type range is T6.5 to Y2, and allowed filters are MKO J, WISE W2
            - *filippazzo2015* (default): Abs Mag/SpT relation from Filippazzo et al. (2015). 
              Allowed spectral type range is M6 to T9, and allowed filters are 2MASS J and WISE W2.
            - *faherty2016*: Abs Mag/SpT relation for field dwarfs from `Faherty et al. (2016) <http://adsabs.harvard.edu/abs/2016ApJS..225...10F>`_.
              Allowed spectral type range is M6 to T9, and allowed filters are 2MASS J, H, Ks and WISE W1, W2, W3
            - *faherty2016-group*: Abs Mag/SpT relation for "group" dwarfs from `Faherty et al. (2016) <http://adsabs.harvard.edu/abs/2016ApJS..225...10F>`_.
              Allowed spectral type range is M7 to L7, and allowed filters are 2MASS J, H, Ks and WISE W1, W2, W3
            - *faherty2016-young*: Abs Mag/SpT relation for "young" dwarfs from `Faherty et al. (2016) <http://adsabs.harvard.edu/abs/2016ApJS..225...10F>`_.
              Allowed spectral type range is M7 to L7, and allowed filters are 2MASS J, H, Ks and WISE W1, W2, W3

        :param unc: uncertainty of ``spt`` (default = 0)
        :param nsamples: number of Monte Carlo samples for error computation (default = 100)

    :Output: 
    
        2 element tuple providing the absolute magnitude and its uncertainty

    :Example:
        >>> import splat
        >>> print splat.typeToMag('L3', '2MASS J')
            (12.730064813273996, 0.4)
        >>> print splat.typeToMag(21, 'MKO K', ref = 'burgasser')
            (10.705292820099999, 0.26)
        >>> print splat.typeToMag(24, '2MASS J', ref = 'faherty')
            Invalid filter given for Abs Mag/SpT relation from Faherty et al. (2012)
            (nan, nan)
        >>> print splat.typeToMag('M0', '2MASS H', ref = 'dupuy')
            Spectral Type is out of range for Abs Mag/SpT relation from Dupuy & Liu (2012) Abs Mag/SpT relation
            (nan, nan)
    """

#Keywords alternatives
    ref = kwargs.get('reference', ref)
    ref = kwargs.get('set', ref)
    unc = kwargs.get('uncertainty', unc)
    unc = kwargs.get('error', unc)

#Convert spectral type string to number
    if isinstance(spt,str):
        sptn = typeToNum(spt, uncertainty=unc)
    elif isinstance(spt,int) or isinstance(spt,float):
        sptn = copy.deepcopy(spt)
    else:
        raise ValueError('\nInput spectral type {} must be a string, float or int'.format(spt))

# check that you can use the proscribed relation and filter
    filtcheck = checkFilterName(filt,verbose=verbose)
    if filtcheck == False: return numpy.nan,numpy.nan
    else: filt=filtcheck

    refcheck = checkAbsMag(ref,filt=filt,verbose=verbose)
    if refcheck == False: return numpy.nan,numpy.nan
    else: ref=refcheck

    sptoffset = SPT_ABSMAG_RELATIONS[ref]['sptoffset']
    coeff = SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['coeff']
    rng = SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['range']
    fitunc = SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['fitunc']
    refstring = 'Absolute {}/SpT relation from {}'.format(filt,shortRef(SPT_ABSMAG_RELATIONS[ref]['bibcode']))
    if verbose: print('\nUsing {}'.format(refstring))

# compute magnitude if its in the right spectral type range
    if (rng[0] <= sptn <= rng[1]):
        abs_mag = numpy.polyval(coeff, sptn-sptoffset)
        abs_mag_error = fitunc
        if unc > 0.:
            vals = numpy.polyval(coeff, numpy.random.normal(sptn - sptoffset, unc, nsamples))
#            abs_mag = numpy.nanmean(vals)
            abs_mag_error = (numpy.nanstd(vals)**2+fitunc**2)**0.5
        return abs_mag, abs_mag_error
    else:
        if verbose: sys.stderr.write('\nSpectral Type {} is out of range for {}\n'.format(typeToNum(sptn),refstring))
        return numpy.nan, numpy.nan


def typeToMag(spt, filt, uncertainty=0.,reference='filippazzo2015',verbose=False,nsamples=100,mask=True,mask_value=numpy.nan,**kwargs):
    """
    :Purpose: 

    Takes a spectral type and a filter, and returns the expected absolute magnitude based on empirical relations

    :Required Inputs: 

        :param spt: string or integer of the spectral type
        :param filter: filter for which to retrieve absolute magnitude, which must be defined for the given reference set.
        You can check what filters are available by printing splat.SPT_ABSMAG_RELATIONS[reference]['filters'].keys(), where reference is, e.g., 'filippazzo2015'

    :Optional Inputs: 

        :param reference: Abs Mag/SpT relation used to compute the absolute magnitude (also 'ref' and 'set'). These are defined in splat.SPT_ABSMAG_RELATIONS and are currently as follows:

            - *dahn2002*: Abs Mag/SpT relation from `Dahn et al. (2002) <http://adsabs.harvard.edu/abs/2002AJ....124.1170D>`_
              Allowed spectral type range is M7 to L8, and allowed filters are 2MASS J
            - *hawley2002*: Abs Mag/SpT relation from `Hawley et al. (2002) <http://adsabs.harvard.edu/abs/2002AJ....123.3409H>`_
              Allowed spectral type range is M0 to T6, and allowed filters are 2MASS J
            - *cruz2003*: Abs Mag/SpT relation from `Cruz et al. (2003) <http://adsabs.harvard.edu/abs/2003AJ....126.2421C>`_
              Allowed spectral type range is M6 to L8, and allowed filters are 2MASS J
            - *tinney2003*: Abs Mag/SpT relation from `Tinney et al. (2003) <http://adsabs.harvard.edu/abs/2003AJ....126..975T>`_.
              Allowed spectral type range is L0 to T7.5, and allowed filters are Cousins I, UKIRT Z, J, K and 2MASS J, Ks
            - *burgasser2007*: Abs Mag/SpT relation from `Burgasser (2007) <http://adsabs.harvard.edu/abs/2007ApJ...659..655B>`_.
              Allowed spectral type range is L0 to T8, and allowed filters are MKO K.
            - *looper2008*: Abs Mag/SpT relation from `Looper et al. (2008) <http://adsabs.harvard.edu/abs/2008ApJ...685.1183L>`_.
              Allowed spectral type range is L0 to T8, and allowed filters are 2MASS J, H, Ks.
            - *dupuy2012*: Abs Mag/SpT relation from `Dupuy & Liu (2012) <http://adsabs.harvard.edu/abs/2012ApJS..201...19D>`_.
              Allowed spectral type range is M6 to T9, and allowed filters are MKO Y, J, H,K, LP, 2MASS J, H, Ks, and WISE W1, W2.
            - *faherty2012*: Abs Mag/SpT relation from `Faherty et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...752...56F>`_.
              Allowed spectral type range is L0 to T8, and allowed filters are MKO J, H, K.
            - *tinney2014*: Abs Mag/SpT relation from `Tinney et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...796...39T>`_.
              Allowed spectral type range is T6.5 to Y2, and allowed filters are MKO J, WISE W2
            - *filippazzo2015* (default): Abs Mag/SpT relation from Filippazzo et al. (2015). 
              Allowed spectral type range is M6 to T9, and allowed filters are 2MASS J and WISE W2.
            - *faherty2016*: Abs Mag/SpT relation for field dwarfs from `Faherty et al. (2016) <http://adsabs.harvard.edu/abs/2016ApJS..225...10F>`_.
              Allowed spectral type range is M6 to T9, and allowed filters are 2MASS J, H, Ks and WISE W1, W2, W3
            - *faherty2016-group*: Abs Mag/SpT relation for "group" dwarfs from `Faherty et al. (2016) <http://adsabs.harvard.edu/abs/2016ApJS..225...10F>`_.
              Allowed spectral type range is M7 to L7, and allowed filters are 2MASS J, H, Ks and WISE W1, W2, W3
            - *faherty2016-young*: Abs Mag/SpT relation for "young" dwarfs from `Faherty et al. (2016) <http://adsabs.harvard.edu/abs/2016ApJS..225...10F>`_.
              Allowed spectral type range is M7 to L7, and allowed filters are 2MASS J, H, Ks and WISE W1, W2, W3

        :param uncertainty: uncertainty of ``spt`` (default = 0)
        :param nsamples: number of Monte Carlo samples for error computation (default = 100)

    :Output: 
    
        2 element tuple providing the absolute magnitude and its uncertainty

    :Example:
        >>> import splat
        >>> print splat.typeToMag('L3', '2MASS J')
            (12.730064813273996, 0.4)
        >>> print splat.typeToMag(21, 'MKO K', ref = 'burgasser')
            (10.705292820099999, 0.26)
        >>> print splat.typeToMag(24, '2MASS J', ref = 'faherty')
            Invalid filter given for Abs Mag/SpT relation from Faherty et al. (2012)
            (nan, nan)
        >>> print splat.typeToMag('M0', '2MASS H', ref = 'dupuy')
            Spectral Type is out of range for Abs Mag/SpT relation from Dupuy & Liu (2012) Abs Mag/SpT relation
            (nan, nan)
    """

# Keywords alternatives
    for f in ['unc','spt_e','error']:
        if f in list(kwargs.keys()):
            uncertainty = kwargs.get(f,uncertainty)
    unc = copy.deepcopy(uncertainty)
    for f in ['ref','set','method','model','relation']:
        if f in list(kwargs.keys()):
            reference = kwargs.get(f,reference)
    ref = copy.deepcopy(reference)

# Check and convert spectral type variable
#    spt_type = type(spt)
    sptn = copy.deepcopy(spt)
    if isinstance(sptn,str): sptn = [sptn]
    try:
        sptn = list(sptn)
    except:
        sptn = [sptn]
    if isinstance(sptn[0],str):
        sptn = [typeToNum(s) for s in sptn]

    try:
        sptn = numpy.array(sptn)
    except:
        raise ValueError('\nInput spectral type {} must be a string, float, int, list or numpy array'.format(spt))

# Check uncertainties
    uncn = copy.deepcopy(unc)
    if not isinstance(uncn,list) and not isinstance(uncn,numpy.ndarray): uncn = [uncn]
    if len(uncn) == 1:
        uncn = numpy.zeros(len(sptn))+float(uncn[0])

    try:
        uncn = numpy.array(uncn)
    except:
        raise ValueError('\nInput spectral type uncertainty {} must be a float, int, list or numpy array'.format(unc))
    uncn = numpy.abs(uncn)

# check that you can use the proscribed relation and filter
    filtcheck = checkFilterName(filt,verbose=verbose)
    if filtcheck == False: 
        if verbose: print('\nDid not recognize filter {}'.format(filt))
        return numpy.nan,numpy.nan
    else: filt=filtcheck

    refcheck = checkAbsMag(ref,filt=filt,verbose=verbose)
    if refcheck == False: 
        if verbose: print('\nDid not recognize relation {} or filter {} not in this relation'.format(ref,filt))
        return numpy.nan,numpy.nan
    else: ref=refcheck

# read in relevant information
    refstring = 'Absolute {}/SpT relation from {}'.format(filt,SPT_ABSMAG_RELATIONS[ref]['bibcode'])
    if verbose: print('\nUsing {}'.format(refstring))

# polynomial method
    if SPT_ABSMAG_RELATIONS[ref]['method'] == 'polynomial':

# compute values
        abs_mag = numpy.polyval(SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['coeff'], sptn-SPT_ABSMAG_RELATIONS[ref]['sptoffset'])
        abs_mag_error = numpy.zeros(len(sptn))+SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['fitunc']
# mask out absolute magnitudes if they are outside spectral type range
        if mask == True:
            abs_mag[numpy.logical_or(sptn<SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['range'][0],sptn>SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['range'][1])] = mask_value
            abs_mag_error[numpy.logical_or(sptn<SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['range'][0],sptn>SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['range'][1])] = mask_value
            if verbose: print('{} values are outside relation range'.format(len(abs_mag[numpy.logical_or(sptn<SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['range'][0],sptn>SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['range'][1])])))
# perform monte carlo error estimate (slow)
        if numpy.nanmin(uncn) > 0.:
            for i,u in enumerate(uncn):
                if abs_mag[i] != mask_value and abs_mag[i] != numpy.nan:
                    vals = numpy.polyval(SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['coeff'], numpy.random.normal(sptn[i] - SPT_ABSMAG_RELATIONS[ref]['sptoffset'], uncn, nsamples))
                    abs_mag_error[i] = (numpy.nanstd(vals)**2+SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['fitunc']**2)**0.5

# interpolation method
    elif SPT_ABSMAG_RELATIONS[ref]['method'] == 'interpolate':
        rng = [numpy.nanmin(SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['spt']),numpy.nanmax(SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['spt'])]

# compute values
        f = interp1d(SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['spt'],SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['values'],bounds_error=False,fill_value=numpy.nan)
        abs_mag = f(sptn)

# compute uncertainties
        if 'rms' in list(SPT_ABSMAG_RELATIONS[ref]['filters'][filt].keys()):
            fe = interp1d(SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['spt'],SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['rms'],bounds_error=False,fill_value=numpy.nan)
            abs_mag_error = fe(sptn)
        elif 'scatter' in list(SPT_ABSMAG_RELATIONS[ref]['filters'][filt].keys()):
            abs_mag_error = numpy.zeros(len(sptn))+SPT_ABSMAG_RELATIONS[ref]['filters'][filt]['scatter']
        elif 'scatter' in list(SPT_ABSMAG_RELATIONS[ref].keys()):
            abs_mag_error = numpy.zeros(len(sptn))+SPT_ABSMAG_RELATIONS[ref]['scatter']
        else: 
            abs_mag_error = numpy.zeros(len(sptn))
        
# mask out absolute magnitudes if they are outside spectral type range
        if mask == True:
            abs_mag[numpy.logical_or(sptn<rng[0],sptn>rng[1])] = mask_value
            abs_mag_error[numpy.logical_or(sptn<rng[0],sptn>rng[1])] = mask_value
            if verbose: print('{} values are outside relation range'.format(len(abs_mag[numpy.logical_or(sptn<rng[0],sptn>rng[1])])))
# perform monte carlo error estimate (slow)
        if numpy.nanmin(uncn) > 0. and nsamples>1:
            print('doing error calc')
            for i,u in enumerate(uncn):
                if absmag[i] != mask_value and absmag[i] != numpy.nan:
                    vals = f(numpy.random.normal(sptn[i], uncn, nsamples))
                    abs_mag_error[i] = (numpy.nanstd(vals)**2+abs_mag_error[i]**2)**0.5
    else:
        raise ValueError('Unknown method {} for {}'.format(SPT_ABSMAG_RELATIONS[ref]['method'],refstring))


# # compute absolute magnitudes
#     abs_mag = numpy.polyval(coeff, sptn-sptoffset)
#     abs_mag_error = numpy.zeros(len(sptn))+fitunc

# # mask out absolute magnitudes if they are outside spectral type range
#     if mask == True:
#         abs_mag[numpy.logical_or(sptn<rng[0],sptn>rng[1])] = mask_value
#         abs_mag_error[numpy.logical_or(sptn<rng[0],sptn>rng[1])] = mask_value
#         if verbose: print('{} values are outside relation range'.format(len(abs_mag[numpy.logical_or(sptn<rng[0],sptn>rng[1])])))

# # perform monte carlo error estimate (slow)
#     if numpy.nanmin(uncn) > 0.:
#         for i,u in enumerate(uncn):
#             if absmag[i] != mask_value and absmag[i] != numpy.nan:
#                 vals = numpy.polyval(coeff, numpy.random.normal(sptn[i] - sptoffset, uncn, nsamples))
# #            abs_mag = numpy.nanmean(vals)
#                 abs_mag_error[i] = (numpy.nanstd(vals)**2+fitunc**2)**0.5

# return values in same dimension as input
    if len(sptn) == 1: return float(abs_mag[0]),float(abs_mag_error[0])
    else: return abs_mag, abs_mag_error


def typeToTeff(var, uncertainty=[0.], reference='stephens09', nsamples=100, reverse=False, string=False,verbose=False, mask=True, mask_value=numpy.nan, **kwargs):
    '''
    :Purpose: 

    Returns an effective temperature (Teff) and its uncertainty for a given spectral type (or vice versa) based on an empirical relation

    :Required Inputs:

        :param inp: A single or array of either spectral types or effective temperatures (reverse). 
        If spectral types, these can be ints, floats or strings from 0 (K0) and 49.0 (Y9).
        If temperatures, these can be ints or floats and assumed to be in units of Kelvin.
        Note: you must set reverse=True to convert temperature to spectral type.

    :Optional Inputs:

        :param uncertainty: uncertainty of spectral type/temperature (default = 0.001; also 'unc', 'spt_e')
        :param reference: Teff/SpT relation used to compute the effective temperature (also 'set'). Options are:

        - *stephens* (default): Teff/SpT relation from `Stephens et al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...702..154S>`_.
          Allowed spectral type range is M6 to T8 and uses alternate coefficients for L3 to T8.
        - *golimowski*: Teff/SpT relation from `Golimowski et al. (2004) <http://adsabs.harvard.edu/abs/2004AJ....127.3516G>`_.
          Allowed spectral type range is M6 to T8.
        - *looper*: Teff/SpT relation from `Looper et al. (2008) <http://adsabs.harvard.edu/abs/2008ApJ...685.1183L>`_.
          Allowed spectral type range is L0 to T8.
        - *marocco*: Teff/SpT relation from `Marocco et al. (2013) <http://adsabs.harvard.edu/abs/2013AJ....146..161M>`_.
          Allowed spectral type range is M7 to T8.
        - *filippazzo*: Teff/SpT relation from Filippazzo et al. (2015) <http://adsabs.harvard.edu/abs/2015ApJ...810..158F>`_. Allowed spectral type range is M6 to T9.
        - *faherty*: Teff/SpT relation from Faherty et al. (2016) <http://adsabs.harvard.edu/abs/2016ApJS..225...10F>`_. 
        This relation is defined for normal dwarfs (M7 < SpT < T8), young dwarfs (``-young`` and ``-young2``, M7 < SpT < L7),
        and young dwarfs in groups (``-group``, M7 < SpT < L7)
        - *dupuy*: Teff/SpT relation from Dupuy et al. (2017).  
        This relation is defined for Saumon & Marley (2008) models (``-saumon``, L1.5 < SpT < T5) and Lyon models (``-lyon``, M7 < SpT < T5)

        :param reverse: set to True to convert effective temperature to spectral type (default=False)
        :param nsamples: number of samples to use in Monte Carlo error estimation (default=100)

    :Output:

        A 2-element tuple containing the Teff/SpT and its uncertainty

    :Example:
        >>> import splat
        >>> print splat.typeToTeff(20)
            (2233.4796740905499, 100.00007874571999)
        >>> print splat.typeToTeff(20, unc = 0.3, ref = 'golimowski')
            (2305.7500497902788, 127.62548366132124)
    '''
# Keywords alternatives
    for f in ['unc','spt_e','error']:
        if f in list(kwargs.keys()):
            uncertainty = kwargs.get(f,uncertainty)
    unc = copy.deepcopy(uncertainty)
    for f in ['ref','set','method','model','relation']:
        if f in list(kwargs.keys()):
            reference = kwargs.get(f,reference)
    ref = checkEmpiricalRelation(reference.lower().replace(' ',''),splat.SPT_TEFF_RELATIONS)
    if ref == False:
        print('\nSpT/Teff relation from {} has not be integrated into SPLAT\n\n'.format(reference))
        return numpy.nan, numpy.nan
    if verbose==True: print('\nUsing the SpT/Teff relation from {}\n'.format(ref))

# Check and convert input variable
    inp = copy.deepcopy(var)
    if isUnit(inp): inp = inp.value
    if type(inp) in [int,float,str,numpy.float64]: inp = [inp]
    try:
        inp = list(inp)
    except:
        raise ValueError('\nInput variable {} must be a string, float, int, list, or numpy array'.format(inp))

# Convert spectral type string to number
    if isinstance(inp[0],str):
        inp = [typeToNum(i) for i in inp]

# Check and convert uncertainty variable
    if type(unc) in [int,float]: unc = [unc]
    try:
        unc = list(unc)
    except:
        if verbose==True: print('\nInput uncertainty {} must be a string, float, int, list, or numpy array; ignoring'.format(unc))
        unc = list(numpy.zeros(len(inp)))
    while len(unc) < len(inp): unc.append(unc[-1])

# some special relations
    if 'stephens' in ref.lower() and 'alt' in ref.lower(): ref = 'stephens-alt'
    if 'faherty' in ref.lower() and 'young2' in ref.lower(): ref = 'faherty-young2'
    elif 'faherty' in ref.lower() and 'young' in ref.lower(): ref = 'faherty-young'
    elif 'faherty' in ref.lower() and 'group' in ref.lower(): ref = 'faherty-group'
    else: pass
    if 'dupuy' in ref.lower() and 'saumon' in ref.lower(): ref = 'dupuy-saumon'
    elif 'dupuy' in ref.lower() and 'lyon' in ref.lower(): ref = 'dupuy-lyon'
    else: pass

# check that you can use the proscribed relation and filter
    refcheck = checkDict(ref,SPT_TEFF_RELATIONS,replace=[[' ','_'],['-','_']],verbose=verbose)
    if refcheck == False: 
        raise ValueError(print('\nDid not recognize relation reference {}; try {}'.format(reference,list(SPT_TEFF_RELATIONS.keys()))))
    else: ref=refcheck

# read in relevant information
#    refstr = shortRef(SPT_TEFF_RELATIONS[ref]['bibcode'])
    refstr = SPT_TEFF_RELATIONS[ref]['reference']
    refstring = 'Teff/SpT relation from {}'.format(refstr)
    if verbose: print('\nUsing {}'.format(refstring))

# convert teff into spt
    if reverse == True:
# polynomial method      
        if SPT_TEFF_RELATIONS[ref]['method'] == 'polynomial':
            sptoffset = SPT_TEFF_RELATIONS[ref]['sptoffset']
            fitunc = SPT_TEFF_RELATIONS[ref]['fitunc']
            rng = numpy.array(SPT_TEFF_RELATIONS[ref]['range'])
            coeff = numpy.array(SPT_TEFF_RELATIONS[ref]['coeff'])

            x = numpy.linspace(rng[0]-sptoffset,rng[1]-sptoffset,nsamples)
            rng = numpy.array([numpy.nanmin(numpy.polyval(coeff,x)),numpy.nanmax(numpy.polyval(coeff,x))])
            f = interp1d(numpy.polyval(coeff,x),x,bounds_error=False)
            spto = f(inp)+sptoffset
            spto_e = 0.5*numpy.absolute(f(numpy.array(inp)+fitunc)-f(numpy.array(inp)-fitunc))

# interpolation method
        elif SPT_TEFF_RELATIONS[ref]['method'] == 'interpolate':
            rng = [numpy.nanmin(SPT_TEFF_RELATIONS[ref]['values']),numpy.nanmax(SPT_TEFF_RELATIONS[ref]['values'])]
            xspt = SPT_TEFF_RELATIONS[ref]['spt']
            xval = SPT_TEFF_RELATIONS[ref]['values']
            xrms = SPT_TEFF_RELATIONS[ref]['rms']
            sptoffset = 0.

            f = interp1d(xval,xspt,bounds_error=False)
            n = interp1d(xval,xrms,bounds_error=False)
            spto = f(inp)
            fitunc = numpy.nanmedian(n(inp))
# estimate a relation uncertainty from average scatter for all measures
            spto_e = 0.5*numpy.absolute(f(numpy.array(inp)+n(inp))-f(numpy.array(inp)-n(inp)))

# FAIL: estimate a single relation uncertainty from first measure (for single temperature)
        if numpy.isfinite(spto_e).all() == False:
            spto_e = numpy.zeros(len(spto))+0.5*numpy.absolute(f(inp[0]+fitunc)-f(inp[0]-fitunc))
            if verbose==True: print('using sysunc from first measure')
# FAIL: estimate a single relation uncertainty from middle of relation (for many temperatues)
        if numpy.isfinite(spto_e).all() == False:
            spto_e = numpy.zeros(len(spto))+0.5*numpy.absolute(f(numpy.nanmedian(inp)+fitunc)-f(numpy.nanmedian(inp)-fitunc))      
            if verbose==True: print('using sysunc from middle')
# FAIL: just use a 0.5 subtype error
        if numpy.isfinite(spto_e).all() == False: 
            spto_e = numpy.zeros(len(spto))+0.5
            if verbose==True: print('using fixed sysunc')
            
# fold in uncertainty of measurement along with relation error
        if unc[0] > 0.:
            for i,t in enumerate(inp):
                x = numpy.random.normal(t,unc[i],nsamples)
                vals = f(x)+sptoffset
                spto_e[i] = (spto_e[i]**2+numpy.nanstd(vals)**2)**0.5
# mask out of range
        if mask == True:
            spto[numpy.logical_or(inp<rng[0],inp>rng[1])] = mask_value
            spto_e[numpy.logical_or(inp<rng[0],inp>rng[1])] = mask_value
            if verbose: print('{} values are outside relation range'.format(len(spto[numpy.logical_or(inp<rng[0],inp>rng[1])])))

# convert to strings if desired
        if string == True: spto = [typeToNum(s) for s in spto]

# return values
        if len(inp) == 1: return spto[0],spto_e[0]
        else: return list(spto),list(spto_e)

        
#         vals = f(x)+sptoffset
#             if 'stephens' in ref.lower():
#                 if numpy.min(numpy.polyval(coeff_alt,[r-sptoffset for r in range_alt])) <= teff <= numpy.max(numpy.polyval(coeff_alt,[r-sptoffset for r in range_alt])):
#                     x = numpy.linspace(range_alt[1]-sptoffset,range_alt[0]-sptoffset,nsamples)
#                     f = interp1d(numpy.polyval(coeff_alt,x),x,bounds_error=False)
#                     spto = float(f(teff))+sptoffset
#                     if kwargs.get('string',False) == True: spto = splat.typeToNum(spto)
#                     x = numpy.random.normal(teff,teff_e,nsamples)+sptoffset
#                     vals = f(x)
# # assuming an at least 0.5 spectral type uncertainty
#             spto_e = (numpy.nanstd(vals)**2+0.5**2)**0.5
#             return spto, spto_e
#         else:
#             if verbose: sys.stderr.write('\nTeff is out of range for {:s} Teff/SpT relation\n'.format(reference))
#             return numpy.nan, numpy.nan

# convert spt into teff
    else:
# polynomial method
        if SPT_TEFF_RELATIONS[ref]['method'] == 'polynomial':
            sptoffset = SPT_TEFF_RELATIONS[ref]['sptoffset']
            rng = numpy.array(SPT_TEFF_RELATIONS[ref]['range'])
            fitunc = SPT_TEFF_RELATIONS[ref]['fitunc']
            coeff = numpy.array(SPT_TEFF_RELATIONS[ref]['coeff'])

# compute values
            teff = numpy.polyval(coeff, numpy.array(inp)-sptoffset)
            teff_e = numpy.zeros(len(inp))+fitunc
# mask out temperatures if they are outside spectral type range
            if mask == True:
                teff[numpy.logical_or(inp<rng[0],inp>rng[1])] = mask_value
                teff_e[numpy.logical_or(inp<rng[0],inp>rng[1])] = mask_value
                if verbose == True: print('{} values are outside relation range'.format(len(teff[numpy.logical_or(inp<rng[0],inp>rng[1])])))
# perform monte carlo error estimate (slow)
            if numpy.nanmin(unc) > 0.:
                for i,err in enumerate(unc):
                    if teff[i] != mask_value and teff[i] != numpy.nan:
                        vals = numpy.polyval(coeff, numpy.random.normal(inp[i]-sptoffset, err, nsamples))
                        teff_e[i] = (numpy.nanstd(vals)**2+fitunc**2)**0.5

# interpolation method
        elif SPT_TEFF_RELATIONS[ref]['method'] == 'interpolate':
            rng = [numpy.nanmin(SPT_TEFF_RELATIONS[ref]['spt']),numpy.nanmax(SPT_TEFF_RELATIONS[ref]['spt'])]
            xspt = SPT_TEFF_RELATIONS[ref]['spt']
            xval = SPT_TEFF_RELATIONS[ref]['values']
            xrms = SPT_TEFF_RELATIONS[ref]['rms']
            sptoffset = 0.

# compute values
            f = interp1d(xspt,xval,bounds_error=False,fill_value=numpy.nan)
            fe = interp1d(xspt,xrms,bounds_error=False,fill_value=numpy.nan)
            teff = f(inp)
            teff_e = fe(inp)
# mask out absolute magnitudes if they are outside spectral type range
            if mask == True:
                teff[numpy.logical_or(inp<rng[0],inp>rng[1])] = mask_value
                teff_e[numpy.logical_or(inp<rng[0],inp>rng[1])] = mask_value
                if verbose: print('{} values are outside relation range'.format(len(teff[numpy.logical_or(inp<rng[0],inp>rng[1])])))
# perform monte carlo error estimate (slow)
            if numpy.nanmin(unc) > 0.:
                for i,err in enumerate(unc):
                    if teff[i] != mask_value and teff[i] != numpy.nan:
                        vals = f(numpy.random.normal(inp[i], err, nsamples))
                        teff_e[i] = (numpy.nanstd(vals)**2+teff_e[i]**2)**0.5
        else:
            raise ValueError('Unknown method {} for {}'.format(SPT_TEFF_RELATIONS[ref]['method'],refstring))


# return values
        if len(inp) == 1: return teff[0]*u.K,teff_e[0]*u.K
        else: return teff*u.K,teff_e*u.K


def redden(sp, **kwargs):
    '''
    Description:
      Redden a spectrum based on an either Mie theory or a standard interstellar profile
      using Cardelli, Clayton, and Mathis (1989 ApJ. 345, 245)

    **Usage**

       >>> import splat
       >>> sp = splat.Spectrum(10001)                   # read in a source
       >>> spr = splat.redden(sp,av=5.,rv=3.2)          # redden to equivalent of AV=5

    **Note**
      This routine is still in beta form; only the CCM89 currently works

    '''
    w = sp.wave.value                           # assuming in microns!
    av = kwargs.get('av',0.0)


    if kwargs.get('mie',False):                 # NOT CURRENTLY FUNCTIONING
        a = kwargs.get('a',10.)                 # grain size
        n = kwargs.get('n',1.33)                # complex index of refraction
        x = 2*numpy.pi*a/w
        x0 = 2.*numpy.pi*a/0.55                 # for V-band
        qabs = -4.*x*((n**2-1)/(n**2+2)).imag
        qsca = (8./3.)*(x**4)*(((n**2-1)/(n**2+2))**2).real
#        tau = numpy.pi*(a**2)*(qabs+qsca)
        tau = 1.5*(qabs+qsca)/a    # for constant mass
        qabs0 = -4.*x0*((n**2-1)/(n**2+2)).imag
        qsca0 = (8./3.)*(x0**4)*(((n**2-1)/(n**2+2))**2).real
#        tau0 = numpy.pi*(a**2)*(qabs0+qsca0)
        tau0 = 1.5*(qabs0+qsca0)/a    # for constant mass
        scale = (10.**(-0.4*av))
        absfrac = scale*numpy.exp(numpy.max(tau)-tau)
    else:
        x = 1./w
        a = 0.574*(x**1.61)
        b = -0.527*(x**1.61)
        rv = kwargs.get('rv',3.1)
        absfrac = 10.**(-0.4*av*(a+b/rv))

    if kwargs.get('normalize',False):
        absfrac = absfrac/numpy.median(absfrac)

#    print(tau0, min(tau), max(tau), max(absfrac), min(absfrac))
    spabs = splat.Spectrum(wave=w,flux=absfrac)
    return sp*spabs


def typeToLbol(*args,**kwargs):
    return typeToLuminosity(*args,**kwargs)

def typeToLuminosity(spt, uncertainty=0.,reference='filippazzo2015',verbose=False,nsamples=100,reverse=False,**kwargs):
    """
    :Purpose: 

    Takes a spectral type and returns the expected scaled log luminosity (log Lbol/Lsun) based on empirical relations

    :Required Inputs: 

        :param spt: string or integer of the spectral type

    :Optional Inputs: 

        :param reference: log Lbol/SpT relation reference (also 'ref' and 'set'). These are defined in splat.SPT_LBOL_RELATIONS and are currently as follows:

            - *filippazzo2015* (default): Lbol/SpT relation from `Filippazzo et al. (2015) <http://adsabs.harvard.edu/abs/2013Sci...341.1492D>`_
              Allowed spectral type range is M6 to T9

        :param uncertainty: uncertainty of ``spt`` (default = 0)
        :param reverse: apply reverse approach: given BC, infer spectral type
        :param nsamples: number of Monte Carlo samples for error computation (default = 100)

    :Output: 
    
        2 element tuple providing the absolute magnitude and its uncertainty

    :Example:
        >>> import splat
        >>> print splat.typeToLuminosity('L3')
    """

# Keywords alternatives
    for f in ['unc','spt_e','error']:
        if f in list(kwargs.keys()):
            uncertainty = kwargs.get(f,uncertainty)
    unc = copy.deepcopy(uncertainty)
    for f in ['ref','set','method','model','relation']:
        if f in list(kwargs.keys()):
            reference = kwargs.get(f,reference)
    ref = checkEmpiricalRelation(reference.lower().replace(' ',''),splat.SPT_LBOL_RELATIONS)
    if ref == False:
        print('\nSpT/Lbol relation from {} has not be integrated into SPLAT\n\n'.format(reference))
        return numpy.nan, numpy.nan
    refstring = 'Luminosity/SpT relation for from {}'.format(shortRef(SPT_LBOL_RELATIONS[ref]['bibcode']))
    if verbose: print('\nUsing {}'.format(refstring))

# normal approach: SpT -> Lbol
    if reverse == False:

#Convert spectral type string to number
        if isinstance(spt,str):
            sptn = typeToNum(spt, uncertainty=unc)
        elif isinstance(spt,int) or isinstance(spt,float):
            sptn = copy.deepcopy(spt)
        else:
            raise ValueError('\nInput spectral type {} must be a string, float or int'.format(spt))

# polynomial method
        if SPT_LBOL_RELATIONS[ref]['method'] == 'polynomial':
            rng = SPT_LBOL_RELATIONS[ref]['range']
            if (rng[0] <= sptn <= rng[1]):
                lbol = numpy.polyval(SPT_LBOL_RELATIONS[ref]['coeff'], sptn-SPT_LBOL_RELATIONS[ref]['sptoffset'])
                lbol_error = SPT_LBOL_RELATIONS[ref]['fitunc']
                if unc > 0.:
                    vals = numpy.polyval(SPT_LBOL_RELATIONS[ref]['coeff'], numpy.random.normal(sptn - SPT_LBOL_RELATIONS[ref]['sptoffset'], unc, nsamples))
                    lbol_error = (numpy.nanstd(vals)**2+lbol_error**2)**0.5
                return lbol, lbol_error
            else:
                if verbose: sys.stderr.write('\nSpectral type {} is out of range for {}'.format(typeToNum(sptn),refstring))
                return numpy.nan, numpy.nan

# interpolation method
        elif SPT_LBOL_RELATIONS[ref]['method'] == 'interpolate':
            rng = [numpy.nanmin(SPT_LBOL_RELATIONS[ref]['spt']),numpy.nanmax(SPT_LBOL_RELATIONS[ref]['spt'])]
            if (rng[0] <= sptn <= rng[1]):
                f = interp1d(SPT_LBOL_RELATIONS[ref]['spt'],SPT_LBOL_RELATIONS[ref]['bc'])
                fe = interp1d(SPT_LBOL_RELATIONS[ref]['spt'],SPT_LBOL_RELATIONS[ref]['rms'])
                lbol = float(f(sptn))
                lbol_error = float(fe(sptn))
                if unc > 0.:
                    vals = f(numpy.random.normal(sptn, unc, nsamples))
                    lbol_error = (numpy.nanstd(vals)**2+bc_error**2)**0.5
                return lbol, lbol_error
            else:
                if verbose: sys.stderr.write('\nSpectral type {} is out of range for {}'.format(typeToNum(sptn),refstring))
                return numpy.nan, numpy.nan
        else:
            raise ValueError('Unknown method {} for {}'.format(SPT_LBOL_RELATIONS[ref]['method'],refstring))

# reverse approach: Lbol -> SpT
    else:
        if not isinstance(spt,float): raise ValueError('Running this in reverse you need to provide a log luminosity value instead of {}'.format(spt))
        if SPT_LBOL_RELATIONS[ref]['method'] == 'polynomial':
            rng = SPT_LBOL_RELATIONS[ref]['range']
            x = numpy.linspace(rng[0],rng[1],nsamples)
            y = numpy.polyval(SPT_LBOL_RELATIONS[ref]['coeff'], x-SPT_LBOL_RELATIONS[ref]['sptoffset'])
            f = interp1d(y,x)
            try:
                lbol = float(f(spt))
            except:
                if verbose: print('\nlog luminosity value {} is outside the range expected for {}'.format(spt,refstring))
                return numpy.nan, numpy.nan
            vals = []
            for i in range(nsamples):
                ye = y+numpy.random.normal(0,SPT_LBOL_RELATIONS[ref]['fitunc'])
                f = interp1d(ye,x)
                try:
                    vals.append(f(numpy.random.normal(spt,unc)))
                except:
                    pass
            lbol_error = numpy.nanstd(vals)
            return lbol, lbol_error
        elif SPT_LBOL_RELATIONS[ref]['method'] == 'interpolate':
            f = interp1d(SPT_LBOL_RELATIONS[ref]['bc'],SPT_LBOL_RELATIONS[ref]['spt'])
            try:
                lbol = f(spt)
            except:
                if verbose: print('\nlog luminosity value {} is outside the range expected for {}'.format(spt,refstring))
                return numpy.nan, numpy.nan
            vals = []
            for i in range(nsamples):
                y = numpy.random.normal(SPT_LBOL_RELATIONS[ref]['bc'],SPT_LBOL_RELATIONS[ref]['rms'])
                f = interp1d(y,SPT_LBOL_RELATIONS[ref]['spt'])
                try:
                    vals.append(f(numpy.random.normal(spt,unc)))
                except:
                    pass
            lbol_error = numpy.nanstd(vals)
            return lbol, lbol_error
        else:
            raise ValueError('Unknown method {} for {}'.format(SPT_LBOL_RELATIONS[ref]['method'],refstring))



def typeToBC(spt, filt, uncertainty=0.,reference='filippazzo2015',verbose=False,nsamples=100,reverse=False,**kwargs):
    """
    :Purpose: 

    Takes a spectral type and a filter, and returns the expected bolometric correction BC = M_bol - M_filter

    :Required Inputs: 

        :param spt: string or integer of the spectral type
        :param filter: filter for which to retrieve absolute magnitude, which must be defined for the given reference set.
        You can check what filters are available by printing splat.SPT_BC_RELATIONS[reference]['filters'].keys(), where reference is, e.g., 'filippazzo2015'

    :Optional Inputs: 

        :param reference: Abs Mag/SpT relation used to compute the absolute magnitude (also 'ref' and 'set'). These are defined in splat.SPT_BC_RELATIONS and are currently as follows:

            - *liu2010*: BC/SpT relation from `Liu et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJ...722..311L>`_
              Allowed spectral type range is M6 to T8.5, and allowed filters are MKO J, H, K
            - *dupuy2013*: BC/SpT relation from `Dupuy & Kraus (2013) <http://adsabs.harvard.edu/abs/2013Sci...341.1492D>`_
              Allowed spectral type range is T8 to Y0.5, and allowed filters are MKO Y, J, H
            - *filippazzo2015* (default): BC/SpT relation from `Filippazzo et al. (2015) <http://adsabs.harvard.edu/abs/2013Sci...341.1492D>`_
              Allowed spectral type range is M6 to T8/9, and allowed filters are 2MASS J, Ks
            - *filippazzo2015-young*: BC/SpT relation for young sources from `Filippazzo et al. (2015) <http://adsabs.harvard.edu/abs/2013Sci...341.1492D>`_
              Allowed spectral type range is M7 to T8, and allowed filters are 2MASS J, Ks

        :param uncertainty: uncertainty of ``spt`` (default = 0)
        :param reverse: apply reverse approach: given BC, infer spectral type
        :param nsamples: number of Monte Carlo samples for error computation (default = 100)

    :Output: 
    
        2 element tuple providing the absolute magnitude and its uncertainty

    :Example:
        >>> import splat
        >>> print splat.typeToBC('L3', '2MASS J')
    """

# Keywords alternatives
    for f in ['unc','spt_e','error']:
        if f in list(kwargs.keys()):
            uncertainty = kwargs.get(f,uncertainty)
    unc = copy.deepcopy(uncertainty)
    for f in ['ref','set','method','model','relation']:
        if f in list(kwargs.keys()):
            reference = kwargs.get(f,reference)
    ref = copy.deepcopy(reference)


# check that you can use the proscribed relation and filter
    filtcheck = checkFilterName(filt,verbose=verbose)
    if filtcheck == False: return numpy.nan,numpy.nan
    else: filt=filtcheck

    refcheck = checkBC(ref,filt=filt,verbose=verbose)
    if refcheck == False: return numpy.nan,numpy.nan
    else: ref=refcheck

    if verbose: 
        refstring = 'BC/SpT relation for filter {} from {}'.format(filt,shortRef(SPT_BC_RELATIONS[ref]['bibcode']))
        print('\nUsing {}'.format(refstring))

# normal approach: SpT -> BC
    if reverse == False:

#Convert spectral type string to number
        if isinstance(spt,str):
            sptn = typeToNum(spt, uncertainty=unc)
        elif isinstance(spt,int) or isinstance(spt,float):
            sptn = copy.deepcopy(spt)
        else:
            raise ValueError('\nInput spectral type {} must be a string, float or int'.format(spt))

# polynomial method
        if SPT_BC_RELATIONS[ref]['method'] == 'polynomial':
            rng = SPT_BC_RELATIONS[ref]['filters'][filt]['range']
            if (rng[0] <= sptn <= rng[1]):
                bc = numpy.polyval(SPT_BC_RELATIONS[ref]['filters'][filt]['coeff'], sptn-SPT_BC_RELATIONS[ref]['sptoffset'])
                bc_error = SPT_BC_RELATIONS[ref]['filters'][filt]['fitunc']
                if unc > 0.:
                    vals = numpy.polyval(SPT_BC_RELATIONS[ref]['filters'][filt]['coeff'], numpy.random.normal(sptn - SPT_BC_RELATIONS[ref]['sptoffset'], unc, nsamples))
                    bc_error = (numpy.nanstd(vals)**2+SPT_BC_RELATIONS[ref]['filters'][filt]['fitunc']**2)**0.5
                return bc, bc_error
            else:
                if verbose: sys.stderr.write('\nSpectral type {} is out of range for {}'.format(typeToNum(sptn),refstring))
                return numpy.nan, numpy.nan

# interpolation method
        elif SPT_BC_RELATIONS[ref]['method'] == 'interpolate':
            rng = [numpy.nanmin(SPT_BC_RELATIONS[ref]['filters'][filt]['spt']),numpy.nanmax(SPT_BC_RELATIONS[ref]['filters'][filt]['spt'])]
            if (rng[0] <= sptn <= rng[1]):
                f = interp1d(SPT_BC_RELATIONS[ref]['filters'][filt]['spt'],SPT_BC_RELATIONS[ref]['filters'][filt]['values'])
                fe = interp1d(SPT_BC_RELATIONS[ref]['filters'][filt]['spt'],SPT_BC_RELATIONS[ref]['filters'][filt]['rms'])
                bc = float(f(sptn))
                bc_error = float(fe(sptn))
                if unc > 0.:
                    vals = f(numpy.random.normal(sptn, unc, nsamples))
                    bc_error = (numpy.nanstd(vals)**2+bc_error**2)**0.5
                return bc, bc_error
            else:
                if verbose: sys.stderr.write('\nSpectral type {} is out of range for {}'.format(typeToNum(sptn),refstring))
                return numpy.nan, numpy.nan
        else:
            raise ValueError('Unknown method {} for {}'.format(SPT_BC_RELATIONS[ref]['method'],refstring))

# reverse approach: BC -> SpT
    else:
        if not isinstance(spt,float): raise ValueError('Running this in reverse you need to provide a BC value instead of {}'.format(spt))
        if SPT_BC_RELATIONS[ref]['method'] == 'polynomial':
            rng = SPT_BC_RELATIONS[ref]['filters'][filt]['range']
            x = numpy.linspace(rng[0],rng[1],nsamples)
            y = numpy.polyval(SPT_BC_RELATIONS[ref]['filters'][filt]['coeff'], x-SPT_BC_RELATIONS[ref]['sptoffset'])
            f = interp1d(y,x)
            try:
                bc = float(f(spt))
            except:
                if verbose: print('\nBC value {} is outside the range expected for {}'.format(spt,refstring))
                return numpy.nan, numpy.nan
            vals = []
            for i in range(nsamples):
                ye = y+numpy.random.normal(0,SPT_BC_RELATIONS[ref]['filters'][filt]['fitunc'])
                f = interp1d(ye,x)
                try:
                    vals.append(f(numpy.random.normal(spt,unc)))
                except:
                    pass
            bc_error = numpy.nanstd(vals)
            return bc, bc_error
        elif SPT_BC_RELATIONS[ref]['method'] == 'interpolate':
            f = interp1d(SPT_BC_RELATIONS[ref]['filters'][filt]['bc'],SPT_BC_RELATIONS[ref]['filters'][filt]['spt'])
            try:
                bc = f(spt)
            except:
                if verbose: print('\nBC value {} is outside the range expected for {}'.format(spt,refstring))
                return numpy.nan, numpy.nan
            vals = []
            for i in range(nsamples):
                y = numpy.random.normal(SPT_BC_RELATIONS[ref]['filters'][filt]['bc'],SPT_BC_RELATIONS[ref]['filters'][filt]['rms'])
                f = interp1d(y,SPT_BC_RELATIONS[ref]['filters'][filt]['spt'])
                try:
                    vals.append(f(numpy.random.normal(spt,unc)))
                except:
                    pass
            bc_error = numpy.nanstd(vals)
            return bc, bc_error
        else:
            raise ValueError('Unknown method {} for {}'.format(SPT_BC_RELATIONS[ref]['method'],refstring))


def typeToChi(spt,spt_e = 0.,set='douglas2014',nsamples=100):
    '''
    :Purpose: 

        Computes the chi correction factor to convert from Halpha equivalent width to LHa/Lbol
        TO BE COMPLETED

    :Required Inputs: 

        :param spt: spectral type, either string or numerical

    :Optional Inputs: 

        :param spt_e=0.: uncertainty in spectral type; if zero, not monte carlo is done
        :param set='douglas2014': reference for parameters to use for calculatino
        :param nsamples=100: number of Monte Carlo samples for error computation

    :Output: 

        Estimated chi correction factor and uncertainty

    :Example:

        TBD

    '''
    pass

def magToMass(magnitude,set='mann2019',mag_e = 0,nsamples=100):
    '''
    :Purpose: 

        Computes the mass based on a given absolute magnitude and magnitude/mass relation
        TO BE COMPLETED

    :Required Inputs: 

        :param mag: magnitude to use in computation; should be in the filter defined by the relation

    :Optional Inputs: 

        :param set='douglas2014': reference for parameters to use for calculation
        :param mag_e=0.: uncertainty in magnitude; if zero, not monte carlo is done
        :param nsamples=100: number of Monte Carlo samples for error computation

    :Output: 

        Estimated mass in solar mass units and uncertainty

    :Example:

        TBD

    '''
    pass

#####################################################
######   SUMMARIES OF EMPIRICAL RELATIONS   #########
#####################################################



def metallicity(sp,nsamples=100,ref='rojas',output='all',verbose=False,**kwargs):
    '''
    :Purpose: 

        THIS IS OUT OF DATE

        Metallicity measurement using Na I and Ca I lines and H2O-K2 index as described in 
        `Rojas-Ayala et al.(2012) <http://adsabs.harvard.edu/abs/2012ApJ...748...93R>`_
        Requires higher resolution data than SpeX prism mode. 
        Includes calls to measureIndexSet()_ and measureEWSet()_

    :Required Inputs: 

        :param sp: Spectrum class object, which should contain wave, flux and noise array elements

    :Optional Inputs: 

        :param nsamples=100: number of Monte Carlo samples for error computation

    :Output: 

        Estimated metallicity and uncertainty

    :Example:

        TBD

    .. _measureIndexSet() : api.html#splat.core.measureIndexSet
    .. _measureEWSet() : api.html#splat.core.measureEWSet

    '''

    allowed_refs = ['rojas2012','terrien2012','mann2013','mann2014','newton2014']

    if 'rojas' in ref.lower():
        reference = 'Rojas et al. (2012)'
        bibcode = '2012ApJ...748...93R'
        result = {'reference': reference, 'bibcode': bibcode, 'measures': {}}
        if verbose==True: print('Computing [Fe/H] from {} ({}), valid for M0-M5, ??? < [Fe/H] < ???'.format(reference,bibcode))

        coeff_feh = [-1.039,0.092,0.119]
        coeff_feh_e = [0.17,0.023,0.033]
        feh_unc = 0.100
        coeff_mh = [-0.731,0.066,0.083]
        coeff_mh_e = [0.12,0.016,0.023]
        mh_unc = 0.100

        ews = splat.measureEWSet(sp,ref='rojas2012')
        nai,nai_e = ews['nai']['ew'].to(u.Angstrom).value,ews['nai']['ew_unc'].to(u.Angstrom).value
        cai,cai_e = ews['cai']['ew'].to(u.Angstrom).value,ews['cai']['ew_unc'].to(u.Angstrom).value
        h2ok2,h2ok2_e = splat.measureIndexSet(sp, ref='rojas')['H2O-K2']

        feh = coeff_feh[0]+coeff_feh[1]*nai/h2ok2+coeff_feh[2]*cai/h2ok2
        feh_e = numpy.sqrt(feh_unc**2+numpy.nanstd(
            coeff_feh[0]+\
            coeff_feh[1]*(numpy.random.normal(nai,nai_e,nsamples)/numpy.random.normal(h2ok2,h2ok2_e,nsamples))+\
            coeff_feh[2]*(numpy.random.normal(cai,cai_e,nsamples)/numpy.random.normal(h2ok2,h2ok2_e,nsamples)))**2)
        result['measures']['feh'] = [feh,feh_e]
        if verbose==True: print('K-band [Fe/H] = {:.2f}+/-{:.2f}'.format(feh,feh_e))

        mh = coeff_mh[0]+coeff_mh[1]*nai/h2ok2+coeff_mh[2]*cai/h2ok2
        mh_e = numpy.sqrt(feh_unc**2+numpy.nanstd(
            coeff_mh[0]+\
            coeff_mh[1]*(numpy.random.normal(nai,nai_e,nsamples)/numpy.random.normal(h2ok2,h2ok2_e,nsamples))+\
            coeff_mh[2]*(numpy.random.normal(cai,cai_e,nsamples)/numpy.random.normal(h2ok2,h2ok2_e,nsamples)))**2)
        result['measures']['mh'] = [mh,mh_e]
        if verbose==True: print('K-band [M/H] = {:.2f}+/-{:.2f}'.format(mh,mh_e))


# Terrien et al. 2012
    elif 'terrien' in ref.lower():
        reference = 'Terrrien et al. (2012)'
        bibcode = '2012ApJ...747L..38T'
        result = {'reference': reference, 'bibcode': bibcode, 'measures': {}}
        if verbose==True: print('Computing [Fe/H] from {} ({}), valid for M0--M5, -0.25 < [Fe/H] < +0.3'.format(reference,bibcode))

        coeff_feh_h = [0.340,0.407,0.436,-1.485]
        feh_k_unc = 0.12
        coeff_feh_k = [0.132,0.083,-0.403,-0.616]
        feh_h_unc = 0.12

        ews = splat.measureEWSet(sp,ref='terrien2012')
        caih = ews['cai-h1']['ew'].to(u.Angstrom).value+ews['cai-h2']['ew'].to(u.Angstrom).value
        caih_e = numpy.sqrt((ews['cai-h1']['ew_unc'].to(u.Angstrom).value)**2+(ews['cai-h1']['ew_unc'].to(u.Angstrom).value)**2)
        ind = splat.measureIndexSet(sp, ref='covey')
        h2oh,h2oh_e =ind['H2O-H']
        h2ok,h2ok_e = ind['H2O-K']

        feh_k = coeff_feh_k[0]*ews['nai-k']['ew'].to(u.Angstrom).value+\
            coeff_feh_k[1]*ews['cai-k']['ew'].to(u.Angstrom).value+\
            coeff_feh_k[2]*ind['H2O-K'][0]+coeff_feh_k[3]
        feh_k_e = numpy.sqrt(feh_k_unc**2+numpy.nanstd(
            coeff_feh_k[0]*numpy.random.normal(ews['nai-k']['ew'].to(u.Angstrom).value,ews['nai-k']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_k[1]*numpy.random.normal(ews['cai-k']['ew'].to(u.Angstrom).value,ews['cai-k']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_k[2]*numpy.random.normal(ind['H2O-K'][0],ind['H2O-K'][1],nsamples)+coeff_feh_k[3])**2)
        result['measures']['feh_k'] = [feh_k,feh_k_e]
        if verbose==True: print('K-band [Fe/H] = {:.2f}+/-{:.2f}'.format(feh_k,feh_k_e))

        feh_h = coeff_feh_h[0]*caih+\
            coeff_feh_h[1]*ews['ki-h']['ew'].to(u.Angstrom).value+\
            coeff_feh_h[2]*ind['H2O-H'][0]+coeff_feh_h[3]
        feh_h_e = numpy.sqrt(feh_h_unc**2+numpy.nanstd(
            coeff_feh_h[0]*numpy.random.normal(caih,caih_e,nsamples)+\
            coeff_feh_h[1]*numpy.random.normal(ews['ki-h']['ew'].to(u.Angstrom).value,ews['ki-h']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_h[2]*numpy.random.normal(ind['H2O-H'][0],ind['H2O-H'][1],nsamples)+coeff_feh_h[3])**2)
        result['measures']['feh_h'] = [feh_h,feh_h_e]
        if verbose==True: print('H-band [Fe/H] = {:.2f}+/-{:.2f}'.format(feh_h,feh_h_e))


# Mann et al. 2013
    elif 'mann' in ref.lower() and '14' not in ref.lower():
        reference = 'Mann et al. (2013)'
        bibcode = '2013AJ....145...52M'
        result = {'reference': reference, 'bibcode': bibcode, 'measures': {}}
        if verbose==True: print('Computing [Fe/H] from {} ({}), valid for K5--M5, -1.04 < [Fe/H] < +0.56'.format(reference,bibcode))

        coeff_feh_j = [0.29,0.21,0.26,-0.26,-0.190,-1.03]
        feh_unc_j = 0.07
        coeff_mh_j = [0.32,0.46,0.076,1.213,-1.97]
        mh_unc_j = 0.08
        coeff_feh_h = [0.40,0.51,-0.28,-1.460,0.71]
        feh_unc_h = 0.07
        coeff_mh_h = [0.38,0.40,0.41,0.194,-0.76]
        mh_unc_h = 0.06
        coeff_feh_k = [0.19,0.069,0.083,0.218,-1.55]
        feh_unc_k = 0.06
        coeff_mh_k = [0.12,0.086,0.13,0.245,-1.18]
        mh_unc_k = 0.05

        ews = splat.measureEWSet(sp,ref='mann2013')
        ind = splat.measureIndexSet(sp, ref='mann')
        h2oj,h2oj_e =ind['H2O-J']
        ind = splat.measureIndexSet(sp, ref='covey')
        h2oh,h2oh_e =ind['H2O-H']
        ind = splat.measureIndexSet(sp, ref='rojas')
        h2ok,h2ok_e =ind['H2O-K2']


        feh_j = coeff_feh_j[0]*ews['f10']['ew'].to(u.Angstrom).value+\
            coeff_feh_j[1]*ews['f09']['ew'].to(u.Angstrom).value+\
            coeff_feh_j[2]*ews['f12']['ew'].to(u.Angstrom).value+\
            coeff_feh_j[3]*ews['f13']['ew'].to(u.Angstrom).value+\
            coeff_feh_j[4]*h2oj+coeff_feh_j[5]
        feh_j_e = numpy.sqrt(feh_unc_j**2+numpy.nanstd(\
            coeff_feh_j[0]*numpy.random.normal(ews['f10']['ew'].to(u.Angstrom).value,ews['f10']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_j[1]*numpy.random.normal(ews['f09']['ew'].to(u.Angstrom).value,ews['f09']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_j[2]*numpy.random.normal(ews['f12']['ew'].to(u.Angstrom).value,ews['f12']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_j[3]*numpy.random.normal(ews['f13']['ew'].to(u.Angstrom).value,ews['f13']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_j[4]*numpy.random.normal(h2oj,h2oj_e,nsamples)+coeff_feh_j[5])**2)
        result['measures']['feh_j'] = [feh_j,feh_j_e]
        if verbose==True: print('J-band [Fe/H] = {:.2f}+/-{:.2f}'.format(feh_j,feh_j_e))

        mh_j = coeff_mh_j[0]*ews['f10']['ew'].to(u.Angstrom).value+\
            coeff_mh_j[1]*ews['f11']['ew'].to(u.Angstrom).value+\
            coeff_mh_j[2]*ews['f09']['ew'].to(u.Angstrom).value+\
            coeff_mh_j[3]*h2oj+coeff_mh_j[4]
        mh_j_e = numpy.sqrt(mh_unc_j**2+numpy.nanstd(\
            coeff_mh_j[0]*numpy.random.normal(ews['f10']['ew'].to(u.Angstrom).value,ews['f10']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_mh_j[1]*numpy.random.normal(ews['f11']['ew'].to(u.Angstrom).value,ews['f11']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_mh_j[2]*numpy.random.normal(ews['f09']['ew'].to(u.Angstrom).value,ews['f09']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_mh_j[3]*numpy.random.normal(h2oj,h2oj_e,nsamples)+coeff_mh_j[4])**2)
        result['measures']['mh_j'] = [mh_j,mh_j_e]
        if verbose==True: print('J-band [M/H] = {:.2f}+/-{:.2f}'.format(mh_j,mh_j_e))

        feh_h = coeff_feh_h[0]*ews['f17']['ew'].to(u.Angstrom).value+\
            coeff_feh_h[1]*ews['f14']['ew'].to(u.Angstrom).value+\
            coeff_feh_h[2]*ews['f18']['ew'].to(u.Angstrom).value+\
            coeff_feh_h[3]*h2oh+coeff_feh_h[4]
        feh_h_e = numpy.sqrt(feh_unc_h**2+numpy.nanstd(\
            coeff_feh_h[0]*numpy.random.normal(ews['f17']['ew'].to(u.Angstrom).value,ews['f17']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_h[1]*numpy.random.normal(ews['f14']['ew'].to(u.Angstrom).value,ews['f14']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_h[2]*numpy.random.normal(ews['f18']['ew'].to(u.Angstrom).value,ews['f18']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_h[3]*numpy.random.normal(h2oh,h2oh_e,nsamples)+coeff_feh_h[4])**2)
        result['measures']['feh_h'] = [feh_h,feh_h_e]
        if verbose==True: print('H-band [Fe/H] = {:.2f}+/-{:.2f}'.format(feh_h,feh_h_e))

        mh_h = coeff_mh_h[0]*ews['f17']['ew'].to(u.Angstrom).value+\
            coeff_mh_h[1]*ews['f16']['ew'].to(u.Angstrom).value+\
            coeff_mh_h[2]*ews['f15']['ew'].to(u.Angstrom).value+\
            coeff_mh_h[3]*h2oh+coeff_mh_h[4]
        mh_h_e = numpy.sqrt(mh_unc_j**2+numpy.nanstd(\
            coeff_mh_h[0]*numpy.random.normal(ews['f17']['ew'].to(u.Angstrom).value,ews['f17']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_mh_h[1]*numpy.random.normal(ews['f16']['ew'].to(u.Angstrom).value,ews['f16']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_mh_h[2]*numpy.random.normal(ews['f15']['ew'].to(u.Angstrom).value,ews['f15']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_mh_h[3]*numpy.random.normal(h2oh,h2oh_e,nsamples)+coeff_mh_h[4])**2)
        result['measures']['mh_h'] = [mh_h,mh_h_e]
        if verbose==True: print('H-band [M/H] = {:.2f}+/-{:.2f}'.format(mh_h,mh_h_e))

        feh_k = coeff_feh_k[0]*ews['f19']['ew'].to(u.Angstrom).value+\
            coeff_feh_k[1]*ews['f22']['ew'].to(u.Angstrom).value+\
            coeff_feh_k[2]*ews['f20']['ew'].to(u.Angstrom).value+\
            coeff_feh_k[3]*h2ok+coeff_feh_k[4]
        feh_k_e = numpy.sqrt(feh_unc_h**2+numpy.nanstd(\
            coeff_feh_k[0]*numpy.random.normal(ews['f19']['ew'].to(u.Angstrom).value,ews['f19']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_k[1]*numpy.random.normal(ews['f22']['ew'].to(u.Angstrom).value,ews['f22']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_k[2]*numpy.random.normal(ews['f20']['ew'].to(u.Angstrom).value,ews['f20']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh_k[3]*numpy.random.normal(h2ok,h2ok_e,nsamples)+coeff_feh_k[4])**2)
        result['measures']['feh_k'] = [feh_k,feh_k_e]
        if verbose==True: print('K-band [Fe/H] = {:.2f}+/-{:.2f}'.format(feh_k,feh_k_e))

        mh_k = coeff_mh_k[0]*ews['f19']['ew'].to(u.Angstrom).value+\
            coeff_mh_k[1]*ews['f22']['ew'].to(u.Angstrom).value+\
            coeff_mh_k[2]*ews['f21']['ew'].to(u.Angstrom).value+\
            coeff_mh_k[3]*h2ok+coeff_mh_k[4]
        mh_k_e = numpy.sqrt(mh_unc_j**2+numpy.nanstd(\
            coeff_mh_k[0]*numpy.random.normal(ews['f19']['ew'].to(u.Angstrom).value,ews['f19']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_mh_k[1]*numpy.random.normal(ews['f22']['ew'].to(u.Angstrom).value,ews['f22']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_mh_k[2]*numpy.random.normal(ews['f21']['ew'].to(u.Angstrom).value,ews['f21']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_mh_k[3]*numpy.random.normal(h2ok,h2ok_e,nsamples)+coeff_mh_k[4])**2)
        result['measures']['mh_k'] = [mh_k,mh_k_e]
        if verbose==True: print('K-band [M/H] = {:.2f}+/-{:.2f}'.format(mh_k,mh_k_e))


# Newton et al. 2014
    elif 'newton' in ref.lower():
        reference = 'Newton et al. (2014)'
        bibcode = '2014AJ....147...20N'
        result = {'reference': reference, 'bibcode': bibcode, 'measures': {}}
        if verbose==True: print('Computing [Fe/H] from {} ({}), valid for M1--M5, -1.0 < [Fe/H] < +0.35'.format(reference,bibcode))

        coeff_feh = [0.596,-0.0392,-1.96]
        feh_unc = 0.12

        ews = splat.measureEWSet(sp,ref='rojas2012')

        feh = coeff_feh[0]*ews['nai']['ew'].to(u.Angstrom).value+\
            coeff_feh[1]*ews['nai']['ew'].to(u.Angstrom).value**2+coeff_feh[2]
        nais = numpy.random.normal(ews['nai']['ew'].to(u.Angstrom).value,ews['nai']['ew_unc'].to(u.Angstrom).value,nsamples)
        feh_e = numpy.sqrt(feh_unc**2+numpy.nanstd(coeff_feh[0]*nais+coeff_feh[1]*nais**2+coeff_feh[2])**2)
        result['measures']['feh'] = [feh,feh_e]
        if verbose==True: print('K-band [Fe/H] = {:.2f}+/-{:.2f}'.format(feh,feh_e))

# Mann et al. 2014
    elif 'mann' in ref.lower() and '14' in ref.lower():
        reference = 'Mann et al. (2014)'
        bibcode = '2014AJ....147..160M'
        result = {'reference': reference, 'bibcode': bibcode, 'measures': {}}
        if verbose==True: print('Computing [Fe/H] from {} ({}), valid for M4.5--M9.5, -0.58 < [Fe/H] < +0.56'.format(reference,bibcode))

        coeff_feh = [0.131,0.210,-3.07,1.341]
        feh_unc = 0.07

        ews = splat.measureEWSet(sp,ref='mann2014')
        ind = splat.measureIndexSet(sp, ref='rojas')
#        h2ok,h2ok_e =ind['H2O-K2']

        feh = coeff_feh[0]*ews['nai']['ew'].to(u.Angstrom).value+\
            coeff_feh[1]*ews['cai']['ew'].to(u.Angstrom).value+\
            coeff_feh[2]*ind['H2O-K2'][0]+coeff_feh[3]
        feh_e = numpy.sqrt(feh_unc**2+numpy.nanstd(\
            coeff_feh[0]*numpy.random.normal(ews['nai']['ew'].to(u.Angstrom).value,ews['nai']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh[1]*numpy.random.normal(ews['cai']['ew'].to(u.Angstrom).value,ews['cai']['ew_unc'].to(u.Angstrom).value,nsamples)+\
            coeff_feh[2]*numpy.random.normal(ind['H2O-K2'][0],ind['H2O-K2'][1],nsamples)+coeff_feh[3])**2)
        result['measures']['feh'] = [feh,feh_e]
        if verbose==True: print('K-band [Fe/H] = {:.2f}+/-{:.2f}'.format(feh,feh_e))

    else:
        raise ValueError('Reference {} not recognized as one of the metallicity relations currently in SPLAT; try {}'.format(ref,allowed_refs))

# return
    if output=='feh': 
        if 'rojas' in ref.lower() or 'newton' in ref.lower() or ('mann' in ref.lower() and '14' in ref.lower()): 
            return result['measures']['feh']
        elif 'mann' in ref.lower(): 
            if verbose==True: print('Returning {} K-band [Fe/H] measurement and uncertainty'.format(result['reference']))
            return result['measures']['feh_k']
        elif 'terrien' in ref.lower(): 
            if verbose==True: print('Returning {} H-band [Fe/H] measurement and uncertainty'.format(result['reference']))
            return result['measures']['feh_h']
        else: return result

    elif output=='mh': 
        if 'rojas' in ref.lower(): return result['measures']['mh']
        elif 'mann' in ref.lower() and '14' not in ref.lower(): 
            if verbose==True: print('Returning {} K-band [M/H] measurement and uncertainty'.format(result['reference']))
            return result['measures']['mh_k']
        else: return result

    else: return result
        

#####################################################
######   SUMMARIES OF EMPIRICAL RELATIONS   #########
#####################################################


def info_indices(*args):
    if len(args)>0:
        if isinstance(args[0],list) == True: sets = args[0]
        else: sets = list(args)
        print('\nChecking index sets {}:\n'.format(sets))
    else:
        sets = list(INDEX_SETS.keys())
        print('\nIndex sets currently installed in SPLAT:\n')
    flag=0

    for ref in sets:
        refch = checkDict(ref,INDEX_SETS)
        if refch!=False: 
            flag = 1
            print('{}: {}'.format(refch,INDEX_SETS[refch]['bibcode']))
            for i in INDEX_SETS[refch]['indices']:
                if INDEX_SETS[refch]['indices'][i]['method']=='value':
                    print('\t{}: [{}-{}]'.format(
                        i,(INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        ))
                elif INDEX_SETS[refch]['indices'][i]['method']=='ratio':
                    print('\t{}: [{}-{}] / [{}-{}]'.format(
                        i,(INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        ))
                elif INDEX_SETS[refch]['indices'][i]['method']=='doubleratio':
                    print('\t{}: [{}-{}]/[{}-{}] / [{}-{}]/[{}-{}]'.format(
                        i,(INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[1],
                        ))
                elif INDEX_SETS[refch]['indices'][i]['method']=='line':
                    print('\t{}: 0.5 ([{}-{}] + [{}-{}]) / [{}-{}]'.format(
                        i,(INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[1]
                        ))
                elif INDEX_SETS[refch]['indices'][i]['method']=='sumnum':
                    print('\t{}: ([{}-{}] + [{}-{}]) / [{}-{}]'.format(
                        i,(INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[1]
                        ))
                elif INDEX_SETS[refch]['indices'][i]['method']=='inverse_line':
                    print('\t{}: 2 [{}-{}] / ([{}-{}] + [{}-{}])'.format(
                        i,(INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[1]
                        ))
                elif INDEX_SETS[refch]['indices'][i]['method']=='sumdenum':
                    print('\t{}: [{}-{}] / ([{}-{}] + [{}-{}])'.format(
                        i,(INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[1]
                        ))
                elif INDEX_SETS[refch]['indices'][i]['method']=='change':
                    print('\t{}: 2 ([{}-{}]-[{}-{}]) / ([{}-{}] + [{}-{}])'.format(
                        i,(INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        ))
                elif INDEX_SETS[refch]['indices'][i]['method']=='allers':
                    print('\t{}: ([{}-{}] - [{}-{}]) / ([{}-{}] - [{}-{}]) + ([{}-{}] - [{}-{}]) / ([{}-{}] - [{}-{}]) + '.format(
                        i,(INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][0].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][2].value)[1],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[0],
                        (INDEX_SETS[refch]['indices'][i]['ranges'][1].value)[1],
                        ))
                else:
                    print('\t{}: dunno'.format(i))
    if flag==0: print('Index set(s) {} are not present in SPLAT'.format(sets))
    return


def info_classifyByIndex(*args):
    if len(args)>0:
        if isinstance(args[0],list) == True: sets = args[0]
        else: sets = list(args)
        print('\nChecking index sets {}:\n'.format(sets))
    else:
        sets = list(INDEX_SETS.keys())
        print('\nIndex sets currently installed in SPLAT:\n')
    flag=0

    for ref in sets:
        refch = checkDict(ref,INDEX_SETS)
        if refch!=False: 
            flag = 1