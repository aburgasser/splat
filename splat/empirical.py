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
#from astropy import units as u            # standard units
#from astropy import constants as const        # physical constants in SI units
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# splat functions
from splat.initialize import *
from splat.utilities import *
from splat.photometry import filterMag
from splat.core import classifyByIndex, Spectrum
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



def typeToColor(spt,color, **kwargs):
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

#Keywords
    nsamples = kwargs.get('nsamples', 100)
    ref = kwargs.get('ref', 'skrzypek')
    ref = kwargs.get('set', ref)
    unc = kwargs.get('unc', 0.)

# what empirical relations are used:
    empirical_sets = {
        'skrzypek': {
            'reference': 'Skrzypek et al. (2015)',
            'bibcode': '2015A%26A...574A..78S',
            'rng': [15,38],
            'filters': ['i','z','y','j','h','k','w1','w2'],
            'scatter': 0.07,
            'values': { 
                'i-z': [0.91,1.45,1.77,1.93,1.99,2.01,2.02,2.04,2.1,2.2,2.33,2.51,2.71,2.93,3.15,3.36,3.55,3.7,3.82,3.9,3.95,3.98,4.01,4.08], \
                'z-y': [0.47,0.6,0.7,0.77,0.82,0.86,0.88,0.9,0.92,0.94,0.97,1.0,1.04,1.09,1.16,1.23,1.33,1.43,1.55,1.68,1.81,1.96,2.11,2.26], \
                'y-j': [0.55,0.67,0.78,0.87,0.96,1.04,1.11,1.18,1.23,1.27,1.31,1.33,1.35,1.21,1.2,1.19,1.19,1.18,1.18,1.17,1.16,1.16,1.15,1.15], \
                'j-h': [0.45,0.53,0.56,0.58,0.6,0.63,0.67,0.73,0.79,0.86,0.91,0.96,0.97,0.96,0.9,0.8,0.65,0.46,0.25,0.02,-0.19,-0.35,-0.43,-0.36], \
                'h-k': [0.32,0.39,0.44,0.47,0.51,0.54,0.58,0.63,0.67,0.71,0.74,0.75,0.75,0.71,0.65,0.56,0.45,0.31,0.16,0.01,-0.11,-0.19,-0.2,-0.09], \
                'k-w1': [0.11,0.22,0.25,0.26,0.27,0.29,0.33,0.4,0.48,0.56,0.65,0.72,0.77,0.79,0.79,0.76,0.71,0.65,0.59,0.55,0.54,0.59,0.7,0.9], \
                'w1-w2': [0.17,0.21,0.24,0.26,0.27,0.27,0.28,0.28,0.29,0.3,0.32,0.36,0.41,0.48,0.57,0.68,0.82,0.99,1.19,1.43,1.7,2.02,2.38,2.79]}
        }
    }

#Convert spectral type string to number
    if isinstance(spt,str):
        spt = typeToNum(spt, uncertainty=unc)
    else:
        spt = copy.deepcopy(spt)

    if ref.lower() in list(empirical_sets.keys()):
        reference = empirical_sets[ref.lower()]['reference']
        rng = empirical_sets[ref.lower()]['rng']
        filters = empirical_sets[ref.lower()]['filters']
        values = empirical_sets[ref.lower()]['values']
        scatter = empirical_sets[ref.lower()]['scatter']

    else:
        sys.stderr.write('\nColor set from {} has not be intergrated into SPLAT\n\n'.format(ref))
        return numpy.nan, numpy.nan
    if kwargs.get('verbose',False):
        print('\nUsing the SpT/color trends from {}\n'.format(reference))

# spectral type array
    if (rng[0] <= spt <= rng[1]):

# fill in extra colors - a little inefficient right now  
        if color.lower() not in list(values.keys()):
#            basecolors = 
            c1 = (color.lower()).split('-')[0]


            for i in numpy.arange(len(list(tmpval.keys()))):
                for a in list(tmpval.keys()):
                    for b in list(tmpval.keys()):
                        f1 = a.split('-')
                        f2 = b.split('-')
                        if f1[-1] == f2[0]:
                            k = '{}-{}'.format(f1[0],f2[-1])
                            if k not in list(values.keys()):
                                values[k] = [sum(x) for x in zip(tmpval[a], tmpval[b])]

        if color.lower() in list(values.keys()):
            f = interp1d(numpy.arange(rng[0],rng[1]+1),values[color.lower()],bounds_error=False,fill_value=0.)
            if (unc > 0.):
                vals = [f(x) for x in numpy.random.normal(spt, unc, nsamples)]
                return float(f(spt)), (numpy.nanstd(vals)**2+scatter**2)**0.5
            else:
                return float(f(spt)), scatter
        else:
            sys.stderr.write('\n Color {} is not in reference set for {}\n\n'.format(color,reference))
            return numpy.nan, numpy.nan

    else:
        sys.stderr.write('\n Spectral type {} is outside the range for reference set {}\n\n'.format(typeToNum(spt),reference))
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


def typeToMag(spt, filt, unc=0.,ref='filippazzo2015',verbose=False,nsamples=100,**kwargs):
    """
    :Purpose: 

    Takes a spectral type and a filter, and returns the expected absolute magnitude based on empirical relations

    :Required Inputs: 

        :param spt: string or integer of the spectral type
        :param filter: filter for which to retrieve absolute magnitude, which must be defined for the given reference set.
        You can check what filters are available by printing splat.ABSMAG_SETS[reference]['filters'].keys(), where reference is, e.g., 'filippazzo2015'

    :Optional Inputs: 

        :param reference: Abs Mag/SpT relation used to compute the absolute magnitude (also 'ref' and 'set'). These are defined in splat.ABSMAG_SETS and are currently as follows:

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

    sptoffset = ABSMAG_SETS[ref]['sptoffset']
    coeff = ABSMAG_SETS[ref]['filters'][filt]['coeff']
    rng = ABSMAG_SETS[ref]['filters'][filt]['range']
    fitunc = ABSMAG_SETS[ref]['filters'][filt]['fitunc']
    refstring = 'Absolute {}/SpT relation from {}'.format(filt,shortRef(ABSMAG_SETS[ref]['bibcode']))
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



def typeToTeff(inp, uncertainty=0.001, ref='stephens',nsamples=100, reverse=False, **kwargs):
    '''
    :Purpose: 

    Returns an effective temperature (Teff) and its uncertainty for a given spectral type (or vice versa) based on an empirical relation

    :Required Inputs:

        :param inp: either the spectral type or effective temperature. 
        If a spectral type, this can be a number or a string from 0 (K0) and 49.0 (Y9).
        If a temperature, this is assumed to be in Kelvin; you must set reverse=True for this option

    :Optional Inputs:

        :param uncertainty: uncertainty of spectral type/temperature (default = 0.001; also 'unc', 'spt_e')
        :param ref: Teff/SpT relation used to compute the effective temperature (also 'set'). Options are:

        - *stephens* (default): Teff/SpT relation from `Stephens et al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...702..154S>`_.
          Allowed spectral type range is M6 to T8 and uses alternate coefficients for L3 to T8.
        - *golimowski*: Teff/SpT relation from `Golimowski et al. (2004) <http://adsabs.harvard.edu/abs/2004AJ....127.3516G>`_.
          Allowed spectral type range is M6 to T8.
        - *looper*: Teff/SpT relation from `Looper et al. (2008) <http://adsabs.harvard.edu/abs/2008ApJ...685.1183L>`_.
          Allowed spectral type range is L0 to T8.
        - *marocco*: Teff/SpT relation from `Marocco et al. (2013) <http://adsabs.harvard.edu/abs/2013AJ....146..161M>`_.
          Allowed spectral type range is M7 to T8.
        - *filippazzo*: Teff/SpT relation from Filippazzo et al. (2015). Allowed spectral type range is M6 to T9.

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
# keywords
    verbose = kwargs.get('verbose',False)
    nsamples = kwargs.get('nsamples',100)
    unc = kwargs.get('uncertainty',0.001)
    unc = kwargs.get('unc',unc)
    unc = kwargs.get('spt_e',unc)
    ref = kwargs.get('ref','stephens2009')
    ref = kwargs.get('reference',ref)
    ref = kwargs.get('set',ref)
    ref = kwargs.get('method',ref)

# convert spectral type string to number
    if (type(inp) == str):
        spt = typeToNum(inp,uncertainty=unc)
    else:
        spt = copy.deepcopy(inp)

#    if spt < 20. and 'marocco' not in ref.lower():
#        ref='stephens2009'

# choose among possible options

# Golimowski et al. (2004, AJ, 127, 3516)
    if ('golimowski' in ref.lower()):
        reference = 'Teff/SpT relation from Golimowski et al. (2004)'
        sptoffset = 10.
        coeff = [9.5373e-4,-9.8598e-2,4.0323,-8.3099e1,9.0951e2,-5.1287e3,1.4322e4]
        sptrange = [16.,38.]
        fitunc = 124.

# Looper et al. (2008, ApJ, 685, 1183)
    elif ('looper' in ref.lower()):
        reference = 'Teff/SpT relation from Looper et al. (2008)'
        sptoffset = 20.
        coeff = [9.084e-4,-4.255e-2,6.414e-1,-3.101,1.950,-108.094,2319.92]
        sptrange = [20.,38.]
        fitunc = 87.

# Stephens et al. (2009, ApJ, 702, 1545); using OPT/IR relation for M6-T8
# plus alternate coefficients for L3-T8
    elif ('stephens' in ref.lower()):
        reference = 'Teff/SpT relation from Stephens et al. (2009)'
        sptoffset = 10.
        coeff = [-0.0025492,0.17667,-4.4727,54.67,-467.26,4400.]
        sptrange = [16.,38.]
        fitunc = 100.
        coeff_alt = [-0.011997,1.2315,-50.472,1031.9,-10560.,44898.]
        range_alt = [23.,38.]

# Marocco et al. (2013, AJ, 146, 161)
    elif ('marocco' in ref.lower()):
        reference = 'Teff/SpT relation from Marocco et al. (2013)'
        sptoffset = 10.
        coeff = [7.4211e-5,-8.43736e-3,3.90319e-1,-9.46896,129.141,-975.953,3561.47,-1613.82]
        sptrange = [17.,38.]
        fitunc = 140.

    elif ('filippazzo' in ref.lower()):
        reference = 'Teff/SpT relation from Filippazzo et al. (2015)'
        sptoffset = 10.
        coeff = [1.546e-4, -1.606e-2, 6.318e-1, -1.191e1, 1.155e2, -7.005e2, 4.747e3]
        sptrange = [16., 39.]
        fitunc = 113.

    elif ('faherty' in ref.lower()):
        sptoffset = 10.
        if kwargs.get('young',False) == True:
            sptrange = [17., 27.]
            reference = 'Teff/SpT young relation from Faherty et al. (2016)'
            coeff = [1.330,-6.68637e1,1.23542e3,-1.00688e4,3.27664e4]
            fitunc = 180.
        elif kwargs.get('young2',False) == True:
            sptrange = [17., 27.]
            reference = 'Teff/SpT young2 relation from Faherty et al. (2016)'
            coeff = [9.106e-4,-1.016e-1,4.578,-1.066e2,1.360e3,-9.183e3,2.795e4]
            fitunc = 198.
        elif kwargs.get('group',False) == True:
            sptrange = [17., 27.]
            reference = 'Teff/SpT group relation from Faherty et al. (2016)'
            coeff = [7.383e0,-3.44522e2,4.87986e3]
            fitunc = 172.
        else:
            sptrange = [17., 38.]
            reference = 'Teff/SpT field relation from Faherty et al. (2016)'
            coeff = [1.546e-4,-1.606e-2,6.318e-1,-1.191e1,1.155e2,-7.005e2,4.747e3]
            fitunc = 113.

    elif ('dupuy' in ref.lower()):
        sptoffset = 10.
        if kwargs.get('saumon',False) == True or kwargs.get('sm08',False) == True:
            reference = 'Teff/SpT relation from Dupuy & Liu (2017) with Saumon & Marley (2008) models'
            coeff = [6.001,-284.52,4544.3]
            sptrange = [21.5, 35.]
            fitunc = 80.
        else:
            reference = 'Teff/SpT relation from Dupuy & Liu (2017) with Lyon models'
            coeff = [4.582,-238.03,4251.0]
            sptrange = [17., 35.]
            fitunc = 90.

    else:
        sys.stderr.write('\nInvalid Teff/SpT relation given ({})\n'.format(ref))
        return numpy.nan, numpy.nan

    if reverse == True:

        teff = spt
        teff_e = unc
        if numpy.min(numpy.polyval(coeff,[r-sptoffset for r in sptrange])) <= teff <= numpy.max(numpy.polyval(coeff,[r-sptoffset for r in sptrange])):
            x = numpy.linspace(sptrange[1]-sptoffset,sptrange[0]-sptoffset,nsamples)
            f = interp1d(numpy.polyval(coeff,x),x,bounds_error=False)
            spto = float(f(teff))+sptoffset
            x = numpy.random.normal(teff,teff_e,nsamples)
            vals = f(x)+sptoffset
            if 'stephens' in ref.lower():
                if numpy.min(numpy.polyval(coeff_alt,[r-sptoffset for r in range_alt])) <= teff <= numpy.max(numpy.polyval(coeff_alt,[r-sptoffset for r in range_alt])):
                    x = numpy.linspace(range_alt[1]-sptoffset,range_alt[0]-sptoffset,nsamples)
                    f = interp1d(numpy.polyval(coeff_alt,x),x,bounds_error=False)
                    spto = float(f(teff))+sptoffset
                    if kwargs.get('string',False) == True: spto = splat.typeToNum(spto)
                    x = numpy.random.normal(teff,teff_e,nsamples)+sptoffset
                    vals = f(x)
# assuming an at least 0.5 spectral type uncertainty
            spto_e = (numpy.nanstd(vals)**2+0.5**2)**0.5
            return spto, spto_e
        else:
            if verbose: sys.stderr.write('\nTeff is out of range for {:s} Teff/SpT relation\n'.format(reference))
            return numpy.nan, numpy.nan

    else:

        if (sptrange[0] <= spt <= sptrange[1]):
            teff = numpy.polyval(coeff,spt-sptoffset)
            x = numpy.random.normal(spt,unc,nsamples)
            x = x[numpy.where(numpy.logical_and(x >= sptrange[0],x <= sptrange[1]))]
            vals = numpy.polyval(coeff,x-sptoffset)
            if ('stephens' in ref.lower()):
                if (range_alt[0] <= spt <= range_alt[1]):
                    teff = numpy.polyval(coeff_alt,spt-sptoffset)
                    vals = numpy.polyval(coeff_alt,x-sptoffset)
    #        teff = numpy.nanmean(vals)
            teff_e = (numpy.nanstd(vals)**2+fitunc**2)**0.5
            return teff, teff_e
        else:
            if verbose: sys.stderr.write('\nSpectral Type is out of range for {:s} Teff/SpT relation\n'.format(reference))
            return numpy.nan, numpy.nan


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
    spabs = Spectrum(wave=w,flux=absfrac)
    return sp*spabs


def typeToLbol(*args,**kwargs):
    return typeToLuminosity(*args,**kwargs)

def typeToLuminosity(spt, unc=0.,ref='filippazzo2015',verbose=False,nsamples=100,reverse=False,**kwargs):
    """
    :Purpose: 

    Takes a spectral type and returns the expected scaled log luminosity (log Lbol/Lsun) based on empirical relations

    :Required Inputs: 

        :param spt: string or integer of the spectral type

    :Optional Inputs: 

        :param reference: log Lbol/SpT relation reference (also 'ref' and 'set'). These are defined in splat.LBOL_SETS and are currently as follows:

            - *filippazzo2015* (default): Lbol/SpT relation from `Filippazzo et al. (2015) <http://adsabs.harvard.edu/abs/2013Sci...341.1492D>`_
              Allowed spectral type range is M6 to T9

        :param unc: uncertainty of ``spt`` (default = 0)
        :param reverse: apply reverse approach: given BC, infer spectral type
        :param nsamples: number of Monte Carlo samples for error computation (default = 100)

    :Output: 
    
        2 element tuple providing the absolute magnitude and its uncertainty

    :Example:
        >>> import splat
        >>> print splat.typeToLuminosity('L3')
    """

#Keywords alternatives
    ref = kwargs.get('reference', ref)
    ref = kwargs.get('set', ref)
    unc = kwargs.get('uncertainty', unc)
    unc = kwargs.get('error', unc)


# check that you can use the proscribed relation and filter
    refcheck = checkEmpiricalRelation(ref,LBOL_SETS,verbose=verbose)
    if refcheck == False: return numpy.nan,numpy.nan
    else: ref=refcheck

    refstring = 'Luminosity/SpT relation for from {}'.format(shortRef(LBOL_SETS[ref]['bibcode']))
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
        if LBOL_SETS[ref]['method'] == 'polynomial':
            rng = LBOL_SETS[ref]['range']
            if (rng[0] <= sptn <= rng[1]):
                lbol = numpy.polyval(LBOL_SETS[ref]['coeff'], sptn-LBOL_SETS[ref]['sptoffset'])
                lbol_error = LBOL_SETS[ref]['fitunc']
                if unc > 0.:
                    vals = numpy.polyval(LBOL_SETS[ref]['coeff'], numpy.random.normal(sptn - LBOL_SETS[ref]['sptoffset'], unc, nsamples))
                    lbol_error = (numpy.nanstd(vals)**2+lbol_error**2)**0.5
                return lbol, lbol_error
            else:
                if verbose: sys.stderr.write('\nSpectral type {} is out of range for {}'.format(typeToNum(sptn),refstring))
                return numpy.nan, numpy.nan

# interpolation method
        elif LBOL_SETS[ref]['method'] == 'interpolate':
            rng = [numpy.nanmin(LBOL_SETS[ref]['spt']),numpy.nanmax(LBOL_SETS[ref]['spt'])]
            if (rng[0] <= sptn <= rng[1]):
                f = interp1d(LBOL_SETS[ref]['spt'],LBOL_SETS[ref]['bc'])
                fe = interp1d(LBOL_SETS[ref]['spt'],LBOL_SETS[ref]['rms'])
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
            raise ValueError('Unknown method {} for {}'.format(LBOL_SETS[ref]['method'],refstring))

# reverse approach: Lbol -> SpT
    else:
        if not isinstance(spt,float): raise ValueError('Running this in reverse you need to provide a log luminosity value instead of {}'.format(spt))
        if LBOL_SETS[ref]['method'] == 'polynomial':
            rng = LBOL_SETS[ref]['range']
            x = numpy.linspace(rng[0],rng[1],nsamples)
            y = numpy.polyval(LBOL_SETS[ref]['coeff'], x-LBOL_SETS[ref]['sptoffset'])
            f = interp1d(y,x)
            try:
                lbol = float(f(spt))
            except:
                if verbose: print('\nlog luminosity value {} is outside the range expected for {}'.format(spt,refstring))
                return numpy.nan, numpy.nan
            vals = []
            for i in range(nsamples):
                ye = y+numpy.random.normal(0,LBOL_SETS[ref]['fitunc'])
                f = interp1d(ye,x)
                try:
                    vals.append(f(numpy.random.normal(spt,unc)))
                except:
                    pass
            lbol_error = numpy.nanstd(vals)
            return lbol, lbol_error
        elif BC_SETS[ref]['method'] == 'interpolate':
            f = interp1d(LBOL_SETS[ref]['bc'],LBOL_SETS[ref]['spt'])
            try:
                lbol = f(spt)
            except:
                if verbose: print('\nlog luminosity value {} is outside the range expected for {}'.format(spt,refstring))
                return numpy.nan, numpy.nan
            vals = []
            for i in range(nsamples):
                y = numpy.random.normal(LBOL_SETS[ref]['bc'],LBOL_SETS[ref]['rms'])
                f = interp1d(y,LBOL_SETS[ref]['spt'])
                try:
                    vals.append(f(numpy.random.normal(spt,unc)))
                except:
                    pass
            lbol_error = numpy.nanstd(vals)
            return lbol, lbol_error
        else:
            raise ValueError('Unknown method {} for {}'.format(LBOL_SETS[ref]['method'],refstring))



def typeToBC(spt, filt, unc=0.,ref='filippazzo2015',verbose=False,nsamples=100,reverse=False,**kwargs):
    """
    :Purpose: 

    Takes a spectral type and a filter, and returns the expected bolometric correction BC = M_bol - M_filter

    :Required Inputs: 

        :param spt: string or integer of the spectral type
        :param filter: filter for which to retrieve absolute magnitude, which must be defined for the given reference set.
        You can check what filters are available by printing splat.ABSMAG_SETS[reference]['filters'].keys(), where reference is, e.g., 'filippazzo2015'

    :Optional Inputs: 

        :param reference: Abs Mag/SpT relation used to compute the absolute magnitude (also 'ref' and 'set'). These are defined in splat.BC_SETS and are currently as follows:

            - *liu2010*: BC/SpT relation from `Liu et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJ...722..311L>`_
              Allowed spectral type range is M6 to T8.5, and allowed filters are MKO J, H, K
            - *dupuy2013*: BC/SpT relation from `Dupuy & Kraus (2013) <http://adsabs.harvard.edu/abs/2013Sci...341.1492D>`_
              Allowed spectral type range is T8 to Y0.5, and allowed filters are MKO Y, J, H
            - *filippazzo2015* (default): BC/SpT relation from `Filippazzo et al. (2015) <http://adsabs.harvard.edu/abs/2013Sci...341.1492D>`_
              Allowed spectral type range is M6 to T8/9, and allowed filters are 2MASS J, Ks
            - *filippazzo2015-young*: BC/SpT relation for young sources from `Filippazzo et al. (2015) <http://adsabs.harvard.edu/abs/2013Sci...341.1492D>`_
              Allowed spectral type range is M7 to T8, and allowed filters are 2MASS J, Ks

        :param unc: uncertainty of ``spt`` (default = 0)
        :param reverse: apply reverse approach: given BC, infer spectral type
        :param nsamples: number of Monte Carlo samples for error computation (default = 100)

    :Output: 
    
        2 element tuple providing the absolute magnitude and its uncertainty

    :Example:
        >>> import splat
        >>> print splat.typeToBC('L3', '2MASS J')
    """

#Keywords alternatives
    ref = kwargs.get('reference', ref)
    ref = kwargs.get('set', ref)
    unc = kwargs.get('uncertainty', unc)
    unc = kwargs.get('error', unc)


# check that you can use the proscribed relation and filter
    filtcheck = checkFilterName(filt,verbose=verbose)
    if filtcheck == False: return numpy.nan,numpy.nan
    else: filt=filtcheck

    refcheck = checkBC(ref,filt=filt,verbose=verbose)
    if refcheck == False: return numpy.nan,numpy.nan
    else: ref=refcheck

    refstring = 'BC/SpT relation for filter {} from {}'.format(filt,shortRef(BC_SETS[ref]['bibcode']))
    if verbose: print('\nUsing {}'.format(refstring))

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
        if BC_SETS[ref]['method'] == 'polynomial':
            rng = BC_SETS[ref]['filters'][filt]['range']
            if (rng[0] <= sptn <= rng[1]):
                bc = numpy.polyval(BC_SETS[ref]['filters'][filt]['coeff'], sptn-BC_SETS[ref]['sptoffset'])
                bc_error = BC_SETS[ref]['filters'][filt]['fitunc']
                if unc > 0.:
                    vals = numpy.polyval(BC_SETS[ref]['filters'][filt]['coeff'], numpy.random.normal(sptn - BC_SETS[ref]['sptoffset'], unc, nsamples))
                    bc_error = (numpy.nanstd(vals)**2+BC_SETS[ref]['filters'][filt]['fitunc']**2)**0.5
                return bc, bc_error
            else:
                if verbose: sys.stderr.write('\nSpectral type {} is out of range for {}'.format(typeToNum(sptn),refstring))
                return numpy.nan, numpy.nan

# interpolation method
        elif BC_SETS[ref]['method'] == 'interpolate':
            rng = [numpy.nanmin(BC_SETS[ref]['filters'][filt]['spt']),numpy.nanmax(BC_SETS[ref]['filters'][filt]['spt'])]
            if (rng[0] <= sptn <= rng[1]):
                f = interp1d(BC_SETS[ref]['filters'][filt]['spt'],BC_SETS[ref]['filters'][filt]['bc'])
                fe = interp1d(BC_SETS[ref]['filters'][filt]['spt'],BC_SETS[ref]['filters'][filt]['rms'])
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
            raise ValueError('Unknown method {} for {}'.format(BC_SETS[ref]['method'],refstring))

# reverse approach: BC -> SpT
    else:
        if not isinstance(spt,float): raise ValueError('Running this in reverse you need to provide a BC value instead of {}'.format(spt))
        if BC_SETS[ref]['method'] == 'polynomial':
            rng = BC_SETS[ref]['filters'][filt]['range']
            x = numpy.linspace(rng[0],rng[1],nsamples)
            y = numpy.polyval(BC_SETS[ref]['filters'][filt]['coeff'], x-BC_SETS[ref]['sptoffset'])
            f = interp1d(y,x)
            try:
                bc = float(f(spt))
            except:
                if verbose: print('\nBC value {} is outside the range expected for {}'.format(spt,refstring))
                return numpy.nan, numpy.nan
            vals = []
            for i in range(nsamples):
                ye = y+numpy.random.normal(0,BC_SETS[ref]['filters'][filt]['fitunc'])
                f = interp1d(ye,x)
                try:
                    vals.append(f(numpy.random.normal(spt,unc)))
                except:
                    pass
            bc_error = numpy.nanstd(vals)
            return bc, bc_error
        elif BC_SETS[ref]['method'] == 'interpolate':
            f = interp1d(BC_SETS[ref]['filters'][filt]['bc'],BC_SETS[ref]['filters'][filt]['spt'])
            try:
                bc = f(spt)
            except:
                if verbose: print('\nBC value {} is outside the range expected for {}'.format(spt,refstring))
                return numpy.nan, numpy.nan
            vals = []
            for i in range(nsamples):
                y = numpy.random.normal(BC_SETS[ref]['filters'][filt]['bc'],BC_SETS[ref]['filters'][filt]['rms'])
                f = interp1d(y,BC_SETS[ref]['filters'][filt]['spt'])
                try:
                    vals.append(f(numpy.random.normal(spt,unc)))
                except:
                    pass
            bc_error = numpy.nanstd(vals)
            return bc, bc_error
        else:
            raise ValueError('Unknown method {} for {}'.format(BC_SETS[ref]['method'],refstring))





