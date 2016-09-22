from __future__ import print_function, division

"""
.. note::
         Using a suite of evolutionary models, this code translates 
         between the following brown dwarf parameters: mass, age, 
         temperature, radius, surface gravity, and luminosity. We allow 
         the user to choose a set of evolutionary model 
         (Baraffe, Burrows, or Saumon) and two parameters, then output
         the rest of the interpolated parameters. 
"""

# Standard library imports.
import copy
import requests
import sys

# Related third party imports.
from astropy import units as u
from astropy.cosmology import Planck15, z_at_value
from astropy.io import ascii
import pandas
from math import isnan
import matplotlib.pyplot as plt
import numpy
from scipy.interpolate import interp1d 
import scipy.integrate as integrate 
from splat import SPLAT_PATH, SPLAT_URL
EVOLUTIONARY_MODEL_FOLDER = '/reference/EvolutionaryModels/'
EMODELS = ['baraffe','burrows','saumon']
EPARAMETERS = ['mass','age','temperature','gravity','luminosity','radius']
EPARAMETER_UNITS = {'mass': u.solMass, 'age': u.Gyr, 'temperature': u.K, 'gravity': u.dex(u.cm / u.s**2),\
    'luminosity': u.dex(u.solLum), 'radius': u.solRad}

# change the command prompt
sys.ps1 = 'splat evolve> '


###############################################################################
###############################################################################
def loadEvolModel(*model,**kwargs):
    '''
    :Purpose: Reads in the evolutionary model parameters for the models listed below, which are used to interpolate parameters in `modelParameters()`_. 

    .. _`modelParameters()` : api.html#splat_evolve.modelParameters

    Available models are:

        - **burrows** : Models from `Burrows et al. (2001) <http://adsabs.harvard.edu/abs/2001RvMP...73..719B>`_ for 1 Myr < age < 10 Gyr, 0.005 Msol < mass < 0.2 Msol, and solar metallicity
        - **baraffe** : Models from `Baraffe et al. (2003) <http://adsabs.harvard.edu/abs/2003A&A...402..701B>`_ for 1 Myr < age < 10 Gyr, 0.005 Msol < mass < 0.1 Msol, and solar metallicity (COND dust prescription)
        - **saumon** : Models from `Saumon et al. (2003) <http://adsabs.harvard.edu/abs/2008ApJ...689.1327S>`_ for 3 Myr < age < 10 Gyr, 0.002 Msol < mass < 0.085 Msol, although mass and age ranges vary as the maximum temperature for the models is 2500 K. For these models there are additional options:
        
            - **metallicity** = `solar`, `+0.3`, or `-0.3`
            - **cloud** =  `cloud-free`, `hybrid`, `f2` (sub- and super-solar metallicities are only cloud-free)

    Parameter units (in astropy convention) are: 

        - `masses`: Solar masses
        - `ages`: Gyr
        - `temperature`: K
        - `gravity`: log10 of cm/s/s
        - `luminosity`: log10 of Solar luminosities
        - `radius`: Solar radii

    Models are contained in SPLAT's reference/EvolutionaryModels folder.

    Required Inputs:

    :param: model: string of the name of the evolutionary model set to be used; can be `baraffe` (default), `burrows`, or `saumon`

    Optional Inputs:

    :param: metallicity: for Saumon models, this is the metallicity assumed, and can be a string or integer.  Allowed values are 0 (or `solar` = default), -0.3 (or `subsolar`) or 0.3 (or `supersolar`)
    :param: cloud: for Saumon models, this is the desired cloud prescription, and is a string:

        - no clouds: `cloud` = `nocloud`, `cloud-free` or `nc` (default)
        - hybrid cloud prescription: `cloud` = `hybrid`
        - f2 cloud prescription: `cloud` = `f2`

    Output: 

    Dictionary containing keywords mass, age, temperature, luminosity, gravity, and radius, each linked to the evolutionary parameters retrieved. 

    :Example:
    >>> import splat
    >>> p = splat.loadEvolModel('saumon',metallicity=-0.3,cloud='nc')
    You are using saumon's models.
    >>> for k in list(p.keys()): print('{}: {}'.format(k, p[k][12]))
    age: 0.15
    mass: [ 0.002  0.003  0.004  0.005  0.006  0.007  0.008  0.009  0.01   0.011
      0.012  0.013  0.014  0.015  0.016  0.017  0.018  0.019  0.02   0.022
      0.024  0.026  0.028  0.03   0.033  0.035  0.038  0.04   0.043  0.045
      0.048  0.05   0.053]
    temperature: [  353.   418.   471.   523.   585.   642.   695.   748.   806.   893.
      1146.  1228.  1114.  1113.  1148.  1183.  1227.  1270.  1316.  1402.
      1489.  1572.  1654.  1739.  1853.  1930.  2030.  2096.  2187.  2240.
      2316.  2362.  2426.] 
    gravity: [ 3.576  3.746  3.871  3.972  4.056  4.128  4.191  4.246  4.296  4.335
      4.337  4.368  4.437  4.479  4.512  4.543  4.571  4.597  4.621  4.665
      4.704  4.74   4.772  4.8    4.839  4.861  4.892  4.909  4.931  4.947
      4.966  4.978  4.996]
    luminosity: [-6.691 -6.393 -6.185 -6.006 -5.815 -5.658 -5.527 -5.404 -5.277 -5.098
     -4.628 -4.505 -4.709 -4.724 -4.675 -4.627 -4.568 -4.51  -4.45  -4.342
     -4.24  -4.146 -4.058 -3.969 -3.856 -3.781 -3.69  -3.628 -3.546 -3.5   -3.432
     -3.393 -3.34 ]
    radius: [ 0.1206  0.1214  0.1214  0.1209  0.1202  0.1195  0.1189  0.1182  0.1178
      0.1181  0.123   0.1235  0.1184  0.1167  0.1161  0.1154  0.1151  0.1148
      0.1146  0.1142  0.1139  0.1139  0.1138  0.1141  0.1144  0.115   0.1155
      0.1163  0.1174  0.118   0.1193  0.12    0.121 ]
    '''

# check model
    try: model = model[0].lower()
    except TypeError: raise TypeError("Model must be a string.")
    except IndexError: model = 'baraffe'
    finally: print("You are using " + model + "'s models.")
    assert model in EMODELS, "\nModel {} not in allowed model sets; please use: {}\n".format(model,' '.join(EMODELS))

##################### BARAFFE OR BURROWS MODEL #########################
    if model == 'baraffe' or model == 'burrows':
        ages = ['0.001', '0.005', '0.010', '0.050', '0.100',
                '0.120', '0.500', '1.000', '5.000', '10.000']

        if model == 'baraffe': 
            prefix = 'Baraffe/cond_'
        else: 
            prefix = 'Burrows/b97_'

########################### SAUMON MODEL ##############################
    else:

# set metallicity
        metallicity = kwargs.get('z',False)
        metallicity = kwargs.get('metallicity',metallicity)
        if metallicity == False:
            metallicity = 'solar'
        if isinstance(metallicity,int) or isinstance(metallicity,float):
            metallicity = '{:.1f}'.format(metallicity)
        if metallicity.lower() == 'solar' or metallicity == '0.0':
            Z = 'z0'
        elif metallicity == '0.3' or metallicity == '+0.3' or metallicity == 'supersolar':
            Z = 'z+0.3'
        elif metallicity == '-0.3' or metallicity == 'subsolar':
            Z = 'z-0.3'
        else:
            raise ValueError('\nMetallicity for Saumon model must be 0.0 (solar), +0.3 or -0.3, not {}\n'.format(metallicity))
# set cloud treatment
        cloud = kwargs.get('cloud',False)
        cloud = kwargs.get('clouds',cloud)
        cloud = kwargs.get('cld',cloud)
        if cloud == False:
            cloud = 'nc'
        if isinstance(cloud,int) or isinstance(cloud,float):
            cloud = 'f{:1d}'.format(int(cloud))
        if cloud.lower() == 'hybrid':
            C = 'hybrid'
        elif cloud.lower() == 'f2':
            C = 'f2'
        elif cloud.lower() == 'nc' or cloud.lower() == 'nocloud' or cloud.lower() == 'noclouds' or cloud.lower() == 'cloud-free':
            C = 'nc'
        if metallicity=='-0.3' or metallicity=='0.3':
            C = 'nc'
        else:
            raise ValueError('\nCould not recognize cloud choice for Saumon model: must be cloud-free, hybrid or f2, not {}\n'.format(cloud))

        ages = ['0.003','0.004','0.006','0.008','0.010','0.015','0.020',
                '0.030','0.040','0.060','0.080','0.100','0.150','0.200',
                '0.300','0.400','0.600','0.800','1.000','1.500','2.000',
                '3.000','4.000','6.000','8.000','10.000']

        prefix = 'Saumon/sau08_{:s}_{:s}_'.format(Z,C)

#######################################################################

# read in parameters
    mparam = {}
    for ep in EPARAMETERS:
        mparam[ep] = []

    for i,age in enumerate(ages):
        mfile = prefix+'{:05d}.txt'.format(int(float(age)*1000.))
        try:
            dp=pandas.read_csv(SPLAT_PATH+EVOLUTIONARY_MODEL_FOLDER+mfile,comment='#',sep=',',header=0)
            for ep in EPARAMETERS:
                mparam[ep].append(dp[ep].values)

# this is done in case models are not local - NOTE: currently just throwing an error
        except:
            raise ValueError('Could not find model file {} locally; aborting'.format(mfile))
#            try:
#                print('Could not read in model file {} locally; trying online'.format(mfile))
#                data =ascii.read(requests.get(SPLAT_URL+EVOLUTIONARY_MODEL_FOLDER+mfile).content,comment='#',delimiter='\t')
#            except: 
#                raise ValueError('Could not find model file {} locally or online; aborting'.format(mfile))

    mparam['age'] = [float(i) for i in ages]

    return mparam



def modelParametersSingle(*args, **kwargs):
    '''
    :Purpose: Driver function for modelParameters_, performs actual interpolation of evolutionary models. See SPLAT API for `modelParameters()`_ for details.

    .. _`modelParameters()` : api.html#splat_evolve.modelParameters

    '''

    keywords = list(kwargs.keys())

# check that model is passed correctly
    try: model = args[0]
    except IndexError: 
        model = loadEvolModel('baraffe')
        print('\nWarning: model error; using Baraffe models by default\n')

# retool models to allow for logarithmic interpolation
    lmodel = copy.deepcopy(model)
# strip off units
    lmodel['age'] = [numpy.log10(m) for m in lmodel['age']]
    for i in range(len(lmodel['age'])):
        lmodel['mass'][i] = [numpy.log10(m) for m in lmodel['mass'][i]]
        lmodel['temperature'][i] = [numpy.log10(m) for m in lmodel['temperature'][i]]
        lmodel['radius'][i] = [numpy.log10(m) for m in lmodel['radius'][i]]

# prep output parameters
    params = {}
    for e in EPARAMETERS:
        params[e] = 0.

    for e in EPARAMETERS:
        if e in keywords:
            try: f = float(kwargs[e])
            except: raise ValueError('\nInput paramter {} must be a single number, not {}\n'.format(e,kwargs[e]))
            finally: params[e] = f

    input_type = 'mass_age'
    Ag, Ma, Te, Le, Ge, Re, P = [],[],[],[],[],[],[]


############### UNKNOWN MASS AND AGE - INTERPOLATE AGE FROM OTHER PARAMETERS #################
# for each age, interpolate mass as a function of first parameter and then second parameter as a function of mass
# and obtain second parameter as a function of age; then interpolate the model ages as a function of 
# the second parameter and evaluate for known parameter to get age
###############################################################################
    if (params['mass'] == 0.) and (params['age'] == 0.):

        input_type = 'two_params'
        if params['temperature'] != 0.:
            P.append(['temperature', numpy.log10(params['temperature'])])
        if params['gravity'] != 0.:
            P.append(['gravity', params['gravity']])
        if params['radius'] != 0.:
            P.append(['radius', numpy.log10(params['radius'])])
        if params['luminosity'] != 0.:
            P.append(['luminosity', params['luminosity']])

        for i,age in enumerate(lmodel['age']):
            if min(lmodel[P[0][0]][i]) <= P[0][1] <= max(lmodel[P[0][0]][i]) \
                and min(lmodel[P[1][0]][i]) <= P[1][1] <= max(lmodel[P[1][0]][i]):
                Ag.append(age)
                f = interp1d(lmodel[P[0][0]][i], lmodel['mass'][i])
                Ma = f(P[0][1])
                f = interp1d(lmodel['mass'][i], lmodel[P[1][0]][i])
                Ge.append(f(Ma))

        try: 
            f = interp1d(Ge, Ag)
            params['age'] = 10.**f(P[1][1])
        except: params['age'] = float('nan')

        Ge, Ag, Ma = [], [], []

################ UNKNOWN AGE BUT KNOWN MASS AND ONE OTHER PARAMETER ###########
# interpolate second parameter as a function of mass for each of the age models and evaluate for known mass
# interpolate the model ages as a fucntion of these parameters and evaluate for known parameter
###############################################################################
    if params['age'] == 0. and params['mass'] != 0. and \
        not isnan(params['mass']):

        if input_type != 'two_params': 
            input_type = 'one_param'
            if params['temperature'] != 0.:
                P.append(['temperature', numpy.log10(params['temperature'])])
            elif params['gravity'] != 0.:
                P.append(['gravity', params['gravity']])
            elif params['radius'] != 0.:
                P.append(['radius', numpy.log10(params['radius'])])
            elif params['luminosity'] != 0.:
                P.append(['luminosity', numpy.log10(params['luminosity'])])
            else:
                for k in list(params.keys()):
                    print('{}: {}'.format(k,params[k]))
                print(P)
                raise ValueError('\nProblem with one_param interpolation\n')

        for i,age in enumerate(lmodel['age']):
            if min(lmodel['mass'][i]) <= numpy.log10(params['mass']) <= max(lmodel['mass'][i]):
                Ag.append(age)
                f = interp1d(lmodel['mass'][i], lmodel[P[0][0]][i])
                Ge.append(f(numpy.log10(params['mass'])))

        try: 
            f = interp1d(Ge, Ag)
            params['age'] = 10.**f(P[0][1])
        except: 
            print('\nFailed in age + parameter determination\n')
            params['age'] = float('nan')

        Ge, Ag = [], []


################ KNOWN AGE BUT UNKNOWN MASS AND ONE OTHER PARAMETER ###########
# generate mass as function of second parameter interpolated between two closest age models
# evaluate mass(parameter) (resulting in both mass and age as knowns)
###############################################################################

    if params['age'] != 0. and params['mass'] == 0. and \
        not isnan(params['age']):

        if input_type != 'two_params' and input_type != 'one_param': 
            input_type = 'one_param'
            if params['temperature'] != 0.:
                P.append(['temperature', numpy.log10(params['temperature'])])
            elif params['gravity'] != 0.:
                P.append(['gravity', params['gravity']])
            elif params['radius'] != 0.:
                P.append(['radius', numpy.log10(params['radius'])])
            elif params['luminosity'] != 0.:
                P.append(['luminosity', numpy.log10(params['luminosity'])])
            else:
                for k in list(params.keys()):
                    print('{}: {}'.format(k,params[k]))
                print(P)
                raise ValueError('\nProblem with one_param interpolation\n')

        if numpy.log10(params['age']) < numpy.min(lmodel['age']) or \
            numpy.log10(params['age']) > numpy.max(lmodel['age']):
                print('\nAge of {} is outside range of models, {} to {}\n'.format(params['age'],10.**numpy.min(lmodel['age']),10**numpy.max(lmodel['age'])))
                params['mass'] = numpy.nan

        else:
            adiff = [numpy.log10(params['age'])-a for a in lmodel['age']]
            ai = numpy.argmin(numpy.abs(adiff))
            if adiff[ai] < 0:
                ai-=1
            for i,m in enumerate(lmodel['mass'][ai]):
                if m in lmodel['mass'][ai+1]:
                    Ma.append(m)
                    aj = numpy.argmin(numpy.abs([a-m for a in lmodel['mass'][ai+1]]))
                    vals = [lmodel[P[0][0]][ai][i],lmodel[P[0][0]][ai+1][aj]]
                    f = interp1d(lmodel['age'][ai:ai+2],vals)
                    Ge.append(f(numpy.log10(params['age'])))
            try:
                f = interp1d(Ge, Ma)
                params['mass'] = 10.**f(P[0][1])
            except:
                print('\nFailed in mass + parameter determination\n')
                params['mass'] = numpy.nan

        Ma, Ge = [],[]


###################### KNOWN MASS AND AGE #####################################
# generate parameters as a function of mass interpolated between two closest age models
# evaluate parameters(mass)
###############################################################################
    if params['mass'] != 0. and params['age'] != 0. and \
        not isnan(params['age']) and not isnan(params['mass']):

        for i,age in enumerate(lmodel['age']):
            if min(lmodel['mass'][i]) <= numpy.log10(params['mass']) \
                                              <= max(lmodel['mass'][i]):
                Ag.append(age)
                f =interp1d(lmodel['mass'][i],lmodel['temperature'][i])
                Te.append(f(numpy.log10(params['mass'])))
                f = interp1d(lmodel['mass'][i],lmodel['luminosity'][i])
                Le.append(f(numpy.log10(params['mass'])))
                f = interp1d(lmodel['mass'][i],lmodel['gravity'][i])
                Ge.append(f(numpy.log10(params['mass'])))
                f = interp1d(lmodel['mass'][i],lmodel['radius'][i])
                Re.append(f(numpy.log10(params['mass'])))
  
        if params['temperature'] == 0.:
            try:
                f = interp1d(Ag, Te)
                params['temperature'] = 10.**f(numpy.log10(params['age']))
            except: 
                params['temperature'] = numpy.nan
        if params['luminosity'] == 0.:
            try: 
                f = interp1d(Ag, Le)
                params['luminosity'] = f(numpy.log10(params['age'])).item(0)
            except: 
                params['luminosity'] = numpy.nan
        if params['gravity'] == 0.:
            try: 
                f = interp1d(Ag, Ge) 
                params['gravity'] = f(numpy.log10(params['age'])).item(0)        
            except: 
                params['gravity'] = numpy.nan
        if params['radius'] == 0.:
            try: 
                f = interp1d(Ag, Re)
                params['radius'] = 10.**f(numpy.log10(params['age']))
            except: 
                params['radius'] = numpy.nan
  
        return params

# something failed	  
    else:
        for e in EPARAMETERS:
            params[e] = numpy.nan
        print('\nParameter set is not covered by models\n')
        return params
      


###############################################################################

def modelParameters(*model,**kwargs):
    '''
    :Purpose: Retrieves the evolutionary model parameters given two of the following parameters: mass, age, temperature, luminosity, gravity, or radius. The inputs can be individual values or arrays.  Using the input parameters, the associated evolutionary model parameters are computed through log-linear interpolation of the original model grid. Parameters that fall outside the grid return nan.

    Required Inputs:

    :param: model: Either a string of the name of the evolutionary model set, which can be one of `baraffe` (default), `burrows`, or `saumon`; or a dictionary output from `loadEvolModel()`_ containing model parameters. 
    
    and two (2) of the following:

    :param: mass: input value of list of values for mass (can also be `masses` or `m`)
    :param: age: input value of list of values for age (can also be `ages`, `time` or `a`)
    :param: temperature: input value of list of values for temperature (can also be `temperatures`, `teff`, `temp` or `t`)
    :param: gravity: input value of list of values for gravity (can also be `gravities`, `grav`, `logg` or `g`)
    :param: luminosity: input value of list of values for luminosity (can also be `luminosities`, `lum`, `lbol` or `l`)
    :param: radius: input value of list of values for radius (can also be `radii`, `rad` and `r`)

    .. _`loadEvolModel()` : api.html#splat_evolve.loadEvolModel

    Optional Inputs:

    :param: Parameters for `loadEvolModel()`_ may also be used.

    Output: 

    Dictionary containing keywords mass, age, temperature, luminosity, gravity, and radius, each linked to the evolutionary parameters retrieved. 


    :Example:
    >>> import splat, numpy
    >>> masses = numpy.random.uniform(0.01,0.1,20)
    >>> ages = numpy.random.uniform(0.01,10,20)
    >>> p = splat.modelParameters('baraffe',mass=masses,age=ages)
    You are using baraffe's models.
    >>> print(p.temperature)
    [ 2502.90132332  2818.85920306  1002.64227134  1330.37273021  1192.86976417
      500.45609068  2604.99966013  1017.03307609  1774.18267474  1675.12181635
      2682.9697321   2512.45223777   346.41152614  2066.19972036   843.28528456
      2264.93051445  2767.85660557   348.84214986   922.87030167  2669.27152307] K    
    '''

# read in model
    try: model = model[0]
    except IndexError: 
        if kwargs.get('model',False)==False:
            model = 'baraffe'
        else: model=kwargs.get('model')
    if type(model) is not dict: model = loadEvolModel(model,**kwargs)

    keywords = list(kwargs.keys())

# do some key word replacement
    mkwargs = {}
    for e in EPARAMETERS:
        if e in keywords:
            mkwargs[e] = kwargs[e]
    if 'temperature' not in keywords:
        if 't' in keywords:
            mkwargs['temperature'] = kwargs['t']
        if 'teff' in keywords:
            mkwargs['temperature'] = kwargs['teff']
        if 'temp' in keywords:
            mkwargs['temperature'] = kwargs['temp']
    if 'gravity' not in keywords:
        if 'g' in keywords:
            mkwargs['gravity'] = kwargs['g']
        if 'logg' in keywords:
            mkwargs['gravity'] = kwargs['logg']
        if 'grav' in keywords:
            mkwargs['gravity'] = kwargs['grav']
    if 'mass' not in keywords:
        if 'm' in keywords:
            mkwargs['mass'] = kwargs['m']
    if 'age' not in keywords:
        if 'time' in keywords:
            mkwargs['age'] = kwargs['time']
        if 'a' in keywords:
            mkwargs['age'] = kwargs['a']
    if 'radius' not in keywords:
        if 'r' in keywords:
            mkwargs['radius'] = kwargs['r']
        if 'rad' in keywords:
            mkwargs['radius'] = kwargs['rad']
    if 'luminosity' not in keywords:
        if 'l' in keywords:
            mkwargs['luminosity'] = kwargs['l']
        if 'lum' in keywords:
            mkwargs['luminosity'] = kwargs['lum']
        if 'lbol' in keywords:
            mkwargs['luminosity'] = kwargs['lbol']


# determine length of input arrays
    inparams = {}
    outparams = {}
    pkeys = list(mkwargs.keys())
    for p in EPARAMETERS:
        outparams[p] = []
        if p in pkeys:
            if isinstance(mkwargs[p],float) or isinstance(mkwargs[p],int):
                mkwargs[p] = [mkwargs[p]]
            numberValues = len(mkwargs[p])

# now loop through each parameter set to determine remaining parameters
    for i in range(numberValues):
        for p in pkeys:
            inparams[p] = mkwargs[p][i]
        par = modelParametersSingle(model,**inparams)
        for p in EPARAMETERS:
            outparams[p].append(par[p])


# remove lists if only one parameter set is being calculated
    if len(outparams['temperature']) == 1:
        for e in EPARAMETERS:
            outparams[e] = outparams[e][0]

# add units
    for e in EPARAMETERS:
        outparams[e] *= EPARAMETER_UNITS[e]

    return outparams




def plotModelParameters(parameters,xparam,yparam,**kwargs):
    '''
    :Purpose: Plots pairs of physical star parameters and optionally compares to evolutionary model tracks. 

    Required Inputs:

    :param: parameters: dictionary or nested set of two arrays containing parameters to be plotted. For dictionary, keywords should include the `xparameter` and `yparameter` strings to be plotted. Values associated with keywords can be single numbers or arrays
    :param: xparam: string corresponding to the key in the `parameters` dictionary to be plot as the x (independent) variable. 
    :param: yparam: string corresponding to the key in the `parameters` dictionary to be plot as the y (dependent) variable. 


    Optional Inputs:

    .. _`loadEvolModel()` : api.html#splat_evolve.loadEvolModel

    :param: showmodel: set to True to overplot evolutionary model tracks from `model` (default = True)
    :param: model: either a string of the name of the evolutionary model set, one of  `baraffe` (default), `burrows`, or `saumon`; or a dictionary output from `loadEvolModel()`_ containing model parameters. 
    :param: tracks: string indicating what model tracks to show; can either be `mass` (default) or `age`


    :param: file: name of file to output plot (`output` can also be used)
    :param: show: set to True to show the plot onscreen (default = True)
    :param: figsize: a two-element array defining the figure size (default = [8,6])

    :param: color: color of data symbols (default = 'blue')
    :param: marker: matplotlib marker type for data symbols (default = 'o')
    :param: xlabel: string overriding the x-axis label (default = parameter name and unit)
    :param: ylabel: string overriding the y-axis label (default = parameter name and unit)
    :param: title: string specifying plot title (no title by default)
    :param: tight: set to True to tighten plot to focus on the data points (default = True)

    Output: 

    A matplotlib plot object. Optionally, can also show plot on screen or output plot to a file.

    :Example:
    >>> import splat, numpy
    >>> age_samp = 10.**numpy.random.normal(numpy.log10(1.),0.3,50)
    >>> mass_samp = numpy.random.uniform(0.001,0.1,50)
    >>> p = splat.modelParameters('baraffe',age=age_samp,mass=mass_samp)
    >>> splat.plotModelParameters(p,'age','temperature',showmodels=True,model='baraffe',show=True)
    [plot of temperature vs age for 50 data points with baraffe models overplotted]    
    '''
# check inputs
    if type(parameters) is not dict:
        if len(parameters) != 2:
            raise ValueError('\nInput parameters should be a dictionary or two-element list\n')
        else:
            param = {xparam: parameters[0], yparam: parameters[1]}
    else:
        param = copy.deepcopy(parameters)

    keys = list(param.keys())
    if xparam not in keys:
        raise ValueError('\nCould not find parameter {} in input dictionary\n'.format(xparam))
    if yparam not in keys:
        raise ValueError('\nCould not find parameter {} in input dictionary\n'.format(yparam))

    if isinstance(param[xparam],list) == False:
        param[xparam] = [param[xparam]]
    if isinstance(param[yparam],list) == False:
        param[yparam] = [param[yparam]]


# sort flags
    if xparam=='age' or xparam=='time' or xparam=='a':
        xmparam = 'age'
        xlogflag = True
    elif xparam=='mass' or xparam=='m':
        xmparam = 'mass'
        xlogflag = True
    elif xparam=='temperature' or xparam=='teff' or xparam=='t':
        xmparam = 'temperature'
        xlogflag = True
    elif xparam=='radius' or xparam=='r':
        xmparam = 'radius'
        xlogflag = True
    elif xparam=='gravity' or xparam=='logg' or xparam=='g':
        xmparam = 'gravity'
        xlogflag = False
    elif xparam=='luminosity' or xparam=='lbol' or xparam=='l':
        xmparam = 'luminosity'
        xlogflag = False
    else:
        raise ValueError('\nx-axis parameter {} is not one that can be plotted'.format(xparam))

    if yparam=='age' or yparam=='time' or yparam=='a':
        ymparam = 'age'
        ylogflag = True
    elif yparam=='mass' or yparam=='m':
        ymparam = 'mass'
        ylogflag = True
    elif yparam=='temperature' or yparam=='teff' or yparam=='t':
        ymparam = 'temperature'
        ylogflag = True
    elif yparam=='radius' or yparam=='r':
        ymparam = 'radius'
        ylogflag = True
    elif yparam=='gravity' or yparam=='logg' or yparam=='g':
        ymparam = 'gravity'
        ylogflag = False
    elif yparam=='luminosity' or yparam=='lbol' or yparam=='l':
        ymparam = 'luminosity'
        ylogflag = False
    else:
        raise ValueError('\ny-axis parameter {} is not one that can be plotted'.format(yparam))

# plot parameters
    plt.close('all')
    plt.figure(figsize=kwargs.get('figsize',[8,6]))
    if xlogflag == True and ylogflag == True:
        plt.loglog(param[xparam],param[yparam],color=kwargs.get('color','blue'),marker=kwargs.get('marker','o'))
    elif xlogflag == False and ylogflag == True:
        plt.semilogy(param[xparam],param[yparam],color=kwargs.get('color','blue'),marker=kwargs.get('marker','o'))
    elif xlogflag == True and ylogflag == False:
        plt.semilogx(param[xparam],param[yparam],color=kwargs.get('color','blue'),marker=kwargs.get('marker','o'))
    else:
        plt.plot(param[xparam],param[yparam],color=kwargs.get('color','blue'),marker=kwargs.get('marker','o'))

# read in models to display   
    if kwargs.get('showmodel',True) != False or kwargs.get('showmodels',True) != False:
        if kwargs.get('model',False) == False:
            model = 'baraffe'
        else:
            model = kwargs.get('model')
        try:
            if type(model) is not dict: model = loadEvolModel(model,**kwargs)
        except:
            print('\nProblem in reading in original models\n')
            kwargs['showmodel'] = False

    if kwargs.get('showmodel',True) != False or kwargs.get('showmodels',True) != False:
        tvals,xvals,yvals = [],[],[]

# models tracks trace mass (by default)
        if kwargs.get('tracks','mass') == 'mass':
            masses = []
            for i in model['mass']:
                masses.extend(i)
            masses.sort()
            tvals = numpy.unique(masses)
            for j,m in enumerate(tvals):
                xx,yy = [],[]
                for i,x in enumerate(model['age']):
                    if m in model['mass'][i]:
                        if xmparam != 'age':
                            xx.append(numpy.array(model[xmparam][i])[numpy.where(model['mass'][i]==m)].item(0))
                        else:
                            xx.append(x)
                        if ymparam != 'age':
                            yy.append(numpy.array(model[ymparam][i])[numpy.where(model['mass'][i]==m)].item(0))
                        else:
                            yy.append(x)
                    else:
                        xx.append(numpy.nan)
                        yy.append(numpy.nan)                        
                xvals.append(xx)
                yvals.append(yy)

# models tracks trace isochrones
        else:
            tvals = model['age']
# fix to account for unequal lengths of model values
            maxlen = numpy.max([len(a) for a in models['mass']])
            for i,x in enumerate(tvals):
                t = numpy.zeros(maxlen)
                t.fill(numpy.nan)
                if xparam != 'age':
                    t[0:len(model[xparam][i])] = model[xmparam][i]
                else:
                    t.fill(x)
                xvals.append(t.tolist())
                s = numpy.zeros(maxlen)
                s.fill(numpy.nan)
                if yparam != 'age':
                    s[0:len(model[yparam][i])] = model[ymparam][i]
                else:
                    s.fill(x)
                yvals.append(s.tolist())

# plot them
        for i,x in enumerate(xvals):
            if xlogflag == True and ylogflag == True:
                plt.loglog(xvals[i],yvals[i],color='grey')
            elif xlogflag == False and ylogflag == True:
                plt.semilogy(xvals[i],yvals[i],color='grey')
            elif xlogflag == True and ylogflag == False:
                plt.semilogx(xvals[i],yvals[i],color='grey')
            else:
                plt.plot(xvals[i],yvals[i],color='grey')

# add labels
    plt.xlabel(kwargs.get('xlabel','{} ({})'.format(xmparam,EPARAMETER_UNITS[xmparam])))
    plt.ylabel(kwargs.get('ylabel','{} ({})'.format(ymparam,EPARAMETER_UNITS[ymparam])))
    if kwargs.get('title',False) != False:
        plt.title(kwargs.get('title'))
# tighten plot
    if kwargs.get('tight',True) == True:
        xrng = [numpy.nanmin(param[xparam]),numpy.nanmax(param[xparam])]
        if xlogflag ==True:
            xsep = xrng[1]/xrng[0]
            if xsep != 1.:
                plt.xlim([xrng[0]/(xsep**0.1),xrng[1]*(xsep**0.1)])
        else:
            xsep = xrng[1]-xrng[0]
            if xsep != 0.:
                plt.xlim([xrng[0]-0.05*xsep,xrng[1]+0.05*xsep])
        yrng = [numpy.nanmin(param[yparam]),numpy.nanmax(param[yparam])]
        if ylogflag ==True:
            ysep = yrng[1]/yrng[0]
            if ysep != 1.:
                plt.ylim([yrng[0]/(ysep**0.1),yrng[1]*(ysep**0.1)])
        else:
            ysep = yrng[1]-yrng[0]
            if ysep != 0.:
                plt.ylim([yrng[0]-0.05*ysep,yrng[1]+0.05*ysep])

# save the plot or display
    file = kwargs.get('file',False)
    file = kwargs.get('output',file)
    if file != False:
        plt.savefig(file)
    elif kwargs.get('show',True) == True:
        plt.show()
    else:
        pass

    return plt



def generateAges(num,**kwargs):
    '''
    :Purpose: Generates a distribution of ages based on the defined input distribution. 

    Required Inputs:

    :param: num: number of ages to generate

    Optional Inputs:

    :param: age_range: range of ages to draw from (default = [0.1,10.]); can also specify `minage` and `maxage`
    :param: distribution: either a string set to one of the following to define the type of age distribution (or reverse star formation rate) desired:
        * `uniform`: uniform distribution (default) 
        * `exponential`: exponential age distribution, P(t) ~ e\^(beta x t). You can specify the parameters `beta` or `tau` = 1/beta, or set ``distribution`` to `aumer` or `miller`
        * `double_exponential`: double exponential age distribution, P(t) ~ Ae\^(lambda x t) + e\^(beta x t). You can specify the parameters `beta`, `lambda` and `a` or set ``distribution`` to `aumer_double` (default parameters)
        * `cosmic` or `rujopakarn`: cosmic age distribution with P(t) ~ (1+z(t))\^alpha, where z is the redshift, which is converted to time using the Planck 2015 cosmology. You can specify the parameter `alpha` or set ``distribution`` to `rujopakarn` (default parameters)
        * `peaked`: age distribution that peaks at some early time, written in the form P(t) ~ (t-t0)/(t\^2+t1\^2)\^2. You can specify the parameters `t0` and `t1` or set ``distribution`` to `aumer_peaked` or `just_peaked`
        * `aumer` or `aumer_exponential`: exponential age distribution with parameters from Aumer & Binney (2009): beta = 0.117
        * `aumer_double`: double exponential age distribution with parameters from Aumer & Binney (2009): beta = 0.348, lambda = 2.0, a = 1.e-8
        * `aumer_peaked`: peaked age distribution with parameters from Aumer & Binney (2009): t0 = XXX, t1 = XXX
        * `just` or `just_exponential: exponential age distribution with parameters from Just & Jahriess (2010): beta = 0.125
        * `just_peaked_a`: peaked age distribution with parameters from Just & Jahriess (2010) Model A: t0 = 5.6, t1 = 8.2
        * `just_peaked` or `just_peaked_b`: peaked age distribution with parameters from Just & Jahriess (2010) Model B: t0 = 1.13, t1 = 7.8
        * `miller`: exponential age distribution with parameters from Miller & Scalo (1979): beta = max age / 2
        * `rujopakarn`: cosmic age distribution with parameters from Rujopakarn et al. (2010): beta = max age / 2
        * `input`: user specified age distribution or star formation history; ``input`` must be set to a 2 x N array specifying age and distribution
    :param: distribution can also be set to a 2 x N array specifying an age distribution or star formation history; the first vector should be the ages for the function and the second vector the distribution function
    :param: parameters: dictionary containing the parameters for the age distribution/star formation model being used; options include:
        * `alpha`: power law factor for cosmic age distribution
        * `beta`: power factor in exponential age distribution; positive beta implies a star formation rate that decreases with time
        * `lambda`: second power factor in double exponential age distribution; positive lambda implies a star formation rate that decreases with time
        * `a`: relative scale factor for second exponential in double exponential age distribution
        * `tau`: 1/beta scale factor in exponential age distribution
        * `t0` and `t1`: parameters for peaked age distribution
    :param: sfh: set to True if distribution is a star formation history rather than an age distribution (default = False)

    Output: 

    An array of ages drawn from the desired distribution in units of Gyr

    :Example:
    >>> import splat
    >>> import matplotlib.pyplot as plt
    >>> ages = splat.generateAges(10000,distribution='aumer',age_range=[0.3,8.0])
    >>> plt.hist(ages)
    [histogram of ages in range 0.3-8.0 Gyr]    
    '''

# initial parameters
    distribution = kwargs.get('distribution','uniform')
    mn = kwargs.get('minage',0.1)
    mx = kwargs.get('maxage',10.)
    sfh = kwargs.get('sfh',False)
    age_range = kwargs.get('age_range',[mn,mx])
# protective offset
    if age_range[0] == age_range[1]:
        age_range[1]+=0.0001
# set default parameters
    if kwargs.get('parameters',False) == False:
        parameters = {}
    else:
        parameters = kwargs['parameters']
    if 'beta' not in list(parameters.keys()):
        parameters['beta'] = 1.0
    if 'tau' not in list(parameters.keys()):
        parameters['tau'] = 1./parameters['beta']
    if 'alpha' not in list(parameters.keys()):
        parameters['alpha'] = 3.5
    if 'lambda' not in list(parameters.keys()):
        parameters['lambda'] = 2.0
    if 'a' not in list(parameters.keys()):
        parameters['a'] = 1.e-8
    if 't0' not in list(parameters.keys()):
        parameters['t0'] = 1.13
    if 't1' not in list(parameters.keys()):
        parameters['t1'] = 7.8

# 
# exponential
    if distribution.lower() == 'exponential' or distribution.lower() == 'aumer' or distribution.lower() == 'miller':
        print('using exponential distribution')
        if distribution.lower() == 'aumer':
            parameters['beta'] = 0.117
        if distribution.lower() == 'miller':
            parameters['beta'] = 0.5*numpy.max(age_range)
        if distribution.lower() == 'just':
            parameters['beta'] = 0.125

# use CDF sampling
        if parameters['beta'] != 0.:
            x = numpy.linspace(numpy.min(age_range),numpy.max(age_range),num=10000)
            y = numpy.exp(parameters['beta']*x)
            y -= numpy.min(y)
            y /= numpy.max(y)
            f = interp1d(y,x)
            ages = f(numpy.random.uniform(size=num))
        else:
            ages = numpy.random.uniform(numpy.min(age_range), numpy.max(age_range), size=num)

# double exponential
    elif distribution.lower() == 'double_exponential' or distribution.lower() == 'aumer_double':
        print('using double exponential distribution')
        if distribution.lower() == 'aumer_double':
            parameters['beta'] = 0.348
            parameters['lambda'] = 2.0
            parameters['a'] = 1.e-8

# use CDF sampling
        x = numpy.linspace(numpy.min(age_range),numpy.max(age_range),num=10000)
        y = parameters['a']*numpy.exp(parameters['lambda']*x) + numpy.exp(parameters['beta']*x)
        y -= numpy.min(y)
        y /= numpy.max(y)
        f = interp1d(y,x)
        ages = f(numpy.random.uniform(size=num))

# peaked distribution
    elif distribution.lower() == 'peaked' or distribution.lower() == 'just_peaked' or distribution.lower() == 'just_peaked_a' or distribution.lower() == 'just_peaked_b' or distribution.lower() == 'aumer_peaked':
        print('using peaked distribution')
# Aumer & Binney 2009
        if distribution.lower() == 'aumer_peaked':
            parameters['t0'] = 0.
            parameters['t1'] = 7.23
# Just & Jahriess 2010 Model A
        if distribution.lower() == 'just_peaked_a':
            parameters['t0'] = 5.6
            parameters['t1'] = 8.2
            sfh = True
# Just & Jahriess 2010 Model B (default)
        if distribution.lower() == 'just_peaked' or distribution.lower() == 'just_peaked_b':
            parameters['t0'] = 1.13
            parameters['t1'] = 7.8
            sfh = True

# generate CDF by integration and then do CDF sampling
# note that function is slightly different for the two forms
        x = numpy.linspace(numpy.min(age_range),numpy.max(age_range),num=10000)
        if 'just' in distribution:
            y = (x+parameters['t0'])/((x**2+parameters['t1']**2)**2)
#            print(2./3.*(t0**2+0.75*t1**2)**0.5 - 2./3.*t0)
        else:
            y = (14.-x+parameters['t0'])/(((14.-x)**2+parameters['t1']**2)**2)
#            print(14.-2./3.*(t0**2+0.75*t1**2)**0.5 - 2./3.*t0)
        yc = numpy.cumsum(y)
        yc -= numpy.min(yc)
        yc /= numpy.max(yc)
        f = interp1d(yc,x)
        ages = f(numpy.random.uniform(size=num))

# cosmic star formation rate
    elif distribution.lower() == 'cosmic' or distribution.lower() == 'rujopakarn': 
        print('using cosmic SFH distribution')
        if distribution.lower() == 'rujopakarn': 
            parameters['alpha'] = 3.5

        cosmo = Planck15    # in case we want to change later
        zrng = [z_at_value(cosmo.lookback_time,numpy.min(age_range)*u.Gyr),z_at_value(cosmo.lookback_time,numpy.max(age_range)*u.Gyr)]
# use CDF sampling
        x = numpy.linspace(numpy.min(zrng),numpy.max(zrng),num=10000)
        y = (x+1.)**parameters['alpha']
        y -= numpy.min(y)
        y /= numpy.max(y)
        f = interp1d(y,x)
        z = f(numpy.random.uniform(size=num))
        ages = cosmo.lookback_time(z)

# uniform distribution (default)
    else:
        print('using uniform distribution')
        ages = numpy.random.uniform(numpy.min(age_range), numpy.max(age_range), size=num)

    if sfh:
        print('reversing ages (SFH)')
        ages = numpy.max(ages)-ages

    return ages



def generateMasses(num,**kwargs):
    '''
    :Purpose: Generates a distribution of masses based on the defined input distribution. 

    Required Inputs:

    :param: num: number of masses to generate

    Optional Inputs:

    :param: mass_range: range of masses to draw from (default = [0.01,0.1]); can also specify ``minmass`` and ``maxmass``
    :param: distribution: can be a string set to one of the following to define the type of mass distribution to sample:
        * `uniform`: uniform distribution (default) 
        * `powerlaw` or `power-law`: single power-law distribution, P(M) ~ M\^-alpha. You must specify the parameter `alpha` or set ``distribution`` to TBD
        * `broken-powerlaw' or `broken-power-law: a broken power-law distribution; segments are specified by the parameters `alpha` (N array of numbers) for the slopes and `ranges` (N array of 2-element arrays) for the ranges over which these slopes occur; if the `scales` parameter is also included, the power-law segments are scaled by these factors; otherwise, the segments are forced to be continuous. You can also set ``distribution`` to `kroupa`
        * 'lognormal` or `log-normal`: log normal distribution, P(M) ~ exp(-0.5*(M-M0)\^2/sigmaM^2). You must specify the parameters `M0` and `sigmaM` or set ``distribution`` to `chabrier` (default parameters)
        * `kroupa`: broken power-law distribution with parameters from Kroupa et al. (XXXX): XXXX
        * `chabrier`: lognormal distribution with parameters from Chabrier et al. (XXX): XXXXX
    :param: distribution can also be set to a 2 x N array specifying the mass distribution; the first vector should be the masses for the distribution function and the second vector the distribution function itself
    :param: parameters: dictionary containing the parameters for the age distribution/star formation model being used; options include:
        * `alpha`: exponent for power-law distribution, or array of numbers giving power-law factors for broken power-law distribution
        * `range`: array of 2-element arrays specifying the masses (in units of solar masses) over which the broken-law slopes are defined
        * `scales`: array of numbers specifying relative scaling between the segments in the broken-law distribution
        * `M0` and `sigmaM: parameters for lognormal distribution in units of solar masses

    Output: 

    An array of masses drawn from the desired distribution in units of solar masses

    :Example:
    >>> import splat
    >>> import matplotlib.pyplot as plt
    >>> masses = splat.generateMasses(10000,distribution='power-law',parameters={'alpha': 0.5},mass_range=[0.01,0.08])
    }
    >>> plt.hist(masses)
    [histogram of masses in range 0.01-0.08 solar masses]    
    '''
    pass

    

def generatePopulation(**kwargs):

    parameters = {}

# draw ages - DONE
    age_kwargs = kwargs.get('age_parameters',{})
    parameters['age'] = generateAges(num,**age_kwargs)

# draw masses
    mass_kwargs = kwargs.get('mass_parameters',{})
#    parameters['mass'] = generateMasses(num,**age_kwargs)

# extract evolutionary model parameters
# NEED TO DEAL WITH OUT OF RANGE VALUES
    model_kwargs = kwargs.get('model_parameters',{})
    mp = modelParameters(mass=parameters['mass'],age=parameters['age'],**age_kwargs)
    parameters['gravity'] = mp['gravity']
    parameters['luminosity'] = mp['luminosity']
    parameters['radius'] = mp['radius']
    parameters['temperature'] = mp['temperature']

# determine spectral types from teff
# NOTE: NEED TO ALLOW FOR DIFFERENT RELATIONSHIPS AND AUTOMATICALLY DITCH NAN REGIONS
# ALSO DEAL WITH BAD VALUES
# COULD ALSO DO THIS WITH LUMINOSITIES
    sp = numpy.linspace(16,38,100)
    tf = numpy.array([splat.typeToTeff(spi)[0] for spi in sp])
    f = interp1d(sp,tf)
    parameters['spt'] = f(parameters['temperature'])

# add binary companions if desired

# assign binary orbital properties if desired

# assign sky positions if desired

# assign distances based on density profile if desired

# assign absolute, systemic and apparent magnitudes if desired

# assign age-dependent kinematics if desired

# assign proper and radial motions if desired

# assign apparent binary properties - current projected separation, astrometric offset, primary & secondary RV offsets - if desired

# assign metallicities (?) if desired

# visualize output?

    return parameters



###############################################################################
###################### TESTING FUNCTIONS #####################################
###############################################################################


def test_readmodel():
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

def test_evolve_basic():
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

def test_evolve_accuracy(modelname='baraffe', parameter='temperature',metallicity=0.,clouds='nocloud'):
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


def test_evolve_accuracy_plotting(xparam,yparam,modelname='baraffe',metallicity=0.,clouds='nocloud',**kwargs):
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

def test_ages(num,**kwargs):
    ages = generateAges(num,**kwargs)
    plt.hist(ages)
    plt.ylabel('Number')
    plt.xlabel('Age')
    plt.show()
    return True


if __name__ == '__main__':
#    test_readmodel()
#    test_evolve_basic()
#    test_evolve_accuracy(modelname='saumon',parameter='luminosity')
#    test_evolve_accuracy_plotting('age','temperature',modelname='baraffe',file='/Users/adam/projects/splat/code/testing/test_evolve_plotting.eps')
    test_ages(100000,distribution='exponential',minage=0.1,maxage=12.,parameters={'beta': -0.5})
