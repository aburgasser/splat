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
from astropy.io import ascii
import pandas
from math import isnan
import matplotlib.pyplot as plt
import numpy
from scipy.interpolate import interp1d 
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
    :Purpose: Driver function for modelParameters_, performs actual interpolation of evolutionary models. See SPLAT API for `modelParameters`_ for details.

    .. _`modelParameters` : api.html#splat_evolve.modelParameters

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
    :Purpose: Plots pairs of physical star parameters and optionally compares to evolutionary models. 

    Required Inputs:

    :param: parameters: dictionary or nested set of two arrays containing parameters to be plotted. For dictionary, keywords should include the `xparameter` and `yparameter` strings to be plotted. Values associated with keywords can be single numbers or arrays
    :param: xparameter: string corresponding to the key in the `parameters` dictionary to be plot as the x (independent) variable. 
    :param: yparameter: string corresponding to the key in the `parameters` dictionary to be plot as the y (dependent) variable. 


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

    A matplotlib plot object. Optinally, can also show plot on screen or output plot to a file.

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
    plt = plotModelParameters(p,xparam,yparam,model=model,**kwargs)
    return True

if __name__ == '__main__':
#    test_readmodel()
#    test_evolve_basic()
#    test_evolve_accuracy(modelname='saumon',parameter='luminosity')
    test_evolve_accuracy_plotting('age','temperature',modelname='baraffe',file='/Users/adam/projects/splat/code/testing/test_evolve_plotting.eps')

