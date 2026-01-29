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

# imports: internal
import copy
import glob
import os
import requests

# imports: external
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.constants as constants
from astropy.cosmology import Planck15, z_at_value
from astropy.io import ascii
import pandas
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy
from scipy.interpolate import griddata, interp1d
import scipy.stats as stats

# imports: splat
from splat.initialize import *
from splat.utilities import *
from splat.citations import shortRef


#################################
#                               #
# Evolutionary Model routines   #
#                               #
#################################

def loadEvolModel(modelname='baraffe2003',allparameters=False,verbose=False,**kwargs):
    '''
    Reads in the evolutionary model parameters for SPLAT-definied models, which are used 
    to interpolate parameters in `modelParameters()`_. See `modelInfo()`_ for which models are available.

    Parameters
    ----------

    modelname : str, default = 'baraffe2003'
        Name of the evolutionary model set to be used, must be in EVOLUTIONARY_MODELS dictionary

    allparameters : bool, default = False
        If True, returns all evolutionary model points (no constraints made)

    verbose : bool, default = False
        If True, return verbose output

    It is also possible to provide parameters to downselect model parameters (see example)

    Returns
    -------

    dict
        data structure containing evolutionary model information and a pandas dataframe for model parameters


    Examples
    --------

    A typical usage is a single call with the model name

    >>> import splat.evolve as spev
    >>> p = spev.loadEvolModel('morley2024',verbose=True)
    
    You are using evolutionary models from Morley et al. (2024)
    Setting parameter metallicity to default value 0.0
    Setting parameter cloud to default value hybrid-grav
    
    >>> p
    
    {'file': 'morley2024.csv',
     'citation': 'Morley et al. (2024)',
     'bibcode': '2024ApJ...975...59M',
     'url': 'https://zenodo.org/records/12735103',
     'altname': ['morley',
      'morley24',
      'diamondback',
      'diamondback24',
      'dback',
      'dback24',
      'mor24'],
     'default': {'metallicity': 0.0, 'cloud': 'hybrid-grav'},
     'n_points': 2985,
     'parameters':         mass      age  luminosity  temperature  gravity  radius  metallicity  \
     0     0.0005   0.0020     -5.6889       530.82    2.678   1.650          0.0   
     1     0.0005   0.0023     -5.7588       514.57    2.694   1.621          0.0   
     2     0.0005   0.0027     -5.8415       495.60    2.712   1.588          0.0   
     3     0.0005   0.0031     -5.9144       479.16    2.726   1.562          0.0   
     4     0.0005   0.0036     -5.9902       462.60    2.741   1.536          0.0   
     ...      ...      ...         ...          ...      ...     ...          ...   
     2980  0.0800   7.2000     -3.6800      2191.22    5.337   0.978          0.0   
     2981  0.0800   8.3000     -3.6801      2191.22    5.337   0.978          0.0   
     2982  0.0800   9.7000     -3.6810      2190.15    5.337   0.978          0.0   
     2983  0.0800  11.0000     -3.6813      2189.90    5.337   0.978          0.0   
     2984  0.0800  13.0000     -3.6815      2189.90    5.337   0.978          0.0   
     
                 cloud  
     0     hybrid-grav  
     1     hybrid-grav  
     2     hybrid-grav  
     3     hybrid-grav  
     4     hybrid-grav  
     ...           ...  
     2980  hybrid-grav  
     2981  hybrid-grav  
     2982  hybrid-grav  
     2983  hybrid-grav  
     2984  hybrid-grav  

    
    Note that default pre-selection was applied; you can specify downselect through model parameters;
    this example chooses all possible clouds but narrows metallicty range to 0 to 0.5

    >>> p = spev.loadEvolModel('morley2024',cloud='*',metallicity=[0,0.5],verbose=True)
    
    You are using evolutionary models from Morley et al. (2024)
    Downselecting on metallicity using constraint [0, 0.5]
    
    >>> len(p['parameters'])
    
    18155

    See Also
    --------

    `modelParameters()`_ : Determine evolutionary model parameters given two constraints
    `modelInfo()`_ : List information on available models

.. _`modelParameters()` : api.html#splat_evolve.modelParameters
.. _`modelInfo()` : api.html#splat_evolve.modelInfo

    '''

# check model
    model = checkDict(modelname,EVOLUTIONARY_MODELS)
    if isinstance(model,bool): 
        raise ValueError('\nDid not recognize model name {}; try {}'.format(modelname,list(EVOLUTIONARY_MODELS.keys())))
    if verbose==True: print('You are using evolutionary models from {}'.format(EVOLUTIONARY_MODELS[model]['citation']))

# read in full models
    dp = pandas.read_csv(os.path.join(EVOLUTIONARY_MODEL_FOLDER,EVOLUTIONARY_MODELS[model]['file']),comment='#',sep=',',header=0)
    dpsel = copy.deepcopy(dp)

# downselect models as desired
    if allparameters == False:
        dpselbck = copy.deepcopy(dp)
        for k in list(kwargs.keys()):
            if k in list(dp.columns) and kwargs[k]!='*':
                if verbose==True: print('Downselecting on {} using constraint {}'.format(k,kwargs[k]))
                if isinstance(kwargs[k],list) == True:
                    if len(kwargs[k]) == 2 and isNumber(kwargs[k][0])==True and isNumber(kwargs[k][1])==True:
                        dpsel = dpsel[dpsel[k] >= numpy.nanmin(kwargs[k])]
                        dpsel = dpsel[dpsel[k] <= numpy.nanmax(kwargs[k])]
                    else:
                        for i,x in enumerate(kwargs[k]):
                            if i == 0: dptmp = dpsel[dpsel[k] == x]
                            dptmp = pandas.concat([dptmp,dpsel[dpsel[k] == x]])
                        dpsel = copy.deepcopy(dptmp)
                else:
                    dpsel = dpsel[dpsel[k] == kwargs[k]]
                if len(dpsel) == 0:
                    if verbose==True: print('Warning: this constraint eliminated all models; ignoring')
                    dpsel = copy.deepcopy(dpselbck)
                else: dpselbck = copy.deepcopy(dpsel)
        dpsel.reset_index(inplace=True,drop=True)

# apply defaults if parameter not supplied
        dpselbck = copy.deepcopy(dpsel)
        for k in list(EVOLUTIONARY_MODELS[model]['default'].keys()):
            if k not in list(kwargs.keys()): 
                if verbose==True: print('Setting parameter {} to default value {}'.format(k,EVOLUTIONARY_MODELS[model]['default'][k]))
                dpsel = dpsel[dpsel[k] == EVOLUTIONARY_MODELS[model]['default'][k]]
            if len(dpsel) == 0:
                if verbose==True: print('Warning: this constraint eliminated all models; ignoring')
                dpsel = copy.deepcopy(dpselbck)
            else: dpselbck = copy.deepcopy(dpsel)
        dpsel.reset_index(inplace=True,drop=True)

# convert into dictionary
    mparam = copy.deepcopy(EVOLUTIONARY_MODELS[model])
    mparam['n_points'] = len(dpsel)
    mparam['parameters'] = dpsel
    return mparam



def modelInfo(modelname,file='',verbose=False):
    '''
    Returns information on a model or set of models

    Parameters
    ----------

    modelname : str
        Name of the evolutionary model set to be used, must be in EVOLUTIONARY_MODELS dictionary

    file : str, default = ''
        Set to the full path file name where you wnat to save the model information

    verbose : bool, default = False
        If True, return verbose output

    Returns
    -------

    str
        string containing model in formation, either to screen or to specified file


    Examples
    --------

    >>> import splat.evolve as spev
    >>> p = spev.modelInfo('morley2024')

    Evolutionary models morley2024
    Morley et al. (2024) 
    Download the original models at https://zenodo.org/records/12735103
    Parameter values:
        mass: 0.001 to 0.080
        age: 0.002 to 13.000
        luminosity: -7.811 to -2.432
        temperature: 200.040 to 2499.880
        gravity: 2.566 to 5.509
        radius: 0.741 to 3.170
        metallicity: -0.500 to 0.500
        cloud: hybrid,hybrid-grav,nc

    See Also
    --------

    `modelParameters()`_ : Determine evolutionary model parameters given two constraints
    `loadEvolModel()`_ : Loads the parameter of an evolutionary model

.. _`modelParameters()` : api.html#splat_evolve.modelParameters
.. _`loadEvolModel()` : api.html#splat_evolve.loadEvolModel

    '''

# check model
    if not isinstance(modelname,str):
        print('\nInput parameter should be a string (name of model); you entered {}'.format(modelname))
        return

    model = checkDict(modelname,EVOLUTIONARY_MODELS)
    if isinstance(model,bool): 
        print('\nDid not recognize model name {}; try {}'.format(modelname,list(EVOLUTIONARY_MODELS.keys())))
        return

# read in model information
    par = loadEvolModel(model,allparameters=True)

# prepare output string
    output = 'Evolutionary models {}'.format(modelname)
    cite = '\n '
    if 'citation' in list(par.keys()): output+='\n{}'.format(par['citation'])
    if 'bibcode' in list(par.keys()): output+=' {}'.format(shortRef(par['bibcode']))
    if 'url' in list(par.keys()): output+='\nDownload the original models at {}'.format(EVOLUTIONARY_MODELS[model]['url'])
    output+='\nParameter values:'
    dp = par['parameters']
    for k in list(dp.keys()):
        if isNumber(dp.loc[0,k]):
            output+='\n\t{}: {:.3f} to {:.3f}'.format(k,numpy.nanmin(dp[k]),numpy.nanmax(dp[k]))
        else:
            vals = list(set(list(dp[k])))
            vals.sort()
            x = vals[0]
            for v in vals[1:]: x+=',{}'.format(v)
            output+='\n\t{}: {}'.format(k,x)
    if file != '':
        f = open(file,'w')
        f.print(output)
        f.close()
        if verbose==True: ('Model information written to {}'.format(file))
    else:
        print(output)
    return 





def _modelParametersSingle(orig_model,verbose=False,**kwargs):
    '''
    :Purpose: Driver function for modelParameters_, performs actual interpolation of evolutionary models. See SPLAT API for `modelParameters()`_ for details.

    .. _`modelParameters()` : api.html#splat_evolve.modelParameters

    '''

    logparam = ['age','mass','temperature','radius']
    computeparam = ['temperature','gravity','radius','luminosity']
    keywords = list(kwargs.keys())

# check that model is passed correctly - should be a pandas dataframe
    model = copy.deepcopy(orig_model)
    if not isinstance(model,pandas.core.frame.DataFrame):
        if isinstance(model,dict):
            if 'parameters' in list(model.keys()): model = model['parameters']
            if not isinstance(model,pandas.core.frame.DataFrame):
                try: model = pandas.DataFrame(model)
                except: pass
        if isinstance(model,str):
            model = loadEvolModel(orig_model,**kwargs)['parameters']
    if not isinstance(model,pandas.core.frame.DataFrame):
        raise ValueError('Could not process input model variable, which is type {}\n{}'.format(type(orig_model),orig_model))

#     try: model = args[0]
#     except IndexError: 
#         model = loadEvolModel('baraffe03')
#         print('\nWarning: using Baraffe et al. (2003) models by default\n')

# rescale some model parameters to allow for logarithmic interpolation
    lmodel = copy.deepcopy(model)
# convert to logarithmic values
    for k in logparam:
        lmodel[k] = numpy.array([numpy.log10(m) for m in lmodel[k]])
    ages = list(set(lmodel['age']))
    ages.sort()
    ages = numpy.array(ages)

# prep output parameters
    params = {}
    for k in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):
        params[k] = numpy.nan
        if k in keywords:
            try: f = float(kwargs[k])
            except: raise ValueError('\nInput parameter {} must be a single number, not {}\n'.format(k,kwargs[k]))
            finally: params[k] = f

# set some to log
    for k in logparam:
        if k in list(params.keys()): params[k] = numpy.log10(params[k])

# check that you have at least 2 non-nan values
    nparam = 0
    for k in params.keys():
        if numpy.isnan(params[k]) == False: nparam=nparam+1
    if nparam < 2: 
        if verbose==True: print('Warning: need at least 2 finite parameters; you passed {}'.format(params))
        return params

    input_type = 'mass_age'
    Ag, Ma, Te, Le, Ge, Re, P = [],[],[],[],[],[],[]


    if kwargs.get('debug',False) == True: print(lmodel)

############### UNKNOWN MASS AND AGE - INTERPOLATE AGE FROM OTHER PARAMETERS #################
# for each age, interpolate mass as a function of first parameter 
# then interpolate second parameter as a function of mass
# obtain second parameter as a function of age; then interpolate the model ages as a function of 
# the second parameter and evaluate for known parameter to get age
###############################################################################


# REVISED METHOD USING GRIDDATA    
    if numpy.isnan(params['mass']) == True and numpy.isnan(params['age']) == True:

        input_type = 'two_params'
        for k in computeparam:
            if numpy.isnan(params[k]) == False: P.append([k, params[k]])

# create a grid to extract the mass and age
        points = []
        Ag, Ma = [], []

        for i,age in enumerate(ages):
            lmodelpar1 = lmodel[P[0][0]][lmodel['age'] == age]
            lmodelpar2 = lmodel[P[1][0]][lmodel['age'] == age]
            lmodelmass = list(lmodel['mass'][lmodel['age'] == age])
            if numpy.nanmin(lmodelpar1) <= P[0][1] <= numpy.nanmax(lmodelpar1) \
                and numpy.nanmin(lmodelpar2) <= P[1][1] <= numpy.nanmax(lmodelpar2):
                for j,m in enumerate(lmodelpar1): 
                    points.append((m,lmodelpar2.ioc[j]))
                    Ag.append(age)
                Ma.extend(lmodelmass)

        if kwargs.get('debug',False) == True:
            pts = numpy.array(points)
            print('\n')
            print(pts,pts.shape)
            print('\n')
            print(Ag,len(Ag))
            print('\n')
            print(Ma,len(Ma))
            print('\n')
            print(P[0][1],P[1][1])

        try: 
            params['age'] = griddata(numpy.array(points),Ag,numpy.array((P[0][1],P[1][1])),method='linear')[0]
        except: 
            if kwargs.get('verbose',False) == True: print('\nFailed in 2 parameter determination\n')
            params['age'] = numpy.nan
        try: 
            params['mass'] = griddata(numpy.array(points),Ma,numpy.array((P[0][1],P[1][1])),method='linear')[0]
        except: 
            if kwargs.get('verbose',False) == True: print('\nFailed in 2 parameter determination\n')
            params['mass'] = numpy.nan


        if kwargs.get('debug',False) == True: print('\nMass and Age unknown; determined age to be {} and mass to be {}'.format(params['age'],params['mass']))
        Ge, Ag, Ma = [], [], []


################ UNKNOWN AGE BUT KNOWN MASS AND ONE OTHER PARAMETER ###########
# interpolate second parameter as a function of mass for each of the age models and evaluate for known mass
# interpolate the model ages as a fucntion of these parameters and evaluate for known parameter
###############################################################################
#    if params['age'] == 0. and params['mass'] != 0. and not numpy.isnan(params['mass']):
    if numpy.isnan(params['age']) == True and numpy.isnan(params['mass']) == False:

        if input_type != 'two_params': 
            input_type = 'one_param'
            for k in computeparam:
                if numpy.isnan(params[k]) == False: P.append([k, params[k]])
            if len(P) == 0 or len(P) > 1:
                for k in list(params.keys()):
                    print('{}: {}'.format(k,params[k]))
                print(P)
                raise ValueError('\nProblem with one_param interpolation\n')

        for i,age in enumerate(ages):
            lmodelmass = list(lmodel['mass'][lmodel['age'] == age])
            lmodelpar1 = lmodel[P[0][0]][lmodel['age'] == age]
            if numpy.nanmin(lmodelmass) <= params['mass'] <= numpy.nanmax(lmodelmass):
                Ag.append(age)
                f = interp1d(lmodelmass, list(lmodelpar1))
                Ge.append(f(params['mass']))

        try: 
            f = interp1d(Ge, Ag)
            params['age'] = f(P[0][1])
        except: 
            if kwargs.get('verbose',False) == True: print('\nFailed in age + parameter determination\n')
            params['age'] = numpy.nan

        if kwargs.get('debug',False) == True: print('\nMass known and Age unknown; determined age to be {}'.format(10.**params['age']))

        Ge, Ag = [], []


################ KNOWN AGE BUT UNKNOWN MASS AND ONE OTHER PARAMETER ###########
# generate mass as function of second parameter interpolated between two closest age models
# evaluate mass(parameter) (resulting in both mass and age as knowns)
###############################################################################

#    if params['age'] != 0. and params['mass'] == 0. and not numpy.isnan(params['age']):
    if numpy.isnan(params['age']) == False and numpy.isnan(params['mass']) == True:
        if kwargs.get('debug',False) == True: print(params)

        if input_type != 'two_params' and input_type != 'one_param': 
            input_type = 'one_param'
            for k in computeparam:
                if numpy.isnan(params[k]) == False: P.append([k, params[k]])
            if len(P) == 0 or len(P) > 1:
                for k in list(params.keys()):
                    print('{}: {}'.format(k,params[k]))
                print(P)
                raise ValueError('\nProblem with one_param interpolation\n')


        if params['age'] < numpy.nanmin(ages) or params['age'] > numpy.nanmax(ages):
                if kwargs.get('verbose',False)==True: print('\nAge of {} is outside range of models, {} to {}\n'.format(10.**params['age'],10.**numpy.min(ages),10**numpy.max(ages)))
                params['mass'] = numpy.nan

        else:
            adiff = [params['age']-a for a in ages]
            ai = numpy.argmin(numpy.abs(adiff))
            if adiff[ai] < 0: ai=ai-1
            aii = numpy.nanmin([ai+1,len(ages)-1])
            lmodelmass1 = list(lmodel['mass'][lmodel['age'] == ages[ai]])
            lmodelmass2 = list(lmodel['mass'][lmodel['age'] == ages[aii]])
            lmodelpar1 = lmodel[P[0][0]][lmodel['age'] == ages[ai]]
            lmodelpar2 = lmodel[P[0][0]][lmodel['age'] == ages[aii]]
            for i,m in enumerate(lmodelmass1):
                if m in lmodelmass2:
                    Ma.append(m)
                    aj = numpy.argmin(numpy.abs([a-m for a in lmodelmass2]))
                    vals = [lmodelpar1.iloc[i],lmodelpar2.iloc[aj]]
                    f = interp1d([ages[ai],ages[aii]],vals)
                    Ge.append(f(params['age']))
            try:
                f = interp1d(Ge, Ma)
                params['mass'] = f(P[0][1])
            except:
                if kwargs.get('verbose',False) == True: print('\nFailed in mass + parameter determination\n')
                params['mass'] = numpy.nan

        if kwargs.get('debug',False) == True: print('\nMass unknown and Age known; determined mass to be {}'.format(10.**params['mass']))

        Ma, Ge = [],[]


###################### KNOWN MASS AND AGE #####################################
# FINALLY use a simple grid interpolation
###############################################################################
    if numpy.isnan(params['age']) == False and numpy.isnan(params['mass']) == False:

        if kwargs.get('debug',False) == True: print(params)
  
        for k in computeparam:
            if numpy.isnan(params[k]) == True:
                params[k] = griddata(numpy.transpose([lmodel['mass'],lmodel['age']]),lmodel[k],(params['mass'],params['age']),method='linear')
#                if numpy.isfinite(params[k]) == True: 
#                    params[k] = params[k]*EVOLUTIONARY_MODEL_PARAMETERS[k]['unit']
                if isinstance(params[k],numpy.ndarray): params[k] = float(params[k])

        if kwargs.get('debug',False) == True: print('\nDetermined parameters: {}'.format(params))

# return to linear
        for k in logparam: params[k] = 10**(params[k])

# assign parameters that are constant throughout modelset
        for k in list(model.keys()):
            samp = list(set(list(model[k])))
            if len(samp)==1: params[k] = samp[0]

        return params


# something failed	  
    else:
#        print(params)
        # for e in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):
        #     params[e] = numpy.nan
        if kwargs.get('verbose',False) == True: 
            print('\nParameter set is not covered by input model\n')
        return params
      




def modelParameters(modelname='baraffe2003',verbose=False,debug=False,alt=False,**kwargs):
    '''
    :Purpose: Retrieves the evolutionary model parameters given two of the following parameters: mass, age, temperature, luminosity, gravity, or radius. The inputs can be individual values or arrays.  Using the input parameters, the associated evolutionary model parameters are computed through log-linear interpolation of the original model grid. Parameters that fall outside the grid return nan.

    Required Inputs:

    :param: modelname: Either a string of the name of the evolutionary model set, which can be one of `baraffe` (default), `burrows`, or `saumon`; or a dictionary output from `loadEvolModel()`_ containing model parameters. 
    
    and two (2) of the following:

    :param: mass: input value or list of values for mass (can also be `masses` or `m`)
    :param: age: input value or list of values for age (can also be `ages`, `time` or `a`)
    :param: temperature: input value or list of values for temperature (can also be `temperatures`, `teff`, `temp` or `t`)
    :param: gravity: input value or list of values for gravity (can also be `gravities`, `grav`, `logg` or `g`)
    :param: luminosity: input value or list of values for luminosity (can also be `luminosities`, `lum`, `lbol` or `l`)
    :param: radius: input value or list of values for radius (can also be `radii`, `rad` and `r`)

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
    for k in ['model','set','ref','reference','mset']: modelname=kwargs.get(k,modelname)
    if isinstance(modelname,str): 
        tmp = checkDict(modelname,EVOLUTIONARY_MODELS)
        if isinstance(tmp,bool): raise ValueError('Model set {} not in SPLAT; try {}'.format(modelname,list(EVOLUTIONARY_MODELS.keys())))
        modelname = copy.deepcopy(tmp)
# fix for downselect
        loadkwargs = {}        
        for k in EVOLUTIONARY_MODELS[modelname]['default']:
            if k in list(kwargs.keys()): loadkwargs[k] = kwargs[k]
        model = loadEvolModel(modelname,verbose=verbose,**loadkwargs)
    else: 
        model = copy.deepcopy(modelname)
        modelname = 'user-supplied'
# set up model parameters as pandas dataframe (might be provided)
    modelpars = copy.deepcopy(model)
    if isinstance(model,dict):
        if 'parameters' in list(model.keys()): modelpars = model['parameters']
        if 'mass' in list(model.keys()) and 'age' in list(model.keys()):
            try: modelpars = pandas.DataFrame(model)
            except: pass
        if isinstance(modelpars,dict):
            if 'mass' in list(modelpars.keys()) and 'age' in list(modelpars.keys()):
                try: modelpars = pandas.DataFrame(modelpars)
                except: pass
    if not isinstance(modelpars,pandas.core.frame.DataFrame):
        raise ValueError('Something went wrong in loading in models; you provided {}'.format(modelname))
    if verbose == True: print('Using model {}'.format(modelname))

# do some search parameter name replacement - parameters should be defined in EVOLUTIONARY_MODEL_PARAMETERS
    keywords = list(kwargs.keys())
    mkwargs = {}
    for k in keywords:
        tmp = checkDict(k,EVOLUTIONARY_MODEL_PARAMETERS)
        if isinstance(tmp,bool) == False: mkwargs[tmp] = kwargs[k]
        tmp = checkDict(k[:-1],EVOLUTIONARY_MODEL_PARAMETERS) # remove s
        if isinstance(tmp,bool) == False: mkwargs[tmp] = kwargs[k]

# do some model parameter replacement - parameters should be defined in EVOLUTIONARY_MODEL_PARAMETERS
    for k in list(modelpars.keys()):
        tmp = checkDict(k,EVOLUTIONARY_MODEL_PARAMETERS)
        if isinstance(tmp,bool) == False: modelpars.rename({k:tmp},inplace=True)
        tmp = checkDict(k[:1],EVOLUTIONARY_MODEL_PARAMETERS) # remove s
        if isinstance(tmp,bool) == False and k in list(modelpars.columns): modelpars.rename({k:tmp},inplace=True)

# # add in default options
#     for k in list(modelpars.keys()):
#         if k in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()) and k not in list(mkwargs.keys()):
#             mkwargs[k] = EVOLUTIONARY_MODEL_PARAMETERS


# downselect models based on discrete values
    # for k in list(model.keys()):
    #     if isinstance(model.loc[k,0],str):

    #     if k in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):
    #         if isinstance(EVOLUTIONARY_MODEL_PARAMETERS[k]

# determine length of input arrays and save off units 
    pkeys = list(mkwargs.keys())
#    punits = {}
    if len(pkeys) < 2:
        raise ValueError('\nNeed at least two parameters provided; only {} given: {}'.format(len(pkeys),pkeys))
    for p in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):
        if p in pkeys:
#            punits[p] = EVOLUTIONARY_MODEL_PARAMETERS[p]['unit']
#            print(p,type(mkwargs[p]))
            if isUnit(mkwargs[p]):
#                punits[p] = mkwargs[p].unit
                mkwargs[p] = mkwargs[p].value
            if not isinstance(mkwargs[p],list) and not isinstance(mkwargs[p],numpy.ndarray):
                mkwargs[p] = [mkwargs[p]]
            if isUnit(mkwargs[p]):
 #               punits[p] = mkwargs[p].unit
                mkwargs[p] = mkwargs[p].value
            if isUnit(mkwargs[p][0]):
                mkwargs[p] = [x*value for x in mkwargs[p]]
            
# fill in extra parameters if the number provided is unequal (i.e., default to "most" parameters)            
    numberValues = numpy.nanmax([len(mkwargs[p]) for p in list(mkwargs.keys())])
    for p in list(mkwargs.keys()):
        while len(mkwargs[p]) < numberValues: mkwargs[p].append(mkwargs[p][-1])

# now loop through each parameter set to determine remaining parameters
    inparams = {'alt':alt,'debug':debug,'verbose':verbose}
    outparams = {}
    for i in range(numberValues):
        for p in pkeys: inparams[p] = mkwargs[p][i]
        par = _modelParametersSingle(modelpars,**inparams)
        if i==0: 
            for p in list(par.keys()): outparams[p] = [par[p]]
        else: 
            for p in list(outparams.keys()): outparams[p].append(par[p])
    if debug==True: print(outparams)

# remove lists if only one parameter set is being calculated
    kk = list(outparams.keys())
    if len(outparams[kk[0]]) == 1:
        for k in kk: outparams[k] = outparams[k][0]
            

# add back units
    for k in kk: 
        if k in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):
            if EVOLUTIONARY_MODEL_PARAMETERS[p]['unit'] != None and isUnit(outparams[k])==False and isNumber(outparams[k])==True:
                outparams[k] *= EVOLUTIONARY_MODEL_PARAMETERS[p]['unit']

    return outparams




def plotModelParameters(parameters,xparam,yparam,verbose=False,**kwargs):
    '''
    :Purpose: 

        Plots pairs of physical star parameters and optionally compares to evolutionary model tracks. 

    :Required Inputs:

        :param: parameters: dictionary or nested set of two arrays containing parameters to be plotted. For dictionary, keywords should include the `xparameter` and `yparameter` strings to be plotted. Values associated with keywords can be single numbers or arrays
        :param: xparam: string corresponding to the key in the `parameters` dictionary to be plot as the x (independent) variable. 
        :param: yparam: string corresponding to the key in the `parameters` dictionary to be plot as the y (dependent) variable. 

    :Optional Inputs:

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

    :Output: 

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
    noplot=False
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
        noplot = True 
        if verbose==True: print('\nx-axis parameter {} is not one that can be plotted'.format(xparam))

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
        noplot = True 
        if verbose==True: print('\ny-axis parameter {} is not one that can be plotted'.format(xparam))

# can't plot - just return
    if noplot==True:
        if verbose==True: print('\nUnable to plot; returning')
        return None

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
            model = 'baraffe03'
        else:
            model = kwargs.get('model')
        try:
            if type(model) is not dict: model = loadEvolModel(model,**kwargs)
        except:
            print('\nProblem in reading in original models\n')
            kwargs['showmodel'] = False

# models tracks (mass by default)
        tracks = kwargs.get('tracks',['mass'])
        if isinstance(tracks,bool):
            if tracks == True: 
                tracks = ['mass','age']
            else:
                tracks = ['']
        if tracks == 'all': tracks = ['mass','age']
        if not isinstance(tracks,list): tracks = [tracks]
        tracks = [t.lower() for t in tracks]

        if 'mass' in tracks:
            xvals, yvals, masses = [], [], []
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

            for i,x in enumerate(xvals):
                if xlogflag == True and ylogflag == True:
                    plt.loglog(xvals[i],yvals[i],color='grey',linestyle='-')
                elif xlogflag == False and ylogflag == True:
                    plt.semilogy(xvals[i],yvals[i],color='grey',linestyle='-')
                elif xlogflag == True and ylogflag == False:
                    plt.semilogx(xvals[i],yvals[i],color='grey',linestyle='-')
                else:
                    plt.plot(xvals[i],yvals[i],color='grey',linestyle='-')

# models tracks trace isochrones
        if 'age' in tracks:
            xvals, yvals = [], []
            tvals = model['age']
# fix to account for unequal lengths of model values
            maxlen = numpy.max([len(a) for a in model['mass']])
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
                    plt.loglog(xvals[i],yvals[i],color='grey',linestyle='--')
                elif xlogflag == False and ylogflag == True:
                    plt.semilogy(xvals[i],yvals[i],color='grey',linestyle='--')
                elif xlogflag == True and ylogflag == False:
                    plt.semilogx(xvals[i],yvals[i],color='grey',linestyle='--')
                else:
                    plt.plot(xvals[i],yvals[i],color='grey',linestyle='--')

# add labels
    plt.xlabel(kwargs.get('xlabel','{} ({})'.format(xmparam,EVOLUTIONARY_MODEL_PARAMETERS[xmparam]['unit'])))
    plt.ylabel(kwargs.get('ylabel','{} ({})'.format(ymparam,EVOLUTIONARY_MODEL_PARAMETERS[yparam]['unit'])))
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


##########################
# old code to be removed #
##########################

def loadEvolModel_old(*model,**kwargs):
    '''
    :Purpose: Reads in the evolutionary model parameters for the models listed below, which are used to interpolate parameters in `modelParameters()`_. 

    .. _`modelParameters()` : api.html#splat_evolve.modelParameters

    Available models are:

        - **baraffe** : Models from `Baraffe et al. (2003) <http://adsabs.harvard.edu/abs/2003A&A...402..701B>`_ for 1 Myr < age < 10 Gyr, 0.005 Msol < mass < 0.1 Msol, and solar metallicity (COND dust prescription)
        - **baraffe15** : Models from `Baraffe et al. (2015) <http://adsabs.harvard.edu/abs/2015A&A...577A..42B>`_ for 1 Myr < age < 10 Gyr, 0.01 Msol < mass < 1.4 Msol, and solar metallicity
        - **burrows** : Models from `Burrows et al. (2001) <http://adsabs.harvard.edu/abs/2001RvMP...73..719B>`_ for 1 Myr < age < 10 Gyr, 0.005 Msol < mass < 0.2 Msol, and solar metallicity
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

    Models are contained in SPLAT's EvolutionaryModels folder.

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
    except TypeError: raise TypeError('\nInput parameter must be a string')
    except IndexError: model = 'baraffe15'

    m = checkEvolutionaryModelName(model)
    if m == False: raise ValueError('\nDid not recognize model name {}; try {}'.format(model,list(EVOLUTIONARY_MODELS.keys())))
    model = m
    if kwargs.get('verbose',False): print('You are using evolutionary models from {}'.format(EVOLUTIONARY_MODELS[model]['name']))

# read in models
    files = glob.glob(os.path.join(EVOLUTIONARY_MODEL_FOLDER,model,'{}*.txt'.format(model)))
    if model == 'saumon08':

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
        if metallicity=='-0.3' or metallicity=='0.3':
            C = 'nc'
        if isinstance(cloud,int) or isinstance(cloud,float):
            cloud = 'f{:1d}'.format(int(cloud))
        if cloud.lower() == 'hybrid':
            C = 'hybrid'
        elif cloud.lower() == 'f2':
            C = 'f2'
        elif cloud.lower() == 'nc' or cloud.lower() == 'nocloud' or cloud.lower() == 'noclouds' or cloud.lower() == 'cloud-free':
            C = 'nc'
        else:
            raise ValueError('\nCould not recognize cloud choice for Saumon model: must be cloud-free, hybrid or f2, not {}\n'.format(cloud))

        files = glob.glob(os.path.join(EVOLUTIONARY_MODEL_FOLDER,model,'/{}_{}_{}*.txt'.format(model,Z,C)))

# read in parameters
    ages = [(f.split('_')[-1]).replace('.txt','') for f in files]

    mparam = {}
    mparam['name'] = model
    for ep in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):
        mparam[ep] = []


    for f in files:
        try:
            dp=pandas.read_csv(os.path.normpath(f),comment='#',sep=',',header=0)
            for ep in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):
                mparam[ep].append(dp[ep].values)



# this is done in case models are not local - NOTE: currently just throwing an error
        except:
            raise ValueError('Could not find model file {} locally; aborting'.format(f))
#            try:
#                print('Could not read in model file {} locally; trying online'.format(mfile))
#                data =ascii.read(requests.get(SPLAT_URL+EVOLUTIONARY_MODEL_FOLDER+mfile).content,comment='#',delimiter='\t')
#            except: 
#                raise ValueError('Could not find model file {} locally or online; aborting'.format(mfile))

    mparam['age'] = [float(i)/1000. for i in ages]

    return mparam



def _modelParametersSingle_old(*args, **kwargs):
    '''
    :Purpose: Driver function for modelParameters_, performs actual interpolation of evolutionary models. See SPLAT API for `modelParameters()`_ for details.

    .. _`modelParameters()` : api.html#splat_evolve.modelParameters

    '''

    keywords = list(kwargs.keys())

# check that model is passed correctly
    try: model = args[0]
    except IndexError: 
        model = loadEvolModel('baraffe03')
        print('\nWarning: using Baraffe et al. (2003) models by default\n')

# retool models to allow for logarithmic interpolation
    lmodel = copy.deepcopy(model)
# convert to logarithmic values
    lmodel['age'] = [numpy.log10(m) for m in lmodel['age']]
    for i in range(len(lmodel['age'])):
        lmodel['mass'][i] = [numpy.log10(m) for m in lmodel['mass'][i]]
        lmodel['temperature'][i] = [numpy.log10(m) for m in lmodel['temperature'][i]]
        lmodel['radius'][i] = [numpy.log10(m) for m in lmodel['radius'][i]]

# prep output parameters
    params = {}
#    for e in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):

    for e in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):
        params[e] = 0.
        if e in keywords:
            try: f = float(kwargs[e])
            except: raise ValueError('\nInput parameter {} must be a single number, not {}\n'.format(e,kwargs[e]))
            finally: params[e] = f

    input_type = 'mass_age'
    Ag, Ma, Te, Le, Ge, Re, P = [],[],[],[],[],[],[]


    if kwargs.get('debug',False) == True: print(lmodel)

############### UNKNOWN MASS AND AGE - INTERPOLATE AGE FROM OTHER PARAMETERS #################
# for each age, interpolate mass as a function of first parameter and then second parameter as a function of mass
# and obtain second parameter as a function of age; then interpolate the model ages as a function of 
# the second parameter and evaluate for known parameter to get age
###############################################################################

# REVISED METHOD USING GRIDDATA    
    if (params['mass'] == 0.) and (params['age'] == 0.) and kwargs.get('alt',False) == False:

        Ma1, Ma2 = [], []
        input_type = 'two_params'
        if params['temperature'] != 0.:
            P.append(['temperature', numpy.log10(params['temperature'])])
        if params['gravity'] != 0.:
            P.append(['gravity', params['gravity']])
        if params['radius'] != 0.:
            P.append(['radius', numpy.log10(params['radius'])])
        if params['luminosity'] != 0.:
            P.append(['luminosity', params['luminosity']])

# create a grid to produce extract the mass and age
        points = []
        Ag, Ma = [], []

        for i,age in enumerate(lmodel['age']):
            if numpy.nanmin(lmodel[P[0][0]][i]) <= P[0][1] <= numpy.nanmax(lmodel[P[0][0]][i]) \
                and numpy.nanmin(lmodel[P[1][0]][i]) <= P[1][1] <= numpy.nanmax(lmodel[P[1][0]][i]):
                for j,m in enumerate(lmodel[P[0][0]][i]): 
                    points.append((lmodel[P[0][0]][i][j],lmodel[P[1][0]][i][j]))
                    Ag.append(age)
                Ma.extend(lmodel['mass'][i])

        if kwargs.get('debug',False) == True:
            pts = numpy.array(points)
            print('\n')
            print(pts,pts.shape)
            print('\n')
            print(Ag,len(Ag))
            print('\n')
            print(Ma,len(Ma))
            print('\n')
            print(P[0][1],P[1][1])

        try: 
            params['age'] = 10.**(griddata(numpy.array(points),Ag,numpy.array((P[0][1],P[1][1])),method='linear')[0])
        except: 
            params['age'] = float('nan')
        try: 
            params['mass'] = 10.**(griddata(numpy.array(points),Ma,numpy.array((P[0][1],P[1][1])),method='linear')[0])
        except: 
            params['mass'] = float('nan')

        if kwargs.get('debug',False) == True: print('\nMass and Age unknown; determined age to be {} and mass to be {}'.format(params['age'],params['mass']))
        Ge, Ag, Ma = [], [], []


# OLD WAY OF DOING THIS

    if (params['mass'] == 0.) and (params['age'] == 0.) and kwargs.get('alt',False) == True:

        input_type = 'two_params'
        if params['temperature'] != 0.:
            P.append(['temperature', numpy.log10(params['temperature'])])
        if params['gravity'] != 0.:
            P.append(['gravity', params['gravity']])
        if params['radius'] != 0.:
            P.append(['radius', numpy.log10(params['radius'])])
        if params['luminosity'] != 0.:
            P.append(['luminosity', params['luminosity']])

# compute for each age the appropriate mass for each parameter
        for i,age in enumerate(lmodel['age']):
            if numpy.nanmin(lmodel[P[0][0]][i]) <= P[0][1] <= numpy.nanmax(lmodel[P[0][0]][i]) \
                and numpy.nanmin(lmodel[P[1][0]][i]) <= P[1][1] <= numpy.nanmax(lmodel[P[1][0]][i]):
                Ag.append(age)
                f = interp1d(lmodel[P[0][0]][i], lmodel['mass'][i])
                Ma = f(P[0][1])
                f = interp1d(lmodel['mass'][i], lmodel[P[1][0]][i])
                Ge.append(f(Ma))

        try: 
            f = interp1d(Ge, Ag)
            params['age'] = 10.**f(P[1][1])
        except: 
            params['age'] = float('nan')
            if kwargs.get('debug',False) == True: 
                print('\nFailed:\nP = {}\nGe = {}\nAg = {}'.format(P,Ge,Ag))
                for a in Ag:
                    i = numpy.argmin(numpy.abs(numpy.array(lmodel['age'])-a))
                    print(a,i,numpy.nanmin(lmodel[P[0][0]][i]),P[0][1],numpy.nanmax(lmodel[P[0][0]][i]),numpy.nanmin(lmodel[P[1][0]][i]),P[1][1],numpy.nanmax(lmodel[P[1][0]][i]))
                    f = interp1d(lmodel[P[0][0]][i], lmodel['mass'][i])
                    Ma = f(P[0][1])
                    print(lmodel[P[0][0]][i],lmodel['mass'][i],Ma)
                    f = interp1d(lmodel['mass'][i], lmodel[P[1][0]][i])
                    print(lmodel[P[1][0]][i],lmodel['mass'][i],f(Ma))


        if kwargs.get('debug',False) == True: print('\nMass and Age unknown; determined age to be {}'.format(params['age']))
        Ge, Ag, Ma = [], [], []

################ UNKNOWN AGE BUT KNOWN MASS AND ONE OTHER PARAMETER ###########
# interpolate second parameter as a function of mass for each of the age models and evaluate for known mass
# interpolate the model ages as a fucntion of these parameters and evaluate for known parameter
###############################################################################
    if params['age'] == 0. and params['mass'] != 0. and not numpy.isnan(params['mass']):

        if input_type != 'two_params': 
            input_type = 'one_param'
            if params['temperature'] != 0.:
                P.append(['temperature', numpy.log10(params['temperature'])])
            elif params['gravity'] != 0.:
                P.append(['gravity', params['gravity']])
            elif params['radius'] != 0.:
                P.append(['radius', numpy.log10(params['radius'])])
            elif params['luminosity'] != 0.:
                P.append(['luminosity', params['luminosity']])
            else:
                for k in list(params.keys()):
                    print('{}: {}'.format(k,params[k]))
                print(P)
                raise ValueError('\nProblem with one_param interpolation\n')

        for i,age in enumerate(lmodel['age']):
            if numpy.nanmin(lmodel['mass'][i]) <= numpy.log10(params['mass']) <= numpy.nanmax(lmodel['mass'][i]):
                Ag.append(age)
                f = interp1d(lmodel['mass'][i], lmodel[P[0][0]][i])
                Ge.append(f(numpy.log10(params['mass'])))

        try: 
            f = interp1d(Ge, Ag)
            params['age'] = 10.**f(P[0][1])
        except: 
            if kwargs.get('verbose',False) == True: print('\nFailed in age + parameter determination\n')
            params['age'] = float('nan')

        if kwargs.get('debug',False) == True: print('\nMass known and Age unknown; determined age to be {}'.format(params['age']))

        Ge, Ag = [], []


################ UNKNOWN AGE BUT KNOWN MASS AND ONE OTHER PARAMETER ###########
# generate mass as function of second parameter interpolated between two closest age models
# evaluate mass(parameter) (resulting in both mass and age as knowns)
###############################################################################

    if params['age'] != 0. and params['mass'] == 0. and not numpy.isnan(params['age']):

        if kwargs.get('debug',False) == True: print(params)

        if input_type != 'two_params' and input_type != 'one_param': 
            input_type = 'one_param'
            if params['temperature'] != 0.:
                P.append(['temperature', numpy.log10(params['temperature'])])
            elif params['gravity'] != 0.:
                P.append(['gravity', params['gravity']])
            elif params['radius'] != 0.:
                P.append(['radius', numpy.log10(params['radius'])])
            elif params['luminosity'] != 0.:
                P.append(['luminosity', params['luminosity']])
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
            if adiff[ai] < 0: ai=ai-1
            aii = numpy.nanmin([ai+1,len(lmodel['mass'])-1])
            for i,m in enumerate(lmodel['mass'][ai]):
                if m in lmodel['mass'][aii]:
                    Ma.append(m)
                    aj = numpy.argmin(numpy.abs([a-m for a in lmodel['mass'][aii]]))
                    vals = [lmodel[P[0][0]][ai][i],lmodel[P[0][0]][aii][aj]]
                    f = interp1d([lmodel['age'][ai],lmodel['age'][aii]],vals)
                    Ge.append(f(numpy.log10(params['age'])))
            try:
                f = interp1d(Ge, Ma)
                params['mass'] = 10.**f(P[0][1])
            except:
                if kwargs.get('verbose',False) == True: print('\nFailed in mass + parameter determination\n')
                params['mass'] = numpy.nan

        if kwargs.get('debug',False) == True: print('\nMass unknown and Age known; determined mass to be {}'.format(params['mass']))

        Ma, Ge = [],[]


###################### KNOWN MASS AND AGE #####################################
# generate parameters as a function of mass interpolated between two closest age models
# evaluate parameters(mass)
###############################################################################
    if params['mass'] != 0. and params['age'] != 0. and \
        not numpy.isnan(params['age']) and not numpy.isnan(params['mass']):

        if kwargs.get('debug',False) == True: print(params)

        for i,age in enumerate(lmodel['age']):
            if numpy.nanmin(lmodel['mass'][i]) <= numpy.log10(params['mass']) \
                                              <= numpy.nanmax(lmodel['mass'][i]):
                Ag.append(age)
                f =interp1d(lmodel['mass'][i],lmodel['temperature'][i])
                Te.append(float(f(numpy.log10(params['mass']))))
                f = interp1d(lmodel['mass'][i],lmodel['luminosity'][i])
                Le.append(float(f(numpy.log10(params['mass']))))
                f = interp1d(lmodel['mass'][i],lmodel['gravity'][i])
                Ge.append(float(f(numpy.log10(params['mass']))))
                f = interp1d(lmodel['mass'][i],lmodel['radius'][i])
                Re.append(float(f(numpy.log10(params['mass']))))
  
        if params['temperature'] == 0.:
            try:
                f = interp1d(Ag, Te)
                params['temperature'] = 10.**float(f(numpy.log10(params['age'])))
            except: 
                params['temperature'] = numpy.nan
        if params['luminosity'] == 0.:
            try: 
                f = interp1d(Ag, Le)
                params['luminosity'] = float(f(numpy.log10(params['age'])).item(0))
            except: 
                params['luminosity'] = numpy.nan
        if params['gravity'] == 0.:
            try: 
                f = interp1d(Ag, Ge) 
                params['gravity'] = float(f(numpy.log10(params['age'])).item(0))      
            except: 
                params['gravity'] = numpy.nan
        if params['radius'] == 0.:
            try: 
                f = interp1d(Ag, Re)
                params['radius'] = 10.**float(f(numpy.log10(params['age'])))
            except: 
                params['radius'] = numpy.nan
  
        if kwargs.get('debug',False) == True: print('\nDetermined parameters: {}'.format(params))

        return params


# something failed    
    else:
#        print(params)
        for e in list(EVOLUTIONARY_MODEL_PARAMETERS.keys()):
            params[e] = numpy.nan
        if kwargs.get('verbose',False) == True: 
            print('\nParameter set is not covered by model {}\n'.format(model['name']))
        return params
