"""
.. note::
         These are the spectral modeling functions for SPLAT 
"""

import astropy
import bdevopar
import copy
from datetime import datetime
import os
import sys
import urllib2
#import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import stats
from scipy.integrate import trapz        # for numerical integration
from scipy.interpolate import griddata
import numpy
from astropy.io import ascii            # for reading in spreadsheet
from astropy.table import Table
import splat
import triangle

SPECTRAL_MODEL_FOLDER = '/reference/SpectralModels/'
MODEL_PARAMETER_NAMES = ['teff','logg','z','fsed','cld','kzz','slit']
MODEL_PARAMETERS = {'teff': 1000.0,'logg': 5.0,'z': 0.0,'fsed':'nc','cld':'nc','kzz':'eq','slit':0.5}
DEFINED_MODEL_SET = ['BTSettl2008','burrows06','morley12','morley14','saumon12','drift']
DEFINED_MODEL_NAME = ['BT-Settled (2008)','Burrows (2006)','Morley (2012)','Morley (2014)','Saumon (2012)','Drift (2008)']
TMPFILENAME = 'splattmpfile'
TEN_PARSEC = 443344480.     # ten parsecs in solar radii

# change the command prompt
sys.ps1 = 'splat model> '

#set the SPLAT PATH, either from set environment variable or from sys.path
#SPLAT_PATH = './'
#if os.environ.get('SPLAT_PATH') != None:
#    SPLAT_PATH = os.environ['SPLAT_PATH']
#else:
#    checkpath = ['splat' in r for r in sys.path]
#    if max(checkpath):
#        SPLAT_PATH = sys.path[checkpath.index(max(checkpath))]




def loadInterpolatedModel_NEW(*args,**kwargs):
# path to model
#    kwargs['path'] = kwargs.get('path',SPLAT_PATH+SPECTRAL_MODEL_FOLDER)
#    if not os.path.exists(kwargs['path']):
#        kwargs['remote'] = True
#        kwargs['path'] = SPLAT_URL+SPECTRAL_MODEL_FOLDER        
#    kwargs['set'] = kwargs.get('set','BTSettl2008')
#    kwargs['model'] = True
#    for ms in MODEL_PARAMETER_NAMES:
#        kwargs[ms] = kwargs.get(ms,MODEL_PARAMETERS[ms])

# first get model parameters
    pfile = 'parameters_new.txt'
#    parameters = loadModelParameters(**kwargs)

# insert a switch to go between local and online here

    print 'Running new version'
# read in parameters of available models
    folder = splat.checkLocal(splat.SPLAT_PATH+SPECTRAL_MODEL_FOLDER+kwargs['set']+'/')
    if folder=='':
        raise NameError('\n\nCould not locate spectral model folder {} locally\n'.format(SPECTRAL_MODEL_FOLDER+kwargs['set']+'/'))
    pfile = splat.checkLocal(folder+pfile)
    if pfile=='':
        raise NameError('\n\nCould not locate parameter list in folder {} locally\n'.format(folder))
    parameters = ascii.read(pfile)
#    numpy.genfromtxt(folder+pfile, comments='#', unpack=False, \
#        missing_values = ('NaN','nan'), filling_values = (numpy.nan)).transpose()

 
# check that given parameters are in range
    for ms in MODEL_PARAMETER_NAMES[0:3]:
        if (float(kwargs[ms]) < min(parameters[ms]) or float(kwargs[ms]) > max(parameters[ms])):
            raise NameError('\n\nInput value for {} = {} out of range for model set {}\n'.format(ms,kwargs[ms],kwargs['set']))
    for ms in MODEL_PARAMETER_NAMES[3:7]:
        if (kwargs[ms] not in parameters[ms]):
            raise NameError('\n\nInput value for {} = {} not one of the options for model set {}\n'.format(ms,kwargs[ms],kwargs['set']))

# now identify grid points around input parameters
# first set up a mask for digital parameters
    mask = numpy.ones(len(parameters[MODEL_PARAMETER_NAMES[0]]))
    for ms in MODEL_PARAMETER_NAMES[3:7]:
        m = [1 if a == kwargs[ms] else 0 for a in parameters[ms]]
        mask = mask*m
    
# identify grid points around input parameters
# 3x3 grid for teff, logg, z
# note interpolation and model ranges are separate
    dist = numpy.zeros(len(parameters[MODEL_PARAMETER_NAMES[0]]))
    for ms in MODEL_PARAMETER_NAMES[0:3]:
# first get step size
        ps = list(set(parameters[ms]))
        ps.sort()
        pps = numpy.abs(ps-ps[int(0.5*len(ps))])
        pps.sort()
        step = pps[1]
        dist = dist + ((kwargs[ms]-parameters[ms])/step)**2

# apply digital constraints
    ddist = dist/mask    

# find "closest" models
    mvals = {}
    for ms in MODEL_PARAMETER_NAMES[0:3]:
        mvals[ms] = [x for (y,x) in sorted(zip(ddist,parameters[ms]))]


# this is a guess as to how many models we need to get unique parameter values
    nmodels = 12
    mx,my,mz = numpy.meshgrid(mvals[MODEL_PARAMETER_NAMES[0]][0:nmodels],mvals[MODEL_PARAMETER_NAMES[1]][0:nmodels],mvals[MODEL_PARAMETER_NAMES[2]][0:nmodels])
    mkwargs = kwargs.copy()

    for i,w in enumerate(numpy.zeros(nmodels)):
        for ms in MODEL_PARAMETER_NAMES[0:3]:
            mkwargs[ms] = mvals[ms][i]
        mdl = loadModel(**mkwargs)
        if i == 0:
            mdls = numpy.log10(mdl.flux.value)
        else:
            mdls = numpy.column_stack((mdls,numpy.log10(mdl.flux.value)))

    print mdls[0,:]
    print (float(kwargs[MODEL_PARAMETER_NAMES[0]]),float(kwargs[MODEL_PARAMETER_NAMES[1]]),\
            float(kwargs[MODEL_PARAMETER_NAMES[2]]))
    print (mx.flatten(),my.flatten(),mz.flatten())
    for i in range(nmodels):
        print mvals['teff'][i],mvals['logg'][i],mvals['z'][i]
#
#
#            
## THIS NEXT PART IS BROKEN!
#
#
#
    mflx = numpy.zeros(len(mdl.wave))
    for i,w in enumerate(mflx):
#        print i, (mdls[i],) 
#        mflx[i] = 10.**(griddata((mx.flatten(),my.flatten(),mz.flatten()),val.flatten(),\
#            (float(kwargs['teff']),float(kwargs['logg']),float(kwargs['z'])),'linear'))

        m = mdls[i,:]
        mflx[i] = 10.**(griddata((mx.flatten(),my.flatten(),mz.flatten()),m.flatten(),\
            (float(kwargs[MODEL_PARAMETER_NAMES[0]]),float(kwargs[MODEL_PARAMETER_NAMES[1]]),\
            float(kwargs[MODEL_PARAMETER_NAMES[2]])),'linear'))
    
    return splat.Spectrum(wave=mdl.wave,flux=mflx*mdl.funit,**kwargs)




def loadInterpolatedModel(*args,**kwargs):
    '''
    .. not sure what force does
    :Purpose: Loads interpolated model spectrum based on parameters
    :param set: set of models to use; options include:

        - *'BTSettl2008'*: model set with effective temperature of 400 to 2900 K, surface gravity of 3.5 to 5.5 and metallicity of -3.0 to 0.5 
          from `Allard et al. (2012) <http://adsabs.harvard.edu/abs/2012RSPTA.370.2765A>`_
        - *'burrows06'*: model set with effective temperature of 700 to 2000 K, surface gravity of 4.5 to 5.5, metallicity of -0.5 to 0.5, 
          and sedimentation efficiency of either 0 or 100 from `Burrows et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...640.1063B>`_
        - *'morley12'*: model set with effective temperature of 400 to 1300 K, surface gravity of 4.0 to 5.5, metallicity of 0.0 
          and sedimentation efficiency of 2 to 5 from `Morley et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...756..172M>`_
        - *'morley14'*: model set with effective temperature of 200 to 450 K, surface gravity of 3.0 to 5.0, metallicity of 0.0 
          and sedimentation efficiency of 5 from `Morley et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...787...78M>`_
        - *'saumon12'*: model set with effective temperature of 400 to 1500 K, surface gravity of 3.0 to 5.5 and metallicity of 0.0 
          from `Saumon et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...750...74S>`_
        - *'drift'*: model set with effective temperature of 1700 to 3000 K, surface gravity of 5.0 to 5.5 and metallicity of -3.0 to 0.0 
          from `Witte et al. (2011) <http://adsabs.harvard.edu/abs/2011A%26A...529A..44W>`_
          
    :type set: optional, default = 'BTSettl2008'
    :param local: read in parameter file locally if True
    :type local: optional, default = True
    :param url: string of the url to the SPLAT website
    :type url: optional, default = 'http://pono.ucsd.edu/~adam/splat/'

    :Model Parameters: Below are parameters that define the model:

        - *teff*: effective temperature of the model
        - *logg*: surface gravity of the model
        - *z*: metallicity of the model
        - *fsed*: sedimentation efficiency of the model
        - *cld*: cloud shape function of the model
        - *kzz*: vertical eddy diffusion coefficient of the model
        - *slit*: slit weight of the model
    '''

# attempt to generalize models to extra dimensions
    mkwargs = kwargs.copy()
    mkwargs['force'] = True
    mkwargs['url'] = kwargs.get('url',splat.SPLAT_URL+'/Models/')
    mkwargs['set'] = kwargs.get('set','BTSettl2008')
    mkwargs['model'] = True
    mkwargs['local'] = kwargs.get('local',False)
    for ms in MODEL_PARAMETER_NAMES:
        mkwargs[ms] = kwargs.get(ms,MODEL_PARAMETERS[ms])

# first get model parameters
    parameters = loadModelParameters(**kwargs)
    
# check that given parameters are in range
    for ms in MODEL_PARAMETER_NAMES[0:3]:
        if (float(mkwargs[ms]) < parameters[ms][0] or float(mkwargs[ms]) > parameters[ms][1]):
            raise NameError('\n\nInput value for {} = {} out of range for model set {}\n'.format(ms,mkwargs[ms],mkwargs['set']))
    for ms in MODEL_PARAMETER_NAMES[3:6]:
        if (mkwargs[ms] not in parameters[ms]):
            raise NameError('\n\nInput value for {} = {} not one of the options for model set {}\n'.format(ms,mkwargs[ms],mkwargs['set']))

# identify grid points around input parameters
# 3x3 grid for teff, logg, z
# note interpolation and model ranges are separate
    mrng = []
    rng = []
    for ms in MODEL_PARAMETER_NAMES[0:3]:
        s = float(mkwargs[ms]) - float(mkwargs[ms])%float(parameters[ms][2])
        r = [max(float(parameters[ms][0]),s),min(s+float(parameters[ms][2]),float(parameters[ms][1]))]
        m = copy.deepcopy(r)
#        print s, r, s-float(kwargs[ms])
        if abs(s-float(mkwargs[ms])) < (1.e-3)*float(parameters[ms][2]):
            if float(kwargs[ms])%float(parameters[ms][2])-0.5*float(parameters[ms][2]) < 0.:
                m[1]=m[0]
                r[1] = r[0]+1.e-3*float(parameters[ms][2])
            else:
                m[0] = m[1]
                r[0] = r[1]-(1.-1.e-3)*float(parameters[ms][2])
#        print s, r, m, s-float(kwargs[ms])
        rng.append(r)
        mrng.append(m)
#        print s, r, m
    mx,my,mz = numpy.meshgrid(rng[0],rng[1],rng[2])
    mkwargs0 = mkwargs.copy()

# read in models
# note the complex path is to minimize model reads
    mkwargs['teff'] = mrng[0][0]
    mkwargs['logg'] = mrng[1][0]
    mkwargs['z'] = mrng[2][0]
    md111 = loadModel(**mkwargs)

    mkwargs['z'] = mrng[2][1]
    if (mrng[2][1] != mrng[2][0]):
        md112 = loadModel(**mkwargs)
    else:
        md112 = md111

    mkwargs['logg'] = mrng[1][1]
    mkwargs['z'] = mrng[2][0]
    if (mrng[1][1] != mrng[1][0]):
        md121 = loadModel(**mkwargs)
    else:
        md121 = md111

    mkwargs['z'] = mrng[2][1]
    if (mrng[2][1] != mrng[2][0]):
        md122 = loadModel(**mkwargs)
    else:
        md122 = md121

    mkwargs['teff'] = mrng[0][1]
    mkwargs['logg'] = mrng[1][0]
    mkwargs['z'] = mrng[2][0]
    if (mrng[0][1] != mrng[0][0]):
        md211 = loadModel(**mkwargs)
    else:
        md211 = md111

    mkwargs['z'] = mrng[2][1]
    if (mrng[2][1] != mrng[2][0]):
        md212 = loadModel(**mkwargs)
    else:
        md212 = md112

    mkwargs['logg'] = mrng[1][1]
    mkwargs['z'] = mrng[2][0]
    if (mrng[1][1] != mrng[1][0]):
        md221 = loadModel(**mkwargs)
    else:
        md221 = md211

    mkwargs['z'] = mrng[2][1]
    if (mrng[2][1] != mrng[2][0]):
        md222 = loadModel(**mkwargs)
    else:
        md222 = md221

    mflx = numpy.zeros(len(md111.wave))
    val = numpy.zeros([2,2,2])
    for i,w in enumerate(md111.wave):
        val = numpy.array([ \
            [[numpy.log10(md111.flux.value[i]),numpy.log10(md112.flux.value[i])], \
            [numpy.log10(md121.flux.value[i]),numpy.log10(md122.flux.value[i])]], \
            [[numpy.log10(md211.flux.value[i]),numpy.log10(md212.flux.value[i])], \
            [numpy.log10(md221.flux.value[i]),numpy.log10(md222.flux.value[i])]]])
        mflx[i] = 10.**(griddata((mx.flatten(),my.flatten(),mz.flatten()),val.flatten(),\
            (float(mkwargs0['teff']),float(mkwargs0['logg']),float(mkwargs0['z'])),'linear'))
    
    return splat.Spectrum(wave=md111.wave,flux=mflx*md111.funit,**kwargs)


def loadModel(*args, **kwargs):
    '''
    .. not sure what force does
    :Purpose: Load up a model spectrum based on parameters
    :param set: set of models to use; options include:

        - *'BTSettl2008'*: model set with effective temperature of 400 to 2900 K, surface gravity of 3.5 to 5.5 and metallicity of -3.0 to 0.5 
          from `Allard et al. (2012) <http://adsabs.harvard.edu/abs/2012RSPTA.370.2765A>`_
        - *'burrows06'*: model set with effective temperature of 700 to 2000 K, surface gravity of 4.5 to 5.5, metallicity of -0.5 to 0.5, 
          and sedimentation efficiency of either 0 or 100 from `Burrows et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...640.1063B>`_
        - *'morley12'*: model set with effective temperature of 400 to 1300 K, surface gravity of 4.0 to 5.5, metallicity of 0.0 
          and sedimentation efficiency of 2 to 5 from `Morley et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...756..172M>`_
        - *'morley14'*: model set with effective temperature of 200 to 450 K, surface gravity of 3.0 to 5.0, metallicity of 0.0 
          and sedimentation efficiency of 5 from `Morley et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...787...78M>`_
        - *'saumon12'*: model set with effective temperature of 400 to 1500 K, surface gravity of 3.0 to 5.5 and metallicity of 0.0 
          from `Saumon et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...750...74S>`_
        - *'drift'*: model set with effective temperature of 1700 to 3000 K, surface gravity of 5.0 to 5.5 and metallicity of -3.0 to 0.0 
          from `Witte et al. (2011) <http://adsabs.harvard.edu/abs/2011A%26A...529A..44W>`_
    
    :type set: optional, default = 'BTSettl2008'
    :param local: read in parameter file locally if True
    :type local: optional, default = True
    :param online: read in parameter file online if True
    :type online: optional, default = False
    :param folder: string of the folder name containing the model set
    :type folder: optional, default = ''
    :param filename: string of the filename of the desired model
    :type filename: optional
    :param url: string of the url to the SPLAT website
    :type url: optional, default = 'http://pono.ucsd.edu/~adam/splat/'

    :Model Parameters: Below are parameters that define the model:

        - *teff*: effective temperature of the model
        - *logg*: surface gravity of the model
        - *z*: metallicity of the model
        - *fsed*: sedimentation efficiency of the model
        - *cld*: cloud shape function of the model
        - *kzz*: vertical eddy diffusion coefficient of the model
        - *slit*: slit weight of the model
    '''

# path to model and set local/online
# by default assume models come from local splat directory
    local = kwargs.get('local',True)
    online = kwargs.get('online',not local and splat.checkOnline() != '')
    local = not online
    kwargs['local'] = local
    kwargs['online'] = online
    kwargs['folder'] = kwargs.get('folder','')
    kwargs['model'] = True
    kwargs['force'] = kwargs.get('force',False)
    url = kwargs.get('url',splat.SPLAT_URL)


# a filename has been passed - assume this file is a local file
# and check that the path is correct if its fully provided
# otherwise assume path is inside model set folder
    if (len(args) > 0):
        kwargs['filename'] = args[0]
        if not os.path.exists(kwargs['filename']):
            kwargs['filename'] = kwargs['folder']+os.path.basename(kwargs['filename'])
            if not os.path.exists(kwargs['filename']):
                raise NameError('\nCould not find model file {} or {}'.format(kwargs['filename'],kwargs['folder']+os.path.basename(kwargs['filename'])))
            else:
                return splat.Spectrum(**kwargs)
        else:
            return splat.Spectrum(**kwargs)

# set up the model set
    kwargs['set'] = kwargs.get('set','BTSettl2008')
    if kwargs.get('sed',False):
        kwargs['set'] = 'BTSettl2008'
    kwargs['folder'] = splat.SPLAT_PATH+SPECTRAL_MODEL_FOLDER+kwargs['set']+'/'
        
# preset defaults
    for ms in MODEL_PARAMETER_NAMES:
        kwargs[ms] = kwargs.get(ms,MODEL_PARAMETERS[ms])

# some special defaults
    if kwargs['set'] == 'morley12':
        if kwargs['fsed'] == 'nc':
            kwargs['fsed'] = 'f2'
    if kwargs['set'] == 'morley14':
        if kwargs['fsed'] == 'nc':
            kwargs['fsed'] = 'f5'
        if kwargs['cld'] == 'nc':
            kwargs['cld'] = 'f50'



# check that folder/set is present either locally or online
# if not present locally but present online, switch to this mode
# if not present at either raise error
    folder = splat.checkLocal(kwargs['folder'])
    if folder=='':
        folder = splat.checkOnline(kwargs['folder'])
        if folder=='':
            print '\nCould not find '+kwargs['folder']+' locally or on SPLAT website'
            print '\nAvailable model set options are:'
            for s in DEFINED_MODEL_SET:
                print '\t{}'.format(s)
            raise NameError()
        else:
            kwargs['folder'] = folder
            kwargs['local'] = False
            kwargs['online'] = True
    else:
        kwargs['folder'] = folder

# generate model filename
    kwargs['filename'] = kwargs['folder']+kwargs['set']+'_{:.0f}_{:.1f}_{:.1f}_{}_{}_{}_{:.1f}.txt'.\
        format(float(kwargs['teff']),float(kwargs['logg']),float(kwargs['z'])-0.001,kwargs['fsed'],kwargs['cld'],kwargs['kzz'],float(kwargs['slit']))
    kwargs['name'] = set
    if kwargs.get('sed',False):
        kwargs['filename'] = kwargs['folder']+kwargs['set']+'_{:.0f}_{:.1f}_{:.1f}_nc_nc_eq_sed.txt'.\
            format(float(kwargs['teff']),float(kwargs['logg']),float(kwargs['z'])-0.001)
        kwargs['name'] = set+' SED'

# get model parameters
#        parameters = loadModelParameters(**kwargs)
#        kwargs['path'] = kwargs.get('path',parameters['path'])
# check that given parameters are in range
#        for ms in MODEL_PARAMETER_NAMES[0:3]:
#            if (float(kwargs[ms]) < parameters[ms][0] or float(kwargs[ms]) > parameters[ms][1]):
#                raise NameError('\n\nInput value for {} = {} out of range for model set {}\n'.format(ms,kwargs[ms],kwargs['set']))
#        for ms in MODEL_PARAMETER_NAMES[3:6]:
#            if (kwargs[ms] not in parameters[ms]):
#                raise NameError('\n\nInput value for {} = {} not one of the options for model set {}\n'.format(ms,kwargs[ms],kwargs['set']))


# check if file is present; if so, read it in, otherwise go to interpolated
# locally:
    if kwargs['local']:
        file = splat.checkLocal(kwargs['filename'])
        if file=='':
            if kwargs['force']:
                raise NameError('\nCould not find '+kwargs['filename']+' locally\n\n')
            else:
                return loadInterpolatedModel(**kwargs)
#                kwargs['local']=False
#                kwargs['online']=True
        else:
            try:
                return splat.Spectrum(**kwargs)
            except:
                raise NameError('\nProblem reading in '+kwargs['filename']+' locally\n\n')

# online:
    if kwargs['online']:
        file = splat.checkOnline(kwargs['filename'])
        if file=='':
            if kwargs['force']:
                raise NameError('\nCould not find '+kwargs['filename']+' locally\n\n')
            else:
                return loadInterpolatedModel(**kwargs)
        else:
            try:
                ftype = kwargs['filename'].split('.')[-1]
                tmp = TMPFILENAME+'.'+ftype
                open(os.path.basename(tmp), 'wb').write(urllib2.urlopen(url+kwargs['filename']).read()) 
                kwargs['filename'] = os.path.basename(tmp)
                sp = splat.Spectrum(**kwargs)
                os.remove(os.path.basename(tmp))
                return sp
            except urllib2.URLError:
                raise NameError('\nProblem reading in '+kwargs['filename']+' from SPLAT website\n\n')




def loadModelParameters(**kwargs):
    '''
    .. is cld cloud shape function? kzz the vertical eddy diffusion coefficient?
    :Purpose: Load up model parameters and check model inputs. Parameters include 
                effective temperature, surface gravity (expressed as logg), metallicity, 
                and sedimentation efficiency (for cloudy models only).
    :param parameterFile: name of file containing parameters for spectral models
    :type parameterFile: optional, default = 'parameters.txt'
    :param set: set of models to use; options include:

        - *'BTSettl2008'*: model set with effective temperature of 400 to 2900 K, surface gravity of 3.5 to 5.5 and metallicity of -3.0 to 0.5 
          from `Allard et al. (2012) <http://adsabs.harvard.edu/abs/2012RSPTA.370.2765A>`_
        - *'burrows06'*: model set with effective temperature of 700 to 2000 K, surface gravity of 4.5 to 5.5, metallicity of -0.5 to 0.5, 
          and sedimentation efficiency of either 0 or 100 from `Burrows et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...640.1063B>`_
        - *'morley12'*: model set with effective temperature of 400 to 1300 K, surface gravity of 4.0 to 5.5, metallicity of 0.0 
          and sedimentation efficiency of 2 to 5 from `Morley et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...756..172M>`_
        - *'morley14'*: model set with effective temperature of 200 to 450 K, surface gravity of 3.0 to 5.0, metallicity of 0.0 
          and sedimentation efficiency of 5 from `Morley et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...787...78M>`_
        - *'saumon12'*: model set with effective temperature of 400 to 1500 K, surface gravity of 3.0 to 5.5 and metallicity of 0.0 
          from `Saumon et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...750...74S>`_
        - *'drift'*: model set with effective temperature of 1700 to 3000 K, surface gravity of 5.0 to 5.5 and metallicity of -3.0 to 0.0 
          from `Witte et al. (2011) <http://adsabs.harvard.edu/abs/2011A%26A...529A..44W>`_

    :type set: optional, default = 'BTSettl2008'
    :param online: read in parameter file online if True
    :type online: optional, default = False
    '''
# keyword parameters
    pfile = kwargs.get('parameterFile','parameters.txt')

# legitimate model set?
    if kwargs.get('set',False) not in DEFINED_MODEL_SET:
        raise NameError('\n\nInput model set {} not in defined set of models:\n{}\n'.format(set,DEFINED_MODEL_SET))
    

# read in parameter file - local and not local
    if kwargs.get('online',False):
        try:
            open(os.path.basename(TMPFILENAME), 'wb').write(urllib2.urlopen(splat.SPLAT_URL+SPECTRAL_MODEL_FOLDER+kwargs['set']+'/'+pfile).read())
            p = ascii.read(os.path.basename(TMPFILENAME))
            os.remove(os.path.basename(TMPFILENAME))
        except urllib2.URLError:
            print '\n\nCannot access online models for model set {}\n'.format(set)
#            local = True
    else:            
        if (os.path.exists(pfile) == False):
            pfile = splat.SPLAT_PATH+SPECTRAL_MODEL_FOLDER+kwargs['set']+'/'+os.path.basename(pfile)
            if (os.path.exists(pfile) == False):
                raise NameError('\nCould not find parameter file {}'.format(pfile))
        p = ascii.read(pfile)

# populate output parameter structure
    parameters = {'set': kwargs.get('set'), 'url': splat.SPLAT_URL}
    for ms in MODEL_PARAMETER_NAMES[0:3]:
        if ms in p.colnames:
            parameters[ms] = [float(x) for x in p[ms]]
        else:
            raise ValueError('\n\nModel set {} does not have defined parameter range for {}'.format(set,ms))
    for ms in MODEL_PARAMETER_NAMES[3:6]:
        if ms in p.colnames:
            parameters[ms] = str(p[ms][0]).split(",")
        else:
            raise ValueError('\n\nModel set {} does not have defined parameter list for {}'.format(set,ms))

    return parameters




#### the following codes are in progress

def modelFitGrid(spec, **kwargs):
    '''
    Model fitting code to grid of models
    '''
    print 'This function is not yet implemented'
    pass


def modelFitMCMC(spec, **kwargs):
    '''
    .. still need to add description of emodel
    :Purpose: Uses Markov chain Monte Carlo method to compare an object with models from a 
                given set. Returns the best estimate of the effective temperature, surface 
                gravity, and metallicity. Can also determine the radius of the object by 
                using these estimates. 
    :param spec: Spectrum class object, which should contain wave, flux and 
                  noise array elements.
    :param nsamples: number of Monte Carlo samples
    :type nsamples: optional, default = 1000
    :param initial_cut: the fraction of the initial steps to be discarded. (e.g., if 
                ``initial_cut = 0.2``, the first 20% of the samples are discarded.)
    :type initial_cut: optional, default = 0.1
    :param burn: the same as ``initial_cut``
    :type burn: optional, default = 0.1
    :param set: set of models to use; options include:

        - *'BTSettl2008'*: model set with effective temperature of 400 to 2900 K, surface gravity of 3.5 to 5.5 and metallicity of -3.0 to 0.5 
          from `Allard et al. (2012) <http://adsabs.harvard.edu/abs/2012RSPTA.370.2765A>`_
        - *'burrows06'*: model set with effective temperature of 700 to 2000 K, surface gravity of 4.5 to 5.5, metallicity of -0.5 to 0.5, 
          and sedimentation efficiency of either 0 or 100 from `Burrows et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...640.1063B>`_
        - *'morley12'*: model set with effective temperature of 400 to 1300 K, surface gravity of 4.0 to 5.5, metallicity of 0.0 
          and sedimentation efficiency of 2 to 5 from `Morley et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...756..172M>`_
        - *'morley14'*: model set with effective temperature of 200 to 450 K, surface gravity of 3.0 to 5.0, metallicity of 0.0 
          and sedimentation efficiency of 5 from `Morley et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...787...78M>`_
        - *'saumon12'*: model set with effective temperature of 400 to 1500 K, surface gravity of 3.0 to 5.5 and metallicity of 0.0 
          from `Saumon et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...750...74S>`_
        - *'drift'*: model set with effective temperature of 1700 to 3000 K, surface gravity of 5.0 to 5.5 and metallicity of -3.0 to 0.0 
          from `Witte et al. (2011) <http://adsabs.harvard.edu/abs/2011A%26A...529A..44W>`_
    
    :type set: optional, default = 'BTSettl2008'
    :param model: the same as ``set``
    :type model: optional, default = 'BTSettl2008'
    :param models: the same as ``set``
    :type models: optional, default = 'BTSettl2008'
    :param verbose: give lots of feedback
    :type verbose: optional, default = False
    :param mask_ranges: mask any flux value of ``spec`` by specifying the wavelength range.
                        Must be in microns.
    :type mask_ranges: optional, default = []
    :param mask_telluric: masks certain wavelengths to avoid effects from telluric absorption
    :type mask_telluric: optional, default = False
    :param mask_standard: masks wavelengths below 0.8 and above 2.35 microns
    :type mask_standard: optional, default = True
    :param mask: mask any flux value of ``spec``; has to be an array with length equal as ``spec`` with only 0 (unmask) or 1 (mask).
    :type mask: optional, default = [0, ..., 0] for len(sp1.wave)
    :param radius: calculates and returns radius of object if True
    :type radius: optional
    
    :param filename: filename or filename base for output
    :type filename: optional
    :param filebase: the same as ``filename``
    :type filebase: optional
    :param savestep: indicate when to save data output (e.g. ``savestep = 10`` will save the output every 10 samples)
    :type savestep: optional, default = ``nsamples/10``
    :param dataformat: output data format type
    :type dataformat: optional, default = 'ascii.csv'
    :param initial_guess: array including initial guess of the effective temperature, surface gravity and metallicity of ``spec``.
                            Can also set individual guesses of spectral parameters by using **initial_temperature** or **initial_teff**,
                            **initial_gravity** or **initial_logg**, and **initial_metallicity** or **initial_z**.
    :type initial_guess: optional, default = array of random numbers within allowed ranges
    :param ranges: array of arrays indicating ranges of the effective temperature, surface gravity and metallicity of the model set.
                            Can also set individual ranges of spectral parameters by using **temperature_range** or **teff_range**,
                            **gravity_range** or **logg_range**, and **metallicity_range** or **z_range**.
    :type ranges: optional, default = depends on model set
    :param step_sizes: an array specifying step sizes of spectral parameters. Can also set individual step sizes by using
                        **temperature_step** or **teff_step**, **gravity_step** or **logg_step**, and **metallicity_step** or **z_step**.
    :type step_sizes: optional, default = [50, 0.25, 0.1]
    :param nonmetallicity: if True, sets metallicity = 0
    :type nonmetallicity: optional, default = False
    :param addon: reads in prior calculation and starts from there. Allowed object types are tables, dictionaries and strings.
    :type addon: optional, default = False
    
    :Example:
    >>> import splat
    >>> sp = splat.getSpectrum(shortname='1047+2124')[0]        # T6.5 radio emitter
    >>> spt, spt_e = splat.classifyByStandard(sp,spt=['T2','T8'])
    >>> teff,teff_e = splat.typeToTeff(spt)
    >>> sp.fluxCalibrate('MKO J',splat.typeToMag(spt,'MKO J')[0],absolute=True)
    >>> table = splat.modelFitMCMC(sp, mask_standard=True, initial_guess=[teff, 5.3, 0.], zstep=0.1, nsamples=100, savestep=0, verbose=True)
        Trouble with model BTSettl2008 T=1031.61, logg=5.27, z=-0.02
        At cycle 0: fit = T=1031.61, logg=5.27, z=0.00 with chi2 = 35948.5
        Trouble with model BTSettl2008 T=1031.61, logg=5.27, z=-0.13
        At cycle 1: fit = T=1031.61, logg=5.27, z=0.00 with chi2 = 35948.5 
                                        .
                                        .
                                        .
                            # Skipped a few lines
                                        .
                                        .
                                        .
        Trouble with model BTSettl2008 T=973.89, logg=4.95, z=-0.17
        At cycle 99: fit = T=973.89, logg=4.95, z=0.00 with chi2 = 30569.6 
        <BLANKLINE>
        Number of steps = 170
        <BLANKLINE>
        Best Fit parameters:
        Lowest chi2 value = 29402.3750247 for 169.0 degrees of freedom
        Effective Temperature = 1031.608 (K)
        log Surface Gravity = 5.267 
        Metallicity = 0.000 
        Radius (relative to Sun) from surface fluxes = 0.103 
        <BLANKLINE>
        Median parameters:
        Effective Temperature = 1029.322 + 66.535 - 90.360 (K)
        log Surface Gravity = 5.108 + 0.338 - 0.473 
        Metallicity = 0.000 + 0.000 - 0.000 
        Radius (relative to Sun) from surface fluxes = 0.094 + 0.012 - 0.007 
        <BLANKLINE>
        <BLANKLINE>
        fit_J1047+2124_BTSettl2008
        Quantiles:
        [(0.16, 0.087231370556002871), (0.5, 0.09414839610875167), (0.84, 0.10562967101117798)]
        Quantiles:
        [(0.16, 4.6366512070621884), (0.5, 5.1077094570511488), (0.84, 5.4459108887603094)]
        Quantiles:
        [(0.16, 938.96254520460286), (0.5, 1029.3222563137401), (0.84, 1095.8574021575118)]
        <BLANKLINE>
        Total time elapsed = 0:01:46.340169
    >>> print table
             teff          logg      z       radius         chisqr   
        ------------- ------------- --- --------------- -------------
        1031.60790828 5.26704520744 0.0  0.103152256465 29402.3750247
        1031.60790828 5.26704520744 0.0  0.103152256465 29402.3750247
                  ...           ... ...             ...           ...   # Skipped a few lines
        938.962545205 5.43505121711 0.0  0.125429265207 43836.3720496
        938.962545205 5.43505121711 0.0  0.129294090544 47650.4267022
    '''

# code style note:
#   using the following equivalent terms: teff = temperature, logg = gravity, z = metallicity

# MCMC keywords
    timestart = datetime.now()
    nsample = kwargs.get('nsamples', 1000)
    burn = kwargs.get('initial_cut', 0.1)  # what fraction of the initial steps are to be discarded
    burn = kwargs.get('burn', burn)  # what fraction of the initial steps are to be discarded
    m_set = kwargs.get('set', 'BTSettl2008')
    m_set = kwargs.get('model', m_set)
    m_set = kwargs.get('models', m_set)
    verbose = kwargs.get('verbose', False)
# masking keywords
    mask_ranges = kwargs.get('mask_ranges',[])
    mask_telluric = kwargs.get('mask_telluric',False)
    mask_standard = kwargs.get('mask_standard',True)
    mask = kwargs.get('mask',numpy.zeros(len(spec.wave)))
# plotting and reporting keywords
    showRadius = kwargs.get('radius', spec.fscale == 'Absolute')
    filebase = kwargs.get('filebase', 'fit_'+spec.shortname+'_'+m_set)
    filebase = kwargs.get('filename',filebase)
    kwargs['filebase'] = filebase
    savestep = kwargs.get('savestep', nsample/10)
    dataformat = kwargs.get('dataformat','ascii.csv')
# evolutionary models    
    emodel = kwargs.get('evolutionary_model', 'baraffe')
    emodel = kwargs.get('emodel', emodel)
#    plot = kwargs.get('plot', False)
#    contour = kwargs.get('contour', False)
#    landscape = kwargs.get('landscape', False)
#    xstep = kwargs.get('xstep', 20)
#    ystep = kwargs.get('ystep', 20)

# set mask   
    if (mask_standard == True):
        mask_telluric = True
   
    if mask_telluric == True:
        mask_ranges.append([0.,0.65])        # meant to clear out short wavelengths
        mask_ranges.append([1.35,1.42])
        mask_ranges.append([1.8,1.92])
        mask_ranges.append([2.45,99.]) 
        
    if mask_standard == True:
        mask_ranges.append([0.,0.8])        # standard short cut
        mask_ranges.append([2.35,99.])      # standard long cut
        
    for ranges in mask_ranges:
        mask[numpy.where(((spec.wave.value >= ranges[0]) & (spec.wave.value <= ranges[1])))] = 1
    kwargs['mask_ranges'] = mask_ranges
    kwargs['mask_telluric'] = mask_telluric
    kwargs['mask_standard'] = mask_standard
    
# set the degrees of freedom    
    try:
        slitwidth = spec.slitpixelwidth
    except:
        slitwidth = 3.
    eff_dof = numpy.round((numpy.nansum(mask) / slitwidth) - 1.)

# TBD - LOAD IN ENTIRE MODEL SET

# set ranges for models - input or set by model itself
    rang = splat.loadModelParameters(set = m_set) # Range parameters can fall in
    ranges = kwargs.get('ranges', [rang['teff'][0:2], rang['logg'][0:2], rang['z'][0:2]])
    teff_range = kwargs.get('teff_range',ranges[0])
    teff_range = kwargs.get('temperature_range',teff_range)
    logg_range = kwargs.get('logg_range',ranges[1])
    logg_range = kwargs.get('gravity_range',logg_range)
    z_range = kwargs.get('z_range',ranges[2])
    z_range = kwargs.get('metallicity_range',z_range)

# set initial parameters

    param0_init = kwargs.get('initial_guess',[\
        numpy.random.uniform(teff_range[0],teff_range[1]),\
        numpy.random.uniform(logg_range[0],logg_range[1]),\
        numpy.random.uniform(z_range[0],z_range[1])])
    if len(param0_init) < 3:
        param0_init.append(0.0)
        
    t0 = kwargs.get('initial_temperature',param0_init[0])
    t0 = kwargs.get('initial_teff',t0)
    g0 = kwargs.get('initial_gravity',param0_init[1])
    g0 = kwargs.get('initial_logg',g0)
    z0 = kwargs.get('initial_metallicity',param0_init[2])
    z0 = kwargs.get('initial_z',z0)
    param0 = [t0,g0,z0]

    tstep = kwargs.get('teff_step',50)
    tstep = kwargs.get('temperature_step',tstep)
    gstep = kwargs.get('logg_step',0.25)
    gstep = kwargs.get('gravity_step',gstep)
    zstep = kwargs.get('z_step',0.1)
    zstep = kwargs.get('metallicity_step',zstep)

# catch for no metallicity input - assume not fitting this parameter
    param_step = kwargs.get('step_sizes',[tstep,gstep,zstep])
    if len(param_step) < 3:
        param_step.append[0.0]
        kwargs['nometallicity'] = True

    if kwargs.get('nometallicity',False):
        param_step[2] = 0.
        param0[2] = 0.0

    
# Check that initial guess is within range of models
    if not (teff_range[0] <= param0[0] <= teff_range[1] and \
        logg_range[0] <= param0[1] <= logg_range[1] and \
        z_range[0] <= param0[2] <= z_range[1]):
        sys.stderr.write("\nInitial guess T={}, logg = {} and [M/H] = {} is out of model range;" + \
            "defaulting to a random initial guess in range.".format(param0[0],param0[1],param0[2]))
        param0 = param0_init
        if param0[2] == 0.:
            param_step[2] = 0.


# read in prior calculation and start from there
    if kwargs.get('addon',False) != False:
        addflag = False
# a table is passed
        if isinstance(kwargs.get('addon'),Table):
            t = kwargs.get('addon')
            addflg = True
# a dictionary is passed
        elif isinstance(kwargs.get('addon'),dict):
            t = Table(kwargs.get('addon'))
            addflg = True
# a filename is passed
        elif isinstance(kwargs.get('addon'),str):
            try:
                p = ascii.read(kwargs.get('addon'))
            except:
                print '\nCould not read in parameter file {}'.format(kwargs.get('addon'))


# initial fit cycle
    try:
        model = splat.loadModel(teff = param0[0], logg = param0[1], z = param0[2], set = m_set)
    except:
        raise ValueError('\nInitial {} model with T = {}, logg = {} and [M/H] = {} did not work; aborting.'.format(m_set,param0[0],param0[1],param0[2]))

    chisqr0,alpha0 = splat.compareSpectra(spec, model, mask_ranges=mask_ranges)
    chisqrs = [chisqr0]    
    params = [param0]
    radii = [TEN_PARSEC*numpy.sqrt(alpha0)]
    for i in range(nsample):
        for j in range(len(param0)):
            if param_step[j] > 0.:          # efficient consideration - if statement or just run a model?
                param1 = copy.deepcopy(param0)
                param1[j] = numpy.random.normal(param1[j],param_step[j])
                try:            
                    model = splat.loadModel(teff = param1[0], logg = param1[1],z = param1[2], set = m_set)
                    chisqr1,alpha1 = splat.compareSpectra(spec, model ,mask_ranges=mask_ranges)  

# Probability that it will jump to this new point; determines if step will be taken
                    h = 1. - stats.f.cdf(chisqr1/chisqr0, eff_dof, eff_dof)
#                    print chisqr1, chisqr0, eff_dof, h
                    if numpy.random.uniform(0,1) < h:
                        param0 = copy.deepcopy(param1)
                        chisqr0 = copy.deepcopy(chisqr1)
                        alpha0 = copy.deepcopy(alpha1)
                
# update list of parameters, chi^2 and radii
                    params.append(param0)
                    chisqrs.append(chisqr0)
                    radii.append(TEN_PARSEC*numpy.sqrt(alpha0))
                    
                except:
                    if verbose:
                        print 'Trouble with model {} T={:.2f}, logg={:.2f}, z={:.2f}'.format(m_set,param1[0],param1[1],param1[2])
                    continue

        if verbose:
            print 'At cycle {}: fit = T={:.2f}, logg={:.2f}, z={:.2f} with chi2 = {:.1f}'.format(i,param0[0],param0[1],param0[2],chisqr0)

# save results iteratively
        if i*savestep != 0:
            if i%savestep == 0:
                t = Table(zip(*params[::-1]),names=['teff','logg','z'])
                if param_step[2] == 0.:
                    del t['z']
                if showRadius:
                    t['radius'] = radii
                t['chisqr'] = chisqrs
                t.write(filebase+'rawdata.dat',format=dataformat)
                reportModelFitResults(spec,t,iterative=True,model_set=m_set,**kwargs)

# Final results
    t = Table(zip(*params[::-1]),names=['teff','logg','z'])
    if param_step[2] == 0. or kwargs.get('nometallicity',False):
        del t['z']
    if showRadius:
        t['radius'] = radii
    t['chisqr'] = chisqrs
# cut first x% of parameters
    s = Table(t[burn*len(t):])

# save data
    s.write(filebase+'rawdata.dat',format=dataformat)
    
    reportModelFitResults(spec,s,iterative=False,model_set=m_set,**kwargs)
    if verbose:
        print '\nTotal time elapsed = {}'.format(datetime.now()-timestart)
    return s


def calcLuminosity(sp, mdl=False, absmags=False, **kwargs):
    '''
    Calculate luminosity from photometry and stitching models
    
    '''

    spec_filters = ['SDSS Z','2MASS J','2MASS H','2MASS KS','MKO J','MKO H','MKO K']
    sed_filters = ['SDSS R','SDSS I','WISE W1','WISE W2','WISE W3','WISE W4','IRAC CH1','IRAC CH2','IRAC CH3','IRAC CH4']
    
    if ~isinstance(absmags,dict):
        raise ValueError('\nAbsolute magnitudes should be a dictionary whose keys are one of the following filters:\n{}'.format(spec_filters+sed_filters))

# read in a model if one is not provided based on classification and temperature
    if mdl == False or 'SED' not in mdl.name:
        spt,spt_unc = splat.classifyByIndex(sp)
        teff,unc = splat.typeToTeff(spt)
        mdl = splat.loadModel(teff=teff,logg=5.0,sed=True)

# prep arrays
    flux = []
    flux_unc = []
    flux_wave = []
    
# steps:
# scale spectrum to absolute magnitude if necessary and integrate flux, varying noise and including variance in abs mag factor
    spcopy = sp
    if spcopy.fscale != 'Absolute':
        scale = []
        scale_unc = []
        for k in absmags.keys():
            if k.upper() in spec_filters:
                m = splat.filterMag(spcopy,k)
                scale.extend(10.**(0.4*(m-absmags[k][0])))
# note: need to add in spectral flux uncertainty as well
                scale_unc.extend(numpy.log(10.)*0.4*absmags[k][1]*scale[-1])
        if len(scale) == 0:
            raise ValueError('\nNo absolute magnitudes provided to scale spectrum; you specified:\n{}'.format(absmags.keys()))
        scl,scl_e = splat.weightedMeanVar(scale,scale_unc,uncertainty=True)
        spcopy.scale(numpy.mean(scl))
        spcopy.fscale = 'Absolute'

# integrate data
# NEED TO INSERT UNCERTAINTY HERE
    flux.extend(trapz(spcopy.flux,spcopy.wave))
    flux_unc.extend(0.)
    flux_wave.extend([numpy.nanmin(spcopy.wave),numpy.nanmax(spcopy.wave)])

# scale segments of models scaled to WISE or IRAC bands if available, include variance in abs mag factor
# PROBLEM: WHAT IF SPECTRAL PIECES OVERLAP?
    for k in absmags.keys():
        if k.upper() in sed_filters:
            filterdat = splat.filterProperties(k.upper())
            mdl.fluxCalibrate(k,absmags[k][0])
            w = numpy.where(mdl.wave.value >= filterdat['lambda_min'] and mdl.wave.value <= filterdat['lambda_max'])
            flux.extend(trapz(mdl.flux[w],mdl.wave[w]))
            flux_unc.extend(2.5*numpy.log(10.)*absmags[k][1]*flux[-1])
            flux_wave.extend([filterdat['lambda_min'],filterdat['lambda_max']])

# match model between these scaled pieces and out to ends and integrate, include variance in abs mag factor(s)
# report log luminosity in solar units and uncertainty
# optional report the various pieces and percentages of whole ()
#
# absmags is a dictionary whose keys are filter names and whose elements are 2-element lists of value and uncertainty        
    
def reportModelFitResults(spec,t,*arg,**kwargs):
    '''
    Reports the result of model fitting parameters
    Input is an astropy Table with columns containing parameters fit, and one column for chi-square values ('chisqr')
    Produces triangle plot, best fit model, statistics of parameters
    and saves raw data if iterative = True
    '''    

    evolFlag = kwargs.get('evol',True)
    emodel = kwargs.get('emodel','Baraffe')
    statsFlag = kwargs.get('stats',True)
    triangleFlag = kwargs.get('triangle',True)
    bestfitFlag = kwargs.get('bestfit',True)
    summaryFlag = kwargs.get('summary',True)
    weights = kwargs.get('weight',None)
    filebase = kwargs.get('filebase','modelfit_results')
    statcolumn = kwargs.get('stat','chisqr')
    mset = kwargs.get('model_set','')
    mset = kwargs.get('mset',mset)
    mask_ranges = kwargs.get('mask_ranges',[])
    sigma = kwargs.get('sigma',1.)

# map some common column names to full descriptive texts
    plotname_assoc = {\
        'teff': r'T$_{eff}$',\
        'logg': r'log g',\
        'z': r'[M/H]',\
        'mass': r'M/M$_{\odot}$',\
        'age': r'$\tau$',\
        'lbol': r'log L$_{bol}$/L$_{\odot}$',\
        'radius': r'R/R$_{\odot}$',\
        'radius_evol': r'R/R$_{\odot}$'}

    format_assoc = {\
        'teff': '.0f',\
        'logg': '.2f',\
        'z': '.2f',\
        'mass': '.3f',\
        'age': '.1f',\
        'lbol': '.2f',\
        'radius': '.3f',\
        'radius_evol': '.3f'}

    descrip_assoc = {\
        'teff': 'Effective Temperature',\
        'logg': 'log Surface Gravity',\
        'z': 'Metallicity',\
        'mass': 'Mass',\
        'age': 'Age',\
        'lbol': 'log Luminosity (relative to Sun)',\
        'radius': 'Radius (relative to Sun) from surface fluxes',\
        'radius_evol': 'Radius (relative to Sun) from evolutionary models'}

    unit_assoc = {\
        'teff': 'K',\
#        'logg': r'cm/s$^2$',\
        'age': 'Gyr'}

    if kwargs.get('iterative',False):
        statsFlag = True
        evolFlag = True
        triangleFlag = False
        bestfitFlag = False
        summaryFlag = False

# check that table has the correct properties
    if isinstance(t,astropy.table.Table) == False:
        raise ValueError('\nInput is not an astropy table')
    if len(t.columns) < 2:
        raise ValueError('\nNeed at least two columns in input table')
    if statcolumn not in t.colnames:
        raise ValueError('\n{} column must be present in input table'.format(statcolumn))

    parameters = t.colnames
    parameters.remove(statcolumn)
 
# get the evolutionary model parameters
# turned off for now
    evolFlag = False
    if evolFlag:
        if 'teff' not in t.colnames or 'logg' not in t.colnames:
            print '\nCannot compare to best fit without teff and logg parameters'

        else:
            values=bdevopar.Parameters(emodel, teff=t['teff'], grav=t['logg'])
            t['age'] = values['age']
            t['lbol'] = values['luminosity']
            t['mass'] = values['mass']
            t['radius_evol'] = values['radius']
            parameters = t.colnames
            parameters.remove(statcolumn)
    
   
# calculate statistics
    if statsFlag:
        if weights == True:
            weights = numpy.exp(0.5*(numpy.nanmin(t[statcolumn])-t[statcolumn]))

        print '\nNumber of steps = {}'.format(len(t))    
        
        print '\nBest Fit parameters:'
        print 'Lowest chi2 value = {} for {} degrees of freedom'.format(numpy.nanmin(t[statcolumn]),spec.dof)
        for p in parameters:
            sort = [x for (y,x) in sorted(zip(t[statcolumn],t[p]))]
            name = p
            if p in descrip_assoc.keys():
                name = descrip_assoc[p]
            unit = ''
            if p in unit_assoc.keys():
                unit = '('+unit_assoc[p]+')'
            print '{} = {:.3f} {}'.format(name,sort[0],unit)

        print '\nMedian parameters:'
        for p in parameters:
            sm, mn, sp = distributionStats(t[p],sigma=sigma,weights=weights)      # +/- 1 sigma
            name = p
            if p in descrip_assoc.keys():
                name = descrip_assoc[p]
            unit = ''
            if p in unit_assoc.keys():
                unit = '('+unit_assoc[p]+')'
            print '{} = {:.3f} + {:.3f} - {:.3f} {}'.format(name,mn,sp-mn,mn-sm,unit)
        print '\n'

        
# best fit model
    if bestfitFlag and mset in DEFINED_MODEL_SET:
# check to make sure at least teff & logg are present
        if 'teff' not in t.colnames or 'logg' not in t.colnames:
            print '\nCannot compare to best fit without teff and logg parameters'

        else:
            t.sort(statcolumn)
            margs = {'set': mset, 'teff': t['teff'][0], 'logg': t['logg'][0]}
            legend = [spec.name,'{} T = {:.0f}, logg =  {:.2f}'.format(mset,margs['teff'],margs['logg']),r'$\chi^2$ = {:.0f}, DOF = {:.0f}'.format(t[statcolumn][0],spec.dof)]
            if 'z' in t.colnames:
                margs['z'] = t['z'][0]
                legend[1]+=', z = {:.2f}'.format(margs['z'])
            model = splat.loadModel(**margs)
            chisqr,alpha = splat.compareSpectra(spec, model ,mask_ranges=mask_ranges)
            model.scale(alpha)

            w = numpy.where(numpy.logical_and(spec.wave.value > 0.9,spec.wave.value < 2.35))
            diff = spec-model
            print filebase
            splat.plotSpectrum(spec,model,diff,uncertainty=True,telluric=True,colors=['k','r','b'], \
                legend=legend,filename=filebase+'bestfit.eps',\
                yrange=[1.1*numpy.nanmin(diff.flux.value[w]),1.25*numpy.nanmax([spec.flux.value[w],model.flux.value[w]])])  
 
# triangle plot of parameters
    if triangleFlag:
        y=[]
        labels = []
        fmt = []
        for p in parameters:
            if min(numpy.isfinite(t[p])) == True and numpy.nanstd(t[p]) != 0.:       # patch to address arrays with NaNs in them
                y.append(t[p])
                if p in plotname_assoc.keys():
                    labels.append(plotname_assoc[p])
                else:
                    labels.append(p)
                if p in unit_assoc.keys():
                    labels[-1] = labels[-1]+' ('+unit_assoc[p]+')'
                if p in format_assoc.keys():
                    fmt.append(format_assoc[p])
                else:
                    fmt.append('.2f')
#        print labels        
        print labels, fmt
        fig = triangle.corner(zip(*y[::-1]), labels=list(reversed(labels)), show_titles=True, quantiles=[0.16,0.5,0.84],cmap=cm.Oranges,title_fmt=list(reversed(fmt)),plot_contours=True)
        fig.savefig(filebase+'parameters.eps')
           
# plain language summary
    if summaryFlag:
        pass
            

def distributionStats(x, q=[0.16,0.5,0.84], weights=None, sigma=None, **kwargs):
    '''
    Find key values along distributions based on quantile steps
    This code is derived almost entirely from triangle.py
    '''

# clean data of nans
    xd = x[~numpy.isnan(x)]

    if q is None and sigma is None:
        sigma = 1.
        
    if sigma is not None:
        q = [stats.norm.cdf(-sigma),0.5,stats.norm.cdf(sigma)]
        
    if weights is None:
        return numpy.percentile(xd, [100. * qi for qi in q])
    else:
        wt = weights[~numpy.isnan(x)]
        idx = numpy.argsort(xd)
        xsorted = xd[idx]
        cdf = numpy.add.accumulate(wt[idx])
        cdf /= cdf[-1]
        return numpy.interp(q, cdf, xsorted).tolist()



if __name__ == '__main__':
    basefolder = '/Users/adam/projects/splat/exercises/ex9/'
    sp = splat.getSpectrum(shortname='1047+2124')[0]        # T6.5 radio emitter
    spt,spt_e = splat.classifyByStandard(sp,spt=['T2','T8'])
    teff,teff_e = splat.typeToTeff(spt)
    sp.fluxCalibrate('MKO J',splat.typeToMag(spt,'MKO J')[0],absolute=True)
    table = modelFitMCMC(sp, mask_standard=True, initial_guess=[teff, 5.3, 0.], zstep=0.1, nsamples=100,savestep=0,filebase=basefolder+'fit1047',verbose=True)

