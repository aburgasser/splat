"""
.. note::
         These are the spectral modeling functions for SPLAT 
"""

import copy
import os
import sys
import urllib2
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import stats
from scipy.interpolate import griddata
import numpy
from astropy.io import ascii            # for reading in spreadsheet
import splat

SPECTRAL_MODEL_FOLDER = '/SpectralModels/'
MODEL_PARAMETER_NAMES = ['teff','logg','z','fsed','cld','kzz','slit']
MODEL_PARAMETERS = {'teff': 1000.0,'logg': 5.0,'z': 0.0,'fsed':'nc','cld':'nc','kzz':'eq','slit':0.5}
DEFINED_MODEL_SET = ['BTSettl2008','burrows06','morley12','morley14','saumon12','drift']
TMPFILENAME = 'splattmpfile'

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
    '''load up a model spectrum based on parameters'''

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
            raise NameError('\nCould not find '+kwargs['folder']+' locally or on SPLAT website\n\n')
        else:
            kwargs['folder'] = folder
            kwargs['local'] = False
            kwargs['online'] = True
    else:
        kwargs['folder'] = folder

# generate model filename
    kwargs['filename'] = kwargs['folder']+kwargs['set']+'_{:.0f}_{:.1f}_{:.1f}_{}_{}_{}_{:.1f}.txt'.\
        format(float(kwargs['teff']),float(kwargs['logg']),float(kwargs['z'])-0.001,kwargs['fsed'],kwargs['cld'],kwargs['kzz'],float(kwargs['slit']))

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
    '''Load up model parameters and check model inputs'''
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
            else:
                p = ascii.read(pfile)

# populate output parameter structure
    parameters = {'set': set, 'url': splat.SPLAT_URL}
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
