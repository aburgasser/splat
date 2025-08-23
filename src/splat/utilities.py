# -*- coding: utf-8 -*-
from __future__ import print_function

"""
.. note::
         These are the utility functions for SPLAT 
"""

# imports: internal
import base64
import copy
import os
import pandas
import re
import requests
import string
import sys

# imports - external
import astropy
from astropy.coordinates import Angle,SkyCoord,EarthLocation      # coordinate conversion
from astropy import units as u            # standard units
from astropy.time import Time            # standard units
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects
import numpy
from scipy import stats
from scipy.interpolate import interp1d,InterpolatedUnivariateSpline
from scipy.integrate import trapezoid as trapz


# code constants
import splat
from splat.initialize import *

# Python 2->3 fix for input
try: input=raw_input
except NameError: pass

# change the command prompt
#sys.ps1 = 'splat util> '


#####################################################
###########   SIMPLE HELPER FUNCTIONS   #############
#####################################################


def isNumber(s):
    '''
    :Purpose: Checks if something is a number.

    :param s: object to be checked
    :type s: required

    :Output: True or False

    :Example:
    >>> import splat
    >>> print splat.isNumber(3)
        True
    >>> print splat.isNumber('hello')
        False
    '''
    s1 = copy.deepcopy(s)
    if isinstance(s1,bool): return False
    if isinstance(s1,u.quantity.Quantity): s1 = s1.value
    if isinstance(s1,float): return (True and not numpy.isnan(s1))
    if isinstance(s1,int): return (True and not numpy.isnan(s1))
    try:
        s1 = float(s1)
        return (True and not numpy.isnan(s1))
    except:
        return False

def isUnit(s):
    '''
    :Purpose: 
        Checks if something is an astropy unit quantity; written in response to the 
        many ways that astropy now codes unit quantities

    :Required Inputs: 
        :param s: quantity to be checked

    :Optional Inputs: 
        None

    :Output: 
        True or False

    :Example:
    >>> import splat
    >>> import astropy.units as u
    >>> print splat.isUnit(3)
        False
    >>> print splat.isUnit(3.*u.s)
        True
    >>> print splat.isUnit(3.*u.s/u.s)
        True
    >>> print splat.isUnit((3.*u.s/u.s).value)
        False
    '''
    return isinstance(s,u.quantity.Quantity) or \
        isinstance(s,u.core.Unit) or \
        isinstance(s,u.core.CompositeUnit) or \
        isinstance(s,u.core.IrreducibleUnit) or \
        isinstance(s,u.core.NamedUnit) or \
        isinstance(s,u.core.PrefixUnit)


def numberList(numstr,sort=False):
    '''
    :Purpose: 

        Convert a string listing of numbers into an array of numbers

    :Required Input:

        :param **numstr**: string indicating number list, e.g., '45,50-67,69,72-90'

    :Optional Input:

        :param **sort**: set to True to sort output list (default = False)

    :Output:

        list of integers specified by string

    :Example:
        >>> import splat
        >>> a = splat.numberList('45,50-67,69,72-90')
        >>> print(a)
            [45, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 69, 
            72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90]
    '''
# check inputs
    if not isinstance(numstr,str): raise ValueError('\nInput to numberList {} must be a string'.format(numstr))

    numlist = []
    tmp1 = numstr.replace(' ','')
    tmp2 = tmp1.split(',')
    for a in tmp2:
        tmp3 = a.split(';')
        for b in tmp3:
            tmp4 = b.split('-')
            if len(tmp4) > 1:
                numlist.extend(list(range(int(tmp4[0]),int(tmp4[1])+1)))
            else:
                numlist.append(int(tmp4[0]))
    
    if sort==True: numlist = sorted(numlist)
    return numlist

def padWhereArray(w,mx):
    '''
    Purpose:

        Pads the output of a numpy.where array to select (if available) one more index spot beyond limits

    '''
    if w[0][0] > 0: w = (numpy.insert(w[0],0,w[0][0]-1),)
    if w[0][-1] < mx: w = (numpy.append(w[0],w[0][-1]+1),)
    return w


def readDictFromFile(file,delim='\t',missing_value=None,data_type=[str],verbose=False,**kwargs):
    '''
    :Purpose:

        Reads a simple text file into a series of key: value pairs and placed into a dictionary;
        allows for assignment of variables to arrays

    :Required Inputs:

        :param: file: string containing full path to file to be read in; this file should be an ascii file with simple delimiters

    :Optional Inputs:

        :param: delim: delimiter to separate keys from values
        :param: value_delim: delimiter to separate values; if not provided, defaults to ``delim``
        :param: data_type: single or list of data type to apply to input data; must be str, int, float or complex
        :param: missing_value: variable to replace missing values (keys without data)
        :param: verbose: set to True to provide verbose feedback

    :Outputs:

        A dictionary of input file parameters

    :Example:

        Assume you have a data file of format:

            this 5
            that 5,6,7,8
            other

        >>> import splat
        >>> readDictFromFile('input.txt',delim=' ',value_delim=',',data_type=[int,float])
            {'this': 5, 'that': [6.0, 7.0, 8.0], 'other': None}
    '''
    list_delim = kwargs.get('list_delim',delim)
    list_delim = kwargs.get('value_delim',list_delim)
    
    if os.path.exists(file) == False:
        raise ValueError('\nFile {} cannot be found'.format(file))
    try:
        with open(file) as f: dat = f.read()
        dat = dat.split('\n')
    except:
        raise ValueError('\nUnable to read in file {} as simple ascii file'.format(file))
    if len(dat) == 0:
        if verbose == True: print('\nNo data found in file {}'.format(file))
        return {}
    if len(dat[0].split(delim)) < 2:
        if verbose == True: print('\nWarning: delimiter {} not found in first line of file {}'.format(file))

# data types
    try:
        dtype = list(data_type)
    except:
        dtype = copy.deepcopy(data_type)
    if not isinstance(dtype,list): dtype = [dtype]
#        if verbose == True: print('\nWarning: could not intepret data type input {}, converting all to strings'.format(data_type))
    while len(dtype) < len(dat): dtype.append(dtype[-1])

# separate and convert        
    output = {}
    for i,line in enumerate(dat):
        if line != '':
            sp = line.split(delim)
            ky = sp[0]
            if len(sp) > 1:
                val = sp[1:]
                if list_delim != delim: val = sp[1].split(list_delim)
                d = dtype[i]
                if d not in [str,int,float,complex]: d = str
                cval = []
                for v in val:
                    try: cval.append(d(v))
                    except: pass
                if len(cval) == 1: cval = cval[0]
            else: cval = missing_value
            output[ky] = cval

    return output
        

def writeDictToFile(data,file,delim='\t',verbose=False,**kwargs):
    '''
    :Purpose:

        Writes the contents of a dictionary to a simple ascii file into a series of key value pairs;
        allows for writing of both individual variables and lists (but not nested dictionaries)

    :Required Inputs:

        :param: data: dictionary to be written out; cannot be a nested dictionary but can contain lists
        :param: file: string containing full path to file to be written

    :Optional Inputs:

        :param: delim: delimiter to separate keys from values
        :param: value_delim: delimiter to separate values; if not provided, defaults to ``delim``
        :param: verbose: set to True to provide verbose feedback

    :Outputs:

        An output file

    :Example:

        >>> import splat
        >>> d = {'this': 5., 'that': [4,6,8], 'other': 'something else'}
        >>> writeDictToFile(d,'/Users/adam//Desktop/temp2.txt',delim='\t',value_delim=',')
            True

        Contents of file will be:

            this    5.0
            that    4,6,8
            other   something else

    '''
    value_delim = kwargs.get('value_delim',delim)
    value_delim = kwargs.get('list_delim',value_delim)

    if isinstance(data,dict) == False: 
        raise ValueError('\nInput data is not a dictionary'.format(file))
    try:
        f = open(file,'w')
    except:
        raise ValueError('\nCould not open file {} for writing'.format(file))

    for k in list(data.keys()):
        line = '{}{}'.format(k,delim)
        val = data[k]
        if isinstance(val,str): val = [val]
        try:
            val = list(val)
        except:
            val = [val]
        line = line+'{}'.format(val[0])
        if len(val) > 1:
            for v in val[1:]: line = line+'{}{}'.format(value_delim,v)
        f.write(line+'\n')
    f.close()

    return True

def directoryTree(folder,verbose=ERROR_CHECKING):
    '''
    :Purpose:

        Finds the lowest level directories within a given folder and returns the full paths for these

    :Required Inputs:

        :param: folder: directory to search

    :Optional Inputs:

        :param: verbose: set to True to provide verbose feedback

    :Outputs:

        A list of directory paths

    :Example:

        >>> import splat
        >>> directoryTree(splat.LIBRARY_PUBLIC_FOLDER)
            ['/Users/adam/projects/splat/code/splat//resources/Data/Public/MAGE/',
             '/Users/adam/projects/splat/code/splat//resources/Data/Public/SPEX-PRISM/',
             '/Users/adam/projects/splat/code/splat//resources/Data/Public/LRIS-RED/']
    '''
    paths = []
    if os.path.exists(folder)==False:
        if verbose==True: print('Warning: folder {} cannot be found'.format(folder))
    else:
        for p,d,r in os.walk(folder):
            if not d: paths.append(os.path.join(p,''))
    return paths
    


#####################################################
################   VARIOUS CHECKS   #################
#####################################################


def checkFile(filename,**kwargs):
    '''
    :Purpose: Checks if a spectrum file exists in the SPLAT's library.
    :param filename: A string containing the spectrum's filename.
    :Example:
       >>> import splat
       >>> spectrum1 = 'spex_prism_1315+2334_110404.fits'
       >>> print spl.checkFile(spectrum1)
       True
       >>> spectrum2 = 'fake_name.fits'
       >>> print spl.checkFile(spectrum2)
       False
    '''
    url = kwargs.get('url',SPLAT_URL)+DATA_FOLDER
    return requests.get(url+filename).status_code == requests.codes.ok



def checkLocal(inputfile):
    '''
    :Purpose: Checks if a file is present locally or within the SPLAT
                code directory
    :Example:
       >>> import splat
       >>> spl.checkLocal('spl.py')
       True  # found the code
       >>> spl.checkLocal('parameters.txt')
       False  # can't find this file
       >>> spl.checkLocal('SpectralModels/BTSettl08/parameters.txt')
       True  # found it
    '''
    if not os.path.exists(os.path.normpath(inputfile)):
        if not os.path.exists(os.path.normpath(inputfile)):
            return ''
        else:
            return inputfile
    else:
        return inputfile



def checkOnline(*args):
    '''
    :Purpose: Checks if SPLAT's URL is accessible from your machine--
                that is, checks if you and the host are online. Alternately
                checks if a given filename is present locally or online
    :Example:
       >>> import splat
       >>> spl.checkOnline()
       True  # SPLAT's URL was detected.
       >>> spl.checkOnline()
       False # SPLAT's URL was not detected.
       >>> spl.checkOnline('SpectralModels/BTSettl08/parameters.txt')
       '' # Could not find this online file.
    '''
    output = False

    if len(args) != 0:
        if 'http://' in args[0]:
            try:
                if requests.get(args[0]).status_code == requests.codes.ok:
                    output = args[0]
            except:
                pass
        else:
            try:
                if requests.get(SPLAT_URL+args[0]).status_code == requests.codes.ok:
                    output = SPLAT_URL+args[0]
            except:
                pass
    else:
        try:
            output = requests.get(SPLAT_URL).status_code == requests.codes.ok
        except:
            pass
    return output



def checkOnlineFile(*args):
    '''
    :Purpose: Checks if SPLAT's URL is accessible from your machine--
                that is, checks if you and the host are online. Alternately
                checks if a given filename is present locally or online
    :Example:
       >>> import splat
       >>> spl.checkOnlineFile('SpectralModels/BTSettl08/parameters.txt')
       '' # Could not find this online file.
       >>> spl.checkOnlineFile()
       '' # SPLAT's URL was not detected; you are not online.
    '''
    if (len(args) != 0):
        if 'http://' in args[0]:
            if requests.get(args[0]).status_code == requests.codes.ok:
                return args[0]
            return ''
        else:
            if requests.get(SPLAT_URL+args[0]).status_code == requests.codes.ok:
                return SPLAT_URL+args[0]
            return ''
    else:
        return requests.get(SPLAT_URL).status_code == requests.codes.ok


def checkDict(ref,refdict,altref='altname',replace=[],verbose=False):
    '''
    Purpose: 
        General usage program to check if a key is present in a dictionary, with the option to look through alternate names

    Required Inputs:
        :param ref: A string containing the reference for lumiosity/SpT relation, should be among the keys and alternate names in refdict
        :param refdict: dictionary containing empirical relation information

    Optional Inputs:
        None

    Output:
        A string containing SPLAT's default name for a given reference set, or False if that reference is not present

    Example:

    >>> import splat
    >>> print(splat.checkDict('filippazzo',splat.SPT_LBOL_RELATIONS))
        filippazzo2015
    >>> print(splat.checkDict('burgasser',splat.SPT_BC_RELATIONS))
        False
    '''
    output = False
    refc = copy.deepcopy(ref)

# check reference    
    if not isinstance(refc,str):
        return output
    if len(replace) > 0:
        for rep in replace:
            if isinstance(rep,list) == True and len(rep) > 0: refc = refc.replace(rep[0],rep[1])
    for k in list(refdict.keys()):
        if refc.lower()==k.lower(): output = k
        if altref in list(refdict[k].keys()):
            if refc.lower() in [x.lower() for x in list(refdict[k][altref])]: output = k
    if output == False:
        if verbose: print('\nCould not find item {} in input dictionary; try: {}'.format(ref,list(refdict.keys())))

    return output


def checkEmpiricalRelation(ref,refdict,verbose=False):
    '''
    Purpose: 
        General checking program for empirical relation dictionaries 

    Required Inputs:
        :param ref: A string containing the reference for lumiosity/SpT relation, should be among the keys and alternate names in refdict
        :param refdict: dictionary containing empirical relation information

    Optional Inputs:
        None

    Output:
        A string containing SPLAT's default name for a given reference set, or False if that reference is not present

    Example:

    >>> import splat
    >>> print(splat.checkEmpiricalRelation('filippazzo',splat.SPT_LBOL_RELATIONS))
        filippazzo2015
    >>> print(splat.checkEmpiricalRelation('burgasser',splat.SPT_BC_RELATIONS))
        False
    '''
    output = False

# check reference    
    if not isinstance(ref,str):
        return output
    for k in list(refdict.keys()):
        if ref.lower()==k.lower() or ref.lower() in refdict[k]['altname']:
            output = k
    if output == False:
        if verbose: print('\nReference {} is not among those present in the reference dictionary; try: {}'.format(ref,list(refdict.keys())))

    return output



def checkInstrument(instrument):
    '''

    Purpose: 
        Checks that an instrument name is one of the available instruments, including a check of alternate names

    Required Inputs:
        :param: instrument: A string containing the instrument name to be checked. This should be one of the instruments in the global parameter splat.initialize.INSTRUMENTS
        
    Optional Inputs:
        None

    Output:
        A string containing SPLAT's default name for a given instrument, or False if that instrument is not present

    Example:

    >>> import splat
    >>> splat.checkInstrument('SPEX PRISM')
        SPEX-PRISM
    >>> splat.checkInstrument('LRIS')
        LRIS-RED
    >>> splat.checkInstrument('somethingelse')
        False
    '''
    return checkDict(instrument,INSTRUMENTS,replace=[['_','-'],[' ','-']])
    # output = False
    # if not isinstance(instrument,str):
    #     return output
    # for k in list(INSTRUMENTS.keys()):
    #     if instrument.upper()==k.upper() or instrument.upper().replace(' ','_').replace('_','-')==k.upper() or instrument.upper() in [a.upper() for a in INSTRUMENTS[k]['altname']]:
    #         output = k
    # return output


def checkFilterName(f,verbose=False):
    '''

    Purpose: 
        Checks that an input filter name is one of the available filters, including a check of alternate names

    Required Inputs:
        :param: filter: A string containing the filter name to be checked. This should be one of the names listed in `splat.FILTERS.keys()` or name alternates
        
    Optional Inputs:
        None

    Output:
        A string containing SPLAT's default name for a given filter, or False if that filter is not present

    Example:

    >>> import splat
    >>> print(splat.checkFilterName('2MASS_KS'))
        2MASS_KS
    >>> print(splat.checkFilterName('2mass k'))
        2MASS_KS
    >>> print(splat.checkFilterName('somethingelse'))
        False
    '''
    output = False
    if not isinstance(f,str):
        return output
    for k in list(FILTERS.keys()):
        if f.lower().replace(' ','_').replace('-','_') == k.lower() or f.lower().replace(' ','_') in [x.lower() for x in FILTERS[k]['altname']]:
            output = k
    if verbose==True and output==False:
        print('\nSPLAT does not contain the filter {}'.format(f))
    return output


def checkSpectralModelName(model):
    '''

    Purpose: 
        Checks that an input model name is one of the available spectral models, including a check of alternate names

    Required Inputs:
        :param: model: A string containing the spectral model to be checked. This should be one of the models listed in `loadModel()`_

    .. _`loadModel()` : api.html#splat_model.loadModel
        
    Optional Inputs:
        None

    Output:
        A string containing SPLAT's default name for a given model set, or False if that model set is not present

    Example:

    >>> import splat
    >>> print(splat.checkSpectralModelName('burrows'))
        burrows06
    >>> print(splat.checkSpectralModelName('allard'))
        BTSettl2008
    >>> print(splat.checkSpectralModelName('somethingelse'))
        False
    '''
    return checkDict(model,SPECTRAL_MODELS)

    # output = False
    # if not isinstance(model,str):
    #     return output
    # for k in list(SPECTRAL_MODELS.keys()):
    #     if model.lower()==k.lower() or model.lower() in SPECTRAL_MODELS[k]['altname']:
    #         output = k
    # return output


def checkEvolutionaryModelName(model):
    '''

    Purpose: 
        Checks that an input model name is one of the available evolutionary models, including a check of alternate names

    Required Inputs:
        :param: model: A string containing the evolutionary model to be checked. This should be one of the models listed in splat.EVOLUTIONARY_MODELS.keys()
        
    Optional Inputs:
        None

    Output:
        A string containing SPLAT's default name for a given model set, or False if that model set is not present

    Example:

    >>> import splat
    >>> print(splat.checkEvolutionaryModelName('burrows'))
        burrows01
    >>> print(splat.checkEvolutionaryModelName('allard'))
        False
    '''
    output = False
    if not isinstance(model,str):
        return output
    for k in list(EVOLUTIONARY_MODELS.keys()):
        if model.lower()==k.lower() or model.lower() in EVOLUTIONARY_MODELS[k]['altname']:
            output = k
    return output


def checkAbsMag(ref,filt='',verbose=False):
    '''

    Purpose: 
        Checks that an input reference name and filter are among the available sets for `typeToMag()`_, 
        including a check of alternate names

    .. _`typeToMag()` : TMP

    Required Inputs:
        :param ref: A string containing the reference for absolute magnitude relation, 
        among the keys and alternate names in splat.SPT_ABSMAG_RELATIONS

    Optional Inputs:
        :param filt: A string containing the filter name, to optionally check if this filter is among those defined in the reference set

    Output:
        A string containing SPLAT's default name for a given reference set, or False if that reference is not present

    Example:

    >>> import splat
    >>> print(splat.checkEvolutionaryModelName('burrows'))
        burrows01
    >>> print(splat.checkEvolutionaryModelName('allard'))
        False
    '''
    output = False

# check reference    
    if not isinstance(ref,str):
        return output
    for k in list(SPT_ABSMAG_RELATIONS.keys()):
        if ref.lower()==k.lower() or ref.lower() in SPT_ABSMAG_RELATIONS[k]['altname']:
            output = k
    if output == False:
        if verbose: print('\nReference {} is not among those used in SPLAT; try: {}'.format(ref,list(SPT_ABSMAG_RELATIONS.keys())))
        return output

# check filter
    if filt != '':
        filt = checkFilterName(filt)
        if filt == False:
            if verbose: print('\nFilter {} is not among the filters used in SPLAT; try: {}'.format(filt,list(FILTERS.keys())))
            return False
        if filt not in list(SPT_ABSMAG_RELATIONS[output]['filters'].keys()):
            if verbose: print('\nFilter {} is not among the filters defined for the {} absolutel magnitude relation; try: {}'.format(filt,output,list(SPT_ABSMAG_RELATIONS[output]['filters'].keys())))
            return False

    return output

def checkBC(ref,filt='',verbose=False):
    '''

    Purpose: 
        Checks that an input reference name and filter are among the available sets for `typeToBC()`_, 
        including a check of alternate names

    .. _`typeToBC()` : TMP

    Required Inputs:
        :param ref: A string containing the reference for absolute magnitude relation, 
        among the keys and alternate names in splat.SPT_BC_RELATIONS

    Optional Inputs:
        :param filt: A string containing the filter name, to optionally check if this filter is among those defined in the reference set

    Output:
        A string containing SPLAT's default name for a given reference set, or False if that reference is not present

    Example:

    >>> import splat
    >>> print(splat.checkBC('filippazzo','2MASS J'))
        filippazzo2015
    >>> print(splat.checkBC('dupuy','2MASS J'))
        False
    '''
    output = False

# check reference    
    if not isinstance(ref,str):
        return output
    for k in list(SPT_BC_RELATIONS.keys()):
        if ref.lower()==k.lower() or ref.lower() in SPT_BC_RELATIONS[k]['altname']:
            output = k
    if output == False:
        if verbose: print('\nReference {} is not among those used in SPLAT; try: {}'.format(ref,list(SPT_BC_RELATIONS.keys())))
        return output

# check filter
    if filt != '':
        filt = checkFilterName(filt)
        if filt == False:
            if verbose: print('\nFilter {} is not among the filters used in SPLAT; try: {}'.format(filt,list(FILTERS.keys())))
            return False
        if filt not in list(SPT_BC_RELATIONS[output]['filters'].keys()):
            if verbose: print('\nFilter {} is not among the filters defined for the {} absolutel magnitude relation; try: {}'.format(filt,output,list(SPT_BC_RELATIONS[output]['filters'].keys())))
            return False

    return output

def checkLbol(ref,verbose=False):
    '''

    Purpose: 
        Checks that an input reference name are among the available sets for `typeToLuminosity()`_, 
        including a check of alternate names

    .. _`typeToLuminosity()` : TMP

    Required Inputs:
        :param ref: A string containing the reference for lumiosity/SpT relation, 
        among the keys and alternate names in splat.SPT_LBOL_RELATIONS

    Optional Inputs:
        None

    Output:
        A string containing SPLAT's default name for a given reference set, or False if that reference is not present

    Example:

    >>> import splat
    >>> print(splat.checkLbol('filippazzo'))
        filippazzo2015
    >>> print(splat.checkBC('burgasser'))
        False
    '''
    output = False

# check reference    
    if not isinstance(ref,str):
        return output
    for k in list(SPT_LBOL_RELATIONS.keys()):
        if ref.lower()==k.lower() or ref.lower() in SPT_LBOL_RELATIONS[k]['altname']:
            output = k
    if output == False:
        if verbose: print('\nReference {} is not among those used in SPLAT; try: {}'.format(ref,list(SPT_LBOL_RELATIONS.keys())))
        return output

    return output



def checkTelescope(location):
    '''
    Purpose: 
        Checks that a location name is one of the telecopes listed in splat.initialize.TELESCOPES, including a check of alternate names

    Required Inputs:
        :param: location: A string containing the telescope/site name to be checked. This should be one of the locations in the global parameter splat.initialize.TELESCOPES
        
    Optional Inputs:
        None

    Output:
        A string containing SPLAT's default name for a given telescope, or False if that telecope is not present

    Example:

    >>> import splat
    >>> print(splat.checkTelescope('keck'))
        KECK
    >>> print(splat.checkTelescope('mauna kea'))
        KECK
    >>> print(splat.checkTelescope('somethingelse'))
        False
    '''
    output = False
    if not isinstance(location,str):
        return output
    for k in list(TELESCOPES.keys()):
        if location.upper().replace(' ','_').replace('-','_')==k.upper() or location.upper().replace(' ','_').replace('-','_') in [a.upper() for a in TELESCOPES[k]['altname']]:
            output = k
    return output

def checkLocation(location):
    '''

    Purpose: 
        Duplicate of checkTelescope()

    '''
    return checkTelescope(location)


#####################################################
##############   SIMPLE CONVERSIONS   ###############
#####################################################

#def caldateToDate(d):
    '''
    :Purpose: Convert from numeric date to calendar date, and vice-versa.
    :param d: A numeric date of the format '20050412', or a date in the
                calendar format '2005 Jun 12'
    :Example:
       >>> import splat
       >>> caldate = splat.dateToCaldate('20050612')
       >>> print caldate
       2005 Jun 12
       >>> date = splat.caldateToDate('2005 June 12')
       >>> print date
       20050612
    '''
#    return properDate(d,output='YYYY MMM DD')


#def dateToCaldate(d):
    '''
    :Purpose: Converts numeric date to calendar date

    :param date: String in the form 'YYYYMMDD'
    :type date: required

    :Output: Date in format YYYY MMM DD

    :Example:
    >>> import splat
    >>> splat.dateToCaldate('19940523')
        1994 May 23
    '''
#    d1 = copy.deepcopy(d)
#    if isNumber(d1): d1 = str(d1)
#    return d1[:4]+' '+MONTHS[int(d1[5:6])-1]+' '+d1[-2:]


def properDate(din,**kwargs):
    '''
    :Purpose: Converts various date formats into a standardized date of YYYY-MM-DD

    :param d: Date to be converted.
    :param format: Optional input format of the following form:
        * 'YYYY-MM-DD': e.g., 2011-04-03 (this is default output)
        * 'YYYYMMDD': e.g., 20110403
        * 'YYMMDD': e.g., 20110403
        * 'MM/DD/YY': e.g., 03/04/11
        * 'MM/DD/YYYY': e.g., 03/04/2011
        * 'YYYY/MM/DD': e.g., 2011/03/04
        * 'DD/MM/YYYY': e.g., 04/03/2011
        * 'DD MMM YYYY': e.g., 04 Mar 2011
        * 'YYYY MMM DD': e.g., 2011 Mar 04
    :type format: Optional, string
    :param output: Format of the output based on the prior list
    :type output: Optional, string

    :Example:
    >>> import splat
    >>> splat.properDate('20030502')
        '2003-05-02'
    >>> splat.properDate('2003/05/02')
        '02-2003-05'
    >>> splat.properDate('2003/05/02',format='YYYY/MM/DD')
        '2003-05-02'
    >>> splat.properDate('2003/05/02',format='YYYY/MM/DD',output='YYYY MMM DD')
        '2003 May 02'

    Note that the default output format can be read into an astropy.time quantity
    >>> import splat
    >>> from astropy.time import Time
    >>> t = Time(splat.properDate('20030502'))
    >>> print(t)
        2003-05-02 00:00:00.000
    '''

    dformat = kwargs.get('format','')
    oformat = kwargs.get('output','YYYY-MM-DD')
    d = copy.deepcopy(din)
    if not isinstance(d,str): d = str(int(d))
    if len(d)==0:
        print('\nCould not determine format of input date {}; please provide a format string\n'.format(din))
        return ''        

# some defaults
    if '/' in d and dformat == '':       # default American style
        if len(d) <= 8:
            dformat = 'MM/DD/YY'
        else:
            dformat = 'MM/DD/YYYY'
    if True in [c.lower() in d.lower() for c in MONTHS] and dformat == '':
        if isNumber(d.replace(' ','')[3]):
            dformat = 'YYYY MMM DD'
        else:
            dformat = 'DD MMM YYYY'
    if 'T' in d and dformat == '':       # default American style
        d = d.split('T')[0]
    if isNumber(d) and dformat == '':
        if len(str(d)) <= 6:
            dformat = 'YYMMDD'
        else:
            dformat = 'YYYYMMDD'            

# no idea
    if dformat == '':
        print('\nCould not determine format of input date {}; please provide a format string\n'.format(din))
        return ''

# case statement for conversion to YYYY-MM-DD
    if dformat == 'YYYYMMDD':
        dp = d[:4]+'-'+d[4:6]+'-'+d[-2:]
    elif dformat == 'YYMMDD':
        if int(d[:2]) > 50:
            dp = '19'+d[:2]+'-'+d[2:4]+'-'+d[-2:]
        else:
            dp = '20'+d[:2]+'-'+d[2:4]+'-'+d[-2:]
    elif dformat == 'MM/DD/YYYY':
        tmp = d.split('/')
        if len(tmp[0]) == 1:
            tmp[0] = '0'+tmp[0]
        if len(tmp[1]) == 1:
            tmp[1] = '0'+tmp[1]
        dp = tmp[2]+'-'+tmp[0]+'-'+tmp[1]
    elif dformat == 'MM/DD/YY':
        tmp = d.split('/')
        if len(tmp[0]) == 1:
            tmp[0] = '0'+tmp[0]
        if len(tmp[1]) == 1:
            tmp[1] = '0'+tmp[1]
        if int(tmp[2]) > 50:
            dp = '19'+tmp[2]+'-'+tmp[0]+'-'+tmp[1]
        else:
            dp = '20'+tmp[2]+'-'+tmp[0]+'-'+tmp[1]
    elif dformat == 'YYYY/MM/DD':
        tmp = d.split('/')
        if len(tmp[2]) == 1:
            tmp[2] = '0'+tmp[2]
        if len(tmp[1]) == 1:
            tmp[1] = '0'+tmp[1]
        dp = tmp[0]+'-'+tmp[1]+'-'+tmp[2]
    elif dformat == 'DD/MM/YYYY':
        tmp = d.split('/')
        if len(tmp[0]) == 1:
            tmp[0] = '0'+tmp[0]
        if len(tmp[1]) == 1:
            tmp[1] = '0'+tmp[1]
        dp = tmp[2]+'-'+tmp[1]+'-'+tmp[0]
    elif dformat == 'DD/MM/YY':
        tmp = d.split('/')
        if len(tmp[0]) == 1:
            tmp[0] = '0'+tmp[0]
        if len(tmp[1]) == 1:
            tmp[1] = '0'+tmp[1]
        if int(tmp[2]) > 50:
            dp = '19'+tmp[2]+'-'+tmp[1]+'-'+tmp[0]
        else:
            dp = '20'+tmp[2]+'-'+tmp[1]+'-'+tmp[0]
    elif dformat == 'DD MMM YYYY':
        tmp = d.split(' ')
        if len(tmp[0]) == 1:
            tmp[0] = '0'+tmp[0]
        for i,c in enumerate(MONTHS):
            if c.lower() == tmp[1].lower():
                mref = str(i+1)
        if len(mref) == 1:
            mref = '0'+mref
        dp = tmp[2]+'-'+mref+'-'+tmp[0]
    elif dformat == 'DD-MMM-YYYY':
        tmp = d.split(' ')
        if len(tmp[0]) == 1:
            tmp[0] = '0'+tmp[0]
        for i,c in enumerate(MONTHS):
            if c.lower() == tmp[1].lower():
                mref = str(i+1)
        if len(mref) == 1:
            mref = '0'+mref
        dp = tmp[2]+'-'+mref+'-'+tmp[0]
    elif dformat == 'YYYY MMM DD':
        tmp = d.split(' ')
        if len(tmp[2]) == 1:
            tmp[2] = '0'+tmp[2]
        for i,c in enumerate(MONTHS):
            if c.lower() == tmp[1].lower():
                mref = str(i+1)
        if len(mref) == 1:
            mref = '0'+mref
        dp = tmp[0]+'-'+mref+'-'+tmp[2]
    elif dformat == 'YYYY-MMM-DD':
        tmp = d.split(' ')
        if len(tmp[2]) == 1:
            tmp[2] = '0'+tmp[2]
        for i,c in enumerate(MONTHS):
            if c.lower() == tmp[1].lower():
                mref = str(i+1)
        if len(mref) == 1:
            mref = '0'+mref
        dp = tmp[0]+'-'+mref+'-'+tmp[2]
    else:
        dp = d

# case statement for conversion from YYYY-MM-DD to desired output format
    if oformat == 'YYYYMMDD':
        df = dp.replace('-','')
    elif oformat == 'YYMMDD':
        df = dp.replace('-','')[2:]
    elif oformat == 'MM/DD/YYYY':
        tmp = dp.split('-')
        df = tmp[1]+'/'+tmp[2]+'/'+tmp[0]
    elif oformat == 'MM/DD/YY':
        tmp = dp.split('-')
        df = tmp[1]+'/'+tmp[2]+'/'+tmp[0][2:]
    elif oformat == 'YYYY/MM/DD':
        tmp = dp.split('-')
        df = tmp[0]+'/'+tmp[1]+'/'+tmp[2]
    elif oformat == 'DD/MM/YYYY':
        tmp = dp.split('-')
        df = tmp[2]+'/'+tmp[1]+'/'+tmp[0]
    elif oformat == 'DD/MM/YY':
        tmp = dp.split('-')
        df = tmp[2]+'/'+tmp[1]+'/'+tmp[0][2:]
    elif oformat == 'DD MMM YYYY':
        tmp = dp.split('-')
        df = tmp[2]+' '+MONTHS[int(tmp[1])-1]+' '+tmp[0]
    elif oformat == 'DD-MMM-YYYY':
        tmp = dp.split('-')
        df = tmp[2]+'-'+MONTHS[int(tmp[1])-1]+'-'+tmp[0]
    elif oformat == 'YYYY MMM DD':
        tmp = dp.split('-')
        df = tmp[0]+' '+MONTHS[int(tmp[1])-1]+' '+tmp[2]
    elif oformat == 'YYYY-MMM-DD':
        tmp = dp.split('-')
        df = tmp[0]+'-'+MONTHS[int(tmp[1])-1]+'-'+tmp[2]
    else:
        df = dp

    return df




def checkKeys(input,parameters,**kwargs):
    '''
    :Purpose: Checks the input kwargs keys against the expected parameters of a function to make sure the right parameters are passed.

    :param input: input dictionary to a function (i.e., kwargs).
    :param parameters: allowed parameters for the function
    :param forcekey: (optional, default = False) if True, raises a Value Error if an incorrect parameter is passed
    '''
    kflag = False
    forcekey = kwargs.get('forcekey',False)
    for k in input.keys():
        if k not in parameters:
            print('\nParameter Warning!\nUnknown input keyword {}'.format(k))
            kflag = True
    if kflag:
        if forcekey:
            raise ValueError('Valid keywords are {}\n'.format(parameters))
        else:
            print('Valid keywords are {}\n'.format(parameters))


def coordinateToDesignation(c,prefix='J',sep='',split='',decimal=False):
    '''
    :Purpose: Converts right ascension and declination into a designation string

    :param c: RA and Dec coordinate to be converted; can be a SkyCoord object with units of degrees,
              a list with RA and Dec in degrees, or a string with RA measured in hour
              angles and Dec in degrees

    :Output: Designation string

    :Example:
    >>> import splat
    >>> from astropy.coordinates import SkyCoord
    >>> c = SkyCoord(238.86, 9.90, unit="deg")
    >>> print splat.coordinateToDesignation(c)
        J15552640+0954000
    >>> print splat.coordinateToDesignation([238.86, 9.90])
        J15552640+0954000
    >>> print splat.coordinateToDesignation('15:55:26.4 +09:54:00.0')
        J15552640+0954000
    '''
# input is ICRS
    # decreplace = ''
    # if decimal==True: decreplace='.'
    if isinstance(c,SkyCoord):
        cc = copy.deepcopy(c)
    else:
        cc = properCoordinates(c)
# input is [RA,Dec] pair in degrees
    output = '{}{}{}{}'.format(prefix, cc.ra.to_string(unit=u.hour, sep=sep, precision=2, pad=True), \
        split , cc.dec.to_string(unit=u.degree, sep=sep, precision=1, alwayssign=True, pad=True))
    if decimal==False: output = output.replace('.','')
    # if sys.version_info.major == 2:
    #     return string.replace('{}{0}{}{1}'.format(prefix,cc.ra.to_string(unit=u.hour, sep=sep, precision=2, pad=True), \
    #     splitstr, cc.dec.to_string(unit=u.degree, sep=sep, precision=1, alwayssign=True, pad=True)),'.',decreplace)
    # else:
    #     return str.replace('{}{0}{}{1}'.format(prefix,cc.ra.to_string(unit=u.hour, sep=sep, precision=2, pad=True), \
    #     splitstr, cc.dec.to_string(unit=u.degree, sep=sep, precision=1, alwayssign=True, pad=True)),'.',decreplace)
    return output


def designationToCoordinate(value, icrs=True, **kwargs):
    '''
    :Purpose: Convert a designation string into a RA, Dec tuple or ICRS SkyCoord objects (default)

    :param value: Designation string with RA measured in hour angles and Dec in degrees
    :type value: required
    :param icrs: returns astropy SkyCoord coordinate in ICRS frame if ``True``
    :type icrs: optional, defualt = True

    :Output: Coordinate, either as [RA, Dec] or SkyCoord object

    :Example:
    >>> import splat
    >>> splat.designationToCoordinate('J1555264+0954120')
    <SkyCoord (ICRS): (ra, dec) in deg
        (238.8585, 9.90333333)>
    '''

# remove any unnecessary symbols or trailing letters
    a = re.sub('[j.:hms]','',str(value).lower())
    a = re.sub(' ','0',a)
    while not splat.isNumber(a[-1]): a=a[:-1]
    fact = 1.
    spl = a.split('+')
    if len(spl) == 1:
        spl = a.split('-')
        fact = -1.
    if len(spl) == 1:
        raise ValueError('Input quantity {} is not a proper designation (missing +/-)'.format(value))
    ra = 15.*float(spl[0][0:2])
    if (len(spl[0]) > 2):
        ra+=15.*float(spl[0][2:4])/60.
    if (len(spl[0]) > 4):
        ra+=15.*float(spl[0][4:6])/3600.
    if (len(spl[0]) > 6):
        ra+=15.*float(spl[0][6:8])/360000.
    dec = float(spl[1][0:2])
    if (len(spl[1]) > 2):
        dec+=float(spl[1][2:4])/60.
    if (len(spl[1]) > 4):
        dec+=float(spl[1][4:6])/3600.
    if (len(spl[1]) > 6):
        dec+=float(spl[1][6:8])/360000.
    dec*=fact
    if icrs == True: return SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    else: return [ra,dec]


def designationToCoordinateString(designation,delimiter=' ',radec_delimiter=' '):
    '''
    :Purpose: 

        Convert a designation string into a coordinate string with delimiters between hour, minute, second, etc.

    :Required Inputs: 

        :param designation: designation, which should be a string of the form 'J12345678+01234567'

    :Optional Inputs: 

        :param: delimiter = ' ': delimiter between coordinate elements
        :param: radec_delimiter = ' ': delimiter between RA and declination substrings

    :Output: 

        coordinate string of the form '12 34 56.78 +01 23 45.67' (depending on delimiters)

    :Example:
    >>> import splat
    >>> splat.designationToCoordinateString('J1555264+0954120')
    15 55 26.4 +09 54 12.0
    >>> splat.designationToCoordinateString('J155526400+095412000',delimiter=':')
    15 55 26.400 +09 54 12.000
    '''
    if not isinstance(designation,str): raise ValueError('Input variable must be a string')

    d = designation.replace('J','').replace('j','').replace('.','')
    dsym = '+'
    tmp = d.split(dsym)
    if len(tmp) != 2: 
        dsym = '-'
        tmp = d.split(dsym)
    if len(tmp) != 2: raise ValueError('problem processing designation string {}'.format(d))
    output = tmp[0][0:2]+delimiter+tmp[0][2:4]+delimiter+tmp[0][4:6]
    if len(tmp[0]) > 6: output = output+'.'+tmp[0][6:]
    output = output+radec_delimiter+dsym+tmp[1][0:2]+delimiter+tmp[1][2:4]+delimiter+tmp[1][4:6]
    if len(tmp[1]) > 6: output = output+'.'+tmp[1][6:]

    return output


def designationToShortName(value):
    '''
    :Purpose: Produce a shortened version of designation

    :param value: Designation string with RA measured in hour angles and Dec in degrees
    :type value: required

    :Output: Shorthand designation string

    :Example:
    >>> import splat
    >>> print splat.designationToShortName('J1555264+0954120')
        J1555+0954
    '''
    if isinstance(value,str):
        a = re.sub('[j.:hms]','',value.lower())
        mrk = '+'
        spl = a.split(mrk)
        if len(spl) == 1:
            mrk = '-'
            spl = a.split(mrk)
        if len(spl) == 2:
            return 'J'+spl[0][0:4]+mrk+spl[1][0:4]
        else:
            return value
    else:
        raise ValueError('\nMust provide a string value for designation\n\n')



def properCoordinates(c,frame='icrs',icrs=True,**kwargs):
    '''
    :Purpose: Converts various coordinate forms to the proper SkyCoord format. Convertible forms include lists and strings.

    :param c: coordinate to be converted. Can be a list (ra, dec) or a string.

    :Example:
    >>> import splat
    >>> print splat.properCoordinates([104.79, 25.06])
        <SkyCoord (ICRS): ra=104.79 deg, dec=25.06 deg>
    >>> print splat.properCoordinates('06:59:09.60 +25:03:36.0')
        <SkyCoord (ICRS): ra=104.79 deg, dec=25.06 deg>
    >>> print splat.properCoordinates('J06590960+2503360')
        <SkyCoord (ICRS): ra=104.79 deg, dec=25.06 deg>
    '''
    if isinstance(c,SkyCoord):
        output = c
    elif isinstance(c,list):
        output = SkyCoord(c[0]*u.deg,c[1]*u.deg,frame=frame)

# input is sexigessimal string - assumed ICRS
    elif isinstance(c,str):
        if c[0] == 'J':
            output = designationToCoordinate(c,**kwargs)
        else:
            output = SkyCoord(c,frame='icrs', unit=(u.hourangle, u.deg))
    else:
        raise ValueError('\nCould not parse input format\n\n')

# add distance
    if kwargs.get('distance',False) != False:
        d = copy.deepcopy(kwargs['distance'])
        if not isUnit(d): d=d*u.pc
        d.to(u.pc)
        output = SkyCoord(output,distance = d)
#        except:
#            print('\nWarning: could not integrate distance {} into coordinate'.format(distance))

# convert to icrs by default
    if icrs == True: return output.icrs
    else: return output



def typeToNum(inp, subclass='dwarf', error='', uncertainty=0., luminosity_class = '', metallicity_class='', age_class = '', color_class='', peculiar=False, verbose=False, **kwargs):
    '''
    :Purpose: 

        Converts between string and numeric spectral types, with the option of specifying the class prefix/suffix and uncertainty tags

    :Required inputs: 

        :param inp: Spectral type to convert. Can convert a number or a string from 0.0 (K0) and 49.0 (Y9).

    :Optional inputs: 

        :param: error = '': flag to indicate magnitude of classification uncertainty; by default ':' for uncertainty > 1 subtypes and '::' for uncertainty > 2 subtype added as suffix to string output. Can also use `err`.
        :param: uncertainty = 0: numerical uncertainty of classification; can also use `unc`
        :param: subclass = 'dwarf': spectral class; options include:

            - *field* or *fld* or *alpha*: object is a field dwarf - no prefix/suffix to string output
            - *sd* or *subdwarf*: object is a subdwarf - 'sd' prefix to string output
            - *dsd* or *d/sd*: object is an intermediate subdwarf - 'd/sd' prefix to string output
            - *esd*: object is an extreme subdwarf - 'esd' prefix to string output
            - *usd*: object is an ultra subdwarf - 'usd' prefix to string output
            - *delta*: object is a extremely low surface gravity dwarf (~1 Myr) - 'delta' suffix to string output
            - *vlg* or *gamma* or *lowg*: object is a low surface gravity dwarf (~10 Myr) - 'gamma' suffix to string output
            - *intg* or *beta*: object is an intermediate surface gravity dwarf (~100 Myr) - 'beta' suffix to string output
            - *giant*: object is a giant with luminosity class III suffix added to string output
            - *subgiant*: object is a subgiant with luminosity class IV suffix added to string output
            - *supergiant*: object is a supergiant with luminosity class I suffix added to string output

        :param: metallicity_class = '': metallicity class of object, traditionally represented by 'sd','d/sd','esd','usd', and added on as prefix to string output. Can also use `lumclass`
        :param: luminosity_class = '': luminosity class of object traditionally represented by roman numerals (e.g., 'III') and added on as suffix to string output. Can also use `lumclass`
        :param: age_class = '': age class of object, traditionally one of 'alpha', 'beta', 'gamma', 'delta' and added on as suffix to string output (see subclass). Can also use 'ageclass'
        :param: color_class: color class of object, traditionally 'b' (for blue) or 'r' (for red), added as prefix to string output. Can also use 'colorclass'
        :param: peculiar = False: Set to True if object is peculiar, which adds a 'pec' suffix to string output
        :param: verbose = False: Set to True to provide more feedback

    :Outputs: 

        The number or string of a spectral type

    :Example:
        >>> import splat
        >>> print splat.typeToNum(30)
            T0.0
        >>> print splat.typeToNum('T0.0')
            30.0
        >>> print splat.typeToNum(27, peculiar = True, uncertainty = 1.2, lumclass = 'II')
            L7.0IIp:
        >>> print splat.typeToNum(50)
            Spectral type number must be between 0 (K0) and 49.0 (Y9)
            nan
    '''
# keywords
    error = kwargs.get('err','')
    uncertainty = kwargs.get('unc',uncertainty)
    luminosity_class = kwargs.get('lumclass',luminosity_class)
    metallicity_class = kwargs.get('z_class',metallicity_class)
    metallicity_class = kwargs.get('metal_class',metallicity_class)
    age_class = kwargs.get('ageclass',age_class)
    colorclass = kwargs.get('colorclass','')
    peculiar = kwargs.get('peculiar',False)
    spletter = 'KMLTY'

# as of 12/18/2017, this only works on individual inputs
    if isinstance(inp,list):
        raise ValueError('\nInput to typeToNum() must be a single element (string or number)')

# convert input into an array
    # output = []
    # var = copy.deepcopy(inp)

    # if not isinstance(var,list): var = [var]
    # if not isinstance(error,list): error = [error]
    # if not isinstance(unc,list): unc = [unc]
    # if not isinstance(subclass,list): subclass = [subclass]
    # if not isinstance(lumclass,list): lumclass = [lumclass]
    # if not isinstance(ageclass,list): ageclass = [ageclass]
    # if not isinstance(colorclass,list): colorclass = [colorclass]
    # if len(error) < len(var):
    #     for i in numpy.arange(len(var)-len(error)): error.append(error[-1])
    # if len(unc) < len(var):
    #     for i in numpy.arange(len(var)-len(unc)): unc.append(unc[-1])

# number -> spectral type
    if isNumber(inp):
        # if len(subclass) < len(var):
        #     for i in numpy.arange(len(var)-len(subclass)): subclass.append(subclass[-1])
        # if len(lumclass) < len(var):
        #     for i in numpy.arange(len(var)-len(lumclass)): lumclass.append(lumclass[-1])
        # for i,l in enumerate(lumclass):
        #     if l != '': lumclass[i]=' '+lumclass[i]
        # if len(ageclass) < len(var):
        #     for i in numpy.arange(len(var)-len(ageclass)): ageclass.append(ageclass[-1])
        # if len(colorclass) < len(var):
        #     for i in numpy.arange(len(var)-len(colorclass)): colorclass.append(colorclass[-1])

        spind = int(abs(inp/10.))
        if spind < 0 or spind >= len(spletter): 
            if verbose: print('Spectral type number must be between 0 ({}0) and {} ({}9)'.format(spletter[0],len(spletter)*10.-1.,spletter[-1]))
            return 'N/A'
        spdec = numpy.around(inp,1)-spind*10.
# deal with subclasses
        if subclass.lower() == 'sd' or subclass.lower() == 'subdwarf': metallicity_class = 'sd'
        if subclass.lower() == 'dsd' or subclass.lower() == 'd/sd': metallicity_class = 'd/sd'
        if subclass.lower() == 'esd': metallicity_class = 'esd'
        if subclass.lower() == 'usd': metallicity_class = 'usd'
        if subclass.lower() == 'giant': luminosity_class = 'III'
        if subclass.lower() == 'subgiant': luminosity_class = 'IV'
        if subclass.lower() == 'supergiant': luminosity_class = 'I'
        if subclass.lower() == 'delta': age_class = 'delta'
        if subclass.lower() == 'vlg' or subclass.lower() == 'vl-g' or subclass.lower() == 'lowg' or subclass.lower() == 'low-g' or subclass.lower() == 'gamma': age_class = 'gamma'
        if subclass.lower() == 'intg' or subclass.lower() == 'int-g' or subclass.lower() == 'beta': age_class = 'beta'
        if uncertainty > 1.: error = ':'
        if uncertainty > 2.: error = '::'
        pstr = ''
        if peculiar == True: pstr = 'p'

        return '{}{}{}{:3.1f}{}{}{}{}'.format(color_class,metallicity_class,spletter[spind],spdec,age_class,luminosity_class,pstr,error)


# spectral type -> number
    elif isinstance(inp,str):
#        output = []
        if (sys.version_info.major == 2):
            inp = string.split(inp,sep='+/-')[0]    # remove +/- sides
        else:
            inp = inp.split('+/-')[0]    # remove +/- sides
        inp = inp.replace('...','').replace(' ','')

        sptype = re.findall('[{}]'.format(spletter),inp.upper())
        outval = 0.

# specialty classes                
        if len(sptype) >= 1:
            ytype = re.findall('[abcd]',inp.split('p')[-1])
            if len(ytype) == 1: age_class = ytype[0]
            if inp.find('pec') != -1:
                 peculiar = True
                 inp.replace('pec','')
            if inp.find('p') != -1:
                 peculiar = True
                 inp.replace('p','')
            if inp.find('alpha') != -1:
                 age_class = 'alpha'
                 inp.replace('alpha','')
            if inp.find('beta') != -1:
                 age_class = 'beta'
                 inp.replace('beta','')
            if inp.find('gamma') != -1:
                 age_class = 'gamma'
                 inp.replace('gamma','')
            if inp.find('delta') != -1:
                 age_class = 'delta'
                 inp.replace('delta','')
            if inp.find('esd') != -1:
                 subclass = 'esd'
                 inp.replace('esd','')
            elif inp.find('usd') != -1:
                 subclass = 'usd'
                 inp.replace('usd','')
            elif inp.find('d/sd') != -1:
                 subclass = 'd/sd'
                 inp.replace('d/sd','')
            elif inp.find('sd') != -1:
                 subclass = 'sd'
                 inp.replace('sd','')
            if inp.count('I') > 0:
                 luminosity_class = ''.join(re.findall('I',inp))
                 inp.replace('I','')
            if inp.count(':') > 0:
                 error = ''.join(re.findall(':',inp))
                 inp.replace(':','')
            if inp[0] == 'b' or inp[0] == 'r':
                 color_class = inp[0]
                 inp.replace('b','')
                 inp.replace('r','')

            outval = spletter.find(sptype[0])*10.
            spind = inp.find(sptype[0])+1
            if spind < len(inp):
                if inp.find('.') < 0:
                    if isNumber(inp[spind]):
                        outval = outval+float(inp[spind])
                else:
                    try:
                        outval = outval+float(inp[spind:spind+3])
                        spind = spind+3
                    except:
                        if verbose: print('\nProblem converting input type {} to a numeric type'.format(inp))
                        outval = numpy.nan
            return outval

        else:
            if verbose: print('\nOnly spectral classes {} are handled by typeToNum'.format(spletter))
            return numpy.nan

# none of the above - return the input
    else:
        if verbose: print('\nWarning: could not recognize format of spectral type {}\n'.format(inp))
        return inp


def UVW(coord,distance,mu,rv,e_distance = 0.,e_mu = [0.,0.],e_rv = 0.,nsamp=100,full=False,verbose=False):
    '''
    THIS FUNCTION NEEDS CLEANING
    '''
    try:
        from uvwxyz.uvwxyz import uvw as uvwcalc
    except:
        raise ValueError('\nMust have installed package uvwxyz to run this module: https://github.com/dr-rodriguez/uvwxyz')
    try:
        c = properCoordinates(coord)
    except:
        raise ValueError('\nCoordinate input {} is in incorrect format'.format(coord))

    if not isinstance(mu,list) and not isinstance(mu,numpy.ndarray): 
        raise ValueError('\nProper motion input {} must be a 2-element list'.format(mu))
    if not isinstance(e_mu,list) and not isinstance(e_mu,numpy.ndarray): 
        raise ValueError('\nProper motion uncertainty input {} must be a 2-element list'.format(e_mu))

    u,v,w = uvwcalc(c.ra.degree,c.dec.degree,distance,mu[0],mu[1],rv)
    if e_distance!=0 or e_mu[0]!=0 or e_mu[1]!=0 or e_rv!=0:
        us,vs,ws = uvwcalc(c.ra.degree,c.dec.degree,numpy.random.normal(distance,e_distance,nsamp),numpy.random.normal(mu[0],e_mu[0],nsamp),numpy.random.normal(mu[1],e_mu[1],nsamp),numpy.random.normal(rv,e_rv,nsamp))
        return [u,numpy.std(us)],[v,numpy.std(vs)],[w,numpy.std(ws)]
    else: return [u,0],[v,0],[w,0]

def lbolToMbol(lbol,err=0.,scale='log',sun_scale=True,reverse=False):
    l0 = 3.0128e28*u.Watt # in watts
    lsun = u.Lsun

# Lbol -> Mbol
    if reverse==False:
        lb = copy.deepcopy(lbol)
        le = copy.deepcopy(err)
        if scale=='linear':
            if not isUnit(lb): 
                if sun_scale==True: lb=lb*lsun
                else: lb=lb*(l0.unit)
            lb = numpy.log10((lb/lsun).decompose())
            if not isUnit(le): 
                if sun_scale==True: le=le*lsun
                else: le=le*(l0.unit)
            le = numpy.log10((le/lsun).decompose())
        mout = -2.5*lb-2.5*numpy.log10((lsun/l0).decompse())
        mout_e = 2.5*le
        if err == 0.:
            return mout
        else:
            return mout,mout_e
    
# Mbol -> Lbol
    else:
        mb = copy.deepcopy(lbol)
        mbe = copy.deepcopy(err)
        lout = l0*10.**(-0.4*mb)
        lout_e = lout*0.4*numpy.log(10)*mbe
        if scale=='linear':
            if err == 0.:
                return lout
            else:
                return lout,lout_e
        else:
            lout_e = ((lout_e/lout).decompose())/numpy.log(10.)
            lout = numpy.log10((lout/lsun).decompose())
            if err == 0.: return lout
            else: return lout.value,lout_e.value


def xyz(coordinate,center='sun',r0=8000*u.pc,z0=25*u.pc,unit=u.pc,**kwargs):
    '''
    :Purpose:

        A "fast" method for converting a coordinate to heliocentric or galactocentric XYZ (cartesian) galaxy coordinates. 
        This assumes a right handed orientation with X from Sun to Galactic center, Y from Sun to the direction of Galactic rotation, and Z from Sun toward Galactic North. 
        Note that the astropy SkyCoord method also provides a way of producing `XYZ equatorial coordinates <http://docs.astropy.org/en/stable/api/astropy.coordinates.CartesianRepresentation.html>`_

    :Required Inputs:

        :param coordinate: A coordinate or list of coordinate variables, something that can be converted to astropy SkyCoord by `splat.properCoordinates()`_

    :Optional Inputs:

        :param distance: If not included in the coordinate variable, the distance to the source in pc (default: None)
        :param center = 'sun': centering of coordinates; by default this is the Sun, but for full galacitic coordindates set to 'galactic'
        :param r0 = 8000 pc: radial distance between Sun and Galactic center 
        :param z0 = 25 pc: vertical distance between Sun and Galactic plane 
        :param unit = astropy.units.pc: preferred unit

    :Outputs:

        A tuple (x,y,z), each of which is an array of x,y,z Galactic coordinates in preferred units

    :Example:

        >>> import splat
        >>> c = splat.properCoordinates('J05591914-1404488',distance=10.2)
        >>> splat.xyz(c)
            (<Quantity -7.442377515807463 pc>, <Quantity -6.2399837133240235 pc>, <Quantity -3.116668119908577 pc>)
        >>> splat.xyz(c,center='galactic')
            (<Quantity 7992.5576224841925 pc>, <Quantity -6.2399837133240235 pc>, <Quantity 21.883331880091422 pc>)

    .. _`splat.properCoordinates() <REF>`        
    
    '''
# check inputs
    if not splat.isUnit(unit): unit = u.pc

    if not isinstance(coordinate,list): c = [coordinate]
    else: c = coordinate
    if not isinstance(c[0],SkyCoord):
        try:
            c = [splat.properCoordinates(cd,**kwargs) for cd in c]
        except:
            raise ValueError('{} is not a proper coordinate'.format(coordinate))

    if not isinstance(kwargs.get('distance',False),bool): distance=kwargs['distance']
    elif str(c[0].distance.unit) != '': distance = [cd.distance for cd in c]
    else:
        raise ValueError('No distance value provided')
    if isinstance(distance,numpy.ndarray): distance = list(distance)
    if not isinstance(distance,list): distance = [distance]
    if splat.isUnit(distance[0]): distance = [float(d.to(unit).value) for d in distance]

    if splat.isUnit(r0): r0 = r0.to(unit).value
    if splat.isUnit(z0): z0 = z0.to(unit).value

    l = [cd.galactic.l.radian for cd in c]
    b = [cd.galactic.b.radian for cd in c]

# make sure arrays are of the same length 
    while len(distance) < len(l): distance.append(distance[-1])
    while len(l) < len(distance): 
        l.append(l[-1])
        b.append(b[-1])

# compute xyz        
    distance = numpy.array(distance)
    l = numpy.array(l)
    b = numpy.array(b)
    x = distance*numpy.cos(l)*numpy.cos(b)
    y = distance*numpy.sin(l)*numpy.cos(b)
    z = distance*numpy.sin(b)

    if center.lower() == 'galactic':
        x = x+r0
        z = z+z0

    if len(distance) == 1:
        return x[0]*unit,y[0]*unit,z[0]*unit
    else:
        return x*unit,y*unit,z*unit


def baryVel(coord,obstime,location='keck',correction='barycenter'):
    '''
    :Purpose: 

        Computes the barycentric or heliocentric velocity in a direction and from a specific Earth location

    :Required Inputs:

        - :param coord: Coordinate of source; should be astropy.coordinates.SkyCoord, but can also be converted from splat.propoCoordinates
        - :param obstime: A date/time, preferred in astropy.time.Time format but can be converted from splat.properDate

    :Optional Inputs:

        - :param location: location on Earth, specified by astropy.coordinates.EarthLocation; string of location; 
        dictionary containing 'ra', 'dec', and 'height'; or array of [ra,dec,height] (default = 'keck')
        - :param correction: type of correction, can be either 'barycentric' or 'heliocentric' (default = 'heliocentric')

    :Output:

        The velocity correction in km/s

    :Example:
        >>> import splat
        >>> coord = splat.properCoordinates('J15104786-2818174')
        >>> print(splat.baryVel(coord,'2017-07-31',location='keck')
            -27.552554878923033 km / s
    '''
# check coordinate
    if not isinstance(coord,SkyCoord):
        try:
            c = properCoordinates(coord)
        except:
            raise ValueError('\nCould not convert coordinate input {} to a SkyCoord'.format(coord))
    else: c = copy.deepcopy(coord)

# check time
    if not isinstance(obstime,Time):
        try:
            t = Time(obstime)
        except:
            raise ValueError('\nCould not convert time input {} into a Time variable'.format(obstime))
    else: t = copy.deepcopy(obstime)

# check location
    if not isinstance(location,EarthLocation):
        if isinstance(location,str):
            loc = checkTelescope(location)
            if loc != False:
                l = EarthLocation.from_geodetic(lat=TELESCOPES[loc]['lat'], lon=TELESCOPES[loc]['lon'], height=TELESCOPES[loc]['height'])
            else:
                try:
                    l =  EarthLocation.of_site(location)
                except:
                    raise ValueError('\nCould not convert location input {} into an EarthLocation; may be offline'.format(location))
        elif isinstance(location,list) or isinstance(location,float):
            try:
                if len(location) == 2:
                    if not isUnit(l[0]): location = [x*u.deg for x in l]
                    l = EarthLocation.from_geodetic(lat=location[0], lon=location[1])
                elif len(location) == 3:
                    if not isUnit(location[0]): 
                        location[0] = l[0]*u.deg
                        location[1] = l[1]*u.deg
                        location[2] = l[2]*u.m
                    l = EarthLocation.from_geodetic(lat=location[0], lon=location[1], height=location[2])
                else:
                    raise ValueError('\nCould not convert location input {} into an EarthLocation'.format(location))
            except:
                raise ValueError('\nCould not convert location input {} into an EarthLocation'.format(location))
        elif isinstance(location,dict):
            try:
                l = EarthLocation.from_geodetic(**location)
            except:
                raise ValueError('\nCould not convert location input {} into an EarthLocation'.format(location))
        else:
            raise ValueError('\nCould not convert location input {} into an EarthLocation'.format(location))
    else: l = copy.deepcopy(location)

# flag if we're not online
    auto_max_age = 14.*u.day
    if checkOnline() == False: auto_max_age = None

# make correction
    if 'bary' in correction.lower():
        return c.radial_velocity_correction(obstime=t, location=l).to(u.km/u.s)  

    elif 'helio' in correction.lower():
        return c.radial_velocity_correction('heliocentric',obstime=t, location=l).to(u.km/u.s)  
    else:
        raise ValueError('\n Could not interpret preferred correction {} '.format(correction))


def lsfRotation(vsini,vsamp,epsilon=0.6):
    '''
    Purpose: 

        Generates a line spread function for rotational broadening, based on Gray (1992) 
        Ported over by Chris Theissen and Adam Burgasser from the IDL routine 
        `lsf_rotate <https://idlastro.gsfc.nasa.gov/ftp/pro/astro/lsf_rotate.pro>`_ writting by W. Landsman

    Required Inputs:  

        :param: **vsini**: vsini of rotation, assumed in units of km/s
        :param: **vsamp**: sampling velocity, assumed in unit of km/s. vsamp must be smaller than vsini or else a delta function is returned

    Optional Inputs:

        :param: **epsilon**: limb darkening parameter based on Gray (1992)

    Output:

        Line spread function kernel with length 2*vsini/vsamp (forced to be odd)

    :Example:
        >>> import splat
        >>> kern = lsfRotation(30.,3.)
        >>> print(kern)
            array([ 0.        ,  0.29053574,  0.44558751,  0.55691445,  0.63343877,
            0.67844111,  0.69330989,  0.67844111,  0.63343877,  0.55691445,
            0.44558751,  0.29053574,  0.        ])
    '''
# limb darkening parameters
    e1 = 2. * (1. - epsilon)
    e2 = numpy.pi * epsilon/2.
    e3 = numpy.pi * (1. - epsilon/3.)

# vsini must be > vsamp - if not, return a delta function
    if vsini <= vsamp:
        print('\nWarning: velocity sampling {} is broader than vsini {}; returning delta function')  
        lsf = numpy.zeros(5)  
        lsf[2] = 1.
        return lsf

# generate LSF
    nsamp = numpy.ceil(2.*vsini/vsamp)
    if nsamp % 2 == 0:
        nsamp+=1
    x = numpy.arange(nsamp)-(nsamp-1.)/2.
    x = x*vsamp/vsini
    x2 = numpy.absolute(1.-x**2)

    return (e1*numpy.sqrt(x2) + e2*x2)/e3

#####################################################
############   STATISTICAL FUNCTIONS   ##############
#####################################################


def distributionStats(x, q=[0.16,0.5,0.84], weights=None, sigma=None, **kwargs):
    '''
    :Purpose: Find key values along distributions based on quantile steps.
              This code is derived almost entirely from triangle.py.
    '''

# clean data of nans
    xd = numpy.array(copy.deepcopy(x))
    xd0 = copy.deepcopy(xd)
    xd = xd[~numpy.isnan(xd)]

    if q is None and sigma is None:
        sigma = 1.
        
    if sigma is not None:
        q = [stats.norm.cdf(-sigma),0.5,stats.norm.cdf(sigma)]
        
    if weights is None:
        return numpy.percentile(xd, [100. * qi for qi in q])
    else:
        wt = numpy.array(copy.deepcopy(weights))
        wt = wt[~numpy.isnan(xd0)]
        idx = numpy.argsort(xd)
        xsorted = xd[idx]
        cdf = numpy.add.accumulate(wt[idx])
#        print(xsorted,cdf,wt[idx],type(xd),type(cdf))
        cdff = [float(c) for c in cdf]
        cdfn = [c/cdff[-1] for c in cdff]
        return numpy.interp(q, cdfn, xsorted).tolist()

def gauss(x,*p):
    '''
    Simple gaussian function for curve fit analysis
    '''

    A,mu,sig,c = p
    return c+A*numpy.exp(-(x-mu)**2/(2*sig**2))


def reMap(x1,y1,x2,nsamp=100,method='fast'):
    '''
    :Purpose: 

        Maps a function y(x) onto a new grid x'. If x' is higher resolution this is done through interpolation; 
        if x' is lower resolution, this is done by integrating over the relevant pixels

    Required Inputs:

        :param x1: x-axis values for original function
        :param y1: y-axis values for original function
        :param x2: x-axis values for output function

    Optional Inputs:

        :param nsamp: Number of samples for stepwise integration if going from high resolution to low resolution

    Output:

        y-axis values for resulting remapped function

    :Example:

    >>> # a coarse way of downsampling spectrum
    >>> import splat, numpy
    >>> sp = splat.Spectrum(file='high_resolution_spectrum.fits')
    >>> w_low = numpy.linspace(numpy.min(sp.wave.value),numpy.max(sp.wave.value),len(sp.wave.value)/10.)
    >>> f_low = splat.integralResample(sp.wave.value,sp.flux.value,w_low)
    >>> n_low = splat.integralResample(sp.wave.value,sp.noise.value,w_low)
    >>> sp.wave = w_low*sp.wave.unit
    >>> sp.flux = f_low*sp.flux.unit
    >>> sp.noise = n_low*sp.noise.unit
    '''

# check inputs
    if x2[0] < x1[0] or x2[-1] > x1[-1]: 
        raise ValueError('\nOutput x range {} to {} must be within input x range {} to {}'.format(x2[0],x2[-1],x1[0],x1[-1]))

# low resolution -> high resolution: interpolation
    if len(x1) <= len(x2): 
        f = interp1d(x1,y1,bounds_error=False,fill_value=0.)
        y2 = f(x2)

# high resolution -> low resolution: integrate
    else: 

# slow flux-preserving method
        if method == 'splat':
            xs = [numpy.max([x1[0],x2[0]-0.5*(x2[1]-x2[0])])]
            for i in range(len(x2)-1): xs.append(x2[i]+0.5*(x2[i+1]-x2[i]))
            xs.append(numpy.min([x2[-1]+0.5*(x2[-1]-x2[-2]),x1[-1]]))

# integral loop
            y2 = []
            for i in range(len(x2)):
                dx = numpy.linspace(xs[i],xs[i+1],nsamp)
                y2.append(trapz(f(dx),x=dx)/trapz(numpy.ones(nsamp),x=dx))

# fast method
        elif method == 'fast':
            baseline = numpy.polynomial.Polynomial.fit(x1, y1, 4)
            ip       = InterpolatedUnivariateSpline(x1, y1/baseline(x1), k=3)
            y2       = baseline(x2)*ip(x2)

    return y2



def integralResample_OLD(xh, yh, xl, nsamp=100):
    '''
    :Purpose: A 1D integral smoothing and resampling function that attempts to preserve total flux. Uses
    scipy.interpolate.interp1d and scipy.integrate.trapz to perform piece-wise integration

    Required Inputs:

    :param xh: x-axis values for "high resolution" data
    :param yh: y-axis values for "high resolution" data
    :param xl: x-axis values for resulting "low resolution" data, must be contained within high resolution and have fewer values

    Optional Inputs:

    :param nsamp: Number of samples for stepwise integration

    Output:

    y-axis values for resulting "low resolution" data

    :Example:
    >>> # a coarse way of downsampling spectrum
    >>> import splat, numpy
    >>> sp = splat.Spectrum(file='high_resolution_spectrum.fits')
    >>> w_low = numpy.linspace(numpy.min(sp.wave.value),numpy.max(sp.wave.value),len(sp.wave.value)/10.)
    >>> f_low = splat.integralResample(sp.wave.value,sp.flux.value,w_low)
    >>> n_low = splat.integralResample(sp.wave.value,sp.noise.value,w_low)
    >>> sp.wave = w_low*sp.wave.unit
    >>> sp.flux = f_low*sp.flux.unit
    >>> sp.noise = n_low*sp.noise.unit
    '''
# check inputs
    if xl[0] < xh[0] or xl[-1] > xh[-1]: raise ValueError('\nLow resolution x range {} to {} must be within high resolution x range {} to {}'.format(xl[0],xl[-1],xh[0],xh[-1]))
    if len(xl) > len(xh): raise ValueError('\nTarget x-axis must be lower resolution than original x-axis')

# set up samples
    xs = [numpy.max([xh[0],xl[0]-0.5*(xl[1]-xl[0])])]
    for i in range(len(xl)-1): xs.append(xl[i]+0.5*(xl[i+1]-xl[i]))
    xs.append(numpy.min([xl[-1]+0.5*(xl[-1]-xl[-2]),xh[-1]]))

    f = interp1d(xh,yh)

# integral loop
    ys = []
    for i in range(len(xl)):
        dx = numpy.linspace(xs[i],xs[i+1],nsamp)
        ys.append(trapz(f(dx),x=dx)/trapz(numpy.ones(nsamp),x=dx))
#    plt.plot(xh,yh,color='k')
#    plt.plot(xl,ys,color='r')

    return ys


def integralResample(xh, yh, xl, nsamp=100,method='fast'):
    '''
    :Purpose: A 1D integral smoothing and resampling function that attempts to preserve total flux. Uses
    scipy.interpolate.interp1d and scipy.integrate.trapz to perform piece-wise integration

    Required Inputs:

    :param xh: x-axis values for "high resolution" data
    :param yh: y-axis values for "high resolution" data
    :param xl: x-axis values for resulting "low resolution" data, must be contained within high resolution and have fewer values

    Optional Inputs:

    :param nsamp: Number of samples for stepwise integration

    Output:

    y-axis values for resulting "low resolution" data

    :Example:
    >>> # a coarse way of downsampling spectrum
    >>> import splat, numpy
    >>> sp = splat.Spectrum(file='high_resolution_spectrum.fits')
    >>> w_low = numpy.linspace(numpy.min(sp.wave.value),numpy.max(sp.wave.value),len(sp.wave.value)/10.)
    >>> f_low = splat.integralResample(sp.wave.value,sp.flux.value,w_low)
    >>> n_low = splat.integralResample(sp.wave.value,sp.noise.value,w_low)
    >>> sp.wave = w_low*sp.wave.unit
    >>> sp.flux = f_low*sp.flux.unit
    >>> sp.noise = n_low*sp.noise.unit
    '''
# check inputs
    if xl[0] < xh[0] or xl[-1] > xh[-1]: raise ValueError('\nLow resolution x range {} to {} must be within high resolution x range {} to {}'.format(xl[0],xl[-1],xh[0],xh[-1]))
    if len(xl) > len(xh): raise ValueError('\nTarget x-axis must be lower resolution than original x-axis')

# set up samples
    if method == 'splat':
        xs = [numpy.max([xh[0],xl[0]-0.5*(xl[1]-xl[0])])]
        for i in range(len(xl)-1): xs.append(xl[i]+0.5*(xl[i+1]-xl[i]))
        xs.append(numpy.min([xl[-1]+0.5*(xl[-1]-xl[-2]),xh[-1]]))

        f = interp1d(xh,yh,bounds_error=False,fill_value=0.)

# integral loop
        ys = []
        for i in range(len(xl)):
            dx = numpy.linspace(xs[i],xs[i+1],nsamp)
            ys.append(trapz(f(dx),x=dx)/trapz(numpy.ones(nsamp),x=dx))
#    plt.plot(xh,yh,color='k')
#    plt.plot(xl,ys,color='r')

# NOTE: THIS METHOD TENDS TO PRODUCE ALL NAN ARRAYS - NEED TO RETHINK THIS
    elif method == 'fast':
#        print(xh,yh)
        baseline = numpy.polynomial.Polynomial.fit(xh, yh, 4)
        ip       = InterpolatedUnivariateSpline(xh, yh/baseline(xh), k=3)
        ys       = baseline(xl)*ip(xl)

    return ys





def randomSphereAngles(num,longitude_range=[0,2*numpy.pi],latitude_range=[-0.5*numpy.pi,0.5*numpy.pi],exclude_longitude_range=[],exclude_latitude_range=[],degrees=False,**kwargs):
    '''
    :Purpose: 

        Draw a set of angles from a uniform spherical distribution, with areal inclusion and exclusion constraints.
        Note that latitude range is assumed to run from -pi/2 to +pi/2

    :Required Input:

        :param num: number of points to draw

    :Optional Input:

        :param: longitude_range = [0,2pi]: range over longitude to draw values
        :param: latitude_range = [-pi,+pi]: range over latitude to draw values
        :param: exclude_longitude_range = []: range of longitudes to exclude values
        :param: exclude_latitude_range = []: range of latitudes to exclude values
        :param: degrees = False: by default, radians are assumed; set to True to convert to degrees (also checks if inclusion/exclusion ranges are in degrees)

    :Output:

        2 arrays of longitudes and latitudes drawn uniformly over select area

    :Example:
        >>> import splat
        >>> splat.randomSphereAngles(10)
            (array([ 2.52679013,  0.85193769,  5.98514797,  0.89943465,  5.36310536,
             5.34344768,  0.01743906,  4.93856229,  0.06508084,  0.5517308 ]),
            array([-0.53399501,  0.04208564,  0.03089855, -0.60445954,  0.55800151,
             0.80119146, -0.19318715,  0.76230148, -0.5935969 , -0.65839849]))
        >>> splat.randomSphereAngles(10,latitude_range=[-10,10],degrees=True)
            (array([  28.55709202,  297.34760719,  152.79525894,   71.08745583,
                     153.56948338,   80.68486463,    7.75479896,  100.8408509 ,
                     356.63091754,   66.16572906]),
             array([ 0.6747939 , -1.00316889, -2.26239023,  9.27397372, -8.96797181,
                     7.34796163, -1.93175289,  3.07888912,  0.69826684, -5.08428339]))
    '''
# check inputs - convert to radians if necessary   
    if degrees==True and numpy.max(numpy.absolute(longitude_range)) > 2.*numpy.pi:
        longitude_range = [l*numpy.pi/180. for l in longitude_range]
    if degrees==True and numpy.max(numpy.absolute(latitude_range)) > numpy.pi:
        latitude_range = [l*numpy.pi/180. for l in latitude_range]

# longitude - uniformly distributed
    longitude = numpy.random.uniform(0,1,num)*(longitude_range[1]-longitude_range[0])+longitude_range[0]

# latitude - distributed by P(phi) = 1/2 cos(phi) for -pi/2 < phi < pi/2
    x = numpy.linspace(latitude_range[0],latitude_range[1],num)
    cdf = 0.5*(numpy.sin(x)+1.)
    cdf = cdf-numpy.nanmin(cdf)
    cdf = cdf/numpy.nanmax(cdf)
    f = interp1d(cdf,x)
    latitude = f(numpy.random.uniform(0,1,num))

# exclude ranges specified
    if len(exclude_longitude_range) > 0:
        if degrees==True and numpy.max(numpy.absolute(exclude_longitude_range)) > 2.*numpy.pi:
            exclude_longitude_range = [l*numpy.pi/180. for l in exclude_longitude_range]
        longex = longitude[longitude<numpy.nanmin(exclude_longitude_range)]
        longex = numpy.concatenate((longex,longitude[longitude>numpy.nanmax(exclude_longitude_range)]))
        while len(longex) < num:
            longitude = numpy.random.uniform(0,1,num)*(longitude_range[1]-longitude_range[0])+longitude_range[0]
            longex = numpy.concatenate((longex,longitude[longitude<numpy.nanmin(exclude_longitude_range)]))
            longex = numpy.concatenate((longex,longitude[longitude>numpy.nanmax(exclude_longitude_range)]))
        longitude = longex[:num]

    if len(exclude_latitude_range) > 0:
        if degrees==True and numpy.max(numpy.absolute(exclude_latitude_range)) > numpy.pi:
            exclude_latitude_range = [l*numpy.pi/180. for l in exclude_latitude_range]
        latex = latitude[latitude<numpy.nanmin(exclude_latitude_range)]
        latex = numpy.concatenate((latex,latitude[latitude>numpy.nanmax(exclude_latitude_range)]))
        while len(latex) < num:
            x = numpy.linspace(latitude_range[0],latitude_range[1],num)
            cdf = 0.5*(numpy.sin(x)+1.)
            cdf = cdf-numpy.nanmin(cdf)
            cdf = cdf/numpy.nanmax(cdf)
            f = interp1d(cdf,x)
            latitude = f(numpy.random.uniform(0,1,num))
            latex = numpy.concatenate((latex,latitude[latitude<numpy.nanmin(exclude_latitude_range)]))
            latex = numpy.concatenate((latex,latitude[latitude>numpy.nanmax(exclude_latitude_range)]))
        latitude = latex[:num]        

# outputs; convert to degrees if desired    
    if degrees==True:
        latitude = latitude*180./numpy.pi
        longitude = longitude*180./numpy.pi
    return longitude, latitude
        

def weightedMeanVar(vals, winp, *args, **kwargs):
    '''
    :Purpose: 

        Computes weighted mean of an array of values through various methods. Returns weighted mean and weighted uncertainty.

    :Required Inputs:

        :param **vals**: array of values
        :param **winp**: array of weights associated with ``vals``

    :Optional Inputs:

        :param **method**: type of weighting to be used. Options include:

            - *default*: (default) ``winp`` is taken to be actual weighting values
            - *uncertainty*: uncertainty weighting, where ``winp`` is the uncertainties of ``vals``
            - *ftest*: ftest weighting, where ``winp`` is the chi squared values of ``vals``

        :param **weight_minimum**: minimum possible weight value (default = 0.)
        :param **dof**: effective degrees of freedom (default = len(vals) - 1)

        .. note:: When using ``ftest`` method, extra ``dof`` value is required

    :Output:

        Weighted mean and uncertainty

    :Example:
        >>> import splat
        >>> splat.weightedMeanVar([3.52, 5.88, 9.03], [0.65, 0.23, 0.19])
            (5.0057009345794379, 4.3809422657000594)
        >>> splat.weightedMeanVar([3.52, 5.88, 9.03], [1.24, 2.09, 2.29], method = 'uncertainty')
            (5.0069199363443841, 4.3914329968409946)
    '''

    method = kwargs.get('method','')
    minwt = kwargs.get('weight_minimum',0.)
    dof = kwargs.get('dof',len(vals)-1)
    if (numpy.nansum(winp) <= 0.):
        weights = numpy.ones(len(vals))
    if isinstance(winp,u.quantity.Quantity):
        winput = winp.value
    else:
        winput = copy.deepcopy(winp)

# uncertainty weighting: input is unceratinties
    if (method == 'uncertainty'):
        weights = [w**(-2) for w in winput]
# ftest weighting: input is chisq values, extra dof value is required
    elif (method == 'ftest'):
# fix issue of chi^2 = 0
        minchi = numpy.nanmin(winput)
        weights = numpy.array([stats.f.pdf(w/minchi,dof,dof) for w in winput])
# just use the input as the weights
    else:
        weights = [w for w in winput]

    weights = weights/numpy.nanmax(weights)
    weights[numpy.where(weights < minwt)] = 0.
    mn = numpy.nansum(vals*weights)/numpy.nansum(weights)
    var = numpy.nansum(weights*(vals-mn)**2)/numpy.nansum(weights)
    if (method == 'uncertainty'):
        var+=numpy.nansum([w**2 for w in winput])/(len(winput)**2)

    return mn,numpy.sqrt(var)


#####################################################
###############   DATABASE HELPERS   ################
#####################################################

def checkDBCoordinates(db,designation_keyword='DESIGNATION',ra_keyword='RA',dec_keyword='DEC',shortname_keyword='SHORTNAME'):
# designation -> ra, dec
    if designation_keyword in list(db.keys()):
        if ra_keyword not in list(db.keys()) or dec_keyword not in list(db.keys()):
            coord = [designationToCoordinate(d) for d in db[designation_keyword]]
            db[ra_keyword] = [c.ra.deg for c in coord]
            db[dec_keyword] = [c.dec.deg for c in coord]
# ra,dec -> designation
    else:
        if ra_keyword not in list(db.keys()) or dec_keyword not in list(db.keys()):
            print('Warning: cannot populate designation column {} without RA column {} and DEC column {}'.format(designation_keyword,ra_keyword,dec_keyword))
        else:
            db[designation_keyword] = [coordinateToDesignation([db[ra_keyword].iloc[i],db[dec_keyword].iloc[i]]) for i in range(len(db))]
# designation -> shortname
    if designation_keyword in list(db.keys()):
        if shortname_keyword not in list(db.keys()):
            db[shortname_keyword] = [designationToShortName(d) for d in db[designation_keyword]]
    return db


def keySearch(keys, key_name='KEY', db=pandas.DataFrame(), verbose=True):
    '''
    Purpose
    -------

    General purpose function that searches a pandas database for values of keys that match in 
    key_name column. Use by `keySource()`_ and `keySpectrum()`_.

    Parameters
    ----------

    keys : int or list
        integer or list of integers corresponding to keys

    key_name = 'KEY' : string [optional]
        name of column to search 

    db = blank pandas DataFrame : pandas DataFrame [optional]
        pandas data Frame to search; code checks that it contains key_name column

    verbose = True : boolean [optional]
        set to True to have program return verbose output

    Outputs
    -------

    pandas DataFrame containing the rows that match input keys, or empty pandas DataFrame

    Example
    -------

        TBD

    Dependencies
    ------------
        None

    '''
# check db is a pandas DataFrame, or try to convert
    if isinstance(db,pandas.DataFrame) == False:
        try: dbconv = pandas.DataFrame(db)
        except: pass
    else: dbconv = copy.deepcopy(db)
    if isinstance(dbconv,pandas.DataFrame) == False:
        if verbose==True: print('Passed database is not a pandas dataframe')
        return pandas.DataFrame()

# check key is in data frame
    if key_name not in list(dbconv.columns):
        raise ValueError('Cannot find key column {} in dataframe'.format(key_name))

# vectorize and make integer
    if isinstance(keys,list) == False: keys = [keys]
    keys = [int(k) for k in keys]

# search dataframe
    sdb = dbconv[[x in keys for x in dbconv[key_name]]]
    if len(sdb) == 0 and verbose==True: print('No sources found with key(s) = {}'.format(*keys))
    return sdb
    

#####################################################
################   CODE MANAGEMENT   ################
#####################################################
#
# Note that all of these should have a checkAccess() flag
#
#####################################################

def codeStats():
    if checkAccess() == False:
        raise ValueError('You do not have sufficient access to run this program\n')

# library statistics - # of total/public spectra, # of total/public sources, # of source papers for public data
    sall = splat.searchLibrary()
    print('Total number of spectra = {} of {} sources'.format(len(sall),len(numpy.unique(numpy.array(sall['SOURCE_KEY'])))))
    s = splat.searchLibrary(public=True)
    print('Total number of public spectra = {} of {} sources'.format(len(s),len(numpy.unique(numpy.array(s['SOURCE_KEY'])))))

# data citations 
    pubs = numpy.unique(numpy.array(sall['BIBCODE']))
    print('Total number of citations for all spectra = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(s['BIBCODE']))
    print('Total number of citations for public spectra = {}'.format(len(pubs)))
    cites = []
    cites_html = []
    for p in pubs:
        try:
            cites_html.append('<li>{} [<a href="{}">NASA ADS</a>]'.format(splat.citations.longRef(str(p)),splat.citations.citeURL(str(p))))
            cites.append('{}'.format(splat.citations.longRef(str(p))))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    cites.sort()
    # with open(DOCS_FOLDER+'_static/citation_list.txt', 'w') as f:
    #     f.write('Data references in SPLAT:\n')
    #     for c in cites:
    #         f.write('{}\n'.format(c))
    # cites_html.sort()
    # with open(DOCS_FOLDER+'_static/citation_list.html', 'w') as f:
    #     f.write('<ul>\n')
    #     for c in cites_html:
    #         f.write('\t{}\n'.format(c))
    #     f.write('</ul>\n')

# source citations 
    pubs = numpy.unique(numpy.array(sall['DISCOVERY_REFERENCE'].replace(numpy.nan,'')))
    print('Total number of citations for all sources = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['OPT_TYPE_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all optical spectral types = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['NIR_TYPE_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all NIR spectral types = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['LIT_TYPE_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all literature spectral types = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['GRAVITY_CLASS_OPTICAL_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all optical gravity types = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['GRAVITY_CLASS_NIR_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all NIR gravity types = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['CLUSTER_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all cluster associations = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['BINARY_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all binary associations = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['SBINARY_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all spectral binary associations = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['COMPANION_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all companion associations = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['SIMBAD_SPT_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all SIMBAD SpTs = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['PARALLEX_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all parallaxes = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['MU_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all proper motions = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['RV_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all RVs = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))
    pubs = numpy.unique(numpy.array(sall['VSINI_REF'].replace(numpy.nan,'')))
    print('Total number of citations for all vsini values = {}'.format(len(pubs)))
    for p in pubs:
        try:
            x = splat.citations.longRef(str(p))
        except:
            print('\tWarning: no bibtex information for citation {}'.format(p))

# histogram of spectral types - all spectra
    sptrng = [16,40]
    xticks = range(sptrng[0],sptrng[1])
    labels = [splat.typeToNum(x) for x in range(sptrng[0],sptrng[1])]
    for i in range(2):
        if i == 0:
            s1 = sall[sall['OBJECT_TYPE'] == 'VLM']
            fname = 'all'
        else:
            s1 = s[s['OBJECT_TYPE'] == 'VLM']
            fname = 'published'
        spex_spts = []
        opt_spts = []
        nir_spts = []
        spts = []
        for i,x in enumerate(s1['SPEX_TYPE']):
            spt = -99.
            if splat.isNumber(splat.typeToNum(x)): 
                sspt = splat.typeToNum(x)
                spex_spts.append(sspt)
                spt = copy.deepcopy(sspt)

            nspt = splat.typeToNum(s1['NIR_TYPE'].iloc[i])
            if splat.isNumber(nspt):
                nir_spts.append(spt)
                if nspt > 28.: spt = copy.deepcopy(nspt)
            ospt = splat.typeToNum(s1['OPT_TYPE'].iloc[i])
            if splat.isNumber(ospt):
                opt_spts.append(spt)
                if ospt < 29.: spt = copy.deepcopy(ospt)
            if spt > 0: spts.append(spt)
    # SpeX type
        sptarr = numpy.array(spex_spts)
        plt.figure(figsize=(14,6))
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= sptrng[0],sptarr < 20))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='green', alpha=0.75)
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= 20,sptarr < 30))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='red', alpha=0.75)
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= 30,sptarr < sptrng[1]))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='b', alpha=0.75)
        plt.xticks(xticks,labels)
        plt.xlabel('SpeX Spectral Type')
        plt.ylabel('log10 Number')
        plt.xlim([sptrng[0]-0.5,sptrng[1]+0.5])
        plt.legend(['M dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= sptrng[0],sptarr < 20))])),'L dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= 20,sptarr < 30))])),'T dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= 30,sptarr < sptrng[1]))]))])
#        plt.savefig(DOCS_FOLDER+'_images/spt_spex_distribution_{}.png'.format(fname))
        plt.show()
        plt.clf()
    # Optical type
        sptarr = numpy.array(opt_spts)
        plt.figure(figsize=(14,6))
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= sptrng[0],sptarr < 20))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='green', alpha=0.75)
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= 20,sptarr < 30))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='red', alpha=0.75)
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= 30,sptarr < sptrng[1]))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='b', alpha=0.75)
        plt.xticks(xticks,labels)
        plt.xlabel('Optical Spectral Type')
        plt.ylabel('log10 Number')
        plt.xlim([sptrng[0]-0.5,sptrng[1]+0.5])
        plt.legend(['M dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= sptrng[0],sptarr < 20))])),'L dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= 20,sptarr < 30))])),'T dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= 30,sptarr < sptrng[1]))]))])
#        plt.savefig(DOCS_FOLDER+'_images/spt_optical_distribution_{}.png'.format(fname))
        plt.show()
        plt.clf()
    # NIR type
        sptarr = numpy.array(nir_spts)
        plt.figure(figsize=(14,6))
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= sptrng[0],sptarr < 20))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='green', alpha=0.75)
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= 20,sptarr < 30))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='red', alpha=0.75)
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= 30,sptarr < sptrng[1]))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='b', alpha=0.75)
        plt.xticks(xticks,labels)
        plt.xlabel('NIR Spectral Type')
        plt.ylabel('log10 Number')
        plt.xlim([sptrng[0]-0.5,sptrng[1]+0.5])
        plt.legend(['M dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= sptrng[0],sptarr < 20))])),'L dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= 20,sptarr < 30))])),'T dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= 30,sptarr < sptrng[1]))]))])
#        plt.savefig(DOCS_FOLDER+'_images/spt_nir_distribution_{}.png'.format(fname))
        plt.show()
        plt.clf()
    # Adopted type
        sptarr = numpy.array(spts)
        plt.figure(figsize=(14,6))
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= sptrng[0],sptarr < 20))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='green', alpha=0.75)
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= 20,sptarr < 30))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='red', alpha=0.75)
        n, bins, patches = plt.hist(sptarr[numpy.where(numpy.logical_and(sptarr >= 30,sptarr < sptrng[1]))], bins=len(range(sptrng[0],sptrng[1])), log=True, range=sptrng, facecolor='b', alpha=0.75)
        plt.xticks(xticks,labels)
        plt.xlabel('Adopted Spectral Type')
        plt.ylabel('log10 Number')
        plt.xlim([sptrng[0]-0.5,sptrng[1]+0.5])
        plt.legend(['M dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= sptrng[0],sptarr < 20))])),'L dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= 20,sptarr < 30))])),'T dwarfs ({} sources)'.format(len(sptarr[numpy.where(numpy.logical_and(sptarr >= 30,sptarr < sptrng[1]))]))])
#        plt.savefig(DOCS_FOLDER+'_images/spt_adopted_distribution_{}.png'.format(fname))
        plt.show()
        plt.clf()

# histogram of S/N

# map sources on sky
    raref = Angle(numpy.linspace(0,359.,360)*u.degree)
    raref.wrap_at(180.*u.degree)
    ra = Angle(list(sall['RA'])*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(list(sall['DEC'])*u.degree)
    rap = Angle(list(s['RA'])*u.degree)
    rap = rap.wrap_at(180*u.degree)
    decp = Angle(list(s['DEC'])*u.degree)
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    p1 = ax.scatter(ra.radian, dec.radian,color='r',alpha=0.5,s=10)
    p2 = ax.scatter(rap.radian, decp.radian,color='k',alpha=0.5, s=5)
#    ur = ax.plot(raref.radian,Angle([67.]*len(raref)*u.degree).radian,'k--')
#    ur = ax.plot(raref.radian,Angle([-50.]*len(raref)*u.degree).radian,'k--')
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    ax.grid(True)
#    ef = matplotlib.patheffects.withStroke(foreground="w", linewidth=3)
#    axis = ax.axis['lat=0']
#    axis.major_ticklabels.set_path_effects([ef])
 #   axis.label.set_path_effects([ef])
    plt.legend([p1,p2],['All Sources ({})'.format(len(sall)),'Published Sources ({})'.format(len(s))],bbox_to_anchor=(1, 1),bbox_transform=plt.gcf().transFigure)
#    fig.savefig(DOCS_FOLDER+'_images/map_all.png')
    plt.show()
    fig.clf()

# map sources on based on spectral type
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    sm = splat.searchLibrary(spt_range=[10,19.9],spt_type='SPEX')
    sm = sm[sm['OBJECT_TYPE'] == 'VLM']
    sl = splat.searchLibrary(spt_range=[20,29.9],spt_type='SPEX')
    sl = sl[sl['OBJECT_TYPE'] == 'VLM']
    st = splat.searchLibrary(spt_range=[30,39.9],spt_type='SPEX')
    st = st[st['OBJECT_TYPE'] == 'VLM']
    ra = Angle(list(sm['RA'])*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(list(sm['DEC'])*u.degree)
    p1 = ax.scatter(ra.radian, dec.radian,color='k',alpha=0.5,s=10)
    ra = Angle(list(sl['RA'])*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(list(sl['DEC'])*u.degree)
    p2 = ax.scatter(ra.radian, dec.radian,color='r',alpha=0.5,s=10)
    ra = Angle(list(st['RA'])*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(list(st['DEC'])*u.degree)
    p3 = ax.scatter(ra.radian, dec.radian,color='b',alpha=0.5,s=10)
    plt.legend([p1,p2,p3],['M dwarfs ({})'.format(len(sm)),'L dwarfs ({})'.format(len(sl)),'T dwarfs ({})'.format(len(st))],bbox_to_anchor=(1, 1),bbox_transform=plt.gcf().transFigure)
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    ax.grid(True)
#    fig.savefig(DOCS_FOLDER+'_images/map_byspt.png')
    plt.show()
    fig.clf()

# map sources on based on young or field
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    sy = splat.searchLibrary(young=True)
    sy = sy[sy['OBJECT_TYPE'] == 'VLM']
    so = splat.searchLibrary()
    so = so[so['OBJECT_TYPE'] == 'VLM']
    ra = Angle(list(so['RA'])*u.degree)
#    ra = ra.wrap_at(180*u.degree)
#    dec = Angle(so['DEC'].filled(numpy.nan)*u.degree)
#    p1 = ax.scatter(ra.radian, dec.radian,color='k',alpha=0.1,s=5)
    ra = Angle(list(sy['RA'])*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(list(sy['DEC'])*u.degree)
    p1 = ax.scatter(ra.radian, dec.radian,color='r',alpha=0.5,s=10)
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    ax.grid(True)
    plt.legend([p1],['Young Sources ({})'.format(len(sy))],bbox_to_anchor=(1, 1),bbox_transform=plt.gcf().transFigure)
#    fig.savefig(DOCS_FOLDER+'_images/map_young.png')
    plt.show()
    fig.clf()

# pie chart of spectrum types
    ot = numpy.unique(numpy.array(sall['OBJECT_TYPE']))
    otypes = 'STAR','GIANT','WD','GALAXY','OTHER'
    sizes = [len(sall[sall['OBJECT_TYPE']==o]) for o in otypes]
    explode = (0.1,0,0,0,0)

    fig, ax = plt.subplots()
    ax.pie(sizes, explode=explode, labels=otypes, autopct='%1.1f%%',
        shadow=True, startangle=90, pctdistance = 0.7)
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
#    fig.savefig(DOCS_FOLDER+'_images/object_othertypes.png')
    plt.show()


def about():
    '''
    Gives basic information about SPLAT code
    '''
    print('\nSPLAT (SpeX Prism Library and Analysis Toolkit)')
    print('\nSPLAT was created by members of the Cool Star Lab:')
    for a in splat.AUTHORS: print('\t'+a)
    print('\nFunding for SPLAT was provided by the National Aeronautics and Space Administration under grant NNX15AI75G')
    print('\nSPLAT can be downloaded at '+splat.GITHUB_URL)
    print('Documentation can be found at '+splat.DOCUMENTATION_URL)
    print('\nIf you use SPLAT, please cite the software paper '+splat.citations.shortRef(splat.BIBCODE))
    print('\nIf you use any of the data or models in SPLAT, you must cite the original references for these')

    return





