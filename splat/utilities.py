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
import re
import requests
import string
import sys

# imports - external
from astropy.coordinates import Angle,SkyCoord      # coordinate conversion
from astropy import units as u            # standard units
import matplotlib.pyplot as plt
import matplotlib.patheffects
import numpy
from scipy import stats

# code constants
from .initialize import *
import splat

# Python 2->3 fix for input
try: input=raw_input
except NameError: pass

# change the command prompt
sys.ps1 = 'splat util> '


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



def checkAccess(**kwargs):
    '''
    :Purpose: Checks if user has access to unpublished spectra in SPLAT library.
    :Example:
       >>> import splat
       >>> print spl.checkAccess()
       True
    :Note: Must have the file .splat_access in your home directory with the correct passcode to use.
    '''
    result = False

    try:
	    home = os.path.expanduser("~")
	    if home == None:
	        home = './'
	    bcode = requests.get(SPLAT_URL+ACCESS_FILE).content
	    lcode = base64.b64encode(open(home+'/'+ACCESS_FILE,'r').read().encode())
	    if (bcode in lcode):        # changed to partial because of EOL variations
	        result = True
    except:
        result = False

    if (kwargs.get('verbose',False) == True):
        if result == True:
            print('You have full access to all SPLAT data')
        else:
            print('You have access only to published data')
    return result


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
    if not os.path.exists(inputfile):
        if not os.path.exists(SPLAT_PATH+inputfile):
            return ''
        else:
            return SPLAT_PATH+inputfile
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
    if (len(args) != 0):
        if 'http://' in args[0]:
            if requests.get(args[0]).status_code == requests.codes.ok:
                return args[0]
            return False
        else:
            if requests.get(SPLAT_URL+args[0]).status_code == requests.codes.ok:
                return SPLAT_URL+args[0]
            return False
    else:
        return requests.get(SPLAT_URL).status_code == requests.codes.ok



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


def checkInstrument(instrument):
    '''

    Purpose: 
        Checks that an instrument name is one of the available instruments, including a check of alternate names

    Required Inputs:
        :param: instrument: A string containing the instrument name to be checked. This should be one of the instruments in the global parameter splat.initialize.INSTRUMENTS
        
    Optional Inputs:
        None

    Output:
        A string containing SPLAT's default name for a given model set, or False if that model set is not present

    Example:

    >>> import splat
    >>> print(splat._checkModelName('burrows'))
        burrows06
    >>> print(splat._checkModelName('allard'))
        BTSettl2008
    >>> print(splat._checkModelName('somethingelse'))
        False
    '''
    output = False
    if not isinstance(instrument,str):
        return output
    for k in list(INSTRUMENTS.keys()):
        if instrument.lower()==k.lower() or instrument.lower() in INSTRUMENTS[k]['altnames']:
            output = k
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
    >>> print(splat._checkModelName('burrows'))
        burrows06
    >>> print(splat._checkModelName('allard'))
        BTSettl2008
    >>> print(splat._checkModelName('somethingelse'))
        False
    '''
    output = False
    if not isinstance(model,str):
        return output
    for k in list(SPECTRAL_MODELS.keys()):
        if model.lower()==k.lower() or model.lower() in SPECTRAL_MODELS[k]['altnames']:
            output = k
    return output



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
    if not isinstance(d,str): d = str(d)

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
        print('\nCould not determine format of input date; please provide a format string\n')
        return d

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


def coordinateToDesignation(c):
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
    if isinstance(c,SkyCoord):
        cc = c
    else:
        cc = properCoordinates(c)
# input is [RA,Dec] pair in degrees
    if sys.version_info.major == 2:
        return string.replace('J{0}{1}'.format(cc.ra.to_string(unit=u.hour, sep='', precision=2, pad=True), \
        cc.dec.to_string(unit=u.degree, sep='', precision=1, alwayssign=True, pad=True)),'.','')
    else:
        return str.replace('J{0}{1}'.format(cc.ra.to_string(unit=u.hour, sep='', precision=2, pad=True), \
        cc.dec.to_string(unit=u.degree, sep='', precision=1, alwayssign=True, pad=True)),'.','')



def designationToCoordinate(value, **kwargs):
    '''
    :Purpose: Convert a designation srtring into a RA, Dec tuple or ICRS SkyCoord objects (default)

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
    icrsflag = kwargs.get('icrs',True)

    a = re.sub('[j.:hms]','',value.lower())
    fact = 1.
    spl = a.split('+')
    if len(spl) == 1:
        spl = a.split('-')
        fact = -1.
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
    if icrsflag:
        return SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    else:
        return [ra,dec]


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



def properCoordinates(c,**kwargs):
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
        output = SkyCoord(c[0]*u.deg,c[1]*u.deg,frame='icrs')
# input is sexigessimal string
    elif isinstance(c,str):
        if c[0] == 'J':
            output = designationToCoordinate(c,**kwargs)
        else:
            output = SkyCoord(c,frame='icrs', unit=(u.hourangle, u.deg))
    else:
        raise ValueError('\nCould not parse input format\n\n')

    if kwargs.get('distance',False) != False:
        d = copy.deepcopy(kwargs['distance'])
        if not isinstance(d,u.quantity.Quantity):
            d*=u.pc
#        try:
        d.to(u.pc)
        print(output.distance)
        output = SkyCoord(output,distance = d)
#        except:
#            print('\nWarning: could not integrate distance {} into coordinate'.format(distance))

    return output



def typeToNum(inp, **kwargs):
    '''
    :Purpose: Converts between string and numeric spectral types, and vise versa.
    :param input: Spectral type to convert. Can convert a number or a string from 0 (K0) and 49.0 (Y9).
    :param error: magnitude of uncertainty. ':' for uncertainty > 1 and '::' for uncertainty > 2.
    :type error: optional, default = ''
    :param uncertainty: uncertainty of spectral type
    :type uncertainty: optional, default = 0
    :param subclass: subclass of object. Options include:

        - *sd*: object is a subdwarf
        - *esd*: object is an extreme subdwarf
        - *usd*: object is an ultra subdwarf

    :type subclass: optional, default = ''
    :param lumclass: luminosity class of object represented by roman numerals
    :type lumclass: optional, default = ''
    :param ageclass: age class of object
    :type ageclass: optional, default = ''
    :param colorclass: color class of object
    :type colorclass: optional, default = ''
    :param peculiar: if object is peculiar or not
    :type peculiar: optional, default = False

    .. not too sure how colorclass and ageclass work

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
    error = kwargs.get('error','')
    unc = kwargs.get('uncertainty',0.)
    subclass = kwargs.get('subclass','')
    lumclass = kwargs.get('lumclass','')
    ageclass = kwargs.get('ageclass','')
    colorclass = kwargs.get('colorclass','')
    peculiar = kwargs.get('peculiar',False)
    spletter = 'KMLTY'
    verbose = kwargs.get('verbose',False)

# convert input into an array
    output = []
    var = copy.deepcopy(inp)

    if not isinstance(var,list): var = [var]
    if not isinstance(error,list): error = [error]
    if not isinstance(unc,list): unc = [unc]
    if not isinstance(subclass,list): subclass = [subclass]
    if not isinstance(lumclass,list): lumclass = [lumclass]
    if not isinstance(ageclass,list): ageclass = [ageclass]
    if not isinstance(colorclass,list): colorclass = [colorclass]
    if len(error) < len(var):
        for i in numpy.arange(len(var)-len(error)): error.append(error[-1])
    if len(unc) < len(var):
        for i in numpy.arange(len(var)-len(unc)): unc.append(unc[-1])

# number -> spectral type
    if (isNumber(var[0])):
        if len(subclass) < len(var):
            for i in numpy.arange(len(var)-len(subclass)): subclass.append(subclass[-1])
        if len(lumclass) < len(var):
            for i in numpy.arange(len(var)-len(lumclass)): lumclass.append(lumclass[-1])
        for i,l in enumerate(lumclass):
            if l != '': lumclass[i]=' '+lumclass[i]
        if len(ageclass) < len(var):
            for i in numpy.arange(len(var)-len(ageclass)): ageclass.append(ageclass[-1])
        if len(colorclass) < len(var):
            for i in numpy.arange(len(var)-len(colorclass)): colorclass.append(colorclass[-1])
        spind = numpy.array([int(abs(x/10.)) for x in var])
        spdec = numpy.array([numpy.around(x,1) for x in var])-spind*10.
        for i,v in enumerate(var):
            pstr = ''
            if (unc[i] > 1.):
                error[i] = ':'
            if (unc[i] > 2.):
                error[i] = '::'
            if (peculiar):
                pstr = 'p'
            if (0 <= spind[i] < len(spletter)):
                output.append(colorclass[i]+subclass[i]+spletter[spind[i]]+'{:3.1f}'.format(spdec[i])+ageclass[i]+lumclass[i]+pstr+error[i])
            else:
                if verbose: print('Spectral type number must be between 0 ({}0) and {} ({}9)'.format(spletter[0],len(spletter)*10.-1.,spletter[-1]))
                output.append('N/A')

# spectral type -> number
    elif isinstance(var[0],str):
        output = []
        for i,v in enumerate(var):
            if (sys.version_info.major == 2):
                v = string.split(v,sep='+')[0]    # remove +/- sides
                v = string.split(v,sep='-')[0]    # remove +/- sides
                v = string.split(v,sep='/')[0]    # remove / in spectral types
            else:
                v = v.split('+')[0]    # remove +/- sides
                v = v.split('-')[0]    # remove +/- sides
                v = v.split('/')[0]    # remove / in spectral types
            v = v.replace('...','').replace(' ','')

            sptype = re.findall('[{}]'.format(spletter),v.upper())
            outval = 0.
            if (len(sptype) == 1):
                outval = spletter.find(sptype[0])*10.
                spind = v.find(sptype[0])+1
                if (spind < len(v)):
                    if (v.find('.') < 0):
                        if (isNumber(v[spind])):
                            outval += float(v[spind])
                    else:
                        try:
                            outval += float(v[spind:spind+3])
                            spind = spind+3
                        except:
                            if verbose: print('\nProblem converting input type {} to a numeric type'.format(v))
                            outval = numpy.nan
                ytype = re.findall('[abcd]',v.split('p')[-1])
                if (len(ytype) == 1):
                    ageclass[i] = ytype[0]
                if (v.find('p') != -1):
                     peculiar = True
                if (v.find('sd') != -1):
                     subclass[i] = 'sd'
                if (v.find('esd') != -1):
                     subclass[i] = 'esd'
                if (v.find('usd') != -1):
                     subclass[i] = 'usd'
                if (v.count('I') > 0):
                     lumclass[i] = ''.join(re.findall('I',v))
                if (v.count(':') > 0):
                     error[i] = ''.join(re.findall(':',v))
                if (v[0] == 'b' or v[0] == 'r'):
                     colorclass[i] = v[0]
            if (len(sptype) != 1):
                if verbose: print('\nOnly spectral classes {} are handled by typeToNum'.format(spletter))
                outval = numpy.nan
            output.append(outval)

# none of the above - return the input
    else:
        if verbose: print('\nWarning: could not recognize format of spectral type {}\n'.format(inp))
        output = var
    if len(output) == 1 and not isinstance(inp,list): 
        return output[0]
    else:
        return output


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
        print(xsorted,cdf,wt[idx],type(xd),type(cdf))
        cdff = [float(c) for c in cdf]
        cdfn = [c/cdff[-1] for c in cdff]
        return numpy.interp(q, cdfn, xsorted).tolist()


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
    except ValueError:
        return False


def weightedMeanVar(vals, winp, *args, **kwargs):
    '''
    :Purpose: Computes weighted mean of an array of values through various methods. Returns weighted mean and weighted uncertainty.
    :param vals: array of values
    :param winp: array of weights associated with ``vals``
    :param method: input type of weights. Default is where ``winp`` is the actual weights of ``vals``. Options include:

        - *uncertainty*: uncertainty weighting, where ``winp`` is the uncertainties of ``vals``
        - *ftest*: ftest weighting, where ``winp`` is the chi squared values of ``vals``

    :type method: optional, default = ''
    :param weight_minimum: minimum possible weight value
    :type weight_minimum: optional, default = 0.
    :param dof: effective degrees of freedom
    :type dof: optional, default = len(vals) - 1

    .. note:: When using ``ftest`` method, extra ``dof`` value is required

    :Example:
        >>> import splat
        >>> print splat.weightedMeanVar([3.52, 5.88, 9.03], [0.65, 0.23, 0.19])
            (5.0057009345794379, 4.3809422657000594)
        >>> print splat.weightedMeanVar([3.52, 5.88, 9.03], [1.24, 2.09, 2.29], method = 'uncertainty')
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
################   CODE MANAGEMENT   ################
#####################################################
#
# Note that all of these should have a checkAccess() flag
#
#####################################################

def codeStats():
# library statistics - # of total/public spectra, # of total/public sources, # of source papers for public data
    if checkAccess() == False:
        raise ValueError('You do not have sufficient access to run this program\n')

    sall = splat.searchLibrary()
    print('Total number of spectra = {} of {} sources'.format(len(sall),len(numpy.unique(numpy.array(sall['SOURCE_KEY'])))))
    s = splat.searchLibrary(public=True)
    print('Total number of public spectra = {} of {} sources'.format(len(s),len(numpy.unique(numpy.array(s['SOURCE_KEY'])))))

# data citations 
    pubs = numpy.unique(numpy.array(s['DATA_REFERENCE']))
    print('Total number of citations for public spectra = {}'.format(len(pubs)))
    cites = []
    cites_html = []
    for p in pubs:
        try:
            cites_html.append('<li>{} [<a href="{}">NASA ADS</a>]'.format(splat.citations.longRef(str(p)),splat.citations.citeURL(str(p))))
            cites.append('{}'.format(splat.citations.longRef(str(p))))
        except:
            print('Warning: no bibtex information for citation {}'.format(p))
    cites.sort()
    with open(SPLAT_PATH+DOCS_FOLDER+'_static/citation_list.txt', 'w') as f:
        f.write('Data references in SPLAT:\n')
        for c in cites:
            f.write('{}\n'.format(c))
    cites_html.sort()
    with open(SPLAT_PATH+DOCS_FOLDER+'_static/citation_list.html', 'w') as f:
        f.write('<ul>\n')
        for c in cites_html:
            f.write('\t{}\n'.format(c))
        f.write('</ul>\n')

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
            nspt = splat.typeToNum(s1['NIR_TYPE'][i])
            if splat.isNumber(nspt):
                nir_spts.append(spt)
                if nspt > 28.: spt = copy.deepcopy(nspt)
            ospt = splat.typeToNum(s1['OPT_TYPE'][i])
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
        plt.savefig(SPLAT_PATH+DOCS_FOLDER+'_images/spt_spex_distribution_{}.png'.format(fname))
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
        plt.savefig(SPLAT_PATH+DOCS_FOLDER+'_images/spt_optical_distribution_{}.png'.format(fname))
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
        plt.savefig(SPLAT_PATH+DOCS_FOLDER+'_images/spt_nir_distribution_{}.png'.format(fname))
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
        plt.savefig(SPLAT_PATH+DOCS_FOLDER+'_images/spt_adopted_distribution_{}.png'.format(fname))
        plt.clf()

# histogram of S/N

# map sources on sky
    raref = Angle(numpy.linspace(0,359.,360)*u.degree)
    raref.wrap_at(180.*u.degree)
    ra = Angle(sall['RA'].filled(numpy.nan)*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(sall['DEC'].filled(numpy.nan)*u.degree)
    rap = Angle(s['RA'].filled(numpy.nan)*u.degree)
    rap = rap.wrap_at(180*u.degree)
    decp = Angle(s['DEC'].filled(numpy.nan)*u.degree)
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
    fig.savefig(SPLAT_PATH+DOCS_FOLDER+'_images/map_all.png')
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
    ra = Angle(sm['RA'].filled(numpy.nan)*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(sm['DEC'].filled(numpy.nan)*u.degree)
    p1 = ax.scatter(ra.radian, dec.radian,color='k',alpha=0.5,s=10)
    ra = Angle(sl['RA'].filled(numpy.nan)*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(sl['DEC'].filled(numpy.nan)*u.degree)
    p2 = ax.scatter(ra.radian, dec.radian,color='r',alpha=0.5,s=10)
    ra = Angle(st['RA'].filled(numpy.nan)*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(st['DEC'].filled(numpy.nan)*u.degree)
    p3 = ax.scatter(ra.radian, dec.radian,color='b',alpha=0.5,s=10)
    plt.legend([p1,p2,p3],['M dwarfs ({})'.format(len(sm)),'L dwarfs ({})'.format(len(sl)),'T dwarfs ({})'.format(len(st))],bbox_to_anchor=(1, 1),bbox_transform=plt.gcf().transFigure)
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    ax.grid(True)
    fig.savefig(SPLAT_PATH+DOCS_FOLDER+'_images/map_byspt.png')
    fig.clf()

# map sources on based on young or field
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    sy = splat.searchLibrary(young=True)
    sy = sy[sy['OBJECT_TYPE'] == 'VLM']
    so = splat.searchLibrary()
    so = so[so['OBJECT_TYPE'] == 'VLM']
    ra = Angle(so['RA'].filled(numpy.nan)*u.degree)
#    ra = ra.wrap_at(180*u.degree)
#    dec = Angle(so['DEC'].filled(numpy.nan)*u.degree)
#    p1 = ax.scatter(ra.radian, dec.radian,color='k',alpha=0.1,s=5)
    ra = Angle(sy['RA'].filled(numpy.nan)*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = Angle(sy['DEC'].filled(numpy.nan)*u.degree)
    p1 = ax.scatter(ra.radian, dec.radian,color='r',alpha=0.5,s=10)
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    ax.grid(True)
    plt.legend([p1],['Young Sources ({})'.format(len(sy))],bbox_to_anchor=(1, 1),bbox_transform=plt.gcf().transFigure)
    fig.savefig(SPLAT_PATH+DOCS_FOLDER+'_images/map_young.png')
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
    fig.savefig(SPLAT_PATH+DOCS_FOLDER+'_images/object_othertypes.png')


