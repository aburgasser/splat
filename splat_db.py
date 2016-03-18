from __future__ import print_function, division

"""
.. note::
         These are the database functions for SPLAT 
"""

import astropy
#import copy
#from datetime import datetime
import os
import re
import requests
import splat
import sys
#from scipy import stats
import numpy
from astropy.io import ascii, fits            # for reading in spreadsheet
from astropy.table import Table, join            # for reading in table files

DB_FOLDER = '/db/'
DB_ORIGINAL_FILE = 'db_spexprism.txt'
DB_SOURCES_FILE = 'source_data.txt'
DB_SPECTRA_FILE = 'spectral_data.txt'
DB_PHOTOMETRY_FILE = 'photometry_data.txt'
BIBFILE = 'biblibrary.bib'
TMPFILENAME = 'splattmpfile'

# change the command prompt
sys.ps1 = 'splat db> '


def assignSourceID(coordinate,**kwargs):
    '''
    :Purpose:
        Assigns a source ID based on current source database and possible overlap.

    :Note:
        **Currently not functional**

    :Required parameters:
        :param coordinate: astropy Coordinate object with coordinates (RA, Dec) of source, to check against current database

    :Optional parameters:
        :param new: Force the return of a new ID number regardless of search
        :type logical: optional, default = False
        :param validate: Ask for validation from user that a source is new
        :type logical: optional, default = False

    :Output:
        - An id key (long integer) for the source

    '''
    pass


def assignSpectralID(**kwargs):
    '''
    :Purpose:
        Assigns a spectral ID based on current spectral database

    :Note:
        **Currently not functional**

    :Required parameters:
        None

    :Optional parameters:
        None

    :Output:
        - An id key (long integer) for the spectrum

    '''
    pass


def bibTexParser(bib_tex,**kwargs):
    '''
    :Purpose:
        Parses a bibtex segment and returns a dictionary of parameter fields

    :Required parameters:
        :param bib_tex: String containing bibtex data in standard format

    :Optional parameters:
        None

    :Output:
        A dictionary containing the parsed bibtex information

    '''
    bib_dict = {"bib_tex": bib_tex}
    bib_tex.strip('\n')
    # get bib code
    begin = bib_tex.find('{')  
    end = bib_tex.find(',')
    bib_dict["bibcode"] = bib_tex[begin+1:end]
    bib_tex = bib_tex[end+1:]   # remove bib code line
    
    bib_tex =  bib_tex.split(',\n')  # this moght not always work for author lists
    
    for line in bib_tex:
        line = line.strip()
        line = line.replace('{','').replace('}','').replace('\"','').replace('\n','').replace('\t','') 
        line = line.split('=')
        line[0] = line[0].strip().lower()
        line[1] = line[1].strip()
        bib_dict[line[0]] = line[1]

# Journal massaging
    if bib_dict['journal'] == '\\apj':
        bib_dict['journal'] = 'ApJ'
    elif bib_dict['journal'] == '\\apjs':
        bib_dict['journal'] = 'ApJS'
    elif bib_dict['journal'] == '\\aj':
        bib_dict['journal'] = 'AJ'
    elif bib_dict['journal'] == '\\araa':
        bib_dict['journal'] = 'AR&A'
    elif bib_dict['journal'] == '\\aap':
        bib_dict['journal'] = 'A&A'
    elif bib_dict['journal'] == '\\mnras':
        bib_dict['journal'] = 'MNRAS'
    elif bib_dict['journal'] == '\\pasp':
        bib_dict['journal'] = 'PASP'
    elif bib_dict['journal'] == '\\pnas':
        bib_dict['journal'] = 'PNAS'
    else:
        pass


        
    return bib_dict


def shortRef(bib_dict,**kwargs):
    '''
    :Purpose:
        Takes a bibtex dictionary and returns a short (in-line) version of the citation

    :Required parameters:
        :param bib_tex: Dictionary output from bibTexParser, else a bibcode that is fed into bibTexParser

    :Optional parameters:
        None

    :Output:
        A string of the format ``Burgasser, A. J., et al. (2006, ApJ, 710, 1142)``

    '''
    if type(bib_dict) is not dict:
        if type(bib_dict) is str:
            bib_dict = getBibTex(bib_dict,**kwargs)
        else:
            raise NameError('Input to shortRef is neither a bibcode nor a bibTex dictionary')

    authors = bib_dict['author'].split(' and ')
    if len(authors) == 1:
        output = '{}'.format(authors[0].replace('~',' '))
    elif len(authors) == 2:
        output = '{} & {}'.format(authors[0].replace('~',' '),authors[1].replace('~',' '))
#    elif len(a) == 3:
#        output = '{}, {} & {}'.format(a[0].replace('~',' '),a[1].replace('~',' '),a[2].replace('~',' '))
#    else:
#        output = '{}, {}, {}, et al.'.format(a[0].replace('~',' '),a[1].replace('~',' '),a[2].replace('~',' '))
    else:
        output = '{} et al.'.format(authors[0].replace('~',' '))

    return output+' ({}, {}, {}, {})'.format(bib_dict['year'],bib_dict['journal'],bib_dict['volume'],bib_dict['pages'])


def longRef(bib_dict,**kwargs):
    '''
    :Purpose:
        Takes a bibtex dictionary and returns a long (in-line) version of the citation

    :Required parameters:
        :param bib_tex: Dictionary output from bibTexParser, else a bibcode that is fed into bibTexParser

    :Optional parameters:
        None

    :Output:
        A string of the format ``Burgasser, A. J., Cruz, K. L., Cushing, M., et al. SpeX Spectroscopy of Unresolved Very Low Mass Binaries. 
        I. Identification of 17 Candidate Binaries Straddling the L Dwarf/T Dwarf Transition. ApJ 710, 1142 (2010)``

    '''
    if type(bib_dict) is not dict:
        if type(bib_dict) is str:
            bib_dict = getBibTex(bib_dict,**kwargs)
        else:
            raise NameError('Input to shortRef is neither a bibcode nor a bibTex dictionary')

    authors = bib_dict['Author'].split(' and ')
    if len(authors) == 1:
        output = '{}'.format(authors[0].replace('~',' '))
    elif len(authors) == 2:
        output = '{} & {}'.format(authors[0].replace('~',' '),authors[1].replace('~',' '))
    elif len(authors) == 3:
        output = '{}, {} & {}'.format(authors[0].replace('~',' '),authors[1].replace('~',' '),authors[2].replace('~',' '))
    else:
        output = '{}, {}, {}, et al'.format(authors[0].replace('~',' '),authors[1].replace('~',' '),authors[2].replace('~',' '))

    return output+'. {}. {}, {}, {} ({})'.format(bib_dict['title'],bib_dict['journal'],bib_dict['volume'],bib_dict['pages'],bib_dict['year'])



def getBibTex(bibcode,**kwargs):
    '''
    Purpose
        Takes a bibcode and returns a dictionary containing the bibtex information; looks either in internal SPLAT
            or user-supplied bibfile, or seeks online. If nothing found, gives a soft warning and returns False

    :Note:
        **Currently not functional**

    :Required parameters:
        :param bibcode: Bibcode string to look up (e.g., '2014ApJ...787..126L')

    :Optional parameters:
        :param biblibrary: Filename for biblibrary to use in place of SPLAT internal one
        :type string: optional, default = ''
        :param online: If True, go directly online; if False, do not try to go online 
        :type logical: optional, default = null

    :Output:
        - A dictionary containing the bibtex fields, or False if not found

    '''

# go online first if directed to do so
    if kwargs.get('online',False) and checkOnline():
        bib_tex = getBibTexOnline(bibcode)

# read locally first
    else:
        biblibrary = kwargs.get('biblibrary', splat.SPLAT_PATH+DB_FOLDER+BIBFILE)
# check the file
        if not os.path.exists(biblibrary):
            print('Could not find bibtex library {}'.format(biblibrary))
            biblibrary = splat.SPLAT_PATH+DB_FOLDER+BIBFILE

        if not os.path.exists(biblibrary):
            raise NameError('Could not find SPLAT main bibtext library {}; something is wrong'.format(biblibrary))


        with open(biblibrary, 'r') as bib_file:
            text = bib_file.read()
            #print re.search('@[A-Z]+{' + bib_code, bib_file)        
            in_lib = re.search('@[a-z]+{' + bibcode, text)
            if in_lib == None:  
                print('Bibcode {} not in bibtex library {}; checking online'.format(bibcode,biblibrary))
                bib_tex = getBibTexOnline(bibcode)
            else:
                begin = text.find(re.search('@[a-z]+{' + bibcode, text).group(0))
                text = text[begin:]
                end = text.find('\n@')
                bib_tex = text[:end]

    if bib_tex == False:
        return False
    else:
        return bibTexParser(bib_tex)


def getBibTexOnline(bibcode):
    '''
    Purpose
        Takes a bibcode and searches for the bibtex information online through NASA ADS; requires user to be online.
            If successful, returns full bibtex string block; otherwise False.

    :Required parameters:
        :param bibcode: Bibcode string to look up (e.g., '2014ApJ...787..126L')

    :Optional parameters:
        :param bibfile: Filename for bibfile to use in place of SPLAT internal one
        :type string: optional, default = ''
        :param online: If True, go directly online; if False, do not try to go online 
        :type logical: optional, default = null

    :Output:
        - A string block of the basic bibtex information

    '''
    if (checkOnline() == False):
        return False

    url_begin = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode="
    url_end = "&data_type=BIBTEX"
    url = url_begin + bibcode + url_end
    bib_tex = requests.get(url).content
    
    # Check if content is in html which means bad bib_code was given
    if "<HTML>" in bib_tex:
        print('{} is not a valid online bib code.'.format(bibcode))
        return False       
        
    # Cut off extraneous info from website before the bibtex code
    else:
        begin = bib_tex.find('@')
        bib_tex = bib_tex[begin:]
        return bib_tex



def checkFile(filename,**kwargs):
    '''
    :Purpose: Checks if a spectrum file exists in the SPLAT's library.
    :param filename: A string containing the spectrum's filename.
    :Example:
       >>> import splat
       >>> spectrum1 = 'spex_prism_1315+2334_110404.fits'
       >>> print splat.checkFile(spectrum1)
       True
       >>> spectrum2 = 'fake_name.fits'
       >>> print splat.checkFile(spectrum2)
       False
    '''
    url = kwargs.get('url',splat.SPLAT_URL)+DATA_FOLDER
    return requests.get(url+filename).status_code == requests.codes.ok
#    flag = checkOnline()
#    if (flag):
#        try:
#            r = requests.get(url+filename)
#            open(os.path.basename(filename), 'wb').write(r.content)
#            open(os.path.basename(filename), 'wb').write(urllib2.urlopen(url+filename).read())
#        except:
#            flag = False
#    return flag


def checkAccess(**kwargs):
    '''
    :Purpose: Checks if user has access to unpublished spectra in SPLAT library.
    :Example:
       >>> import splat
       >>> print splat.checkAccess()
       True
    :Note: Must have the file .splat_access in your home directory with the correct passcode to use.
    '''
    access_file = '.splat_access'
    result = False

    try:
        home = os.environ.get('HOME')
        if home == None:
            home = './'
        bcode = requests.get(splat.SPLAT_URL+access_file).content
        lcode = base64.b64encode(open(home+'/'+access_file,'r').read())
        if (bcode in lcode):        # changed to partial because of EOL variations
            result = True
    except:
        result = False

    if (kwargs.get('report','') != ''):
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
       >>> splat.checkLocal('splat.py')
       True  # found the code
       >>> splat.checkLocal('parameters.txt')
       False  # can't find this file
       >>> splat.checkLocal('SpectralModels/BTSettl08/parameters.txt')
       True  # found it
    '''
    if not os.path.exists(inputfile):
        if not os.path.exists(splat.SPLAT_PATH+inputfile):
            return ''
        else:
            return splat.SPLAT_PATH+inputfile
    else:
        return inputfile


def checkOnline(*args):
    '''
    :Purpose: Checks if SPLAT's URL is accessible from your machine--
                that is, checks if you and the host are online. Alternately
                checks if a given filename is present locally or online
    :Example:
       >>> import splat
       >>> splat.checkOnline()
       True  # SPLAT's URL was detected.
       >>> splat.checkOnline()
       False # SPLAT's URL was not detected.
       >>> splat.checkOnline('SpectralModels/BTSettl08/parameters.txt')
       '' # Could not find this online file.
    '''
    if (len(args) != 0):
        if 'http://' in args[0]:
            if requests.get(args[0]).status_code == requests.codes.ok:
                return args[0]
            return ''
        else:
            if requests.get(splat.SPLAT_URL+args[0]).status_code == requests.codes.ok:
                return splat.SPLAT_URL+args[0]
            return ''
    else:
        return requests.get(splat.SPLAT_URL).status_code == requests.codes.ok


def fetchDatabase(*args, **kwargs):
    '''
    :Purpose: Get the SpeX Database from either online repository or local drive
    '''
    filename = DB_ORIGINAL_FILE
    if len(args) > 0:
        filename = args[0]
    kwargs['filename'] = kwargs.get('filename',filename)
    kwargs['filename'] = kwargs.get('file',kwargs['filename'])
    kwargs['folder'] = kwargs.get('folder',DB_FOLDER)
    url = kwargs.get('url',splat.SPLAT_URL)+kwargs['folder']
    local = kwargs.get('local',True)
    online = kwargs.get('online',not local and checkOnline())
    local = not online
    kwargs['local'] = local
    kwargs['online'] = online
    kwargs['model'] = True

# check that folder/set is present either locally or online
# if not present locally but present online, switch to this mode
# if not present at either raise error
    folder = checkLocal(kwargs['folder'])
    if folder=='':
        folder = checkOnline(kwargs['folder'])
        if folder=='':
            raise NameError('\nCould not find '+kwargs['folder']+' locally or on SPLAT website\n\n')
        else:
            kwargs['folder'] = folder
            kwargs['local'] = False
            kwargs['online'] = True
    else:
        kwargs['folder'] = folder

# locally:
    if kwargs['local']:
        infile = checkLocal(kwargs['filename'])
        if infile=='':
            infile = checkLocal(kwargs['folder']+'/'+kwargs['filename'])
        if infile=='':
            raise NameError('\nCould not find '+kwargs['filename']+' locally\n\n')
        else:
            try:
                data = ascii.read(infile, delimiter='\t',fill_values='-99.',format='tab')
            except:
                raise NameError('\nCould not load {}: this may be a decoding error\n'.format(infile))


# check if file is present; if so, read it in, otherwise go to interpolated
# online:
    if kwargs['online']:
        infile = checkOnline(kwargs['filename'])
        if infile=='':
            infile = checkOnline(kwargs['folder']+'/'+kwargs['filename'])
        if infile=='':
            raise NameError('\nCould not find '+kwargs['filename']+' on the SPLAT website\n\n')
        try:
#            open(os.path.basename(TMPFILENAME), 'wb').write(urllib2.urlopen(url+infile).read())
            open(os.path.basename(TMPFILENAME), 'wb').write(requests.get(url+infile).content)
            kwargs['filename'] = os.path.basename(tmp)
            data = ascii.read(os.path.basename(TMPFILENAME), delimiter='\t',fill_values='-99.',format='tab')
            os.remove(os.path.basename(TMPFILENAME))
        except:
            raise NameError('\nHaving a problem reading in '+kwargs['filename']+' on the SPLAT website\n\n')


# clean up blanks and convert numerical values to numbers
#    if 'RA' in data.keys():
#        data['RA'][numpy.where(data['RA'] == '')] = '0.'
#    data['ran'] = [float(x) for x in data['ra']]
#    if 'DEC' in data.keys():
#        data['DEC'][numpy.where(data['DEC'] == '')] = '0.'
#    data['decn'] = [float(x) for x in data['dec']]
#    if '2MASS_J' in data.keys():
#        data['2MASS_J'][numpy.where(data['2MASS_J'] == '')] = '99.'
#    data['jmagn'] = [float(x) for x in data['jmag']]
#    if '2MASS_H' in data.keys():
#        data['2MASS_H'][numpy.where(data['2MASS_H'] == '')] = '99.'
#    data['hmagn'] = [float(x) for x in data['hmag']]
#    if '2MASS_KS' in data.keys():
#        data['2MASS_KS'][numpy.where(data['2MASS_KS'] == '')] = '99.'
#    data['kmagn'] = [float(x) for x in data['kmag']]
#    if '2MASS_J_UNC' in data.keys():
#        data['2MASS_J_UNC'][numpy.where(data['2MASS_J_UNC'] == '')] = '99.'
#    data['jmag_errorn'] = [float(x) for x in data['jmag_error']]
#    if '2MASS_H_UNC' in data.keys():
#        data['2MASS_H_UNC'][numpy.where(data['2MASS_H_UNC'] == '')] = '99.'
#    data['hmag_errorn'] = [float(x) for x in data['hmag_error']]
#    if '2MASS_KS_UNC' in data.keys():
#        data['2MASS_KS'][numpy.where(data['2MASS_KS_UNC'] == '')] = '99.'
#    data['kmag_errorn'] = [float(x) for x in data['kmag_error']]
#    if 'RESOLUTION' in data.keys():
#        data['RESOLUTION'][numpy.where(data['RESOLUTION'] == '')] = '120'
#    data['resolutionn'] = [float(x) for x in data['resolution']]
#    if 'AIRMASS' in data.keys():
#       data['AIRMASS'][numpy.where(data['AIRMASS'] == '')] = '1.'
#    data['airmassn'] = [float(x) for x in data['airmass']]
#    if 'MEDIAN_SNR' in data.keys():
#        data['MEDIAN_SNR'][numpy.where(data['MEDIAN_SNR'] == '')] = '0'
#    data['median_snrn'] = [float(x) for x in data['median_snr']]

# convert coordinates to SkyCoord format
#    data['skycoords'] = data['ra']
#    s = []
#    for i in numpy.arange(len(data['RA'])):
#        try:        # to deal with a blank string
#            s.append(SkyCoord(ra=float(data['RA'][i])*u.degree,dec=float(data['DEC'][i])*u.degree,frame='icrs'))
#        except:
#            s.append(SkyCoord(ra=0.*u.degree,dec=0.*u.degree,frame='icrs'))
#    data['SKYCOORDS'] = s

# add in RA/Dec (TEMPORARY)
#    ra = []
#    dec = []
#    for x in data['designation']:
#        c = designationToCoordinate(x,ICRS=False)
#        ra.append(c[0])
#        dec.append(c[1])
#    data['ra'] = ra
#    data['dec'] = dec

#    data['YOUNG'] = ['young' in x for x in data['LIBRARY']]
#    data['SUBDWARF'] = ['subdwarf' in x for x in data['LIBRARY']]
#    data['BINARY'] = ['binary' in x for x in data['LIBRARY']]
#    data['SPBINARY'] = ['sbinary' in x for x in data['LIBRARY']]
#    data['BLUE'] = ['blue' in x for x in data['LIBRARY']]
#    data['RED'] = ['red' in x for x in data['LIBRARY']]
#    data['GIANT'] = ['giant' in x for x in data['LIBRARY']]
#    data['WD'] = ['wd' in x for x in data['LIBRARY']]
#    data['STANDARD'] = ['std' in x for x in data['LIBRARY']]
#    data['COMPANION'] = ['companion' in x for x in data['LIBRARY']]

# add in shortnames
#    data['SHORTNAME'] = [splat.designationToShortName(x) for x in data['DESIGNATION']]

# create literature spt
#    data['LIT_TYPE'] = data['OPT_TYPE']
#    w = numpy.where(numpy.logical_and(data['LIT_TYPE'] == '',data['NIR_TYPE'] != ''))
#    data['LIT_TYPE'][w] = data['NIR_TYPE'][w]
#    sptn = [splat.typeToNum(x) for x in data['LIT_TYPE']]
#    w = numpy.where(numpy.logical_and(sptn > 29.,data['NIR_TYPE'] != ''))
#    data['LIT_TYPE'][w] = data['NIR_TYPE'][w]
#    w = numpy.where(numpy.logical_and(data['lit_type'] == '',typeToNum(data['spex_type']) > 17.))
#    data['lit_type'][w] = data['spex_type'][w]

    return data



def getPhotometry(coordinate,**kwargs):
    '''
    Purpose
        Downloads photometry for a source using astroquery (?)

    :Note:
        **Currently not functional**

    :Required parameters:
        :param coordinate: astropy Coordinate object with coordinates (RA, Dec) of source to be searched

    :Optional parameters:
 
        :param radius: radius in arcseconds for matching
        :type float: optional, default = 5.
        :param validate: Ask for validation from user for matched photometry
        :type logical: optional, default = False
        :param 2MASS: Download 2MASS photometry
        :type logical: optional, default = True
        :param SDSS: Download 2MASS photometry
        :type logical: optional, default = True
        :param UKIDSS: Download 2MASS photometry
        :type logical: optional, default = True
        :param WISE: Download 2MASS photometry
        :type logical: optional, default = True
        :param DENIS: Download 2MASS photometry
        :type logical: optional, default = True

    :Output:
        - A table? filename? if new photometry

    '''
    pass




def importSpectra(input,**kwargs):
    '''
    Purpose
        imports a set of spectra into the SPLAT library; requires manager access.

    :Note:
        **Currently not functional**

    :Required parameters:
        :param input: Can be one of the following:
            - A string or array of strings containing spectra filenames, which must either be fits files
                or an ascii files with three columns (tab- or comma-delimited) containing wavelength, flux and noise.
            - A string filename for a spreadsheet (ascii, tab- or comma-delimited) listing the input spectra, one per row.
                At least one column must be ``filename``; the remaining may be one of the optional inputs listed below

    :Optional parameters:
        :param source_key: Source identification key if a source is already in the SPLAT database; use this for a *new* observation of a previously observed source
        :type long: optional, default = null
        :param spectrum_key: Spectrum identification if a spectrum is already in the SPLAT database; use this for a *replacing* a previous observation
        :type long: optional, default = null
        :param name: Source name
        :type string: optional, default = ''
        :param ra: Source Right Ascension (used for searching for catalog information)
        :type double: optional, default = null
        :param dec: Source Declination (used for searching for catalog information)
        :type double: optional, default = null
        :param designation: Source designation (e.g., 'J1234567+0123456'; used for searching for catalog information)
        :type string: optional, default = ''
        :param opt_type: Source optical spectral type (e.g., 'M7')
        :type string: optional, default = ''
        :param opt_type_reference: Bibcode for optical type reference (e.g., '2014ApJ...787..126L')
        :type string: optional, default = ''
        :param nir_type: Source near-infrared spectral type (e.g., 'M7')
        :type string: optional, default = ''
        :param nir_type_reference: Bibcode for optical type reference (e.g., '2014ApJ...787..126L')
        :type string: optional, default = ''
        :param publication_reference: Publication reference bibcode (e.g., '2014ApJ...787..126L')
        :type string: optional, default = ''
        :param observer: Name of observer
        :type string: optional, default = ''
        :param date_observed: Observation date (e.g., '20100425', '2010 Apr 25')
        :type string: optional, default = ''

    :Output:
        - Source DB update file: spreadsheet containing update to source_data.txt, saved locally as UPDATE_source_data.txt
        - Spectral DB update file: spreadsheet containing update to spectral_data.txt, saved locally as UPDATE_spectral_data.txt
        - Photometry DB update file: spreadsheet containing update to photometry_data.txt, saved locally as UPDATE_photometry_data.txt

    '''
    pass


def keySource(keys, **kwargs):
    '''
    :Purpose: Takes a source key and returns a table with the source information
    :param keys: source key or a list of source keys
    :Example:
    >>> import splat
    >>> print splat.keySource(10001)
        SOURCE_KEY           NAME              DESIGNATION    ... NOTE SELECT
        ---------- ------------------------ ----------------- ... ---- ------
             10001 SDSS J000013.54+255418.6 J00001354+2554180 ...        True
    >>> print splat.keySource([10105, 10623])
        SOURCE_KEY          NAME             DESIGNATION    ... NOTE SELECT
        ---------- ---------------------- ----------------- ... ---- ------
             10105 2MASSI J0103320+193536 J01033203+1935361 ...        True
             10623 SDSS J09002368+2539343 J09002368+2539343 ...        True
    >>> print splat.keySource(1000001)
        No sources found with source key 1000001
        False
    '''

# vectorize
    if isinstance(keys,list) == False:
        keys = [keys]

#    sdb = ascii.read(splat.SPLAT_PATH+DB_FOLDER+SOURCES_DB, delimiter='\t',fill_values='-99.',format='tab')
#    sdb = fetchDatabase(splat.SPLAT_PATH+DB_FOLDER+SOURCES_DB)
    sdb = splat.DB_SOURCES
    sdb['SELECT'] = [x in keys for x in sdb['SOURCE_KEY']]

    if sum(sdb['SELECT']) == 0.:
        print('No sources found with source key {}'.format(keys[0]))
        return False
    else:
        db = sdb[:][numpy.where(sdb['SELECT']==1)]
        return db


def keySpectrum(keys, **kwargs):
    '''
    :Purpose: Takes a spectrum key and returns a table with the spectrum and source information
    :param keys: spectrum key or a list of source keys
    :Example:
    >>> import splat
    >>> print splat.keySpectrum(10001)
        DATA_KEY SOURCE_KEY    DATA_FILE     ... COMPANION COMPANION_NAME NOTE_2
        -------- ---------- ---------------- ... --------- -------------- ------
           10001      10443 10001_10443.fits ...
    >>> print splat.keySpectrum([10123, 11298])
        DATA_KEY SOURCE_KEY    DATA_FILE     ... COMPANION COMPANION_NAME NOTE_2
        -------- ---------- ---------------- ... --------- -------------- ------
           11298      10118 11298_10118.fits ...
           10123      10145 10123_10145.fits ...
    >>> print splat.keySpectrum(1000001)
        No spectra found with spectrum key 1000001
        False
    '''

# vectorize
    if isinstance(keys,list) == False:
        keys = [keys]

#    sdb = ascii.read(splat.SPLAT_PATH+DB_FOLDER+SPECTRA_DB, delimiter='\t',fill_values='-99.',format='tab')
#    sdb = fetchDatabase(splat.SPLAT_PATH+DB_FOLDER+SPECTRA_DB)
    sdb = splat.DB_SPECTRA
    sdb['SELECT'] = [x in keys for x in sdb['DATA_KEY']]

    if sum(sdb['SELECT']) == 0.:
        print('No spectra found with spectrum key {}'.format(keys[0]))
        return False
    else:
#        s2db = ascii.read(splat.SPLAT_PATH+DB_FOLDER+SOURCES_DB, delimiter='\t',fill_values='-99.',format='tab')
#        s2db = fetchDatabase(splat.SPLAT_PATH+DB_FOLDER+SOURCES_DB)
        s2db = splat.DB_SOURCES
        db = join(sdb[:][numpy.where(sdb['SELECT']==1)],s2db,keys='SOURCE_KEY')
        return db



def searchLibrary(*args, **kwargs):
    '''
    :Purpose: Search the SpeX database to extract the key reference for that Spectrum
    :param output: returns desired output of selected results
    :type output: optional, default = 'all'
    :param logic: search logic, can be and`` or ``or``
    :type logic: optional, default = 'and'
    :param combine: same as logic
    :type combine: optional, default = 'and'
    :Example:
    >>> import splat
    >>> print SearchLibrary(shortname = '2213-2136')
        DATA_KEY SOURCE_KEY    DATA_FILE     ... SHORTNAME  SELECT_2
        -------- ---------- ---------------- ... ---------- --------
           11590      11586 11590_11586.fits ... J2213-2136      1.0
           11127      11586 11127_11586.fits ... J2213-2136      1.0
           10697      11586 10697_11586.fits ... J2213-2136      1.0
           10489      11586 10489_11586.fits ... J2213-2136      1.0
    >>> print SearchLibrary(shortname = '2213-2136', output = 'OBSERVATION_DATE')
        OBSERVATION_DATE
        ----------------
                20110908
                20080829
                20060902
                20051017

    .. note:: Note that this is currently only and AND search - need to figure out how to a full SQL style search
    '''

# program parameters
    ref = kwargs.get('output','all')
    radius = kwargs.get('radius',10.)      # search radius in arcseconds
    classes = ['YOUNG','SUBDWARF','BINARY','SPBINARY','RED','BLUE','GIANT','WD','STANDARD','COMPANION']

# logic of search
    logic = 'and'         # default combination
    logic = kwargs.get('combine',logic).lower()
    logic = kwargs.get('logic',logic).lower()
    if (logic != 'and' and logic != 'or'):
        raise ValueError('\nLogical operator '+logic+' not supported\n\n')

# read in source database and add in shortnames and skycoords
#    source_db = ascii.read(splat.SPLAT_PATH+DB_FOLDER+SOURCES_DB, delimiter='\t', fill_values='-99.', format='tab')
#    source_db = fetchDatabase(SOURCES_DB)
    source_db = splat.DB_SOURCES
    source_db['SHORTNAME'] = [splat.designationToShortName(x) for x in source_db['DESIGNATION']]

# first search by source parameters
    source_db['SELECT'] = numpy.zeros(len(source_db['RA']))
    count = 0.

# search by source key
    idkey = kwargs.get('sourcekey',False)
    idkey = kwargs.get('idkey',idkey)
    idkey = kwargs.get('id',idkey)
    if idkey != False:
        if isinstance(idkey,int):
            idkey = [idkey]
        for s in idkey:
            source_db['SELECT'][numpy.where(source_db['SOURCE_KEY'] == s)] += 1
        count+=1.
# search by name
    if kwargs.get('name',False) != False:
        nm = kwargs['name']
        if isinstance(nm,str):
            nm = [nm]
        for n in nm:
            source_db['SELECT'][numpy.where(source_db['NAME'] == n)] += 1
        count+=1.
# search by shortname
    if kwargs.get('shortname',False) != False:
        sname = kwargs['shortname']
        if isinstance(sname,str):
            sname = [sname]
        for sn in sname:
            if sn[0].lower() != 'j':
                sn = 'J'+sn
            source_db['SELECT'][numpy.where(source_db['SHORTNAME'] == sn)] += 1
        count+=1.
# exclude by shortname
    sname = kwargs.get('excludesource',False)
    sname = kwargs.get('excludeshortname',sname)
    if sname != False and len(sname) > 0:
        if isinstance(sname,str):
            sname = [sname]
        for sn in sname:
            if sn[0].lower() != 'j':
                sn = 'J'+sn
#            t = numpy.sum(source_db['SELECT'][numpy.where(source_db['SHORTNAME'] != sn)])
            source_db['SELECT'][numpy.where(source_db['SHORTNAME'] != sn)] += 1
#            if numpy.sum(source_db['SELECT'][numpy.where(source_db['SHORTNAME'] != sn)]) > t:
#                print('rejected '+sn)
        count+=1.
# search by reference list
    if kwargs.get('reference',False) != False:
        refer = kwargs['reference']
        if isinstance(ref,str):
            refer = [refer]
        for r in refer:
            source_db['SELECT'][numpy.where(source_db['DISCOVERY_REFERENCE'] == r)] += 1
        count+=1.
# search by designation
    if kwargs.get('designation',False) != False:
        desig = kwargs['designation']
        if isinstance(desig,str):
            desig = [desig]
        for d in desig:
            source_db['SELECT'][numpy.where(source_db['DESIGNATION'] == d)] += 1
        count+=1.
# search by coordinate - NOTE: THIS IS VERY SLOW RIGHT NOW
    if kwargs.get('coordinate',False) != False:
        coord = kwargs['coordinate']
        if isinstance(coord,SkyCoord):
            cc = coord
        else:
            cc = properCoordinates(coord)
# calculate skycoords
        s = []
        for i in numpy.arange(len(source_db['RA'])):
            try:        # to deal with a blank string
                s.append(SkyCoord(ra=float(source_db['RA'][i])*u.degree,dec=float(source_db['DEC'][i])*u.degree,frame='icrs'))
            except:
                s.append(SkyCoord(ra=numpy.nan*u.degree,dec=numpy.nan*u.degree,frame='icrs'))
        source_db['SKYCOORDS'] = s
        source_db['SEPARATION'] = [cc.separation(source_db['SKYCOORDS'][i]).arcsecond for i in numpy.arange(len(source_db['SKYCOORDS']))]
        source_db['SELECT'][numpy.where(source_db['SEPARATION'] <= radius)] += 1
        count+=1.

# search by spectral type
    spt_range = kwargs.get('spt_range',False)
    spt_range = kwargs.get('spt',spt_range)
    spt_type = kwargs.get('spt_type','LIT_TYPE')
    if spt_range != False:
        if spt_type not in ['LIT_TYPE','SPEX_TYPE','OPT_TYPE','NIR_TYPE']:
            spt_type = 'LIT_TYPE'
        if not isinstance(spt_range,list):        # one value = only this type
            spt_range = [spt_range,spt_range]
        if isinstance(spt_range[0],str):          # convert to numerical spt
            spt_range = [splat.typeToNum(spt_range[0]),splat.typeToNum(spt_range[1])]
        source_db['SPTN'] = [splat.typeToNum(x) for x in source_db[spt_type]]
        source_db['SELECT'][numpy.where(numpy.logical_and(source_db['SPTN'] >= spt_range[0],source_db['SPTN'] <= spt_range[1]))] += 1
        count+=1.

# search by magnitude range
    if kwargs.get('jmag',False) != False:
        mag = kwargs['jmag']
        if not isinstance(mag,list):        # one value = faint limit
            mag = [0,mag]
        source_db['JMAGN'] = [float('0'+x) for x in source_db['2MASS_J']]
        source_db['SELECT'][numpy.where(numpy.logical_and(source_db['JMAGN'] >= mag[0],source_db['JMAGN'] <= mag[1]))] += 1
        count+=1.
    if kwargs.get('hmag',False) != False:
        mag = kwargs['hmag']
        if not isinstance(mag,list):        # one value = faint limit
            mag = [0,mag]
        source_db['HMAGN'] = [float('0'+x) for x in source_db['2MASS_H']]
        source_db['SELECT'][numpy.where(numpy.logical_and(source_db['HMAGN'] >= mag[0],source_db['HMAGN'] <= mag[1]))] += 1
        count+=1.
    if kwargs.get('kmag',False) != False:
        mag = kwargs['kmag']
        if not isinstance(mag,list):        # one value = faint limit
            mag = [0,mag]
        source_db['KMAGN'] = [float('0'+x) for x in source_db['2MASS_KS']]
        source_db['SELECT'][numpy.where(numpy.logical_and(source_db['KMAGN'] >= mag[0],source_db['KMAGN'] <= mag[1]))] += 1
        count+=1.

# young
    if (kwargs.get('young','') != ''):
        source_db['YOUNG'] = [i != '' for i in source_db['GRAVITY_CLASS_CRUZ']] or [i != '' for i in source_db['GRAVITY_CLASS_ALLERS']]
        source_db['SELECT'][numpy.where(source_db['YOUNG'] == kwargs.get('young'))] += 1
        count+=1.

# specific gravity class
    flag = kwargs.get('gravity_class','')
    flag = kwargs.get('gravity',flag)
    if (flag != ''):
        source_db['SELECT'][numpy.where(source_db['GRAVITY_CLASS_CRUZ'] == flag)] += 1
        source_db['SELECT'][numpy.where(source_db['GRAVITY_CLASS_ALLERS'] == flag)] += 1
        count+=1.

# specific cluster
    if (kwargs.get('cluster','') != '' and isinstance(kwargs.get('cluster'),str)):
        source_db['CLUSTER_FLAG'] = [i.lower() == kwargs.get('cluster').lower() for i in source_db['CLUSTER']]
        source_db['SELECT'][numpy.where(source_db['CLUSTER_FLAG'] == True)] += 1
        count+=1.

# giant
    if (kwargs.get('giant','') != ''):
#        kwargs['vlm'] = False
        source_db['GIANT'] = [i != '' for i in source_db['LUMINOSITY_CLASS']]
        source_db['SELECT'][numpy.where(source_db['GIANT'] == kwargs.get('giant'))] += 1
        count+=1.

# luminosity class
    if (kwargs.get('giant_class','') != ''):
        if 'GIANT' not in source_db.keys():
            source_db['GIANT'] = [i != '' for i in source_db['LUMINOSITY_CLASS']]
        source_db['GIANT_FLAG'] = [i.lower() == kwargs.get('giant_class').lower() for i in source_db['GIANT']]
        source_db['SELECT'][numpy.where(source_db['GIANT_FLAG'] == True)] += 1
        count+=1.

# subdwarf
    if (kwargs.get('subdwarf','') != ''):
        source_db['SUBDWARF'] = [i != '' for i in source_db['METALLICITY_CLASS']]
        source_db['SELECT'][numpy.where(source_db['SUBDWARF'] == kwargs.get('subdwarf'))] += 1
        count+=1.

# metallicity class
    if (kwargs.get('subdwarf_class','') != ''):
        source_db['SD_FLAG'] = [i.lower() == kwargs.get('subdwarf_class').lower() for i in source_db['METALLICITY_CLASS']]
        source_db['SELECT'][numpy.where(source_db['SD_FLAG'] == True)] += 1
        count+=1.

# red
    if (kwargs.get('red','') != ''):
        source_db['RED'] = ['red' in i for i in source_db['LIBRARY']]
        source_db['SELECT'][numpy.where(source_db['RED'] == kwargs.get('red'))] += 1
        count+=1.

# blue
    if (kwargs.get('blue','') != ''):
        source_db['BLUE'] = ['blue' in i for i in source_db['LIBRARY']]
        source_db['SELECT'][numpy.where(source_db['BLUE'] == kwargs.get('blue'))] += 1
        count+=1.

# binaries
    if (kwargs.get('binary','') != ''):
        source_db['BINARY_FLAG'] = [i == 'Y' for i in source_db['BINARY']]
        source_db['SELECT'][numpy.where(source_db['BINARY_FLAG'] == kwargs.get('binary'))] += 1
        count+=1.

# spectral binaries
    if (kwargs.get('sbinary','') != ''):
        source_db['SBINARY_FLAG'] = [i == 'Y' for i in source_db['SBINARY']]
        source_db['SELECT'][numpy.where(source_db['SBINARY_FLAG'] == kwargs.get('sbinary'))] += 1
        count+=1.

# companions
    if (kwargs.get('companion','') != ''):
        source_db['COMPANION_FLAG'] = [i != '' for i in source_db['COMPANION_NAME']]
        source_db['SELECT'][numpy.where(source_db['COMPANION_FLAG'] == kwargs.get('companion'))] += 1
        count+=1.

# white dwarfs
    if (kwargs.get('wd','') != ''):
        kwargs['vlm'] = False
        source_db['WHITEDWARF'] = [i == 'WD' for i in source_db['OBJECT_TYPE']]
        source_db['SELECT'][numpy.where(source_db['WHITEDWARF'] == kwargs.get('wd'))] += 1
        count+=1.

# galaxies
    if (kwargs.get('galaxy','') != ''):
        kwargs['vlm'] = False
        source_db['GALAXY'] = [i == 'GAL' for i in source_db['OBJECT_TYPE']]
        source_db['SELECT'][numpy.where(source_db['GALAXY'] == kwargs.get('galaxy'))] += 1
        count+=1.

# carbon stars
    if (kwargs.get('carbon','') != ''):
        kwargs['vlm'] = False
        source_db['CARBON'] = [i == 'C' for i in source_db['OBJECT_TYPE']]
        source_db['SELECT'][numpy.where(source_db['CARBON'] == kwargs.get('carbon'))] += 1
        count+=1.

# VLM dwarfs by default
    if (kwargs.get('vlm','') != ''):
        if (kwargs.get('vlm') == True):
            source_db['SELECT'][numpy.where(source_db['OBJECT_TYPE'] == 'VLM')] += 1
            count+=1.
        if (kwargs.get('vlm') == False):
            source_db['SELECT'][numpy.where(source_db['OBJECT_TYPE'] != 'VLM')] += 1
            count+=1.

# select source keys
    if (count > 0):
        if (logic == 'and'):
            source_db['SELECT'] = numpy.floor(source_db['SELECT']/count)
        elif (logic == 'or'):
            source_db['SELECT'] = numpy.ceil(source_db['SELECT']/count)

        source_keys = source_db['SOURCE_KEY'][numpy.where(source_db['SELECT']==1)]
# no selection made on sources - choose everything
    else:
        source_keys = source_db['SOURCE_KEY']


# read in spectral database
#    spectral_db = ascii.read(splat.SPLAT_PATH+DB_FOLDER+SPECTRA_DB, delimiter='\t',fill_values='-99.',format='tab')
#    spectral_db = fetchDatabase(splat.SPLAT_PATH+DB_FOLDER+SPECTRA_DB)
    spectral_db = splat.DB_SPECTRA
    spectral_db['SELECT'] = numpy.zeros(len(spectral_db['DATA_KEY']))
    count = 0.

    spectral_db['SOURCE_SELECT'] = [x in source_keys for x in spectral_db['SOURCE_KEY']]

# search by filename
    file = kwargs.get('file','')
    file = kwargs.get('filename',file)
    if (file != ''):
        if isinstance(file,str):
            file = [file]
        for f in file:
            spectral_db['SELECT'][numpy.where(spectral_db['DATA_FILE'] == f)] += 1
        count+=1.

# exclude by data key
    if kwargs.get('excludekey',False) != False:
        exkey = kwargs['excludekey']
        if len(exkey) > 0:
            if isinstance(exkey,str):
                exkey = [exkey]
            for f in exkey:
                spectral_db['SELECT'][numpy.where(spectral_db['DATA_KEY'] != f)] += 1
            count+=1.

# exclude by filename
    if kwargs.get('excludefile',False) != False:
        file = kwargs['excludefile']
        if len(file) > 0:
            if isinstance(file,str):
                file = [file]
            for f in file:
                spectral_db['SELECT'][numpy.where(spectral_db['DATA_FILE'] != f)] += 1
            count+=1.

# search by observation date range
    if kwargs.get('date',False) != False:
        date = kwargs['date']
        if isinstance(date,str) or isinstance(date,long) or isinstance(date,float) or isinstance(date,int):
            date = [float(date),float(date)]
        elif isinstance(date,list):
            date = [float(date[0]),float(date[-1])]
        else:
            raise ValueError('\nCould not parse date input {}\n\n'.format(date))
        spectral_db['DATEN'] = [float(x) for x in spectral_db['OBSERVATION_DATE']]
        spectral_db['SELECT'][numpy.where(numpy.logical_and(spectral_db['DATEN'] >= date[0],spectral_db['DATEN'] <= date[1]))] += 1
        count+=1.

# search by S/N range
    if kwargs.get('snr',False) != False:
        snr = kwargs['snr']
        if not isinstance(snr,list):        # one value = minimum S/N
            snr = [float(snr),1.e9]
        spectral_db['SNRN'] = [float('0'+x) for x in spectral_db['MEDIAN_SNR']]
        spectral_db['SELECT'][numpy.where(numpy.logical_and(spectral_db['SNRN'] >= snr[0],spectral_db['SNRN'] <= snr[1]))] += 1
        count+=1.

# combine selection logically
    if (count > 0):
        if (logic == 'and'):
            spectral_db['SELECT'] = numpy.floor(spectral_db['SELECT']/count)
        else:
            spectral_db['SELECT'] = numpy.ceil(spectral_db['SELECT']/count)

    else:
        spectral_db['SELECT'] = numpy.ones(len(spectral_db['DATA_KEY']))

# limit access to public data for most users
    if (not checkAccess() or kwargs.get('published',False) or kwargs.get('public',False)):
        spectral_db['SELECT'][numpy.where(spectral_db['PUBLISHED'] != 'Y')] = 0.

# no matches
    if len(spectral_db[:][numpy.where(numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==1))]) == 0:
        if not kwargs.get('silent',False):
            print('No spectra in the SPL database match the selection criteria')
        return Table()
    else:

# merge databases
#        print(numpy.sum(spectral_db['SELECT']), numpy.sum(spectral_db['SOURCE_SELECT']))
#        print(spectral_db[:][numpy.where(spectral_db['SELECT']==1)])
#        print(spectral_db['SELECT'][numpy.where(spectral_db['SOURCE_SELECT']==1)])
#        print(len(spectral_db[:][numpy.where(numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==1))]))
        db = join(spectral_db[:][numpy.where(numpy.logical_and(spectral_db['SELECT']==1,spectral_db['SOURCE_SELECT']==1))],source_db,keys='SOURCE_KEY')

        if (ref == 'all'):
            return db
        else:
            return db[ref]




# main testing of program
if __name__ == '__main__':
    basefolder = '/Users/adam/projects/splat/exercises/ex9/'
    sp = splat.getSpectrum(shortname='1047+2124')[0]        # T6.5 radio emitter
    spt,spt_e = splat.classifyByStandard(sp,spt=['T2','T8'])
    teff,teff_e = splat.typeToTeff(spt)
    sp.fluxCalibrate('MKO J',splat.typeToMag(spt,'MKO J')[0],absolute=True)
    table = modelFitMCMC(sp, mask_standard=True, initial_guess=[teff, 5.3, 0.], zstep=0.1, nsamples=100,savestep=0,filebase=basefolder+'fit1047',verbose=True)

