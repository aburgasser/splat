# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
.. note::
         These are the database functions for SPLAT 
"""

# imports: internal
import base64
import copy
import csv
import glob
import os
import re
import requests
from shutil import copyfile
import time

# imports: external
import astropy
import numpy
import pandas
from astropy.io import ascii, fits            # for reading in spreadsheet
from astropy.table import Column, Table, join, vstack           # for reading in table files
from astropy.time import Time            # for reading in table files
from astropy.coordinates import Angle,SkyCoord,Galactic,BarycentricTrueEcliptic,match_coordinates_sky,search_around_sky
from astropy import units as u            # standard units
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.nist import Nist
from astroquery.xmatch import XMatch
#from astroquery.gaia import Gaia

# splat requirements
import splat
import splat.plot as splot
from splat.initialize import *
from splat.utilities import *
from splat.empirical import estimateDistance, typeToColor
from splat.core import classifyByStandard, searchLibrary
#from splat import DB_SOURCES, DB_SPECTRA
#import splat as spl

# Python 2->3 fix for input
try: input=raw_input
except NameError: pass

# set timeout limits to 1 minute
Simbad.TIMEOUT = 60
Vizier.TIMEOUT = 60
Nist.TIMEOUT = 60
XMatch.TIMEOUT = 180


##########################################################
###########   DATABASE QUERY, ACCESS & TOOLS   ###########
##########################################################

def prepDB(db_init,raCol='RA',decCol='DEC',desigCol='DESIGNATION',shortnameCol='SHORTNAME',coordinateCol='COORDINATES',force=False,verbose=False):
    '''
    Purpose
    -------

    Populates RA, DEC, DESIGNATION, COORDINATE and SHORTNAME columns if not present
    Requires RA and DEC, or DESIGNATION to be present
    '''
    db = copy.deepcopy(db_init)
    if raCol not in list(db.columns) or decCol not in list(db.columns): 
        if desigCol not in list(db.columns):
            raise ValueError('Database must have columns {} and {}, or {}'.format(raCol,decCol,desigCol))
        else:
            db[coordinateCol] = [designationToCoordinate(d) for d in db[desigCol]]
            if raCol not in list(db.columns) or decCol not in list(db.columns): 
                db[raCol] = [c.ra.degree for c in db[coordinateCol]]
                db[decCol] = [c.dec.degree for c in db[coordinateCol]]
    if desigCol not in list(db.columns):
        db[desigCol] = [coordinateToDesignation([db[raCol].iloc[i],db[decCol].iloc[i]]) for i in range(len(db))]
    if shortnameCol not in list(db.columns):
        db[shortnameCol] = [designationToShortName(d) for d in db[desigCol]]
# force COORDINATES, RA, DEC if desired
    if force == True:
        db[coordinateCol] = [designationToCoordinate(d) for d in db[desigCol]]
        db[raCol] = [c.ra.degree for c in db[coordinateCol]]
        db[decCol] = [c.dec.degree for c in db[coordinateCol]]
        db[shortnameCol] = [designationToShortName(d) for d in db[desigCol]]
# remove extra space in designation string
    db[desigCol] = [x.strip() for x in db[desigCol]]

    return db   

def longCoordList(dp,imax=10000):
    '''
    Purpose
    -------

    Creates a SkyCoord array of more than 10,000 coordinates
    '''
    for i in range(int(len(dp)/imax)):
        c0 = SkyCoord(ra=dp['RA'].iloc[(i*imax):((i+1)*imax)]*u.degree,dec=dp['DEC'].iloc[(i*imax):((i+1)*imax)]*u.degree)
        if i==0: c = c0
        else: c = c0.insert(0,c)
    c0 = SkyCoord(ra=dp['RA'].iloc[((i+1)*imax):]*u.degree,dec=dp['DEC'].iloc[((i+1)*imax):]*u.degree)
    c = c0.insert(0,c)
    return(c)
    
def catXMatch(dp1,dp2,sep=5*u.arcsec,imax = 10000,merge=False):
    '''
    Purpose
    -------

    Cross-matches two catalogs and return only those sources in common within the specified separation 
    '''
    sep = sep.to(u.arcsec)
    dp1p = prepDB(dp1)
    dp2p = prepDB(dp2)
    if len(dp1p)>imax: c1 = longCoordList(dp1p)
    else: c1 = SkyCoord(ra=dp1p['RA']*u.degree,dec=dp1p['DEC']*u.degree)
    if len(dp2p)>imax: c2 = longCoordList(dp2p)
    else: c2 = SkyCoord(ra=dp2p['RA']*u.degree,dec=dp2p['DEC']*u.degree)
    idx, sep2d, sep3d = match_coordinates_sky(c1,c2)
    dp1p['separation'] = [s.to(u.arcsec).value for s in sep2d]
    dp1pc = dp1p[dp1p['separation']<sep.value]
    idx, sep2d, sep3d = match_coordinates_sky(c2,c1)
    dp2p['separation'] = [s.to(u.arcsec).value for s in sep2d]
    dp2pc = dp2p[dp2p['separation']<sep.value]

    dp1pc.reset_index(inplace=True,drop=True)
    dp2pc.reset_index(inplace=True,drop=True)
    
    return dp1pc,dp2pc

def sourceXMatch(c,dp2,sep=5*u.arcsec,imax = 10000):
    '''
    Purpose
    -------

    Cross-matches a coordinate with a catalog and returns only those sources within the specified separation 
    '''
    sep = sep.to(u.arcsec)
    if not isinstance(c,SkyCoord): c = properCoordinates(c)
    dp2p = prepDB(dp2)
    if len(dp2p)>imax: c2 = longCoordList(dp2p)
    else: c2 = SkyCoord(ra=dp2p['RA']*u.degree,dec=dp2p['DEC']*u.degree)
    d2d = c.separation(c2).arcsec
    dp2p['separation'] = d2d
    dp2pc = dp2p[dp2p['separation']<sep.value]
    dp2pc.reset_index(inplace=True,drop=True)
    
    return dp2pc

def fetchDatabase(*args, **kwargs):
    '''
    Purpose
    -------

    Get the SpeX Database from either online repository or local drive
    '''
    filename = 'db_spexprism.txt'   # temporary original database file for backwards compatability
    if len(args) > 0:
        filename = args[0]
    kwargs['filename'] = kwargs.get('filename',filename)
    kwargs['filename'] = kwargs.get('file',kwargs['filename'])
    kwargs['folder'] = kwargs.get('folder',DB_FOLDER)
    url = kwargs.get('url',SPLAT_URL)+kwargs['folder']
    local = kwargs.get('local',True)
    online = kwargs.get('online',not local and checkOnline())
    local = not online
    kwargs['local'] = local
    kwargs['online'] = online
    kwargs['model'] = True

# determine format of file    
    delimiter = kwargs.get('delimiter','')
    fmt = kwargs.get('format','')
    fmt = kwargs.get('fmt',fmt)
    if delimiter == ',' or delimiter == 'comma' or delimiter == 'csv' or kwargs.get('comma',False) == True or ('.csv' in kwargs['filename']):
        delimiter = ','
        fmt = 'csv'
    if delimiter == '\t' or delimiter == 'tab' or kwargs.get('tab',False) == True or ('.txt' in kwargs['filename']):
        delimiter = '\t'
        fmt = 'tab'
    if fmt == '':
        raise NameError('\nCould not determine the file format of '+kwargs['filename']+'; please specify using format or delimiter keywords\n\n')


# check that folder/set is present either locally or online
# if not present locally but present online, switch to this mode
# if not present at either raise error
    folder = checkLocal(kwargs['folder'])
    if folder=='':
        folder = checkOnlineFile(kwargs['folder'])
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
#        print('Reading local')
        infile = checkLocal(kwargs['filename'])
        if infile=='':
            infile = checkLocal(os.path.join(kwargs['folder'],kwargs['filename']))
        if infile=='':
            raise NameError('\nCould not find '+kwargs['filename']+' locally\n\n')
        else:
            try:
                data = ascii.read(os.path.normpath(infile), delimiter=delimiter,fill_values='-99.',format=fmt)
#                data = ascii.read(infile, delimiter='\t',fill_values='-99.',format='tab')
            except:
                raise NameError('\nCould not load {}: this may be a decoding error\n'.format(infile))


# check if file is present; if so, read it in, otherwise go to interpolated
# online:
    if kwargs['online']:
#        print('Reading online')
        infile = checkOnlineFile(kwargs['filename'])
        if infile=='':
            infile = checkOnlineFile(os.path.join(kwargs['folder'],kwargs['filename']))
        if infile=='':
            raise NameError('\nCould not find '+kwargs['filename']+' on the SPLAT website\n\n')
        try:
#            open(os.path.basename(TMPFILENAME), 'wb').write(urllib2.urlopen(url+infile).read())
            open(os.path.basename(TMPFILENAME), 'wb').write(requests.get(url+infile).content)
            kwargs['filename'] = os.path.basename(tmp)
            data = ascii.read(os.path.basename(TMPFILENAME), delimiter=delimiter,fill_values='-99.',format=fmt)
            os.remove(os.path.basename(TMPFILENAME))
        except:
            raise NameError('\nHaving a problem reading in '+kwargs['filename']+' on the SPLAT website\n\n')

    return data


#####################################################
###########  ADDING NEW SPECTRA TO SPLAT   ##########
#####################################################

def addUserSpectra(folder='./',instrument='SPEX-PRISM',mode='update',repeat='retain',radius_repeat=10.*u.arcsec,input_file='input.txt',search_str='*.fits',sources_data_file=DB_SOURCES_FILE,spectra_data_file=DB_SPECTRA_FILE,verbose=True,*args):
    '''
    :Purpose:

        Adds in local spectral data to the underlying SPLAT library
        This program is currently UNDER DEVELOPMENT

    '''
# program constants
    optional_spectra_columns = ['PUBLISHED','DATA_BIBCODE','PROGRAM_PI','OBSERVATION_DATE','OBSERVATION_MJD','OBSERVATION_TIME','OBSERVER','AIRMASS']
    optional_sources_columns = ['NAME','DESIGNATION','RA','DEC','COORDINATES','DISCOVERY_REF','SPT','SPT_REF','SPT_OPT','SPT_OPT_REF','SPT_NIR','SPT_NIR_REF','SPT_LIT','SPT_LIT_REF','LUMINOSITY_CLASS','METALLICITY_CLASS','GRAVITY_CLASS_OPTICAL','GRAVITY_CLASS_OPTICAL_REF','GRAVITY_CLASS_NIR','GRAVITY_CLASS_NIR_REF','CLUSTER','CLUSTER_REF','BINARY','BINARY_TYPE','BINARY_REF','SBINARY','SBINARY_REF','COMPANION_NAME','COMPANION_REF']
    header_spectra_columns = {
        'OBSERVATION_DATE': ['OBS_DATE','OBS-DATE','UT-DATE'],
        'OBSERVATION_TIME': ['OBS_TIME','OBS-TIME','UT-TIME'],
        'OBSERVER': [],
        'AIRMASS': ['Z'],
        'SLIT': ['APERTURE'],
        'DISPERSER': ['GRATING','GRISM','DISPERSE'],
        }
    header_sources_columns = {
        'NAME': ['OBJECT','SOURCE','TARGET'],
        'RA': ['RA-D','RADEG'],
        'DEC': ['DEC-D','DECDEG'],
        }
    dataset_number_factor = 1e6
    now = time.localtime()
    nowstr = str(now.tm_year)+str(now.tm_mon)+str(now.tm_mday)

    if len(args) > 0:
        folder = args[0]
    if len(args) > 1:
        instrument = args[1]

##### STOPPED HERE #####

# check if this has already been read in
#    if folder in DATA_FOLDERS:
#        n = 

# check instrument
    inst = checkInstrument(instrument)
    if inst != False: instrument = inst

# check mode and repeat
    mode_labels = ['new','append','refresh','update']
    if mode.lower() not in mode_labels:
        if verbose==True: print('\nDo not recognize mode = {}; should be one of {}; reverting to update'.format(mode,mode_labels))
        mode = 'update'
    repeat_labels = ['replace','assert','retain','keep']
    if repeat.lower() not in repeat_labels:
        if verbose==True: print('\nDo not recognize repeat = {};  should be one of {}; reverting to retain'.format(repeat,repeat_labels))
        repeat = 'retain'

# check the folder is correctly specified
    if not os.path.exists(folder):
        print('\nCould not find folder {} in local directory structure; skipping')
        return

# check if spectra data file is present; if not, you'll need to generate a new one
    if spectra_data_file not in os.listdir(folder):
        if verbose == True: print('\nCannot find spectral data file {}; generating a new one from input files'.format(spectra_data_file))
        mode = 'new'

# STAGE 1: SET UP A NEW FOLDER OF DATA
    if mode.lower() == 'new':

# check if input file is in place; if not, make one
        if input_file not in os.listdir(folder):
            files = glob.glob(os.path.join(folder,search_str))
            files = [os.path.basename(f) for f in files]
            for f in [input_file,sources_data_file,spectra_data_file]:
                if f in files: files.remove(f)

            # turn into preliminary input.txt file
            input_db = pandas.DataFrame()
            input_db['DATA_FILE'] = files
            input_db['INSTRUMENT'] = [instrument]*len(files)
            if '.txt' in input_file: input_db.to_csv(os.path.join(folder,input_file),sep='\t',index=False)
            elif '.csv' in input_file: input_db.to_csv(os.path.join(folder,input_file),sep=',',index=False)
            elif '.xls' in input_file: input_db.to_excel(os.path.join(folder,input_file),index=False)
            else: raise ValueError('\nDo not recognize file format for {}'.format(input_file))

# prompt to continue?

# read in input file and start building spectral database
        if '.txt' in input_file: input_db = pandas.read_csv(os.path.join(folder,input_file),delimiter='\t')
        elif '.csv' in input_file: input_db = pandas.read_csv(os.path.join(folder,input_file),delimiter=',')
        elif '.xls' in input_file: input_db = pandas.read_excel(os.path.join(folder,input_file))
        else: raise ValueError('\nDo not recognize file format for input file {}'.format(input_file))

# capitalize all columns
        for c in list(input_db.columns):
            if c.upper() not in list(input_db.columns):
                input_db[c.upper()] = input_db[c]
                del input_db[c]

        # adjust instrument
        syn = ['INST','INSTR']
        if 'INSTRUMENT' not in list(input_db.columns): 
            for s in syn:
                if s in list(input_db.columns): 
                    input_db['INSTRUMENT'] = input_db[s]
                    del input_db[s]
        if 'INSTRUMENT' not in list(input_db.columns): 
            input_db['INSTRUMENT'] = [instrument]*len(input_db)
        for i,inst in enumerate(input_db['INSTRUMENT']):
            inst = checkInstrument(inst)
            if inst != False: input_db['INSTRUMENT'].iloc[i] = inst

        # adjust filename
        syn = ['FILE','FILENAME','FILE_NAME']
        if 'DATA_FILE' not in list(input_db.columns):
            for s in syn:
                if s in list(input_db.columns): 
                    input_db['DATA_FILE'] = input_db[s]
                    del input_db[s]

        # establish source and spectra data frames
        sources_db = pandas.DataFrame()
        spectra_db = pandas.DataFrame()

        # prep keys
        n = len(DATA_FOLDERS)
        keys = numpy.arange(len(input_db))+n*dataset_number_factor+1
        sources_db['SOURCE_KEY'] = [int(k) for k in keys]
        spectra_db['DATA_KEY'] = sources_db['SOURCE_KEY']
        spectra_db['SOURCE_KEY'] = sources_db['SOURCE_KEY']

        # required spectral information
        spectra_db['DATA_FILE'] = input_db['DATA_FILE']
        spectra_db['INSTRUMENT'] = input_db['INSTRUMENT']
        spectra_db['DATA_ENTRY'] = [nowstr]*len(input_db)

        # add in optional columns from input 
        for c in optional_spectra_columns:
            if c in list(input_db.columns): spectra_db[c] = input_db[c]
        for c in optional_sources_columns:
            if c in list(input_db.columns): sources_db[c] = input_db[c]
        for c in list(input_db.columns):
            if c not in optional_spectra_columns and c not in optional_sources_columns and c not in list(spectra_db.columns): spectra_db[c] = input_db[c]

# write out the source and spectra folders
        if '.txt' in sources_data_file: sources_db.to_csv(os.path.join(folder,sources_data_file),sep='\t',index=False)
        elif '.csv' in sources_data_file: sources_db.to_csv(os.path.join(folder,sources_data_file),sep=',',index=False)
        elif '.xls' in sources_data_file: sources_db.to_excel(os.path.join(folder,sources_data_file),index=False)
        else: raise ValueError('\nDo not recognize file format for {}'.format(sources_data_file))

        if '.txt' in spectra_data_file: spectra_db.to_csv(os.path.join(folder,spectra_data_file),sep='\t',index=False)
        elif '.csv' in spectra_data_file: spectra_db.to_csv(os.path.join(folder,spectra_data_file),sep=',',index=False)
        elif '.xls' in spectra_data_file: spectra_db.to_excel(os.path.join(folder,spectra_data_file),index=False)
        else: raise ValueError('\nDo not recognize file format for {}'.format(spectra_data_file))

# STAGE 1: SET UP A NEW FOLDER OF DATA
    if mode.lower() == 'new' or mode.lower() == 'append':
        pass

    return



#####################################################
###########   ACCESSING ONLINE CATALOGS   ###########
#####################################################


def getVizierName(catalog,output=False):
    '''
    Wrapper function for Vizier.find_catalogs to help search for Vizier catalog names
    '''
    catalog_list = Vizier.find_catalogs(catalog)
    output={}
    for k,v in catalog_list.items(): 
        print('{}: {}'.format(k,v.description))
        output[k] = v.description
    if output==True: return output
    else: return

def getVizierCatalog(catalog,catnum=0,verbose=True,limit=-1,return_pandas=True):
    '''
    Wrapper function for Vizier.get_catalogs which returns a single whole catalog
    '''
    if limit<0: Vizier.ROW_LIMIT = -1
    else: Vizier.ROW_LIMIT = int(limit)
    catalogs = Vizier.get_catalogs(catalog)
    Vizier.ROW_LIMIT = 50
    if len(catalogs)<1: 
        raise ValueError('Catalog {} is not in Vizier; use getVizierName() to check for catalog ID'.format(catalog))
    catnum = int(numpy.min([catnum,len(catalogs)]))
    if verbose==True: 
        print('{} catalog(s) identified, returning the {}th one: {}'.format(len(catalogs),catnum,list(catalogs.keys())[catnum]))
    if return_pandas==True:
        return catalogs[catnum].to_pandas()
    else:
        return catalogs[catnum]


def queryVizier(coordinate,**kwargs):
    '''
    see `splat.database.getPhotometry()`_

    .. _`splat.database.getPhotometry()` : api.html#splat.database.getPhotometry
    '''
    return getPhotometry(coordinate,**kwargs)

# NOTE: THIS IS NOT PROPERLY PASSING ON THE KEYWORDS

def getPhotometry(coordinate,return_pandas=True,catalog='2MASS',radius=30.*u.arcsec,sort='sep',limit=-1,info=False,nearest=False,verbose=False,**kwargs):
    '''
    Purpose
        Downloads photometry for a single source coordinate using astroquery.
        If you are getting data on multiple sources, it is preferable to use `splat.database.queryXMatch()`_

    .. _`splat.database.queryXMatch()` : api.html#splat.database.queryXMatch

    Required Inputs:
        :param: coordinate: Either an astropy SkyCoord or a variable that can be converted into a SkyCoord using `properCoordinates()`_

    .. _`properCoordinates()` : api.html#properCoordinates
        
    Optional Inputs:
        :param radius: Search radius, nominally in arcseconds although this can be changed by passing an astropy.unit quantity (default = 30 arcseconds)
        :param catalog: Catalog to query, which can be set to the Vizier catalog identifier code or to one of the following preset catalogs:

            * '2MASS' (or set ``2MASS``=True): the 2MASS All-Sky Catalog of Point Sources (`Cutri et al. 2003 <http://adsabs.harvard.edu/abs/2003yCat.2246....0C>`_), Vizier id II/246
            * 'SDSS' (or set ``SDSS``=True): the The SDSS Photometric Catalog, Release 9 (`Adelman-McCarthy et al. 2012 <http://adsabs.harvard.edu/abs/2012ApJS..203...21A>`_), Vizier id V/139
            * 'WISE' (or set ``WISE``=True): the WISE All-Sky Data Release (`Cutri et al. 2012 <http://adsabs.harvard.edu/abs/2012yCat.2311....0C>`_), Vizier id II/311
            * 'ALLWISE' (or set ``ALLWISE``=True): the AllWISE Data Release (`Cutri et al. 2014 <http://adsabs.harvard.edu/abs/2014yCat.2328....0C>`_), Vizier id II/328
            * 'VISTA' (or set ``VISTA``=True): the VIKING catalogue data release 1 (`Edge et al. 2013 <http://adsabs.harvard.edu/abs/2013Msngr.154...32E>`_), Vizier id II/329
            * 'CFHTLAS' (or set ``CFHTLAS``=True): the CFHTLS Survey (T0007 release) by (`Hudelot et al. 2012 <http://adsabs.harvard.edu/abs/2012yCat.2317....0H>`_), Vizier id II/317
            * 'DENIS' (or set ``DENIS``=True): the DENIS DR3 (DENIS Consortium 2005), Vizier id B/denis/denis
            * 'UKIDSS' (or set ``UKIDSS``=True): the UKIDSS-DR8 LAS, GCS and DXS Surveys (`Lawrence et al. 2012 <http://adsabs.harvard.edu/abs/2007MNRAS.379.1599L>`_), Vizier id II/314
            * 'LEHPM' (or set ``LEHPM``=True): the Liverpool-Edinburgh High Proper Motion Catalogue (`Pokorny et al. 2004 <http://adsabs.harvard.edu/abs/2004A&A...421..763P>`_), Vizier id J/A+A/421/763
            * 'SIPS' (or set ``SIPS``=True): the Southern Infrared Proper Motion Survey (`Deacon et al 2005 <http://adsabs.harvard.edu/abs/2005A&A...435..363D>`_), Vizier id J/A+A/435/363
            * 'UCAC4' (or set ``UCAC4``=True): the UCAC4 Catalogue (`Zacharias et al. 2012 <http://adsabs.harvard.edu/abs/2012yCat.1322....0Z>`_), Vizier id I/322A
            * 'USNOB' (or set ``USNO``=True): the USNO-B1.0 Catalog (`Monet et al. 2003 <http://adsabs.harvard.edu/abs/2003AJ....125..984M>`_), Vizier id I/284
            * 'LSPM' (or set ``LSPM``=True): the LSPM-North Catalog (`Lepine et al. 2005 <http://adsabs.harvard.edu/abs/2005AJ....129.1483L>`_), Vizier id I/298
            * 'GAIA-DR1': the GAIA DR1 Catalog (`Gaia Collaboration et al. 2016 <http://adsabs.harvard.edu/abs/2016yCat.1337....0G>`_), Vizier id I/337
            * 'GAIA' or 'GAIA-DR2' (or set ``GAIA``=True): the GAIA DR2 Catalog (REF TBD), Vizier id I/345/gaia2

        :param: sort: String specifying the parameter to sort the returned SIMBAD table by; by default this is the offset from the input coordinate (default = 'sep')
        :param: nearest: Set to True to return on the single nearest source to coordinate (default = False)
        :param: return_pandas: Return a pandas table as opposed to an astropy Table (default = True)
        :param: verbose: Give feedback (default = False)

    Output:
        An astropy or pandas Table that contains data from the Vizier query, or a blank Table if no sources are found

    Example:

    >>> import splat
    >>> import splat.database as spdb
    >>> from astropy import units as u
    >>> c = splat.properCoordinates('J053625-064302')
    >>> v = spdb.querySimbad(c,catalog='SDSS',radius=15.*u.arcsec)
    >>> print(v)
      _r    _RAJ2000   _DEJ2000  mode q_mode  cl ... r_E_ g_J_ r_F_ i_N_  sep  
     arcs     deg        deg                     ... mag  mag  mag  mag   arcs 
    ------ ---------- ---------- ---- ------ --- ... ---- ---- ---- ---- ------
     7.860  84.105967  -6.715966    1          3 ...   --   --   --   --  7.860
    14.088  84.108113  -6.717206    1          6 ...   --   --   --   -- 14.088
    14.283  84.102528  -6.720843    1      +   6 ...   --   --   --   -- 14.283
    16.784  84.099524  -6.717878    1          3 ...   --   --   --   -- 16.784
    22.309  84.097988  -6.718049    1      +   6 ...   --   --   --   -- 22.309
    23.843  84.100079  -6.711999    1      +   6 ...   --   --   --   -- 23.843
    27.022  84.107504  -6.723965    1      +   3 ...   --   --   --   -- 27.022

    '''

# check if online
    if not checkOnline():
        print('\nYou are currently not online; cannot do a Vizier query')
        return Table()

    VIZIER_REF = {
        'SDSS': {'altname': [], 'catalog': u'V/147/sdss12'},
        '2MASS': {'altname': [], 'catalog': u'II/246/out'},
        'USNO': {'altname': ['USNOB','USNO-B','USNOB1.0','USNO-B1.0'], 'catalog': u'I/284/out'},
        'LSPM': {'altname': ['LSPM-N','LSPM-NORTH'], 'catalog': u'I/298/lspm_n'},
        'WISE': {'altname': [], 'catalog': u'II/311/wise'},
        'ALLWISE': {'altname': [], 'catalog': u'II/328/allwise'},
        'CATWISE': {'altname': [], 'catalog': u'II/365/catwise'},
        'UKIDSS': {'altname': [], 'catalog': u'II/314'},
        'CFHT': {'altname': ['CFHTLAS'], 'catalog': u'II/317/sample'},
        'UCAC': {'altname': [], 'catalog': u'I/322A/out'},
        'VISTA': {'altname': [], 'catalog': u'II/329/urat1'},
        'GAIA-DR1': {'altname': ['GAIA1','GAIADR1'], 'catalog': u'II/337/gaia'},
        'GAIA-DR2': {'altname': ['GAIA2','GAIADR2'], 'catalog': u'I/345/gaia2'},
        'GAIA-EDR3': {'altname': ['GAIA','GAIA3','GAIADR3'], 'catalog': u'I/350/gaiaedr3'},
        'PANSTARRS': {'altname': ['PAN-STARRS','PS1'], 'catalog': u'II/349/ps1'},
        'DENIS': {'altname': [], 'catalog': u'B/denis'},
        'LEHPM': {'altname': [], 'catalog': u'J/A+A/421/763'},
        'LEPINE': {'altname': ['LEPINE-MDWARFS'], 'catalog': u'J/AJ/142/138/Mdwarfs'},
        'SIPS': {'altname': [], 'catalog': u'J/A+A/435/363'},
        'MOVERS': {'altname': [], 'catalog': u'J/AJ/151/41'},
        'LATEMOVERS': {'altname': ['LATE-MOVERS'], 'catalog': u'J/AJ/153/92'},
        'GLIESE': {'altname': ['GJ'], 'catalog': u'J/PASP/122/885/table1'},
        'DESHPANDE2013': {'altname': ['DESHPANDE-2013','APOGEE_UCD'], 'catalog': u'J/AJ/146/156/table1'},
        'DITTMAN2014': {'altname': ['DITTMAN-2014','DITTMAN-PARALLAX','DIT16'], 'catalog': u'J/ApJ/784/156/table2'},
        'NEWTON2016': {'altname': ['NEWTON-2016','NEW16'], 'catalog': u'J/ApJ/821/93/table1'},
        'KIRKPATRICK2016': {'altname': ['KIRKPATRICK-2016','ALLWISE-MOTION','KIR16'], 'catalog': u'J/ApJS/224/36/motionobj'},
    }

# give a summary of the built-in catalogs
    if info==True:
        print('Currently available input catalogs:')
        for k in list(VIZIER_REF.keys()):
            line = '\t{}: '.format(k)
            if len(VIZIER_REF[k]['altname'])>0:
                line=line+'(or'
                for a in VIZIER_REF[k]['altname']: line=line+' {}'.format(a)
                line=line+') '
            print(line+'Vizier reference: {}'.format(str(VIZIER_REF[k]['catalog'])))
            catsp = str(VIZIER_REF[k]['catalog']).split('/')
            ctref = catsp[0] 
            for ct in catsp[1:-1]: ctref=ctref+'/'+ct
            print('\tURL = https://cdsarc.unistra.fr/viz-bin/cat/{}\n'.format(ctref))
        return

    for c in list(VIZIER_REF.keys()): 
        if kwargs.get(c,False): catalog = c

# is catalog one of pre-defined ones?
    for c in list(VIZIER_REF.keys()): 
        if kwargs.get(c,False): catalog = c

    cat = checkDict(catalog,VIZIER_REF)
    if cat == False: cat = catalog
    else: cat = VIZIER_REF[cat]['catalog']

# parameters
    if not isUnit(radius): radius = radius*u.arcsec

# convert coordinate if necessary
    if not isinstance(coordinate,SkyCoord):
        try:
            c = properCoordinates(coordinate)
        except:
            print('\n{} is not a proper coordinate'.format(coordinate))
            return numpy.nan
    else:
        c = copy.deepcopy(coordinate)

# search Vizier, sort by separation        
    v = Vizier(columns=["**", "+_r"], catalog=cat)
    if limit<0: v.ROW_LIMIT = -1
    else: v.ROW_LIMIT = int(limit)
    t_vizier = v.query_region(c,radius=radius)
    tv = Table()
    if len(t_vizier) > 0:
        for k in list(t_vizier.keys()):
            if cat in k: tv = t_vizier[k]
    else:
        tv = Table()

    if len(tv)==0: 
        if return_pandas==True: return pandas.DataFrame()
        else: return tv

# sorting
    tv['sep'] = tv['_r']
    if len(tv) > 1:
        sortparam = kwargs.get('sort','sep')
        if sortparam in list(tv.keys()):
            tv.sort(sortparam)
        else:
            if verbose:
                print('\nCannot find sorting keyword {}; try using {}\n'.format(sort,list(tv.keys())))

# return only nearest
#    print(kwargs.get('nearest',False))
    if nearest == True:
#        tv = tv[0]
        while len(tv) > 1:
            tv.remove_row(1)
#        print(tv)

# reformat to convert binary ascii data to text
    for s in list(tv.keys()):
        if isinstance(tv[s][0],bytes) == True or isinstance(tv[s][0],numpy.bytes_)  == True:
            tmp = [x.decode() for x in tv[s]]
            tv.remove_column(s)
            tv[s] = tmp

# convert to pandas if desired
    if return_pandas==True:
        tv = tv.to_pandas()
        fix = list(tv.dtypes[tv.dtypes=='object'].keys())
        if len(fix) > 0: 
            for f in fix:
                tv[f] = tv[f].str.decode('utf8')

    return tv



def querySimbad(variable,radius=30.*u.arcsec,sort='sep',reject_type='',nearest=False,iscoordinate=False,isname=False,clean=False,return_pandas=True,verbose=False,**kwargs):
    '''
    Purpose
        Queries SIMBAD using astroquery for a single source
        If you are getting data on multiple sources, it is preferable to use `splat.database.queryXMatch()`_

    Required Inputs:
        :param: variable: Either an astropy SkyCoord object containing position of a source, a variable that can be converted into a SkyCoord using `spl.properCoordinates()`_, or a string name for a source.
        
    Optional Inputs:
        :param: radius: Search radius, nominally in arcseconds although can be set by assigning and astropy.unit value (default = 30 arcseconds)
        :param: sort: String specifying the parameter to sort the returned SIMBAD table by; by default this is the offset from the input coordinate (default = 'sep')
        :param: reject_type: Set to string or list of strings to filter out object types not desired. Useful for crowded fields (default = None)
        :param: nearest: Set to True to return on the single nearest source to coordinate (default = False)
        :param: iscoordinate: Specifies that input is a coordinate of some kind (default = False)
        :param: isname: Specifies that input is a name of some kind (default = False)
        :param: clean: Set to True to clean the SIMBAD output and reassign to a predefined set of parameters (default = True)
        :param: return_pandas: Return a pandas table as opposed to an astropy Table (default = True)
        :param: verbose: Give lots of feedback (default = False)

    Output:
        An astropy or pandas Table that contains data from the SIMBAD search, or a blank Table if no sources found

    Example:

    >>> import splat
    >>> import splat.database as spdb
    >>> from astropy import units as u
    >>> c = splat.properCoordinates('J053625-064302')
    >>> q = spdb.querySimbad(c,radius=15.*u.arcsec,reject_type='**')
    >>> print(q)
              NAME          OBJECT_TYPE     OFFSET    ... K_2MASS K_2MASS_E
    ----------------------- ----------- ------------- ... ------- ---------
               BD-06  1253B        Star  4.8443894429 ...                  
                [SST2010] 3        Star 5.74624887682 ...   18.36       0.1
                BD-06  1253         Ae* 7.74205447776 ...   5.947     0.024
               BD-06  1253A          ** 7.75783861347 ...                  
    2MASS J05362590-0643020     brownD* 13.4818185612 ...  12.772     0.026
    2MASS J05362577-0642541        Star  13.983717577 ...                  


    .. _`splat.database.queryXMatch()` : api.html#splat.database.queryXMatch
    .. _`spl.properCoordinates()` : api.html#spl.properCoordinates
    '''

# check that online
    if not checkOnline():
        print('\nYou are currently not online; cannot do a SIMBAD query')
        return Table()

# parameters 
    if not isUnit(radius): radius=radius*u.arcsec

# check if this is a coordinate query
    if isinstance(variable,SkyCoord):
        c = copy.deepcopy(variable)
        iscoordinate = True
    elif not isname:
        try:
            c = properCoordinates(variable)
            iscoordinate = True
# this is probably a name
        except:
            isname = True
    else:
        if isinstance(variable,bytes):
            c = variable.decode()
        else:
            c = str(variable)

# prep Simbad search
    sb = Simbad()
    votfields = ['otype','parallax','sptype','propermotions','rot','rvz_radvel','rvz_error',\
    'rvz_bibcode','fluxdata(B)','fluxdata(V)','fluxdata(R)','fluxdata(I)','fluxdata(g)','fluxdata(r)',\
    'fluxdata(i)','fluxdata(z)','fluxdata(J)','fluxdata(H)','fluxdata(K)']
    for v in votfields:
        sb.add_votable_fields(v)

# search SIMBAD by coordinate
    if iscoordinate:
        t_sim = sb.query_region(c,radius=radius)
        if not isinstance(t_sim,Table):
            if verbose:
                print('\nNo sources found; returning empty Table\n')
            return Table()

# if more than one source, sort the results by separation
        sep = [c.separation(SkyCoord(str(t_sim['RA'][lp]),str(t_sim['DEC'][lp]),unit=(u.hourangle,u.degree))).arcsecond for lp in numpy.arange(len(t_sim))]
        t_sim['sep'] = sep

# search SIMBAD by name
    elif isname:
        t_sim = sb.query_object(c)
        if not isinstance(t_sim,Table):
            if verbose:
                print('\nNo sources found; returning empty Table\n')
            return Table()
        t_sim['sep'] = numpy.zeros(len(t_sim['RA']))

    else:
        raise ValueError('problem!')

# sort results by separation by default
    if sort in list(t_sim.keys()):
        t_sim.sort(sort)
    else:
        if verbose:
            print('\nCannot sort by {}; try keywords {}\n'.format(sort,list(t_sim.keys())))


# reject object types not wanted
    if reject_type != '':
        rej = reject_type.split(',')
        for r in rej:
            w = numpy.array([str(r) not in str(o) for o in t_sim['OTYPE']])
            if len(w) > 0:
                t_sim = t_sim[w]

# trim to single source if nearest flag is set
    if iscoordinate and nearest==True:
        while len(t_sim)>1:
            t_sim.remove_row(1) 

# clean up the columns    
    if clean == True and len(t_sim) > 0:
        t_src = Table()

# reformat to convert binary ascii data to text
        for s in list(t_sim.keys()):
            if isinstance(t_sim[s][0],bytes) == True or isinstance(t_sim[s][0],numpy.bytes_)  == True:
                tmp = [x.decode() for x in t_sim[s]]
                t_sim.remove_column(s)
                t_sim[s] = tmp

#        if not isinstance(t_sim['MAIN_ID'][0],str):
        t_src['NAME'] = [x.replace('  ',' ') for x in t_sim['MAIN_ID']]
#        else: 
#            t_src['NAME'] = t_sim['MAIN_ID']
#        if not isinstance(t_sim['OTYPE'][0],str):
        t_src['OBJECT_TYPE'] = [x.replace('  ',' ') for x in t_sim['OTYPE']]
#        else:
#            t_src['OBJECT_TYPE'] = t_sim['OTYPE']
        t_src['OFFSET'] = t_sim['sep']
#        if not isinstance(t_sim['SP_TYPE'][0],str):
        t_src['LIT_SPT'] = [x.replace(' ','') for x in t_sim['SP_TYPE']]
#        else:
#            t_src['LIT_SPT'] = t_sim['SP_TYPE']
#        if not isinstance(t_sim['SP_BIBCODE'][0],str):
        t_src['LIT_SPT_REF'] = [x.replace(' ','') for x in t_sim['SP_BIBCODE']]
#        else: 
#            t_src['LIT_SPT_REF'] = t_sim['SP_BIBCODE']
        t_src['DESIGNATION'] = ['J{}{}'.format(t_sim['RA'][i],t_sim['DEC'][i]).replace(' ','').replace('.','') for i in range(len(t_sim))] 
        t_src['RA'] = numpy.zeros(len(t_sim))
        t_src['DEC'] = numpy.zeros(len(t_sim))
        for i in range(len(t_sim)):
            c2 = properCoordinates(t_src['DESIGNATION'][i])
            t_src['RA'][i] = c2.ra.value
            t_src['DEC'][i] = c2.dec.value
        t_src['PARALLAX'] = [str(p).replace('--','') for p in t_sim['PLX_VALUE']]
        t_src['PARALLAX_E'] = [str(p).replace('--','') for p in t_sim['PLX_ERROR']]
#        if not isinstance(t_sim['PLX_BIBCODE'][0],str):
        t_src['PARALLEX_REF'] = [x.replace(' ','') for x in t_sim['PLX_BIBCODE']]
#        else:
#            t_src['PARALLEX_REF'] = t_sim['PLX_BIBCODE']
        t_src['MU_RA'] = [str(p).replace('--','') for p in t_sim['PMRA']]
        t_src['MU_DEC'] = [str(p).replace('--','') for p in t_sim['PMDEC']]
        t_src['MU'] = numpy.zeros(len(t_sim))
        for i in range(len(t_sim)):
            if t_src['MU_RA'][i] != '':
                t_src['MU'][i] = (float(t_src['MU_RA'][i])**2+float(t_src['MU_DEC'][i])**2)**0.5
        t_src['MU_E'] = [str(p).replace('--','') for p in t_sim['PM_ERR_MAJA']]
#        if not isinstance(t_sim['PM_BIBCODE'][0],str):
        t_src['MU_REF'] = [x.replace(' ','') for x in t_sim['PM_BIBCODE']]
#        else:
#            t_src['MU_REF'] = t_sim['PM_BIBCODE']
        t_src['RV'] = [str(p).replace('--','') for p in t_sim['RVZ_RADVEL']]
        t_src['RV_E'] = [str(p).replace('--','') for p in t_sim['RVZ_ERROR']]
#        if not isinstance(t_sim['RVZ_BIBCODE'][0],str):
        t_src['RV_REF'] = [x.replace(' ','') for x in t_sim['RVZ_BIBCODE']]
#        else:
#            t_src['RV_REF'] = t_sim['RVZ_BIBCODE']
        t_src['VSINI'] = [str(p).replace('--','') for p in t_sim['ROT_Vsini']]
        t_src['VSINI_E'] = [str(p).replace('--','') for p in t_sim['ROT_err']]
#        if not isinstance(t_sim['ROT_bibcode'][0],str):
        t_src['VSINI_REF'] = [x.replace(' ','') for x in t_sim['ROT_bibcode']]
#        else:
#            t_src['VSINI_REF'] = t_sim['ROT_bibcode']
        t_src['J_2MASS'] = [str(p).replace('--','') for p in t_sim['FLUX_J']]
        t_src['J_2MASS_E'] = [str(p).replace('--','') for p in t_sim['FLUX_ERROR_J']]
        t_src['H_2MASS'] = [str(p).replace('--','') for p in t_sim['FLUX_H']]
        t_src['H_2MASS_E'] = [str(p).replace('--','') for p in t_sim['FLUX_ERROR_H']]
        t_src['K_2MASS'] = [str(p).replace('--','') for p in t_sim['FLUX_K']]
        t_src['K_2MASS_E'] = [str(p).replace('--','') for p in t_sim['FLUX_ERROR_K']]
    else:
        t_src = t_sim.copy()

# convert to pandas if desired
    if return_pandas==True:
        t_src = t_src.to_pandas()
#        fix = list(t_src.dtypes[t_src.dtypes=='object'].keys())
#        if len(fix) > 0: 
#            for f in fix:
#                t_src[f] = t_src[f].str.decode('utf8')

    return t_src



def _querySimbad2(t_src,designation='DESIGNATION',**kwargs):
    '''
    Purpose
        Internal function that queries Simbad and populates data for source table.

    :Note:
        **this program is in beta testing; bugs/errors are likely**

    :Required parameters:
        :param table: an astropy Table object, requires the presence of DESIGNATION column

    :Optional parameters:
        :param simbad_radius = 30 arcseconds: circular radius to search for sources (note: must be an angular quantity)
        :param export = '': filename to which to export resulting table to; if equal to a null string then no expoer is made. Note that a populated table is returned in either case
        :param closest = False: return only the closest source to given coordinate
    '''    

# parameters 
    simbad_radius = kwargs.get('simbad_radius',30.*u.arcsec)
    verbose = kwargs.get('verbose',True)
# checks
    if designation not in t_src.keys():
        raise NameError('\nDesigation column {} is required for input table to querySimbad\n'.format(designation))
    if 'SIMBAD_SEP' not in t_src.keys():
        t_src['SIMBAD_SEP'] = Column(numpy.zeros(len(t_src)),dtype='float')
# must be online
    if not checkOnline():
        print('\nYou are currently not online so cannot query Simbad\n')
        return t_src

# if necessary, populate columns that are expected for source database
    for c in list(DB_SOURCES.keys()):
        if c not in t_src.keys():
            t_src[c] = Column([' '*50 for des in t_src['DESIGNATION']],dtype='str')

# prep Simbad search
    sb = Simbad()
    votfields = ['otype','parallax','sptype','propermotions','rot','rvz_radvel','rvz_error',\
    'rvz_bibcode','fluxdata(B)','fluxdata(V)','fluxdata(R)','fluxdata(I)','fluxdata(g)','fluxdata(r)',\
    'fluxdata(i)','fluxdata(z)','fluxdata(J)','fluxdata(H)','fluxdata(K)']
    for v in votfields:
        sb.add_votable_fields(v)

# search by source
    for i,des in enumerate(t_src['DESIGNATION']):
        print(i,des)
        c = designationToCoordinate(des)
        try:
            t_sim = sb.query_region(c,radius=simbad_radius)
        except:
            t_sim = None
# source found in query
        if isinstance(t_sim,Table):
# many sources found
#            if len(t_sim) >= 1:      # take the closest position
            if verbose:
                print('\nSource {} Designation = {} {} match(es)'.format(i+1,des,len(t_sim)))
                print(t_sim)

            sep = [c.separation(SkyCoord(str(t_sim['RA'][lp]),str(t_sim['DEC'][lp]),unit=(u.hourangle,u.degree))).arcsecond for lp in numpy.arange(len(t_sim))]
            t_sim['sep'] = sep
            t_sim.sort('sep')
            if len(t_sim) > 1:
                while len(t_sim)>1:
                    t_sim.remove_row(1) 
# one source found
#            else:
#                t_sim['sep'] = [c.separation(SkyCoord(str(t_sim['RA'][0]),str(t_sim['DEC'][0]),unit=(u.hourangle,u.degree))).arcsecond]

# fill in information
            t_src['SIMBAD_NAME'][i] = t_sim['MAIN_ID'][0]
            t_src['NAME'][i] = t_src['SIMBAD_NAME'][i]
            t_src['SIMBAD_OTYPE'][i] = t_sim['OTYPE'][0]
            if not isinstance(t_sim['SP_TYPE'][0],str):
                t_sim['SP_TYPE'][0] = t_sim['SP_TYPE'][0].decode()
            spt = t_sim['SP_TYPE'][0]
            spt.replace(' ','').replace('--','')
            t_src['SIMBAD_SPT'][i] = spt
            t_src['SIMBAD_SPT_REF'][i] = t_sim['SP_BIBCODE'][0]
            t_src['SIMBAD_SEP'][i] = t_sim['sep'][0]
            if spt != '':
                t_src['LIT_TYPE'][i] = t_src['SIMBAD_SPT'][i]
                t_src['LIT_TYPE_REF'][i] = t_src['SIMBAD_SPT_REF'][i]
            t_src['DESIGNATION'][i] = 'J{}{}'.format(t_sim['RA'][0],t_sim['DEC'][0]).replace(' ','').replace('.','')
            coord = properCoordinates(t_src['DESIGNATION'][i])
            t_src['RA'][i] = coord.ra.value
            t_src['DEC'][i] = coord.dec.value
            t_src['OBJECT_TYPE'][i] = 'VLM'
            if 'I' in t_sim['SP_TYPE'][0] and 'V' not in t_sim['SP_TYPE'][0]:
                t_src['LUMINOSITY_CLASS'][i] = 'I{}'.format(t_sim['SP_TYPE'][0].split('I',1)[1])
                t_src['OBJECT_TYPE'][i] = 'GIANT'
            if 'VI' in t_sim['SP_TYPE'][0] or 'sd' in t_sim['SP_TYPE'][0]:
                t_src['METALLICITY_CLASS'][i] = '{}sd'.format(t_sim['SP_TYPE'][0].split('sd',1)[0])
            t_src['PARALLAX'][i] = str(t_sim['PLX_VALUE'][0]).replace('--','')
            t_src['PARALLAX_E'][i] = str(t_sim['PLX_ERROR'][0]).replace('--','')
            if isinstance(t_sim['PLX_BIBCODE'][0],str):
                t_src['PARALLEX_REF'][i] = str(t_sim['PLX_BIBCODE'][0]).replace('--','')
            else:
                t_src['PARALLEX_REF'][i] = t_sim['PLX_BIBCODE'][0].decode()
            t_src['MU_RA'][i] = str(t_sim['PMRA'][0]).replace('--','')
            t_src['MU_DEC'][i] = str(t_sim['PMDEC'][0]).replace('--','')
#                try:            # this is in case MU is not present
            t_src['MU'][i] = (float('{}0'.format(t_src['MU_RA'][i]))**2+float('{}0'.format(t_src['MU_DEC'][i]))**2)**0.5
            t_src['MU_E'][i] = str(t_sim['PM_ERR_MAJA'][0]).replace('--','')
#                except:
#                    pass
            t_src['MU_REF'][i] = t_sim['PM_BIBCODE'][0]
            t_src['RV'][i] = str(t_sim['RVZ_RADVEL'][0]).replace('--','')
            t_src['RV_E'][i] = str(t_sim['RVZ_ERROR'][0]).replace('--','')
            t_src['RV_REF'][i] = t_sim['RVZ_BIBCODE'][0]
            t_src['VSINI'][i] = str(t_sim['ROT_Vsini'][0]).replace('--','')
            t_src['VSINI_E'][i] = str(t_sim['ROT_err'][0]).replace('--','')
            t_src['VSINI_REF'][i] = t_sim['ROT_bibcode'][0]
            if isinstance(t_sim['FLUX_J'][0],str):
                t_src['J_2MASS'][i] = t_sim['FLUX_J'][0].replace('--','')
            else:
                t_src['J_2MASS'][i] = t_sim['FLUX_J'][0]
            if isinstance(t_sim['FLUX_ERROR_J'][0],str):
                t_src['J_2MASS_E'][i] = t_sim['FLUX_ERROR_J'][0].replace('--','')
            else:
                t_src['J_2MASS_E'][i] = t_sim['FLUX_ERROR_J'][0]
            if isinstance(t_sim['FLUX_H'][0],str):
                t_src['H_2MASS'][i] = t_sim['FLUX_H'][0].replace('--','')
            else:
                t_src['H_2MASS'][i] = t_sim['FLUX_H'][0]
            if isinstance(t_sim['FLUX_ERROR_H'][0],str):
                t_src['H_2MASS_E'][i] = t_sim['FLUX_ERROR_H'][0].replace('--','')
            else:
                t_src['H_2MASS_E'][i] = t_sim['FLUX_ERROR_H'][0]
            if isinstance(t_sim['FLUX_K'][0],str):
                t_src['KS_2MASS'][i] = t_sim['FLUX_K'][0].replace('--','')
            else:
                t_src['KS_2MASS'][i] = t_sim['FLUX_K'][0]
            if isinstance(t_sim['FLUX_ERROR_K'][0],str):
                t_src['KS_2MASS_E'][i] = t_sim['FLUX_ERROR_K'][0].replace('--','')
            else:
                t_src['KS_2MASS_E'][i] = t_sim['FLUX_ERROR_K'][0]

    return


# query the NIST database

def queryNist(element,wave_range,clean=['Observed'],noclean=False,verbose=True,wavelength_type='vacuum'):
# check inputs
    if not isinstance(element,str):
        raise ValueError('\nElement input must be a string like "K I", not {}'.format(element))
    if len(element.strip().split(' ')) == 1:
        element = element+' I'
    if len(element.strip().split(' ')) != 2:
        raise ValueError('\nElement input must be a string like "K I", not {}'.format(element))
    if not isUnit(wave_range[0]): wave_range = [w*u.micron for w in wave_range]  

    t = Nist.query(wave_range[0],wave_range[1],linename=element,energy_level_unit='eV',wavelength_type=wavelength_type)
    if noclean == False:
        for m in clean:
            t = t[~t[m].mask]
    if len(t) == 0 and verbose == True: print('\nNo lines found; check element, wavelength range, or set noclean=True')
    return(t)



def queryXMatch(db,radius=30.*u.arcsec,catalog='2MASS',file='',desigCol='DESIGNATION',raCol='RA',decCol='DEC',verbose=False,clean=True,drop_repeats=True,use_select_columns=False,select_columns=[],prefix=None,info=False,debug=False,*args):
    '''
    Purpose
        Queries databases in the XXX XMatch service (REF), including SIMBAD
        This is the preferred manner for extracting data for large numbers of sources

    Required Inputs:
        :param: db: a pandas Dataframe (FUTURE: astropy Table, dict, or file name for csv, txt or xls file). 
        This must contain column(s) for designation (specified in `desigCol`) and/or RA (specified in `raCol`) and DEC (specified in `decCol`)

    .. _`spl.properCoordinates()` : api.html#spl.properCoordinates
        
    Optional Inputs:
        :param radius: Search radius, nominally in arcseconds although can be set by assigning and astropy.unit value (default = 30 arcseconds)
        :param desigCol: column in db that specifies source designations ('Jhhmmss[.]s±ddmmss[.]s')
        :param raCol: column in db that specifies source RAs (in degrees)
        :param decCol: column in db that specifies source DECs (in degrees)
        :param catalog: Database to query, which can be set one of the follow presets or any catalog listed in astroquery.xmatch.XMatch.get_available_tables():

            * 'SIMBAD' (or set ``SIMBAD``=True): query SIMBAD (coordinate search only)
            * '2MASS' (or set ``2MASS``=True): query the 2MASS All-Sky Catalog of Point Sources (`Cutri et al. 2003 <http://adsabs.harvard.edu/abs/2003yCat.2246....0C>`_), Vizier id II/246
            * 'SDSS' (or set ``SDSS``=True): query the SDSS Photometric Catalog, Release 12 (NEED REF), Vizier id V/147
            * 'SDSS9' (or set ``SDSS``=True): query the SDSS Photometric Catalog, Release 9 (`Adelman-McCarthy et al. 2012 <http://adsabs.harvard.edu/abs/2012ApJS..203...21A>`_), Vizier id V/147
            * 'ALLWISE' (or set ``ALLWISE``=True): the AllWISE Data Release (`Cutri et al. 2014 <http://adsabs.harvard.edu/abs/2014yCat.2328....0C>`_), Vizier id II/328
            * 'DENIS' (or set ``DENIS``=True): the DENIS DR3 (DENIS Consortium 2005), Vizier id B/denis/denis
            * 'GAIA-DR1': the GAIA DR1 Catalog (`Gaia Collaboration et al. 2016 <http://adsabs.harvard.edu/abs/2016yCat.1337....0G>`_), Vizier id I/337
            * 'GAIA' or 'GAIA-DR2' (or set ``GAIA``=True): the GAIA DR2 Catalog (REF TBD), Vizier id I/345/gaia2, accessed using astroquery.gaia
            
        :param nearest: Set to True to return only the single nearest source to each coordinate (default = True)
        :param clean: Set to True to clean the SIMBAD output and reassign to a predefined set of parameters (default = True)
        :param file: Write the output to a csv or xlsx file (default = '' or not saved)
        :param verbose: Give lots of feedback (default = False)

        :param sort: String specifying the parameter to sort the returned SIMBAD table by; by default this is the offset from the input coordinate (default = 'sep')
        :param return_pandas: Return a pandas table as opposed to an astropy Table (default = True)
        :param reject_type: Set to string or list of strings to filter out object types not desired. Useful for crowded fields (default = None)

    Output:
        A pandas Dataframe that contains data from the search, or a blank frame if no sources found

    Example:

    >>> import splat
    >>> from astropy import units as u
    >>> c = spl.properCoordinates('J053625-064302')
    >>> q = spl.querySimbad(c,radius=15.*u.arcsec,reject_type='**')
    >>> print(q)
              NAME          OBJECT_TYPE     OFFSET    ... K_2MASS K_2MASS_E
    ----------------------- ----------- ------------- ... ------- ---------
               BD-06  1253B        Star  4.8443894429 ...                  
                [SST2010] 3        Star 5.74624887682 ...   18.36       0.1
                BD-06  1253         Ae* 7.74205447776 ...   5.947     0.024
               BD-06  1253A          ** 7.75783861347 ...                  
    2MASS J05362590-0643020     brownD* 13.4818185612 ...  12.772     0.026
    2MASS J05362577-0642541        Star  13.983717577 ...                  

    '''
    callloop = 5

# pre-defined catalogs
    XMATCH_CATALOGS = {
        'SIMBAD': {'altname': [],'vref': u'simbad', 'select_columns': ['main_id','ra','dec','main_type','sp_type','plx','pmra','pmdec','radvel','B', 'V', 'R', 'J', 'H', 'K', 'u', 'g', 'r', 'i', 'z']},\
        '2MASS': {'altname': [],'vref': u'vizier:II/246/out', 'select_columns': ['2MASS','RAJ2000','DEJ2000','Jmag','e_Jmag','Hmag','e_Hmag','Kmag','e_Kmag','MeasureJD']},\
        'DENIS': {'altname': [],'vref': u'vizier:B/denis/denis', 'select_columns': ['DENIS','RAJ2000','DEJ2000','Imag','e_Imag','Jmag','e_Jmag','Kmag','e_Kmag','Obs.JD']},\
        'SDSS': {'altname': ['SDSS16'],'vref': u'vizier:V/154/sdss16', 'select_columns': ['SDSS16','RAdeg','DEdeg','umag','e_umag','gmag','e_gmag','rmag','e_rmag','imag','e_imag','zmag','e_zmag','pmRA','e_pmRA','pmDE','e_pmDE','ObsDate','objID','SpObjID','spInst','spType','spCl','subCl','MJD']},\
        'SDSS12': {'altname': ['SDSS12'],'vref': u'vizier:V/147/sdss12', 'select_columns': ['SDSS12','RAdeg','DEdeg','umag','e_umag','gmag','e_gmag','rmag','e_rmag','imag','e_imag','zmag','e_zmag','pmRA','e_pmRA','pmDE','e_pmDE','ObsDate','objID','SpObjID','spType','spCl']},\
        'SDSS9': {'altname': [],'vref': u'vizier:V/139/sdss9', 'select_columns': ['SDSS9','RAdeg','DEdeg','umag','e_umag','gmag','e_gmag','rmag','e_rmag','imag','e_imag','zmag','e_zmag','pmRA','e_pmRA','pmDE','e_pmDE','ObsDate','objID','SpObjID','spType','spCl']},\
#        'CATWISE': {'altname': ['CAT'],'vref': u'vizier:II/365/catwise', 'select_columns': ['objID','RA_ICRS','DE_ICRS','Name','MJD','pmRA','e_pmRA','pmDE','e_pmDE','W1mproPM','e_W1mproPM','W2mproPM','e_W2mproPM']},\
        'ALLWISE': {'altname': [],'vref': u'vizier:II/328/allwise', 'select_columns': ['AllWISE','RAJ2000','DEJ2000','W1mag','e_W1mag','W2mag','e_W2mag','W3mag','e_W3mag','W4mag','e_W4mag','pmRA','e_pmRA','pmDE','e_pmDE','ID']},\
        'GAIA-DR1': {'altname': ['GAIADR1','GAIA1'],'vref': u'vizier:I/337/gaia', 'select_columns': ['source_id','ra','dec','ref_epoch','phot_g_mean_mag','phot_g_mean_flux','phot_g_mean_flux_error','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error']},\
        'GAIA-DR2': {'altname': ['GAIADR2','GAIA2'],'vref': u'vizier:I/345/gaia2', 'select_columns': ['source_id','ra','dec','phot_g_mean_mag','phot_g_mean_flux','phot_g_mean_flux_error','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error']},\
        'GAIA-EDR3': {'altname': ['GAIA-DR3','GAIAEDR3','GAIA3','GAIA'],'vref': u'vizier:I/350/gaiaedr3', 'select_columns': ['source_id','ra','dec','phot_g_mean_mag','phot_g_mean_flux','phot_g_mean_flux_error','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error']},\
        'PANSTARRS': {'altname': ['PAN-STARRS','PS1'], 'vref': u'vizier:II/349/ps1', 'select_columns': ['objID','RAJ2000','DEJ2000','Epoch','gmag','e_gmag','rmag','e_rmag','imag','e_imag','zmag','e_zmag','ymag','e_ymag']},
        'UKIDSS': {'altname': ['UKIDSS-LAS','UKIDSS-LAS9','UKIDSS-DR9','UKIDSS-LAS-DR9'], 'vref': u'vizier:II/319/las9', 'select_columns': ['JName','RAJ2000','DEJ2000','Epoch','yAperMag3','yAperMag3Err','j_1AperMag3','j_1AperMag3Err','hAperMag3','hAperMag3Err','kAperMag3','kAperMag3Err','mergedClass']},
# not yet integrated
#        'WISE': {'altname': ['WISE'],'vref': u'vizier:II/311/wise', 'select_columns': ['AllWISE','RAJ2000','DEJ2000','W1mag','e_W1mag','W2mag','e_W2mag','W3mag','e_W3mag','W4mag','e_W4mag','pmRA','e_pmRA','pmDE','e_pmDE','ID']},\
#        'UCAC': {'altname': ['UCAC'],'vref': u'vizier:II/322A/las9', 'select_columns': ['AllWISE','RAJ2000','DEJ2000','W1mag','e_W1mag','W2mag','e_W2mag','W3mag','e_W3mag','W4mag','e_W4mag','pmRA','e_pmRA','pmDE','e_pmDE','ID']},\
#        'MOVERS': {'altname': ['MOVERS'],'vref': u'vizier:J/AJ/151/41/movers', 'select_columns': ['AllWISE','RAJ2000','DEJ2000','W1mag','e_W1mag','W2mag','e_W2mag','W3mag','e_W3mag','W4mag','e_W4mag','pmRA','e_pmRA','pmDE','e_pmDE','ID']},\
#        'LATEMOVERS': {'altname': ['LATEMOVERS','LATE-MOVERS'],'vref': u'vizier:J/AJ/153/92/lmovers', 'select_columns': ['AllWISE','RAJ2000','DEJ2000','W1mag','e_W1mag','W2mag','e_W2mag','W3mag','e_W3mag','W4mag','e_W4mag','pmRA','e_pmRA','pmDE','e_pmDE','ID']},\
#        'WISE': {'vref': u'II/311', 'select_columns': 
#        'VISTA': {'vref': u'II/329', 'select_columns': 
#        'CFHT': {'vref': u'II/317', 'select_columns': 
#        'LEHPM': {'vref': u'J/A+A/421/763', 'select_columns': 
#        'SIPS': {'vref': u'J/A+A/435/363', 'select_columns': 
#        'UCAC': {'vref': u'I/340/ucac5', 'select_columns': 
#        'USNO': {'vref': u'I/284', 'select_columns': 
#        'LSPM': {'vref': u'I/298', 'select_columns': 
        }

# give a summary of the built-in catalogs
    if info==True:
        print('Currently available input catalogs:')
        for k in list(XMATCH_CATALOGS.keys()):
            line = '\t{}: '.format(k)
            if len(XMATCH_CATALOGS[k]['altname'])>0:
                line=line+'(or'
                for a in XMATCH_CATALOGS[k]['altname']: line=line+' {}'.format(a)
                line=line+') '
            print(line+'Vizier reference: {}'.format(str(XMATCH_CATALOGS[k]['vref'])))
            if 'vizier:' in str(XMATCH_CATALOGS[k]['vref']):
                catsp = str(XMATCH_CATALOGS[k]['vref']).split('/')
                ctref = catsp[0].replace('vizier:','') 
                for ct in catsp[1:-1]: ctref=ctref+'/'+ct
                print('\tVizier URL = https://cdsarc.unistra.fr/viz-bin/cat/{}\n'.format(ctref))
            else: print()

        return

# check db has DESIGNATION and fill in columns
#    print(db.columns,raCol in list(db.columns),decCol in list(db.columns))
    if desigCol not in list(db.columns) or raCol not in list(db.columns) or decCol not in list(db.columns):
        db = prepDB(db,raCol=raCol,decCol=decCol,desigCol=desigCol)
    if desigCol not in list(db.columns):
        raise ValueError('\nInput database must have at least the designation column {}; this one has {}'.format(desigCol,db.columns))

# add RA and DEC if needed
    # if raCol not in list(db.columns) or decCol not in list(db.columns):
    #     db['COORDINATES'] = [designationToCoordinate(d) for d in db[desigCol]]
    #     db[raCol] = [c.ra.degree for c in db['COORDINATES']]
    #     db[decCol] = [c.dec.degree for c in db['COORDINATES']]
    basecols = [desigCol,raCol,decCol]
    if not isUnit(radius): radius = radius*u.arcsec
        
# define catalog
    if len(args) > 0: catalog = args[0]
    cat = checkDict(catalog,XMATCH_CATALOGS)
    if cat == False: 
        cat = catalog.upper()
        vref = 'vizier:'+catalog
    else: 
        vref = XMATCH_CATALOGS[cat]['vref']
#    if catalog.upper() in list(XMATCH_CATALOGS.keys()):
#        cat = catalog.upper()
#        vref = XMATCH_CATALOGS[cat]['vref']
        if use_select_columns == True and len(XMATCH_CATALOGS[cat]['select_columns']) > 0: 
            select_columns = XMATCH_CATALOGS[cat]['select_columns']
#        else: select_columns = []

# check that catalog is there
    if XMatch.is_table_available(vref) == False:
        print('\n{} is not one of the catalogs in astroquery.xmatch; try using queryVizer()'.format(catalog))
        return db
    if prefix == None: prefix = cat

# use XMatch
    t = Table()
    t = t.from_pandas(db[basecols])
    t_match = XMatch.query(t,vref,radius,colRA1=raCol,colDec1=decCol,columns=["**", "+_r"])
    db_match = t_match.to_pandas()
    if debug==True: 
        print('Found {} matches'.format(len(db_match)))
        print(db_match.iloc[0])

# reject repeats if desired
    if drop_repeats == True:
        db_match.drop_duplicates(subset=desigCol,keep='first',inplace=True)
        db_match.reset_index(drop=True,inplace=True)

# constrain columns and rename
    if len(select_columns)>0:
        if len(select_columns) == 0: 
            newcols = list(db_match.columns)
        else:
            newcols = copy.deepcopy(basecols)
            newcols.append('angDist')
            newcols.extend(select_columns)
# check that all columns are present
        ncdup = copy.deepcopy(newcols)
        for s in ncdup:
            if s not in list(db_match.columns): 
                print('Warning: could not find column named {}'.format(s))
                newcols.remove(s)
        if len(newcols) > 0: db_match = db_match[newcols]

# rename columns
    if prefix != None:
        rename = {}
        for c in list(db_match.columns): rename[c] = prefix+'_'+c
        for c in list(basecols): rename[c] = c
        db_match = db_match.rename(index=str,columns=rename)

# merge and drop redundant columns
    db_merge = pandas.merge(db,db_match,how='left',on=desigCol,suffixes=('','_DROP'))
    for c in list(db_merge.columns):
        if '_DROP' in c: del db_merge[c]

    if debug==True: 
        print(db_merge.iloc[0])
            

# save out
    if file != '':
        if file.split('.')[-1] == 'csv' or file.split('.')[-1] == 'txt':   
            db_merge.to_csv(file,index=False)
        elif file.split('.')[-1] == 'xls' or file.split('.')[-1] == 'xlsx': 
            db_merge.to_excel(file,index=False)
        else:
            print('\nWarning: did not know how to save to {}; not saving'.format(file))
                
    return db_merge




#####################################################
###########   ADDING SPECTRA TO LIBRARY   ###########
#####################################################



def importSpectra(*args,**kwargs):
    '''
    Purpose
        imports a set of spectra into the SPLAT library; requires manager access.

    :Note:
        **this program is in beta testing; bugs/errors are likely**

    :Optional parameters:
        :param data_folder = "./": Full path to folder containing data; by default this is the current directory
        :param review_folder = "./review/": Full path to folder in which review materials will be kept; by default a new folder ``review`` will be created inside the data_folder
        :param spreadsheet = "": Filename for a spreadsheet (ascii, tab- or comma-delimited) listing the input spectra, one per row. At least one column must be named ``filename`` or ``file`` that contains the name of the data file; the following columns are also recommended:

            * ``designation``: source desigation; e.g., ``J15420830-2621138`` (strongly recommended)
            * ``ra`` and ``dec``: Right Ascension and declination in decimal format (only needed if no designation column provided)
            * ``name``: source name, designation will be used if not provided
            * ``type``, ``opt_type``, ``nir_type``: spectral type of source (string); ``type`` will default to ``lit_type``
            * ``date`` or ``observation_date``: date of observation in format YYYYMMDD
            * ``slit``: slit width used (for computing resolution)
            * ``airmass``: airmass of observation
            * ``observer``: last name of primary observer
            * ``bibcode``: bibcode of data reference

    :Output:
        - Source DB update file: spreadsheet containing update to source_data.txt, saved in review folder as source_data.txt
        - Spectral DB update file: spreadsheet containing update to spectral_data.txt, saved locally as UPDATE_spectral_data.txt
        - Photometry DB update file: spreadsheet containing update to photometry_data.txt, saved locally as UPDATE_photometry_data.txt

    '''
# check user access
    if checkAccess() == False:
        print('\nSpectra may only be imported into library by designated manager or while online; please email {}'.format(SPLAT_EMAIL))
        return

# check online
#    if spl.checkOnline() == False:
#        print('\nWarning! You are not currently online so you will not be able to retrieve SIMBAD and Vizier data\n')

# set up optional inputs
    simbad_radius = kwargs.get('simbad_radius',60.*u.arcsec)
    if not isUnit(simbad_radius): simbad_radius=simbad_radius*u.arcsec

    vizier_radius = kwargs.get('vizier_radius',30.*u.arcsec)
    if not isUnit(vizier_radius): vizier_radius=vizier_radius*u.arcsec

    data_folder = kwargs.get('data_folder','./')
    data_folder = kwargs.get('dfolder',data_folder)
    data_folder = kwargs.get('folder',data_folder)
    if data_folder[-1] != '/':
        data_folder+='/'
    review_folder = kwargs.get('review_folder','{}/review/'.format(data_folder))
    review_folder = kwargs.get('rfolder',review_folder)
    if review_folder[-1] != '/':
        review_folder+='/'
    spreadsheet = kwargs.get('spreadsheet','')
    spreadsheet = kwargs.get('sheet',spreadsheet)
    spreadsheet = kwargs.get('entry',spreadsheet)
    instrument = kwargs.get('instrument','UNKNOWN')
    verbose = kwargs.get('verbose',True)

# make sure relevant files and folders are in place
    if not os.path.exists(review_folder):
        try:
            os.makedirs(review_folder)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
#        raise NameError('\nCannot find review folder {}'.format(review_folder))
    if not os.path.exists(data_folder):
        raise NameError('\nCannot find data folder {}'.format(data_folder))
    if not os.path.exists('{}/published'.format(review_folder)):
        try:
            os.makedirs('{}/published'.format(review_folder))
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    if not os.path.exists('{}/unpublished'.format(review_folder)):
        try:
            os.makedirs('{}/unpublished'.format(review_folder))
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

# if spreadsheet is given, use this to generate list of files
    if spreadsheet != '':
        try:
            t_input = fetchDatabase(spreadsheet)        
        except:
            try:
                t_input = fetchDatabase(data_folder+spreadsheet)        
            except:
                raise NameError('\nCould not find spreadsheet {} in local or data directories\n'.format(spreadsheet))
        tkeys = list(t_input.keys())
        if 'FILENAME' in tkeys:
            files = t_input['FILENAME']
        elif 'FILE' in tkeys:
            files = t_input['FILE']
        elif 'FILES' in tkeys:
            files = t_input['FILES']
        else:
            raise NameError('\nSpreadsheet {} does not have a column named filename; aborting\n'.format(spreadsheet))
        if data_folder not in files[0]:
            files = [data_folder+f for f in files]

# otherwise search for *.fits and *.txt files in data folder
    else:
        files = glob.glob(os.path.normpath(data_folder+'*.fits'))+glob.glob(os.path.normpath(data_folder+'*.txt'))
        if len(files) == 0:
            raise NameError('\nNo spectral files in {}\n'.format(data_folder))

# what instrument is this?
    s = splat.Spectrum(filename=files[0])
    if 'INSTRUME' in list(s.header.keys()):
        instrument = s.header['INSTRUME'].replace(' ','').upper()
    if 'INSTR' in list(s.header.keys()):
        instrument = s.header['INSTR'].replace(' ','').upper()
        if 'MODENAME' in list(s.header.keys()):
            instrument+=' {}'.format(s.header['MODENAME'].replace(' ','').upper())

    if instrument.upper().replace(' ','_') in list(INSTRUMENTS.keys()):
        instrument_info = INSTRUMENTS[instrument.upper().replace(' ','_')]
    else:
        instrument_info = {'instrument_name': instrument, 'resolution': 0.*u.arcsec, 'slitwidth': 0.}

# prep tables containing information
    t_spec = Table()
    for c in list(DB_SPECTRA.keys()):
        t_spec[c] = Column([' '*200 for f in files],dtype='str')
    t_src = Table()
    for c in list(DB_SOURCES.keys()):
        t_src[c] = Column([' '*200 for f in files],dtype='str')
    source_id0 = numpy.max(DB_SOURCES['SOURCE_KEY'])
    spectrum_id0 = numpy.max(DB_SPECTRA['DATA_KEY'])

# read in files into Spectrum objects
    if verbose: print('\nReading in {} files from {}'.format(len(files),data_folder))
#    splist = []
    t_spec['DATA_FILE'] = Column(files,dtype='str')
    t_spec['SPECTRUM'] = [splat.Spectrum(filename=f) for f in files]
    t_spec['INSTRUMENT'] = [instrument_info['instrument_name'] for f in files]
#    for f in files:
#        splist.append()

# populate spec array
    if verbose: print('\nGenerating initial input tables')
    t_spec['SOURCE_KEY'] = Column(numpy.arange(len(files))+source_id0+1,dtype='int')
    t_spec['DATA_KEY'] = Column(numpy.arange(len(files))+spectrum_id0+1,dtype='int')
#    t_spec['SPECTRUM'] = [sp for sp in splist]
    t_spec['QUALITY_FLAG'] = Column(['OK' for f in t_spec['DATA_FILE']],dtype='str')
    t_spec['PUBLISHED'] = Column(['N' for f in t_spec['DATA_FILE']],dtype='str')
#  measurements
    t_spec['MEDIAN_SNR'] = Column([sp.computeSN() for sp in t_spec['SPECTRUM']],dtype='float')
    t_spec['SPEX_TYPE'] = Column([classifyByStandard(sp,string=True,method=kwargs.get('method','kirkpatrick'),mask_telluric=True)[0] for sp in t_spec['SPECTRUM']],dtype='str')
    t_spec['SPEX_GRAVITY_CLASSIFICATION'] = Column([classifyGravity(sp,string=True) for sp in t_spec['SPECTRUM']],dtype='str')
# populate spectral data table from fits file header
    for i,sp in enumerate(t_spec['SPECTRUM']):
        if 'DATE_OBS' in list(sp.header.keys()):
            t_spec['OBSERVATION_DATE'][i] = sp.header['DATE_OBS'].replace('-','')
            t_spec['JULIAN_DATE'][i] = Time(sp.header['DATE_OBS']).mjd
        if 'DATE' in list(sp.header.keys()):
            t_spec['OBSERVATION_DATE'][i] = sp.header['DATE'].replace('-','')
            if verbose: print(i,t_spec['OBSERVATION_DATE'][i],properDate(t_spec['OBSERVATION_DATE'][i],output='YYYYMMDD'))
            t_spec['JULIAN_DATE'][i] = Time(sp.header['DATE']).mjd
        if 'TIME_OBS' in list(sp.header.keys()):
            t_spec['OBSERVATION_TIME'][i] = sp.header['TIME_OBS'].replace(':',' ')
        if 'MJD_OBS' in list(sp.header.keys()):
            t_spec['JULIAN_DATE'][i] = sp.header['MJD_OBS']
        if 'OBSERVER' in list(sp.header.keys()):
            t_spec['OBSERVER'][i] = sp.header['OBSERVER']
        if 'RESOLUTION' in list(sp.header.keys()):
            t_spec['RESOLUTION'][i] = sp.header['RESOLUTION']
        elif 'RES' in list(sp.header.keys()):
            t_spec['RESOLUTION'][i] = sp.header['RES']
        elif 'SLITW' in list(sp.header.keys()):
            t_spec['RESOLUTION'][i] = instrument_info['resolution']*(instrument_info['slitwidth'].value)/sp.header['SLITW']
        elif 'SLTW_ARC' in list(sp.header.keys()):
            t_spec['RESOLUTION'][i] = instrument_info['resolution']*(instrument_info['slitwidth'].value)/sp.header['SLTW_ARC']
        if 'AIRMASS' in list(sp.header.keys()):
            t_spec['AIRMASS'][i] = sp.header['AIRMASS']
        if 'VERSION' in list(sp.header.keys()):
            v = sp.header['VERSION']
            t_spec['REDUCTION_SPEXTOOL_VERSION'][i] = 'v{}'.format(v.split('v')[-1])
# populate spectral data table from spreadsheet 
    if spreadsheet != '':
#        if 'FILENAME' in tkeys:
#            t_spec['DATA_FILE'] = t_input['FILENAME']
        if 'DATE' in tkeys:
            t_spec['OBSERVATION_DATE'] = [properDate(str(a),output='YYYYMMDD') for a in t_input['DATE']]
#            for a in t_input['DATE']:
#                print(a,spl.properDate(str(a)),Time(spl.properDate(str(a),output='YYYY-MM-DD')),Time(spl.properDate(str(a),output='YYYY-MM-DD')).mjd)
            t_spec['JULIAN_DATE'] = [Time(properDate(str(a),output='YYYY-MM-DD')).mjd for a in t_input['DATE']]
        if 'RESOLUTION' in tkeys:
            t_spec['RESOLUTION'] = [r for r in t_input['RESOLUTION']]
# CHANGE THIS TO BE INSTRUMENT SPECIFIC
        if 'SLIT' in tkeys:
            t_spec['RESOLUTION'] = [t_spec['RESOLUTION']*(instrument_info['slitwidth'].value)/float(s) for s in t_input['SLIT']]
        if 'AIRMASS' in tkeys:
            t_spec['AIRMASS'] = t_input['AIRMASS']
        if 'OBSERVER' in tkeys:
            t_spec['OBSERVER'] = t_input['OBSERVER']
        if 'BIBCODE' in tkeys:
            t_spec['BIBCODE'] = t_input['BIBCODE']
            for i,ref in enumerate(t_spec['BIBCODE']):
                if ref != '':
                    t_spec['PUBLISHED'][i] = 'Y'

#    for c in splist[0].header.keys():
#        if c != 'HISTORY':
#            print('{} {}'.format(c,splist[0].header[c]))

    t_src['SOURCE_KEY'] = t_spec['SOURCE_KEY']
    t_src['GRAVITY_CLASS_NIR'] = t_spec['SPEX_GRAVITY_CLASSIFICATION']
    t_src['GRAVITY_CLASS_NIR_REF'] = Column(['SPL' for sp in t_spec['SPECTRUM']],dtype='str')
    t_spec['COMPARISON_SPECTRUM'] = [STDS_DWARF_SPEX[spt] for spt in t_spec['SPEX_TYPE']]
    t_spec['COMPARISON_TEXT'] = [' '*200 for spt in t_spec['SPEX_TYPE']]
    for i,spt in enumerate(t_spec['SPEX_TYPE']):
        t_spec['COMPARISON_TEXT'][i] = '{} standard'.format(spt)

# determine coordinates as best as possible
    for i,sp in enumerate(t_spec['SPECTRUM']):
#        if i == 0:
#            for k in list(sp.header.keys()):
#                print(k,sp.header[k])
        if 'TCS_RA' in list(sp.header.keys()) and 'TCS_DEC' in list(sp.header.keys()):
            sp.header['RA'] = sp.header['TCS_RA']
            sp.header['DEC'] = sp.header['TCS_DEC']
            sp.header['RA'] = sp.header['RA'].replace('+','')
        if t_src['DESIGNATION'][i].strip() == '' and 'RA' in list(sp.header.keys()) and 'DEC' in list(sp.header.keys()):
            if sp.header['RA'] != '' and sp.header['DEC'] != '':
                t_src['DESIGNATION'][i] = 'J{}+{}'.format(sp.header['RA'].replace('+',''),sp.header['DEC']).replace(':','').replace('.','').replace('+-','-').replace('++','+').replace('J+','J').replace(' ','')
#            print('DETERMINED DESIGNATION {} FROM RA/DEC'.format(t_src['DESIGNATION'][i]))
        if t_src['RA'][i].strip() == '' and t_src['DESIGNATION'][i].strip() != '':
            coord = properCoordinates(t_src['DESIGNATION'][i])
            t_src['RA'][i] = coord.ra.value
            t_src['DEC'][i] = coord.dec.value
#            print('DETERMINED RA/DEC FROM DESIGNATION {}'.format(t_src['DESIGNATION'][i]))
#    print(t_src['DESIGNATION'],t_src['RA'],t_src['DEC'])
# populate source data table from spreadsheet
    if spreadsheet != '':
        if 'DESIGNATION' in tkeys:
            t_src['DESIGNATION'] = t_input['DESIGNATION']
            t_src['NAME'] = t_src['DESIGNATION']
# may want to check how we overrule fits file headers
            coord = [properCoordinates(s) for s in t_src['DESIGNATION']]
            t_src['RA'] = [c.ra.value for c in coord]
            t_src['DEC'] = [c.dec.value for c in coord]
        if 'NAME' in tkeys:
            t_src['NAME'] = t_input['NAME']
        if 'RA' in tkeys and 'DEC' in tkeys:
            if isNumber(t_input['RA'][0]):
                t_src['RA'] = t_input['RA']
                t_src['DEC'] = t_input['DEC']
        if 'TYPE' in tkeys:
            t_src['LIT_TYPE'] = t_input['TYPE']
        if 'OPT_TYPE' in tkeys:
            t_src['OPT_TYPE'] = t_input['OPT_TYPE']
        if 'NIR_TYPE' in tkeys:
            t_src['NIR_TYPE'] = t_input['NIR_TYPE']
        if 'J' in tkeys:
            t_src['J_2MASS'] = t_input['J']
        if 'J_E' in tkeys:
            t_src['J_2MASS_E'] = t_input['J_E']
        if 'H' in tkeys:
            t_src['H_2MASS'] = t_input['H']
        if 'H_E' in tkeys:
            t_src['H_2MASS_E'] = t_input['H_E']
        if 'K' in tkeys:
            t_src['KS_2MASS'] = t_input['K']
        if 'KS' in tkeys:
            t_src['KS_2MASS'] = t_input['KS']
        if 'K_E' in tkeys:
            t_src['KS_2MASS_E'] = t_input['K_E']
        if 'KS_E' in tkeys:
            t_src['KS_2MASS_E'] = t_input['KS_E']

#    for c in DB_SOURCES.keys():
#        if c not in t_src.keys():
#            t_src[c] = Column([' '*50 for sp in splist],dtype='str')        # force string

# transfer spectral types
    for i,t in enumerate(t_src['NIR_TYPE']):
        if t.replace(' ','') == '':
            t_src['NIR_TYPE'][i] = t_spec['SPEX_TYPE'][i]
            t_src['NIR_TYPE_REF'][i] = 'SPL'
        if t_src['LIT_TYPE'][i].replace(' ','') == '':
            t_src['LIT_TYPE'][i] = t_spec['SPEX_TYPE'][i]
            t_src['LIT_TYPE_REF'][i] = 'SPL'


# now do a SIMBAD search for sources based on coordinates
    if kwargs.get('nosimbad',False) == False:
        if verbose:
            print('\nSIMBAD search')
        _querySimbad2(t_src,simbad_radius=simbad_radius)


# fill in missing 2MASS photometry with Vizier query
    if kwargs.get('novizier',False) == False:
        if verbose:
            print('\n2MASS photometry from Vizier')

        if not checkOnline():
            if verbose:
                print('\nCould not perform Vizier search, you are not online')
        else:
            for i,jmag in enumerate(t_src['J_2MASS']):
                if float('{}0'.format(jmag.replace('--',''))) == 0.0:
                    t_vizier = getPhotometry(properCoordinates(t_src['DESIGNATION'][i]),radius=vizier_radius,catalog='2MASS')

        # multiple sources; choose the closest
                    if len(t_vizier) > 0:
                        t_vizier.sort_values('_r')
        #                print(len(t_vizier),t_vizier.keys())
        #                while len(t_vizier)>1:
        #                    t_vizier.remove_row(1) 
                        if verbose:
                            print('\n{}'.format(t_src['DESIGNATION'][i]))
                            print(t_vizier)
                        t_src['DESIGNATION'][i] = 'J{}'.format(t_vizier['_2MASS'][0])
                        t_src['J_2MASS'][i] = str(t_vizier['Jmag'][0]).replace('--','')
                        t_src['J_2MASS_E'][i] = str(t_vizier['e_Jmag'][0]).replace('--','')
                        t_src['H_2MASS'][i] = str(t_vizier['Hmag'][0]).replace('--','')
                        t_src['H_2MASS_E'][i] = str(t_vizier['e_Hmag'][0]).replace('--','')
                        t_src['KS_2MASS'][i] = str(t_vizier['Kmag'][0]).replace('--','')
                        t_src['KS_2MASS_E'][i] = str(t_vizier['e_Kmag'][0]).replace('--','')

    # add in distance if spectral type and magnitude are known
    for i,spt in enumerate(t_src['LIT_TYPE']):
        if spt.replace(' ','') != '' and float('{}0'.format(str(t_src['J_2MASS'][i]).replace('--',''))) != 0.0:
    #            print(spt,t_src['J_2MASS'][i],t_src['J_2MASS_E'][i])
            dist = estimateDistance(spt=spt,filter='2MASS J',mag=float(t_src['J_2MASS'][i]))
            if not numpy.isnan(dist[0]):
                t_src['DISTANCE_PHOT'][i] = dist[0]
                t_src['DISTANCE_PHOT_E'][i] = dist[1]
                t_src['DISTANCE'][i] = dist[0]
                t_src['DISTANCE_E'][i] = dist[1]
        if float('{}0'.format(str(t_src['PARALLAX'][i]).replace('--',''))) != 0.0 and float('{}0'.format(str(t_src['PARALLAX_E'][i]).replace('--',''))) != 0.0 :
            t_src['DISTANCE'][i] = 1000./float(t_src['PARALLAX'][i])
            t_src['DISTANCE_E'][i] = float(t_src['DISTANCE'][i])*float(t_src['PARALLAX_E'][i])/float(t_src['PARALLAX'][i])
    # compute vtan
        if float('{}0'.format(str(t_src['MU'][i]).replace('--',''))) != 0.0 and float('{}0'.format(str(t_src['DISTANCE'][i]).replace('--',''))) != 0.0:
            t_src['VTAN'][i] = 4.74*float(t_src['DISTANCE'][i])*float(t_src['MU'][i])/1000.

    # clear up zeros
        if float('{}0'.format(str(t_src['J_2MASS'][i]).replace('--',''))) == 0.0:
            t_src['J_2MASS'][i] = ''
            t_src['J_2MASS_E'][i] = ''
        if float('{}0'.format(str(t_src['H_2MASS'][i]).replace('--',''))) == 0.0:
            t_src['H_2MASS'][i] = ''
            t_src['H_2MASS_E'][i] = ''
        if float('{}0'.format(str(t_src['KS_2MASS'][i]).replace('--',''))) == 0.0:
            t_src['KS_2MASS'][i] = ''
            t_src['KS_2MASS_E'][i] = ''
        if float('{}0'.format(str(t_src['PARALLAX'][i]).replace('--',''))) == 0.0:
            t_src['PARALLAX'][i] = ''
            t_src['PARALLAX_E'][i] = ''
        if float('{}0'.format(str(t_src['MU'][i]).replace('--',''))) == 0.0:
            t_src['MU'][i] = ''
            t_src['MU_E'][i] = ''
            t_src['MU_RA'][i] = ''
            t_src['MU_DEC'][i] = ''
        if float('{}0'.format(str(t_src['RV'][i]).replace('--',''))) == 0.0:
            t_src['RV'][i] = ''
            t_src['RV_E'][i] = ''
        if float('{}0'.format(str(t_src['VSINI'][i]).replace('--',''))) == 0.0:
            t_src['VSINI'][i] = ''
            t_src['VSINI_E'][i] = ''
        if float('{}0'.format(str(t_src['SIMBAD_SEP'][i]).replace('--',''))) == 0.0:
            t_src['SIMBAD_SEP'][i] = ''
        if t_src['GRAVITY_CLASS_NIR'][i] == '':
            t_src['GRAVITY_CLASS_NIR_REF'][i] = ''

    # compute J-K excess and color extremity
        if spt.replace(' ','') != '' and float('{}0'.format(str(t_src['J_2MASS'][i]).replace('--',''))) != 0.0 and float('{}0'.format(str(t_src['KS_2MASS'][i]).replace('--',''))) != 0.0:
            t_src['JK_EXCESS'][i] = float(t_src['J_2MASS'][i])-float(t_src['KS_2MASS'][i])-typeToColor(spt,'J-K')[0]
            if t_src['JK_EXCESS'][i] == numpy.nan or t_src['JK_EXCESS'][i] == '' or t_src['JK_EXCESS'][i] == 'nan':
                t_src['JK_EXCESS'][i] = ''
            elif float(t_src['JK_EXCESS'][i]) > 0.3:
                t_src['COLOR_EXTREMITY'][i] == 'RED'
            elif float(t_src['JK_EXCESS'][i]) < -0.3:
                t_src['COLOR_EXTREMITY'][i] == 'BLUE'
            else:
                pass


# check for previous entries
    t_src['SHORTNAME'] = [designationToShortName(d) for d in t_src['DESIGNATION']]
    if 'SHORTNAME' not in list(DB_SOURCES.keys()):
        DB_SOURCES['SHORTNAME'] = [designationToShortName(d) for d in DB_SOURCES['DESIGNATION']]
    for i,des in enumerate(t_src['DESIGNATION']):

# check if shortnames line up
        if t_src['SHORTNAME'][i] in DB_SOURCES['SHORTNAME']:
            for c in list(t_src.keys()):
                t_src[c][i] = DB_SOURCES[c][numpy.where(DB_SOURCES['SHORTNAME'] == t_src['SHORTNAME'][i])][0]
            t_spec['SOURCE_KEY'][i] = t_src['SOURCE_KEY'][i]

# check if SIMBAD names line up
        elif t_src['SIMBAD_NAME'][i] != '' and t_src['SIMBAD_NAME'][i] in DB_SOURCES['SIMBAD_NAME']:
            for c in t_src.keys():
                if t_src[c][i] == '':
                    t_src[c][i] = DB_SOURCES[c][numpy.where(DB_SOURCES['SIMBAD_NAME'] == t_src['SIMBAD_NAME'][i])][0]
            t_spec['SOURCE_KEY'][i] = t_src['SOURCE_KEY'][i]

        else:
            pass

# check to see if prior spectrum was taken on the same date (possible redundancy)
        matchlib = searchLibrary(idkey=t_src['SOURCE_KEY'][i],date=t_spec['OBSERVATION_DATE'][i])
# previous observation on this date found - retain in case this is a better spectrum
        if len(matchlib) > 0.:
            mkey = matchlib['DATA_KEY'][0]
            if verbose:
                print('Previous spectrum found in library for data key {}'.format(mkey))
            t_spec['COMPARISON_SPECTRUM'][i] = splat.Spectrum(int(mkey))
            t_spec['COMPARISON_TEXT'][i] = 'repeat spectrum: {}'.format(mkey)
# no previous observation on this date - retain the spectrum with the highest S/N
        else:
            matchlib = searchLibrary(idkey=t_src['SOURCE_KEY'][i])
            if len(matchlib) > 0:
                matchlib.sort('MEDIAN_SNR')
                matchlib.reverse()
                t_spec['COMPARISON_SPECTRUM'][i] = splat.Spectrum(int(matchlib['DATA_KEY'][0]))
                t_spec['COMPARISON_TEXT'][i] = 'alternate spectrum: {} taken on {}'.format(matchlib['DATA_KEY'][0],matchlib['OBSERVATION_DATE'][0])
#                print(matchlib['DATA_KEY'][0])
#                print(t_spec['COMPARISON_TEXT'][i])


# generate check plots
    legend = []
    for i,sp in enumerate(t_spec['SPECTRUM']):
        legend.extend(['Data Key: {} Source Key: {}\n{}'.format(t_spec['DATA_KEY'][i],t_spec['SOURCE_KEY'][i],t_spec['SPECTRUM'][i].name),'{} {}'.format(t_spec['COMPARISON_SPECTRUM'][i].name,t_spec['COMPARISON_TEXT'][i])])
    for s in t_spec['COMPARISON_SPECTRUM']: print(s)
    splot.plotBatch([s for s in t_spec['SPECTRUM']],comparisons=[s for s in t_spec['COMPARISON_SPECTRUM']],normalize=True,output=review_folder+'/review_plots.pdf',legend=legend,noise=True,telluric=True)


# output database updates
    if 'SHORTNAME' in t_src.keys():
        t_src.remove_column('SHORTNAME')
    if 'SELECT' in t_src.keys():
        t_src.remove_column('SELECT')
    if 'SELECT' in t_spec.keys():
        t_spec.remove_column('SELECT')   
    if 'SOURCE_SELECT' in t_spec.keys():
        t_spec.remove_column('SOURCE_SELECT')
    if 'SPECTRUM' in t_spec.keys():
        t_spec.remove_column('SPECTRUM')
    if 'COMPARISON_SPECTRUM' in t_spec.keys():
        t_spec.remove_column('COMPARISON_SPECTRUM')
    if 'COMPARISON_TEXT' in t_spec.keys():
        t_spec.remove_column('COMPARISON_TEXT')
#    for i in numpy.arange(len(t_spec['NOTE'])):
#        t_spec['NOTE'][i] = compdict[str(t_spec['DATA_KEY'][i])]['comparison_type']
    t_src.write(review_folder+'/source_update.csv',format='ascii.csv')
    t_spec.write(review_folder+'/spectrum_update.csv',format='ascii.csv')

# open up windows to review spreadsheets
# NOTE: WOULD LIKE TO MAKE THIS AUTOMATICALLY OPEN FILE
#    app = QtGui.QApplication(sys.argv)
#    window = Window(10, 5)
#    window.resize(640, 480)
#    window.show()
#    app.exec_()

    print('\nSpectral plots and update speadsheets now available in {}'.format(review_folder))
    response = input('Please review and edit, and press any key when you are finished...\n')


# NEXT STEP - MOVE FILES TO APPROPRIATE PLACES, UPDATE MAIN DATABASES
# source db
    t_src = fetchDatabase(review_folder+'/source_update.csv',csv=True)
#    if 'SIMBAD_SEP' in t_src.keys():
#        t_src.remove_column('SIMBAD_SEP')

#    for col in t_src.colnames:
#        tmp = t_src[col].astype(DB_SOURCES[col].dtype)
#        t_src.replace_column(col,tmp)
#    t_merge = vstack([DB_SOURCES,t_src])
#    t_merge.sort('SOURCE_KEY')
#    if 'SHORTNAME' in t_merge.keys():
#        t_merge.remove_column('SHORTNAME')
#    if 'SELECT' in t_merge.keys():
#        t_merge.remove_column('SELECT')
#    t_merge.write(review_folder+DB_SOURCES_FILE,format='ascii.tab')

# spectrum db
    t_spec = fetchDatabase(review_folder+'/spectrum_update.csv',csv=True)

# move files
    for i,file in enumerate(t_spec['DATA_FILE']):
        t_spec['DATA_FILE'][i] = '{}_{}.fits'.format(t_spec['DATA_KEY'][i],t_spec['SOURCE_KEY'][i])
#        print(file[-4:],t_spec['DATA_FILE'][i])
        if file[-4:] == 'fits':
            if t_spec['PUBLISHED'][i] == 'Y':
                copyfile(file,'{}/published/{}'.format(review_folder,t_spec['DATA_FILE'][i]))
#                if verbose:
#                    print('Moved {} to {}/published/'.format(t_spec['DATA_FILE'][i],review_folder))
            else:
                copyfile(file,'{}/unpublished/{}'.format(review_folder,t_spec['DATA_FILE'][i]))
#                if verbose:
#                    print('Moved {} to {}/unpublished/'.format(t_spec['DATA_FILE'][i],review_folder))
        else:
#            print(data_folder+file)
            sp = splat.Spectrum(file=file)
            if t_spec['PUBLISHED'][i] == 'Y':
                sp.export('{}/published/{}'.format(review_folder,t_spec['DATA_FILE'][i]))
#                if verbose:
#                    print('Moved {} to {}/published/'.format(t_spec['DATA_FILE'][i],review_folder))
            else:
                sp.export('{}/unpublished/{}'.format(review_folder,t_spec['DATA_FILE'][i]))
#                if verbose:
#                    print('Moved {} to {}/unpublished/'.format(t_spec['DATA_FILE'][i],review_folder))

# save off updated spectrum update file
    t_spec.write(review_folder+'/spectrum_update.csv',format='ascii.csv')

# merge and export - THIS WASN'T WORKING
#    for col in t_spec.colnames:
#        print(col,DB_SPECTRA[col].dtype)
#        tmp = t_spec[col].astype(DB_SPECTRA[col].dtype)
#        t_spec.replace_column(col,tmp)
#    t_merge = vstack([DB_SPECTRA,t_spec])
#    t_merge.sort('DATA_KEY')
#    if 'SHORTNAME' in t_merge.keys():
#        t_merge.remove_column('SHORTNAME')
#    if 'SELECT' in t_merge.keys():
#        t_merge.remove_column('SELECT')
#    if 'SOURCE_SELECT' in t_merge.keys():
#        t_merge.remove_column('SOURCE_SELECT')
#    if 'DATEN' in t_merge.keys():
#        t_merge.remove_column('DATEN')
#    t_merge.write(review_folder+DB_SPECTRA_FILE,format='ascii.tab')

    if verbose:
        print('\nDatabases updated; be sure to add these to primary databases in {}'.format(DB_FOLDER))
        print('and to move spectral files from {}/published and {}/unpublished/ to {}\n'.format(review_folder,review_folder,DATA_FOLDER))

    return


