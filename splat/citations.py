# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
.. note::
         These are the database functions for SPLAT 
"""

# imports: internal
import os
import re
import requests

# imports: external
import numpy

# splat functions and codes
from splat.initialize import *
from splat.utilities import *
import splat


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
    if 'journal' in list(bib_dict.keys()):
        if bib_dict['journal'] == '\\apj':
            bib_dict['journal'] = 'ApJ'
        elif bib_dict['journal'] == '\\apjl':
            bib_dict['journal'] = 'ApJ Letters'
        elif bib_dict['journal'] == '\\apjs':
            bib_dict['journal'] = 'ApJS'
        elif bib_dict['journal'] == '\\aj':
            bib_dict['journal'] = 'AJ'
        elif bib_dict['journal'] == '\\actaa':
            bib_dict['journal'] = 'Acta Astronomica'
        elif bib_dict['journal'] == '\\araa':
            bib_dict['journal'] = 'AR&A'
        elif bib_dict['journal'] == '\\aap':
            bib_dict['journal'] = 'A&A'
        elif bib_dict['journal'] == '\\icarus':
            bib_dict['journal'] = 'Icarus'
        elif bib_dict['journal'] == '\\mnras':
            bib_dict['journal'] = 'MNRAS'
        elif bib_dict['journal'] == '\\nat':
            bib_dict['journal'] = 'Nature'
        elif bib_dict['journal'] == '\\pasp':
            bib_dict['journal'] = 'PASP'
        elif bib_dict['journal'] == '\\solphys':
            bib_dict['journal'] = 'Solar Physics'
        elif bib_dict['journal'] == '\\pnas':
            bib_dict['journal'] = 'PNAS'
        else:
            pass
    else: 
        bib_dict['journal'] = 'UNKNOWN'
        
    return bib_dict


def veryShortRef(bib_dict,**kwargs):
    '''
    :Purpose:
        Takes a bibtex dictionary and returns a short (in-line) version of the citation

    :Required parameters:
        :param bib_tex: Dictionary output from bibTexParser, else a bibcode that is fed into bibTexParser

    :Optional parameters:
        None

    :Output:
        A string of the format ``Burgasser et al. (2006)``

    '''
    if type(bib_dict) is not dict:
        if type(bib_dict) is numpy.str:
            bib_dict = str(bib_dict)
        if type(bib_dict) is str:
            bib_dict = getBibTex(bib_dict,**kwargs)
            if isinstance(bib_dict,dict) == False: return ''
        else:
            if kwargs.get('verbose',False): print('Input to shortRef is neither a bibcode nor a bibTex dictionary')
            return ''

    authors = bib_dict['author'].split(' and ')
    a = authors[0].replace('~',' ').split(' ')
    a = a[0].replace(',','')
    if len(authors) == 1:
        output = a
    else:
        output = '{} et al.'.format(a)

# fill in missing data
    if 'year' not in bib_dict.keys():
        bib_dict['year'] = ''

    return output+' ({})'.format(bib_dict['year'])

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
        if type(bib_dict) is numpy.str:
            bib_dict = str(bib_dict)
        if type(bib_dict) is str:
            bib_dict = getBibTex(bib_dict,**kwargs)
            if isinstance(bib_dict,dict) == False: return ''
        else:
            if kwargs.get('verbose',False): print('Input to shortRef is neither a bibcode nor a bibTex dictionary')
            return ''

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

# fill in missing data
    if 'year' not in bib_dict.keys():
        bib_dict['year'] = ''
    if 'journal' not in bib_dict.keys():
        bib_dict['journal'] = ''
    if 'volume' not in bib_dict.keys():
        bib_dict['volume'] = ''
    if 'pages' not in bib_dict.keys():
        bib_dict['pages'] = ''

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
        if type(bib_dict) is numpy.str:
            bib_dict = str(bib_dict)
        if type(bib_dict) is str:
            bib_dict = getBibTex(bib_dict,**kwargs)
            if isinstance(bib_dict,dict) == False: return ''
        else:
            if kwargs.get('verbose',False): print('Input to longRef is neither a bibcode nor a bibTex dictionary')
            return ''

    authors = bib_dict['author'].split(' and ')
    if len(authors) == 1:
        output = '{}'.format(authors[0].replace('~',' '))
    elif len(authors) == 2:
        output = '{} & {}'.format(authors[0].replace('~',' '),authors[1].replace('~',' '))
    elif len(authors) == 3:
        output = '{}, {} & {}'.format(authors[0].replace('~',' '),authors[1].replace('~',' '),authors[2].replace('~',' '))
    else:
        output = '{}, {}, {}, et al'.format(authors[0].replace('~',' '),authors[1].replace('~',' '),authors[2].replace('~',' '))

# fill in missing data
    if 'year' not in bib_dict.keys():
        bib_dict['year'] = ''
    if 'title' not in bib_dict.keys():
        bib_dict['title'] = ''
    if 'journal' not in bib_dict.keys():
        bib_dict['journal'] = ''
    if 'volume' not in bib_dict.keys():
        bib_dict['volume'] = ''
    if 'pages' not in bib_dict.keys():
        bib_dict['pages'] = ''

    return output+' "{}". {}, {}, {} ({})'.format(bib_dict['title'],bib_dict['journal'],bib_dict['volume'],bib_dict['pages'],bib_dict['year'])


def verylongRef(bib_dict,**kwargs):
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
        if type(bib_dict) is numpy.str:
            bib_dict = str(bib_dict)
        if type(bib_dict) is str:
            bib_dict = getBibTex(bib_dict,**kwargs)
            if isinstance(bib_dict,dict) == False: return ''
        else:
            if kwargs.get('verbose',False): print('Input to verylongRef is neither a bibcode nor a bibTex dictionary')
            return ''

    authors = bib_dict['author'].split(' and ')
    if len(authors) == 1:
        output = '{}'.format(authors[0].replace('~',' '))
    elif len(authors) == 2:
        output = '{} & {}'.format(authors[0].replace('~',' '),authors[1].replace('~',' '))
    else:
        output=''
        for a in authors[:-3]:
            output+='{}, '.format(a.replace('~',' '))
        output+='{} & {}'.format(authors[-2].replace('~',' '),authors[-1].replace('~',' '))

# fill in missing data
    if 'year' not in bib_dict.keys():
        bib_dict['year'] = ''
    if 'title' not in bib_dict.keys():
        bib_dict['title'] = ''
    if 'journal' not in bib_dict.keys():
        bib_dict['journal'] = ''
    if 'volume' not in bib_dict.keys():
        bib_dict['volume'] = ''
    if 'pages' not in bib_dict.keys():
        bib_dict['pages'] = ''

    return output+' "{}". {}, {}, {} ({})'.format(bib_dict['title'],bib_dict['journal'],bib_dict['volume'],bib_dict['pages'],bib_dict['year'])


def citeURL(bib_dict,**kwargs):
    '''
    :Purpose:
        Generate the URL corresponding to a citation, based on the bibcode and NASA ADS syntax

    :Required parameters:
        :param bib_tex: Dictionary output from bibTexParser, else a bibcode that is fed into bibTexParser

    :Optional parameters:
        None

    :Output:
        A string of the format ``Burgasser, A. J., Cruz, K. L., Cushing, M., et al. SpeX Spectroscopy of Unresolved Very Low Mass Binaries. 
        I. Identification of 17 Candidate Binaries Straddling the L Dwarf/T Dwarf Transition. ApJ 710, 1142 (2010)``

    '''
    if type(bib_dict) is not dict:
        if type(bib_dict) is numpy.str:
            bib_dict = str(bib_dict)
        if type(bib_dict) is str:
# assume this is a bibcode
            return 'http://adsabs.harvard.edu/abs/{}'.format(bib_dict)
        else:
            raise NameError('Input to citeURL is neither a bibcode nor a bibTex dictionary')

    else:
        if 'bibcode' in list(bib_dict.keys()):
            return 'http://adsabs.harvard.edu/abs/{}'.format(bib_dict['bibcode'])
        else:
            raise NameError('BibTex dictionary does not contain a bibcode')



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
        biblibrary = kwargs.get('biblibrary', SPLAT_PATH+DB_FOLDER+BIBFILE)
# check the file
        if not os.path.exists(biblibrary):
            if kwargs.get('verbose',True) == True: print('Could not find bibtex library {}'.format(biblibrary))
            biblibrary = SPLAT_PATH+DB_FOLDER+BIBFILE

        if not os.path.exists(biblibrary):
            raise NameError('Could not find SPLAT main bibtext library {}; something is wrong'.format(biblibrary))


        with open(biblibrary, 'r') as bib_file:
            text = bib_file.read()
            #print re.search('@[A-Z]+{' + bib_code, bib_file)        
            in_lib = re.search('@[a-z]+{' + bibcode, text)
            if in_lib == None:  
                if kwargs.get('force',False): return False
                if kwargs.get('verbose',False) == True: print('Bibcode {} not in bibtex library {}; checking online'.format(bibcode,biblibrary))
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
    if not checkOnline():
        return False

    url_begin = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode="
    url_end = "&data_type=BIBTEX"
    url = url_begin + bibcode + url_end
    bib_tex = requests.get(url).content
    
    # Check if content is in html which means bad bib_code was given
    if isinstance(bib_tex,bytes):
        bib_tex = bib_tex.decode()
    if "<HTML>" in bib_tex:
        print('{} is not a valid online bib code.'.format(bibcode))
        return False       
        
    # Cut off extraneous info from website before the bibtex code
    else:
        begin = bib_tex.find('@')
        bib_tex = bib_tex[begin:]
        return bib_tex


