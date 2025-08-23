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




def bibTexParser(bib_input,**kwargs):
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
    bib_dict = {'bib_tex': bib_input}
    bib_tex = bib_input.split(',\n')
    # get bib code
    bib_dict['type'] = bib_tex[0][1:bib_tex[0].find('{')]
    bib_dict['bibcode'] = bib_tex[0][bib_tex[0].find('{')+1:]
#    begin = bib_tex.find('{')  
#    end = bib_tex.find(',')
#    bib_dict['bibcode'] = bib_tex[begin+1:end]
#    bib_tex = bib_tex[end+1:]   # remove bib code line
    
#    bib_tex =  bib_tex.split(',\n')  # this moght not always work for author lists
    
    for line in bib_tex[1:]:
        line = line.strip()
        line = line.replace('{','').replace('}','').replace('\"','').replace('\n@','').replace('\n','').replace('\t','') 
        line = line.split('=')
# TEMPORARY FIX AS THIS LINE ISN'T WORKING RIGHT
        if len(line) > 1:
            line[0] = line[0].strip().lower()
            line[1] = line[1].strip()
            bib_dict[line[0]] = line[1]

# Journal massaging: look up table
    if 'journal' in list(bib_dict.keys()):
        if bib_dict['journal'][1:].lower() in list(JOURNALS_LONGNAMES.keys()): bib_dict['journal'] = JOURNALS_LONGNAMES[bib_dict['journal'][1:].lower()]
        elif bib_dict['journal'][2:].lower() in list(JOURNALS_LONGNAMES.keys()): bib_dict['journal'] = JOURNALS_LONGNAMES[bib_dict['journal'][2:].lower()]
        else: pass

        
    return bib_dict


def veryShortRef(bib_dict,**kwargs):
    '''
    :Purpose:
        Takes a bibtex entry and returns a short (in-line) version of the citation

    :Required parameters:
        :param bib_tex: Dictionary output from bibTexParser, else a bibcode that is fed into bibTexParser

    :Optional parameters:
        None

    :Output:
        A string of the format ``Burgasser et al. (2006)``

    '''
    if type(bib_dict) is not dict:
        if type(bib_dict) is str: bib_dict = getBibTeX(bib_dict,**kwargs)
        if len(bib_dict) == 0: return ''
    if type(bib_dict) is not dict:
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
        if type(bib_dict) is str: 
            bib_dict = getBibTeX(bib_dict,**kwargs)
            if len(bib_dict) == 0: return ''
    if type(bib_dict) is not dict:
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
        if type(bib_dict) is str:
            bib_dict = getBibTeX(bib_dict,**kwargs)
            if len(bib_dict) == 0: return ''
    if type(bib_dict) is not dict:
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


def veryLongRef(bib_dict,**kwargs):
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
            bib_dict = getBibTeX(bib_dict,**kwargs)
            if len(bib_dict) == 0: return ''
    if type(bib_dict) is not dict:
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
        if type(bib_dict) is str:
            bib_dict = str(bib_dict)
        if type(bib_dict) is str:
# assume this is a bibcode
            return '{}{}/abstract'.format(CITATION_URL_BASE,bib_dict)
        else:
            raise NameError('Input to citeURL is neither a bibcode nor a bibTex dictionary')

    else:
        if 'bibcode' in list(bib_dict.keys()):
            return '{}{}/abstract'.format(CITATION_URL_BASE,bib_dict['bibcode'])
        else:
            raise NameError('BibTex dictionary does not contain a bibcode')

def processBiblibrary(biblibrary,verbose=False):
    '''
    Purpose
        Processes a bibtex .bib library (multiple bibtex entries) into a dictionary whose keys
        are the bibcode


    :Required parameters:
        :param biblibrary: .bib file containing the bibtex entries

    :Optional parameters:
        :param: verbose = False: Set to True to provide extensive feedback

    :Output:
        A dictionary containing the bibtex information, organized by bibcode key

    '''    
    if not os.path.exists(os.path.normpath(biblibrary)):
        raise ValueError('Could not find bibtex library file {}'.format(biblibrary))

    with open(os.path.normpath(biblibrary), 'r') as bib_file:
        text = bib_file.read()

# find all of the bibtex codes
    output = {}
    flg = 'upper'
    in_lib = re.search('@[A-Z]+{', text)
    if in_lib==None:
        flg = 'lower'
        in_lib = re.search('@[a-z]+{', text)
        if in_lib==None:
            raise ValueError('Cannot find any bib entries in text {}'.format(text[:1000]))
    while in_lib != None:
        if flg=='upper': in_lib = re.search('@[A-Z]+{', text)
        else: in_lib = re.search('@[a-z]+{', text)
        asc = text[in_lib.start():]
        in_lib = re.search('\n@', asc)
        if in_lib != None:
            text = asc[(in_lib.start()-2):]
            asc = asc[:in_lib.end()]
        p = bibTexParser(asc)
        output[p['bibcode']] = p
    return output


def getBibTeX(bibcode, biblibrary=BIBFILE, online=False, verbose=True):
    '''
    Purpose
    -------
    Takes a bibcode and returns a dictionary containing the bibtex information; 
    looks either in internal SPLAT or user-supplied bibfile, or seeks online. 
    If nothing found, gives a soft warning and returns False

    Parameters
    ----------

    bibcode : str
        Bibcode string to look up (e.g., '2014ApJ...787..126L')

    biblibrary = splat.BIBFILE: str [optional]
        File pointing to a bibtex library file; by default points to internal library

    online = False : bool [optional]
        If True, go directly online; if False, do not try to go online 
        NOTE: CURRENLY SET TO NOT ONLINE DUE TO CHANGE IN ADS API

    verbose = True : bool [optional]
        Set to True to provide feedback

    Outputs
    -------

    dictionary containing bibtex information, or blank dictionary if nothing found

    Example
    -------
    
    TBD

    Dependencies
    ------------
    
    None

    '''

# go online first if directed to do so
    # if online==True and checkOnline():
    #     bib_tex = getBibTeXOnline(bibcode)

# read locally first
#    else:
# check the file
    if not os.path.exists(os.path.normpath(biblibrary)):
        if verbose == True: print('Could not find bibtex library {}'.format(biblibrary))
        biblibrary = os.path.join(CITATION_RESOURCES_FOLDER,BIBFILE)
    if not os.path.exists(os.path.normpath(biblibrary)):
        raise NameError('Could not find SPLAT main bibtext library {}; something is wrong'.format(biblibrary))

# open and read
    bib_tex = {}
    with open(os.path.normpath(biblibrary), 'r') as bib_file:
        text = bib_file.read()
        #print re.search('@[A-Z]+{' + bib_code, bib_file)        
        in_lib = re.search('@[a-z]+{' + bibcode, text)
        if in_lib != None:  
#            if force==True: return bib_tex
#            if kwargs.get('verbose',False) == True: print('Bibcode {} not in bibtex library {}; checking online'.format(bibcode,biblibrary))
#            bib_tex = getBibTeXOnline(bibcode)
#        else:
            begin = text.find(re.search('@[a-z]+{' + bibcode, text).group(0))
            text = text[begin:]
            end = text.find('\n@')
            bib_tex = text[:end]

    if len(bib_tex) == 0: return bib_tex
    else: return bibTexParser(bib_tex)


def getBibTeXOnline(bibcode,verbose=False):
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
        if verbose==True: print('{} is not a valid online bib code.'.format(bibcode))
        return False       
        
    # Cut off extraneous info from website before the bibtex code
    else:
        begin = bib_tex.find('@')
        bib_tex = bib_tex[begin:]
        return bib_tex

def nsfbib(biblibrary,file=''):
    '''
    Purpose
        Takes a biblibrary and generates an NSF-formatted reference list, based on the following requirements:
        "Each reference must include the names of all authors (in the same sequence in which they 
        appear in the publication), the article and journal title, book title, volume number, page numbers, a
        and year of publication. If the proposer has a website 
        address readily available, that information should be included in the citation."

    :Required parameters:

        :param biblibrary: .bib file containing the bibtex entries

    :Optional parameters:
        :param file: Filename to save latex output (default='citations.tex') 

    :Output:
        A latex file containing relevant bibliographic information

    '''

    try:
        cites = processBiblibrary(biblibrary,verbose=True)
    except:
        raise ValueError('Could not parse biblibrary input {}; make sure this a full path to the .bib file'.format(biblibrary))

# process into a set of strings
    output = []
    for c in list(cites.keys()):
        line = cites[c]['author']
        li = line.rsplit(' and ',1)
        line = ', & '.join(li)
        line = line.replace(' and ',', ').replace('~','')
        line = line+' "'+cites[c]['title']+'." '+cites[c]['year']
# article?  
        if cites[c]['type'] == 'article':
            if 'journal' in list(cites[c].keys()): line=line+', '+cites[c]['journal']
            if 'volume' in list(cites[c].keys()): line=line+', '+cites[c]['volume']
            if 'pages' in list(cites[c].keys()): line=line+', '+cites[c]['pages']
# book?  
        if cites[c]['type'] == 'book':
            if 'publisher' in list(cites[c].keys()): line=line+', '+cites[c]['publisher']
# add url
        if 'adsurl' in list(cites[c].keys()): line=line+' ('+cites[c]['adsurl']+')'
        elif 'bdsk-url-1' in list(cites[c].keys()): line=line+' ('+cites[c]['bdsk-url-1']+')'
        else: pass
        output.append(line)
    output.sort()

# save to file if provided
    if file != '':
        try:
            f = open(os.path.normpath(file),'w')
            for o in output: f.write(o+'\n')
            f.close()
            return True
        except:
            print('Warning: problem saving to output file {}'.format(file))

    line = ''
    for o in output: line=line+o+'\n'
    return line



