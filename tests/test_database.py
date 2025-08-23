# -*- coding: utf-8 -*-
from __future__ import print_function

# this is the test function set for splat photometry functions

# imports - internal
import copy
import glob
import os
# # imports - external
import numpy
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units
from astropy import coordinates as coord      # coordinate conversion
from astropy.io import fits
from numpy.testing import assert_allclose

# splat functions and constants
import splat
import splat.database as spdb
from splat.initialize import *

#####################
# TESTING FUNCTIONS #
#####################

def test_fetchDatabase():
# fetchDatabase
    pass

def test_queryVizier():
# queryVizier
# getPhotometry
    pass

def test_querySimbad():
# querySimbad
# _querySimbad2
    pass

def test_importSpectra():
# importSpectra
    pass


# def OLD_test_baseline():
#     basefolder = '/Users/adam/projects/splat/exercises/ex9/'
#     sp = splat.getSpectrum(shortname='1047+2124')[0]        # T6.5 radio emitter
#     spt,spt_e = splat.classifyByStandard(sp,spt=['T2','T8'])
#     teff,teff_e = splat.typeToTeff(spt)
#     sp.fluxCalibrate('MKO J',spl.typeToMag(spt,'MKO J')[0],absolute=True)
#     table = splat.modelFitMCMC(sp, mask_standard=True, initial_guess=[teff, 5.3, 0.], zstep=0.1, nsamples=100,savestep=0,filebase=basefolder+'fit1047',verbose=True)


# def OLD_test_ingest(folder='./',**kwargs):
#     splat.importSpectra(data_folder=folder,**kwargs)

# def OLD_test_combine():
# # source db
#     data_folder = '/Users/adam/projects/splat/adddata/daniella/spex_prism_160218/'
#     review_folder = '/Users/adam/projects/splat/adddata/review/'
#     t_src = spdb.fetchDatabase(review_folder+'source_update.csv',csv=True)
# # convert all t_src columns to the same format in DB_SOURCES
#     for col in t_src.colnames:
#         tmp = t_src[col].astype(DB_SOURCES[col].dtype)
#         t_src.replace_column(col,tmp)
#     t_merge = vstack([DB_SOURCES,t_src])
#     t_merge.sort('SOURCE_KEY')
#     t_merge.write(review_folder+DB_SOURCES_FILE,format='ascii.tab')

# # spectrum db
#     t_spec = fetchDatabase(review_folder+'spectrum_update.csv',csv=True)
# # move files
# # WARNING - ASSUMING THESE ARE FITS; NEED A FIX IF THEY ARE NOT
# # COULD JUST READ IN TO SPECTRUM OBJECT AND OUTPUT AS FITS FILE
#     for i,file in enumerate(t_spec['DATA_FILE']):
#         t_spec['DATA_FILE'][i] = '{}_{}.fits'.format(t_spec['DATA_KEY'][i],t_spec['SOURCE_KEY'][i])
#         if t_spec['PUBLISHED'][i] == 'Y':
#             os.copyfile(data_folder+file,'{}/published/{}'.format(review_folder,t_spec['DATA_FILE'][i]))
#             print('Moved {} to {}/published/'.format(t_spec['DATA_FILE'][i],review_folder))
#         else:
#             os.copyfile(data_folder+file,'{}/unpublished/{}'.format(review_folder,t_spec['DATA_FILE'][i]))
#             print('Moved {} to {}/unpublished/'.format(t_spec['DATA_FILE'][i],review_folder))
# # convert all t_src columns to the same format in DB_SOURCES
#     for col in t_spec.colnames:
#         tmp = t_spec[col].astype(DB_SPECTRA[col].dtype)
#         t_spec.replace_column(col,tmp)
#     t_merge = vstack([DB_SPECTRA,t_spec])
#     t_merge.sort('DATA_KEY')
#     t_merge.write(review_folder+DB_SPECTRA_FILE,format='ascii.tab')
#     print('\nDatabases updated; be sure to move these from {} to {}{}'.format(review_folder,CODE_PATH,DB_FOLDER))
#     print('and to move spectral files from {}/published and {}/unpublished/\n'.format(review_folder,review_folder))




