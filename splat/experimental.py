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