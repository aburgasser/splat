#!/usr/bin/python -tt

from __future__ import absolute_import

import sys
import os
import numpy
import scipy
import astropy
from astropy.io import ascii, fits			# for reading in spreadsheet
from astropy.table import Table, join			# for reading in table files
import matplotlib.pyplot as plt
from scipy.integrate import trapz		# for numerical integration
from scipy.interpolate import interp1d
import re
import urllib2
from astropy import units as u			# standard units
from astropy.coordinates import ICRS, Galactic		# coordinate conversion

