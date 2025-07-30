"""
This subpackage contains modules and packages for testing the SPLAT code
"""
import numpy
import os
import copy
import glob
from astropy import units as u            # standard units
from astropy import constants as const        # physical constants in SI units
from astropy import coordinates as coord      # coordinate conversion
from astropy.io import fits
from numpy.testing import assert_allclose
