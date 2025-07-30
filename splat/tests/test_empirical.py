# -*- coding: utf-8 -*-
from __future__ import print_function

# this is the test function set for splat empirical functions

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
import splat.empirical as spem

# things to test
# isNumber 


#####################
# TESTING FUNCTIONS #
#####################


def test_typeToMag():
    pass

def test_typeToColor():
    pass

def test_typeToTeff():
    t = spem.typeToTeff('M7')
    assert len(t) == 2
    assert isinstance(t[0],u.quantity.Quantity)
    assert isinstance(t[1],u.quantity.Quantity)
    assert isinstance(t[0].value,float)
    assert isinstance(t[1].value,float)
    t = spem.typeToTeff(17)
    assert len(t) == 2
    assert isinstance(t[0],u.quantity.Quantity)
    assert isinstance(t[1],u.quantity.Quantity)
    assert isinstance(t[0].value,float)
    assert isinstance(t[1].value,float)
    t = spem.typeToTeff(['M7','L4','T5'])
    assert len(t) == 2
    assert len(t[0]) == 3
    assert len(t[1]) == 3
    assert isinstance(t[0],u.quantity.Quantity)
    assert isinstance(t[1],u.quantity.Quantity)
    assert isinstance((t[0].value)[0],float)
    assert isinstance((t[0].value)[1],float)
    assert isinstance((t[0].value)[2],float)
    assert isinstance((t[1].value)[0],float)
    assert isinstance((t[1].value)[1],float)
    assert isinstance((t[1].value)[2],float)
    t = spem.typeToTeff([17,24,35])
    assert len(t) == 2
    assert len(t[0]) == 3
    assert len(t[1]) == 3
    assert isinstance(t[0],u.quantity.Quantity)
    assert isinstance(t[1],u.quantity.Quantity)
    assert isinstance((t[0].value)[0],float)
    assert isinstance((t[0].value)[1],float)
    assert isinstance((t[0].value)[2],float)
    assert isinstance((t[1].value)[0],float)
    assert isinstance((t[1].value)[1],float)
    assert isinstance((t[1].value)[2],float)    
    t = spem.typeToTeff(17,uncertainty=0.5)
    assert len(t) == 2
    assert isinstance(t[0],u.quantity.Quantity)
    assert isinstance(t[1],u.quantity.Quantity)
    assert isinstance(t[0].value,float)
    assert isinstance(t[1].value,float)
    t = spem.typeToTeff([17,24,35],uncertainty=0.5)
    assert len(t) == 2
    assert len(t[0]) == 3
    assert len(t[1]) == 3
    assert isinstance(t[0],u.quantity.Quantity)
    assert isinstance(t[1],u.quantity.Quantity)
    assert isinstance((t[0].value)[0],float)
    assert isinstance((t[0].value)[1],float)
    assert isinstance((t[0].value)[2],float)
    assert isinstance((t[1].value)[0],float)
    assert isinstance((t[1].value)[1],float)
    assert isinstance((t[1].value)[2],float)    
    t = spem.typeToTeff([17,24,35],uncertainty=[0.5,1.0,0.8])
    assert len(t) == 2
    assert len(t[0]) == 3
    assert len(t[1]) == 3
    assert isinstance(t[0],u.quantity.Quantity)
    assert isinstance(t[1],u.quantity.Quantity)
    assert isinstance((t[0].value)[0],float)
    assert isinstance((t[0].value)[1],float)
    assert isinstance((t[0].value)[2],float)
    assert isinstance((t[1].value)[0],float)
    assert isinstance((t[1].value)[1],float)
    assert isinstance((t[1].value)[2],float)    
    t = spem.typeToTeff(1875,reverse=True)
    assert len(t) == 2
    assert isinstance(t[0],float)
    assert isinstance(t[1],float)
    t = spem.typeToTeff([1400,1875,2482],reverse=True)
    assert len(t) == 2
    assert isinstance(t[0],list)
    assert isinstance(t[1],list)
    assert isinstance(t[0][0],float)
    assert isinstance(t[0][1],float)
    assert isinstance(t[0][2],float)
    assert isinstance(t[1][0],float)
    assert isinstance(t[1][1],float)
    assert isinstance(t[1][2],float)
    t = spem.typeToTeff(1875,uncertainty=150.,reverse=True)
    assert len(t) == 2
    assert isinstance(t[0],float)
    assert isinstance(t[1],float)
    t = spem.typeToTeff([1400,1875,2482],uncertainty=100.,reverse=True)
    assert len(t) == 2
    assert isinstance(t[0],list)
    assert isinstance(t[1],list)
    assert isinstance(t[0][0],float)
    assert isinstance(t[0][1],float)
    assert isinstance(t[0][2],float)
    assert isinstance(t[1][0],float)
    assert isinstance(t[1][1],float)
    assert isinstance(t[1][2],float)
    t = spem.typeToTeff([1400,1875,2482],uncertainty=[100.,50.,150.],reverse=True)
    assert len(t) == 2
    assert isinstance(t[0],list)
    assert isinstance(t[1],list)
    assert isinstance(t[0][0],float)
    assert isinstance(t[0][1],float)
    assert isinstance(t[0][2],float)
    assert isinstance(t[1][0],float)
    assert isinstance(t[1][1],float)
    assert isinstance(t[1][2],float)
    t = spem.typeToTeff([1400,1875,2482],uncertainty=[100.,50.,150.],reverse=True,string=True)
    assert len(t) == 2
    assert isinstance(t[0],list)
    assert isinstance(t[1],list)
    assert isinstance(t[0][0],str)
    assert isinstance(t[0][1],str)
    assert isinstance(t[0][2],str)
    assert isinstance(t[1][0],float)
    assert isinstance(t[1][1],float)
    assert isinstance(t[1][2],float)

def test_typeToLuminosity():
# typeToMag
    pass

def test_typeToBC():
# typeToMag
    pass

def test_estimateDistance():
# typeToMag
    pass

def test_redden():
# typeToMag
    pass


