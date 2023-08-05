#!/usr/bin/env python
# coding: utf-8
######################################
###          asteroid.py           ###
######################################
from splat.initialize import *
from splat.core import Spectrum
from splat.plot import plotSpectrum
######################################
############# internal ###############
import shutil
import math
import re
import copy
######################################
############# external ###############
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits
from numpy.linalg import eig
from scipy import stats, interpolate
######################################
# ### version warning
# if sys.version_info.major != 3:
#     raise NameError('\The Asteroids package is built for Python3, please be aware of issues while using a an earlier version.\n')



#####################################################
###########   AsteroidSpectrum class   ##############
#####################################################
class AsteroidSpectrum(Spectrum):
    '''
    As extension of the splat.Spectrum class, the AsteroidSpectrum subclass of objects functions similarily
    to splat.Spectrum objects. With the exception of the "flux" measurement being a unitless reflectance,  
    an AsteroidSpectrum object also represents itself as an "asteroid spectrum" to distinguish itself from
    other spectra.
    
    See also:
    splat.Spectrum
    splat.core.Spectrum
    '''
    def __init__(self, *args, **kwargs):
        super().__init__(dimensionless = True, flux_label = 'Reflectance', *args, **kwargs)
        
    # asteroids represent themselves differently when called
    def __repr__(self):
        return '{} asteroid spectrum of {}'.format(self.instrument,self.name)
    
    def plot(self, **kwargs):
        super().plot(**kwargs)



############################################
#######  AsteroidTemplate class   ##########
############################################
class AsteroidTemplate(dict):
    '''
    The AsteroidTemplate class functions like the python native dictionary, with the difference that it takes a
    path to a folder of templates as an argument. On constructon, AsteroidTemplate becomes a dictionary where the
    keys are template names and the values are the corresponding AsteroidSpectrum objects created using the data
    in the template files.
    '''
    def __init__(self, /, template_directory, *args, **kwargs):
        count = 0
        for filename in os.scandir(template_directory):
            if filename.is_file():
                full_path_name = os.path.abspath(filename)
                spectrum_name = os.path.split(full_path_name)[-1].split('.')[0]
                try:
                    sp = AsteroidSpectrum(full_path_name, name = spectrum_name, instrument = 'Template')
                    super().update({spectrum_name.split('_')[-1] : sp})
                except:    
                    try:
                        dp = pd.read_csv(full_path_name, delimiter='\s+')
                        dp.rename(columns={'Wave':'wave','Flux':'flux','Noise':'noise'},inplace=True)
                        sp = AsteroidSpectrum(full_path_name, name = spectrum_name, instrument = 'Template')
                        super().update({spectrum_name.split('_')[-1] : sp})
                    except:
                        print('Asteroid templates failed to load at ' + full_path_name)
                        break
                else:
                    count += 1
        print("Created dictionary with " + str(count) + " classes from " + os.path.split(os.path.normpath(template_directory))[-1])

        
        
### for storing package relevant data
def _initialize_asteroid_template(template_path):
    '''
    :Purpose:
        Stores the asteroid templates in splat/resources/SpectralTemplates/ to a (global) variable
    '''
    return AsteroidTemplate(template_path)


def _load_demeo_pcs(pcs_path):
    '''
    From Demeo's supplementary materials, this function reads a .txt files containing the eigenvectors and means that 
    DeMeo 2009 used to define the Bus-DeMeo Taxonomy. 
  
    Parameters
    ----------
    pcs_path : string
        file path of the item containing some header text, eigenvectors, and mean values of each wavelength channel
        defining the Demeo classification scheme in IR.
    
    Returns
    -------
    Two arrays containing:
        demeo_pc1 : 
            vector for principle component 1
            
        demeo_pc2 : 
            vector for principle component 2
            
        demeo_pc3 : 
            vector for principle component 3 
            
        demeo_pc4 : 
            vector for principle component 4
            
        demeo_pc5 : 
            vector for principle component 5
            
    demeo_chn_means : array
        the mean values used in obtaining the eigenvectors by Demeo
        
    '''
    
    # reading the eigenvectors and mean values txt file
    with open(pcs_path) as f:
        eigv_and_means = str(f.readlines())
    eigv_and_means = eigv_and_means.split()

    # removing the headers
    eigv_and_means = eigv_and_means[11:]

    demeo_ev_channels = []
    demeo_pc1 = []
    demeo_pc2 = []
    demeo_pc3 = []
    demeo_pc4 = []
    demeo_pc5 = []
    demeo_chn_means = []

    for item in eigv_and_means:
        temp = item.split('\\n') # removing the "new line" before each new channel entry
        temp = temp[0].split('\\t') # removing the tab white spaces between values

        demeo_ev_channels.append(float(temp[0][1:])) # optional, allows us to check the order of the channels being handeled corresponding to the PCs 
        demeo_pc1.append(float(temp[1]))
        demeo_pc2.append(float(temp[2]))
        demeo_pc3.append(float(temp[3]))
        demeo_pc4.append(float(temp[4]))
        demeo_pc5.append(float(temp[5]))
       
        # the mean values of Demeo's data, which was used to obtain the PC's eigenvectors
        # this is necessary in calculating the PC scores for each Spectrum
        demeo_chn_means.append(float(temp[6].split('\'')[0])) 
    
    # wavelength channels
    try:
        demeo_ev_channels = np.array(demeo_ev_channels)
        print('Stored principal components from ' + os.path.split(pcs_path)[-1])
    except:
        print('Unable to store principal components from ' + os.path.split(pcs_path)[-1])
                                                                                     
    # means of each channel
    demeo_chn_means = np.array(demeo_chn_means)
    
    return np.array([demeo_pc1,demeo_pc2,demeo_pc3,demeo_pc4,demeo_pc5]),demeo_chn_means



### for extracting names from spex_prism items (after template comparison with classifyFolder)
# def _shorten_result_keys(results_dictionary):
#     shortname_results_dictionary = {}
#     counter = 1
#     for key, value in results_dictionary.items():
#         temp_list = key.split('__')
#         name = temp_list[1]

#         # check if there are multiple spectra of the same object with a classification result already in the dictionary
#         # the shortname of a spectra of some object beyond the first receives an index, beginning at _1, then _2, _3, etc.
#         if(shortname_results_dictionary.get(name) == None):
#             shortname_results_dictionary[name] = value
#             counter = 1
#         else:
#             shortname_results_dictionary[name + '_' + str(counter)] = value
#             counter += 1
#     return shortname_results_dictionary



DEMEO_TEMPLATES_PATH = SPLAT_PATH + '/resources/SpectralTemplates/DemeoAsteroids'
DEMEO_PCS_PATH = SPLAT_PATH + '/resources/SpectralTemplates/DemeoAsteroids/pca/demeo_ir_eigenvectors_and_channel_means.txt'

## These are the words we ignore when reading the fits header 'OBJECT' field
IGNORED_FROM_HEADER = ['asteroid', 'Asteroid']

## Taxonomic template dictionaries
DEMEO_CLASSES = _initialize_asteroid_template(DEMEO_TEMPLATES_PATH)
# FUTURE_TAXON = _initialize_asteroid_template(TAXON_FOLDER_PATH)

# obtaining the principal component eivenvectors of the classification scheme from a .txt file
PC_VECTORS, CHANNEL_MEANS = _load_demeo_pcs(DEMEO_PCS_PATH)



#####################################################
##########    Main Classifying Method     ###########
#####################################################

##########    PCA classifying routine      ##########
def classifyAsteroid(spectrumObject, *args, **kwargs):
    '''
    Classification of an asteroid's spectrum using principal component analysis or direct comparison with templates.
    
    Parameters
    ----------
    spectrumObject : Spectrum object
    
    method : string, default = 'PCA'
        specifies the method used to classify an AsteroidSpectrum. 'PCA': the Bus-DeMeo Taxonomy
    
    scores : bool, default = False
        returns a dictionary instead of only the classification verdict.
        Dictionary :
        {'Name': str(asteroid name), 'Taxon': str(verdict), 'PC Scores': list(PCA scores), 'Slope': float(slope)}
    
    waverange : array of two floats, default = [0.85, 2.45]
        specifies wavelength range to consider [start, end).
    
    norm : float, default = 1.2
        wavelength at which the spline fit is normalized, in microns
    
    plot : bool, default = False
        plots the spectrum data after normalization of spline fit at norm microns 
        and the reconstruction using principal components. verbose = True : plots the spectrum as
        it is transformed by different steps during classification.
    
    print_result : bool, default = False
        Shows the name, verdict, PCA scores and slope of the object that was classified.
        
    verbose : bool, default = False
        prints the data points of the spectrum data after spline fit and normalization 
        at wavelength (µm) defined by the "norm" variable, and the corresponding reconstruction
        using principal components
    
    force_classify : bool, default = False
        asteroid spectrum that does not span the boundaries defined by the "waverange" variable
        will by default be indicated with a verdict of '-1'. Toggling foce_run to True will
        instead try to classify the spectrum regardless of its shorter wavelength range.
        
    Returns
    -------
    verdict : String
        the classification of the Spectrum object
    
    If scores was toggled True : Dictionary
        {'Name': str(asteroid name), 'Taxon': str(verdict), 'PC Scores': list(PCA scores), 'Slope': float(slope)}
    Example
    -------
    >>> ex1_spec = MITHNEOS_to_Spectrum('a053435.sp272.txt')
    >>> classifyAsteroid(ex1_spec, plot = True)
    
    See Also
    --------
    AsteroidSpectrum
    visnir_to_Spectrum
    MITHNEOS_to_Spectrum
        
    '''
    method = kwargs.get('method', 'PCA')
    waverange = kwargs.get('waverange', [0.85, 2.45])
    norm = kwargs.get('norm', 1.2)
    plot = kwargs.get('plot', False)
    print_result = kwargs.get('print_result', False)
    verbose = kwargs.get('verbose', False)
    scores = kwargs.get('scores', False)
    force_classify = kwargs.get('force_classify', False)
    
    # storing a copy of the input splat.Spectrum Object
    sp = copy.deepcopy(spectrumObject)
    
    # preparing the object, cubic spline, dividing out the slope
    sp_reduced, slope_value, slope_line = form_data_pairs(sp, waverange, norm)
    sp_reflectance = sp_reduced[0][1]
    # subtracting the reflectance by mean values from Demeo's taxonomy definition
    subtracted_sp_reflectance = sp_reflectance - CHANNEL_MEANS
    # calculate the PC scores from the PC1- PC5
    pc_scores = calc_scores(PC_VECTORS, subtracted_sp_reflectance)

    
    # if spectrum wavelength range is too short
    if not(sp.wave[0].value <= waverange[0] and sp.wave[-1].value >= waverange[1]):
        if print_result:
            print(sp.name)
            print('\033[1mWarning\033[0m, spectrum ' + sp.name + ' does not span specified waverange ' + str(waverange) + 'µm.')
        
        # not proceeding with PCA, all returns are defaulted to -1 here
        if not force_classify:
            verdict = '-1'
    
            if print_result:
                print('Slope: ' + str(slope_value))
                print('PC Scores:', pc_scores)
                print('Verdict:', verdict)
            
            if scores:
                pc_scores = [-1, -1, -1, -1, -1]
                slope_value = -1
                return {'Name': sp.name, 'Taxon' : verdict, 'PC Scores': pc_scores, 'Slope': slope_value}
            else:
                return verdict
            
        elif force_classify: # apply PCA to short spectrum
            verdict = classifyStep1(pc_scores, slope_value)
        
    if plot:
        xyzip = zip(sp.wave.value, sp.flux.value)
        temp_flux = sp.flux.value
        for wave, flux in xyzip:
            if wave == norm : 
                # attempt to normalize the original spectrum flux @ "norm" µm, if measurement @ normµm is available
                temp_flux = temp_flux/flux
                break
            elif abs(1 - (norm/wave))*100 <= 0.5:
                temp_flux = temp_flux/flux
                break
            
        # recreate spectrum flux using PCA results, WITH SLOPE
        reconstructed_flux = (np.dot(pc_scores, PC_VECTORS) + CHANNEL_MEANS)*slope_line
        # the complete reconstruction, sp_re AsteroidSpectrum, should look like a smoothed version of the original spectrum
        sp_re = AsteroidSpectrum(wave = sp_reduced[0][0], flux = reconstructed_flux)
        # spline fit of the original spectrum (user input), normalized at norm microns WITH SLOPE
        sp_norm_spline_slope = AsteroidSpectrum(wave = sp_reduced[0][0], flux=sp_reflectance*slope_line)
        
        # scaling_factor = reconstructed_flux[0]/temp_flux[0]
        # temp_flux = temp_flux*scaling_factor
        # the original spectrum scaled to display next to the normalized reconstruccted spectrum
        sp_scaled = AsteroidSpectrum(wave = sp.wave.value, flux = temp_flux)
        
        plot_title = 'Classification of ' + sp.name
        # show plot of every step of the classification: original > spline fit > PCA (reconstructed to show accuracy)
        if verbose:
            plot_legend = ['Original normalized to unity at ' + str(norm) + 'µm' ,'Spline fit normalized to unity at ' + str(norm) + 'µm','Reconstruction using components', 'Difference between spline and reconstruction']
            plot_colors = ['tab:blue','tab:red','tab:orange', 'tab:gray']
            plotSpectrum([sp_scaled, sp_norm_spline_slope,sp_re, sp_norm_spline_slope - sp_re], title =plot_title, color =plot_colors, legend =plot_legend) # plots 4 things
        else:
            # plotting only the original and PCA reconstruction using the scores calculated
            plot_legend = ['Original', 'PCA Reconstruction']
            plot_colors = ['tab:blue', 'tab:orange']
            plotSpectrum([sp_scaled, sp_re], title =plot_title, color =plot_colors, legend =plot_legend) # plots 2 things

    # the classifying step(s)    
    verdict = classifyStep1(pc_scores, slope_value)
        
    if verbose:    
        reconstructed_flux = (np.dot(pc_scores, PC_VECTORS) + CHANNEL_MEANS)*slope_line
        print("Flux after normalizing the cubic spline fit at " + str(norm) + ": \n", sp_reflectance*slope_line)
        print("Reconstructed flux: \n", reconstructed_flux)
    
    if print_result:
        print(sp.name)
        print('Slope: ' + str(slope_value))
        print('PC Scores:', pc_scores)
        print('Verdict:', verdict)

    if scores:
        return {'Name': sp.name, 'Taxon' : verdict, 'PC Scores': pc_scores, 'Slope': slope_value}
    else:
        return verdict

# def _PCA_classifying()
#     todo: put _PCA into one method, and _template into another. classifyAsteroid calls either one
#     depending on the "method" kwarg

    
def classifyStep1(pc_scores, slope):
    '''
    Step 1 of the classification process, proceeding to step 2 if necessary.
    
    Parameters
    ----------
    pc_scores : array
        principal component scores of the object
    
    slope : float
        slope of the object
        
    Returns
    -------
    verdict : string
        the classification of the Spectrum object
        
    See Also
    --------
    classifyAsteroid
    calc_scores
    form_data_pairs
    classifyStep2
    
    Notes
    -----
    This function is called internally by other functions and is not intended for the end user.
    
    IR step 1: End members
    
    PCir3 = PCir2 − 0.08  Line 1  
    PCir1 = PCir2 + 0.15  Line 2  
    PCir1 = PCir2 − 0.10  Line 3  
    PCir1 = PCir2 − 0.40  Line 4 
    '''
    PC1, PC2, PC3, PC4, PC5 = pc_scores
    
    line1 = PC2 - 0.08
    line2 = PC2 + 0.15
    line3 = PC2 - 0.10
    line4 = PC2 - 0.40
    
    if(PC3 >= line1):
        on_or_above_line1 = True
    else:
        on_or_above_line1 = False
        
    if(PC1 >= line2):
        on_or_above_line2 = True
    else:
        on_or_above_line2 = False
    
    if(PC1 >= line3):
        on_or_above_line3 = True
    else:
        on_or_above_line3 = False
    
    if(PC1 >= line4):
        on_or_above_line4 = True
    else:
        on_or_above_line4 = False
    
    line_positions = [on_or_above_line1, on_or_above_line2, on_or_above_line3, on_or_above_line4]
    
    # End Members
    if(PC1 >= 0.5):
        return 'V'
    elif((0.29 <= PC1 < 0.5) and (PC5 <= 0.05)):
        return 'Sv, Sr'
    elif((PC2 <= -0.5) and (PC4 >= 0.15) and (-0.40 < PC1 <= 0)):
        return 'O'
    elif((0.25 <= PC2 < 0.5) and (PC5 >= 0.06) and (PC3 >= 0.05)):
        return 'R'
    elif(not on_or_above_line1 and (slope >= 0.24)):
        return 'D'
    elif((PC1 <= -0.4) and (PC2 <= -0.2) and (PC4 >= -0.07) and (slope >= 0.5) and (PC3 >= 0)):
        return 'A'
    elif(PC1 <= -0.4):
        return 'Sa'
    else:
        return classifyStep2(pc_scores, line_positions)


    
def classifyStep2(pc_scores, line_positions):
    '''
    Step 2 of the classification process, proceeding to step 3 if necessary.
    
    Parameters
    ----------
    pc_scores : array
        principal component scores of the object
    
    line_positions : array
        conditions (booleans) corresponding to the position of the object in component space 
        
    Returns
    -------
    verdict : string
        the classification of the Spectrum object
        
    See Also
    --------
    classifyAsteroid
    calc_scores
    form_data_pairs
    classifyStep1
    classifyStep3
    
    Notes
    -----
    This function is called internally by other functions and is not intended for the end user.
    
    IR step 2: S-complex

    PCir3 = PCir2 − 0.08  Line 1  
    PCir1 = PCir2 + 0.15  Line 2  
    PCir1 = PCir2 − 0.10  Line 3  
    PCir1 = PCir2 − 0.40  Line 4 
    '''
    PC1, PC2, PC3, PC4, PC5 = pc_scores
    on_or_above_line1, on_or_above_line2, on_or_above_line3, on_or_above_line4 = line_positions
    
    line3 = PC2 - 0.10
    line4 = PC2 - 0.40
    
    if(PC1 == line3):
        on_line3 = True
    else:
        on_line3 = False
    
    if(PC1 == line4):
        on_line4 = True
    else:
        on_line4 = False
    
    if(on_or_above_line1 and on_or_above_line2): # on or above line 1 and 2
        return 'S, Sr, Sq, Q'
    elif(on_or_above_line1 and (not on_or_above_line2 and on_or_above_line3)): # on or above line 1 and between line 2 and line 3, line 2 is always > line 3
        return 'S, Sq, Q, L, K'
    elif(on_or_above_line1 and ((on_or_above_line4 and not on_or_above_line3) or (on_line3) or (on_line4))): # on or above line 1 and on or between line 3 and 4
         return 'K, L, Sq'
    else:
        return classifyStep3(pc_scores, line_positions)


    
def classifyStep3(pc_scores, line_positions):
    '''
    Step 3 of the classification process, only called of necessary.
    
    Parameters
    ----------
    pc_scores : array
        principal component scores of the object
    
    line_positions : array
        conditions (booleans) corresponding to the position of the object in component space 
        
    Returns
    -------
    verdict : string
        the classification of the Spectrum object
        
    See Also
    --------
    classifyAsteroid
    get_scores
    form_data_pairs
    classifyStep1
    classifyStep2
    
    Notes
    -----
    This function is called internally by other functions and is not intended for the end user.
    
    IR step 3: C- and X-complexes
    
    PCir3 = PCir2 − 0.08  Line 1  
    PCir1 = PCir2 + 0.15  Line 2  
    PCir1 = PCir2 − 0.10  Line 3  
    PCir1 = PCir2 − 0.40  Line 4 
    '''
    PC1, PC2, PC3, PC4, PC5 = pc_scores
    on_or_above_line1, on_or_above_line2, on_or_above_line3, on_or_above_line4 = line_positions
    
    line2 = PC2 + 0.15
    line3 = PC2 - 0.10
    line4 = PC2 - 0.40
    
    if(PC1 == line2):
        on_line2 = True
    else:
        on_line2 = False
    
    if(PC1 == line3):
        on_line3 = True
    else:
        on_line3 = False
    
    if(PC1 == line4):
        on_line4 = True
    else:
        on_line4 = False
    
    if((not on_or_above_line1) and ((on_or_above_line4 and not on_or_above_line3) or (on_line3) or (on_line4))): # below line 1 and on or between line 3 and 4
        return 'X-, C-complexes, L, K, T'
    elif((not on_or_above_line1) and (on_or_above_line3 and not on_or_above_line2) and not on_line3): # below line 1 and between 2 and 3
        return 'X-, C-complexes'
    elif((not on_or_above_line1) and (not on_or_above_line4)):
        return 'C, B, L, Cb, X'
    else:
        return 'Indeterminate'


    
def spline_fit(spec, wave_range = [0.85, 2.45]):
    '''
    Perform cubic spline interpolation, the wave axis is separated into 0.05 micron steps.
    
    Parameters
    ----------
    spec : splat.Spectrum object
        The spectra to be processed before principal component analysis.
    
    wave_range : array of two floats, default = [0.85, 2.45]
        specifies wavelength range to consider [start, end).
    
    
    Returns
    -------
    tuple of variables (xnew, ynew, tck)
        xnew : array
            new wave domain separated by 0.05 micron steps
        
        ynew : array 
            the dependent axis through cubic spline fit
        
        tck : tuple
            A tuple (t,c,k) containing the vector of knots, the B-spline coefficients, and the degree of the spline.
            
    See Also
    --------
    scipy.interpolate.splrep
    
    Notes
    -----
    This function is called internally by other functions and is not intended for the end user.
    
    '''
    
    low_range = wave_range[0]
    high_range = wave_range[1]
        
    # Removing the units
    x = spec.wave/u.micron
    y = spec.flux
    
    tck = interpolate.splrep(x, y, s=0.006) # spline prep, parameters to plot the interpolation 
    # Side note, not all spectrrum data encompasses 0.45 - 2.45 microns. Default k = 3 (cubic), as in source material
    xnew = np.arange(low_range, high_range, 0.05) # separating the wavelengths into 0.05 micron increments, from 0.85 - 2.45
    ynew = interpolate.splev(xnew, tck, der=0) # using the spline prep, obtain the "best fit" reflectance values
    
    return (xnew, ynew, tck)



def norm_splined(xx, yy, tck, norm = 1.2):
    '''
    Normalize the data to unity at some wavelength by division.
    
    Parameters
    ----------
    xx : array
        wave domain separated by 0.05 micron steps
    
    yy : array
        the dependent axis through cubic spline fit
        
    tck : tuple
        a tuple (t,c,k) containing the vector of knots, the B-spline coefficients, and the degree of the spline.
    
    norm : float, default = 1.2
        wavelength at which the spline fit is normalized
    
    Returns
    -------
    xxnew : array
        wave domain separated by 0.05 micron steps
    
    yynew : array
        the dependent axis through cubic spline fit, normalized at "norm" microns
            
    Notes
    -----
    This function is called internally by other functions and is not intended for the end user.
    
    '''
    yy_val = interpolate.splev(np.array([norm]), tck, der=0)
    yynew = yy/yy_val # now, 'norm' micron reflectance value = 1, does not affect slope of the spectrum
    xxnew = xx
    return xxnew, yynew



def lin_regress(xxnew, yynew):
    '''
    Using linear regression, slope of a spectrum is calculated.
    
    Parameters
    ----------
    xxnew : array
        wave domain separated by 0.05 micron steps
    
    yynew : array
        the dependent axis through cubic spline fit, normalized at "norm" microns
    
    Returns
    -------
    gamma : float
        the slope of the normalized, cubic spline-fitted data
            
    Notes
    -----
    This function is called internally by other functions and is not intended for the end user.
    
    '''
    nom = np.sum((xxnew - np.mean(xxnew))*(yynew - np.mean(yynew)))
    denom = np.sum((xxnew - np.mean(xxnew))**2)
    gamma = nom/denom
    return gamma



def form_data_pairs(obj, wave_range = [0.85, 2.45], norm = 1.2):
    '''
    Prepares data for principal component analysis (PCA).
    
    Parameters
    ----------
    obj : Spectrum object
         the spectra to be prepared for PCA
    
    wave_range : array of two floats, default = [0.85, 2.45]
        specifies wavelength range to consider [start, end).
        
    norm : float, default = 1.2
        wavelength at which the spline fit is normalized
    
    Returns
    -------
    tp_list : nested array
        contains two arrays, the 'wave' of the Spectra object and the 'flux' of some Spectrum object
        after cubic spline fit and normalizing, with the slope divided out.
    
    tp_slope : float
        slope of the normalized, cubic spline-fitted data obtained through linear regression
    
    tp_reg_line : array
        flux values of the regression line through the normalized, spline-fitted data where (norm, 1) is a point
    
    Notes
    -----
    This function is called internally by other functions and is not intended for the end user.
    
    '''
    
    tp_list = [] 
    tp_x = obj.wave/u.micron
    tp_y = obj.flux
    # obtain variables from cubic spline fitting
    tp_xx, tp_yy, tp_spline_vars = spline_fit(obj, wave_range)
    # normalize the data to unity at some wavelength define by 'norm'
    tp_xxnew, tp_yynew = norm_splined(tp_xx, tp_yy, tp_spline_vars, norm = norm)
    # calculate the slope of data using linear regression
    tp_slope = lin_regress(tp_xxnew, tp_yynew)

    # having normalized at 'norm' microns, here we translate the regression line to pass through (norm, 1)
    tp_reg_line = tp_xxnew*tp_slope + (1 - tp_slope*norm)
    # dividing out the slope using the line we just defined above
    tp_yynew_no_slope = tp_yynew/tp_reg_line
    # add to list, and return an array of wavelength-(slope removed)reflectance pairs
    tp_list.append((tp_xxnew, tp_yynew_no_slope))
    return np.array(tp_list), tp_slope, tp_reg_line



def MITHNEOS_to_Spectrum(txtPath):
    '''
    Read in a .txt file data from MITHNEOs to a Spectrum object used in classification. as a .txt file contains both V and NIR data, as well as numbers
    indicating what are the omitted rows. This function takes a string for the path of a file and
    returns a splat.Spectrum object constructed with the data in that file
    
    Parameters
    ----------
    txtPath : string
        path of the desired file 
    
    Returns
    -------
    sp : Spectrum object
        Containing only IR and NIR data in that file, where flux is unitless reflectance and wave is in microns
    
    Notes
    -----
    As a .txt file from MITHNEOS likely contains both V and NIR data, as well as numbers
    indicating what are the omitted rows in a specific syntax. This function parse the file for the IR and NIR portion.
    
    '''
    
    data = open(txtPath, 'r')
    
    abspath = os.path.abspath(txtPath)
    objName = os.path.split(abspath)[-1]

    dataString = data.read()

    # MITHNEOS data contains V band as well, NIR data is always separated by an extra new line after the V band
    ir_only_data = dataString.split('\n\n')

    if len(ir_only_data) > 1:
        ir_only_data = ir_only_data[1]
    else:
        ir_only_data = ir_only_data[0]

    dataList = ir_only_data.split("\n")

    # removing entries with 0 in "number of frames used to obtain data", noted as excluded values on MITHNEOS
    pruned_dataList = []
    for entry in dataList:
        if entry is not '': # ignores empty rows
            removed_spaces = entry.split()
            if (int(removed_spaces[-1]) != 0): 
                pruned_dataList.append(removed_spaces) # strings, of each value

    # pruned_dataList is a list of lists, of which the smaller list are rows of the .txt file
    # "Wavelength","reflectance", "sigma1", "# of frames"
    # casting strings to floats and constructing the columns we need, wave, flux, (and uncertainty)

    wavelengths = [] # wavelength
    flx = [] # reflectance
#     sig1 = [] # sigma 1 uncertainty 

    for row in pruned_dataList:
        wavelengths.append(float(row[0]))
        flx.append(float(row[1]))
#         sig1.append(float(row[2]))

    spwave = np.array(wavelengths)
    spflux = np.array(flx)
#     spone_sigma = np.array(sig1)

    sp = AsteroidSpectrum(wave=spwave,flux=spflux, flux_unit=u.dimensionless_unscaled)
    sp.name = objName
    
    return sp



def visnir_to_Spectrum(txtPath, instrument=None):
    '''
    Read in a .txt file data from vis-nir to a Spectrum object used in classification. as a .txt file contains both V and NIR data, as well as numbers
    indicating the omitted rows. This function takes a string for the path of a file and returns an AsteroidSpectrum instance constructed with the data in that file
    
    Parameters
    ----------
    txtPath : string
        path of the desired file 
    
    Returns
    -------
    sp : Spectrum object
        Containing only IR and NIR data in that file, where flux is unitless reflectance and wave is in microns
    
    Notes
    -----
    As a .txt file from MITHNEOS likely contains both V and NIR data, as well as numbers
    indicating what are the omitted rows in a specific syntax. This function parse the file for the IR and NIR portion.
    
    '''
    
    data = open(txtPath, 'r')

    full_path_name = os.path.abspath(txtPath)
    objName = os.path.split(full_path_name)[-1].split('.')[0]

    dataString = data.read()

    # MITHNEOS data contains V band as well, NIR data is always separated by an extra new line after the V band
    ir_only_data = dataString.split('\n \t \t \n')

    if len(ir_only_data) > 1 and len(ir_only_data[1]) != 0:
        ir_only_data = ir_only_data[1]
    else:
        ir_only_data = ir_only_data[0]

    dataList = ir_only_data.split("\n")
    
    wavelengths = [] # wavelength
    flx = [] # reflectance

    for row in dataList:
        tmp = row.split('\t')
        if(len(tmp) == 3):
            wavelengths.append(float(tmp[0]))
            flx.append(float(tmp[1]))
#         sig1.append(float(row[2]))

    spwave = np.array(wavelengths)
    spflux = np.array(flx)
#     spone_sigma = np.array(sig1)

    if instrument is not None:
        sp = AsteroidSpectrum(wave=spwave,flux=spflux, name=objName, instrument=instrument)
    else:
        sp = AsteroidSpectrum(wave=spwave,flux=spflux, name=objName)
    
    return sp



def calc_scores(pcs, subtracted_reflectance):
    '''
    Calculate PC scores of a spectrum given the PC eigenvectors, returns an array of scores. This function calculates
    scores based on given eigenvectors, and does not generate new eigprincipal component eigenvectors from a set of data.
    
    Parameters
    ----------
    pcs : array of arrays
        containing a number (5) of eigenvectors corresponding to the principal components
    
    subtracted_reflectance : array
        containing the reflectance of an asteroid at different wavelengths, with the mean from Demeo's supplemtary
        material subtracted from it.
    
    Returns
    -------
    score_list : array
        of principal component scores for a Spectrum according to the principal components in pcs
        
    See Also
    --------
    load_demeo_pcs
    
    '''
    score_list = []
    for pc in pcs:
        score_list.append(np.dot(pc.T, subtracted_reflectance.T))

    return np.array(score_list)

