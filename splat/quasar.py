##########################################################################
###                        Import Statements                           ###
##########################################################################

#internal imports
import copy

# external imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from tqdm import tqdm
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_generic_continuum

# splat imports
from splat.initialize import *
from splat.core import Spectrum, compareSpectra
from splat.plot import plotSpectrum


##########################################################################
###                      Redshift Functions                            ###
##########################################################################

def find_best_redshift(targ_sp, z_limits, steps, num_sims, template = 'Glikman_2006', bounds = False, MC = False, statistics = False, plotting = False):
    '''
    PURPOSE
    --------------------------------------------------------------
    Runs Monte-Carlo (MC) simulations of find_redshift() to determine the best redshift with error

    PARAMETERS
    --------------------------------------------------------------
    targ_sp: SPLAT Spectrum object
        The target quasar
    z_limits: tuple of floats
        The starting and ending redshift estimate
    steps: int
        The number of steps between the starting and ending redshift estimates
    num_sims: int
        The number of MC simulations to run
    template: string, default = 'Glikman_2006'
        The template used in cross-correlation
    bounds: boolean, default = False
        Whether to plot the emission features and correlation bounds of the target spectrum
    MC: boolean, default = False
        Whether to plot the original target spectrum with the MC spectra
    statistics: boolean, default = False
        Whether to plot the correlation coefficient/chi-squared against redshift estimates for the first MC simulation
    plotting: boolean, default = False
        Whether to plot the target with the template the best redshift
    
    OUTPUTS
    --------------------------------------------------------------
    z_mean: float
        The best redshift estimate given by the MC simulations
    z_std: float
        The associated error on z_mean given by the MC simulations

    EXAMPLE
    --------------------------------------------------------------
    >>> z_mean, z_std = run_mc_sims(targ_sp, (0, 0.5), 100, 5)
    The best redshift estimate is 0.21140939597315436 +/- 0.0

    DEPENDENCIES
    --------------------------------------------------------------
    numpy
    copy
    splat
    astropy
    specutils
    matplotlib
    '''

    #find target emission features and correlation bounds to calculate over; plot if desired
    corr_targ_sp = correct_targ(targ_sp)
    sub_targ_sp = subtract_continuum(corr_targ_sp)
    ft_pixs = find_features(sub_targ_sp)
    corr_bounds = find_corr_bounds(sub_targ_sp, ft_pixs)
    if bounds:
        plot_corr_bounds(corr_targ_sp, sub_targ_sp, ft_pixs, corr_bounds)
    
    #generate MC spectra; plot if desired
    mc_sims = generate_mc(corr_targ_sp, num_sims)
    if MC:
        plot_mc(corr_targ_sp, mc_sims)

    #read in the template
    if template == 'Glikman_2006':
        data_path = SPLAT_PATH + '/resources/SpectralTemplates/GlikmanQuasars'
        temp_df = pd.read_csv(data_path + '/Glikman_2006.csv', names = ['wave', 'flux', 'noise'])
        temp_sp = Spectrum(wave = list(temp_df['wave']), flux = list(temp_df['flux']), noise = list(temp_df['noise']), name = 'template')
    
    #find redshift for each MC simulation
    z_list = []
    for count, mc in enumerate(tqdm(mc_sims)):
        best_z, z_arr, corr_arr, chi_arr = find_redshift(mc, temp_sp, z_limits, steps, corr_bounds)
        z_list.append(best_z)
        #plot statistics if desired
        if count == 0 and statistics:
            plot_stats(z_arr, corr_arr, chi_arr)

    #calculate the best redshift with error
    z_mean = np.mean(z_list)
    z_std = np.std(z_list)
    print('The best redshift estimate is', z_mean, '+/-', z_std)

    #plot results if desired
    if plotting:
        plot_results(corr_targ_sp, temp_sp, z_mean)
    
    return z_mean, z_std


def find_redshift(mc, temp_sp, z_limits, steps, corr_bounds):
    '''
    PURPOSE
    --------------------------------------------------------------
    Estimates the redshift of a target quasar by cross-correlating with a quasar template

    PARAMETERS
    --------------------------------------------------------------
    mc: SPLAT Spectrum object
        An MC simulation of the target quasar
    temp_sp: SPLAT Spectrum object
        The template quasar
    z_limits: tuple of floats
        The starting and ending redshift estimates
    steps: int
        The number of steps between the starting and ending redshift estimates
    corr_bounds: list
        List of pixel bounds over which the cross correlation is calculated
    
    OUTPUTS
    --------------------------------------------------------------
    best_z: float
        The best redshift estimate
    z_arr: list
        List of all tested redshifts
    corr_arr: list
        List of correlation coefficients corresponding to each redshift estimate
    chi_arr: list
        List of chi-squared corresponding to each redshift estimate

    EXAMPLE
    --------------------------------------------------------------
    >>> best_z, z_arr, corr_arr, chi_arr = find_redshift(mc, temp_sp, (0, 0.5), 150, corr_bounds)
    >>> print(best_z)
    0.21476510067114093

    DEPENDENCIES
    --------------------------------------------------------------
    numpy
    copy
    splat
    '''  
        
    #prepare to find the best redshift
    z_arr = np.linspace(z_limits[0], z_limits[1], steps)
    corr_arr = []
    chi_arr = []
    
    #find the best redshift
    for z in z_arr:
        int_targ_sp, int_temp_sp = interpolate_spectra(mc, temp_sp, z)
        int_corr_bounds = revise_corr_bounds(mc, int_targ_sp, corr_bounds)
        corr, chi = calc_corr(int_targ_sp, int_temp_sp, int_corr_bounds)
        corr_arr.append(corr)
        chi_arr.append(chi)
    best_z = z_arr[np.argmax(corr_arr)]
    
    return best_z, z_arr, corr_arr, chi_arr



##########################################################################
###                         Core Functions                             ###
##########################################################################

def correct_targ(targ_sp):
    '''
    PURPOSE
    --------------------------------------------------------------
    Remove telluric regions and smooth target spectrum

    PARAMETERS
    --------------------------------------------------------------
    targ_sp: SPLAT Spectrum object
        The target quasar

    OUTPUTS
    --------------------------------------------------------------
    corr_targ_sp: SPLAT Spectrum object
        The corrected target quasar

    EXAMPLE
    --------------------------------------------------------------
    >>> corr_targ_sp = apply_corrections(targ_sp)
    >>> print(targ_sp.flux)
    [9.81469e-16 9.30530e-16 8.51516e-16 ... 1.19232e-14 1.08384e-13
        2.62537e-14] erg / (micron s cm2)
    >>> print(corr_targ_sp.flux)
    [1.11676769e-16 1.20395136e-16 1.29406886e-16 ... 4.87359296e-17
        4.52625760e-17 4.11283324e-17] erg / (micron s cm2)

    DEPENDENCIES
    --------------------------------------------------------------
    copy
    splat
    '''
    
    #create duplicate to apply corrections to
    targc_sp = copy.deepcopy(targ_sp)
    
    #remove telluric absorption regions
    targc_sp.remove([0, 1.1])
    targc_sp.remove([1.35, 1.4])
    # targc_sp.remove([1.8, 1.93])
    targc_sp.remove([2.56, 2.86])
    targc_sp.remove([3.31, 3.32])
    targc_sp.remove([4.1, 4.5])

    #smooth spectra to emphasize emission features
    targc_sp.smooth(100)
    
    return targc_sp


def subtract_continuum(corr_targ_sp):
    '''
    PURPOSE
    --------------------------------------------------------------
    Subtract the continuum from the target spectrum to more easily locate emission features

    PARAMETERS
    --------------------------------------------------------------
    corr_targ_sp: SPLAT Spectrum object
        The corrected target quasar
    
    OUTPUTS
    --------------------------------------------------------------
    sub_targ_sp: SPLAT Spectrum object
        The continuum-subtracted target quasar

    EXAMPLE
    --------------------------------------------------------------
    >>> sub_targ_sp = subtract_continuum(corr_targ_sp)
    >>> print(corr_targ_sp.flux)
    [1.11676769e-16 1.20395136e-16 1.29406886e-16 ... 4.87359296e-17
        4.52625760e-17 4.11283324e-17] erg / (micron s cm2)
    >>> print(sub_targ_sp.flux)
    [-1.12708952e-16 -1.03971113e-16 -9.49398963e-17 ... -5.44362540e-17
        -5.78565235e-17 -6.19376472e-17] erg / (micron s cm2)
    
    DEPENDENCIES
    --------------------------------------------------------------
    numpy
    specutils (Spectrum1D)
    astropy (units)
    splat (Spectrum)
    '''
    
    #convert target to specutil spectrum
    wave = np.array(corr_targ_sp.wave)
    flux = np.array(corr_targ_sp.flux)
    targ_util = Spectrum1D(spectral_axis = wave * u.um, flux = flux * u.Jy)
    
    #continuum-fit the target and subtract it from the spectrum
    fit = fit_generic_continuum(targ_util)
    continuum = fit(wave * u.um)
    sub_flux = flux - np.array(continuum)
    sub_targ_sp = Spectrum(wave, sub_flux)
    
    return sub_targ_sp


def find_features(sub_targ_sp):
    '''
    PURPOSE
    --------------------------------------------------------------
    Determine the central pixel of the emission features in the target spectrum

    PARAMETERS
    --------------------------------------------------------------
    sub_targ_sp: SPLAT Spectrum object
        The continuum-subtracted target quasar
    
    OUTPUTS
    --------------------------------------------------------------
    ft_pixs: list
        List of the central pixel of the emission features in the target spectrum

    EXAMPLE
    --------------------------------------------------------------
    >>> ft_pixs = find_features(sub_targ_sp)
    >>> print(ft_pixs)
    [684, 1911, 2416, 2985, 3279, 3529, 4598]

    DEPENDENCIES
    --------------------------------------------------------------
    numpy
    '''
    
    #get the wavelengths and flux 
    wave = np.array(sub_targ_sp.wave)
    flux = np.array(sub_targ_sp.flux)
    
    #prepare to find features
    ft_pixs = []
    bright_thresh = 0   #minimum intensity to qualify as emission feature
    bright_pixs = np.where(flux > bright_thresh)[0]   #candidate pixels for emission features
    search_size = 0.1   #one-sided wavelength window to search for emission features
    bright_count = 80   #minimum number of neighboring bright pixels to count as emission feature
    conv = len(wave)/(wave[-1] - wave[0])   #pixel/micron conversion
    
    #locate emission features
    for pix in bright_pixs:
        cand_inten = flux[pix]   #intensity of candidate pixel
        l_bound = int(pix - (search_size*conv))
        r_bound = int(pix + (search_size*conv))
        neighbor_intens = flux[l_bound:r_bound]   #intensity of neighboring pixels

        #if the candidate lies on edges, skip it
        if pix in range(0, 100) or pix in range(len(wave)-100, len(wave)):
            continue
        
        #if the candidate does not have enough good neighbors, skip it
        num_good_neighbors = len(neighbor_intens[neighbor_intens > bright_thresh])   #counts the number of neighbors above the threshold
        if num_good_neighbors < bright_count:
            continue
        
        #if reached, the candidate pixel is a feature. isolate the center by finding the local maximum
        max_neighbor_inten = np.max(neighbor_intens)   #get the highest neighbor intensity
        if np.all(cand_inten >= max_neighbor_inten):   #if the pixel intensity is higher than the neighbors, that pixel is the center
            ft_pixs.append(pix)
    
    return ft_pixs


def find_corr_bounds(targ_sp, ft_pixs):
    '''
    PURPOSE
    --------------------------------------------------------------
    Determine the correlation bounds around each emission feature in the target spectrum

    PARAMETERS
    --------------------------------------------------------------
    targ_sp: SPLAT Spectrum object
        The target quasar
    ft_pixs: list
        List of the central pixel of the emission features in the target spectrum
    
    OUTPUTS
    --------------------------------------------------------------
    corr_bounds: list
        List of pixel bounds over which the cross-correlation is calculated

    EXAMPLE
    --------------------------------------------------------------
    >>> corr_bounds = find_corr_bounds(targ_sp, ft_pixs)
    >>> print(corr_bounds)
    [(581, 786), (1808, 2013), (2313, 2518), (2882, 3087), (3176, 3381), (3426, 3631), (4495, 4700)]
    
    DEPENDENCIES
    --------------------------------------------------------------
    numpy
    '''
    
    #prepare to find correlation bounds
    corr_bounds = []
    wave = np.array(targ_sp.wave)
    
    #establish padding around central pixel to set correlation bounds
    corr_pad_wave = 0.05   #one-sided wavelength padding to set correlation bounds
    conv = len(wave)/(wave[-1] - wave[0])   #pixel/micron conversion
    corr_pad_pix = corr_pad_wave * conv
    
    #find the correlation bounds in pixels
    for pix in ft_pixs:
        left = int(pix - corr_pad_pix)
        right = int(pix + corr_pad_pix)
        corr_bounds.append((left, right))
        
    return corr_bounds


def generate_mc(corr_targ_sp, num_sims):
    '''
    PURPOSE
    --------------------------------------------------------------
    Generate MC target spectra

    PARAMETERS
    --------------------------------------------------------------
    corr_targ_sp: SPLAT Spectrum object
        The corrected target quasar
    num_sims: int
        The number of MC simulations to run
    
    OUTPUTS
    --------------------------------------------------------------
    mc_sims: list
        List of SPLAT Spectrums representing the MC simulations

    EXAMPLE
    --------------------------------------------------------------
    >>> mc_sims = generate_mc(corr_targ_sp, 5)
    >>> print(mc_sims)
    [ spectrum of MC 1,  spectrum of MC 2,  spectrum of MC 3,  spectrum of MC 4,  spectrum of MC 5]
    
    DEPENDENCIES
    --------------------------------------------------------------
    copy
    numpy (random)
    '''

    #prepare to generate MC simulations
    mc_sims = []

    #generate MC simulations
    for i in list(range(num_sims)):
        targc_sp = copy.deepcopy(corr_targ_sp)
        noise = np.random.normal(targc_sp.flux.value, targc_sp.noise.value) * targc_sp.flux_unit
        new_flux = targc_sp.flux + noise
        new_spec = Spectrum(wave = targc_sp.wave, flux = new_flux, noise = targc_sp.noise, name = 'MC {}'.format(i+1))
        mc_sims.append(new_spec)

    return mc_sims


def calc_corr(int_targ_sp, int_temp_sp, int_corr_bounds):
    '''
    PURPOSE
    --------------------------------------------------------------
    Calculate the correlation coefficient between the target and template over the correlation bounds

    PARAMETERS
    --------------------------------------------------------------
    int_targ_sp: SPLAT Spectrum object
        The interpolated target quasar
    int_temp_sp: SPLAT Spectrum object
        The interpolated template quasar
    int_corr_bounds: list
        Revised list of pixel bounds over which the cross correlation is calculated, accounting for interpolation
    
    OUTPUTS
    --------------------------------------------------------------
    corr: float
        The correlation coefficient between the target and template
    chi: float
        The chi-squared between the target and template

    EXAMPLE
    --------------------------------------------------------------
    >>> corr, chi = calc_corr(int_targ_sp, int_temp_sp, int_corr_bounds)
    >>> print(corr, chi)
    (0.9818672176993626, 800718.8997413844)

    DEPENDENCIES
    --------------------------------------------------------------
    numpy
    splat
    '''
    
    #prepare to calculate correlation
    wave = np.array(int_targ_sp.wave)
    targ_flux = np.array(int_targ_sp.flux)
    temp_flux = np.array(int_temp_sp.flux)
    
    #generate mask array for pixels not within the correlation bounds
    mask_pixs = np.ones(len(wave), dtype = bool)
    for bound_num in list(range(len(int_corr_bounds))):
        bound_left = int_corr_bounds[bound_num][0]
        bound_right = int_corr_bounds[bound_num][1]
        for pix in list(range(len(wave))):
            if (pix >= bound_left) and (pix <= bound_right):
                mask_pixs[pix] = 0   
                
    #mask out the pixels not within the correlation bounds and calculate correlation
    targ_flux[mask_pixs] = 0
    temp_flux[mask_pixs] = 0
    corr = np.corrcoef(targ_flux, temp_flux)[0][1]
    
    #calculate the chi-squared
    chi = compareSpectra(int_targ_sp, int_temp_sp)[0].value
    
    return corr, chi



##########################################################################
###                       Interpolation Function                       ###
##########################################################################

def interpolate_spectra(mc, temp_sp, z):
    '''
    PURPOSE
    --------------------------------------------------------------
    Redshift the template spectrum and interpolate it onto the target spectrum

    PARAMETERS
    --------------------------------------------------------------
    mc: SPLAT Spectrum object
        An MC simulation of the target quasar
    temp_sp: SPLAT Spectrum object
        The template quasar
    z: float
        The redshift estimate
    
    OUTPUTS
    --------------------------------------------------------------
    int_targ_sp: SPLAT Spectrum object
        The interpolated target quasar SPLAT spectrum
    int_temp_sp:SPLAT Spectrum object
        The interpolated template quasar SPLAT spectrum with the same wavelength range as the target

    EXAMPLE
    --------------------------------------------------------------
    >>> int_targ_sp, int_temp_sp = interpolate_spectra(targ_sp, temp_sp, 0.5)
    >>> print(int_targ_sp.wave[0], int_targ_sp.wave[-1])
    1.10046 micron 4.09902 micron
    >>> print(int_temp_sp.wave[0], int_temp_sp.wave[-1])
    1.10046 micron 4.09902 micron

    DEPENDENCIES
    --------------------------------------------------------------
    copy
    splat
    '''
    
    #prepare to interpolate
    targc_sp = copy.deepcopy(mc)
    tempc_sp = copy.deepcopy(temp_sp)
    tempc_sp.redshift(z)

    #interpolate the template onto the target to match the wavelength ranges
    wave_min, wave_max = give_new_wave(targc_sp, tempc_sp)
    tempc_sp.trim([wave_min-0.01, wave_max+0.01])
    targc_sp.trim([wave_min, wave_max])
    tempc_sp.toWavelengths(targc_sp.wave)
        
    return targc_sp, tempc_sp


def give_new_wave(targ_sp, temp_sp):
    '''
    PURPOSE
    --------------------------------------------------------------
    Give the interpolated wavelength range endpoints

    PARAMETERS
    --------------------------------------------------------------
    targ_sp: SPLAT Spectrum object
        The target quasar
    temp_sp: SPLAT Spectrum object
        The template quasar
    
    OUTPUTS
    --------------------------------------------------------------
    max_lows: float
        The low end of the interpolated wavelength range
    min_highs: float
        The high end of the interpolated wavelength range

    EXAMPLE
    --------------------------------------------------------------
    >>> give_new_wave(targ_sp, temp_sp)
    (0.80777, 3.50954)

    DEPENDENCIES
    --------------------------------------------------------------
    None
    '''
    
    targ_low = targ_sp.wave[0].value
    temp_low = temp_sp.wave[0].value
    max_lows = max(targ_low, temp_low)
    
    targ_high = targ_sp.wave[-1].value
    temp_high = temp_sp.wave[-1].value
    min_highs = min(targ_high, temp_high)
    
    return max_lows, min_highs


def revise_corr_bounds(mc, int_targ_sp, corr_bounds):
    '''
    PURPOSE
    --------------------------------------------------------------
    Adjust the correlation bounds to account for the new interpolated wavelength range

    PARAMETERS
    --------------------------------------------------------------
    mc: SPLAT Spectrum object
        An MC simulation of the target quasar
    int_targ_sp: SPLAT Spectrum object
        The interpolated target quasar
    corr_bounds: list
        List of pixel bounds over which the cross correlation is calculated
    
    OUTPUTS
    --------------------------------------------------------------
    int_corr_bounds: list
        Revised list of pixel bounds over which the cross correlation is calculated, accounting for interpolation

    EXAMPLE
    --------------------------------------------------------------
    >>> int_corr_bounds = revise_corr_bounds(mc, int_targ_sp, corr_bounds)
    >>> print(int_corr_bounds)
    [(581, 786), (1808, 2013), (2313, 2518), (2882, 3087), (3176, 3381), (3426, 3631), (4495, 4700)]

    DEPENDENCIES
    --------------------------------------------------------------
    numpy
    '''
    
    #change correlation bounds from pixels to wavelength
    wave = np.array(mc.wave)
    corr_bounds_wave = [(wave[left_pix], wave[right_pix]) for (left_pix, right_pix) in corr_bounds]
    
    #get the wavelength limits of the interpolated target
    wave_min = int_targ_sp.wave.value[0]
    wave_max = int_targ_sp.wave.value[-1]
    
    #determine indicies of all correlation bounds that are out of range after interpolation
    remove_inds = []
    for bound_num in list(range(len(corr_bounds_wave))):
        left_wave = corr_bounds_wave[bound_num][0]
        right_wave = corr_bounds_wave[bound_num][1]
        if not(left_wave > wave_min and right_wave > wave_min and left_wave < wave_max and right_wave < wave_max):
            remove_inds.append(bound_num)
    
    #remove out of range correlation bounds
    for ind in sorted(remove_inds, reverse = True):
        del corr_bounds[ind]
        
    return corr_bounds



##########################################################################
###                          Plotting Functions                        ###
##########################################################################


def plot_corr_bounds(corr_targ_sp, sub_targ_sp, ft_pixs, corr_bounds):
    '''
    PURPOSE
    --------------------------------------------------------------
    Plot the original target, continuum-subtracted target, emission feature centers, and correlation bounds

    PARAMETERS
    --------------------------------------------------------------
    corr_targ_sp: SPLAT Spectrum object
        The corrected target quasar
    sub_targ_sp: SPLAT Spectrum object
        The continuum-subtracted quasar
    ft_pixs: list
        List of the central pixel of the emission features in the target spectrum
    corr_bounds: list
        List of pixel bounds over which the cross-correlation is calculated

    OUTPUTS
    --------------------------------------------------------------
    None (matploblib figure)

    EXAMPLE
    --------------------------------------------------------------
    >>> plot_corr_bounds(corr_targ_sp, sub_targ_sp, ft_pixs, corr_bounds)

    DEPENDENCIES
    --------------------------------------------------------------
    None
    '''
    
    #plot the original target and continuum-subtracted target
    plt.figure()
    plt.title('Target with Correlation Bounds')
    plt.plot(corr_targ_sp.wave, corr_targ_sp.flux, label = 'Original')
    plt.plot(sub_targ_sp.wave, sub_targ_sp.flux, label = 'Continuum-Subtracted')
    plt.plot(sub_targ_sp.wave, sub_targ_sp.wave.value*0, 'k--')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    
    #plot the emission features
    wave = np.array(corr_targ_sp.wave)
    y_min, y_max = plt.gca().get_ylim()
    plt.vlines(wave[ft_pixs], y_min, y_max, color = 'red', label = 'Features')
    
    #plot the correlation bounds
    plt.vlines([wave[left] for (left, right) in corr_bounds], y_min, y_max, color = 'green', label = 'Correlation Bounds')
    plt.vlines([wave[right] for (left, right) in corr_bounds], y_min, y_max, color = 'green')
    
    plt.legend()
    plt.show()


def plot_mc(corr_targ_sp, mc_sims):
    '''
    PURPOSE
    --------------------------------------------------------------
    Plot the original target with the MC spectra

    PARAMETERS
    --------------------------------------------------------------
    corr_targ_sp: SPLAT Spectrum object
        The corrected target quasar
    mc_sims: list
        List of SPLAT Spectrums representing the MC simulations

    OUTPUTS
    --------------------------------------------------------------
    None (matploblib figure)

    EXAMPLE
    --------------------------------------------------------------
    >>> plot_mc(corr_targ_sp, mc_sims)

    DEPENDENCIES
    --------------------------------------------------------------
    None
    '''

    plt.figure()
    plt.title('MC Simulations')
    plt.plot(corr_targ_sp.wave, corr_targ_sp.flux, label = 'original')
    for count, mc in enumerate(mc_sims):
        plt.plot(mc.wave, mc.flux, label = 'MC %i' %(count+1))
    plt.plot(corr_targ_sp.wave, corr_targ_sp.wave.value*0, 'k--')
    plt.legend()
    plt.show()


def plot_results(targ_sp, temp_sp, best_z):
    '''
    PURPOSE
    --------------------------------------------------------------
    Plot the target spectrum and the template spectrum at the found redshift

    PARAMETERS
    --------------------------------------------------------------
    targ_sp: SPLAT Spectrum object
        The original target quasar
    temp_sp: SPLAT Spectrum object
        The template quasar
    best_z: float
        The estimated redshift
    
    OUTPUTS
    --------------------------------------------------------------
    None (matplotlib figure)

    EXAMPLE
    --------------------------------------------------------------
    >>> plot_results(targ_sp, temp_sp, best_z)

    DEPENDENCIES
    --------------------------------------------------------------
    matplotlib
    copy
    '''

    #apply redshift and make the test and target wavelength ranges match
    est_targ_sp, est_temp_sp,  = interpolate_spectra(targ_sp, temp_sp, best_z)
    est_temp_sp.name = 'template at estimated z: {}'.format(best_z)
    sp_list = [est_targ_sp, est_temp_sp]
    
    #redshift the template to the true z
    # true_targ_sp, true_temp_sp,  = interpolate_spectra(targ_sp, temp_sp, targ_dict[targ_sp.name]['true z'])
    # true_temp_sp.name = 'template at true z: {}'.format(targ_dict[targ_sp.name]['true z'])
    # sp_list.append(true_temp_sp)
    
    #plot all spectra
    plotSpectrum(sp_list, multiplot = True, layout = [4, 1])


def plot_stats(z_arr, corr_arr, chi_arr):
    '''
    PURPOSE
    --------------------------------------------------------------
    Plot the correlation coefficient/chi-squared against redshift estimates for the first MC simulation

    PARAMETERS
    --------------------------------------------------------------
    z_arr: list
        List of all tested redshifts
    corr_arr: list
        List of correlation coefficients corresponding to each redshift estimate
    chi_arr: list
        List of chi-squared corresponding to each redshift estimate
    
    OUTPUTS
    --------------------------------------------------------------
    None (matplotlib figure)

    EXAMPLE
    --------------------------------------------------------------
    >>> plot_stats(z_arr, corr_arr, chi_arr)

    DEPENDENCIES
    --------------------------------------------------------------
    matplotlib
    '''
    
    #plot the correlation coefficient vs. redshift
    plt.figure(figsize = (10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(z_arr, corr_arr)
    plt.title('Correlation Coefficient vs. Redshift')
    plt.xlabel('Redshift')
    plt.ylabel('Correlation Coefficient')

    #plot the chi-squared vs. redshift 
    plt.subplot(1, 2, 2)
    plt.plot(z_arr, chi_arr)
    plt.title('Chi-Squared vs. Redshift')
    plt.xlabel('Redshift')
    plt.ylabel('Chi-Squared')
    
    plt.tight_layout(pad = 4)
    plt.show()