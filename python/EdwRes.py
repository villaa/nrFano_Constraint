import numpy as np
from functools import partial

def get_heatRes(sig0, a, E_keV):
    """return the heat resolution at energy E_keV.  sig0, E_keV assumed to be in units of keV."""
    
    # see eqn (5) in 2004 NIMA Edelweiss paper
    sigH = np.sqrt(sig0**2 + (a*E_keV)**2)
    
    # multiply by 2.355 to get FWHM
    return sigH

def get_heatRes_func(FWHM0, FWHM122, aH=None):
    """returns a resolution function given the FWHM values at 0 keV and 122 keV"""
    
    # convert from FWHM to sigma
    # or maybe aH is calculated from the FWHM values?
    sig0 = FWHM0 #/2.355
    sig122 = FWHM122 #/2.355
    
    # calculate aH, which is unitless
    if aH is None:
        aH = np.sqrt((sig122**2 - sig0**2)/122**2)

    #print ("aH is: ", aH)
    
    # create function
    return partial(get_heatRes, sig0, aH)

def get_ionRes_func(FWHM_center, FWHM_guard, FWHM122):
    FWHM0 = np.sqrt(FWHM_center**2 + FWHM_guard**2)
    
    return get_heatRes_func(FWHM0, FWHM122)

def Q_avg(E_keV):
    return 0.16*np.power(E_keV,0.18)

def get_sig_gamma(sigI, sigH, V, E_keV):
    return ((1+V/3)/E_keV)*np.sqrt((sigI(E_keV)/2.355)**2 + (sigH(E_keV)/2.355)**2)

def get_sig_gamma_func(FWHM_center, FWHM_guard, FWHM122_ion, FWHM0_heat, FWHM122_heat, V, aH=None):
    # get the ionization resolution function
    sigI = get_ionRes_func(FWHM_center, FWHM_guard, FWHM122_ion)
    
    # get the heat resolution function
    sigH = get_heatRes_func(FWHM0_heat, FWHM122_heat, aH)
    
    return partial(get_sig_gamma, sigI, sigH, V)

def get_sig_nuc(sigI, sigH, V, E_keV):
    return ((1+V/3)/E_keV)*np.sqrt((sigI(E_keV)/2.355)**2 + (sigH(E_keV)/2.355)**2)

def get_sig_nuc_func(FWHM_center, FWHM_guard, FWHM122_ion, FWHM0_heat, FWHM122_heat, V, aH=None):
    # get the ionization resolution function
    sigI = get_ionRes_func(FWHM_center, FWHM_guard, FWHM122_ion)
    
    # get the heat resolution function
    sigH = get_heatRes_func(FWHM0_heat, FWHM122_heat, aH)
    
    return partial(get_sig_gamma, sigI, sigH, V)
