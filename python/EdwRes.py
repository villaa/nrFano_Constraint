import numpy as np
from functools import partial
FWHM_to_SIG = 1 / (2*np.sqrt(2*np.log(2)))

def get_heatRes(sig0, a, E_keV):
    """return the heat resolution (1 sigma) at energy E_keV.  sig0, E_keV assumed to be in units of keV."""
    
    # see eqn (5) in 2004 NIMA Edelweiss paper
    sigH = np.sqrt(sig0**2 + (a*E_keV)**2)
    
    # multiply by 2.355 to get FWHM
    return sigH

def get_heatRes_func(FWHM0, FWHM122, aH=None):
    """returns a resolution function given the FWHM values at 0 keV and 122 keV"""
    
    # convert from FWHM to sigma
    # note that aH is calculated from the FWHM values
    # in the Edelweiss paper, so the aH values they report
    # are 2.355 times larger than values of a calculated using eq 5
    sig0 = FWHM0 * FWHM_to_SIG
    sig122 = FWHM122 * FWHM_to_SIG
    
    # calculate aH, which is unitless
    if aH is None:
        aH = np.sqrt((sig122**2 - sig0**2)/122**2)

    print ("aH is: ", aH)
    
    # create function
    return partial(get_heatRes, sig0, aH)

def get_ionRes_func(FWHM_center, FWHM_guard, FWHM122):
    FWHM0 = np.sqrt(FWHM_center**2 + FWHM_guard**2)
    
    return get_heatRes_func(FWHM0, FWHM122)

def Q_avg(E_keV):
    return 0.16*np.power(E_keV,0.18)

def get_sig_gamma(sigI, sigH, V, E_keV):
    return ((1+V/3)/E_keV)*np.sqrt((sigI(E_keV))**2 + (sigH(E_keV))**2)

def get_sig_neutron(sigI, sigH, V, Er_keV):
    E_keVee_I = np.multiply(Q_avg(Er_keV), Er_keV)
    E_keVee_H = np.multiply((1+(V/3.0)*Q_avg(Er_keV))/(1+(V/3.0)), Er_keV)
    # we're pretty sure Edelweiss uses the correct (above) conversion
    # and not the incorrect (below) conversion
    #E_keVee_H = np.multiply(Q_avg(Er_keV), Er_keV)

    a = np.multiply(1+(V/3)*Q_avg(Er_keV), sigI(E_keVee_I))
    b = np.multiply((1+V/3)*Q_avg(Er_keV), sigH(E_keVee_H))

    return (1/Er_keV)*np.sqrt(a**2 + b**2)

def get_sig_gamma_func(FWHM_center, FWHM_guard, FWHM122_ion, FWHM0_heat, FWHM122_heat, V, aH=None):
    # get the ionization resolution function
    sigI = get_ionRes_func(FWHM_center, FWHM_guard, FWHM122_ion)
    
    # get the heat resolution function
    sigH = get_heatRes_func(FWHM0_heat, FWHM122_heat, aH)
    
    return partial(get_sig_gamma, sigI, sigH, V)


def get_sig_nuc_func(FWHM_center, FWHM_guard, FWHM122_ion, FWHM0_heat, FWHM122_heat, V, aH=None):
    # get the ionization resolution function
    sigI = get_ionRes_func(FWHM_center, FWHM_guard, FWHM122_ion)
    
    # get the heat resolution function
    sigH = get_heatRes_func(FWHM0_heat, FWHM122_heat, aH)
    
    return partial(get_sig_neutron, sigI, sigH, V)

def get_sig_nuc_func_fit(FWHM_center, FWHM_guard, FWHM122_ion, FWHM0_heat, FWHM122_heat, V, aH=None, C = None):
    def fit_func(E_keVee):
        sig_nuc_func = get_sig_nuc_func(FWHM_center, FWHM_guard, FWHM122_ion, FWHM0_heat, FWHM122_heat, V, aH)

        return np.sqrt(np.pow(sig_nuc_func(E_keVee),2) + np.pow(C,2))

    return fit_func
