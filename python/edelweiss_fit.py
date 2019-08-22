import numpy as np
import lmfit as lmf

from EdwRes import *
from prob_dist import *

# lmfit needs a residuals function
# par_dict is a dictionary of the form
# {'V': 4.0, 'eps_eV': 3.0}
def residualNR(params, x, data, eps_data, par_dict):
    ion_center_0keV = params['ion_center_0keV']
    ion_guard_0keV = params['ion_guard_0keV']
    ion_122keV = params['ion_122keV']
    heat_0keV = params['heat_0keV']
    heat_122keV = params['heat_122keV']
    aH = params['aH']
    C = params['C']
    m = params['m']

    fit_func = get_sig_nuc_func_alt(ion_center_0keV, ion_guard_0keV, ion_122keV, heat_0keV, heat_122keV, par_dict, aH, C, m)
    
    model = fit_func(x)

    return (data-model) / eps_data

def residualER(params, x, data, eps_data, par_dict):
    ion_center_0keV = params['ion_center_0keV']
    ion_guard_0keV = params['ion_guard_0keV']
    ion_122keV = params['ion_122keV']
    heat_0keV = params['heat_0keV']
    heat_122keV = params['heat_122keV']
    aH = params['aH']

    fit_func = get_sig_gamma_func(ion_center_0keV, ion_guard_0keV, ion_122keV, heat_0keV, heat_122keV, par_dict, aH)
    
    model = fit_func(x)
    
    return (data-model) / eps_data

def edelweiss_fit(pars, ER_data, NR_data):
    # GGA3 parameters from Edelweiss tables
    ion_center_0keV = 1.3
    ion_guard_0keV = 1.5
    heat_0keV = 0.4
    ion_122keV = 3.1 
    heat_122keV = 2.7
    aH = 0.0157
    pars = {
        'V' : 4.0,
        'eps_eV' : 3.0
        }

    # first Edelweiss fits the ER band to get aH
    paramsER = lmf.Parameters()
    paramsER.add('ion_center_0keV', value=ion_center_0keV, vary=False)
    paramsER.add('ion_guard_0keV', value=ion_guard_0keV, vary=False)
    paramsER.add('ion_122keV', value=ion_122keV, vary=False)
    paramsER.add('heat_0keV', value=heat_0keV, vary=False)
    paramsER.add('heat_122keV', value=heat_122keV, vary=False)
    paramsER.add('aH', value=0.01638)

    ER_band_fit = edelweiss_fit_ER(paramsER, pars, \
                                   ER_data['Erecoil'], ER_data['sigma'], ER_data['sigma_err'], \
                                   residualER)

    # then Edelweiss fits the NR band to get C (and we get m)
    # aH is fixed to the value determined in the ER band fit
    paramsNR = lmf.Parameters()
    paramsNR.add('ion_center_0keV', value=ion_center_0keV, vary=False)
    paramsNR.add('ion_guard_0keV', value=ion_guard_0keV, vary=False)
    paramsNR.add('ion_122keV', value=ion_122keV, vary=False)
    paramsNR.add('heat_0keV', value=heat_0keV, vary=False)
    paramsNR.add('heat_122keV', value=heat_122keV, vary=False)
    paramsNR.add('aH', value=0.03808, vary=False)
    paramsNR.add('C', value=0.04)
    paramsNR.add('m', value=0)
    
    NR_band_fit = edelweiss_fit_NR()