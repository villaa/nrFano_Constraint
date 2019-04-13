import numpy as np
import matplotlib.pyplot as plt
import resfuncRead as rfr
import pandas as pd
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from data_check import * 
from bin_data_check import *

f = 1.5
fp = .75

k = 0.157 # from lindhard
q = 1.602*10**-19 #electron charge 
V = 4.0 # voltage Bias
eps = .0033 #keV

#for detector 1
p_alpha = 0.155393
p_beta = 9.60343*10**(-11)
p_gamma = 0.000506287
q_alpha = 0.166004
q_beta = 0.0023371
q_gamma = 9.52576*10**(-5)


def Yield_Er(Er_True,fano):
    Yield_er = []
    ER = []
    E_true =[]
    Eq_mean = []
    Ep_mean = []
    sigp_mean = []
    sigq_mean = []
 
    Er_True = np.sort(Er_True)
    
    for Eer in Er_True:

        #For mean number of electron hole pairs 
        N_mean = Eer/eps 

        E_p_mean = Eer+(V/1000)*N_mean
        E_q_mean = N_mean*eps
        sig_p_mean = np.sqrt(p_alpha + p_beta*E_p_mean + p_gamma*(E_p_mean**2) + (V/1000)**2*fano*N_mean)
        sig_q_mean = np.sqrt(q_alpha + q_beta*E_q_mean + q_gamma*(E_q_mean**2) + eps**2*fano*N_mean)
        

        #for fano varried number of electron hole pairs. 
        N_var = np.random.normal(0,np.sqrt(fano *N_mean)) + N_mean

        E_p_var = Eer + (V/1000)*N_var
        E_q_var = N_var*eps #with Fano
        sig_p_var = np.sqrt(p_alpha + p_beta*E_p_var + p_gamma*(E_p_var**2))  
        sig_q_var = np.sqrt(q_alpha + q_beta*E_q_var + q_gamma*(E_q_var**2))

        # find 'measured' Ep and Eq using N_var
        Ep_smear = np.random.normal(E_p_var,sig_p_var) #smearing 'measured' phonon energy with res that contains fano factor 
        Eq_smear = np.random.normal(E_q_var,sig_q_var) #smearing 'measured' charge energy with res that contains fano factor 
        

        Er = Ep_smear - (V/eps/1000)*Eq_smear # measured recoil energy 
     
        Yield = Eq_smear/Er # "measured" Ionization Yield
        
        #append stuff 
        E_true.append(Eer)
        Ep_mean.append(E_p_mean) # used for 'expected'
        Eq_mean.append(E_q_mean) #used for 'expected'
        sigp_mean.append(sig_p_mean) #used for 'expected'
        sigq_mean.append(sig_q_mean) #used for 'expected
        ER.append(Er)
        Yield_er.append(Yield)

    df = pd.DataFrame({'E_true':E_true,'E_measured':ER,'Yield':Yield_er,'Ep_mean':Ep_mean,'Eq_mean':Eq_mean,'sigp_mean':sigp_mean,'sigq_mean':sigq_mean})
    
    return df

