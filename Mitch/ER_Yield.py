
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
  
 
    Eer = np.sort(Er_True)
    
    N_mean = Eer/eps 
    sigN = np.sqrt(N_mean * fano) 

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



    df = pd.DataFrame({'EQ_smear':Eq_smear,'EP_smear':Ep_smear,'N_mean':N_mean,'sig_N':sigN,'E_true':Eer,'E_measured':Er,'Yield':Yield,'Ep_mean':E_p_mean,'Eq_mean':E_q_mean,'sigp_mean':sig_p_mean,'sigq_mean':sig_q_mean})
    
    return df
