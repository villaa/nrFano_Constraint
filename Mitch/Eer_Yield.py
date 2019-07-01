
import numpy as np
import matplotlib.pyplot as plt
import resfuncRead as rfr
import pandas as pd
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from band_check import * 
from bin_data_check import *


import sys
sys.path.insert(0,'/Users/Mitch 1/Desktop/nrFano_Constraint/python')
import lindhard as lind
lpar = lind.getLindhardPars('Ge',True)


f = 1.5
fp = .75

k = 0.157 # from lindhard
q = 1.602*10**-19 #electron charge 
V = 4.0 # voltage Bias
eps = .0033 #keV

#for detector 1
#p_alpha = 0.155393
p_alpha = 0.05
p_beta = 9.60343*10**(-11)
p_gamma = 0.000506287
#q_alpha = 0.166004
q_alpha = 0.1
q_beta = 0.0023371
q_gamma = 9.52576*10**(-5)





def Yield_Er(Er_True,fano):
    
    
 
    Eer = np.sort(Er_True)


    '''For generating data for version 1 of our simulation. Here, the fano factor is accounted for in the resolutions.'''

    N_mean = Eer/eps

    sigN = np.sqrt(N_mean * fano) 

    Ep_true_v1 = Eer+(V/1000)*N_mean
    Eq_true_v1 = N_mean*eps


    sig_p_v1= np.sqrt(p_alpha + p_beta*Ep_true_v1 + p_gamma*(Ep_true_v1**2)+ (V/1000)**2*sigN**2)
    sig_q_v1 = np.sqrt(q_alpha + q_beta*Eq_true_v1 + q_gamma*(Eq_true_v1**2)+ eps**2*sigN**2)

    Ep_measured_v1 = np.random.normal(Ep_true_v1,sig_p_v1)
    Eq_measured_v1 = np.random.normal(Eq_true_v1,sig_q_v1)

    Er_v1 = Ep_measured_v1 - (V/eps/1000)*Eq_measured_v1

    Yield_v1 = Eq_measured_v1/Er_v1

    Er_v1 = np.sort(Er_v1)


    df_v1 = pd.DataFrame({'E_true':Eer,'N_mean':N_mean,'Ep_true':Ep_true_v1,'Eq_true':Eq_true_v1,'sigp':sig_p_v1,'sigq':sig_q_v1,'Ep_measured':Ep_measured_v1,'Eq_measured':Eq_measured_v1,'E_measured':Er_v1,'Yield':Yield_v1})


 
    '''These resolutions are for the dependent exoected (v2) pdf there is no variation in the resolution. For the data, there will be. So these will not be used for generating (V2) data'''
    sig_p_expected= np.sqrt(p_alpha + p_beta*Ep_true_v1 + p_gamma*(Ep_true_v1**2))
    sig_q_expected = np.sqrt(q_alpha + q_beta*Eq_true_v1 + q_gamma*(Eq_true_v1**2))


    '''This section is for generating measured data v2 fano factor is accounted for in the number of electron hole pairs created for a single recoil energy'''    #for fano varried number of electron hole pairs. 

    N_var = np.random.normal(0,sigN) + N_mean
    #N_var =  N_mean

    Ep_true_v2 = Eer+ (V/1000)*N_var
    Eq_true_v2 = N_var*eps #with Fano


    sig_p_v2 = np.sqrt(p_alpha + p_beta*Ep_true_v2 + p_gamma*(Ep_true_v2**2))  
    sig_q_v2 = np.sqrt(q_alpha + q_beta*Eq_true_v2 + q_gamma*(Eq_true_v2**2))

    # find 'measured' Ep and Eq using N_var
    Ep_measured_v2 = np.random.normal(Ep_true_v2,sig_p_v2) #smearing 'measured' phonon energy with res that contains fano factor 
    Eq_measured_v2 = np.random.normal(Eq_true_v2,sig_q_v2) #smearing 'measured' charge energy with res that contains fano factor 


    Er_v2 = Ep_measured_v2 - (V/eps/1000)*Eq_measured_v2 # measured recoil energy 

    Yield_v2 = Eq_measured_v2/Er_v2 # "measured" Ionization Yield

    Er_v2 = np.sort(Er_v2)


    df_v2 = pd.DataFrame({'E_true':Eer,'N_mean':N_mean,'N_var':N_var, 'sig_N':sigN, 'Ep_true':Ep_true_v2,'Eq_true':Eq_true_v2,'sigp':sig_p_v2,'sigq':sig_q_v2,'Ep_measured':Ep_measured_v2,'Eq_measured':Eq_measured_v2,'E_measured':Er_v2,'Yield':Yield_v2,'sigp_expected':sig_p_expected,'sigq_expected':sig_q_expected})


    return df_v1,df_v2




