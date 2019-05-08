
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
p_alpha = 0.155393
p_beta = 9.60343*10**(-11)
p_gamma = 0.000506287
q_alpha = 0.166004
q_beta = 0.0023371
q_gamma = 9.52576*10**(-5)





def Yield_NR(Er_True):
    
    
 
    Enr = np.sort(Er_True)

    Y = ynr_mu(Enr)
    print(Y)

    a = 0.16
    b = 0.18 
    c = 0.04
    A = Enr**(1-b)
    B = Enr**(1-2*b)
    C = Enr**(1-3*b)

    fano = c**2/((eps*a/(A))+(2*V*a**2/(B*1000))+(2*(V/1000)**2*a**3/(eps*C))) #edelweiss fano


    N_mean = Y*Enr/eps

    sigN = np.sqrt(N_mean * fano) 

    E_p_mean = Enr+(V/1000)*N_mean
    E_q_mean = N_mean*eps

    '''These resolutions are for the independent(v1) pdf'''
    sig_p_v1= np.sqrt(p_alpha + p_beta*E_p_mean + p_gamma*(E_p_mean**2)+ (V/1000)**2*fano*N_mean)
    sig_q_v1 = np.sqrt(q_alpha + q_beta*E_q_mean + q_gamma*(E_q_mean**2)+ eps**2*fano*N_mean)

    '''These resolutions are for the dependent(v2) pdf'''
    sig_p_v2= np.sqrt(p_alpha + p_beta*E_p_mean + p_gamma*(E_p_mean**2))
    sig_q_v2 = np.sqrt(q_alpha + q_beta*E_q_mean + q_gamma*(E_q_mean**2))

 

    '''This section is for generating measured data'''    #for fano varried number of electron hole pairs. 
    N_var = np.random.normal(0,np.sqrt(fano *N_mean)) + N_mean
    #N_var =  N_mean

    E_p_var = Enr+ (V/1000)*N_var
    E_q_var = N_var*eps #with Fano


    sig_p_var = np.sqrt(p_alpha + p_beta*E_p_var + p_gamma*(E_p_var**2))  
    sig_q_var = np.sqrt(q_alpha + q_beta*E_q_var + q_gamma*(E_q_var**2))

    # find 'measured' Ep and Eq using N_var
    Ep_smear = np.random.normal(E_p_var,sig_p_var) #smearing 'measured' phonon energy with res that contains fano factor 
    Eq_smear = np.random.normal(E_q_var,sig_q_var) #smearing 'measured' charge energy with res that contains fano factor 


    Er = Ep_smear - (V/eps/1000)*Eq_smear # measured recoil energy 

    Yield = Eq_smear/Er # "measured" Ionization Yield

    Er_m = np.sort(Er)


    df = pd.DataFrame({'E_measured':Er_m,'E_true':Enr,'N_mean':N_mean,'sig_N':sigN,'Yield':Yield,'Ep_mean':E_p_mean,'Eq_mean':E_q_mean,'sigp_v1':sig_p_v1,'sigq_v1':sig_q_v1,'sigp_v2':sig_p_v2,'sigq_v2':sig_q_v2})
    df2 = pd.DataFrame({'Yield':Yield,'ER':Er,'EP':Ep_smear,'EQ':Eq_smear,'E_sort':Er_m,'EP_var': E_p_var,'EQ_var':E_q_var})

    return df,df2




