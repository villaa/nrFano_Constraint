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
q_beta = 0.00233716

q_gamma = 9.52576*10**(-5)


def Yield_NR(Er_true):
    
    import sys
    sys.path.insert(0,'/Users/Mitch 1/Desktop/nrFano_Constraint/python')
    import lindhard as lind
    lpar = lind.getLindhardPars('Ge',True)
    
    ER = []
    Yield = []
    sigp_nr = []
    sigq_nr = []
    EPnr = []
    EQnr = []
    Enr_true = []
    

    for Enr in Er_true:

        #Enr = np.random.choice(E1nr)
        Enr_true.append(Enr)
        

        a = 0.16
        b = 0.18 
        c = 0.04
        A = Enr**(1-b)
        B = Enr**(1-2*b)
        C = Enr**(1-3*b)


        #Yield Edelweiss 
        #Y = a*Enr**b

         #use the "calculated" value of k
        #Y = lind.getLindhard(lpar)
        Y = ynr_mu(Enr)


        #EdelWeiss Fano

        F = c**2/((eps*a/(A))+(2*V*a**2/(B*1000))+(2*(V/1000)**2*a**3/(eps*C)))

        #Constant Fano
        #F = 10

        #Neh = Y(Enr*1000)*Enr/eps #number of electron-hole pairs. 
        Neh = Y*Enr/eps #number of electron-hole pairs. 
        N_eh = Neh + np.random.normal(0,np.sqrt(F*Neh))

        Ptnr = Enr+(V/1000)*N_eh #central value of Pt
        Qnr = N_eh*eps
        
        Ptnr1 = Enr + (V/1000)*Neh # mean Neh, not varied. 
        Qnr1 = Neh*eps #mean Neh 


       #With Fano Factor. 
        sig_p = np.sqrt(p_alpha + p_beta*Ptnr1 + p_gamma*(Ptnr1**2))  #Phonon uncertainty (energy dependent)
        sig_q = np.sqrt(q_alpha + q_beta*Qnr1 + q_gamma*(Qnr1**2)) #Charge uncertainty 
        sigp_nr.append(sig_p)
        sigq_nr.append(sig_q)

        Fnr = np.random.normal(0.0,sig_p) #random sample assuming phonon variance 
        Fq = np.random.normal(0.0,sig_q) #random sampel assuming charge variance 


        Ptnr1 = Ptnr + Fnr 
        Qnr1 = Qnr + Fq 
        Ernr = Ptnr1 - (V/(eps*1000))*Qnr1 #Meaured ER
        yield1 = Qnr1 / Ernr  #Measured Yieldf

        EPnr.append(Ptnr1)
        EQnr.append(Qnr1)

        #Store Stuff 
        ER.append(Ernr)
       
        Yield.append(yield1)
        
    df = pd.DataFrame({'E_true':Enr_true,'E_measured':ER,'Yield':Yield,'EP':EPnr,'EQ':EQnr,'Sigp':sigp_nr,'Sigq':sigq_nr})
     
        
    return  df



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

