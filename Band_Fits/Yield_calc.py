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
    EQ = []
    EQ_true = []
    EP = []
    EP_true = []
    sig_p = []
    sig_q = []
    sig_p0 = []
    Er_fano = fano 


    QER = []
 
    Er_True = np.sort(Er_True)
    
    
    for Eer in Er_True:


        E_true.append(Eer) #sorted array of 'true' recoil energies. 

 
        #calculate number of electron hole pairs. 
    
        N_e = Eer/eps #mean number of electron hole pairs producde. 
        N_er = np.random.normal(0,np.sqrt(Er_fano *N_e)) + N_e  # 0.13 is the fano factor for germainium 
        
        #calculate Phonon energy 
        Pter = Eer+(V/1000)*N_er #phonon energy, with fano 
        
        Pter_noF = Eer+(V/1000)*N_e # phonon energy, no fano 
        
        #calculate Charge energy
        
        Qer1 = N_e*eps #no Fano
        Qer1 = N_er*eps #with Fano
        
        
        #resolutions
        
        sig_p1 = np.sqrt(p_alpha + p_beta*Pter + p_gamma*(Pter**2)) #phonon res, with fano 
        sig_p.append(sig_p1)
        
        sig_p_noF = np.sqrt(p_alpha + p_beta*Pter_noF + p_gamma*(Pter_noF**2)) #phonon res, no fano
        sig_p0.append(sig_p_noF)
    
        sig_q1 = np.sqrt(q_alpha + q_beta*Eer + q_gamma*(Eer**2)) #Charge uncertainty
        sig_q.append(sig_q1)
        


        
        # Find 'Measured' phonon energy. 

        Pter1 = np.random.normal(Pter,sig_p1) #smearing 'measured' phonon energy with res that contains fano factor 
        EP.append(Pter1)
        EP_true.append(Pter_noF) # no smear, no fano factor. 
      
        # Find 'Measured' charge energy.
        Qer = np.random.normal(Eer,sig_q1) #smearing 'measured' charge energy with res that contains fano factor 
        EQ.append(Qer)
        EQ_true.append(Eer) 
       
        
     

        Erer = Pter1 - (V/eps/1000)*Qer # measured recoil energy 
        ER.append(Erer)

       # Y_er = yer_mu(Erer)
        
        Yield2 = Qer/Erer # "measured" Ionization Yield
        Yield_er.append(Yield2)
    
    df = pd.DataFrame({'E_true':E_true,'E_measured':ER,'Yield':Yield_er,'EP':EP,'EP_true':EP_true,'EQ_true':EQ_true,'EQ':EQ,'Sigp':sig_p,'Sigq':sig_q,'Sigp_noF':sig_p0})
    
    return df

