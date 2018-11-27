import numpy as np
import matplotlib.pyplot as plt
import resfuncRead as rfr
import pandas as pd
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from data_check import * 
from bin_data_check import *

ER = []
Yield = []

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

#for Detector 2 

'''p_alpha = 0.082
p_beta = 0.0108
p_gamma = 0

q_alpha = 0.126
q_beta = 0.0066
q_gamma = 7.8*10**-5
'''
def Yield_NR(N):
    
    import sys
    sys.path.insert(0,'/Users/Mitch 1/Desktop/nrFano_Constraint/python')
    import lindhard as lind
    lpar = lind.getLindhardPars('Ge',True)

    
    #Er = np.random.uniform(10,200,N)
    Er = np.random.exponential(40,np.uint32(N*0.3))
    Er = Er[Er<=150]
    Er = Er[Er>=10]

    for i in np.arange(0,N):

        Enr = np.random.choice(Er)
        

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


       #With Fano Factor. 
        sig_p = np.sqrt(p_alpha + p_beta*Ptnr + p_gamma*(Ptnr**2))  #Phonon uncertainty (energy dependent)
        sig_q = np.sqrt(q_alpha + q_beta*Qnr + q_gamma*(Qnr**2)) #Charge uncertainty 


        Fnr = np.random.normal(0.0,sig_p) #random sample assuming phonon variance 
        Fq = np.random.normal(0.0,sig_q) #random sampel assuming charge variance 


        Ptnr1 = Ptnr + Fnr 
        Qnr1 = Qnr + Fq 
        Ernr = Ptnr1 - (V/(eps*1000))*Qnr1 #Meaured ER
        yield1 = Qnr1 / Ernr  #Measured Yieldf



        #Store Stuff 
        ER.append(Ernr)
        Yield.append(yield1)
     
        
    return ER, Yield 



def Yield_Er(N):
    Yield_er = []
    ERer = []
    EQ = []
    bins  = np.array([10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110])
    E1er = np.random.uniform(10,150,N)#from anthony, Er's are close enough to randomly distributed. 
    #E1er = (bins[:-1] + bins[1:]) / 2

    QER = []

    for i in np.arange(N):

        Eer = np.random.choice(E1er) #randomly sample from Energy dist 
        #Eer = 40

        Pter = (1+(V/eps/1000))*Eer

        #More physical 
        ''' 
        N_e = Eer/eps 
        N_er = np.random.normal(0,np.sqrt(0.13*N_e)) + N_e  # 0.13 is the fano factor for germainium 
        #N_er = N_e
        Pter = Eer+(V/1000)*N_er'''
        
        sig_pee = np.sqrt(p_alpha + p_beta*Pter + p_gamma*(Pter**2)) #Phonon uncertainty 
        sig_qee = np.sqrt(q_alpha + q_beta*Eer + q_gamma*(Eer**2)) #Charge uncertainty

        #For Linear Fit Model
        
        
        
        
        '''
        alpha1 = 0.406241
        beta1 =0.00899969
        sig_qee = alpha1 + beta1*Eer   #Charge resolution 
        
        alpha2 = 0.125107
        beta2 = 0.0223536
        sig_pee = alpha2 + beta2*Pter1 #Phonon resolution '''
        
        
        
        
        Fer = np.random.normal(0.0,sig_pee) #Random energy assuming phonon variance 
        Fqe = np.random.normal(0.0,sig_qee) #Random assuming charge variane


        
        Pter = Pter + Fer
        #Qer = N_er*eps
        Qer = Eer
        Qer = Qer + Fqe
        EQ.append(Qer)
     

        Erer = Pter - (V/eps/1000)*Qer
        ERer.append(Erer)

        Y_er = yer_mu(Erer)
        
        Yield2 = Qer/Erer +(Y_er-1)
        Yield_er.append(Yield2)
        
    return ERer, Yield_er, EQ