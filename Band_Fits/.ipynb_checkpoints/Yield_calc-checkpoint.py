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

#for Detector 2 

'''p_alpha = 0.082
p_beta = 0.0108
p_gamma = 0

q_alpha = 0.126
q_beta = 0.0066
q_gamma = 7.8*10**-5
'''
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
    
    #Er = np.random.uniform(10,200,N)
    #Er = np.random.exponential(40,np.uint32(N*0.3))
    #Er = Er[Er<=150]
    #Er = Er[Er>=10]
    #E1nr = E_true
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


       #With Fano Factor. 
        sig_p = np.sqrt(p_alpha + p_beta*Ptnr + p_gamma*(Ptnr**2))  #Phonon uncertainty (energy dependent)
        sig_q = np.sqrt(q_alpha + q_beta*Qnr + q_gamma*(Qnr**2)) #Charge uncertainty 
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
     
        
    return  Enr_true, ER, Yield, EQnr,EPnr,sigp_nr,sigq_nr 



def Yield_Er(Er_True):
    Yield_er = []
    ERer = []
    Er_true =[]
    EQ = []
    EP = []
    sig_p = []
    sig_q = []
   # bins  = np.array([10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110])
    #E1er = np.random.uniform(10,150,N)#from anthony, Er's are close enough to randomly distributed.
    #Eer = np.random.choice(E1er,N)
    #E1er = (bins[:-1] + bins[1:]) / 2
    #E1er = np.array([10.7,25.2,40.3,75.2])
    E1er = Er_True
    print("Bin center is:", E1er)

    QER = []
 
    for Eer in Er_True:


        Er_true.append(Eer)
        #Eer = 40

       # Pter = (1+(V/eps/1000))*Eer

        #More physical 
        
        N_e = Eer/eps 
        Er_fano = 0#0.13
        N_er = np.random.normal(0,np.sqrt(Er_fano *N_e)) + N_e  # 0.13 is the fano factor for germainium 
        #N_er = N_e
        Pter = Eer+(V/1000)*N_er
        
        sig_pee = np.sqrt(p_alpha + p_beta*Pter + p_gamma*(Pter**2)) #Phonon uncertainty 
        sig_qee = np.sqrt(q_alpha + q_beta*Eer + q_gamma*(Eer**2)) #Charge uncertainty
        
        sig_p.append(sig_pee)
        sig_q.append(sig_qee)
        
        #For Linear Fit Model
        '''
        alpha1 = 0.406241
        beta1 =0.00899969
        sig_qee = alpha1 + beta1*Eer   #Charge resolution 
        
        alpha2 = 0.125107
        beta2 = 0.0223536
        sig_pee = alpha2 + beta2*Pter1 #Phonon resolution '''
        
        
        
               

        Pter1 = np.random.normal(Pter,sig_pee)
      
        Qer = np.random.normal(Eer,sig_qee)

        EQ.append(Qer)
        EP.append(Pter1)
     

        Erer = Pter1 - (V/eps/1000)*Qer
        ERer.append(Erer)

       # Y_er = yer_mu(Erer)
        
        Yield2 = Qer/Erer
        Yield_er.append(Yield2)
        
    return Er_true,ERer, Yield_er, EQ, EP, sig_p, sig_q


