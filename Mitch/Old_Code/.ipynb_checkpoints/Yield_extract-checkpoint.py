import numpy as np
import matplotlib.pyplot as plt
import resfuncRead as rfr
import pandas as pd
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit



k = 0.157 # from lindhard
q = 1.602*10**-19 #electron charge 
V = 4.0 # voltage Bias
eps = .0033 #keV

#for detector 1.

p_alpha = 0.155393
p_beta = 0 #9.60343*10**(-11)
p_gamma = 0.000506287

q_alpha = 0.166004
q_beta = 0 #0.00233716
q_gamma = 9.52576*10**(-5)
'--------------------------'
ER = []
Yield = []
PtNr = []
QR = []
sigQ = []
sigP = []
sigQ1 = []
sigP1 = []
ENR = []
U1 = []
U2 = []
N_e_h = []
QNR=[]
PT1 = []
fano = []
'---------------------------'



def Yield_NR(N,Er):
    import sys
    sys.path.insert(0,'/Users/Mitch 1/Desktop/nrFano_Constraint/python')
    import lindhard as lind
    lpar = lind.getLindhardPars('Ge',True)

    
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
        Y = lind.getLindhard(lpar)


        #EdelWeiss Fano

        F = c**2/((eps*a/(A))+(2*V*a**2/(B*1000))+(2*(V/1000)**2*a**3/(eps*C)))


        #Constant Fano
        #F = 0

        Neh = Y(Enr*1000)*Enr/eps #number of electron-hole pairs. 

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

        yield1 = Qnr1 / Ernr #Measured Yield 



        #Derivatives for Fano
        U_1 = Ptnr/(Ptnr-Qnr1*(V/(eps*1000)))**2 #derivative with respect to Qnr1
        U_2 = -Qnr1/(Ptnr1 -(Qnr1*V/(eps*1000)))**2      #derivative with respect to Ptnf



        #Store Shit 
        ENR.append(Enr)
        N_e_h.append(Neh)
        sigQ.append(sig_q)
        sigP.append(sig_p)
        PT1.append(Ptnr1)
        QNR.append(Qnr1)
        ER.append(Ernr)
        Yield.append(yield1)
        U1.append(U_1)
        U2.append(U_2)
        fano.append(F)

        
    return ENR,N_e_h,sigQ,sigP,PT1,QNR,ER,Yield,U1,U2,fano


def Yield_ER(N):
    
    Yield_er = []
    ERer = []
    #x = np.arange(0,2000,0.1)
    #E1er = 82*np.exp(-x/304) #For ER from Kennedy Thesis
    E1er = np.random.uniform(0,200,N) #from anthony, Er's are close enough to randomly distributed. 

    QER = []

    for i in np.arange(N):

        Eer = np.random.choice(E1er) #randomly sample from Energy dist 

        Pter = (1+(V/eps/1000))*Eer

        N_e = Eer/eps 
        N_er = np.random.normal(0,np.sqrt(0.1*N_e)) + N_e


        sig_pee = np.sqrt(p_alpha + p_beta*Pter + p_gamma*(Pter**2)) #Phonon uncertainty 
        sig_qee = np.sqrt(q_alpha + q_beta*Eer + q_gamma*(Eer**2)) #Charge uncertainty 

        Fer = np.random.normal(0.0,sig_pee) #Random energy assuming phonon variance 
        Fqe = np.random.normal(0.0,sig_qee) #Random assuming charge variane


        Pter = Eer+(V/1000)*N_er
        Pter = Pter + Fer
        Qer = N_er*eps
        Qer = Qer + Fqe
        QER.append(Qer)

        Erer = Pter - (V/eps/1000)*Qer
        ERer.append(Erer)

        Yield2 = Qer/Erer
        Yield_er.append(Yield2)
        
    return Yield_er, ERer 