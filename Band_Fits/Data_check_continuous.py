
import resfuncRead as rfr
import numpy as np
from Yield_calc import * 
import pandas as pd 
from data_check import *
from matplotlib.ticker import FuncFormatter
from prob_dist import * 
from Dist_check import * 
from scipy import integrate 
from Hist_plot import * 


from tabulate import tabulate

    
        
def continuous_containment(df,s,band_func):
    
    V = 4
    eps = 0.0033
    expected = []

   # print("before sort Er_true is:",Er_true)
    recoil_type, upper,lower = band_func(df.E_true,s)

    #Data = np.vstack((Er,Yield,upper,lower,EP,EQ,sigp,sigq,Er_true)).T
    #Data1 = Data[np.argsort(Data[:, 0])]
    
    df1 = pd.DataFrame({'upper':upper,'lower':lower})
    
    df = pd.concat([df,df1],axis=1)



    E = np.array(df.E_measured)
    Er_true = np.array(df.E_true)
    x = np.array(df.Yield) # Yield 
    y = np.array(df.upper) #upper bound
    z = np.array(df.lower) #lower bound
    P = np.array(df.EP) #EPe
    P_true = np.array(df.EP_true)
    Q = np.array(df.EQ) #EQ
    Q_true = np.array(df.EQ_true)
    Sp = np.array(df.Sigp) #Sigma_p
    Sq = np.array(df.Sigq)
    Sp_noF = np.array(df.Sigp_noF) #Sigma_q


            #print(x)

            #look at distributions graphically f
    k= (V/eps/1000)
    u = np.arange(0,2,0.002) #electron recoils 
           # u = np.linspace(0.1,0.5,1000) #for nuclear recoils. 




            #v = np.linspace(np.mean(z),np.mean(y),1000)
            #g = np.trapz(prob,v)
            
    for a,b,c,d,e,f in zip(P_true,Q_true,Sp_noF,Sq,z,y): 
        
        g = integrate.quad(lambda x: dist_check(x,a,b,c,d,k),np.mean(e),np.mean(f) )
        H =g[0]*100
           
        expected.append(H)


          

    return expected, Er_true


