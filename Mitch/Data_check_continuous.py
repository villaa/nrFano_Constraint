import resfuncRead as rfr
import numpy as np 
import pandas as pd 
from band_check import *
from matplotlib.ticker import FuncFormatter
import sys
sys.path.append('../python')
import prob_dist as PD
from Dist_check import * 
from scipy import integrate 
from Hist_plot import * 



    
        
def continuous_containment(df,s,band_type):
    
    V = 4
    eps = 0.0033
    expected_v1 = []
    expected_v2= [] 
   # print("before sort Er_true is:",Er_true)
    recoil_type, upper,lower = band_type(df.E_true,s)
    
    df1 = pd.DataFrame({'upper':upper,'lower':lower})
    
    df = pd.concat([df,df1],axis=1)


    E = np.array(df.E_measured)
    Er_true = np.array(df.E_true)
    x = np.array(df.Yield) # Yield 
    y = np.array(df.upper) #upper bound
    z = np.array(df.lower) #lower bound
    Ep_mean = np.array(df.Ep_mean)
    Eq_mean = np.array(df.Eq_mean)

    '''For Independent PDF'''
    Sp_v1 = np.array(df.sigp_v1) #Sigma_p
    Sq_v1 = np.array(df.sigq_v1)

    '''For Dependent pdf '''
    Sp_v2 = np.array(df.sigp_v2) #Sigma_p
    Sq_v2 = np.array(df.sigq_v2)
    SN = np.array(df.sig_N)
   



    N_mean = np.array(df.N_mean) #Sigma_q


    k= (V/eps/1000)
    u = np.arange(0,2,0.002) #electron recoils 
  # u = np.linspace(0.1,0.5,1000) #for nuclear recoils. 

    '''Calculate Expected Independent Containment Fraction'''
    for a,b,c,d,e,f in zip(Ep_mean,Eq_mean,Sp_v1,Sq_v1,z,y):
      f = integrate.quad(lambda x: dist_check_v1(x,a,b,c,d,k),np.mean(e),np.mean(f))

      h = f[0]*100

      expected_v1.append(h)



    '''Calculate Expected Dependent Containment Fraction'''
    for a,b,c,d,e,f,g in zip(Er_true,N_mean,Sp_v2,Sq_v2,SN,z,y):
        
        g = integrate.quad(lambda x: dist_check_v2(x,a,b,c,d,e),np.mean(f),np.mean(g))

        H =g[0]*100
           
        expected_v2.append(H)




      

    return expected_v1, expected_v2, Er_true