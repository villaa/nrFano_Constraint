
import resfuncRead as rfr
import numpy as np 
import pandas as pd 
from data_check import *
from matplotlib.ticker import FuncFormatter
from prob_dist import * 
from Dist_check import * 
from scipy import integrate 
from Hist_plot import * 


from tabulate import tabulate

    
        
def continuous_containment(df,s,band_type):
    
    V = 4
    eps = 0.0033
    expected = []
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
    Sp_mean = np.array(df.sigp_mean) #Sigma_p
    Sq_mean = np.array(df.sigq_mean)
   
    k= (V/eps/1000)
    u = np.arange(0,2,0.002) #electron recoils 
  # u = np.linspace(0.1,0.5,1000) #for nuclear recoils. 

            
    for a,b,c,d,e,f in zip(Ep_mean,Eq_mean,Sp_mean,Sq_mean,z,y): 
        
        g = integrate.quad(lambda x: dist_check(x,a,b,c,d,k),np.mean(e),np.mean(f) )
        H =g[0]*100
           
        expected.append(H)


          

    return expected, Er_true


