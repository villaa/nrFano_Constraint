
#%%
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
import numpy as np
import matplotlib.pyplot as plt
import resfuncRead as rfr
from scipy.stats import kde
from data_check import *
from Band_plots import *
from ER_Yield import Yield_Er 
from NR_Yield import Yield_Nr
from bin_data_check import * 
from Dist_check import *  
from Data_check_continuous import * 


#%%
from prob_dist import * 
from Dist_check import * 

N_mean = 40000
SN = 476
Sq = 2
Sp = 7
Er = 149
V = 0.004
e = 0.003

Yield = np.linspace(0,2,1000)


'''#parameters
N_mean = 10
SN = 1
Sq = 1
Sp = 1
Er = 10
V = 1
e = 1'''

prob = ratio_dist_fano(Yield,Er,N_mean,Sp,Sq,SN,V,e)


#%%
plt.plot(Yield,prob)
plt.show()




