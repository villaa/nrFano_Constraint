
#%%
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
import seaborn as sns 
import matplotlib as mpl
mpl.rcParams.update(_VSCode_defaultMatplotlib_Params)
sns.set_style(rc = {'figure.facecolor':'white'})

import numpy as np
import matplotlib.pyplot as plt
import resfuncRead as rfr
from scipy.stats import kde
from band_check import *
from Band_plots import *
from Eer_Yield import Yield_Er
from NR_Yield import Yield_NR
from bin_data_check_v1 import * 
from bin_data_check_v2 import * 
from Dist_check import *  
from Data_check_continuous import *
from EpEq_space import * 
 


#%%
N =100000
s = 1
#fano =0.13
fano = 'EDW' # 'known' fano factor for electron recoils. 
Fano = 0 #electron recoil fano factor 

bins = np.array([10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110,])
bins_cont = np.linspace(10,150,N)

bins_cont = np.random.choice(bins_cont,N)

Eer = np.random.choice(bins,N)

#df_er_v1,df_er_v2= Yield_Er(Eer,Fano) #Electron Recoil Band with fano (BINNED)
#df_er_v1_count,df_er_v2_count = Yield_Er(bins_cont,Fano) #indep

df_nr_v1,df_nr_v2= Yield_NR(Eer) #Nuclear Recoil Band with fano (BINNED)
df_nr_v1_count,df_nr_v2_count = Yield_NR(bins_cont) #NR band for continuos dist of energies 


#%%
NR_band_plot(df_nr_v2,1)
#NR_band_plot(df_nr_v2_count,2)

#ER_band_plot(df_er_v1,1)
#ER_band_plot(df_er_v2,2)


#%%
'''For Nuclear Recoils'''
expected_v1, expected_v2, Er_true = continuous_containment(df_nr_v1_count,df_nr_v2_count,s,band_nr)
'''For Electron Recoils '''
#expected_v1, expected_v2, Er_true = continuous_containment(df_er_v1_count,df_er_v2_count,s,band_er)


#%%

cut_idx = 'E_true' # True energy
#cut_idx = 'E_measured' # Measured Energy
er = np.arange(0,2,0.002) #electron recoils 
nr = np.linspace(0.1,0.5,1000) #for nuclear recoils. 

'''For Nuclear Recoils'''
df1,bincenters1 = bin_check_v1(df_nr_v1,1,band_nr,bins,cut_idx,expected_v1,Er_true,fano,nr)
#df2,bincenters2 = bin_check_v2(df_nr_v2,1,band_nr,bins,cut_idx,expected_v2,Er_true,fano,nr)
'''For Electron Recoils'''
#df1,bincenters1 = bin_check_v1(df_er_v1,1,band_er,bins,cut_idx,expected_v1,Er_true,Fano,er)
#df2,bincenters2 = bin_check_v2(df_er_v2,1,band_er,bins,cut_idx,expected_v2,Er_true,Fano,er)







#%%
'''Continous EQ_EP Space Splot '''
plt.figure(figsize=(9.0,8.0))
plt.scatter(df2_nr.EP,df2_nr.EQ,s = 10,label = 'Nuclear Recoils')
plt.scatter(df2_er.EP,df2_er.EQ,color = 'black', s = 10,label = 'Electron Recoils')

plt.title('Simulated $E_P$ $E_Q$ Space',size = 16)
plt.xlabel('$E_P$ [keV]',size = 14)
plt.ylabel('$E_Q$ [keV]',size = 14)
plt.grid(linestyle='-', linewidth=0.35)
plt.legend()

plt.savefig('/Users/Mitch 1/Desktop/Thesis_Plots/EP_EQ_Space.png')

plt.show()




#%%
'''Fano Factor as a Function of ER (edel) plot'''
plt.figure(figsize=(9.0,8.0))
plt.plot(df_cont.E_true,df_cont.Fano)
plt.title('Nuclear Recoil Fano Factor',size = '18')
plt.xlabel('Recoil Energy [keV]',size = '16')
plt.ylabel('Fano Factor',size = '16')
plt.grid(linestyle='-', linewidth=0.35)
plt.savefig('/Users/Mitch 1/Desktop/Thesis_Plots/Fano_Factor.png')
plt.show() 

#%%


