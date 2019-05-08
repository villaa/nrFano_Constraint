
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
from ER_Yield import Yield_Er
from NR_Yield import Yield_NR
from bin_data_check import * 
from Dist_check import *  
from Data_check_continuous import *


#%%
N = 50000
s = 1
fano =5 # 'known' fano factor for electron recoils. 

bins = np.array([10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110])
bins_cont = np.linspace(10,150,N)

Eer = np.random.choice(bins,N)

df,DF= Yield_Er(Eer,fano) #Electron Recoil Band with fano (BINNED)
df_count,df2 = Yield_Er(bins_cont,fano) #indep

#df,DF= Yield_NR(Eer) #Nuclear Recoil Band with fano (BINNED)
#df_count,df2 = Yield_NR(bins_cont) #NR band for continuos dist of energies 


#%%
expected_v1, expected_v2, Er_true = continuous_containment(df_count,s,band_er)

#%%

cut_idx = 'E_true' # True energy
#cut_idx = 'E_measured' # Measured Energy
er = np.arange(0,2,0.002) #electron recoils 
nr = np.linspace(0.1,0.5,1000) #for nuclear recoils. 

df,bincenters = bin_check(df,1,band_er,bins,cut_idx,expected_v1,expected_v2,Er_true,fano,er)






#%%
'''
p_alpha = 0.155393
p_beta = 9.60343*10**(-11)
p_gamma = 0.000506287
q_alpha = 0.166004
q_beta = 0.0023371
q_gamma = 9.52576*10**(-5)


def sigmaq(Eq):
    return np.sqrt(q_alpha + q_beta*Eq + q_gamma*(Eq**2))
def sigq_der(Eq):
    return 0.5*np.power(q_alpha + q_beta*Eq + q_gamma*(Eq**2),-0.5)*(q_beta + 2*q_gamma*(Eq))

def sigmap(Ep):
    return np.sqrt(p_alpha + p_beta*Ep + p_gamma*(Ep**2))
def sigp_der(Ep):
    return 0.5*np.power(p_alpha + p_beta*Ep + p_gamma*(Ep**2),-0.5)*(p_beta + 2*p_gamma*(Ep))

Eq = np.linspace(10,170,100000)
Ep = np.linspace(20,250,100000)




plt.figure(figsize=(8, 6))
#plt.plot(Eq,sigq_der(Eq))
plt.plot(Ep,df_count.sig_p_var,label = 'Varried N')
plt.plot(Ep,df_count.sigp_mean,label = 'Varried SigP')
plt.title("simgap compare ")
plt.ylabel('$\sigma_q^|$')
plt.xlabel('$E_q$')
plt.legend()
plt.grid(True)


plt.show()
'''

#%%
'''
m,c,r,p,se1= stats.linregress(df_count.EP_smear,df_count.EQ_smear)

cm1lab="$"+('y=%2.2fx+%2.2f, r^2=%1.2f'%(m,c,r**2))+"$"

plt.figure(figsize=(9.0,8.0))
plt.plot(df_count.EP_smear, m*df_count.EP_smear+c,'r--',label = cm1lab)
plt.scatter(df_count.EP_smear,df_count.EQ_smear,s= 2.5,label = 'Data')

plt.title('$E_P$ vs. $E_Q$')
plt.xlabel('$E_P$')
plt.ylabel('$E_Q$')
plt.grid(linestyle='-', linewidth=0.35)
plt.legend()

plt.savefig('/Users/Mitch 1/Desktop/EP_EQ_FIT.png')

plt.show()
'''

