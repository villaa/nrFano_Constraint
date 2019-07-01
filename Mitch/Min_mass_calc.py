'''This file calculated the minimum accessible mass based on the fano factor'''
#%%
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
import seaborn as sns 
import matplotlib as mpl
mpl.rcParams.update(_VSCode_defaultMatplotlib_Params)
sns.set_style(rc = {'figure.facecolor':'white'})


from ER_Yield import Yield_Er
from NR_Yield import Yield_NR
from Intersection_loc import * 
from WIMP_Mass import * 

import numpy as np
import matplotlib.pyplot as plt






#%%

N =100000
 
Fano = 0.13 #electron recoil fano factor 
NR_Fano = np.linspace(0,100,11)
bins = np.array([1,3,5,8,12,15])
bins1 = np.array([0,1,3,5,8,12,15])
WIMP_M = []
ER_Min = []
intersection = []

for x in NR_Fano:  

    Eer = np.random.choice(bins,N)

    df_er,DF_er= Yield_Er(Eer,Fano) #Electron Recoil Band with fano (BINNED)

    df,DF_nr= Yield_NR(Eer,x) #Nuclear Recoil Band with fano (BINNED)



  #  bins = np.array([0,5,10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110])
    V = 4
    eps = 0.0033

    Intersection = Intersection_loc(DF_er,DF_nr,'E_true',bins1,x) 


    EP = Intersection.x 
    EQ = Intersection.y 

    Er_min = EP - (V/eps/1000)*EQ
    ER_Min.append(Er_min)

    W_mass = WIMP_MASS(Er_min)
    WIMP_M.append(W_mass)


    print('Minimum Mass for Fano = ',x,'is: ',W_mass, '[GeV]')



#%%

m,b,r,p,se= stats.linregress(NR_Fano,WIMP_M)

line = m*NR_Fano + b 

equation = 'y = ' + str(round(m,4)) + 'x' ' + ' + str(round(b,4))


plt.figure(figsize=(9.0,8.0))

#plt.scatter(NR_Fano, ER_Min)
plt.scatter(NR_Fano, WIMP_M,s = 35)
plt.plot(NR_Fano, line,'--',color = 'black')

plt.title('Minimum Accessable Mass',size = 14)
plt.xlabel('Nuclear Recoil Fano Factor',size = 14)
#plt.ylabel('Crossover Energy [keV]')
plt.ylabel('WIMP Mass [GeV]',size = 14)
plt.grid(linestyle='-', linewidth=0.35)
plt.text(70,4.6,equation)


plt.savefig('/Users/Mitch 1/Desktop/Min_Mass_Fano.png')



#%%

m1,b1,r1,p1,se1= stats.linregress(NR_Fano,ER_Min)

line1 = m1*NR_Fano + b1 

equation = 'y = ' + str(round(m,4)) + 'x' ' + ' + str(round(b,4))


plt.figure(figsize=(9.0,8.0))

plt.scatter(NR_Fano, ER_Min)
#plt.scatter(NR_Fano, intersction,s = 35)
plt.plot(NR_Fano, line1,'--',color = 'black')

plt.title('Crossover Energy vs. Fano Factor',size = 16)
plt.xlabel('Nuclear Recoil Fano Factor',size = 14)
plt.ylabel('Crossover Energy [keV]',size = 14)
plt.grid(linestyle='-', linewidth=0.35)
plt.text(70,3.8,equation)


plt.savefig('/Users/Mitch 1/Desktop/Min_Er_Fano.png')

#%%
