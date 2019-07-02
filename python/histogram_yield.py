import numpy as np
from scipy.special import erf
import math 
from scipy.integrate import quad
import resfuncRead as rfr
import scipy.optimize as so
import pandas as pds





def QEr_Ebin(Q, Ernr, bins=[5, 10, 20, 30, 40, 50, 70,150]):

    
    #create a dataframe
    nr_df = pds.DataFrame(data={'yield':Q, 'energy':Ernr})

    #bin the data
    nr_df['binned'] = pds.cut(nr_df['energy'],bins)

    #print stats in each bin
    s = nr_df.groupby(pds.cut(nr_df['energy'], bins=bins)).size()
    #print (s)

    #create list of vectors for histogrammin'
    bindf = nr_df.groupby(pds.cut(nr_df['energy'], bins=bins))['yield'].apply(list)
    #print(bindf)

    return bindf 

def QEr_Qhist(bindf, qbins=np.linspace(0,0.6,40)):

 
    xcq = (qbins[:-1] + qbins[1:]) / 2

    qhistos = np.zeros((np.shape(qbins)[0]-1,0))
    qerrs = np.zeros((np.shape(qbins)[0]-1,0))
   
    for i,Qv in enumerate(bindf):
      n,nx = np.histogram(Qv,bins=qbins)
      n = np.reshape(n,(np.shape(n)[0],1))
      qhistos = np.append(qhistos,n,axis=1)
      qerrs = np.append(qerrs,np.sqrt(n),axis=1)
      qerrs[qerrs==0]=1

    return qhistos,qerrs
