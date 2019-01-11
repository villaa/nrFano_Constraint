import resfuncRead as rfr
import numpy as np
from Yield_calc import * 
import pandas as pd 
from data_check import *
from matplotlib.ticker import FuncFormatter




def Corr(EP,EQ): 
    
    Data = np.vstack((EP,EQ)).T
    Data1 = Data[np.argsort(Data[:, 0])]

    df = pd.DataFrame(Data1)
    Corr1 = []
    bins  = [10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110]
    df1 = pd.cut(df[0],bins, labels = False)
    bin_spacing = [10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110]
    for q in np.arange(0,8):
        x = np.array(df[0].loc[df1 == q])
        y = np.array(df[1].loc[df1 == q])

        Corr = np.corrcoef(x,y)
        
        Corr1.append(Corr[0,1])
        

        #print('The correlation coefficent for the',bin_spacing[q],'keV bin is:',Corr[0,1])

  
    #Brute Force (Check)

    Coeff = []
    for q in np.arange(0,9):
        x = np.array(df[0].loc[df1 == q])
        y = np.array(df[1].loc[df1 == q])

        S_XX = np.sum(x**2) - (np.sum(x))**2/(len(x))
        S_YY = np.sum(y**2) - (np.sum(y))**2/(len(y))
        S_XY = np.sum(x*y) - ((np.sum(x))*(np.sum(y)))/len(x) 

        p = S_XY/(np.sqrt((S_XX*S_YY)))

        Coeff.append(p)

    #print(Coeff)

    print('CORRELATION COEFFICENT FOR EQ AND EP')
    print('---------------------------------------------------------------')
    print('Bin(keV)','\t', "Numpy Coefficient", '\t', 'Brute Force Check')
    print('---------------------------------------------------------------')
    bin_spacing = '10-13.4','13.4-18.1','18.1-24.5','24.5-33.1','33.1-44.8','44.8-60.6','60.6-80.2','80.2-110','110-150'


    for x,y,z in zip(bin_spacing,Corr1,Coeff):
        print(x, '\t','\t', '{0:1.2f}'.format(y), '\t','\t','\t','{0:1.2f}'.format(z))