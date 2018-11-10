
import resfuncRead as rfr
import numpy as np
from Yield_calc import * 
import pandas as pd 
from data_check import *

from tabulate import tabulate




def bin_check(Yield,Er,s,band_func):
    
    
    Percent = []
    
    recoil_type, upper,lower = band_func(Er,s)

    Data = np.vstack((Er,Yield,upper,lower)).T
    Data1 = Data[np.argsort(Data[:, 0])]
    
    df = pd.DataFrame(Data1)

    bins  = [10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110,150]

    df1 = pd.cut(df[0],bins, labels = False)

    for q in np.arange(0,9):

    

     
        x = np.array(df[1].loc[df1 == q])
        y = np.array(df[2].loc[df1 == q])
        z = np.array(df[3].loc[df1 == q])


        up,down,N = compare(x,y,z)

        percent = 100*(N - (up+down))/N

        Percent.append(percent)
        
    
    print('--------------------------------------------')  
    print(s,"SIGMA",recoil_type, "RECOIL BAND")
    print('--------------------------------------------')  
    bin_spacing = '10-13.4','13.4-18.1','18.1-24.5','24.5-33.1','33.1-44.8','44.8-60.6','60.6-80.2','80.2-110','110-150'

    print("Bin Spacing (keV)", '\t', "Percent in band")     #table column headings
    print("--------------", '\t', '\t' "-------------------")


    for x,y in zip(bin_spacing,Percent):
        print(x, '\t','\t', '{0:1.2f}'.format(y))
            
    print('--------------------------------------------')
    
    