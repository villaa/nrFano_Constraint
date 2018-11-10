import resfuncRead as rfr
import numpy as np
from Yield_calc import * 
import pandas as pd 

from tabulate import tabulate

bpar_er = rfr.getBandFunc('data/erband_R133') #reads in band data/fit for er data
bpar_nr = rfr.getBandFunc('data/nrband_R133')

ynr_mu = rfr.makeBFunc(bpar_nr[1]['mu'])
ynr_muv = np.vectorize(ynr_mu)
ynr_sig = rfr.makeBFunc(bpar_nr[1]['sig'],True)
ynr_sigv = np.vectorize(ynr_sig)

yer_mu = rfr.makeBFunc(bpar_er[1]['mu']) # sets average fit from 1st (2nd) col in data table. 
yer_muv = np.vectorize(yer_mu) #puts mean data for er in 1D array 
yer_sig = rfr.makeBFunc(bpar_er[1]['sig'],True) #sets uncertainty 
yer_sigv = np.vectorize(yer_sig)


def bin_check_NR(s,N,ER,Yield):
    
    Percent = []


    #for nuclear recoils 
    upper = ynr_mu(ER)+s*ynr_sigv(ER)
    lower = ynr_mu(ER)-s*ynr_sigv(ER)

    Data = np.vstack((ER,Yield,upper,lower)).T
    Data1 = Data[np.argsort(Data[:, 0])]
    
    df = pd.DataFrame(Data1)

    bins  = [10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110,150]

    df1 = pd.cut(df[0],bins, labels = False)

    for q in np.arange(0,9): # q is the bin number. 

        count = []
        num = []
        count1 = []
        num1 = []


        x = np.array(df[1].loc[df1 == q]) #yield
        y = np.array(df[2].loc[df1 == q]) #upper band
        z = np.array(df[3].loc[df1 == q]) #lower band 


        for i,j in zip(x,y):
            if i > j:
                n = 1
            else: 
                n = 0 

            count.append(n)
        #print(count)

        for i in count: 
            if i == 1:
                num.append(i) #number of times the data point is outside upper band


        for i,j in zip(x,z):
            if i < j:
                n = 1
            else: 
                n = 0 

            count1.append(n)
        #print(count)

        for i in count: 
            if i == 1:
                num1.append(i) #number of times the data point is outside lower band


        N= len(x)
        percent = 100*(N - np.abs(len(num)+len(num1)))/N

        Percent.append(percent)
        
    
    print('--------------------------------------------')  
    print(s,"SIGMA NUCLEAR RECOIL BAND")
    print('--------------------------------------------')  
    bin_spacing = '10-13.4','13.4-18.1','18.1-24.5','24.5-33.1','33.1-44.8','44.8-60.6','60.6-80.2','80.2-110','110-150'

    print("Bin Spacing (keV)", '\t', "Percent in band")     #table column headings
    print("--------------", '\t', '\t' "-------------------")


    for x,y in zip(bin_spacing,Percent):
        print(x, '\t','\t', '{0:1.2f}'.format(y))
            
    print('--------------------------------------------')       
            