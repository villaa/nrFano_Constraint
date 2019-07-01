import resfuncRead as rfr
import numpy as np 
import pandas as pd 
from band_check import *
from prob_dist import * 
import matplotlib.pyplot as plt

import scipy.stats as stats



def EpEQ_bin(df1,df2,cut_idx,bins): 

    EQ_er_edge = []
    EP_er_Center = []
    EQ_nr_edge = []
    EP_nr_Center = []


    
    df1['bin'] = pd.cut(df1[cut_idx],bins)
    df2['bin'] = pd.cut(df2[cut_idx],bins)
    

    for bin_name, bin_data in df1.groupby('bin'):

        
        if len(bin_data) > 1: 

            EQ_er = np.array(bin_data.EQ)
            EP_er = np.array(bin_data.EP)
            #EQ_nr = np.array(df2.EQ)
            #EP_nr = np.array(df2.EP) 



            EP_er_center = np.mean(EP_er)
            EP_er_Center.append(EP_er_center)

            x = np.mean(EQ_er)-2*np.std(EQ_er)
            EQ_er_edge.append(x)


    for bin_name, bin_data1 in df2.groupby('bin'):

        
        if len(bin_data1) > 1: 


            EQ_nr = np.array(bin_data1.EQ)
            EP_nr = np.array(bin_data1.EP)
  

          

            Ep_nr_center = np.mean(EP_nr)
            EP_nr_Center.append(Ep_nr_center)

            y = np.mean(EQ_nr)+2*np.std(EQ_nr)
            EQ_nr_edge.append(y)

    return EP_er_Center, EQ_er_edge, EP_nr_Center, EQ_nr_edge



