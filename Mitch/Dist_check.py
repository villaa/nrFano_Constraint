import numpy as np 
import matplotlib.pyplot as plt
import resfuncRead as rfr
from scipy.stats import gaussian_kde
from prob_dist import * 


def dist_check(Yield,EP,EQ,sigP,sigQ,k):
    V = 4 #Voltage 
    eps = 0.003
    mu_p = np.mean(EP)
    sig_p = np.mean(sigP)
    #sig_p = sigP
    mu_q = np.mean(EQ)
    sig_q = np.mean(sigQ)
    #sig_q = sigQ
    
    res_p = mu_p/sig_p
    res_q = mu_q/sig_q
    r = sig_p/sig_q

    prob = ratio_dist(Yield, res_p, res_q, r, k)

    return prob

    #print(res_p,res_q,r)

def dist_check_fano(Yield,E,N_mean,Sp_mean,Sq_mean,SN):
    V = 4 

    N_mean = np.mean(N_mean)    
    Er = np.mean(E)
    Sp = np.mean(Sp_mean)
    Sq = np.mean(Sq_mean)
    SN = np.mean(SN)

    prob = ratio_dist_fano(Yield,Er,N_mean,Sp,Sq,SN,V/1000)
    
    return prob 