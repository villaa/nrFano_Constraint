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

    

    print(res_p,res_q,r)

    
    prob = ratio_dist(Yield, res_p, res_q, r, k)
    
    return prob 