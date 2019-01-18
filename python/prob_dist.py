import numpy as np

# returns the probability of z, where z is the yield
# z = Eq/(Ep - k*Eq)
# this distribution assumes that Eq and Ep are independent, 
# Eq has mean mu_q and width sig_q and is gaussian
# Ep has mean mu_p and width sig_p and is gaussian
# res_p = mu_p/sig_p
# res_q = mu_q/sig_q
# r = sig_p/sig_q
# k is defined as e*voltage/(energy needed to create one e/h pair)
def ratio_dist(z, res_p, res_q, r, k):
    F1 = np.exp(-0.5*(res_q**2 + res_p**2)) / (np.pi*(r*z**2 + (1/r)*(1+k*z)*2))
    
    G11 = () / ()
