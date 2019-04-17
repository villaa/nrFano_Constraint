import numpy as np
from scipy.special import erf

# returns the probability of z, where z is the yield
# z = Eq/(Ep - k*Eq)
# this distribution assumes that Eq and Ep are independent, 
# Eq has mean mu_q and width sig_q and is gaussian
# Ep has mean mu_p and width sig_p and is gaussian
# res_p = mu_p/sig_p
# res_q = mu_q/sig_q
# r = sig_p/sig_q
# k is defined as e*voltage/(energy needed to create one e/h pair)
# so for a Si detector at e.g. 4V, k ~ 4/3 
def ratio_dist(z, res_p, res_q, r, k):
    F1 = np.exp(-0.5*(res_q**2 + res_p**2)) / (np.pi*(r*z**2 + (1/r)*(1+k*z)*2))
    
    G11 = (r*(z*res_q*r + (1+k*z)*res_p)) / (np.sqrt(2*np.pi)*np.power(z**2 * r**2 + (1+k*z)**2, 3/2))
    G12 = np.exp(-(z*res_p*r - (1+k*z)*res_q)**2 / (2*(z**2 * r**2 + (1+k*z)**2)))
    G13 = erf((z*res_q + (1+k*z)*res_p/r) / np.sqrt(2*(z**2 + (1+k*z)**2 / r**2)))

    return F1 + G11*G12*G13

# x is the yield, Eq/Er
# Er is the recoil energy, in units of keV
# meanN is the mean number of e/h pairs created given Er
# sdP is the standard deviation of the phonon signal, units of ??
# sdQ is the standard deviation of the charge signal, units of ??
# sdN is the standard deviation of the number of electron-hole pairs, unitless
# V is the voltage across the detector, in units of kV??
def ratio_dist_fano(x, Er, meanN, sdP, sdQ, sdN, V):
    e = 1 #charge of electron
    k = (sdP**2)*(sdQ**2)+(V**2)*(sdQ**2)*(sdN**2)+(e**2)*(sdN**2)*(sdP**2)

    def g(x):
        ans = (1/math.sqrt(math.pi))+x*(np.exp(x**2))*math.erf(x)
        return ans
        
    def A(x):
        ans = ((((x*(V/e)+1)*sdQ)**2)+((x*sdP)**2)+((e*sdN)**2))/(2*k)
        return ans
        
    def B(x):
        ans = ((V/e)*(sdQ**2)*(Er*x+e*meanN)+x*e*meanN*(((V*sdQ/e)**2)+(sdP**2))+Er*((sdQ**2)+((e*sdN)**2)))/(k)
        return ans

    C = ((((meanN*V+Er)*sdQ)**2)+(((meanN*sdP)**2)+((Er*sdN)**2))*(e**2))/(2*k)

    ans = (1/(2*math.sqrt(math.pi*k)))*(1/A(x))*np.exp(-C)*g((B(x))/(2*math.sqrt(A(x))))
    return ans