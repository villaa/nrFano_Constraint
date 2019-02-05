import resfuncRead as rfr
import numpy as np

bpar_er = rfr.getBandFunc('data/erband_R133') #reads in band data/fit for er data
bpar_nr = rfr.getBandFunc('data/nrband_R133')

ynr_mu = rfr.makeBFunc(bpar_nr[1]['mu'])
ynr_muv = np.vectorize(ynr_mu)
ynr_sig = rfr.makeBFunc(bpar_nr[1]['sig'],True)
ynr_sigv = np.vectorize(ynr_sig)

yer_mu = rfr.makeBFunc(bpar_er[1]['mu']) # sets average fit from 1st (2nd) col in data table. 
yer_muv = np.vectorize(yer_mu) #puts mean data for er in 1D array 
yer_sig = rfr.makeBFunc(bpar_er[1]['sig'],True) #sets uncertainty 
yer_sigv = np.vectorize(yer_sig) #puts uncertainty into 1D array


def band_nr(Er,s):
    upper = ynr_mu(Er)+s*ynr_sigv(Er)
    lower = ynr_mu(Er)-s*ynr_sigv(Er)
    
    return "NUCLEAR",upper,lower
    
def band_er(Er,s):
    #upper = yer_mu(Er)+s*yer_sigv(Er)
    #lower = yer_mu(Er)-s*yer_sigv(Er) 

    upper = 1+s*yer_sigv(Er)
    lower = 1-s*yer_sigv(Er)
    #print("S should be:",s)
    #print("The True Er should be:",Er)
    #print("in function upper value is:",upper)
    return "ELECTRON",upper, lower 



def compare(Yield,upper,lower):
    up = np.sum(Yield>upper)
    down = np.sum(Yield<lower)
    N = len(Yield)
    
    return up,down,N


def band_check(Yield,Er,s,band_func):
    
    recoil_type, upper,lower = band_func(Er,s)
    
    up,down,N = compare(Yield,upper,lower)
    
    percent = 100*(N - np.abs(up+down))/N
    
    
    print('For',N,'number of data points')
    print('Perecnt of data in', recoil_type, 'for',s,'simga =',percent)
    print('Number of data points outside bands is:', up+down)
    print('-----------------------------------------------------------')

    

    
