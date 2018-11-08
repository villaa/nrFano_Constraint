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


def band_check_nr(a,s,X):
    upper = ynr_mu(X)+s*ynr_sigv(X)
    lower = ynr_mu(X)-s*ynr_sigv(X)
    count = []
    count1 = [] 
    num = []
    num1 = [] 
    

    for x,y in zip(a, upper): 
        if x >= y:
            n = 1
        else:
            n = 0 
        count.append(n)
    
    for x,y in zip(a,lower):
        if x <= y:
            n = 1
        else:
            n = 0
        count1.append(n)       
        
        
    for i in count: 
        if i == 1:
            #print(x)
            num.append(i)     

    
    for i in count1: 
        if i == 1:
            #print(x)
            num1.append(i)     
            
    N= len(a)
    percent = 100*(N - np.abs(len(num)+len(num1)))/N
    print('For',N,'number of data points')
    print('Perecnt of data in nuclear band for',s,'simga =',percent)
    print('Number of data points outside bands is:',len(num)+len(num1))
    print('-----------------------------------------------------------')
    
    
def band_check_er(b,s,Y):
    upper = yer_mu(Y)+s*yer_sigv(Y)
    lower = yer_mu(Y)-s*yer_sigv(Y)
    counte = []
    counte1 = [] 
    nume = []
    nume1 = [] 
    

    for x,y in zip(b, upper): 
        if x >= y:
            n = 1
        else:
            n = 0 
        counte.append(n)
    
    for x,y in zip(b,lower):
        if x <= y:
            n = 1
        else:
            n = 0
        counte1.append(n)       
        
        
    for i in counte: 
        if i == 1:
            #print(x)
            nume.append(i)     

    
    for i in counte1: 
        if i == 1:
            #print(x)
            nume1.append(i)     
            
    N = len(b)
    percent = 100*(N - np.abs(len(nume)+len(nume1)))/N
    print('For',N,'number of data points')
    print('Perecnt of data in electron band for',s,'simga =',percent)
    print('Number of data points outside bands is:',len(nume)+len(nume1))
    print('-----------------------------------------------------------')