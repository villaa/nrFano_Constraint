import numpy as np 
import matplotlib.pyplot as plt
import resfuncRead as rfr


bpar_er = rfr.getBandFunc('data/erband_R133') #reads in band data/fit for er data
bpar_nr = rfr.getBandFunc('data/nrband_R133') #reads in band data/fit for nr data 

yer_mu = rfr.makeBFunc(bpar_er[1]['mu']) # sets average fit from 1st (2nd) col in data table. 
yer_muv = np.vectorize(yer_mu) #puts mean data for er in 1D array 
yer_sig = rfr.makeBFunc(bpar_er[1]['sig'],True) #sets uncertainty 
yer_sigv = np.vectorize(yer_sig) #puts uncertainty into 1D array

#following does the same but for nuclear recoils. 
ynr_mu = rfr.makeBFunc(bpar_nr[1]['mu'])
ynr_muv = np.vectorize(ynr_mu)
ynr_sig = rfr.makeBFunc(bpar_nr[1]['sig'],True)
ynr_sigv = np.vectorize(ynr_sig)



def NR_ER_Band(Enr,Yield,Erer,Yield_er):

    X = Enr
    X = np.sort(X)
    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    
    
    ax1.plot(Enr,Yield,'+',color='b',linewidth=2,markersize=3)
    
    ax1.plot(X,ynr_muv(X),color='yellow',linestyle='--',label = 'Mean')
    ax1.plot(X,ynr_mu(X)+3*ynr_sig(X),'k-')
    ax1.plot(X,ynr_mu(X)-3*ynr_sig(X),'k-')

    

    Y = Erer
    Y = np.sort(Y)
    ax1.plot(Erer,Yield_er,'+',color='k',linewidth=2,markersize=3)
    
    ax1.plot(Y,yer_muv(Y),color='yellow',linestyle='--')
    ax1.plot(Y,yer_mu(Y)+3*yer_sig(Y),'k-',label = '3 $\sigma$')
    ax1.plot(Y,yer_mu(Y)-3*yer_sig(Y),'k-')


    ax1.plot()
    ax1.set_ylim(0,1.8)
    ax1.set_xlim(0,200)
    ax1.set_xlabel('$E_R$(keV)',size = '18')
    ax1.set_ylabel('Yield',size = '18')
    ax1.set_title('Bands: Edelweiss Fano', size = '20')
    ax1.tick_params(axis='both', labelsize = '20')
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='-')
    
    ax1.legend(loc=1,prop={'size':12})
    plt.show()