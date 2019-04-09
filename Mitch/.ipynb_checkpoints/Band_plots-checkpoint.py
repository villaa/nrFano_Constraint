import numpy as np 
import matplotlib.pyplot as plt
import resfuncRead as rfr
from scipy.stats import gaussian_kde

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


def NR_band_plot(df):
    
    Enr = df.E_measured
    a = df.Yield
    X = df.E_measured
    X = np.sort(X)
    
    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    xy = np.vstack([Enr,a])
    z = gaussian_kde(xy)(xy)
    
    
    #ax1.plot(Enr,a,'+',color='b',linewidth=2,markersize=3)
    cm = plt.cm.get_cmap('nipy_spectral')
    cax = ax1.scatter(Enr, a, c=z, s=10, cmap = cm)
    
    ax1.plot(X,ynr_muv(X),color='yellow',linestyle='--',label = 'Mean')
    ax1.plot(X,ynr_mu(X)+1*ynr_sig(X),'r-',label = '1 $\sigma$')
    ax1.plot(X,ynr_mu(X)-1*ynr_sig(X),'r-')
    
    

    ax1.plot(X,ynr_mu(X)+2*ynr_sig(X),'m-',label = '2 $\sigma$')
    ax1.plot(X,ynr_mu(X)-2*ynr_sig(X),'m-')

    ax1.plot(X,ynr_mu(X)+3*ynr_sig(X),'k-',label = '3 $\sigma$')
    ax1.plot(X,ynr_mu(X)-3*ynr_sig(X),'k-')

    ax1.plot()
    ax1.set_ylim(0.1,.5)
    ax1.set_xlim(0,160)
    
    
    ax1.set_xlabel('$E_R$(keV)',size = '18')
    ax1.set_ylabel('Yield',size = '18')
    ax1.set_title('Data Check NR Band C = 0.04', size = '20')
    ax1.tick_params(axis='both', labelsize = '20')
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='-')
    

    ax1.legend(loc=1,prop={'size':12})
    plt.savefig('figures/ENr_band_edelweissfano1.png')
    plt.show()
    
    
def ER_band_plot(Erer,b,N):
    import time
    start_time = time.time()
    
    Y = Erer
    Y = np.sort(Y)
    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    xy = np.vstack([Erer,b])
    z = gaussian_kde(xy)(xy)

    #ax1.plot(Erer,b,'+',color='b',linewidth=2,markersize=3,alpha = 0.5)
    cm = plt.cm.get_cmap('nipy_spectral')
    cax = ax1.scatter(Erer, b, c=z, s=10, cmap = cm)
    #ax1.hist2d(a, b, (50, 50), cmap=plt.cm.jet)
    
    
     
    '''   
    ax1.plot(Y,yer_muv(Y),color='blue',linestyle='--',label = 'Mean')
    
    ax1.plot(Y,yer_mu(Y)+1*yer_sig(Y),'r-',label = '1 $\sigma$')
    ax1.plot(Y,yer_mu(Y)-1*yer_sig(Y),'r-')

    ax1.plot(Y,yer_mu(Y)+2*yer_sig(Y),'m-',label = '2 $\sigma$')
    ax1.plot(Y,yer_mu(Y)-2*yer_sig(Y),'m-')

    ax1.plot(Y,yer_mu(Y)+3*yer_sig(Y),'k-',label = '3 $\sigma$')
    ax1.plot(Y,yer_mu(Y)-3*yer_sig(Y),'k-')
    '''
    
    
    
    
    x = np.ones(len(Y))
    
    ax1.plot(Y,x,color='blue',linestyle='--',label = 'Mean')
    
    ax1.plot(Y,1+1*yer_sig(Y),'r-',label = '1 $\sigma$')
    ax1.plot(Y,1-1*yer_sig(Y),'r-')

    ax1.plot(Y,1+2*yer_sig(Y),'m-',label = '2 $\sigma$')
    ax1.plot(Y,1-2*yer_sig(Y),'m-')

    ax1.plot(Y,1+3*yer_sig(Y),'k-',label = '3 $\sigma$')
    ax1.plot(Y,1-3*yer_sig(Y),'k-')
    
    
    
    
    
    
    
    
    
    
    plt.axvline(10, color='r', linestyle='-')

    ax1.plot()
    ax1.set_ylim(0.7,1.3)
    ax1.set_xlim(0,165)
    
    
    ax1.set_xlabel('$E_R$(keV)',size = '18')
    ax1.set_ylabel('Yield',size = '18')
    ax1.set_title('Data Check ER Band', size = '20')
    ax1.tick_params(axis='both', labelsize = '20')
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='-')
    

    ax1.legend(loc=1,prop={'size':12})
    #fig.colorbar(cax)
    #plt.savefig('figures/ERer_Band_')
    plt.show()
    
    print("My program took", (time.time() - start_time)/60, " mintues to run")
    

def NR_ER_plot(Enr,Erer,a,b):
    
    
    # Nuclear Recoil Band 
    X = Enr
    X = np.sort(X)
    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    
    
    ax1.plot(Enr,a,'+',color='b',linewidth=2,markersize=3)
    
    ax1.plot(X,ynr_muv(X),color='yellow',linestyle='--',label = 'Mean')
    ax1.plot(X,ynr_mu(X)+1*ynr_sig(X),'r-',label = '1 $\sigma$')
    ax1.plot(X,ynr_mu(X)-1*ynr_sig(X),'r-')
    
    

    ax1.plot(X,ynr_mu(X)+2*ynr_sig(X),'m-',label = '2 $\sigma$')
    ax1.plot(X,ynr_mu(X)-2*ynr_sig(X),'m-')

    ax1.plot(X,ynr_mu(X)+3*ynr_sig(X),'k-',label = '3 $\sigma$')
    ax1.plot(X,ynr_mu(X)-3*ynr_sig(X),'k-')

    
    # Electron recoil band 
    Y = Erer
    Y = np.sort(Y)

    ax1.plot(Erer,b,'+',color='k',linewidth=2,markersize=3)
    
    ax1.plot(Y,yer_muv(Y),color='yellow',linestyle='--')
    
    ax1.plot(Y,yer_mu(Y)+1*yer_sig(Y),'r-')
    ax1.plot(Y,yer_mu(Y)-1*yer_sig(Y),'r-')

    ax1.plot(Y,yer_mu(Y)+2*yer_sig(Y),'m-')
    ax1.plot(Y,yer_mu(Y)-2*yer_sig(Y),'m-')

    ax1.plot(Y,yer_mu(Y)+3*yer_sig(Y),'k-')
    ax1.plot(Y,yer_mu(Y)-3*yer_sig(Y),'k-')
    
    
    ax1.set_xlabel('$E_R$(keV)',size = '18')
    ax1.set_ylabel('Yield',size = '18')
    ax1.set_title('NR_ER Yield Plot', size = '20')
    ax1.tick_params(axis='both', labelsize = '20')
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='-')
    ax1.set_xlim(0,150)
    ax1.set_ylim(-0.3,1.5)
    ax1.legend(loc=4,prop={'size':12})
    
    plt.savefig('figures/NRER_Band_fits_fano_shifted.png')
    plt.show()
    
    
def EP_EQ_plot(EP,EQ):
    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    xy = np.vstack([EP,EQ])
    z = gaussian_kde(xy)(xy)

    #ax1.plot(Erer,b,'+',color='b',linewidth=2,markersize=3,alpha = 0.5)
    cm = plt.cm.get_cmap('nipy_spectral')
    cax = ax1.scatter(EP, EQ, c=z, s=10, cmap = cm)

    ax1.set_xlabel('$E_P$ (keV)',size = '18')
    ax1.set_ylabel('$E_Q$ (keV)',size = '18')
    ax1.set_title('EQ - EP Space', size = '20')
    ax1.tick_params(axis='both', labelsize = '20')
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='-')
    #fig.colorbar()
    
    plt.savefig('figures/EQ_EP_space')
    