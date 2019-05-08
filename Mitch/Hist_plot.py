import numpy as np
import matplotlib.pyplot as plt

def hist_plot(data,prob,prob1,array,bins,fano,recoil_type):
    
    #for x in zip(bin_names):
    
    bin_center = '{0:1.2f}'.format(bins.mid)


    y2,bine = np.histogram(data,bins = 10000, density=True)
    bine1 =  bine#0.5*(bine[1:]+bine[:-1])

    mu = np.mean(data)
    sigma = np.std(data)
    Gaussian = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bine1 - mu)**2 / (2 * sigma**2) )

    
    xmin=min(data)
    xmax=max(data)

    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    
    ax1.hist(data,density=True,histtype = "step",color ='red',label = 'Data',linewidth="3")
#    ax1.hist(data,density=True,label = 'Data',linewidth="3")

    #ax1.plot(bine,Gaussian,linewidth=1, color='r',linestyle = '--',label = 'Gaussian')
    ax1.plot(array,prob,color ='black',label = 'Expected_dep')
    ax1.plot(array,prob1,color ='blue',linestyle = '--',label = 'Expected_indep')


    y,bined = np.histogram(data, normed = True)
    bin_center = 0.5*(bined[1:]+bined[:-1])


    ax1.set_xlabel('Yield',size = '18')
    ax1.set_ylabel('Count',size = '18')
    ax1.set_title('Electron Recoil Yield Distribution' +str(bins) + ' kev, '  'fano = ' +str(fano))
   # ax1.set_xlim(0.7,1.3)

    ax1.set_ylim(bottom = 0.1)
    ax1.set_yscale('log')
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='--')
    ax1.legend()
    ax1.tick_params(axis='both', labelsize = '20')
    #plt.savefig('figures/Dist_fits/nofano'+ str(bin_center)+'.png')
    plt.savefig('/Users/Mitch 1/Desktop/plots/fano = '+str(fano)+'bin = '+str(bins)+'.png')

    plt.show()

    return y,bin_center

def hist_plot_simple(data,label):
    
    
    xmin=min(data)
    xmax=max(data)

    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    
    ax1.hist(data,density=False,histtype = "step",color ='red',label = 'Data',linewidth="3")

    ax1.set_xlabel('Energy [keV]',size = '18')
    ax1.set_ylabel('Count',size = '18')
    ax1.set_title(str(label)) 
    #ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(bottom = 0.1)
    #ax1.set_yscale('log')
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='--')
    ax1.legend(loc=1,prop={'size':18})
    ax1.tick_params(axis='both', labelsize = '20')
    
    #plt.savefig('figures/true'+str(label)+'.png')

    plt.show()

    