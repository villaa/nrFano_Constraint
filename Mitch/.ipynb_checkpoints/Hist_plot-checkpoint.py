import numpy as np
import matplotlib.pyplot as plt

def hist_plot(data,prob,array,bins,fano):
    
    #for x in zip(bin_names):
    
    bin_center = '{0:1.2f}'.format(bins.mid)

    
    xmin=min(data)
    xmax=max(data)

    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    
    ax1.hist(data,normed=True,histtype = "step",color ='red',label = 'Data',linewidth="3")
    ax1.plot(array,prob,color ='black',label = 'Expected')



    ax1.set_xlabel('Yield',size = '18')
    ax1.set_ylabel('Count',size = '18')
    ax1.set_title("Electron Recoil Yield Distribution " +str(bins) + ' kev, '  'fano = ' +str(fano))
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(bottom = 0.1)
    ax1.set_yscale('log')
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='--')
    ax1.legend(loc=1,prop={'size':18})
    ax1.tick_params(axis='both', labelsize = '20')
    #plt.savefig('figures/Dist_fits/nofano'+ str(bin_center)+'.png')
    plt.savefig('figures/Dist_fits/'+ str(bin_center)+'.png')


    plt.show()