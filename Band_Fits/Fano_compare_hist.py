import numpy as np
import matplotlib.pyplot as plt

def hist_plot_compare(data,data1,prob,prob1,array,bins):
    
    #for x in zip(bin_names):
    
    bin_center = '{0:1.2f}'.format(bins.mid)

    
    xmin=min(data)
    xmax=max(data)

    fig,axes = plt.subplots(2,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    ax2 = axes
    
    ax1.hist(data,normed=True,histtype = "step",color ='red',label = 'Data',linewidth="3")
    ax1.plot(array,prob,color ='black',label = 'Expected')
    
    ax2.hist(data1,normed=True,histtype = "step",color ='red',label = 'Data',linewidth="3")
    ax2.plot(array,prob1,color ='black',label = 'Expected')


    ax1.set_xlabel('Yield',size = '18')
    ax1.set_ylabel('Count',size = '18')
    ax1.set_title("Electron Recoil Yield Distribution With Fano " + str(bins) + ' keV')
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(bottom = 0.1)
    ax1.set_yscale('log')
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='--')
    ax1.legend(loc=1,prop={'size':18})
    ax1.tick_params(axis='both', labelsize = '20')
    #plt.savefig('figures/Dist_fits/'+ str(bin_center)+'.png')

    
    ax2.set_xlabel('Yield',size = '18')
    ax2.set_ylabel('Count',size = '18')
    ax2.set_title("Electron Recoil Yield Distribution No Fano " + str(bins) + ' keV')
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(bottom = 0.1)
    ax2.set_yscale('log')
    ax2.grid(True)
    ax2.yaxis.grid(True,which='minor',linestyle='--')
    ax2.legend(loc=1,prop={'size':18})
    ax2.tick_params(axis='both', labelsize = '20')
    #plt.savefig('figures/Dist_fits/'+ str(bin_center)+'.png')

    plt.show()