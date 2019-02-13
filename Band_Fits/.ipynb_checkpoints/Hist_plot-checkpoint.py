import numpy as np
import matplotlib.pyplot as plt

def hist_plot(data,prob,array,bins):
    
    #for x in zip(bin_names):
    
    bin_center = '{0:1.2f}'.format(bins.mid)

    
    xmin=min(data)
    xmax=max(data)

    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes
    
    n_res,n_resx = np.histogram(data,normed=True)
    xresc = (n_resx[:-1] + n_resx[1:]) / 2
    #step, = ax1.step(xresc,n_res,color ='red',label = 'Data')
    ax1.hist(data,normed=True,histtype = "step",color ='red',label = 'Data',linewidth="3")
    ax1.plot(array,prob,color ='black',label = 'Expected')




    ax1.set_xlabel('Yield',size = '18')
    ax1.set_ylabel('Count',size = '18')
    ax1.set_title("Electron Recoil Yield Distribution " + str(bin_center) + ' keV')
    ax1.set_xlim(xmin, xmax)
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='--')
    ax1.legend(loc=1,prop={'size':18})
    ax1.tick_params(axis='both', labelsize = '20')


    plt.show()