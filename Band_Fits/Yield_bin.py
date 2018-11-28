
import resfuncRead as rfr
import numpy as np
from Yield_calc import * 
import pandas as pd 
from data_check import *

import numpy as np
import matplotlib.pyplot as plt
import resfuncRead as rfr
import pandas as pd
from scipy.stats import norm
from scipy import stats
import matplotlib.mlab as mlab
plt.rcParams.update({'font.size': 30})
import seaborn as sns




def Hist(data):
    
    
    df = data

    bins  = [10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110]

    df1 = pd.cut(df[0],bins, labels = False)
    
    


    for i in np.arange(0,8):
        
        x = np.array(df[1].loc[df1 == i])
        print('Amount of data in '+ str(bins[i])+ ' keV each bin is:',len(x))
        
        n_res,n_resx = np.histogram(x,60)
        bin = (n_resx[:-1] + n_resx[1:]) / 2
        
        fig,ax1 = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
        
        
        #sns.set_style('darkgrid')
        #sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8})
        
        sns.distplot(x,bins = 50,color="blue",fit=stats.norm, kde=False,hist_kws={"histtype": "step", "linewidth": 3,"alpha": 1, "color": "r"})
        #sns.distplot(x,color="blue",fit=stats.gennorm,norm_hist = True, kde=False,)


            
        
        ax1.set_yscale('log')
        ax1.set_xlabel('Yield',size = '18')
        #ax1.set_ylabel('Counts',size = '1')
        ax1.set_title('Yield for the ' + str(bins[i])+ 'keV Bin',size ='20')
        ax1.grid(True)
        ax1.yaxis.grid(True,which='minor',linestyle='--')
        #ax1.legend(loc=1,prop={'size':18})
        ax1.tick_params(axis='both', labelsize = '20')
        
        plt.show()
        
     
        
        
        

        
     
    
        '''   
         x = np.array(df[1].loc[df1 == i])
    
        xmin=min(x)
        xmax=max(x)
        amax= 30
        #print(np.std(val_stat))
        #print(xmin,xmax)
        
        
        #xt = plt.xticks()[0]    
        #lnspc = np.linspace(xmin, xmax, len(x))  
        
        
        n_res,n_resx = np.histogram(x,300,range=(xmin,xmax))
        xresc = (n_resx[:-1] + n_resx[1:]) / 2
        
        m, s = stats.norm.fit(x) # get mean and standard deviation  
        pdf_g = stats.norm.pdf(xresc, m, s) # get theoretical values in our interval

   

        #set up a 1d plot
        fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
        ax1 = axes

        #X = np.arange(-amax,amax,0.01) 
        #step, = ax1.step(xresc,n_res, where='mid',color='r', linestyle='-', label=Yield, linewidth=3)
        ax1.step(xresc,n_res,where='mid',color='r', linestyle='-')
        ax1.plot(xresc, pdf_g, label="Norm") #plot of fit
        ymin=0
        ymax=max(n_res)


        ax1.set_yscale('log')
        #ax1.set_xlim(xmin, xmax) #in pairs
        ax1.set_ylim(ymin,700)
        ax1.set_xlabel('Yield',size = '18')
        ax1.set_ylabel('Counts',size = '1')
        ax1.set_title("Yield ",size ='20')
        ax1.grid(True)
        ax1.yaxis.grid(True,which='minor',linestyle='--')
        #ax1.legend(loc=1,prop={'size':18})
        ax1.tick_params(axis='both', labelsize = '20')
        
        #plt.plot(lnspc, pdf_g, label="Norm") #plot of fit 
        #ax1.plot(lnspc, pdf_g, label="Norm") #plot of fit
        
        for axis in ['top','bottom','left','right']:
          ax1.spines[axis].set_linewidth(2)

        plt.tight_layout()
        plt.show()'''
