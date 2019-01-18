
import resfuncRead as rfr
import numpy as np
from Yield_calc import * 
import pandas as pd 
from data_check import *
from matplotlib.ticker import FuncFormatter



from tabulate import tabulate

    
        
def bin_check(Yield,Er,s,band_func):
    
    
    Percent1 = []
    Percent2 = []
    Percent3 = []
    Error = [] 
    sig_high = []
    sig_low = []
    
    recoil_type, upper,lower = band_func(Er,s)

    Data = np.vstack((Er,Yield,upper,lower)).T
    Data1 = Data[np.argsort(Data[:, 0])]
    
    df = pd.DataFrame(Data1)

    bins  = [10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110]

    df1 = pd.cut(df[0],bins, labels = False)

    for q in np.arange(0,8):

    
     
        x = np.array(df[1].loc[df1 == q])
        y = np.array(df[2].loc[df1 == q])
        z = np.array(df[3].loc[df1 == q])


        up,down,N = compare(x,y,z) # up and down are the number of data points OUTSIDE the bands. 

        percent1 = 100*(N - (up+down))/N
        Percent1.append(percent1)
        
        #Error
        p_up = N-up/N 
        #q_up = 1-p_up 
        #Error_up = np.sqrt((N*p_up*q_up)) 
        
        p_down = N-down/N
        #q_down = 1-p_down
        #Error_down = np.sqrt((N*p_down*q_down)) 
        
        #error = np.sqrt((up/N)**2*(Error_up)**2 + (down/N)**2*(Error_down)**2)
        
        
        
        n = (N-up-down)
        p = n/N
        
        sig_d = np.sqrt((N*p*(1-p)))/N
        Error.append(sig_d*100)
            
        
        
        #looking at symmetry 
        percent2 = 100*(N - (2*up))/N
        percent3 = 100*(N - (2*down))/N
        
        
        #Error for symmetry 
        
        n1 = N-up
        n2 = N-down
        
        sig_up = np.sqrt((2/N)**2*(N*((n1/N)*(1-(n1/N)))))
        sig_down = np.sqrt((2/N)**2*(N*((n2/N)*(1-(n2/N)))))
        
        sig_high.append(100*sig_up)
        sig_low.append(100*sig_down)
        
        
        
        
        Percent2.append(percent2)
        Percent3.append(percent3)
        
        print(up,down,N)
    
    print('--------------------------------------------')  
    print(s,"SIGMA",recoil_type, "RECOIL BAND")
    print('--------------------------------------------')  
    bin_spacing = '10-13.4','13.4-18.1','18.1-24.5','24.5-33.1','33.1-44.8','44.8-60.6','60.6-80.2','80.2-110','110-150'

    print("Bin Spacing (keV)", '\t', "Percent in band", '\t', "Percent from high", '\t', "Percent from low")     #table column headings
    print("--------------", '\t', '\t' "-------------------", '\t' "-------------------", '\t' "-------------------")


    for x,y,e,z,q,t,l in zip(bin_spacing,Percent1,Error,Percent2,Percent3,sig_high,sig_low):
        print(x, '\t','\t', '{0:1.2f}'.format(y), '\t','±','{0:1.2f}'.format(e),'%', '\t', '{0:1.2f}'.format(z), '\t','±','{0:1.2f}'.format(t),'\t', '{0:1.2f}'.format(q), '\t','±','{0:1.2f}'.format(l))
            
    print('--------------------------------------------')
    bin  = [10,13.4,18.1,24.5,33.1,44.8,60.6,80.2]
    
    
    
    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes

    ax1.scatter(bin,Percent1,label = "Percent in Band")
    plt.errorbar(bin,Percent1,yerr=Error,fmt ='o',label = 'Error', ecolor = 'purple', Linestyle = 'None', capsize=5, capthick=0.5)
    plt.axhline(68, color='r', linestyle='--',Label = "68%")
    ax1.set_xlabel('Bins',size = '18')
    ax1.set_ylabel('Percent of Data in Bin',size = '18')
    ax1.set_title('Error in 1$\sigma$ Data Check', size = '20')
    plt.xticks(bin)
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='-')
    ax1.legend(loc=1,prop={'size':12})
    #plt.savefig('figures/ENr_Error_edelweissfano.png')
    plt.show()
    

    return df