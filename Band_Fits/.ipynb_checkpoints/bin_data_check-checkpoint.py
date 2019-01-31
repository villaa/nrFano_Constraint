
import resfuncRead as rfr
import numpy as np
from Yield_calc import * 
import pandas as pd 
from data_check import *
from matplotlib.ticker import FuncFormatter
from prob_dist import * 
from Dist_check import * 
from scipy import integrate 



from tabulate import tabulate

    
        
def bin_check(Yield,Er,Er_true,s,band_func,EP,EQ,sigp,sigq):
    
    V = 4
    eps = 0.0033
    Percent1 = []
    Percent2 = []
    Percent3 = []
    Error = [] 
    sig_high = []
    sig_low = []
    
   # print("before sort Er_true is:",Er_true)
    recoil_type, upper,lower = band_func(Er_true,s)

    Data = np.vstack((Er,Yield,upper,lower,EP,EQ,sigp,sigq,Er_true)).T
    Data1 = Data[np.argsort(Data[:, 0])]
    
    df = pd.DataFrame(Data1)

    #bins  = [10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110] #use 8
    bins = [9,16,25.2,40.3,75.2]

    df1 = pd.cut(df[0],bins, labels = False)

    for q in np.arange(0,4):

    
        E = np.array(df[0].loc[df1 == q])
        x = np.array(df[1].loc[df1 == q]) # Yield 
        y = np.array(df[2].loc[df1 == q]) #upper bound
        z = np.array(df[3].loc[df1 == q]) #lower bound
        P = np.array(df[4].loc[df1 == q]) #EP
        Q = np.array(df[5].loc[df1 == q]) #EQ
        Sp = np.array(df[6].loc[df1 == q]) #Sigma_p
        Sq = np.array(df[7].loc[df1 == q]) #Sigma_q
        Er_true = np.array(df[8].loc[df1 == q])
        
        #print(x)
        
        #look at distributions graphically 
        k= (V/eps/1000)
        u = np.arange(0,2,0.002)
        
        prob = dist_check(u,P,Q,Sp,Sq,k) #amy's defined PDF 
        
        plt.figure()
        plt.plot(u,prob)
        plt.hist(x,normed = True)
        plt.xlim(0.75,1.5)
        plt.show()
        
  
        up,down,N = compare(x,y,z) # up and down are the number of data points OUTSIDE the bands. 

        percent1 = 100*(N - (up+down))/N
        Percent1.append(percent1)
        
        #Error
        p_up = N-up/N 
        #q_up = 1-p_up 
        #Error_up = np.sqrt((N*p_up*q_up)) 
        
        p_down = N-down/N
        
       
        print("Er is:",Er_true)
        print("number of oddball entries",np.sum(Er_true > 12),len(Er_true))
        print("Yield is:",x,"Yield is centered at:",np.mean(x))
        print("lower is: ", z)
        print("Upper is:", y)
        
        
        #v = np.linspace(np.mean(z),np.mean(y),1000)
       # g = np.trapz(prob,v)
        g = integrate.quad(lambda x: dist_check(x,P,Q,Sp,Sq,k),np.mean(z),np.mean(y) )
        print('Area under curve is:',g)
        
        
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
        
        #print(up,down,N)
    
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


