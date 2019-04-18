import resfuncRead as rfr
import numpy as np 
import pandas as pd 
from data_check import *
from matplotlib.ticker import FuncFormatter
from prob_dist import * 
from Dist_check import * 
from scipy import integrate 
from Hist_plot import * 


from tabulate import tabulate

    
        
def bin_check(df,s,band_func,bins,cut_idx,expected,Er_true,fano):
    
    V = 4
    eps = 0.0033
    Percent1 = []
    Percent2 = []
    Percent3 = []
    Error = [] 
    sig_high = []
    sig_low = []
    expected1 = []
    bin_names = []
    bincenters = []
    
    
   # print("before sort Er_true is:",Er_true)
    recoil_type, upper,lower = band_func(df.E_true,s)

    #Data = np.vstack((Er,Yield,upper,lower,EP,EQ,sigp,sigq,Er_true)).T
    #Data1 = Data[np.argsort(Data[:, 0])]
    
    df1 = pd.DataFrame({'upper':upper,'lower':lower})
    
    df = pd.concat([df,df1],axis=1)

    #bins  = [10,13.4,18.1,24.5,33.1,44.8,60.6,80.2,110] #use 8
   # bins = [9,16,25.2,40.3,75.2]

    df['bin'] = pd.cut(df[cut_idx],bins)

    
   # return df1
        


    #for q in np.arange(0,8):
    for bin_name, bin_data in df.groupby('bin'):

        
        
        if len(bin_data) > 1: 
            
  
            bin_names.append(bin_name)
  

            E = np.array(bin_data.E_measured)
            E_true = np.array(bin_data.E_true)
            Yield = np.array(bin_data.Yield) # Yield 
            upper_bound = np.array(bin_data.upper) #upper bound
            lower_bound = np.array(bin_data.lower) #lower bound
            Ep_mean = np.array(bin_data.Ep_mean)
            Eq_mean = np.array(bin_data.Eq_mean) #EQ
            Sp_mean = np.array(bin_data.sigp_mean) #Sigma_p
            Sq_mean = np.array(bin_data.sigq_mean)
            N_mean = np.array(bin_data.N_mean) #Sigma_q
            SN = np.array(bin_data.sig_N)
             

            bincenters.append(np.mean(E_true))
            #print(x)

            #look at distributions graphically f
            k= (V/eps/1000)
            u = np.arange(0,2,0.002) #electron recoils 
           # u = np.linspace(0.1,0.5,1000) #for nuclear recoils. 

            #prob = dist_check(u,Ep_mean,Eq_mean,Sp_mean,Sq_mean,k) #amy's defined PDF 
            prob = dist_check_fano(u,E,N_mean,Sp_mean,Sq_mean,SN)
            hist_plot(Yield,prob,u,bin_name,fano)  
            
        
            up,down,N = compare(Yield,upper_bound,lower_bound) # up and down are the number of data points OUTSIDE the bands. 

            percent1 = 100*(N - (up+down))/N
            Percent1.append(percent1)

        
            p_up = N-up/N 
            p_down = N-down/N


            #g = integrate.quad(lambda x: dist_check(x,Ep_mean,Eq_mean,Sp_mean,Sq_mean,k),np.mean(upper_bound),np.mean(lower_bound) )
            g = integrate.quad(lambda x: dist_check_fano(x,E,N_mean,Sp_mean,Sq_mean,SN),np.mean(upper_bound),np.mean(lower_bound))
            H =g[0]*100
            #H = g*500
           # print('Area under curve is:',g)
            expected1.append(H)


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
    #'10-13.4','13.4-18.1','18.1-24.5','24.5-33.1','33.1-44.8','44.8-60.6','60.6-80.2','80.2-110','110-150'

    print("Bin Spacing (keV)", '\t', "Percent in band",'\t', "Expected",'\t', '\t', "Percent from high", '\t', "Percent from low")     #table column headings
    print("--------------", '\t', '\t' "-------------------", '\t' "-------------------", '\t' "-------------------", '\t' "-------------------")
    bin_spacing = np.array(bin_names).astype(str)

    for x,y,e,h,z,q,t,l in zip(bin_spacing,Percent1,Error,expected1,Percent2,Percent3,sig_high,sig_low):
        print(x, '\t','\t', '{0:1.2f}'.format(y), '\t','±','{0:1.2f}'.format(e),'%', '\t','{0:1.2f}'.format(np.abs(h)), '\t','\t','\t','{0:1.2f}'.format(z), '\t','±','{0:1.2f}'.format(t),'\t', '{0:1.2f}'.format(q), '\t','±','{0:1.2f}'.format(l))
            
    print('--------------------------------------------')
    
    
    #bin = [10.7,25.2,40.3]
    bin_center = [x.mid for x in bin_names]
    print(bin_center)
    
    
    #f = np.linspace(11,96,100)
    
    #print(len(Percent1),len(bin))
    
    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes

    #ax1.scatter(bincenters,Percent1,label = "Simulated")
    #plt.plot(bincenters,expected1,label = "Expected")
    plt.plot(Er_true,expected,label = "Expected",color ='blue')
    plt.errorbar(bincenters,Percent1,yerr=Error,fmt ='o',label = 'Simulated', ecolor = 'purple', Linestyle = 'None', capsize=5, capthick=0.5)
    plt.axhline(68, color='r', linestyle='--',Label = "68%")
    ax1.set_xlabel('Recoil Energy [keV]',size = '18')
    ax1.set_ylabel('Percent of Data Contained in Band',size = '18')
    ax1.set_title('1$\sigma$ Containment Fraction for Electron Recoils Fano =' + str(fano) , size = '15')
    plt.xticks(bins)
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='-')
    ax1.legend(loc=1,prop={'size':12})
    plt.savefig('Notes/Dist_fits/Eer_Error_Fano =' + str(fano)+'.png')
    plt.show()
    

    return df,bincenters