import resfuncRead as rfr
import numpy as np 
import pandas as pd 
from band_check import *
from matplotlib.ticker import FuncFormatter
from prob_dist import * 

import sys
sys.path.append('../python')
import prob_dist as PD

from Dist_check import * 
from scipy import integrate 
from Hist_plot import * 
from Stat_analysis import * 
import pylab 
import scipy.stats as stats

    
        
def bin_check(df,s,band_func,bins,cut_idx,expected_v1,expected_v2,Er_true,fano,u):
    
    V = 4
    eps = 0.0033
    Percent1 = []
    Percent2 = []
    Percent3 = []
    Error = [] 
    sig_high = []
    sig_low = []
    expected11 = []
    bin_names = []
    bincenters = []
    
    
   # print("before sort Er_true is:",Er_true)
    recoil_type, upper,lower = band_func(df.E_true,s)
    df1 = pd.DataFrame({'upper':upper,'lower':lower})
    df = pd.concat([df,df1],axis=1)
    df['bin'] = pd.cut(df[cut_idx],bins)

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
            N_mean = np.array(bin_data.N_mean)

            '''For Independent pdf'''
            Sp_v1 = np.array(bin_data.sigp_v1) #Sigma_p
            Sq_v1 = np.array(bin_data.sigq_v1)
             #Sigma_q
            '''For Dependent pdf'''
            Sp_v2 = np.array(bin_data.sigp_v2) #Sigma_p
            Sq_v2 = np.array(bin_data.sigq_v2)            
            SN = np.array(bin_data.sig_N)
             


            bincenters.append(np.mean(E_true))
            bin_center = bin_name.mid

        
            k= (V/eps/1000)
            pdf_v1 = dist_check_v1(u,Ep_mean,Eq_mean,Sp_v1,Sq_v1,k) #amy's defined PDF 
            pdf_v2 = dist_check_v2(u,E_true,N_mean,Sp_v2,Sq_v2,SN)

            #hist_plot(Yield,prob,u,bin_name,fano)  #for new distribution 
            hist_plot(Yield,pdf_v2,pdf_v1,u,bin_name,fano,recoil_type) #for normal Distribution 

            #stat_analysis(Yield,pdf_v2,bin_center)
    
            #g = integrate.quad(lambda x: dist_check(x,Ep_mean,Eq_mean,Sp_mean,Sq_mean,k),np.mean(upper_bound),np.mean(lower_bound) )
            g = integrate.quad(lambda x: dist_check_v2(x,E_true,N_mean,Sp_v2,Sq_v2,SN),np.mean(upper_bound),np.mean(lower_bound))
            H =g[0]*100
            expected11.append(H) #binned (for data points)



            up,down,N = compare(Yield,upper_bound,lower_bound) # up and down are the number of data points OUTSIDE the bands. 
            percent1 = 100*(N - (up+down))/N
            Percent1.append(percent1)


            n = (N-up-down)
            p = n/N
            sig_d = np.sqrt((N*p*(1-p)))/N
            Error.append(sig_d*100) # error in expected containment (data) ss
                                                                                   
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
    print("Bin Spacing (keV)", '\t', "Percent in band",'\t', "Expected",'\t', '\t', "Percent from high", '\t', "Percent from low")     #table column headings
    print("--------------", '\t', '\t' "-------------------", '\t' "-------------------", '\t' "-------------------", '\t' "-------------------")
    bin_spacing = np.array(bin_names).astype(str)

    for x,y,e,h,z,q,t,l in zip(bin_spacing,Percent1,Error,expected11,Percent2,Percent3,sig_high,sig_low):
        print(x, '\t','\t', '{0:1.2f}'.format(y), '\t','±','{0:1.2f}'.format(e),'%', '\t','{0:1.2f}'.format(np.abs(h)), '\t','\t','\t','{0:1.2f}'.format(z), '\t','±','{0:1.2f}'.format(t),'\t', '{0:1.2f}'.format(q), '\t','±','{0:1.2f}'.format(l))
            
    print('--------------------------------------------')
    
    fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    ax1 = axes


    plt.plot(Er_true,expected_v1,label = "Indp_Expected",color ='magenta',linestyle = '--')
    plt.plot(Er_true,expected_v2,label = "Dep_Expected",color ='blue')
    plt.errorbar(bincenters,Percent1,yerr=Error,fmt ='o',label = 'Data', ecolor = 'purple', Linestyle = 'None', capsize=5, capthick=0.5)
    plt.axhline(68, color='r', linestyle='--',Label = "68%")
    ax1.set_xlabel('Recoil Energy [keV]',size = '18')
    ax1.set_ylabel('Percent of Data Contained in Band',size = '18')
    #ax1.set_title('1$\sigma$ Containment Fraction for Electron Recoils Fano =' + str(fano) , size = '15')
    ax1.set_title('1$\sigma$ Containment Fraction for Electron Recoils fano = '+str(fano), size = '15')
    #ax1.set_title('1$\sigma$ Containment Fraction for Nuclear Recoils Edelweiss Fano', size = '15')    
    plt.xticks(bins)
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='-')
    ax1.legend(loc=1,prop={'size':12})
    plt.savefig('/Users/Mitch 1/Desktop/ERContainment_f=5.png')
    #plt.savefig('/Users/Mitch 1/Desktop/NRContainment.png')
    #plt.savefig('Notes/Dist_fits/Eer_Error_Fano =' + str(fano)+'.png')
    plt.show()
    

    return df,bincenters