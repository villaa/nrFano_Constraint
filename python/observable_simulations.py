import numpy as np
import EdwRes as er
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

def simQEr(label='GGA3',V=4.0,aH=None,C=None,F=0.0,highstats=True):

    #get detector resolution 
    sigHv,sigIv,sigQerv,sigH_NRv,sigI_NRv,sigQnrv = er.getEdw_det_res(label,V, \
    'data/edw_res_data.txt',aH,C) 

    eps = 3.0/1000 #keV per pair, I usually use 3.3 for the numerator, but Edw. uses 3.
    #print(sigQnrv)

    #yield models
    a=0.16
    b=0.18
    Qbar = lambda Er: a*Er**b

    #include a nominal Fano factor
    #F=0.0 #for NRs the factor is probably much higher than this

    if not highstats:
      print('lowstats !')
      f = h5py.File("data/k100_252Cf_shield_Edw_NRs.h5","r")
    else:
      f = h5py.File("data/k100_252Cf_shield_Edw_NRs_large.h5","r")

    #get the data variables
    nr_energies = np.asarray(f['nr_Fano/nr_energies'])
    nr_hits = np.asarray(f['nr_Fano/nr_hits'])
    
    Enr = nr_energies*1000 #initial energies are in MeV
    Enr_ss = nr_energies[nr_hits==1]*1000 #initial energies are in MeV

    Enr_sum = np.sum(Enr,1)
    Enr_ss_sum = np.sum(Enr_ss,1)
    
    #step 1
    EIhit_av = Qbar(Enr)*Enr
    EIhit_av_ss = Qbar(Enr_ss)*Enr_ss
    
    #step 2
    Nhit_av = EIhit_av/eps
    Nhit_av_ss = EIhit_av_ss/eps
    
    #step 3
    Nhit = np.around(np.random.normal(Nhit_av,np.sqrt(F*Nhit_av))).astype(np.float)
    Nhit_ss = np.around(np.random.normal(Nhit_av_ss,np.sqrt(F*Nhit_av_ss))).astype(np.float)
    
    #step 4
    EHhit = (Enr + Nhit*V/1000.0)/(1+(V/(1000*eps)))
    EHhit_ss = (Enr_ss + Nhit_ss*V/1000.0)/(1+(V/(1000*eps)))
    
    #step 5
    EIhit = eps*Nhit
    EIhit_ss = eps*Nhit_ss
    
    #step 6
    EI = np.sum(EIhit,1)
    EI_ss = np.sum(EIhit_ss,1)
    
    #step 7
    EH = np.sum(EHhit,1)
    EH_ss = np.sum(EHhit_ss,1)
    
    #step 8
    if C is not None:
      EI = EI + np.random.normal(0.0,np.sqrt(sigIv(EI) + C**2*Enr_sum**2))
      EI_ss = EI_ss + np.random.normal(0.0,np.sqrt(sigIv(EI_ss)+C**2*Enr_ss_sum**2))
    else:
      EI = EI + np.random.normal(0.0,sigIv(EI))
      EI_ss = EI_ss + np.random.normal(0.0,sigIv(EI_ss))
    
    #step 9
    EH = EH + np.random.normal(0.0,sigHv(EH))
    EH_ss = EH_ss + np.random.normal(0.0,sigHv(EH_ss))
    
    #step 10
    Ernr = (1+(V/(1000*eps)))*EH - (V/(1000*eps))*EI
    Ernr_ss = (1+(V/(1000*eps)))*EH_ss - (V/(1000*eps))*EI_ss
    
    #step 11
    Q = EI/Ernr
    Q_ss = EI_ss/Ernr_ss

    f.close()

    return Q,Ernr,Q_ss,Ernr_ss
