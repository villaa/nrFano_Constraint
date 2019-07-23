import numpy as np
from scipy.special import erf
import math 
from scipy.integrate import quad
import resfuncRead as rfr
import scipy.optimize as so
import pandas as pds
import lmfit as lmf
import nrfano_stats as nfs
import copy





def QEr_Ebin(Q, Ernr, bins=[5, 10, 20, 30, 40, 50, 70,150],silent=False):

    
    #create a dataframe
    nr_df = pds.DataFrame(data={'yield':Q, 'energy':Ernr})

    #bin the data
    nr_df['binned'] = pds.cut(nr_df['energy'],bins)

    #print stats in each bin
    s = nr_df.groupby(pds.cut(nr_df['energy'], bins=bins)).size()
    if not silent:
      print (s)

    #create list of vectors for histogrammin'
    bindf = nr_df.groupby(pds.cut(nr_df['energy'], bins=bins))['yield'].apply(list)
    #print(bindf)

    return bindf 

def QEr_Qhist(bindf, qbins=np.linspace(0,0.6,40)):

 
    xcq = (qbins[:-1] + qbins[1:]) / 2

    qhistos = np.zeros((np.shape(qbins)[0]-1,0))
    qerrs = np.zeros((np.shape(qbins)[0]-1,0))

    #get errors for 0-20 here
    fcerrs = nfs.largestErr_fast()
   
    for i,Qv in enumerate(bindf):
      n,nx = np.histogram(Qv,bins=qbins)
      n = np.reshape(n,(np.shape(n)[0],1))
      qhistos = np.append(qhistos,n,axis=1)
      qerrs0 = np.sqrt(n)
      qerrs0[n<=20] = fcerrs[n[n<=20]]
      qerrs = np.append(qerrs,qerrs0,axis=1)
      #qerrs[n<=20] = fcerrs[n[n<=20]] 
      #use gaussian errors
      #qerrs = np.append(qerrs,np.sqrt(n),axis=1)
      #qerrs[qerrs==0]=1
      #use FC poissonian errors
      #qerrs = np.append(qerrs,nfs.largestErr(n),axis=1)

    return qhistos,qerrs

def QEr_Qfit(qhistos,qerrs, qbins=np.linspace(0,0.6,40),damps=0.1,dmu=1.0,dsig=0.1,silent=False):

    xcq = (qbins[:-1] + qbins[1:]) / 2

    qamps = np.zeros((np.shape(qhistos)[1],))
    qampserrs = np.zeros((np.shape(qhistos)[1],))
    qmus = np.zeros((np.shape(qhistos)[1],))
    qmuerrs = np.zeros((np.shape(qhistos)[1],))
    qsigs = np.zeros((np.shape(qhistos)[1],))
    qsigerrs = np.zeros((np.shape(qhistos)[1],))

    startamps = damps*np.ones((np.shape(qhistos)[1],)) 
    startmus = dmu*np.ones((np.shape(qhistos)[1],)) 
    startsigs = dsig*np.ones((np.shape(qhistos)[1],)) 

    for i,h in enumerate(qhistos[0,:]):
      if not silent:
        print('fitting {}'.format(i))

      qsum = np.sum(qhistos[:,i])
      #do it with lmfit
      params = lmf.Parameters()
      params.add('amp', value=startamps[i])
      params.add('mean', value=startmus[i])
      params.add('sig', value=startsigs[i])
      lmfout = lmf.minimize(gauss_residual, params, args=(xcq, qhistos[:,i]/qsum, qerrs[:,i]/qsum))
      #print(lmf.fit_report(lmfout))
      if not silent:
        print('lmfit results')
        print(lmf.report_fit(lmfout.params))
      qamps[i] = lmfout.params['amp'].value
      qmus[i] = lmfout.params['mean'].value
      qsigs[i] = lmfout.params['sig'].value
 
      #somtimes covariance doesn't exist (if bad fit)
      try: 
        qampserrs[i] = np.sqrt(lmfout.covar[0,0])
        qmuerrs[i] = np.sqrt(lmfout.covar[1,1])
        qsigerrs[i] = np.sqrt(lmfout.covar[2,2])
      except:
        print('bad fit')
        qampserrs[i] = -1 
        qmuerrs[i] = -1 
        qsigerrs[i] = -1 

    if not silent:
      print(qsigs)
      print(qsigerrs)

    return qamps,qampserrs,qmus,qmuerrs,qsigs,qsigerrs

def gauss_residual(params, x, data, eps_data):
    amp = params['amp']
    mean = params['mean']
    sig = params['sig']


    model = amp * np.exp(-(x-mean)**2/(2*sig**2))

    return (data-model) / eps_data

#a bin-centering correction
def bc_corr(E,sig,n=10,Qbar=lambda E: 0.16*E**(0.18)):


    Ev = np.linspace(np.amin(E),np.amax(E),n) 
    Ec = (Ev[:-1] + Ev[1:]) / 2

    #create a dataframe
    df = pds.DataFrame(data={'energy':E})

    s = df.groupby(pds.cut(df['energy'], bins=Ev)).size()[:]
    s = np.asarray(s)
    s = s/np.sum(s)
    #print(s)

    #make a function
    func = {} 
    for i,Er in enumerate(Ec):
      idx = 'f{}'.format(i)
      idx = str(idx)
      prefac = copy.copy(s[i])
      #func[idx] = lambda Q,prefac=prefac,Er=Er: prefac*(1/np.sqrt(2*np.pi*sig(Er)**2))*np.exp(-(Qbar(Er)-Q)**2/(2*sig(Er)**2))
      func[idx] = (lambda prefac,Er: (lambda Q: prefac*(1/np.sqrt(2*np.pi*sig(Er)**2))*np.exp(-(Qbar(Er)-Q)**2/(2*sig(Er)**2))))(prefac,Er)

    fmean = lambda Q: Q*fsum(Q,func)
    f2mom = lambda Q: Q**2*fsum(Q,func)

    #print(quad(fsum,-10,10,args=(func,)))
    mean = quad(fmean,-2,2)[0]
    tmom = quad(f2mom,-2,2)[0]
    #print(mean)
    #print(tmom)

    #print(np.sqrt(tmom-mean**2))

    return (np.sqrt(tmom-mean**2))

def fsum(x,f):

    fsum=0
    for fi in f:
      fsum+=f[fi](x)

    return fsum
