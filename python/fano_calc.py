import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py
#warnings.resetwarnings()
from scipy.integrate import quad
import resfuncRead as rfr
import scipy.optimize as so
import prob_dist as pd
import os


def writeFano(file='fanoout.h5'):


  Ef = np.linspace(10,100,120)
  F = np.zeros(np.shape(Ef))

  for i,E in enumerate(Ef):
    print('calculating {} out of {}'.format(i+1,np.shape(Ef)[0]))
    lowsig = pd.sigroot(0.001,E)

    findF = lambda F,Er,C: pd.sigroot(F,Er)**2 - lowsig**2 -C**2

    F[i] = so.brentq(findF,0.001,200,rtol=0.001,maxiter=100,args=(E,0.035,))
    print('at energy E = {} keV NRFano is {}'.format(E,F[i]))

  of = h5py.File(file, 'w')
  d = Ef 
  #hits dataset
  dset_hits = of.create_dataset("nr_Fano_Extracted/nr_energies", np.shape(d), dtype=np.dtype('float64').type, compression="gzip", compression_opts=9)
  dset_hits[...] = d
  d = F 
  #hits dataset
  dset_hits = of.create_dataset("nr_Fano_Extracted/fano", np.shape(d), dtype=np.dtype('float64').type, compression="gzip", compression_opts=9)
  dset_hits[...] = d

  of.close()
  return

def calcQWidth(n,F=10,V=4,eps=(3/1000),alpha=(1/100),Qbar=lambda x: 0.16*x**0.18,aH=0.035,path='./'):

  Er = np.linspace(7,100,n)
  emin = np.min(Er)
  emax = np.max(Er)
  #n = np.shape(Er)[0]
  epslabel = eps*1000
  rtype = ''
  if(Qbar(10)<0.8):
    rtype='nr'
  else:
    rtype='er'

  emins = '{:01.1f}'.format(emin)
  ns = '{:04.0f}'.format(n)
  Fs = '{:03.0f}'.format(F)
  Vs = '{:2.1f}'.format(V)
  epss = '{:1.1f}'.format(epslabel)
  alphas = '{:1.3f}'.format(alpha)
  aHs = '{:1.3f}'.format(aH)
 
  if(aH==0.035):
    filename='EdwYieldWidths-emin{}-n{}-F{}-V{}-eps{}-alpha{}-type{}.h5'.format(emins,ns,Fs,Vs,epss,alphas,rtype)
  else:
    filename='EdwYieldWidths-emin{}-n{}-F{}-V{}-eps{}-alpha{}-aH{}-type{}.h5'.format(emins,ns,Fs,Vs,epss,alphas,aHs,rtype)

  out=[]
  if(os.path.exists('{}data/{}'.format(path,filename))):
    #just open it and return the array 
    f = h5py.File('{}data/{}'.format(path,filename),"r")
    out = np.asarray(f['sigma'])
    

  else:
    #have to compute everything and store the result
    out = np.zeros(np.shape(Er))
    for i,E in enumerate(Er):
      print('calculating {} of {} points (for filename {})'.format(i+1,n,filename))
      out[i] = pd.sigrootEdw(F,E,V,eps,alpha,Qbar,aH) 

    f = h5py.File('{}data/{}'.format(path,filename),"w")
    dset = f.create_dataset('sigma',np.shape(out),dtype=np.dtype('float64').type,compression="gzip",compression_opts=9)
    dset[...] = out
    dset = f.create_dataset('Er',np.shape(out),dtype=np.dtype('float64').type,compression="gzip",compression_opts=9)
    dset[...] = Er 

  #f.close() 
  return (out,Er)

def RWCalc(filename='test.h5',det='GGA3',band='ER',F=0.00001,V=4.0,alpha=(1/10000.0),aH=0.035,Erv=None,sigv=None,erase=False):

  #n=10
  #Er = np.linspace(7,100,n)
  #emin = np.min(Er)
  #emax = np.max(Er)

  #emins = '{:01.1f}'.format(emin)
  #ns = '{:04.0f}'.format(n)
  Fs = '{:03.0f}'.format(F)
  Vs = '{:2.1f}'.format(V)
  alphas = '{:1.3E}'.format(alpha)
  aHs = '{:1.3f}'.format(aH)
 
  path='{}/{}/{}/{}/{}/{}/'.format(det,band,Vs,alphas,aHs,Fs)

  print(path)

  #check for path
  f = h5py.File(filename,'a')
  exEr = path+'Er' in f
  exsig = path+'sig' in f
  print(exEr)
 

  #make some vector
  if exEr&exsig&~erase:
    Er = np.asarray(f[path+'Er'])
    sig = np.asarray(f[path+'sig'])
  else:
    Er = np.zeros((0,))
    sig = np.zeros((0,))

  #add in the data supplied
  if (Erv is not None)&(sigv is not None):
    Er = np.append(Er,Erv)
    sig = np.append(sig,sigv)

  if exEr&exsig:
    del f[path+'Er']
    del f[path+'sig']

  #sort the array
  idxEr = np.argsort(Er)
  Er = Er[idxEr]
  sig = sig[idxEr]

  Er,uidx = np.unique(Er,return_index=True)
  sig = sig[uidx]

  dset = f.create_dataset(path+'Er',np.shape(Er),dtype=np.dtype('float64').type, \
  compression="gzip",compression_opts=9)
  dset[...] = Er
  dset = f.create_dataset(path+'sig',np.shape(Er),dtype=np.dtype('float64').type, \
  compression="gzip",compression_opts=9)
  dset[...] = sig

  f.close()

  return (Er,sig)

def storeQWidth(n,filename='test.h5',det='GGA3',band='ER',F=0.00001,V=4.0,alpha=(1/10000.0),aH=0.035,erase=False):

  Er = np.linspace(7,100,n)
  emin = np.min(Er)
  emax = np.max(Er)

  Fs = '{:03.0f}'.format(F)
  Vs = '{:2.1f}'.format(V)
  alphas = '{:1.3E}'.format(alpha)
  aHs = '{:1.3f}'.format(aH)

  (Er_stored,sig_stored) = RWCalc(filename,det,band,F,V,alpha,aH)
  n_stored = np.shape(Er_stored)[0]

  print(Er)
  print(Er_stored)
  print(sig_stored)

  #calculate density and overlap
  if n_stored>0:
    emin_stored = np.min(Er_stored)
    emax_stored = np.max(Er_stored)
    ovr = (emax_stored-emin_stored)/(emax-emin)
  else:
    emin_stored = 0 
    emax_stored = 0 
    ovr = 0

  if ((emax_stored-emin_stored)>0)&((emax-emin)>0):
    den = (n_stored/(emax_stored-emin_stored))/(n/(emax-emin))
  else: 
    den = 0

  print(ovr)
  print(den)

  if erase:
    E_needed = Er
  else:
    idx_needed = ~np.isin(Er,Er_stored)
    E_needed = Er[idx_needed]

  print(E_needed)

  sigcalc = np.zeros(np.shape(E_needed))
  for i,E in enumerate(E_needed):
    print('Calculating with sigmomEdw for E = {:3.2f} keV'.format(E))
    sigcalc[i] = pd.sigmomEdw(E,band=band,label=det,F=F,V=V,aH=aH,alpha=alpha)
    print(sigcalc[i])
     
  print(E_needed)
  print(sigcalc)
  (Er_new,sig_new) = RWCalc(filename,det,band,F,V,alpha,aH,Erv=E_needed,sigv=sigcalc,erase=erase)
  return (Er_new,sig_new)
