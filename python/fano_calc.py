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
 
  if(aH==0.035):
    filename='EdwYieldWidths-emin{}-n{}-F{}-V{}-eps{}-alpha{}-type{}.h5'.format(emin,n,F,V,epslabel,alpha,rtype)
  else:
    filename='EdwYieldWidths-emin{}-n{}-F{}-V{}-eps{}-alpha{}-aH{}-type{}.h5'.format(emin,n,F,V,epslabel,alpha,aH,rtype)

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
