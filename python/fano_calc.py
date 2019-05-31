import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py
#warnings.resetwarnings()
from scipy.integrate import quad
import resfuncRead as rfr
import scipy.optimize as so
import prob_dist as pd


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

