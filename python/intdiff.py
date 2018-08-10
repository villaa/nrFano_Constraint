#this file gives some functions for computation of Lindhard function 
import numpy as np


#try to make a version of f(t^1/2) because I'll need it
#return a callable
def getft12(version='s2'):

  if version=='s2':
    f = lambda x: 0.343
  elif version=='TF_approx' :
    f = lambda x: 1.43*x**0.35
  else:
    raise ValueError('getft12: invalid request for stopping funtion')

  return f
