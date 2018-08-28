#this file gives some functions for computation of Lindhard function 
import numpy as np
import dataPython as dp #my text file library
import scipy.interpolate as inter
import scipy.integrate as integrate


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
#try to get the Thomas-Fermi potential in several limits 
#return a callable
def getphi0(version='LT',file=None):

  p0pr = -1.5880464 #seen N-MISC-18-001

  if version=='LT':
    f = lambda x: 1 + p0pr*x + (4/3)*x**(3/2) + (2/5)*p0pr*x**(5/2) \
         + (1/3)*x**3 + (3/70)*p0pr**2*x**(7/2) + (2/15)*p0pr*x**4 \
         + (4/63)*((7/6)-(1/16)*p0pr**3)*x**(9/2)
  elif version=='HT' :
    f = lambda x: 144/x**3 
  elif version=='numeric':
    #get the data
    if(file==None) :
      raise ValueError('getphi0: data file does not exist')
      return None

    data = dp.getXYdata(file)
    #print(data.keys())

    #convert to numpy arrays
    data['xx']= np.asarray(data['xx'])
    data['yy']= np.asarray(data['yy'])

    #print(np.min(data['yy']))

    #spline fit
    f = inter.InterpolatedUnivariateSpline (data['xx'], data['yy'], k=3)
  else:
    raise ValueError('getphi0: invalid request for stopping funtion')

  return f
#try to get the Thomas-Fermi potential's derivative 
#return a callable
def getgradphi0(version='LT',file=None):

  f = getphi0(version,'data/phi0_NACI_format_mod.txt')

  #make a grid of x, and calculate the derivative on the grid
  dx=0.001
  X  = np.arange(0.001,1000,dx)
  y = np.gradient(f(X))

  #spline fit
  fpr = inter.InterpolatedUnivariateSpline (X, y, k=3)

  return fpr
#calculate the function g(xi) see N-MISC-18-002 pg 21
def g(xi=1,version='LT'):

  #get u and derivative
  f = getphi0(version,'data/phi0_NACI_format_mod.txt')
  fpr = getgradphi0(version,'data/phi0_NACI_format_mod.txt')

  #make a callable for integrand
  integrand = lambda x: np.cos(x)*(f(xi/np.cos(x)) - (xi/np.cos(x))*fpr(xi/np.cos(x)))

  #integrate
  result = integrate.quad(integrand, 0.0, np.pi/2.0,epsrel=0.01)

  return result
