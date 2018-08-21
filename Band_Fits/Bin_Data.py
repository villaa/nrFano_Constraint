import numpy as np
import matplotlib.pyplot as plt

import matplotlib.mlab as mlab
import math 
from Read_Data import Read_File
from scipy import stats
import pandas as pd
import matplotlib.mlab as mlab
from scipy.stats import norm
from scipy.optimize import curve_fit

Enr,Yield = Read_File("test_data.txt") #reads in data file. 

#Enr = np.arange(0,100,10)
#Yield = np.arange(0,0.4,0.04)

combined = np.vstack((Enr, Yield)).T\

#print(combined)


a,b,c,d,e,f,g = np.array_split(combined,7) # about 50keV each 

#print(test)


df = pd.DataFrame(data=combined)

#print(df)

bins  = [0,50,100,150,200,250,300,400]

df1 = pd.cut(df[0],bins)


print(df1)

hist = np.histogram2d(Enr,Yield,bins = [0,50,100,150,200,250,300,400])


df2 = pd.DataFrame.as_matrix(df1)

#print(df2)
'''
plt.hist(test[0])
plt.figure()
plt.hist(test[1])
plt.show()
plt.figure()
plt.hist(test[2])
plt.show()



test1 = test[0]
test2 = test[1]
test3 = test[2]
test4 = test[3]
test5 = test[4]
test6 = test[5]
test7 = test[6]




n,bins = np.histogram(b[:,1],425)

def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

popt, pcov = curve_fit(gaussian, b[:,1],bins)
    
plt.figure()
plt.hist(b[:,1],bins = 7)
plt.plot(b[:,1],gaussian(b[:,1], *popt))

plt.show()


mu,sigma = norm.fit(b[:,1])

n,bins = np.histogram(b[:,1],425)

y = mlab.normpdf(bins, mu, sigma)

plt.figure()
plt.hist(b[:,1],bins =7,normed = True)
plt.plot(bins,y, 'r--', linewidth = 2)
plt.show()
'''



