#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 12:34:40 2018

@author: Mitch
"""
"""
This is the cleaned up file.
"""

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



mu,sigma = norm.fit(b[:,1])
n,bins = np.histogram(b[:,1],50)
y = mlab.normpdf(bins, mu, sigma)

plt.figure()
plt.hist(b[:,1],bins =50,normed = True)
plt.plot(bins,y, 'r--', linewidth = 2)
plt.show()

print("The STD is:", sigma)




