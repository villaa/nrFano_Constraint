#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 22:30:08 2018

@author: Mitch
"""

import numpy as np
import matplotlib.pyplot as plt

d
k = 0.157 # from lindhard
q = 1.602*10**-19 #electron charge 
V = 2.0 # voltage Bias
sigt = 0.03
sigq = 0.01 
eps = 0.036

x =np.arange(0,200) # keV 

Er = 93*np.exp(-x/29) #For NR
#Er = np.random.normal(0, 0.3, 200)

ER = []
Yield = []
Ran = []

for i in np.arange(0,200):
    
    E = np.random.choice(Er)
    Ran.append(E)
    

    Y = 0.3

    EQ = Y*Er 
    Et = (EQ*q*V)/eps + Er 
    

    Fq = 1/np.sqrt(2*3.14*sigq**2)*np.exp(-(x-0.3*Er)**2/2*sigq**2)# equation 5.23
   # FQ = Fq.n
    
    Ft = 1/np.sqrt(2*3.14*sigt*82)*np.exp(-(x-(1+0.3*(q*V/0.036))*Er)**2/2*sigt**2) #equation 5.2


    qimean = EQ + Fq 

    ptnf = Et + Ft 

    ER1 = ptnf - qimean*(q*V/eps)
    ER.append(ER1)

    yield1 = qimean / ER1 
    Yield.append(yield1)

plt.scatter(ER,Yield)
plt.ylim(0,0.5)


