#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'Dependent_Ratio'))
	print(os.getcwd())
except:
	pass

#%%
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as st
from scipy import special as sp


#%%
# Q-Q Plot Definition
def QQplot(data, **kwargs):
    mean = kwargs.get('mean', np.mean(data))
    std = kwargs.get('std', math.sqrt(np.mean((data-np.mean(data))**2)))
    datalbl = kwargs.get('lbl', "Data")
    data.sort()
    Z = (data-mean)/std
    n = len(data)
    L = np.arange(1,n + 1,1) / (n + 1)
    Qnorm = math.sqrt(2)*sp.erfinv(2*L-1)
    MaxZ = max(abs(Z))
    MaxQ = max(abs(Qnorm))
    Max = max(MaxZ, MaxQ)
    plt.scatter(Qnorm, Z, marker = ".")
    plt.plot([-Max, Max],[-Max, Max], color = "black")
    plt.title("Q-Q Plot")
    plt.xlabel("Theoretical Normal Quantiles")
    plt.ylabel(datalbl + " Quantiles")
    plt.show()

#parameters for X and Y
meanX = 5
sdX = 1
meanY = 5
sdY = 1

N = 10000

k = 5000
count = 0
smpl_mean = []

while (count < k):
    np.random.seed((k+1)*count)
    #Z=X/Y
    Xdata = np.random.normal(meanX,sdX,N)
    Ydata = np.random.normal(meanY,sdY,N)
    Zdata = Xdata / Ydata
    smpl_mean.append(np.mean(Zdata))
    count += 1

mean_smpl_mean = np.mean(smpl_mean)
stdev_smpl_mean = np.sqrt(np.mean((smpl_mean-mean_smpl_mean)**2))

print("Mean = ", mean_smpl_mean)
print("stdev = ", stdev_smpl_mean)

density = True
bins = 100
hist, bin_edges = np.histogram(smpl_mean, bins, density=density)

f = (1/(stdev_smpl_mean*np.sqrt(2*np.pi)))*np.exp(-(bin_edges[0:-1]-mean_smpl_mean)*(bin_edges[0:-1]-mean_smpl_mean)/(2*(stdev_smpl_mean*stdev_smpl_mean)))

plt.figure()
plt.plot(bin_edges[0:-1], (hist), drawstyle = "steps-mid", color="red", label = "Simulated Data")
plt.plot(bin_edges[0:-1], f, color="black", label = "Normal Plot")
plt.title("Histogram of Sample Data")
plt.xlabel("Value of Z")
plt.ylabel("Probability")
plt.grid()
plt.legend()
plt.show()

QQplot(smpl_mean)

#st.shapiro(smpl_mean)


#%%



