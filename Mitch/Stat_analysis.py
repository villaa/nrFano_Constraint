from scipy import stats
import matplotlib.pyplot as plt
import pylab
import numpy as np 


def stat_analysis(Yield,prob,bin_center):

    plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
    stats.probplot(Yield, dist="norm", plot=pylab)
    pylab.title('Q-Q plot '+ str(bin_center)+' keV')
    plt.grid(True)
    plt.savefig('figures/Q-Q_plot '+str(bin_center)+' keV.png')
    pylab.show()

    G = stats.skew(Yield)
    K = stats.kurtosis(Yield,fisher=False)
    W = stats.shapiro(Yield)


    percent_diff = abs(K-3)/(3)

    #chisq = np.sum((Yield-prob)**2)/np.sqrt(len(Yield))

    



    print('skew is: ',G)
    print('kurtosis is: ', K )
    print('Shapero-Wilk results: ', W )
    print('diff = ', percent_diff)
    #print("chi-Squared = ", chisq)