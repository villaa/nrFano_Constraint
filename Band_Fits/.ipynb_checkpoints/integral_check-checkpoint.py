import numpy as np


def integrat(f, a, b, N):
    x = np.linspace(a+(b-a)/(2*N), b-(b-a)/(2*N), N)
    fx = f 
    area = np.sum(fx)*(b-a)/N
    return area

