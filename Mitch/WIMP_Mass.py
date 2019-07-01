'''Calculate the mass of a WIMP given its recoil energy '''

import numpy as np 

def WIMP_MASS(ER): 

    V_max = 780000 #[m/s] 
    #M_G = 1.205*10**-25 #Mass of Germanium Nucleus [kg]
    M_G = 72.61 * 1.66*10**-27 #Mass of Germanium Nucleus [kg]


    Er = ER*(1.602*10**-16) #Recoil Energy in Joules [J] from keV! 



    M_X = (M_G*Er)/(V_max*np.sqrt(2*M_G*Er) - Er) #Mass in Kg 


    M_X = M_X*(2.998*10**8)**2*((1/1.62)*10**16) #convert mass from kg to keV

    M_X =M_X*10**-6 #keV to GeV

    return M_X 