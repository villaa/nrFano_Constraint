# Electron Recoil Band Distribution 

## Theory 

In the current model, the yield for the electron recoil band is being simulated assuming that the charge energy $E_q$ and the phonon energy $E_p$ are independent.

$$E_p = E_{er} + VN_{eh}$$
Where 
$$N_{eh} = \frac{E_{er}}{\epsilon}$$

To make sure that the measured values for $E_p$ and $E_q$ are generated independently, we randomly sample from distributions for $E_p$ and $E_q$ that were generated using the phonon and charge resolutions from cdms:

$$E_{p_{m}} = random.normal(E_p,\sigma_p)$$
$$E_{q_{m}} = random.normal(E_{er},\sigma_q)$$
Where $E_{p_m}$ and $E_{q_m}$ are the *measured* values.
\
The measured recoil energy is then: 

$$Er_{m} = E_{p_m} - \frac{E_{q_m}V}{\epsilon} $$
 
\
The Yield is then: 

$$Y = \frac{E_{q_m}}{Er_m}$$



## Expected Distribution (CDMS)


When fitting the data from CDMS, the collaboration assumed that the distribution for electron recoils is normally distributed:

<p align = "center" >
    <img src="figures/CDMS_fit.png" width="350" height="300" >
</p>

The bands that CDMS use are generated under this assumption. However, if we look at the containment fractions for the $1\sigma$ CDMS bands with simulated data, we see that the assumption of normality isn't correct:

<p align ="center" >
    <img src="figures/ER_containment1 .png" width ="500" height="200"  >
</p>

One can see that there is two much data in the $1\sigma$ recoil band. If the yield distribution was normal, we would see that the percent of data in the band for each bin is $68\%$.

## Expected Distribution (Ratio of Gaussians)


