


# Recoil Band Simulation


### Electron Recoils 
 To simulate the bands, we need an initial **simulated** recoil energy., or "true" recoil energy. The energy is found by randomly drawning from a uniform distribution. This energy, $E_{er}$, is used to calculate the true phonon energy $pt$ and the number of electron-hole pairs produced : 

 $$ N_{e/h} = \frac{E_{er}}{\epsilon}$$
Where $\epsilon$ is the average amount of energy to create 1 electron-hole pair. 




To find the measured values for the phonon energy and charge energy, we need to account for the detector resolutions. The charge resolution is $\sigma_Q$  and the phonon resolution is $\sigma_P$.
For detector 1, Run 33.


(From note: High Threshold Analysis, 325, Dan Jardin  )

$$ \sigma_Q = \sqrt{\alpha + \beta E_{er} + \gamma E_{er}}$$
$$\alpha = 0.166, 
\beta = 0.0023,
\gamma = 9.515*10^{-5}\\
$$
$$ \sigma_P = \sqrt{\alpha_p + \beta_p pt + \gamma_p pt^2} $$
$$
\alpha_p = 0.155,
\beta_p = 9.6*10^{-11},
\gamma_p = 0.0005
$$


Similar to how we find $N_{e/h}$, we smear the charge and phonon energy by the detector resolutions to find the measured quantities $Pter$ and $Q_{er}$. 

$$ pt = E_{er} + \frac{V}{\epsilon} N_{e/h}\\
q_{er} = N_{e/h}\epsilon\\$$
$$
Pter = pt + randnormal(0,\sigma_P)\\
Q_{e/r} = q_{er} + randnormal(0,\sigma_Q)$$

The measured recoil energy can now be calculated. 
$$ Er_{er} = Pter + VN_{e/h}$$
We can now calculated the measured Yield, $Y$, for the detector. 
$$ Y = \frac{Q_{er}}{Er_{er}}
$$

Here is a plot of the simulated electron recoil data. The band fits shown are 1,2,3 sigma bands. The fits for the bands 

![Electron Recoil Band](Band_fit/figures/ERer_Band_fits.png)