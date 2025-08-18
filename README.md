# OC-sim
OC-sim is a group of functions designed to generate Observed-Calculated (O-C) diagram data for pulsating stars experiencing circular reflex motion, and inferring a rate of period change by fitting a parabola to the O-C diagram. It first simulates light curves of a variable star orbitted by a companion, assuming circular orbits and no transits. These lightcurves are fed into [Pyriod](https://github.com/keatonb/Pyriod/), which fits sinusoidal models to the simulated lightcurves, and subsequently returns time delays measured from phase variations (with uncertainties). This represents the O-C data, and a constant rate of period change can be calculated from a parabolic fit.

This code was written as part of the published work "Gravitational Influence from Planets on the Measured Rates of Period Change of Pulsating White Dwarfs" by Ling Xuan Yao, Keaton J. Bell, & Andrew Dublin, 2025 (The Astrophysical Journal, in press).

The main function is `analyze_reflex_OC`

```python
from OC_sim import analyze_reflex_OC
```

It takes many arguments that change the details of the calculation:

```
Args:
    nobs: number of continuous light curves to generate. (integer; default: 50)
    spacingbetween: time spacing between each epoch of continuous observation (astropy time units; default: 28 days)
    duration: length of each individual light curve (astropy time units; default: 1 day)
    cadence: spacing between observations within a light curve (astropy time units; default: 60 seconds)
    centertime: whether to subtract mean of times to center time=0, otherwise first time is 0 (boolean; default: True)
    pulsationperiod: period of the intrinsic sinusoidal signal (with astropy time units; default: 215 s)
    m1: mass of the variable source (with astropy mass units; default: 0.6 solar masses)
    m2: mass of body orbiting (with astropy mass units; default: 1 Jupiter mass)
    separation: distance between variable source and orbiting body (with astropy distance units; default: 27.25 AU)
    inclination: inclination angle of orbit relative to observer (float radians or astropy angle units; default: pi/2)
    orbitalphase: phase of the variable star in orbit at time=0, measured counter-clockwise from the first quadrant 
                  (float radians or astropy angle units; default: 3pi/2... closest to observer)
    pulsationamplitude: relative fractional semi-amplitude of pulsational variability (float; default: 0.05)
    noiselevel: scale of point-by-point Gaussian noise (float; default: 0.01)
    refine_period: whether to redo the analysis with an updated reference period based on initial linear fit (bool; 
                   default: True)
    plot: whether to display a plot of the final O-C diagram (bool; default: True)
Returns:
    A pandas DataFrame with measured rate of period change based on parabolic fit to O-C diagram.
```

Here is an example output from `analyze_reflex_OC(nobs = 580)`:

<img width="400" height="300" alt="github" src="https://github.com/user-attachments/assets/52c4013b-2a39-4152-81ed-8c880bfb99fa" />

```
reference_period    214.99999629964532 s
Pdot                         5.29195e-15
Pdot_err                     3.11952e-17
reducedchi2                  1.00642e+00
```

Thank you [@keatonb](https://github.com/keatonb) for contributing to this project.
