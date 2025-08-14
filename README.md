# OC-sim
OC-sim is a group of functions designed to generate Observed-Calculated (O-C) diagram data for pulsating stars experiencing circular reflex motion, and inferring a rate of period change by fitting a parabola to the O-C diagram. It first simulates lightcurves of a white dwarf orbitted by a smaller mass object (i.e. a planet), assuming circular orbits with no transit dipping. These lightcurves are fed into [Pyriod](https://github.com/keatonb/Pyriod/tree/master?tab=readme-ov-file), which fits sinusoidal signals to simulated lightcurves, and subsequently return phase data and its corresponding uncertainties. The phase data and uncertainties are used to plot O-C diagrams with a polynomial fit overlayed on phase data, where the rate of period change is calculated from the polynomial fit.



Here is an example output:

<img width="400" height="300" alt="github" src="https://github.com/user-attachments/assets/52c4013b-2a39-4152-81ed-8c880bfb99fa" />





Thank you [@keatonb](https://github.com/keatonb) for contributing to this project!
