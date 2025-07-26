# OC-sim
OC-sim is a group of functions designed to generate Observed-Calculated (O-C) diagrams for pulsating white dwarfs and rate of period change analysis. It first simulates lightcurves of a white dwarf orbitted by a smaller mass object (i.e. a planet), assuming circular orbits with no transit dipping. These lightcurves are fed into [Pyriod](https://github.com/keatonb/Pyriod/tree/master?tab=readme-ov-file), which fits sinusoidal signals to simulated lightcurves, and subsequently return phase data and its corresponding uncertainties. The phase data and uncertainties are used to plot O-C diagrams with a polynomial fit overlayed on phase data, where the rate of period change is calculated from the polynomial fit.



Thank you @keatonb for contributing to this project
