import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.options.display.float_format = '{:.5e}'.format # Display small number with scientific notation
import lightkurve as lk
from Pyriod import Pyriod
from astropy import constants as c, units as u
from astropy.coordinates import Angle

def timefunction(nobs = 50, spacingbetween = 28*u.day, duration = 1*u.day, cadence = 60*u.s):
    """Generate timestamps for a series of simulated light curves.

    Simulates a total of "nobs" individual light curves separated in time by "spacingbetween".
    Each light curve has a duration of "duration" with individual data points sampled every "cadence".
    Default values will produce timestamps for 50 light curves, separated by 28 days, each sampled
    every 60 seconds for a durations of 1 day.
    
    Args:
        nobs: number of continuous light curves to generate. (integer; default: 50)
        spacingbetween: time spacing between each epoch of continuous observation (astropy time units; default: 28 days)
        duration: length of each individual light curve (astropy time units; default: 1 day)
        cadence: spacing between observations within a light curve (astropy time units; default: 60 seconds)
    
    Returns:
        A list of time arrays for each individual light curve (with astropy units of days)
    """
    # Each light curve has same duration and cadence
    eachlc = np.arange(int((duration/cadence).decompose().value))*cadence
    # Repeat light curves nobs times every spacingbetween
    times = [(i*spacingbetween+eachlc).to(u.d) for i in range(nobs)]
    return times

def fluxfunction(times,pulsationperiod = 215*u.s, m1=0.6*c.M_sun, m2 = c.M_jup, separation = 27.25*c.au, 
                 inclination = np.pi/2, orbitalphase = 3*np.pi/2, centerphase = True,
                 pulsationamplitude = 0.05, noiselevel = 0.01): 
    """Generate light curve flux values for a pulsating star with reflex motion at given times.

    For a list of time arrays (individual light curves with astropy units), simulate data for an intrinsically coherent 
    sinusoidal variable experiencing time modulaton from circular reflex motion about a center of mass.
    
    Args:
        times: list of timestamp arrays for each light curve to be generated (list of arrays with astropy time units)
        pulsationperiod: period of the intrinsic sinusoidal signal (with astropy time units; default: 215 s)
        m1: mass of the variable source (with astropy mass units; default: 0.6 solar masses)
        m2: mass of body orbiting (with astropy mass units; default: 1 Jupiter mass)
        separation: distance between variable source and orbiting body (with astropy distance units; default: 27.25 AU)
        inclination: inclination angle of orbit relative to observer (float radians or astropy angle units; default: pi/2)
        orbitalphase: phase of the variable star in orbit at a specified time (see centerphase), measured counter-clockwise 
                      from the first quadrant (float radians or astropy angle units; default: 3pi/2... closest to observer)
        centerphase: if True, orbitalphase is phase at average timestep, otherwise orbitalphase is phase at time t=0 (bool; default: True)
        pulsationamplitude: relative fractional semi-amplitude of pulsational variability (float; default: 0.05)
        noiselevel: scale of point-by-point Gaussian noise (float; default: 0.01)
    Returns:
        A list of lightkurve.LightCurve objects for each individual light curve.
    """
    if type(inclination) is not u.quantity.Quantity:
        inclination *= u.rad
    if type(orbitalphase) is not u.quantity.Quantity:
        orbitalphase *= u.rad
    # pulsation angular frequency
    omega1 = (2*np.pi/pulsationperiod)*u.rad
    # Calculate orbital period from Kepler's 3rd law
    P = np.sqrt((4*np.pi**2)*((separation)**3)/(c.G*(m1+m2)))  
    omega2 = (2*np.pi/P)*u.rad #orbital angular frequency
    # Compute distance of m1 from center of mass
    com = separation*m2/(m1+m2)
    # in units of light-travel-time
    tprime = (com/c.c)*np.sin(inclination)
    # Is the orbital phase at the observation midpoint?
    referencetime = 0
    if centerphase:
        referencetime = np.mean(times*u.d)
    # Make empty list to store light curves
    lightcurves = []
    # Loop over arrays of timestamps and generate lightcurves for each
    for timestamps in times:
        reflex = tprime*np.sin(omega2*(timestamps-referencetime)+orbitalphase) #time delays
        noise = noiselevel*np.random.randn(len(timestamps)) # Gaussian measurement noise
        flux = pulsationamplitude*np.sin(omega1*(timestamps-reflex)+np.pi*u.rad) + noise + 1 # Normalized to 1
        lightcurves.append(lk.LightCurve(time=timestamps, flux = flux, flux_err=noiselevel).normalize())
    return lightcurves

def OC(lightcurves, referenceperiod):
    """Measure Observed-Calculated of a sinusoidal signal across a series of light curves.
    
    Measure variations in signal arrival time (Observed) compared to a constant reference period (Calculated)
    using phase variations measured with the Pyriod frequency analysis code. One timing data point is collected
    for each light curve provided in a list.
    
    Args:
        lightcurves: list of lightkurve.LightCurve objects, for O-C to be measured from each.
        referenceperiod: reference period for Calculated times (float in units of seconds, or with astropy units)

    Returns:
        pandas DataFrame with average epoch of each light curve (average time divided by reference period), O-C 
        time delay, and O-C uncertainty.

    Notes:
        A correction is attempted when the phase of the signal crosses the 0-2pi boundary. Assumes that signal
        coverage is frequent enought that large jumps in O-C phase (> pi) are not expected.
    """
    # convert reference period to frequency in units of microHz
    if type(referenceperiod) is not u.quantity.Quantity:
        referenceperiod *= u.s
    referencefrequency = (1/referenceperiod).to(u.uHz) # microHz for Pyriod

    # Prepare to collect results
    nlightcurves = len(lightcurves) # number of light curves in list
    meantime = np.zeros(nlightcurves) # array to store phase variations
    phase = np.zeros(nlightcurves) # array to store phase variations
    phaseerr = np.zeros(nlightcurves) # array to store phase errors
    for i in range(nlightcurves): 
        # Make Pyriod instance to measure phase variations
        # Disable GUI functionality (unneeded and slow)
        meantime[i] = np.mean(lightcurves[i].time.value)
        pyriod = Pyriod(lightcurves[i],gui=False,rescale_covar=True, freq_unit = "muHz")  
        pyriod.add_signal(referencefrequency.value,fixfreq=True) # Do not vary frequency in the fit
        pyriod.fit_model() # Fitting sinusoidal model to time series data
        phase[i] = pyriod.fitvalues['phase'].values[0] # record best-fit phase
        phaseerr[i] = pyriod.fitvalues['phaseerr'].values[0] # record uncertainty
    # Format results of O-C
    oc = -np.array(phase)*referenceperiod # Unprocessed O-C signal
    ocerr = np.array(phaseerr)*referenceperiod
    # Apply phase boundary correction
    ind = np.diff(oc) > referenceperiod/2 # difference between consecutive cells are greater than half the pulsation period
    oc[1:] = oc[1:]-np.cumsum(ind)*referenceperiod # Subtracts cycles where previous condition (OC>pulsationperiod/2) is fulfilled
    oc -= oc[0] #start at O-C=0
    epoch = (meantime*u.d/referenceperiod).decompose() # Converting time to cycles
    return pd.DataFrame({"Epoch":epoch,"OC": oc,"OC_err": ocerr})

def analyze_reflex_OC(nobs = 50, spacingbetween = 28*u.day, duration = 1*u.day, cadence = 60*u.s,
                      pulsationperiod = 215*u.s, m1=0.6*c.M_sun, m2 = c.M_jup, 
                      separation = 27.25*c.au, inclination = np.pi/2,  orbitalphase = 3*np.pi/2, centerphase = True,
                      pulsationamplitude = 0.05, noiselevel = 0.01, refine_period = True, plot = True):
    """Full O-C analysis of time series data for a variable star experiencing reflex motion.

    Fits a parabola to the data and report a P-dot (constant rate of period change) even though the similated data
    has sinusoidal O-C variations. Code developed for the project "Gravitational Influence from Planets on the Measured 
    Rates of Period Change of Pulsating White Dwarfs" by Yao, Bell, & Dublin (submitted).

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
    """
    # Generate structure timestamps for measuring signal phase across multiple epochs
    times = timefunction(nobs, spacingbetween, duration, cadence)
    # Generate light curves with signals varied by reflex motion
    lcs = fluxfunction(times, pulsationperiod, m1, m2, separation, inclination, orbitalphase, 
                       centerphase, pulsationamplitude, noiselevel)
    # Generate O-C data using intrinsic period as reference period
    referenceperiod = pulsationperiod
    oc = OC(lcs, referenceperiod=referenceperiod)
    # Fit parabola to the OC data
    fit,cov = np.polyfit(oc["Epoch"], oc["OC"], 2, full=False, cov=True) # Fitting 2nd degree polynomial
    # If refine_period, update reference period based on linear term and redo
    if refine_period:
        referenceperiod = pulsationperiod + fit[1]*u.s
        # Recompute OC with new reference period
        oc = OC(lcs, referenceperiod=referenceperiod)
        fit,cov = np.polyfit(oc["Epoch"], oc["OC"], 2, full=False, cov=True) # 2nd degree polynomial
    # Calculate P-dot
    pdot= 2*fit[0]/referenceperiod.value # Calculates Pdot 
    pdoterr = 2*np.sqrt(cov[0,0])/referenceperiod.value # Calculates error on pdot using variance from covariance matrix
    # Determine reduced Chi^2 of the fit
    z = np.poly1d(fit) # Function from fit
    redchi2 = np.sum(((oc["OC"] - z(oc["Epoch"]))**2.) / (oc["OC_err"]**2.))/(len(oc)-3)
    # Prepare output
    output = pd.Series({"reference_period":referenceperiod,
                        "Pdot":pdot, 
                        "Pdot_err":pdoterr,
                        "reducedchi2":redchi2})
    if plot: # Plot results
        fig, ax = plt.subplots(figsize=(4,3),tight_layout=True)
        # Plot OC data
        ax.errorbar(oc["Epoch"]/1e6, oc["OC"], yerr = oc["OC_err"], linestyle='', color='0.3')
        # Plot best-fit polynomial
        ax.plot(oc["Epoch"]/1e6,z(oc["Epoch"]),color='red') # Plotting O-C fit
        # Axis labels
        ax.set_xlabel(r'Epoch ($10^6$ cycles)')
        ax.set_ylabel('$O-C$ (s)')
        # Top axis in years
        def epoch2year(x): # conversion functions for axes
            return x * 1e6 * referenceperiod.to(u.yr).value
        def year2epoch(x):
            return x / (referenceperiod.to(u.yr).value * 1e6)
        yrax = ax.secondary_xaxis('top', functions=(epoch2year, year2epoch))
        yrax.set_xlabel('Years')
        plt.show()
    return output #return results