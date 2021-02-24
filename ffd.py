import numpy as np

def FFD(ED, TotDur=1., Lum=30., fluxerr=0., dur=[], logY=True):
    '''
    Generate a cumulative Flare Frequency Distribution figure, and 
    approximate uncertainties in both dimensions.
    
    Not a complicated task, just tedious.
    
    Parameters
    ----------
    ED : array of Equiv Dur's, need to include a luminosity!
    TotDur : total duration of observations, in days
    Lum : the luminosity of the star
    fluxerr : the average flux errors (in relative flux units!)
    dur : array of flare durations.
    
    Returns
    -------
    ffd_x, ffd_y, ffd_xerr, ffd_yerr
    
    X coordinate always assumed to be log_10(Energy)
    Y coordinate is log_10(N/Day) by default, but optionally is N/Day
    
    '''
    
    ss = np.argsort(np.array(ED))[::-1]
    ffd_x = np.log10(ED[ss]) + Lum
    
    Num = np.arange(1, len(ffd_x)+1)
    ffd_y = Num / TotDur
    
    # approximate the Poisson Y errors using Gehrels (1986), eqn 7, assuming S=1
    Perror = np.sqrt(Num + 0.75) + 1.0
    ffd_yerr = Perror / TotDur
    
    if logY:
        # transform FFD Y and Y Error into log10
        ffd_yerr = np.abs(ffd_yerr / np.log(10.) / ffd_y)
        ffd_y = np.log10(ffd_y)


    if len(dur)==len(ffd_x):
        # compute X uncertainties for FFD
        # Make basic S/N equation, assume error = 1/SN
        S2N = ED / np.sqrt(ED + (fluxerr * dur * 86400.))
        ffd_xerr = np.abs((1./S2N) / np.log(10.) / ED) # convert to log
    else:
        # not particularly meaningful, but an interesting shape. NOT reccomended
        print('Warning: Durations not set. Making bad assumptions about the FFD X Error!')
        ffd_xerr = (np.sqrt(10**(ffd_x-np.nanmax(ffd_x)))/10**(ffd_x-np.nanmax(ffd_x)))/np.nanmax(ffd_x)

    return ffd_x, ffd_y, ffd_xerr, ffd_yerr