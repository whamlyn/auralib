"""
auralib module containing various window functions for signal processing.

Author:   Wes Hamlyn
Created:  24-Mar-2020
Last Mod: 24-Mar-2020
"""

import numpy as np

def papoulis(nsamp, normed=True):
    """
    Calculate a time-domain Papoulis (a.k.a Bohman) window

    Input:
    ------
        nsamp = number of samples to include in the window
    
    Output:
    -------
        w = window amplitudes
        
    Reference:
    ----------
        https://prod-ng.sandia.gov/techlib-noauth/access-control.cgi/2017/174042.pdf
        page: 84
    """
    
    tmin = -0.5
    tmax = 0.5
    dt = (tmax-tmin)/(nsamp-1)

    t = np.arange(nsamp)*dt + tmin
    
    w = (np.pi**2)/4 * (1-2*np.abs(t))*np.cos(2*np.pi*np.abs(t)) + \
        np.pi/4*np.sin(2*np.pi*np.abs(t))
    
    if normed:
        w = w / w.max()
        
    return w


def cosine(nsamp, taper_len):
    """
    Calculate a time-domain window with cosine tapers at each end
    
    nsamp = total length of window operator
    nsamp_taper = length of cosine taper in samples at ends of window
    """
    
    # build cosine ramp from zero to one
    x = np.linspace(-np.pi, 0.0, taper_len)
    costaper = 0.5*(np.cos(x)+1)
    
    # make window with ones everwhere and add cosine tapers to ends
    w = np.ones(nsamp)
    w[0:taper_len] = costaper
    w[-taper_len:] = costaper[-1::-1]
    
    return w


