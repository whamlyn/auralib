# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 14:00:20 2024

@author: whamlyn
"""

import numpy as np
import matplotlib.pyplot as plt
import auralib as aura
import pandas as pd

def butterworth_lp(freq, fc, n, return_freq=False):
    """
    Low pass butterworth filter
    fc = high corner frequency
    n = frequency rolloff (# of poles)
	
    Written by: Wes Hamlyn
    Created:    3-Jan-2015
    Modified:   1-Dec-2016
    """
    
    w = 2*np.pi*freq
    wc = 2*np.pi*fc

    filt = np.sqrt(1/(1+(w/wc)**(2.0*n)))
    
    
    if return_freq:
        return filt, freq
    else:
        return filt
    
    
    
def butterworth_hp(freq, fc, n, return_freq=False):
    """
    High pass butterworth filter
    fc = low corner frequency
    n = frequency rolloff (# of poles)
	
    Written by: Wes Hamlyn
    Created:    3-Jan-2015
    Modified:   1-Dec-2016
    """


    w = 2*np.pi*freq
    wc = 2*np.pi*fc
    
    filt = np.sqrt(-1/(1 + (w/wc)**(2.0*n)) + 1)
    
    if return_freq:
        return filt, freq
    else:
        return filt


def butterworth_bp(freq, fc1, n1, fc2, n2, return_freq=False):
    """
    Bandpass butterworth filter
    fc1 = low corner frequency
    n1  = low frequency rolloff (# of poles)
    fc2 = high corner frequency
    n2  = high frequency rolloff (# of poles)
	
    Written by: Wes Hamlyn
    Created:    3-Jan-2015
    Modified:   1-Dec-2016
    """

    w = 2*np.pi*freq
    
    wc1 = 2*np.pi*fc1
    wc2 = 2*np.pi*fc2
    
    filt_lp = np.sqrt(1/(1+(w/wc2)**(2.0*n2)))
    filt_hp = np.sqrt(-1/(1 + (w/wc1)**(2.0*n1)) + 1)
    
    filt = filt_lp * filt_hp

    if return_freq:
        return filt, freq
    else:
        return filt