"""
AuraQI module for various types of time and frequency domain filters.

Author:   Wes Hamlyn
Created:  3-Jan-20156
Last Mod: 3-Jan-2016
"""

import numpy as np


def phaserot(data, deg):
    """
    apply phase rotation of n degrees to a 1D array
    """
    
    from scipy.signal import hilbert
    deg = float(deg)
    
    deg = deg*np.pi/180.0
    datah = hilbert(data)
    datah = np.imag(datah)
    data = np.cos(deg)*data - np.sin(deg)*datah
    
    return data
    



def bw_lp(nyq, nsamp, fc, n):
    """
    Low pass butterworth filter
    fc = high corner frequency
    n = frequency rolloff (# of poles)
    """
    
    freq = np.linspace(-nyq, nyq, nsamp)
    w = 2*np.pi*freq
    wc = 2*np.pi*fc

    filt = np.sqrt(1/(1+(w/wc)**(2.0*n)))
    
    return filt, freq
    
    
    
def bw_hp(nyq, nsamp, fc, n):
    """
    High pass butterworth filter
    fc = low corner frequency
    n = frequency rolloff (# of poles)
    """
    
    freq = np.linspace(-nyq, nyq, nsamp)
    w = 2*np.pi*freq
    wc = 2*np.pi*fc
    
    filt = np.sqrt(-1/(1 + (w/wc)**(2.0*n)) + 1)
    
    return filt, freq



def bw_bp(nyq, nsamp, fc1, n1, fc2, n2):
    """
    Bandpass butterworth filter
    fc1 = low corner frequency
    n1  = low frequency rolloff (# of poles)
    fc2 = high corner frequency
    n2  = high frequency rolloff (# of poles)
    """
    
    freq = np.linspace(-nyq, nyq, nsamp)
    w = 2*np.pi*freq
    
    wc1 = 2*np.pi*fc1
    wc2 = 2*np.pi*fc2
    
    filt_lp = np.sqrt(1/(1+(w/wc2)**(2.0*n2)))
    filt_hp = np.sqrt(-1/(1 + (w/wc1)**(2.0*n1)) + 1)
    
    filt = filt_lp * filt_hp
    
    return filt, freq




def fkspec(tdata, dx, dy):
    """
    Calculates the FK Spectrum of a 2D array.
    
    ToDo:
    - Add padding to power of 2
    - Add detrending
    - Add cosine taper
    
    Written by: Wes Hamlyn
    Created:    8-Dec-2016
    Modified:   8-Dec-2016
    """
    
    tdata = np.array(tdata)

    ntrc, nsamp = tdata.shape
    
    df = np.fft.fft2(tdata)
    df = np.fft.fftshift(df)

    aspec = np.sqrt(df.imag**2.0 + df.real**2.0)
    pspec = np.arctan2(df.imag, df.real)
    
    nyqx = 0.5/dx
    nyqy = 0.5/dy
    
    kx = np.linspace(-nyqx, nyqx, ntrc)
    ky = np.linspace(-nyqy, nyqy, nsamp)
    
    return aspec, pspec, kx, ky
    
    
def plot_fkspec(ax, aspec, dx, dy):
    """
    Plots the calculated FK Spectrum to a figure axis.
    """
    
    if ax != 0:
        nyqx = 0.5/dx
        nyqy = 0.5/dy
        
        hdl = ax.imshow(np.log10(aspec.T), extent=[-nyqx, nyqx, -nyqy, nyqy])
        
        ax.set_ylim([0, nyqy])
        ax.set_xlabel('Kx (cyc/m)')
        ax.set_ylabel('Freq (Hz)')
        ax.set_aspect('auto')
        
        fig = ax.get_figure()
        fig.colorbar(hdl, ax=ax)
        
        return hdl

def fk_get_dip(ax):
    
    fig = ax.get_figure()
    
    x, y = fig.ginput(1)
    
    M = x/y
    Y = M*x + 0
    

def KLT(a):
    """
    Returns Karhunen Loeve Transform of the input and the transformation matrix and eigenval
    
    Ex:
    import numpy as np
    a  = np.array([[1,2,4],[2,3,10]])
    
    kk,m = KLT(a)
    print(kk)
    print(m)
    
    # to check, the following should return the original a
    print(np.dot(kk.T,m).T)
        
    """
    a = a - np.mean(a)
    val, vec = np.linalg.eig(np.cov(a))
    klt = np.dot(vec, a)
    return klt, vec, val
    
    

    
    
    
    
    
    
    
    
    
    
    