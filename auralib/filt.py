"""
AuraQI module for various types of time and frequency domain filters.

Author:   Wes Hamlyn
Created:  3-Jan-2015
Last Mod: 1-Dec-2016

"""

from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
#from ._utils import _maybe_get_pandas_wrapper

def apply_filter(tdata, filt):
    """
    Convenience function to apply a frequency domain filter operator to time
    series data (i.e. seismic traces).  Note that both filter and time-series
    must have the same number of samples.
	
    Written by: Wes Hamlyn
    Created:    30-Nov-2016
    Modified:   19-Sep-2017
    """
    
    tdataf = np.fft.fft(tdata)
    aspec = np.abs(tdataf)
    pspec = np.arctan2(tdataf.imag, tdataf.real)
    
    aspec2 = (aspec*filt)
    tdataf2 = aspec2*np.exp(pspec*1j)
    tdata2 = np.fft.ifft(tdataf2)
    tdata2 = tdata2.real
    
    return tdata2


def ormsby(f1, f2, f3, f4, dt, nsamp):
    """
    Design a frequency domain trapezoidal bandpass filter
    
    Usage:
    ------
    filter = ormsby(f1, f2, f3, f4, dt, nsamp)
    
    f1: Low truncation frequency in Hz
    f2: Low cut frequency in Hz
    f3: High cut frequency in Hz
    f4: High truncation frequency in Hz
    dt: sample rate in seconds
    nsamp: length of filter in samples
    
    TO DO:
	1) Build cosine taper into low and high frequency filter ramps to avoid
	   sharp edges at filter corners.
	
    Written by: Wes Hamlyn
    Created:    30-Nov-2016
    Modified:   1-Dec-2016
	"""
        
    f1 = float(f1)
    f2 = float(f2)
    f3 = float(f3)
    f4 = float(f4)
    dt = float(dt)
    nsamp = int(nsamp)
    
    # Calculate slope and y-int for low frequency ramp
    if f1 == f2:
        pass
    else:
        M1 = 1/(f2-f1)
        b1 = -M1*f1
    
    # Calculate slope and y-int for high frequency ramp
    if f3 == f4:
        pass
    else:
        M2 = -1/(f4-f3)
        b2 = -M2*f4
    
    # Initialize frequency and filter arrays
    freq = np.fft.fftfreq(nsamp, dt)
    freq = np.fft.fftshift(freq)
    filt = np.zeros(nsamp)
    
    # Build low frequency filter ramp
    idx = np.nonzero((np.abs(freq)>=f1) & (np.abs(freq)<f2))
    if f1 == f2:
        filt[idx] = 0
    else:
        filt[idx] = M1*np.abs(freq)[idx]+b1
    
    # Build central filter flat
    idx = np.nonzero((np.abs(freq)>=f2) & (np.abs(freq)<=f3))
    filt[idx] = 1.0
    
    # Build high frequency filter ramp
    idx = np.nonzero((np.abs(freq)>f3) & (np.abs(freq)<=f4))
    if f3 == f4:
        filt[idx] = 0
    else:
        filt[idx] = M2*np.abs(freq)[idx]+b2
    
    # Un-shift the frequencies
    filt = np.fft.ifftshift(filt)
    
    return filt

def butterworth_lp(nyq, nsamp, fc, n, return_freq=False):
    """
    Low pass butterworth filter
    fc = high corner frequency
    n = frequency rolloff (# of poles)
	
    Written by: Wes Hamlyn
    Created:    3-Jan-2015
    Modified:   1-Dec-2016
    """
    dt = 0.5/nyq
    
    freq = np.fft.fftshift(np.fft.fftfreq(nsamp, dt))
    w = 2*np.pi*freq
    wc = 2*np.pi*fc

    filt = np.sqrt(1/(1+(w/wc)**(2.0*n)))
    
    # Un-shift the frequencies
    filt = np.fft.ifftshift(filt)
    freq = np.fft.ifftshift(freq)
    
    if return_freq:
        return filt, freq
    else:
        return filt
    
    
    
def butterworth_hp(nyq, nsamp, fc, n, return_freq=False):
    """
    High pass butterworth filter
    fc = low corner frequency
    n = frequency rolloff (# of poles)
	
    Written by: Wes Hamlyn
    Created:    3-Jan-2015
    Modified:   1-Dec-2016
    """
    dt = 0.5/nyq
    
    freq = np.fft.fftshift(np.fft.fftfreq(nsamp, dt))
    w = 2*np.pi*freq
    wc = 2*np.pi*fc
    
    filt = np.sqrt(-1/(1 + (w/wc)**(2.0*n)) + 1)
    
    # Un-shift the frequencies
    filt = np.fft.ifftshift(filt)
    freq = np.fft.ifftshift(freq)
    
    if return_freq:
        return filt, freq
    else:
        return filt


def butterworth_bp(nyq, nsamp, fc1, n1, fc2, n2, return_freq=True):
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
    dt = 0.5/nyq
    
    freq = np.fft.fftshift(np.fft.fftfreq(nsamp, dt))
    w = 2*np.pi*freq
    
    wc1 = 2*np.pi*fc1
    wc2 = 2*np.pi*fc2
    
    filt_lp = np.sqrt(1/(1+(w/wc2)**(2.0*n2)))
    filt_hp = np.sqrt(-1/(1 + (w/wc1)**(2.0*n1)) + 1)
    
    filt = filt_lp * filt_hp
    
    # Un-shift the frequencies
    filt = np.fft.ifftshift(filt)
    freq = np.fft.ifftshift(freq)

    if return_freq:
        return filt, freq
    else:
        return filt

    
def phaserot(data, deg):
    """
    Apply phase rotation of n degrees to a 1D array
	
    Written by: Wes Hamlyn
    Created:    8-Dec-2015
    Modified:   8-Dec-2015
    """
    
    from scipy.signal import hilbert
    deg = float(deg)
    
    deg = deg*np.pi/180.0
    datah = hilbert(data)
    datah = np.imag(datah)
    data = np.cos(deg)*data - np.sin(deg)*datah
    
    return data


def ampspec(tdata, dt, ax=-1, scale='amp', ls='k-', lw=1, label=None):
    """
    Convenience function for plotting amplitude, power, or dB spectrum 
    displays to an existing matplotlib axis.  Operates on 1D and 2D data arrays.
    If using 2D arrays, first dimension must correspond to trace number, second
    dimension must correspond to sample number on a trace.
    
    Written by: Wes Hamlyn
    Created:    30-Nov-2016
    Modified:   1-Dec-2016
    """
    
    
    if np.ndim(tdata) == 2:
        ntrc = tdata.shape[0]
        nsamp = tdata.shape[1]
    
    elif np.ndim(tdata) == 1:
        ntrc = 1
        nsamp = len(tdata)
        
    tdataf = np.fft.fftshift(np.fft.fft(tdata))
    
    if ntrc > 1:
        aspec = np.abs(tdataf).sum(axis=0) / ntrc
    else:
        aspec = np.abs(tdataf)
    
    pwrspec = aspec ** 2.0
    db = 20.0 * np.log10(aspec/np.max(aspec))

    flbl = np.fft.fftfreq(nsamp, dt)
    flbl = np.fft.fftshift(flbl)
    
    if ax==-1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Spectral Analysis Window')
    
    if scale == 'amp':
        hdl = ax.plot(flbl, aspec, ls, lw=lw, label=label)
        ax.set_ylabel('Amplitude')
    
    if scale == 'logamp':
        hdl = ax.semilogy(flbl, aspec, ls, lw=lw, label=label)
        ax.set_ylabel('Amplitude')
        ax.grid(True, which='minor', axis='y')

    elif scale == 'power':
        hdl = ax.plot(flbl, pwrspec, ls, lw=lw, label=label)
        ax.set_ylabel('Power')
    
    elif scale == 'db':
        hdl = ax.plot(flbl, db, ls, lw=lw, label=label)
        ax.set_ylabel('dB Down')

    ax.set_xlabel('Frequency (Hz)')
    ax.set_xlim([0, 0.5/dt])
    
    return hdl, ax


def fkspec(tdata, dx, dy):
    """
    Calculates the FK Spectrum of a 2D array.
    
    ToDo:
    - Add padding to power of 2
    - Add detrending
    - Add cosine taper
    
    Written by: Wes Hamlyn
    Created:    8-Dec-2015
    Modified:   8-Dec-2015
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
	
    Written by: Wes Hamlyn
    Created:    8-Dec-2015
    Modified:   8-Dec-2015
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


def KLT(a):
    """
    Returns Karhunen Loeve Transform of the input and the transformation matrix 
    and eigenvalues.
	
	*** IN DEVELOPMENT ***
    
    Ex:
    import numpy as np
    a  = np.array([[1,2,4],[2,3,10]])
    
    kk,m = KLT(a)
    print(kk)
    print(m)
    
    # to check, the following should return the original a
    print(np.dot(kk.T,m).T)
    
	Written by: Wes Hamlyn
    Created:    8-Dec-2015
    Modified:   8-Dec-2015
    """
	
    a = a - np.mean(a)
    val, vec = np.linalg.eig(np.cov(a))
    klt = np.dot(vec, a)
    return klt, vec, val
    
    
def hpfilter(X, lamb=20000):
    """
    Hodrick-Prescott filter
    Copied verbatim from statsmodels package for convenience to use in auralib.
    Edited to remove dependence on pandas (expect numpy array, WesH: Jan, 2019)

    Parameters
    ----------
    X : array-like
        The 1d ndarray timeseries to filter of length (nobs,) or (nobs,1)
    lamb : float
        The Hodrick-Prescott smoothing parameter. A value of 1600 is
        suggested for quarterly data. Ravn and Uhlig suggest using a value
        of 6.25 (1600/4**4) for annual data and 129600 (1600*3**4) for monthly
        data.

    Returns
    -------
    cycle : array
        The estimated cycle in the data given lamb.
    trend : array
        The estimated trend in the data given lamb.

    Examples
    ---------
    >>> import statsmodels.api as sm
    >>> import pandas as pd
    >>> dta = sm.datasets.macrodata.load_pandas().data
    >>> index = pd.DatetimeIndex(start='1959Q1', end='2009Q4', freq='Q')
    >>> dta.set_index(index, inplace=True)

    >>> cycle, trend = sm.tsa.filters.hpfilter(dta.realgdp, 1600)
    >>> gdp_decomp = dta[['realgdp']]
    >>> gdp_decomp["cycle"] = cycle
    >>> gdp_decomp["trend"] = trend

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> gdp_decomp[["realgdp", "trend"]]["2000-03-31":].plot(ax=ax,
    ...                                                      fontsize=16);
    >>> plt.show()

    .. plot:: plots/hpf_plot.py

    Notes
    -----
    The HP filter removes a smooth trend, `T`, from the data `X`. by solving

    min sum((X[t] - T[t])**2 + lamb*((T[t+1] - T[t]) - (T[t] - T[t-1]))**2)
     T   t

    Here we implemented the HP filter as a ridge-regression rule using
    scipy.sparse. In this sense, the solution can be written as

    T = inv(I - lamb*K'K)X

    where I is a nobs x nobs identity matrix, and K is a (nobs-2) x nobs matrix
    such that

    K[i,j] = 1 if i == j or i == j + 2
    K[i,j] = -2 if i == j + 1
    K[i,j] = 0 otherwise

    See Also
    --------
    statsmodels.tsa.filters.bk_filter.bkfilter
    statsmodels.tsa.filters.cf_filter.cffilter
    statsmodels.tsa.seasonal.seasonal_decompose

    References
    ----------
    Hodrick, R.J, and E. C. Prescott. 1980. "Postwar U.S. Business Cycles: An
        Empricial Investigation." `Carnegie Mellon University discussion
        paper no. 451`.
    Ravn, M.O and H. Uhlig. 2002. "Notes On Adjusted the Hodrick-Prescott
        Filter for the Frequency of Observations." `The Review of Economics and
        Statistics`, 84(2), 371-80.
    """
	
    #_pandas_wrapper = _maybe_get_pandas_wrapper(X)
    X = np.asarray(X, float)
    if X.ndim > 1:
        X = X.squeeze()
    nobs = len(X)
    I = sparse.eye(nobs,nobs)
    offsets = np.array([0,1,2])
    data = np.repeat([[1.],[-2.],[1.]], nobs, axis=1)
    K = sparse.dia_matrix((data, offsets), shape=(nobs-2,nobs))

    use_umfpack = True
    trend = spsolve(I+lamb*K.T.dot(K), X, use_umfpack=use_umfpack)

    cycle = X-trend
    #if _pandas_wrapper is not None:
    #    return _pandas_wrapper(cycle), _pandas_wrapper(trend)
    return cycle, trend
 
   
    
    
    
    
    
    
    
    