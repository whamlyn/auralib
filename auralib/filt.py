"""
auralib module for various types of time and frequency domain filters.

Author:   Wes Hamlyn
Created:  3-Jan-2015
Last Mod: 2-Sept-2026

"""

from __future__ import absolute_import
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.ndimage import gaussian_filter1d

#from ._utils import _maybe_get_pandas_wrapper


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


def apply_filter(tdata, filt, detrend=None):
    """
    Convenience function to apply a frequency domain filter operator to time
    series data (i.e. seismic traces).  Note that both filter and time-series
    must have the same number of samples.
	
    Written by: Wes Hamlyn
    Created:    30-Nov-2016
    Modified:   19-Sep-2017
    """
    
    if detrend == 'mean':
        trend = np.nanmean(tdata)
    elif detrend == 'linear':
        trend = np.linspace(tdata[0], tdata[-1], len(tdata))
    else:
        trend = 0.0

    tdata = tdata - trend
    
    #win = cosine(len(tdata), 51)
    #tdata = tdata*win
    
    tdataf = np.fft.fft(tdata)
    aspec = np.abs(tdataf)
    pspec = np.arctan2(tdataf.imag, tdataf.real)
    
    aspec2 = (aspec*filt)
    tdataf2 = aspec2*np.exp(pspec*1j)
    tdata2 = np.fft.ifft(tdataf2)
    tdata2 = tdata2.real
    
    tdata2 = tdata2 + trend
    
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


def butterworth_bp(nyq, nsamp, fc1, n1, fc2, n2, return_freq=False):
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
    filt = filt / np.max(filt)
    
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
    
    rad = np.deg2rad(deg)
    datah = hilbert(data)
    datah = np.imag(datah)
    data = np.cos(rad)*data - np.sin(rad)*datah
    
    return data


def ampspec(tdata, dt, ax=-1, scale='amp', smoother=0.0, *args, **kwargs):
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
    
    flbl = np.fft.fftfreq(nsamp, dt)
    flbl = np.fft.fftshift(flbl)
    
    # Sometimes the spectrum may be noisy/spiky, this is a simple moving
    # average smoother to make the trends in the spectrum a bit more clear
    if smoother > 0:
        aspec_sm = []
        for each in flbl:
            fmin = each - smoother/2
            fmax = each + smoother/2
            mask = (flbl>=fmin) & (flbl<=fmax)
            aspec_sm.append(np.mean(aspec[mask]))
        
        aspec_sm = np.array(aspec_sm)
        aspec = aspec_sm
        
    pwrspec = aspec ** 2.0
    db = 20.0 * np.log10(aspec/np.max(aspec))


    
    if ax==-1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Spectral Analysis Window')
    
    if scale == 'amp':
        hdl = ax.plot(flbl, aspec, *args, **kwargs)
        ax.set_ylabel('Amplitude')
        
    if scale == 'amp_norm':
        hdl = ax.plot(flbl, aspec/np.max(aspec), *args, **kwargs)
        ax.set_ylabel('Normalized Amplitude')
    
    elif scale == 'logamp':
        hdl = ax.semilogy(flbl, aspec, *args, **kwargs)
        ax.set_ylabel('Amplitude')
        ax.grid(True, which='minor', axis='y')

    elif scale == 'power':
        hdl = ax.plot(flbl, pwrspec, *args, **kwargs)
        ax.set_ylabel('Power')
    
    elif scale == 'db':
        hdl = ax.plot(flbl, db, *args, **kwargs)
        ax.set_ylabel('dB Down')

    ax.set_xlabel('Frequency (Hz)')
    ax.set_xlim([0, 0.5/dt])
    
    return hdl, ax


def fkspec(tdata, dx=1, dy=1):
    """
    Calculates the FK Spectrum of a 2D array.
    
    ToDo:
    - Add padding to power of 2
    - Add detrending
    - Add cosine taper
    
    Written by: Wes Hamlyn
    Created:    8-Dec-2015
    Modified:   8-Dec-2022
    """
    
    tdata = np.array(tdata)

    ntrc, nsamp = tdata.shape
    
    df = np.fft.fft2(tdata)
    df = np.fft.fftshift(df)

    aspec = np.abs(df)
    pspec = np.angle(df)
    
    kx = np.fft.fftshift(np.fft.fftfreq(ntrc, dx))
    ky = np.fft.fftshift(np.fft.fftfreq(nsamp, dy))
    
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
 
   
def fx_noise_attenuation_v1(
    data,
    window_samples=256,
    overlap=0.5,
    rank=4,
    spatial_window=None,
    dt=None,
    fmin=None,
    fmax=None,
    strength=1.0,
):
    """
    Attenuate spatially incoherent noise using FX rank reduction.

    Parameters
    ----------
    data : array_like, shape (n_traces, n_samples)
        Input seismic section. The first dimension is trace number and
        the second dimension is time/depth sample number.

    window_samples : int, default=256
        Length of each overlapping time window, in samples.

    overlap : float, default=0.5
        Fractional overlap between adjacent time windows. Must satisfy
        0 <= overlap < 1.

    rank : int, default=4
        Number of singular components retained in each frequency slice.
        Smaller values produce stronger attenuation but increase the risk
        of damaging crossing dips or complex signal.

    spatial_window : int or None, default=None
        Number of rows in the Hankel trajectory matrix. If None, roughly
        half the number of traces is used.

    dt : float or None, default=None
        Sample interval in seconds. Required only when fmin or fmax is set.

    fmin, fmax : float or None
        Minimum and maximum frequencies, in Hz, to filter. Frequencies
        outside this range are passed unchanged.

    strength : float, default=1.0
        Blending strength:

            0.0 -> unchanged input
            1.0 -> fully filtered result

        Intermediate values provide more conservative attenuation.

    Returns
    -------
    filtered : ndarray, shape (n_traces, n_samples)
        Noise-attenuated seismic section.

    noise : ndarray, shape (n_traces, n_samples)
        Estimated noise, calculated as data - filtered.

    Notes
    -----
    NaN values are temporarily replaced by zero during processing and
    restored to NaN in the returned arrays.
    """
    x = np.asarray(data)

    if x.ndim != 2:
        raise ValueError(
            "data must have shape (n_traces, n_samples)."
        )

    if not np.issubdtype(x.dtype, np.number):
        raise TypeError("data must contain numeric values.")

    if not 0.0 <= overlap < 1.0:
        raise ValueError("overlap must satisfy 0 <= overlap < 1.")

    if rank < 1:
        raise ValueError("rank must be at least 1.")

    if not 0.0 <= strength <= 1.0:
        raise ValueError("strength must be between 0 and 1.")

    if (fmin is not None or fmax is not None) and dt is None:
        raise ValueError("dt is required when fmin or fmax is specified.")

    if dt is not None and dt <= 0:
        raise ValueError("dt must be positive.")

    n_traces, n_samples = x.shape

    if n_traces < 3:
        raise ValueError("At least three traces are required.")

    if n_samples < 2:
        raise ValueError("At least two samples per trace are required.")

    # Work in floating point and preserve the locations of missing samples.
    output_dtype = np.result_type(x.dtype, np.float64)
    x = x.astype(output_dtype, copy=True)

    nan_mask = ~np.isfinite(x)
    x[nan_mask] = 0.0

    window_samples = int(min(max(window_samples, 2), n_samples))
    hop = max(1, int(round(window_samples * (1.0 - overlap))))

    if spatial_window is None:
        spatial_window = n_traces // 2 + 1

    spatial_window = int(spatial_window)
    spatial_window = max(2, min(spatial_window, n_traces - 1))

    # H has shape (L, K), with L + K - 1 = n_traces.
    L = spatial_window
    K = n_traces - L + 1
    max_rank = min(L, K)
    rank = min(int(rank), max_rank)

    # Padding lets the Hann windows reconstruct the ends of the section.
    pad = window_samples // 2
    pad_mode = "reflect" if n_samples > 1 else "edge"

    padded = np.pad(
        x,
        ((0, 0), (pad, pad)),
        mode=pad_mode,
    )

    padded_samples = padded.shape[1]

    # Add enough padding for the final complete window.
    if padded_samples <= window_samples:
        n_windows = 1
    else:
        n_windows = int(
            np.ceil((padded_samples - window_samples) / hop)
        ) + 1

    required_samples = (n_windows - 1) * hop + window_samples
    extra = required_samples - padded_samples

    if extra > 0:
        padded = np.pad(
            padded,
            ((0, 0), (0, extra)),
            mode=pad_mode,
        )

    result = np.zeros_like(padded, dtype=output_dtype)
    weights = np.zeros(padded.shape[1], dtype=np.float64)

    # Square-root Hann gives analysis and synthesis windows whose product
    # is a conventional Hann overlap-add window.
    analysis_window = np.sqrt(np.hanning(window_samples))

    if np.all(analysis_window == 0):
        analysis_window = np.ones(window_samples)

    frequencies = (
        np.fft.rfftfreq(window_samples, d=dt)
        if dt is not None
        else None
    )

    def filter_frequency_slice(spatial_vector):
        """
        Low-rank Hankel reconstruction of one complex frequency slice.
        """
        # Construct Hankel trajectory matrix:
        # H[i, j] = spatial_vector[i + j]
        index = (
            np.arange(L)[:, None]
            + np.arange(K)[None, :]
        )
        H = spatial_vector[index]

        U, singular_values, Vh = np.linalg.svd(
            H,
            full_matrices=False,
        )

        H_filtered = (
            U[:, :rank] * singular_values[:rank]
        ) @ Vh[:rank, :]

        # Convert the Hankel matrix back to a spatial vector by averaging
        # along anti-diagonals.
        reconstructed = np.zeros(
            n_traces,
            dtype=np.complex128,
        )
        counts = np.zeros(n_traces, dtype=np.float64)

        for i in range(L):
            reconstructed[i:i + K] += H_filtered[i, :]
            counts[i:i + K] += 1.0

        return reconstructed / counts

    for window_index in range(n_windows):
        start = window_index * hop
        stop = start + window_samples

        segment = padded[:, start:stop]
        windowed_segment = segment * analysis_window[None, :]

        spectrum = np.fft.rfft(
            windowed_segment,
            axis=1,
        )

        filtered_spectrum = spectrum.copy()

        for frequency_index in range(spectrum.shape[1]):
            apply_filter = True

            if frequencies is not None:
                frequency = frequencies[frequency_index]

                if fmin is not None and frequency < fmin:
                    apply_filter = False

                if fmax is not None and frequency > fmax:
                    apply_filter = False

            # Leave DC unchanged because rank reduction of the zero-frequency
            # slice can introduce an unnecessary spatial bias.
            if frequency_index == 0:
                apply_filter = False

            if apply_filter:
                predicted = filter_frequency_slice(
                    spectrum[:, frequency_index]
                )

                filtered_spectrum[:, frequency_index] = (
                    (1.0 - strength)
                    * spectrum[:, frequency_index]
                    + strength * predicted
                )

        filtered_segment = np.fft.irfft(
            filtered_spectrum,
            n=window_samples,
            axis=1,
        )

        # Synthesis window.
        filtered_segment *= analysis_window[None, :]

        result[:, start:stop] += filtered_segment
        weights[start:stop] += analysis_window**2

    valid_weights = weights > np.finfo(float).eps
    result[:, valid_weights] /= weights[valid_weights][None, :]

    # Remove the padding.
    filtered = result[:, pad:pad + n_samples]

    noise = x - filtered

    filtered[nan_mask] = np.nan
    noise[nan_mask] = np.nan

    return filtered, noise


# Fx and Fxy noise attenuation updates

def _window_starts(n, window_length, overlap):
    """
    Return window starting indices, including a final window aligned
    with the end of the data.
    """
    if not 0.0 <= overlap < 1.0:
        raise ValueError("overlap must satisfy 0 <= overlap < 1.")

    window_length = min(int(window_length), n)

    if window_length < 1:
        raise ValueError("window_length must be at least 1.")

    if window_length == n:
        return [0]

    step = max(
        1,
        int(round(window_length * (1.0 - overlap))),
    )

    starts = list(range(0, n - window_length + 1, step))
    final_start = n - window_length

    if starts[-1] != final_start:
        starts.append(final_start)

    return starts


def _positive_hann(n):
    """
    Construct a Hann-like taper with nonzero endpoints.

    Nonzero endpoints prevent samples on the outside boundaries from
    receiving zero weight during overlap-add reconstruction.
    """
    if n <= 1:
        return np.ones(n, dtype=np.float64)

    return np.hanning(n + 2)[1:-1]


def _frequency_selection(
    window_samples,
    dt,
    fmin,
    fmax,
):
    """
    Return frequencies and a Boolean mask identifying frequencies
    that should be filtered.
    """
    if (fmin is not None or fmax is not None) and dt is None:
        raise ValueError(
            "dt is required when fmin or fmax is specified."
        )

    if dt is not None and dt <= 0:
        raise ValueError("dt must be positive.")

    if fmin is not None and fmin < 0:
        raise ValueError("fmin cannot be negative.")

    if fmax is not None and fmax < 0:
        raise ValueError("fmax cannot be negative.")

    if (
        fmin is not None
        and fmax is not None
        and fmin > fmax
    ):
        raise ValueError("fmin cannot be greater than fmax.")

    if dt is None:
        frequencies = np.arange(
            window_samples // 2 + 1,
            dtype=np.float64,
        )
    else:
        frequencies = np.fft.rfftfreq(
            window_samples,
            d=dt,
        )

    selected = np.ones(
        frequencies.shape,
        dtype=bool,
    )

    # Leave DC unchanged.
    selected[0] = False

    if fmin is not None:
        selected &= frequencies >= fmin

    if fmax is not None:
        selected &= frequencies <= fmax

    return frequencies, selected


def fx_noise_attenuation_v2(
    data,
    window_samples=256,
    time_overlap=0.5,
    rank=4,
    embedding_length=None,
    dt=None,
    fmin=None,
    fmax=None,
    strength=1.0,
):
    """
    Attenuate spatially incoherent noise using FX Hankel rank reduction.

    Parameters
    ----------
    data : array_like, shape (n_traces, n_samples)
        Input seismic section.

    window_samples : int, default=256
        Length of each time window, in samples.

    time_overlap : float, default=0.5
        Fractional overlap between adjacent time windows.

    rank : int, default=4
        Number of singular components retained in each Hankel matrix.

    embedding_length : int or None, default=None
        Number of rows in the Hankel trajectory matrix.

        This is not a local spatial aperture. All traces contribute to
        the Hankel reconstruction.

        If None, approximately half the number of traces is used.

    dt : float or None, default=None
        Temporal sample interval in seconds.

    fmin, fmax : float or None
        Frequency range, in Hz, to filter. Frequencies outside this
        range are passed unchanged.

    strength : float, default=1.0
        Filter blending strength:

            0.0 = unchanged input
            1.0 = fully rank-reduced output

    Returns
    -------
    filtered : ndarray
        Filtered section with shape (n_traces, n_samples).

    noise : ndarray
        Estimated noise, calculated as data - filtered.
    """
    input_data = np.asarray(data)

    if input_data.ndim != 2:
        raise ValueError(
            "data must have shape (n_traces, n_samples)."
        )

    if not np.issubdtype(input_data.dtype, np.number):
        raise TypeError("data must contain numeric values.")

    if not 0.0 <= time_overlap < 1.0:
        raise ValueError(
            "time_overlap must satisfy 0 <= time_overlap < 1."
        )

    if not 0.0 <= strength <= 1.0:
        raise ValueError(
            "strength must satisfy 0 <= strength <= 1."
        )

    if rank < 1:
        raise ValueError("rank must be at least 1.")

    n_traces, n_samples = input_data.shape

    if n_traces < 3:
        raise ValueError("At least three traces are required.")

    if n_samples < 2:
        raise ValueError(
            "At least two samples per trace are required."
        )

    output_dtype = np.result_type(
        input_data.dtype,
        np.float64,
    )

    working_data = input_data.astype(
        output_dtype,
        copy=True,
    )

    invalid_mask = ~np.isfinite(working_data)
    working_data[invalid_mask] = 0.0

    window_samples = min(
        max(int(window_samples), 2),
        n_samples,
    )

    if embedding_length is None:
        embedding_length = n_traces // 2 + 1

    embedding_length = int(embedding_length)

    if not 2 <= embedding_length <= n_traces - 1:
        raise ValueError(
            "embedding_length must be between 2 and "
            "n_traces - 1."
        )

    n_hankel_rows = embedding_length
    n_hankel_columns = (
        n_traces - n_hankel_rows + 1
    )

    maximum_rank = min(
        n_hankel_rows,
        n_hankel_columns,
    )

    rank = min(int(rank), maximum_rank)

    _, selected_frequencies = _frequency_selection(
        window_samples=window_samples,
        dt=dt,
        fmin=fmin,
        fmax=fmax,
    )

    time_starts = _window_starts(
        n=n_samples,
        window_length=window_samples,
        overlap=time_overlap,
    )

    # Square-root taper for analysis and synthesis.
    time_taper = np.sqrt(
        _positive_hann(window_samples)
    )

    filtered_sum = np.zeros_like(
        working_data,
        dtype=output_dtype,
    )

    weight_sum = np.zeros(
        n_samples,
        dtype=np.float64,
    )

    hankel_indices = (
        np.arange(n_hankel_rows)[:, None]
        + np.arange(n_hankel_columns)[None, :]
    )

    def reconstruct_frequency_slice(spatial_vector):
        """
        Perform low-rank Hankel reconstruction of a complex spatial
        frequency vector.
        """
        hankel_matrix = spatial_vector[hankel_indices]

        U, singular_values, Vh = np.linalg.svd(
            hankel_matrix,
            full_matrices=False,
        )

        filtered_hankel = (
            U[:, :rank] * singular_values[:rank]
        ) @ Vh[:rank, :]

        reconstructed = np.zeros(
            n_traces,
            dtype=np.complex128,
        )

        counts = np.zeros(
            n_traces,
            dtype=np.float64,
        )

        # Average the anti-diagonals of the filtered Hankel matrix.
        for row in range(n_hankel_rows):
            reconstructed[
                row:row + n_hankel_columns
            ] += filtered_hankel[row, :]

            counts[
                row:row + n_hankel_columns
            ] += 1.0

        return reconstructed / counts

    for time_start in time_starts:
        time_stop = time_start + window_samples

        segment = working_data[:, time_start:time_stop]

        analysis_segment = (
            segment * time_taper[None, :]
        )

        spectrum = np.fft.rfft(
            analysis_segment,
            axis=1,
        )

        filtered_spectrum = spectrum.copy()

        for frequency_index in np.flatnonzero(
            selected_frequencies
        ):
            original_slice = spectrum[:, frequency_index]

            predicted_slice = reconstruct_frequency_slice(
                original_slice
            )

            filtered_spectrum[:, frequency_index] = (
                (1.0 - strength) * original_slice
                + strength * predicted_slice
            )

        filtered_segment = np.fft.irfft(
            filtered_spectrum,
            n=window_samples,
            axis=1,
        )

        # Synthesis taper.
        filtered_segment *= time_taper[None, :]

        filtered_sum[:, time_start:time_stop] += (
            filtered_segment
        )

        weight_sum[time_start:time_stop] += (
            time_taper**2
        )

    filtered = filtered_sum / weight_sum[None, :]

    noise = working_data - filtered

    filtered[invalid_mask] = np.nan
    noise[invalid_mask] = np.nan

    return filtered, noise


def fxy_noise_attenuation_v2(
    data,
    window_samples=256,
    time_overlap=0.5,
    spatial_window=(16, 16),
    spatial_overlap=(0.5, 0.5),
    rank=4,
    dt=None,
    fmin=None,
    fmax=None,
    strength=1.0,
):
    """
    Attenuate spatially incoherent noise using local FXY eigenimage
    rank reduction.

    Parameters
    ----------
    data : array_like
        Input volume with shape:

            (n_inline, n_crossline, n_samples)

    window_samples : int, default=256
        Length of each time window, in samples.

    time_overlap : float, default=0.5
        Fractional overlap between time windows.

    spatial_window : tuple of int, default=(16, 16)
        Local spatial aperture:

            (n_inline_window, n_crossline_window)

        If None, the complete inline-crossline aperture is used.

    spatial_overlap : float or tuple of float, default=(0.5, 0.5)
        Fractional overlap between local spatial windows.

        A scalar applies the same overlap in both directions.

    rank : int, default=4
        Number of singular components retained in each local
        inline-crossline frequency slice.

    dt : float or None, default=None
        Temporal sample interval in seconds.

    fmin, fmax : float or None
        Frequency range, in Hz, to filter.

    strength : float, default=1.0
        Filter blending strength:

            0.0 = unchanged input
            1.0 = fully rank-reduced output

    Returns
    -------
    filtered : ndarray
        Filtered volume with the same shape as data.

    noise : ndarray
        Estimated noise:

            noise = data - filtered

    Notes
    -----
    The method assumes the data are regularly sampled along the inline
    and crossline directions.

    Each local frequency slice has shape:

        (inline_window, crossline_window)

    and is reconstructed using a truncated complex-valued SVD.
    """
    input_data = np.asarray(data)

    if input_data.ndim != 3:
        raise ValueError(
            "data must have shape "
            "(n_inline, n_crossline, n_samples)."
        )

    if not np.issubdtype(input_data.dtype, np.number):
        raise TypeError("data must contain numeric values.")

    if not 0.0 <= time_overlap < 1.0:
        raise ValueError(
            "time_overlap must satisfy 0 <= time_overlap < 1."
        )

    if not 0.0 <= strength <= 1.0:
        raise ValueError(
            "strength must satisfy 0 <= strength <= 1."
        )

    if rank < 1:
        raise ValueError("rank must be at least 1.")

    n_inline, n_crossline, n_samples = input_data.shape

    if n_inline < 2 or n_crossline < 2:
        raise ValueError(
            "At least two inlines and two crosslines are required."
        )

    if n_samples < 2:
        raise ValueError(
            "At least two samples per trace are required."
        )

    output_dtype = np.result_type(
        input_data.dtype,
        np.float64,
    )

    working_data = input_data.astype(
        output_dtype,
        copy=True,
    )

    invalid_mask = ~np.isfinite(working_data)
    working_data[invalid_mask] = 0.0

    window_samples = min(
        max(int(window_samples), 2),
        n_samples,
    )

    if spatial_window is None:
        inline_window = n_inline
        crossline_window = n_crossline
    else:
        if len(spatial_window) != 2:
            raise ValueError(
                "spatial_window must contain two values."
            )

        inline_window = min(
            max(int(spatial_window[0]), 2),
            n_inline,
        )

        crossline_window = min(
            max(int(spatial_window[1]), 2),
            n_crossline,
        )

    if np.isscalar(spatial_overlap):
        inline_overlap = float(spatial_overlap)
        crossline_overlap = float(spatial_overlap)
    else:
        if len(spatial_overlap) != 2:
            raise ValueError(
                "spatial_overlap must be a scalar or contain "
                "two values."
            )

        inline_overlap = float(spatial_overlap[0])
        crossline_overlap = float(spatial_overlap[1])

    if not 0.0 <= inline_overlap < 1.0:
        raise ValueError(
            "Inline spatial overlap must satisfy 0 <= overlap < 1."
        )

    if not 0.0 <= crossline_overlap < 1.0:
        raise ValueError(
            "Crossline spatial overlap must satisfy "
            "0 <= overlap < 1."
        )

    maximum_rank = min(
        inline_window,
        crossline_window,
    )

    rank = min(int(rank), maximum_rank)

    _, selected_frequencies = _frequency_selection(
        window_samples=window_samples,
        dt=dt,
        fmin=fmin,
        fmax=fmax,
    )

    inline_starts = _window_starts(
        n=n_inline,
        window_length=inline_window,
        overlap=inline_overlap,
    )

    crossline_starts = _window_starts(
        n=n_crossline,
        window_length=crossline_window,
        overlap=crossline_overlap,
    )

    time_starts = _window_starts(
        n=n_samples,
        window_length=window_samples,
        overlap=time_overlap,
    )

    # Square-root analysis/synthesis tapers.
    inline_taper = np.sqrt(
        _positive_hann(inline_window)
    )

    crossline_taper = np.sqrt(
        _positive_hann(crossline_window)
    )

    time_taper = np.sqrt(
        _positive_hann(window_samples)
    )

    spatial_taper = (
        inline_taper[:, None]
        * crossline_taper[None, :]
    )

    analysis_taper = (
        spatial_taper[:, :, None]
        * time_taper[None, None, :]
    )

    filtered_sum = np.zeros_like(
        working_data,
        dtype=output_dtype,
    )

    weight_sum = np.zeros_like(
        working_data,
        dtype=np.float64,
    )

    taper_weight = analysis_taper**2

    for inline_start in inline_starts:
        inline_stop = inline_start + inline_window

        for crossline_start in crossline_starts:
            crossline_stop = (
                crossline_start + crossline_window
            )

            for time_start in time_starts:
                time_stop = time_start + window_samples

                tile = working_data[
                    inline_start:inline_stop,
                    crossline_start:crossline_stop,
                    time_start:time_stop,
                ]

                analysis_tile = tile * analysis_taper

                spectrum = np.fft.rfft(
                    analysis_tile,
                    axis=2,
                )

                filtered_spectrum = spectrum.copy()

                for frequency_index in np.flatnonzero(
                    selected_frequencies
                ):
                    frequency_slice = spectrum[
                        :,
                        :,
                        frequency_index,
                    ]

                    U, singular_values, Vh = np.linalg.svd(
                        frequency_slice,
                        full_matrices=False,
                    )

                    rank_reduced_slice = (
                        U[:, :rank]
                        * singular_values[:rank]
                    ) @ Vh[:rank, :]

                    filtered_spectrum[
                        :,
                        :,
                        frequency_index,
                    ] = (
                        (1.0 - strength) * frequency_slice
                        + strength * rank_reduced_slice
                    )

                filtered_tile = np.fft.irfft(
                    filtered_spectrum,
                    n=window_samples,
                    axis=2,
                )

                # Synthesis taper.
                filtered_tile *= analysis_taper

                filtered_sum[
                    inline_start:inline_stop,
                    crossline_start:crossline_stop,
                    time_start:time_stop,
                ] += filtered_tile

                weight_sum[
                    inline_start:inline_stop,
                    crossline_start:crossline_stop,
                    time_start:time_stop,
                ] += taper_weight

    valid_weights = weight_sum > np.finfo(float).eps

    filtered = np.zeros_like(
        filtered_sum,
        dtype=output_dtype,
    )

    filtered[valid_weights] = (
        filtered_sum[valid_weights]
        / weight_sum[valid_weights]
    )

    # This should rarely be needed because the tapers have nonzero edges.
    filtered[~valid_weights] = working_data[~valid_weights]

    noise = working_data - filtered

    filtered[invalid_mask] = np.nan
    noise[invalid_mask] = np.nan

    return filtered, noise


import numpy as np


def fx_mssa_noise_attenuation(
    data,
    window_samples=256,
    time_overlap=0.5,
    rank=4,
    embedding_length=11,
    n_iterations=1,
    damping_exponent=None,
    dt=None,
    fmin=None,
    fmax=None,
    strength=1.0,
):
    """
    Attenuate spatially incoherent noise using frequency-domain
    multichannel singular spectrum analysis (MSSA).

    Parameters
    ----------
    data : array_like, shape (n_traces, n_samples)
        Input seismic section.

        Axis 0: trace number
        Axis 1: time or depth sample

    window_samples : int, default=256
        Length of each time-analysis window, in samples.

    time_overlap : float, default=0.5
        Fractional overlap between adjacent time windows.

    rank : int, default=4
        Number of signal components retained in each Hankel matrix.

    embedding_length : int or None, default=None
        Number of rows in the Hankel trajectory matrix.

        If None, a nearly square Hankel matrix is used:

            embedding_length = n_traces // 2 + 1

    n_iterations : int, default=1
        Number of Cadzow iterations.

        One iteration performs conventional MSSA rank reduction.
        Additional iterations alternate between low-rank projection and
        Hankel-structure projection.

    damping_exponent : float or None, default=None
        Controls singular-value damping.

        None
            Hard rank truncation. The first `rank` singular values are
            retained unchanged.

        Positive float
            Apply damped rank reduction:

                s_out[i] = s[i] *
                    (1 - (s[rank] / s[i])**damping_exponent)

            where s[rank] is the first rejected singular value.

        Values around 1 to 4 are reasonable starting points.

    dt : float or None, default=None
        Sample interval in seconds. Required if fmin or fmax is used.

    fmin, fmax : float or None
        Frequency range, in Hz, over which MSSA is applied. Frequencies
        outside this range are passed unchanged.

    strength : float, default=1.0
        Blend between the input and MSSA-reconstructed data:

            0.0 = unchanged input
            1.0 = complete MSSA output

    Returns
    -------
    filtered : ndarray, shape (n_traces, n_samples)
        MSSA-filtered seismic section.

    noise : ndarray, shape (n_traces, n_samples)
        Estimated noise:

            noise = data - filtered
    """
    input_data = np.asarray(data)

    if input_data.ndim != 2:
        raise ValueError(
            "data must have shape (n_traces, n_samples)."
        )

    if not np.issubdtype(input_data.dtype, np.number):
        raise TypeError("data must contain numeric values.")

    if not 0.0 <= time_overlap < 1.0:
        raise ValueError(
            "time_overlap must satisfy 0 <= time_overlap < 1."
        )

    if not 0.0 <= strength <= 1.0:
        raise ValueError(
            "strength must satisfy 0 <= strength <= 1."
        )

    if rank < 1:
        raise ValueError("rank must be at least 1.")

    if n_iterations < 1:
        raise ValueError("n_iterations must be at least 1.")

    if damping_exponent is not None:
        if damping_exponent <= 0:
            raise ValueError(
                "damping_exponent must be positive or None."
            )

    if (fmin is not None or fmax is not None) and dt is None:
        raise ValueError(
            "dt is required when fmin or fmax is specified."
        )

    if dt is not None and dt <= 0:
        raise ValueError("dt must be positive.")

    if fmin is not None and fmin < 0:
        raise ValueError("fmin cannot be negative.")

    if fmax is not None and fmax < 0:
        raise ValueError("fmax cannot be negative.")

    if (
        fmin is not None
        and fmax is not None
        and fmin > fmax
    ):
        raise ValueError("fmin cannot exceed fmax.")

    n_traces, n_samples = input_data.shape

    if n_traces < 3:
        raise ValueError("At least three traces are required.")

    if n_samples < 2:
        raise ValueError(
            "At least two samples per trace are required."
        )

    output_dtype = np.result_type(
        input_data.dtype,
        np.float64,
    )

    working_data = input_data.astype(
        output_dtype,
        copy=True,
    )

    invalid_mask = ~np.isfinite(working_data)
    working_data[invalid_mask] = 0.0

    window_samples = min(
        max(int(window_samples), 2),
        n_samples,
    )

    if embedding_length is None:
        embedding_length = n_traces // 2 + 1

    embedding_length = int(embedding_length)

    if not 2 <= embedding_length <= n_traces - 1:
        raise ValueError(
            "embedding_length must lie between 2 and "
            "n_traces - 1."
        )

    n_rows = embedding_length
    n_columns = n_traces - n_rows + 1

    maximum_rank = min(n_rows, n_columns)
    rank = min(int(rank), maximum_rank)

    # Precalculate Hankel indices:
    #
    #     H[i, j] = spatial_vector[i + j]
    #
    hankel_indices = (
        np.arange(n_rows)[:, None]
        + np.arange(n_columns)[None, :]
    )

    def make_hankel(spatial_vector):
        """Embed one spatial frequency vector in a Hankel matrix."""
        return spatial_vector[hankel_indices]

    def diagonal_average(hankel_matrix):
        """
        Convert a Hankel matrix back to a spatial vector by averaging
        along its anti-diagonals.
        """
        reconstructed = np.zeros(
            n_traces,
            dtype=np.complex128,
        )

        counts = np.zeros(
            n_traces,
            dtype=np.float64,
        )

        for row in range(n_rows):
            reconstructed[
                row:row + n_columns
            ] += hankel_matrix[row, :]

            counts[
                row:row + n_columns
            ] += 1.0

        return reconstructed / counts

    def reduce_rank(hankel_matrix):
        """Apply hard or damped singular-value rank reduction."""
        U, singular_values, Vh = np.linalg.svd(
            hankel_matrix,
            full_matrices=False,
        )

        retained_values = singular_values[:rank].copy()

        if (
            damping_exponent is not None
            and rank < singular_values.size
        ):
            first_rejected = singular_values[rank]

            safe_values = np.maximum(
                retained_values,
                np.finfo(float).eps,
            )

            damping_weights = 1.0 - (
                first_rejected / safe_values
            ) ** damping_exponent

            damping_weights = np.clip(
                damping_weights,
                0.0,
                1.0,
            )

            retained_values *= damping_weights

        return (
            U[:, :rank] * retained_values
        ) @ Vh[:rank, :]

    def mssa_reconstruct(spatial_vector):
        """
        Apply one or more Cadzow MSSA iterations to a spatial
        frequency vector.
        """
        estimate = spatial_vector.astype(
            np.complex128,
            copy=True,
        )

        for _ in range(n_iterations):
            hankel_matrix = make_hankel(estimate)

            low_rank_hankel = reduce_rank(
                hankel_matrix
            )

            # Diagonal averaging projects the matrix back onto the set
            # of Hankel-structured matrices.
            estimate = diagonal_average(
                low_rank_hankel
            )

        return estimate

    def window_starts(n, length, overlap):
        """Return complete overlapping-window start locations."""
        if length >= n:
            return [0]

        step = max(
            1,
            int(round(length * (1.0 - overlap))),
        )

        starts = list(
            range(0, n - length + 1, step)
        )

        final_start = n - length

        if starts[-1] != final_start:
            starts.append(final_start)

        return starts

    def positive_hann(length):
        """
        Hann-like taper with nonzero endpoints for stable boundary
        reconstruction.
        """
        if length <= 1:
            return np.ones(
                length,
                dtype=np.float64,
            )

        return np.hanning(length + 2)[1:-1]

    time_starts = window_starts(
        n_samples,
        window_samples,
        time_overlap,
    )

    # The square-root taper is applied during both analysis and synthesis.
    time_taper = np.sqrt(
        positive_hann(window_samples)
    )

    frequencies = np.fft.rfftfreq(
        window_samples,
        d=dt if dt is not None else 1.0,
    )

    selected_frequencies = np.ones(
        frequencies.size,
        dtype=bool,
    )

    # Leave DC unchanged.
    selected_frequencies[0] = False

    if fmin is not None:
        selected_frequencies &= frequencies >= fmin

    if fmax is not None:
        selected_frequencies &= frequencies <= fmax

    filtered_sum = np.zeros_like(
        working_data,
        dtype=output_dtype,
    )

    weight_sum = np.zeros(
        n_samples,
        dtype=np.float64,
    )

    for time_start in time_starts:
        time_stop = time_start + window_samples

        segment = working_data[
            :,
            time_start:time_stop,
        ]

        analysis_segment = (
            segment * time_taper[None, :]
        )

        spectrum = np.fft.rfft(
            analysis_segment,
            axis=1,
        )

        filtered_spectrum = spectrum.copy()

        for frequency_index in np.flatnonzero(
            selected_frequencies
        ):
            original_slice = spectrum[
                :,
                frequency_index,
            ]

            mssa_slice = mssa_reconstruct(
                original_slice
            )

            filtered_spectrum[
                :,
                frequency_index,
            ] = (
                (1.0 - strength) * original_slice
                + strength * mssa_slice
            )

        filtered_segment = np.fft.irfft(
            filtered_spectrum,
            n=window_samples,
            axis=1,
        )

        # Apply the synthesis taper.
        filtered_segment *= time_taper[None, :]

        filtered_sum[
            :,
            time_start:time_stop,
        ] += filtered_segment

        weight_sum[
            time_start:time_stop
        ] += time_taper**2

    filtered = (
        filtered_sum
        / weight_sum[None, :]
    )

    noise = working_data - filtered

    filtered[invalid_mask] = np.nan
    noise[invalid_mask] = np.nan

    return filtered, noise


def fxy_mssa_noise_attenuation(
    data,
    window_samples=256,
    time_overlap=0.5,
    spatial_window=(16, 16),
    spatial_overlap=(0.5, 0.5),
    embedding_shape=11,
    rank=4,
    n_iterations=1,
    damping_exponent=None,
    dt=None,
    fmin=None,
    fmax=None,
    strength=1.0,
):
    """
    Attenuate incoherent noise in a 3D seismic volume using local
    frequency-space multichannel singular spectrum analysis (MSSA).

    Parameters
    ----------
    data : array_like
        Input seismic volume with shape:

            (n_inline, n_crossline, n_samples)

    window_samples : int, default=256
        Length of each time-analysis window, in samples.

    time_overlap : float, default=0.5
        Fractional overlap between adjacent time windows.

    spatial_window : tuple of int or None, default=(16, 16)
        Local spatial processing aperture:

            (inline_window, crossline_window)

        If None, the complete spatial aperture is processed as one tile.

    spatial_overlap : float or tuple of float, default=(0.5, 0.5)
        Fractional overlap between adjacent spatial tiles. A scalar uses
        the same overlap in both directions.

    embedding_shape : tuple of int or None, default=None
        Dimensions of the 2D SSA embedding patch:

            (inline_embedding, crossline_embedding)

        If None, approximately half of each spatial-window dimension is
        used.

    rank : int, default=4
        Number of singular components retained.

    n_iterations : int, default=1
        Number of Cadzow iterations. Each iteration alternates between
        low-rank projection and trajectory-matrix reconstruction.

    damping_exponent : float or None, default=None
        Singular-value damping exponent.

        None:
            Use conventional hard rank truncation.

        Positive number:
            Shrink the retained singular values relative to the first
            rejected singular value. Values between 1 and 4 are useful
            starting points.

    dt : float or None, default=None
        Temporal sample interval in seconds. Required when fmin or fmax
        is specified.

    fmin, fmax : float or None
        Frequency range, in Hz, over which MSSA is applied. Frequencies
        outside this range are passed unchanged.

    strength : float, default=1.0
        Blend between the input and MSSA estimate:

            0.0 = unchanged input
            1.0 = complete MSSA output

    Returns
    -------
    filtered : ndarray
        Filtered volume with the same shape as the input.

    noise : ndarray
        Estimated removed noise:

            noise = data - filtered

    Notes
    -----
    MSSA is applied independently to each temporal-frequency slice
    within each local spatial tile.

    The trajectory matrix has dimensions:

        (Lx * Ly, Kx * Ky)

    where:

        Kx = Nx - Lx + 1
        Ky = Ny - Ly + 1

    Nx and Ny are the local spatial-window dimensions, and Lx and Ly
    are the embedding dimensions.
    """
    input_data = np.asarray(data)

    if input_data.ndim != 3:
        raise ValueError(
            "data must have shape "
            "(n_inline, n_crossline, n_samples)."
        )

    if not np.issubdtype(input_data.dtype, np.number):
        raise TypeError("data must contain numeric values.")

    if not 0.0 <= time_overlap < 1.0:
        raise ValueError(
            "time_overlap must satisfy 0 <= time_overlap < 1."
        )

    if not 0.0 <= strength <= 1.0:
        raise ValueError(
            "strength must satisfy 0 <= strength <= 1."
        )

    if rank < 1:
        raise ValueError("rank must be at least 1.")

    if n_iterations < 1:
        raise ValueError("n_iterations must be at least 1.")

    if damping_exponent is not None:
        if damping_exponent <= 0:
            raise ValueError(
                "damping_exponent must be positive or None."
            )

    if (fmin is not None or fmax is not None) and dt is None:
        raise ValueError(
            "dt is required when fmin or fmax is specified."
        )

    if dt is not None and dt <= 0:
        raise ValueError("dt must be positive.")

    if fmin is not None and fmin < 0:
        raise ValueError("fmin cannot be negative.")

    if fmax is not None and fmax < 0:
        raise ValueError("fmax cannot be negative.")

    if (
        fmin is not None
        and fmax is not None
        and fmin > fmax
    ):
        raise ValueError("fmin cannot exceed fmax.")

    n_inline, n_crossline, n_samples = input_data.shape

    if n_inline < 3 or n_crossline < 3:
        raise ValueError(
            "At least three inlines and three crosslines are required."
        )

    if n_samples < 2:
        raise ValueError(
            "At least two samples per trace are required."
        )

    output_dtype = np.result_type(
        input_data.dtype,
        np.float64,
    )

    working_data = input_data.astype(
        output_dtype,
        copy=True,
    )

    invalid_mask = ~np.isfinite(working_data)
    working_data[invalid_mask] = 0.0

    window_samples = min(
        max(int(window_samples), 2),
        n_samples,
    )

    # --------------------------------------------------------------
    # Parse the local spatial-window dimensions
    # --------------------------------------------------------------
    if spatial_window is None:
        inline_window = n_inline
        crossline_window = n_crossline
    else:
        if len(spatial_window) != 2:
            raise ValueError(
                "spatial_window must contain two integers."
            )

        inline_window = min(
            max(int(spatial_window[0]), 3),
            n_inline,
        )

        crossline_window = min(
            max(int(spatial_window[1]), 3),
            n_crossline,
        )

    # --------------------------------------------------------------
    # Parse spatial overlap
    # --------------------------------------------------------------
    if np.isscalar(spatial_overlap):
        inline_overlap = float(spatial_overlap)
        crossline_overlap = float(spatial_overlap)
    else:
        if len(spatial_overlap) != 2:
            raise ValueError(
                "spatial_overlap must be a scalar or contain "
                "two values."
            )

        inline_overlap = float(spatial_overlap[0])
        crossline_overlap = float(spatial_overlap[1])

    if not 0.0 <= inline_overlap < 1.0:
        raise ValueError(
            "Inline overlap must satisfy 0 <= overlap < 1."
        )

    if not 0.0 <= crossline_overlap < 1.0:
        raise ValueError(
            "Crossline overlap must satisfy 0 <= overlap < 1."
        )

    # --------------------------------------------------------------
    # Parse MSSA embedding dimensions
    # --------------------------------------------------------------
    if embedding_shape is None:
        inline_embedding = inline_window // 2 + 1
        crossline_embedding = crossline_window // 2 + 1
    else:
        if len(embedding_shape) != 2:
            raise ValueError(
                "embedding_shape must contain two integers."
            )

        inline_embedding = int(embedding_shape[0])
        crossline_embedding = int(embedding_shape[1])

    if not 2 <= inline_embedding <= inline_window - 1:
        raise ValueError(
            "The inline embedding dimension must be between "
            "2 and inline_window - 1."
        )

    if not 2 <= crossline_embedding <= crossline_window - 1:
        raise ValueError(
            "The crossline embedding dimension must be between "
            "2 and crossline_window - 1."
        )

    inline_positions = (
        inline_window - inline_embedding + 1
    )

    crossline_positions = (
        crossline_window - crossline_embedding + 1
    )

    trajectory_rows = (
        inline_embedding * crossline_embedding
    )

    trajectory_columns = (
        inline_positions * crossline_positions
    )

    maximum_rank = min(
        trajectory_rows,
        trajectory_columns,
    )

    rank = min(int(rank), maximum_rank)

    # --------------------------------------------------------------
    # Helper functions
    # --------------------------------------------------------------
    def window_starts(n, length, overlap):
        """Return complete overlapping-window start locations."""
        if length >= n:
            return [0]

        step = max(
            1,
            int(round(length * (1.0 - overlap))),
        )

        starts = list(
            range(0, n - length + 1, step)
        )

        final_start = n - length

        if starts[-1] != final_start:
            starts.append(final_start)

        return starts

    def positive_hann(length):
        """Create a Hann-like taper with nonzero endpoints."""
        if length <= 1:
            return np.ones(
                length,
                dtype=np.float64,
            )

        return np.hanning(length + 2)[1:-1]

    def make_trajectory_matrix(frequency_slice):
        """
        Convert one 2D spatial frequency slice into a trajectory
        matrix of overlapping flattened spatial patches.
        """
        trajectory = np.empty(
            (
                trajectory_rows,
                trajectory_columns,
            ),
            dtype=np.complex128,
        )

        column = 0

        for ix in range(inline_positions):
            for iy in range(crossline_positions):
                patch = frequency_slice[
                    ix:ix + inline_embedding,
                    iy:iy + crossline_embedding,
                ]

                trajectory[:, column] = patch.reshape(
                    -1,
                    order="C",
                )

                column += 1

        return trajectory

    def trajectory_average(trajectory):
        """
        Map a trajectory matrix back to its 2D spatial slice by
        averaging all overlapping embedded patches.
        """
        reconstructed = np.zeros(
            (inline_window, crossline_window),
            dtype=np.complex128,
        )

        counts = np.zeros(
            (inline_window, crossline_window),
            dtype=np.float64,
        )

        column = 0

        for ix in range(inline_positions):
            for iy in range(crossline_positions):
                patch = trajectory[:, column].reshape(
                    (
                        inline_embedding,
                        crossline_embedding,
                    ),
                    order="C",
                )

                reconstructed[
                    ix:ix + inline_embedding,
                    iy:iy + crossline_embedding,
                ] += patch

                counts[
                    ix:ix + inline_embedding,
                    iy:iy + crossline_embedding,
                ] += 1.0

                column += 1

        return reconstructed / counts

    def reduce_rank(trajectory):
        """Apply hard or damped singular-value rank reduction."""
        U, singular_values, Vh = np.linalg.svd(
            trajectory,
            full_matrices=False,
        )

        retained_values = singular_values[:rank].copy()

        if (
            damping_exponent is not None
            and rank < singular_values.size
        ):
            first_rejected = singular_values[rank]

            safe_values = np.maximum(
                retained_values,
                np.finfo(float).eps,
            )

            damping_weights = 1.0 - (
                first_rejected / safe_values
            ) ** damping_exponent

            damping_weights = np.clip(
                damping_weights,
                0.0,
                1.0,
            )

            retained_values *= damping_weights

        return (
            U[:, :rank] * retained_values
        ) @ Vh[:rank, :]

    def mssa_reconstruct(frequency_slice):
        """
        Apply one or more MSSA/Cadzow iterations to a local 2D
        frequency slice.
        """
        estimate = frequency_slice.astype(
            np.complex128,
            copy=True,
        )

        for _ in range(n_iterations):
            trajectory = make_trajectory_matrix(
                estimate
            )

            low_rank_trajectory = reduce_rank(
                trajectory
            )

            estimate = trajectory_average(
                low_rank_trajectory
            )

        return estimate

    # --------------------------------------------------------------
    # Window locations
    # --------------------------------------------------------------
    inline_starts = window_starts(
        n_inline,
        inline_window,
        inline_overlap,
    )

    crossline_starts = window_starts(
        n_crossline,
        crossline_window,
        crossline_overlap,
    )

    time_starts = window_starts(
        n_samples,
        window_samples,
        time_overlap,
    )

    # --------------------------------------------------------------
    # Analysis and synthesis tapers
    # --------------------------------------------------------------
    inline_taper = positive_hann(
        inline_window
    )

    crossline_taper = positive_hann(
        crossline_window
    )

    # The time taper is applied during both FFT analysis and synthesis.
    time_taper = np.sqrt(
        positive_hann(window_samples)
    )

    # The spatial taper is used when blending processed spatial tiles.
    spatial_taper = (
        inline_taper[:, None]
        * crossline_taper[None, :]
    )

    synthesis_taper = (
        spatial_taper[:, :, None]
        * time_taper[None, None, :]
    )

    weight_taper = (
        spatial_taper[:, :, None]
        * time_taper[None, None, :]**2
    )

    # --------------------------------------------------------------
    # Frequency selection
    # --------------------------------------------------------------
    frequencies = np.fft.rfftfreq(
        window_samples,
        d=dt if dt is not None else 1.0,
    )

    selected_frequencies = np.ones(
        frequencies.size,
        dtype=bool,
    )

    # Leave DC unchanged.
    selected_frequencies[0] = False

    if fmin is not None:
        selected_frequencies &= frequencies >= fmin

    if fmax is not None:
        selected_frequencies &= frequencies <= fmax

    selected_indices = np.flatnonzero(
        selected_frequencies
    )

    # --------------------------------------------------------------
    # Main processing loop
    # --------------------------------------------------------------
    filtered_sum = np.zeros_like(
        working_data,
        dtype=output_dtype,
    )

    weight_sum = np.zeros_like(
        working_data,
        dtype=np.float64,
    )

    for inline_start in inline_starts:
        inline_stop = inline_start + inline_window

        for crossline_start in crossline_starts:
            crossline_stop = (
                crossline_start + crossline_window
            )

            for time_start in time_starts:
                time_stop = time_start + window_samples

                tile = working_data[
                    inline_start:inline_stop,
                    crossline_start:crossline_stop,
                    time_start:time_stop,
                ]

                # Apply only the temporal analysis taper. Applying a
                # spatial taper before MSSA would alter its low-rank
                # plane-wave structure.
                analysis_tile = (
                    tile * time_taper[None, None, :]
                )

                spectrum = np.fft.rfft(
                    analysis_tile,
                    axis=2,
                )

                filtered_spectrum = spectrum.copy()

                for frequency_index in selected_indices:
                    original_slice = spectrum[
                        :,
                        :,
                        frequency_index,
                    ]

                    mssa_slice = mssa_reconstruct(
                        original_slice
                    )

                    filtered_spectrum[
                        :,
                        :,
                        frequency_index,
                    ] = (
                        (1.0 - strength) * original_slice
                        + strength * mssa_slice
                    )

                filtered_tile = np.fft.irfft(
                    filtered_spectrum,
                    n=window_samples,
                    axis=2,
                )

                filtered_tile *= synthesis_taper

                filtered_sum[
                    inline_start:inline_stop,
                    crossline_start:crossline_stop,
                    time_start:time_stop,
                ] += filtered_tile

                weight_sum[
                    inline_start:inline_stop,
                    crossline_start:crossline_stop,
                    time_start:time_stop,
                ] += weight_taper

    valid_weights = weight_sum > np.finfo(float).eps

    filtered = np.zeros_like(
        filtered_sum,
        dtype=output_dtype,
    )

    filtered[valid_weights] = (
        filtered_sum[valid_weights]
        / weight_sum[valid_weights]
    )

    filtered[~valid_weights] = working_data[~valid_weights]

    noise = working_data - filtered

    filtered[invalid_mask] = np.nan
    noise[invalid_mask] = np.nan

    return filtered, noise



def _smooth_frequency(a: np.ndarray, n_bins: int) -> np.ndarray:
    """Smooth the last dimension with a normalized boxcar."""
    if n_bins <= 1:
        return a.copy()
    if n_bins % 2 == 0:
        n_bins += 1
    pad = n_bins // 2
    kernel = np.ones(n_bins, dtype=float) / n_bins
    padded = np.pad(a, [(0, 0)] * (a.ndim - 1) + [(pad, pad)], mode="edge")
    return np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="valid"), -1, padded)


def _four_corner_taper(
    frequencies: np.ndarray, corners: tuple[float, float, float, float]
) -> np.ndarray:
    """Return a cosine-tapered bandpass defined by f1, f2, f3, f4."""
    f1, f2, f3, f4 = corners
    if not (0 <= f1 <= f2 < f3 <= f4 <= frequencies[-1]):
        raise ValueError("Require 0 <= f1 <= f2 < f3 <= f4 <= Nyquist")

    taper = np.zeros_like(frequencies)
    taper[(frequencies >= f2) & (frequencies <= f3)] = 1.0

    if f2 > f1:
        m = (frequencies > f1) & (frequencies < f2)
        u = (frequencies[m] - f1) / (f2 - f1)
        taper[m] = 0.5 - 0.5 * np.cos(np.pi * u)
    if f4 > f3:
        m = (frequencies > f3) & (frequencies < f4)
        u = (frequencies[m] - f3) / (f4 - f3)
        taper[m] = 0.5 + 0.5 * np.cos(np.pi * u)
    return taper


def time_variant_spectral_whitening(
    data: np.ndarray,
    dt: float,
    *,
    window_length: float = 0.400,
    overlap: float = 0.75,
    freq_corners: tuple[float, float, float, float] = (5.0, 10.0, 80.0, 100.0),
    smooth_hz: float = 8.0,
    alpha: float = 0.5,
    water_level: float = 0.05,
    max_gain: float | None = 5.0,
    mode: str = "ensemble",
) -> np.ndarray:
    """Apply stabilized, partial time-variant spectral whitening.

    Parameters
    ----------
    data : array_like
        One trace ``(n_samples,)`` or traces ``(n_traces, n_samples)``.
    dt : float
        Sample interval in seconds.
    window_length : float, default 0.400
        STFT analysis-window length in seconds.
    overlap : float, default 0.75
        Fractional overlap between adjacent windows, in ``[0, 1)``.
    freq_corners : (f1, f2, f3, f4)
        Cosine-tapered whitening passband in Hz. The response rises from zero
        at f1 to one at f2, stays flat through f3, and falls to zero at f4.
    smooth_hz : float, default 8.0
        Width of frequency smoothing applied to the amplitude spectrum.
    alpha : float, default 0.5
        Whitening strength: 0 retains the spectral envelope and 1 requests
        full flattening within the selected band.
    water_level : float, default 0.05
        Stabilizer as a fraction of the median passband amplitude.
    max_gain : float or None, default 5.0
        Maximum whitening gain. Use None for no explicit cap.
    mode : {'ensemble', 'trace'}, default 'ensemble'
        ``ensemble`` derives one operator per time window from the median
        spectrum of all traces, preserving lateral amplitude relationships
        better. ``trace`` derives a separate operator for every trace.

    Returns
    -------
    ndarray
        Whitened data with the same shape as the input.

    Notes
    -----
    The local operator is

        G(t, f) = [A_ref(t) / (A_smooth(t, f) + epsilon(t))] ** alpha

    followed by a four-corner frequency taper and optional gain clipping.
    Square-root Hann analysis/synthesis windows and explicit overlap
    normalization are used during reconstruction.
    """
    x = np.asarray(data)
    was_1d = x.ndim == 1
    if was_1d:
        x = x[None, :]
    if x.ndim != 2:
        raise ValueError("data must be 1-D or 2-D with time on the final axis")
    if x.shape[1] < 2 or not np.all(np.isfinite(x)):
        raise ValueError("data must contain at least two finite samples per trace")
    if dt <= 0 or window_length <= 0 or smooth_hz <= 0:
        raise ValueError("dt, window_length, and smooth_hz must be positive")
    if not 0 <= overlap < 1:
        raise ValueError("overlap must be in [0, 1)")
    if not 0 <= alpha <= 1:
        raise ValueError("alpha must be in [0, 1]")
    if water_level < 0 or (max_gain is not None and max_gain <= 0):
        raise ValueError("water_level must be nonnegative and max_gain positive")
    if mode not in {"ensemble", "trace"}:
        raise ValueError("mode must be 'ensemble' or 'trace'")

    original_dtype = x.dtype
    x = x.astype(np.float64, copy=False)
    n_traces, n_samples = x.shape
    nwin = max(4, int(round(window_length / dt)))
    if nwin % 2:
        nwin += 1
    hop = max(1, int(round(nwin * (1.0 - overlap))))
    nfft = nwin

    analysis = np.sqrt(np.hanning(nwin))
    frequencies = np.fft.rfftfreq(nfft, dt)
    band_taper = _four_corner_taper(frequencies, freq_corners)
    passband = band_taper >= 0.999
    if not np.any(passband):
        raise ValueError("freq_corners produce no FFT bins in the full passband")
    df = frequencies[1] - frequencies[0]
    smooth_bins = max(1, int(round(smooth_hz / df)))

    # Padding centers the first/last windows on the first/last data samples.
    half = nwin // 2
    xpad = np.pad(x, ((0, 0), (half, half)), mode="reflect")
    out = np.zeros_like(xpad)
    weight = np.zeros(xpad.shape[1], dtype=float)

    starts = list(range(0, xpad.shape[1] - nwin + 1, hop))
    last = xpad.shape[1] - nwin
    if starts[-1] != last:
        starts.append(last)

    tiny = np.finfo(float).tiny
    for start in starts:
        frame = xpad[:, start : start + nwin] * analysis
        spectrum = np.fft.rfft(frame, n=nfft, axis=-1)
        amplitude = np.abs(spectrum)

        if mode == "ensemble":
            envelope = np.median(amplitude, axis=0, keepdims=True)
        else:
            envelope = amplitude
        envelope = _smooth_frequency(envelope, smooth_bins)

        reference = np.median(envelope[..., passband], axis=-1, keepdims=True)
        epsilon = water_level * reference
        gain = (reference / np.maximum(envelope + epsilon, tiny)) ** alpha
        if max_gain is not None:
            gain = np.minimum(gain, max_gain)
        gain *= band_taper

        reconstructed = np.fft.irfft(spectrum * gain, n=nfft, axis=-1)
        out[:, start : start + nwin] += reconstructed * analysis
        weight[start : start + nwin] += analysis**2

    out /= np.maximum(weight, tiny)[None, :]
    out = out[:, half : half + n_samples]

    if np.issubdtype(original_dtype, np.floating):
        out = out.astype(original_dtype, copy=False)
    return out[0] if was_1d else out



def seismic_sharpen(data, sigma=2.0, amount=1.0):
    """
    Sharpen seismic data along the time/sample axis only.

    Parameters
    ----------
    data : ndarray
        Seismic data of arbitrary dimensionality.
        The last dimension must be time/sample.

        Examples:
            (nsamples,)
            (ntraces, nsamples)
            (nilines, nxlines, nsamples)

    sigma : float
        Standard deviation of Gaussian smoothing kernel, in samples.

    amount : float
        Strength of sharpening.
        0.0 = no sharpening
        1.0 = normal sharpening
        >1.0 = stronger sharpening

    Returns
    -------
    sharpened : ndarray
        Sharpened seismic data.
    """

    data = np.asarray(data)

    smooth = gaussian_filter1d(
        data,
        sigma=sigma,
        axis=-1,
        mode="reflect"
    )

    highpass = data - smooth

    sharpened = data + amount * highpass

    return sharpened
