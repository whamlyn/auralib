"""
auralib module to define various types of wavelets.

Author:   Wes Hamlyn
Created:   3-Jan-2016
Last Mod: 17-Aug-2016

"""

import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
import auralib as aura


def ricker(cfreq, phase, dt, wvlt_length):
    """
    Calculate a ricker wavelet
    
    Usage:
    ------
    t, wvlt = wvlt_ricker(cfreq, phase, dt, wvlt_length)
    
    cfreq: central frequency of wavelet in Hz
    phase: wavelet phase in degrees
    dt: sample rate in seconds
    wvlt_length: length of wavelet in seconds
    """
    cfreq = float(cfreq)
    phase = float(phase)
    dt = float(dt)
    wvlt_length = float(wvlt_length)
    
    if aura.utils.iseven(wvlt_length/dt):
        wvlt_length += dt
        
    nsamp = int(wvlt_length/dt)
    t_max = wvlt_length*0.5    
    t = np.linspace(-t_max, t_max, nsamp)
    
    #t = np.linspace(-wvlt_length/2, (wvlt_length-dt)/2, wvlt_length/dt)
    wvlt = (1.0 - 2.0*(np.pi**2)*(cfreq**2)*(t**2)) * np.exp(-(np.pi**2)*(cfreq**2)*(t**2))
    
    if phase != 0:
        phase = phase*np.pi/180.0
        wvlth = signal.hilbert(wvlt)
        wvlth = np.imag(wvlth)
        wvlt = np.cos(phase)*wvlt - np.sin(phase)*wvlth
    
    return t, wvlt


def wvlt_bpass(f1, f2, f3, f4, phase, dt, wvlt_length):
    """
    Calculate a trapezoidal bandpass wavelet
    
    Usage:
    ------
    t, wvlt = wvlt_bpass(f1, f2, f3, f4, phase, dt, wvlt_length)
    
    f1: Low truncation frequency of wavelet in Hz
    f2: Low cut frequency of wavelet in Hz
    f3: High cut frequency of wavelet in Hz
    f4: High truncation frequency of wavelet in Hz
    phase: wavelet phase in degrees
    dt: sample rate in seconds
    wvlt_length: length of wavelet in seconds
    """
    
    from numpy.fft import fft, ifft, fftfreq, fftshift, ifftshift
    
    f1 = float(f1)
    f2 = float(f2)
    f3 = float(f3)
    f4 = float(f4)
    phase = float(phase)
    dt = float(dt)
    wvlt_length = float(wvlt_length)
    
    nsamp = int(wvlt_length/dt + 1)
    
    
    freq = fftfreq(nsamp, dt)
    freq = fftshift(freq)
    
    # Calculate slope and y-int for low frequency ramp
    M1 = 1/(f2-f1)
    b1 = -M1*f1
    
    # Calculate slope and y-int for high frequency ramp
    M2 = -1/(f4-f3)
    b2 = -M2*f4
    
    # Build initial frequency and filter arrays
    freq = fftfreq(nsamp, dt)
    freq = fftshift(freq)
    filt = np.zeros(nsamp)
    
    # Build LF ramp
    idx = np.nonzero((np.abs(freq)>=f1) & (np.abs(freq)<f2))
    filt[idx] = M1*np.abs(freq)[idx]+b1
    
    # Build central filter flat
    idx = np.nonzero((np.abs(freq)>=f2) & (np.abs(freq)<=f3))
    filt[idx] = 1.0
    
    # Build HF ramp
    idx = np.nonzero((np.abs(freq)>f3) & (np.abs(freq)<=f4))
    filt[idx] = M2*np.abs(freq)[idx]+b2
    
    # Unshift the frequencies and convert filter to fourier coefficients
    filt2 = ifftshift(filt)
    Af = filt2*np.exp(np.zeros(filt2.shape)*1j)
    
    # Convert filter to time-domain wavelet
    wvlt = fftshift(ifft(Af))
    wvlt = np.real(wvlt)
    wvlt = wvlt/np.max(np.abs(wvlt)) # normalize wavelet by peak amplitude
    
    # Generate array of wavelet times
    t = np.linspace(-wvlt_length*0.5, wvlt_length*0.5, nsamp)
    
    
    # Apply phase rotation if desired
    if phase != 0:
        phase = phase*np.pi/180.0
        wvlth = signal.hilbert(wvlt)
        wvlth = np.imag(wvlth)
        wvlt = np.cos(phase)*wvlt - np.sin(phase)*wvlth
    
    return t, wvlt


def extract_wavelet(a, dt, wlen_t, sm_freq=10.0, phase=0.0):
    """
    Extract wavelet using frequency domain method. Currently works only on a
    single trace but can be modified.
    """
    
    # calculate frequency axis labels
    freq = np.fft.fftshift(np.fft.fftfreq(len(a), d=dt))
    df = np.median(np.diff(freq))
    
    # transform trace amplitudes to fourier domain and make amplitude spectrum
    af = np.fft.fftshift(np.fft.fft(a))
    aspec = np.abs(af)
    
    # apply smoother to amplitude spectrum so that we're not fitting a spline
    # to a noisy spectrum
    smlen = int(np.round(sm_freq/df))
    aspec_sm = np.convolve(np.ones(smlen), aspec, mode='same')/smlen
    
    # calculate the wavelet length to the nearest sample
    wlen_s = int(np.round(wlen_t/dt))
    
    # resample the smoothed amplitude spectrum to the wavelet length
    f = interp1d(freq, aspec_sm, kind='linear')
    freq2 = np.linspace(freq.min(), freq.max(), wlen_s)
    aspec2 = f(freq2)
    pspec2 = np.zeros(len(aspec2))
    
    # option: set zero hertz amplitude component to zero if desired
    #aspec2[freq2==0.0] = 0.0
    
    # inverse fourier transform the resampled smoothed spectrum to the time
    # domain.
    # Note: For some reason I don't totally understand I need to apply an extra
    # ifftshift to get the wavelet centered on zero and not wrapping around.
    # need to investigate this
    wf = aspec2*np.exp(1j*pspec2)
    wa = np.fft.ifft(np.fft.ifftshift(wf)).real
    wa = np.fft.fftshift(wa)
    
    # make wavelet time sample array
    wt = np.arange(len(wa))*dt
    wt = wt - np.mean(wt)
    
    if phase!=0.0:
        wa = aura.filt.phaserot(wa, phase)
    
    return wt, wa

