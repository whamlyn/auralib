"""
AuraQI module for creating synthetic seismograms.

Author:   Wes Hamlyn
Created:  16-Aug-2016
Last Mod: 17-Aug-2016
"""

import numpy as np
from scipy.interpolate import interp1d
from copy import deepcopy


def d2t(depth, data, td_depth, td_twt, target_twt):
    """
    Convenience function for converting a well log from depth to time given
    time-depth pairs.  This function both converts from depth to time and
    resamples the log data to a regular sampling rate in time.
    """
    
    f1 = interp1d(td_depth, td_twt, kind='linear', fill_value='extrapolate')
    twt_logstep = f1(depth)
        
    f2 = interp1d(twt_logstep, data, kind='linear', fill_value='extrapolate')
    data_twt = f2(target_twt)
    
    return data_twt



def t2d(twt, data, td_depth, td_twt, target_depth):
    """
    Convenience function for converting a well log from depth to time given
    time-depth pairs.  This function both converts from time to depth and
    resamples the log data to a regular sampling rate in depth.
    """
    
    f1 = interp1d(td_twt, td_depth, kind='linear', fill_value='extrapolate')
    depth_timestep = f1(twt)
        
    f2 = interp1d(depth_timestep, data, kind='linear', fill_value='extrapolate')
    data_depth = f2(target_depth)
    
    return data_depth



def backus_length(f_dom, Vs0, mode=1):
    """
    Returns the backus average length for scattering and transmission limits.
    
    L = backus_length(f_dom, Vs0, mode=1)
    
    if mode==1:
        L = Ls = scattering limit (reflection seismic)
        
    if mode==2:
        L = Lt = transmission limit (ray tracing)

    """
    
    Vs0 = float(Vs0)
    f_dom = float(f_dom)
    
    Ls = Vs0/(3.0*f_dom)
    Lt = 2.0*Vs0/f_dom
    
    if mode == 1:
        L = Ls
            
    if mode == 2:
        L = Lt
    
    return L
        


def backus_average(vp, vs, rho, nsamp):
    """
    Computes the backus average from vp, vs, and density logs.
    
    Vpv, Vph, Vsv, Vsh, rho = backus_average(vp, vs, rho, nsamp)
    """
    
    # make sure number of samples is odd
    nsamp = np.round(nsamp)
    if nsamp == 0:
        nsamp = 1
    if nsamp%2.0 == 0:
        nsamp = int(np.ceil(nsamp))
    
    print('Backus averaging over %i samples' % nsamp)
    
    filt = np.ones(nsamp)
    
    # pad top and bottom of input logs to avoid filtering edge effects.  Do
    # this by duplicating the first and last sample in the vp, vs, and rho logs
    # by ~1/2 the Backus filtering length
    pad_length = int((nsamp-1)/2)
        
    vp = np.hstack([np.ones(pad_length)*vp[0], vp, np.ones(pad_length)*vp[-1]])
    vs = np.hstack([np.ones(pad_length)*vs[0], vs, np.ones(pad_length)*vs[-1]])
    rho = np.hstack([np.ones(pad_length)*rho[0], rho, np.ones(pad_length)*rho[-1]])
    
    mu = rho*(vs**2.0)
    K = rho*(vp**2.0) - 4.0/3.0*mu
    lame = K - 2.0/3.0*mu
    
    
    A1 = (4.0*mu*(lame+mu))/(lame+2.0*mu)
    A2 = (1.0/(lame+2.0*mu))**-1.0
    A3 = (lame/(lame+2.0*mu))**2.0
    
    A1 = np.convolve(filt, A1, mode='same')/nsamp
    A2 = np.convolve(filt, A2, mode='same')/nsamp
    A3 = np.convolve(filt, A3, mode='same')/nsamp
    A = A1 + A2*A3
    
    B1 = (2.0*mu*lame)/(lame+2.0*mu)
    B2 = (1.0/(lame+2.0*mu))**-1.0
    B3 = (lame/(lame+2.0*mu))**2.0
    
    B1 = np.convolve(filt, B1, mode='same')/nsamp
    B2 = np.convolve(filt, B2, mode='same')/nsamp
    B3 = np.convolve(filt, B3, mode='same')/nsamp
    B = B1 + B2*B3
    
    C = (1.0/(lame+2.0*mu))**-1.0
    C = np.convolve(filt, C, mode='same')/nsamp
    
    F1 = (1.0/(lame+2.0*mu))**-1.0
    F2 = lame/(lame+2.0*mu)
    
    F1 = np.convolve(filt, F1, mode='same')/nsamp
    F2 = np.convolve(filt, F2, mode='same')/nsamp
    F = F1*F2
    
    D = (1.0/mu)**-1.0
    D = np.convolve(filt, D, mode='same')/nsamp
    
    M = np.convolve(filt, mu, mode='same')/nsamp
    
    rho_b = np.convolve(filt, rho, mode='same')/nsamp
    
    Vpv_b = np.sqrt(C/rho_b)
    Vph_b = np.sqrt(A/rho_b)
    Vsv_b = np.sqrt(D/rho_b)
    Vsh_b = np.sqrt(M/rho_b)
    
    # remove padding added to avoid edge effects
    Vpv_b = Vpv_b[pad_length:-pad_length]
    Vph_b = Vph_b[pad_length:-pad_length]
    Vsv_b = Vsv_b[pad_length:-pad_length]
    Vsh_b = Vsh_b[pad_length:-pad_length]
    rho_b = rho_b[pad_length:-pad_length]
    
    return Vpv_b, Vph_b, Vsv_b, Vsh_b, rho_b
    
    
def plot_wigva(ax, tdata, t, baseline=0, peak=True, trough=True, 
          lcolor='k', pcolor=[0.2, 0.2, 1.0], tcolor=[1.0, 0.2, 0.2]):
    """
    Python function to plot synthetic wiggle traces with variable area fill.
    """
    
    from scipy.interpolate import interp1d
    
    # to make the peak/trough fill more accurate near zero-crossings, increase
    # the sampling of synthetic traces by a factor of 5
    nsamp = len(tdata)
    t2 = np.linspace(t.min(), t.max(), nsamp*5)
    
    # interpolate the new trace
    #t2 = np.arange(t.min(), t.max(), dt2)
    f = interp1d(t, tdata, kind='linear')
    trc2 = f(t2)
    
    if (peak == True) | (trough == True):
        for i in range(0, 1):
            # plot peak fill
            if peak==True:
                ax.fill_betweenx(t2, baseline, trc2, where=trc2>=baseline, 
                                 alpha=0.6, facecolor=pcolor,  edgecolor=pcolor)
            
            # plot trough fill
            if trough==True:
                ax.fill_betweenx(t2, baseline, trc2, where=trc2<=baseline, 
                                 alpha=0.6, facecolor=tcolor,  edgecolor=tcolor)
    
        ax.plot(trc2, t2, lcolor)