import numpy as np
import matplotlib.pyplot as plt
import auralib as aura
from numpy.fft import fftfreq, fft, ifft, fftshift, ifftshift
from scipy.interpolate import interp1d
import scipy as sp


def get_traces_for_matching_filter(basefile, monfile, step):
    buf1 = aura.segy.Segy(basefile)
    buf2 = aura.segy.Segy(monfile)
    
    dt = buf1.bhead['samp_rate']*1e-6
    nsamp = buf1.bhead['num_samp']
    twtr = np.arange(nsamp)*dt
    
    
    tdata1r = []
    tdata2r = []
    trcnum = []
    for i in np.arange(0, buf1.num_traces, step):
        
        tmp1 = buf1.read_tdata(i)
        
        if np.mean(np.abs(tmp1)) > 0.0:
            
            tmp2 = buf2.read_tdata(i)
            
            if np.mean(np.abs(tmp2)) > 0.0:
                tdata1r.append(tmp1)
                tdata2r.append(tmp2)
                trcnum.append(i)
    
    tdata1r = np.array(tdata1r)
    tdata2r = np.array(tdata2r)
    trcnum = np.array(trcnum)
    
    if len(trcnum) > 0:
        print('Found %i live traces. Proceeding to next step...' %
              (len(trcnum)))
        return trcnum, twtr, tdata1r, tdata2r
    
    else:
        print('Failed to find live traces. Terminating execution...')


def calc_ampmatch_operator(tdata1, tdata2, twt):
    
    dt = np.mean(np.diff(twt))
    nsamp = len(twt)
    
    freq = fftshift(fftfreq(nsamp, dt))
    tdata1f = fftshift(fft(tdata1, axis=1))
    tdata2f = fftshift(fft(tdata2, axis=1))
    
    aspec1 = np.abs(tdata1f)
    aspec2 = np.abs(tdata2f)
    
    aspec1_avg = np.mean(aspec1, axis=0)
    aspec2_avg = np.mean(aspec2, axis=0)

    
    aspec_op_raw = aspec1_avg / aspec2_avg 
    f1 = 5.0
    f2 = 115.0
    fidx = np.nonzero((np.abs(freq)>=f1) & (np.abs(freq)<=f2))
    aspec_op = np.ones(len(freq))
    aspec_op[fidx] = aspec_op_raw[fidx]
    
    return freq, aspec_op


def save_ampmatch_operator(opfile, freq, aspec_op):
    
    with open(opfile, 'w') as fd:
        for i in range(len(freq)):
            txt = '%f,%f\n' % (freq[i], aspec_op[i])
            fd.write(txt)



def calc_match_filter(d1, d2, npts=-1):
    """
    Calcualte least squares matching filter to correct for amplitude and
    phase differences.
    
    Inputs:
        d1 = trace from master survey
        d2 = trace from slave survey
        npts = number of samples in matching filter
    
    Outputs:
        a = matching filter operator (time-domain)
    """
    
    if npts == -1:
        npts = len(d1)
    
    # build toeplitz matrix of slave survey trace
    r0 = np.zeros(npts)
    r0[0] = d2[0]
    d2_pad = np.hstack([d2, np.zeros(npts-1)])
    D2 = sp.linalg.toeplitz(d2_pad, r0)
    
    # build colum vector of master matrix reflectivities
    D1 = np.hstack([d1, np.zeros(npts-1)])
    D1 = D1.reshape([-1, 1])
    
    # Calcualte least squares match filter
    A = np.dot(np.dot(np.linalg.inv(np.dot(D2.T, D2)), D2.T), D1)
    a = A.flatten()
    
    return a


def apply_match_filter(a, d2):
    """
    Apply a least squares matching filter operator to a data vector
    
    Inputs:
        a = matching filter operator (time-domain)
        d2 = trace from slave survey
    
    Outputs:
        d2m = trace from slave survey after applying matching filter operator
    """
    
    
    npts = len(a)
    A = a.reshape([-1, 1])
    
    # build toeplitz matrix of slave survey trace
    r0 = np.zeros(npts)
    r0[0] = d2[0]
    d2_pad = np.hstack([d2, np.zeros(npts-1)])
    D2 = sp.linalg.toeplitz(d2_pad, r0)
    
    # Apply matching operator to slave survey trace
    D2m = np.dot(D2, A)

    # Remove extra data due to padding operations
    d2m = D2m.flatten()
    d2m = d2m[0:-(npts-1)]
    
    return d2m










