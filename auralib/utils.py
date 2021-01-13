"""
auralib module for various utility functions.

Author:   Wes Hamlyn
Created:  16-Aug-2016
Last Mod: 17-Aug-2016

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.interpolate import interp1d
import win32clipboard as cb
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms



def iseven(number):
    """
    Convenience function to determine if a value is even.  Returns boolean
    value or numpy array (depending on input)
    """
    
    result = (number % 2) == 0

    return result
    

def isodd(number):
    """
    Convenience function to determine if a value is odd.  Returns boolean
    value or numpy array (depending on input)
    """
    
    result = (number % 2) == 1

    return result


def islive(tdata, tol=0.0):
    """
    Convenience function to verify if a seismic trace is live or dead by 
    examining amplitudes. If the 

    Parameters
    ----------
    tdata : list, tuple, numpy array

    Returns
    -------
    flag : Boolean, True if live trace, False if dead trace
    """
    
    meanval = np.mean(np.abs(tdata))
    
    if meanval > tol:
        flag = True
    else:
        flag = False
    
    return flag


def digitize_poly():
    """
    Convenience function that uses the matplotlib ginput function to capture
    polygon points and return as two separate numpy arrays for the x and y
    click coordinates. 
    """
    
    pts = plt.ginput(n=-1, timeout=60, show_clicks=True,
                     mouse_add=1, mouse_pop=2, mouse_stop=3)
    
    pts = np.array(pts)
    xpoly = pts[:, 0]
    ypoly = pts[:, 1]
    
    return xpoly, ypoly



def inpoly(datax, datay, polyx, polyy):
    """
    Determine if points are inside a polygon
    
    datax, datay = arrays of data coordinates
    polyx, polyy = arrays of polygon verticies
    
    Returns a list of indicies corresponding to points that are inside the 
    polygon.
    
    Note:   There are some issues with this function for points that
            lie along edges and verticies of the polygon.  These specific cases
            need to be tested and debugged.  In general, this function will
            work well enough for manually digitized polygons.
    """
    
    
    p_xmin = polyx.min()
    p_xmax = polyx.max() 
    p_ymin = polyy.min() 
    p_ymax = polyy.max()
    
    nvert = len(polyx)
    polyx = np.hstack([polyx, polyx[0]])
    polyy = np.hstack([polyy, polyy[0]])
    
    idx = []
    for i in range(0, len(datax)):
        count = 0
        
        if (datax[i] <= p_xmin) | (datax[i] >= p_xmax) | \
           (datay[i] <= p_ymin) | (datay[i] >= p_ymax):
            next
                
        for j in range(0, nvert):
            val = (polyx[j+1]-polyx[j]) * (datay[i]-polyy[j]) / \
                    (polyy[j+1]-polyy[j]) + polyx[j]
            if ((polyy[j] > datay[i]) != (polyy[j+1] > datay[i])) & \
                (datax[i] < val):
                count += 1
        
        if isodd(count):
            idx.append(i)
    
    return idx



def confidence_ellipse(x, y, ax, n_std=3.0, ec='k', facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, ec=ec, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
        
    return ax.add_patch(ellipse)


def confidence_ellipse2(mean, cov, ax, n_std=3.0, ec='k', facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    mean : array-like, shape (2, )
        matrix of mean values
    
    cov : array-like, shape (n, )
        covariance matrix

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, ec=ec, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = mean[0]

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = mean[1]

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def get_dist(ax):
    """
    Convenience function to get a distance from a map and plot line in axes
    """
    
    pts = plt.ginput(n=2)
    
    x1 = pts[0][0]
    x2 = pts[1][0]
    y1 = pts[0][1]
    y2 = pts[1][1]
    
    dist = np.sqrt( (x2-x1)**2 + (y2-y1)**2 )
    dist_txt = '%.1f m' % dist
    
    ax.plot([x1, x2], [y1, y2], 'k-')
    ax.plot([x1, x2], [y1, y2], 'ko')
    
    
    ax.text(x2+10, y2, dist_txt)
    
    #plt.draw()
    
    return x1, x2, y1, y2, dist



def nextpow2(value):
    """
    Returns the next power of 2 integer larger than value.
    """
    
    x = value - 1
    npow2 = 1 << x.bit_length()
    
    return npow2



def padzeros(data, pad_length):
    """
    Pads a 1D data array with zeros, typically for input to FFT.
    
    data = input data vector
    pad_length = desired output length
    
    data_pad = data vector padded with zeros
    pad_start = number of samples added to beginning of vector
    pad_end = number of samples added to end of vector
    """
    
    nsamp = len(data)
    padding = pad_length - nsamp
    
    pad_start = int(np.floor(padding*0.5))
    pad_end = int(np.ceil(padding*0.5))

    data_pad = np.pad(data, [pad_start, pad_end], 
                     mode='linear_ramp', end_values=0)
    
    return data_pad, pad_start, pad_end



def clip_seis_amp(tdata, min_clip='None', max_clip='None'):
    """
    Apply amplitude clipping to seismic data traces.
    """
    
    tdata = np.array(tdata)
    
    if min_clip == 'None':
        min_clip = np.min(tdata)
        
    if max_clip == 'None':
        max_clip = np.max(tdata)
    
    idx_min = np.nonzero(tdata <= min_clip)
    idx_max = np.nonzero(tdata >= max_clip)
    
    tdata[idx_min] = min_clip
    tdata[idx_max] = max_clip
    
    return tdata



def smooth_log(log, smooth_window):
    """
    Function to smooth logs using a running average.
    """
	
    # make sure smoothing window is an odd number    
    if smooth_window%2 == 0: # 0 is even number; 1 is odd number
        smooth_window = smooth_window + 1
    
    # pad top and bottom of log with duplicate values to avoid edge effects
    s1 = log[0]
    s2 = log[-1]
    padtop = s1*np.ones(int(smooth_window/2))
    padend = s2*np.ones(int(smooth_window/2))
    logpad = np.hstack([padtop, log, padend])
    
    # design boxcar operator and smooth the log
    filt = np.ones(smooth_window)/smooth_window
    logsm = sp.convolve(logpad, filt, 'valid')

    return logsm


def resample(data, md, md_new, null=-999.25):
    """
    Function to resample a 1D data array containing np.nan (e.g. well log data)
    """
    
    # find nan samples
    idx = np.nonzero(np.isnan(data))
    
    # replace np.nans with a null value
    data_prep = data*1
    data_prep[idx] = null
    
    # make a psuedo-boolean array of zeros ond ones to indicate nan vs live 
    # samples
    bool_null = np.zeros_like(data)
    bool_null[idx] = 1
    
    # interpolate the prepared data array using cubic spline
    f = interp1d(md, data_prep, kind='cubic', bounds_error=False, fill_value=np.nan)
    data_resamp = f(md_new)
    
    # interpolate the pseudo-boolean array using nearest interpoation
    f1 = interp1d(md, bool_null, kind='nearest', bounds_error=False, fill_value=0)
    bool_null2 = f1(md_new)
    
    # replace null values with np.nan
    idx = np.nonzero(bool_null2==1)
    data_resamp[idx] = np.nan
    
    return data_resamp


def get_clipboard():
    """
    Convenience function to copy contents of clipboard to a numpy array. First
    written to copy data from and Excel spreadsheet via a copy/paste type of
    operation.
    
    Assumes all data are numeric and floating point (could easily be improved)
    
    Last modified: 17-May-2019
    """
    
    cb.OpenClipboard()
    data = cb.GetClipboardData()
    cb.CloseClipboard()
    
    tmp = []
    for each in data.split('\r\n'):
        tmp.append(each.split('\t'))
    
    tmp.pop(-1)    
    data2 = np.array(tmp, dtype='float')
    return data2


def split(text):
    """
    Convenience function to parse a string into a list using two or more spaces
    as the delimiter.
    
    Parameters
    ----------
    text : string
        Contains text to be parsed

    Returns
    -------
    textout : list
        list of strings containing the parsed input text
    """
    
    textout = [s.strip() for s in text.split('  ') if s.strip()]
    
    return textout


def calc_rms(x):
    """
    Calculate root mean square value of a numeric array

    Parameters
    ----------
    x : list, tuple, numpy array
        Contains values to compute RMS value of

    Returns
    -------
    rms : float
        RMS value of the input array x

    """
    
    n = len(x)
    rms = np.sqrt(np.sum(x**2)/n)
    return rms


def calc_nrms(x1, x2):
    """
    Calculate the normalized RMS error between two numeric arrays
    
    Parameters
    ----------
    x1, x2 : lists, tuples, numpy arrays
             Contains values to compute NRMS value of

    Returns
    -------
    nrms : float
        NRMS value of the input arrays x1, x2

    """
    
    rms1 = calc_rms(x1)
    rms2 = calc_rms(x2)
    rms12 = calc_rms(x1-x2)
    
    nrms = 2 * rms12/(rms1+rms2)
    
    return nrms


def predictability(x1, x2, as_percent=True):
    """
    Calculate the predictability value between two numeric arrays
    
    Parameters
    ----------
    x1, x2 : lists, tuples, numpy arrays
             Contains values to compute predictability of

    Returns
    -------
    pred : float
        predictability value of the input arrays x1, x2

    """
    
    # compute autocorrelation at lag zero
    ac1 = np.sum(x1*x1)
    ac2 = np.sum(x2*x2)
    
    # compute crosscorrelation at lag zero
    xc12 = np.sum(x1*x2)
    
    # compute predictability
    pred = (xc12*xc12)/(ac1*ac2)
    
    # express as percentage
    if as_percent:
        pred = pred*100.0
    
    return pred
    

