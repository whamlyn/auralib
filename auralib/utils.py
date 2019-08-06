"""
auralib module for various utility functions.

Author:   Wes Hamlyn
Created:  16-Aug-2016
Last Mod: 17-Aug-2016

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import win32clipboard as cb


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
    Returns the next power of 2 larger than value.
    """
    
    npow2 = 1<<(value-1).bit_length()
    
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