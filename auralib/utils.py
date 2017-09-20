"""
AuraQI module for various utility functions.

Author:   Wes Hamlyn
Created:  16-Aug-2016
Last Mod: 17-Aug-2016

Copyright 2016 Wes Hamlyn

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy as np
import matplotlib.pyplot as plt


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
    

def plot_blocky(ax, data, zdata, linespec='b-', lw=1):
    """
    Convenience function for plotting a blocky log.
    
    Ensure that the zdata log has 1 more sample than the data log.
    """
    
    for i in range(0, len(data)):
        ax.plot([data[i], data[i]], [zdata[i], zdata[i+1]], linespec, lw=lw)
        
    for i in range(1, len(data)):
        ax.plot([data[i-1], data[i]], [zdata[i], zdata[i]], linespec, lw=lw)


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