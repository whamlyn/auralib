"""
AuraQI module containing helper functions for reading data from RokDoc exports.

Author:   Wes Hamlyn
Created:  25-Mar-2016
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

def load_horizon_3d(infile, il_min, il_max, xl_min, xl_max):
    """
    Function to load a 3D horizon exported from RokDoc.  This script expects
    columns for Inline, Crossline, UMTX, UTMY, and Attribute value
    """
    
    import numpy as np
    
    il_min = int(il_min)
    il_max = int(il_max)
    xl_min = int(xl_min)
    xl_max = int(xl_max)
    
    #  calculate number of inlines and crosslines
    num_il = il_max - il_min + 1
    num_xl = xl_max - xl_min + 1
    
    #  build a 2D array to store horizon values
    mapdata = np.ones((num_il, num_xl)) * np.nan
    
    #  read the RokDoc horizon export ASCII file    
    buf = np.loadtxt(infile, skiprows=6)
    
    #  map horizon values from the imported columnar arrangement to a 
    #  2D numpy array corresponding to inline/crosslines
    ili = buf[:, 2] - il_min
    xli = buf[:, 3] - xl_min
    ili = np.array(ili, dtype='int')
    xli = np.array(xli, dtype='int')
    
    zval = np.array(buf[:, 4], dtype='float')

    idx = np.nonzero(buf[:, 4] != -999.25)[0]  # ignore horizon nulls

    mapdata[ili[idx], xli[idx]] = zval[idx]
    
    return mapdata


def load_2d_pdf(infile):
    import numpy as np
    
    fd1 = open(infile, 'r')
    
    for i in range(0, 11):
        buf = fd1.readline()
        if i == 0:
            name = buf.split('"')[1].strip()
        if i == 2:
            xattr = buf.split(':')[-1].strip()
    
        if i == 3:
            xmin = float(buf.split(':')[-1].strip())
    
        if i == 4:
            dx = float(buf.split(':')[-1].strip())
    
        if i == 5:
            nx = float(buf.split(':')[-1].strip())
        
        if i == 7:
            yattr = buf.split(':')[-1].strip()
        
        if i == 8:
            ymin = float(buf.split(':')[-1].strip())
            
        if i == 9:
            dy = float(buf.split(':')[-1].strip())
        
        if i == 10:
            ny = float(buf.split(':')[-1].strip())
    
    fd1.close()
    
    data = np.loadtxt(infile, skiprows=13)
    data = np.flipud(data)
    
    pdf = {}
    pdf['name'] = name
    pdf['xattr'] = xattr
    pdf['xmin'] = xmin
    pdf['dx'] = dx
    pdf['nx'] = nx
    pdf['yattr'] = yattr
    pdf['ymin'] = ymin
    pdf['dy'] = dy
    pdf['ny'] = ny
    pdf['data'] = data
    
    return pdf


def plot_2d_pdf(ax, pdf, cont_inc=-1, color='b'):
    """
    Need to use mplot3d axes.  Make sure to import this using the following
    syntax:
    
        from mpl_toolkits.mplot3d import Axes3D
    
    Then be sure the axis is created using:
    
        ax = fig.add_subplot(111, projection='3d')
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    #from mayavi import mlab
    
    name = pdf['name']
    xattr = pdf['xattr']
    xmin = pdf['xmin']
    dx = pdf['dx']
    nx = pdf['nx']
    yattr = pdf['yattr']
    ymin = pdf['ymin']
    dy = pdf['dy']
    ny = pdf['ny']
    data = pdf['data']
    
    idx = np.nonzero(data<0)
    data[idx] = np.nan
    
    #  Calculate max X and Y axis bin centres
    xmax = xmin + dx*nx
    ymax = ymin + dy*ny
    
    #  Build 1D arrays with bin centre coordinates
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    
    #  Build 2D arrays of bin centre coordinates for plotting
    X, Y = np.meshgrid(x, y)
    
    cont_max = np.nanmax(data)
    if cont_inc == -1:
        cont_inc = cont_max / 10.0
    cont_min = cont_inc
    
    print('Contour Min: %f' % cont_min)
    print('Contour Max: %f' % cont_max)
    print('Contour Inc: %f' % cont_inc)
    
    levels = np.arange(cont_min, cont_max+cont_inc, cont_inc)
    cs = ax.contour(X, Y, data, levels=levels, colors=color, linewidth=2.0)
    cs.collections[0].set_label(name)
    
    ax.set_xlabel(xattr)
    ax.set_ylabel(yattr)