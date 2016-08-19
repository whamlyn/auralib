"""
AuraQI module containing helper functions for reading data from RokDoc exports.

Author:   Wes Hamlyn
Created:  25-Mar-2016
Last Mod: 17-Aug-2016
"""



def load_2d_pdf(infile):
    import numpy as np
    
    fd1 = open(infile, 'r')
    
    for i in range(0, 11):
        buf = fd1.readline()
        
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


def plot_2d_pdf(ax, pdf, linecolor='b'):
    """
    Need to use mplot3d axes.  Make sure to import this using the following
    syntax:
    
        from mpl_toolkits.mplot3d import Axes3D
    
    Then be sure the axis is created using:
    
        ax = fig.add_subplot(111, projection='3d')
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    #from mayavi import mlab
    
    xattr = pdf['xattr']
    xmin = pdf['xmin']
    dx = pdf['dx']
    nx = pdf['nx']
    yattr = pdf['yattr']
    ymin = pdf['ymin']
    dy = pdf['dy']
    ny = pdf['ny']
    data = pdf['data']
    
    #  Calculate max X and Y axis bin centres
    xmax = xmin + dx*nx
    ymax = ymin + dy*ny
    
    #  Build 1D arrays with bin centre coordinates
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    
    #  Build 2D arrays of bin centre coordinates for plotting
    X, Y = np.meshgrid(x, y)

    im = ax.plot_surface(X, Y, data, cmap=plt.cm.rainbow, 
                     shade=True, rstride=1, cstride=1)    