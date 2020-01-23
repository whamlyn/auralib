"""
auralib module containing helper functions for reading data from RokDoc exports.

Author:   Wes Hamlyn
Created:  25-Mar-2016
Last Mod: 17-Aug-2016

"""

import numpy as np
from scipy.interpolate import interp1d



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


def load_1d_pdf(infile):
        
    with open(infile, 'r') as fd:
        buf = fd.readlines()
    
    mid_bin_min = float(buf[3].split(':')[1].strip())
    bin_width = float(buf[4].split(':')[1].strip())
    num_bins =  int(buf[5].split(':')[1].strip())
    
    bin_centers = np.arange(num_bins)*bin_width + mid_bin_min
    bin_edges = np.arange(num_bins+1)*bin_width + mid_bin_min - bin_width/2
    
    counts = np.array(buf[8].strip().split(), dtype=np.float32)
    norm_counts = counts / np.max(counts)
    
    f = interp1d(bin_centers, counts, kind='cubic', bounds_error=None, fill_value='extrapolate')
    bin_centers_fine = np.linspace(bin_centers[0], bin_centers[-1], num_bins*10)
    counts_fine = f(bin_centers_fine)
    counts_norm_fine = counts_fine/np.max(counts_fine)
    
    pdf = {'mid_bin_min': mid_bin_min,
           'bin_width': bin_width,
           'num_bins': num_bins,
           'bin_centers': bin_centers,
           'bin_edges': bin_edges,
           'counts': counts, 
           'norm_counts': norm_counts,
           'bin_centers_fine': bin_centers_fine,
           'counts_fine': counts_fine,
           'norm_counts_fine': counts_norm_fine}
    
    return pdf

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


def dump_sgy_from_sxy(sxyfile):
    """
    Dump the segy file path which is associated with a RokDoc index (sxy) file.
    
    This is a first-iteration and may be buggy. No actual knowledge of the sxy
    format has been used in this function, only parsing based on observed 
    patters in several sxy files.
    
    Inputs:
        sxyfile = (string) Path to RokDoc sxy index file

    Outputs:
        sgyfile = (string) Path to sgy file associated with sxy index file
    
    Written by: Wes Hamlyn
    Created:    20-Sep-2017
    Last Mod:   20-Sep-2017
    """
    
    with open(sxyfile, 'rb') as fd:
        fd.seek(0)
        buf = fd.read(1024)

    nbytes = len(buf)

    buf1 = buf.split(b':')
    buf2 = buf1[2].split(b'sgy')[0]

    drive_name = buf1[1].decode(errors='ignore')[-1]
    file_path = buf2.decode(errors='ignore')
    file_ext = 'sgy'

    sgyfile = '' + drive_name + ':' + file_path + file_ext

    return sgyfile



def load_rokdoc_well_markers(infile):
    """
    Function to load well markers exported from RokDoc in ASCII format.
    """
    
    with open(infile, 'r') as fd:
        buf = fd.readlines()
    
    
    marker = []
    well = []
    md = []
    tvdkb = []
    twt = []
    tvdss = []
    x = []
    y = []
    
    for line in buf[5:]:
    
        c1, c2, c3, c4, c5 = line.split("'")
        c6, c7, c8, c9, c10, c11 = c5.strip().split()
        
        marker.append(c2)
        well.append(c4)
        md.append(float(c6))
        tvdkb.append(float(c7))
        twt.append(float(c8))
        tvdss.append(float(c9))
        x.append(float(c10))
        y.append(float(c11))
        
    
    markers = {}
    for each in list(set(well)):
        markers[each] = {}
        
    for i in range(len(marker)):
        cur_well = well[i]
        cur_marker = marker[i]
        cur_md = md[i]
        cur_tvdkb = tvdkb[i]
        cur_tvdss = tvdss[i]
        cur_twt = twt[i]
        cur_x = x[i]
        cur_y = y[i]
        
        markers[cur_well][cur_marker] = {'md': cur_md, 'tvdkb': cur_tvdkb,
                                         'tvdss': cur_tvdss, 'twt': cur_twt,
                                         'x': cur_x, 'y': cur_y}
    
    return markers


def read_rokdoc_fluidset(infile, other1='Other 1', other2='Other 2', other3='Other 3'):
    """
    Convenience function to load a fluid set exported from RokDoc to a Python
    dictionary.
    """
    
    # import required library for parsing xml documents
    import xml.etree.ElementTree as ET
    
    # open the xml fluid set file
    tree = ET.parse(infile)
    root = tree.getroot()
    
    # define fluid type names corresponding to the fluid codes found in the
    # RokDoc fluid set file
    fluidTypeLookup = {'0': 'Water',
                       '1': 'Oil',
                       '2': 'Gas',
                       '3': 'Condensate',
                       '5': 'CO2',
                       '6': 'Other 1',
                       '7': 'Other 2',
                       '8': 'Other 3',
                       '9': 'Heavy Oil'}
    
    
    # Optionally rename the "other" fluid types to have specific names, for
    # example "steam" instead of simply "Other 1"
    fluidTypeLookup['6'] = other1
    fluidTypeLookup['7'] = other2
    fluidTypeLookup['8'] = other3
    
    # initialize an empty dictionary to store the fluid set data
    fset = {}
    
    # read through the fluid set description and parse out relevant information
    # this section of code can likely be expanded as not all perturbations of
    # fluid set files has been investigated. This should work reasonably well
    # for FLAG fluid calculations, no idea how it will work for Batzle-Wang
    # calculated fluid sets.
    desc = root[0].text.split('\n')
    fset['description'] = desc  # store the description as-is for reference
    
    # now do description parsing...
    for line in desc:
        
        if line.split(':')[0].strip() == 'Target depth pressure':
            elem = line.split(':')[1].strip().split()
            fset['Press'] = float(elem[0])
            fset['PressUnit'] = elem[1]
        
        if line.split(':')[0].strip() == 'Target depth temperature':
            elem = line.split(':')[1].strip().split()
            fset['Temp'] = float(elem[0])
            fset['TempUnit'] = elem[1]
        
        if line.split(':')[0].strip() == 'Salinity':
            elem = line.split(':')[1].strip().split()
            fset['Salinity'] = float(elem[0])
            fset['SalinityUnit'] = elem[1]        
        
        if line.split(':')[0].strip() == 'Gas  gravity':
            elem = line.split(':')[1].strip().split()
            fset['GasGravity'] = float(elem[0])
            fset['GasGravityUnit'] = elem[1] 
            
        if line.split(':')[0].strip() == 'Heavy Oil API':
            elem = line.split(':')[1].strip().split()
            fset['HeavyOilAPI'] = float(elem[0])
            fset['HeavyOilAPIUnit'] = elem[1]
            
        if line.split(':')[0].strip() == 'Frequency':
            elem = line.split(':')[1].strip().split()
            fset['Frequency'] = float(elem[0])
            fset['FrequencyUnit'] = elem[1]
            
    # now read and parse the fluid elastic properties
    for elem in root[1]:
        k = float(elem.attrib['api'].replace(',', ''))
        mu = float(elem.attrib['mu'].replace(',', ''))
        vp = float(elem.attrib['vp'].replace(',', ''))
        vs = float(elem.attrib['vs'].replace(',', ''))
        rho = float(elem.attrib['rho'].replace(',', ''))
        
        fluidTypeIndex = elem.attrib['fluidTypeIndex']
        fluidType = fluidTypeLookup[fluidTypeIndex]
        
        fset[fluidType] = {'K': k, 'G': mu, 'R': rho, 'Vp': vp, 'Vs': vs}
    
    # Done!
    
    return fset