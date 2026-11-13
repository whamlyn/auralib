"""
auralib module for various utility functions.

Author:   Wes Hamlyn
Created:  16-Aug-2016
Last Mod: 17-Aug-2016

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.interpolate import interp1d, interp2d, pchip
import win32clipboard as cb
from matplotlib.patches import Ellipse, Polygon
import matplotlib.transforms as transforms
from matplotlib.widgets import PolygonSelector
from matplotlib.backend_bases import LocationEvent
import os
from osgeo import gdal
from pyproj import CRS, Transformer
import math
from pathlib import Path
from typing import Any
import rasterio


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

    pts = plt.ginput(n=-1, timeout=60, show_clicks=False,
                     mouse_add=1, mouse_pop=2, mouse_stop=3)
    
    pts = np.array(pts)
    xpoly = pts[:, 0]
    ypoly = pts[:, 1]

    return xpoly, ypoly


def plot_poly(x, y, ax, **kwargs):
    """
    Convenience funciton to plot a polygon digitized using the above
    digitize_poly() function.
    """
    x = np.array(x)
    y = np.array(y)
    
    xy = np.vstack([x, y]).T
    
    poly1 = plt.Polygon(xy, closed=True, alpha=0.3, **kwargs)
    ax.add_patch(poly1)


def inpoly_old(datax, datay, polyx, polyy):
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
    
    poly = Polygon(np.vstack([polyx, polyy]).T, closed=True)
    mask = poly.contains_points(np.vstack([datax, datay]).T)
    #idx = np.nonzero(mask)

    return mask

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
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
                      facecolor=facecolor, **kwargs)

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



def get_dist(ax, show_azim=False):
    """
    Convenience function to get a distance from a map and plot line in axes
    """

    pts = plt.ginput(n=2)

    x1 = pts[0][0]
    x2 = pts[1][0]
    y1 = pts[0][1]
    y2 = pts[1][1]
    
    dx = x2-x1
    dy = y2-y1
    dist = np.sqrt( dx**2 + dy**2 )

    azim = np.arctan2(dx, dy) * 180/np.pi
    if azim<0:
        azim = azim + 360
    
    if show_azim:
        dist_txt = '%.1f m\n%.1f deg' % (dist, azim)
    else:
        dist_txt = '%.1f m\n%.1f deg' % (dist, azim)

    ax.plot([x1, x2], [y1, y2], 'k-')
    ax.plot([x1, x2], [y1, y2], 'ko')

    xtxt = (x1+x2)/2
    ytxt = (y1+y2)/2
    
    bbox_params = {'boxstyle': 'Round', 'facecolor': 'white',
                   'edgecolor': 'none', 'alpha': 0.7}
    
    ax.text(xtxt, ytxt, dist_txt, bbox=bbox_params, ha='center', va='center', rotation=(90-azim))

    #plt.draw()

    return x1, x2, y1, y2, dist, azim



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
    logsm = np.convolve(logpad, filt, 'valid')

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


def block_log(depth, values, block_boundaries):
    """
    Block a single well log using user-defined depth boundaries.

    Parameters
    ----------
    depth : array_like, shape (m,) or (1, m)
        Depth reference for the original log samples.

    values : array_like, shape (m,) or (1, m)
        Original log values. NaN values are ignored when calculating
        block averages.

    block_boundaries : array_like, shape (n,) or (1, n)
        Depths defining the block boundaries. The boundaries must be
        strictly increasing. There will be n - 1 output blocks.

    Returns
    -------
    blocked_depth : ndarray, shape (n - 1,)
        Midpoint depth of each block.

    blocked_values : ndarray, shape (n - 1,)
        Mean log value within each block. A block containing no valid
        log samples is assigned NaN.

    Notes
    -----
    Each block is defined as

        block_boundaries[i] <= depth < block_boundaries[i + 1]

    except for the final block, which includes its lower and upper
    boundaries.
    """
    depth = np.asarray(depth, dtype=float).reshape(-1)
    values = np.asarray(values, dtype=float).reshape(-1)
    boundaries = np.asarray(
        block_boundaries,
        dtype=float,
    ).reshape(-1)

    if depth.size != values.size:
        raise ValueError(
            "depth and values must contain the same number of samples."
        )

    if boundaries.size < 2:
        raise ValueError(
            "At least two block boundaries are required."
        )

    if not np.all(np.isfinite(depth)):
        raise ValueError("depth cannot contain NaN or infinite values.")

    if not np.all(np.isfinite(boundaries)):
        raise ValueError(
            "block_boundaries cannot contain NaN or infinite values."
        )

    if np.any(np.diff(depth) <= 0):
        raise ValueError("depth must be strictly increasing.")

    if np.any(np.diff(boundaries) <= 0):
        raise ValueError(
            "block_boundaries must be strictly increasing."
        )

    number_of_blocks = boundaries.size - 1

    # Use the midpoint as the representative depth of each block.
    #blocked_depth = 0.5 * ( boundaries[:-1] + boundaries[1:] )
    
    # Use the top of the interval as the representative depth of each block.
    blocked_depth = boundaries[:-1]

    blocked_values = np.full(
        number_of_blocks,
        np.nan,
        dtype=float,
    )

    for i in range(number_of_blocks):
        top = boundaries[i]
        base = boundaries[i + 1]

        # Make every block left-inclusive and right-exclusive.
        # Include the lower boundary of the final block as usual and
        # also include its upper boundary.
        if i < number_of_blocks - 1:
            in_block = (depth >= top) & (depth < base)
        else:
            in_block = (depth >= top) & (depth <= base)

        valid = in_block & np.isfinite(values)

        if np.any(valid):
            blocked_values[i] = np.mean(values[valid])

    return blocked_depth, blocked_values


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


def calc_rms(x, **kwargs):
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
    rms = np.sqrt(np.sum(x**2, **kwargs)/n)
    return rms


def calc_nrms(x1, x2, *args, **kwargs):
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

    rms1 = calc_rms(x1, *args, **kwargs)
    rms2 = calc_rms(x2, *args, **kwargs)
    rms12 = calc_rms(x1-x2, *args, **kwargs)

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


def apply_agc_2d(traces, sample_interval, window_length, normalize=False,
                 epsilon=1e-12):
    """
    Apply RMS automatic gain control independently to multiple seismic traces.

    Parameters
    ----------
    traces : array-like
        Two-dimensional array with shape (m, n), where:
            m = number of traces
            n = number of samples per trace

    sample_interval : float
        Trace sample interval in seconds.

    window_length : float
        AGC window length in seconds.

    epsilon : float, optional
        Minimum RMS amplitude required to apply gain. This prevents
        division by zero and excessive gain in zero-energy intervals.

    Returns
    -------
    agc_traces : numpy.ndarray
        AGC-scaled traces with shape (m, n).

    gain : numpy.ndarray
        Time-varying AGC gain functions with shape (m, n).

    local_rms : numpy.ndarray
        Local RMS amplitudes with shape (m, n).
    """
    traces = np.asarray(traces, dtype=float)

    if traces.ndim != 2:
        raise ValueError(
            "traces must have shape (number_of_traces, number_of_samples)."
        )

    if traces.size == 0:
        raise ValueError("traces cannot be empty.")

    if sample_interval <= 0:
        raise ValueError("sample_interval must be greater than zero.")

    if window_length <= 0:
        raise ValueError("window_length must be greater than zero.")

    number_of_traces, number_of_samples = traces.shape

    # Convert the window length from seconds to samples.
    window_samples = int(round(window_length / sample_interval))
    window_samples = max(window_samples, 1)

    # Use an odd window so it can be centred on each sample.
    if window_samples % 2 == 0:
        window_samples += 1

    # The window does not need to be longer than the trace.
    if window_samples > number_of_samples:
        window_samples = number_of_samples

        # Keep the window odd where possible.
        if window_samples % 2 == 0 and window_samples > 1:
            window_samples -= 1

    boxcar = np.ones(window_samples, dtype=float)

    local_rms = np.empty_like(traces)

    # Number of valid samples in each window. This corrects the
    # normalization near the beginning and end of the traces.
    valid_sample_count = np.convolve(
        np.ones(number_of_samples),
        boxcar,
        mode="same",
    )

    # Calculate the local RMS independently for each trace.
    for trace_index in range(number_of_traces):
        local_sum_square = np.convolve(
            traces[trace_index] ** 2,
            boxcar,
            mode="same",
        )

        local_mean_square = local_sum_square / valid_sample_count
        local_rms[trace_index] = np.sqrt(local_mean_square)

    # Calculate the time-varying gain.
    gain = np.zeros_like(traces)

    valid = local_rms > epsilon
    gain[valid] = 1.0 / local_rms[valid]

    # Apply AGC.
    agc_traces = traces * gain
    
    if normalize:
        nsamp_total = traces.size
        rms_in = np.sqrt(np.sum(traces**2)/(nsamp_total))
        rms_out = np.sqrt(np.sum(agc_traces**2)/(nsamp_total))
        scalar = rms_in/rms_out
        agc_traces = agc_traces*scalar
    
    return agc_traces #, gain, local_rms


# Some functions related to well geometry and domain transforms. These have
# been partially implemented in auralib.well.DevSvy but that really needs some
# more thought to be a properly designed class.
def dev_to_xyz(inc, azm, md, md_step=1.0, method='pchip'):
    """
    Convert deviation survey measured in inclination, azimuth, measured depth
    to x, y, TVDkb using spline interpolant.
    """

    # initialize 1D interpolation functions for inclinations and azimuths
    if method=='pchip':
        finc = pchip(md, inc, extrapolate=True)
        fazm = pchip(md, azm, extrapolate=True)

    else:
        finc = interp1d(md, inc, kind=method, bounds_error=False,
                        fill_value='extrapolate')
        fazm = interp1d(md, azm, kind=method, bounds_error=False,
                        fill_value='extrapolate')

    # interpolate inclinations and azimuths to regularly sampled MD
    mdi = np.arange(md[0], md[-1]+md_step, md_step)
    inci = finc(mdi)
    azmi = fazm(mdi)

    # convert from inc, azm, md to TVDkb, xoffset, and yoffset from tophole
    dmdi = np.diff(mdi)
    dmdi = np.hstack([dmdi, dmdi[-1]])
    dzi = np.cos(np.deg2rad(inci))*dmdi
    dsi = np.sin(np.deg2rad(inci))*dmdi
    dxi = np.sin(np.deg2rad(azmi))*dsi
    dyi = np.cos(np.deg2rad(azmi))*dsi
    tvdkb = np.cumsum(dzi)
    dx = np.cumsum(dxi)
    dy = np.cumsum(dyi)

    return dx, dy, tvdkb, mdi


def md_to_xyz(md, dev, dev_type='angular'):
    """
    Given a deviation survey, calculate the TVDkb from an MD value.

    dev = 2D numpy array with the following convention:
             dev[0, :] = inclination or x-offset from tophole
             dev[1, :] = azimuth  or y-offset from tophole
             dev[2, :] = measured depth along wellbore measured from KB

    Returns either int/float or numpy arrays for each coordinate
    """

    if dev_type=='angular':
        xoff, yoff, tvdkb, mdi = dev_to_xyz(dev[0, :], dev[1, :], dev[2, :])

    elif dev_type=='cartesian':
        xoff = dev[:, 0]
        yoff = dev[:, 1]
        tvdkb = dev[:, 2]

    fx = pchip(mdi, xoff, extrapolate=True)
    fy = pchip(mdi, yoff, extrapolate=True)
    ftvdkb = pchip(mdi, tvdkb, extrapolate=True)

    x = fx(md)
    y = fy(md)
    tvdkb = ftvdkb(md)

    if type(md) in (int, float):
        x = x.tolist()
        y = y.tolist()
        tvdkb = tvdkb.tolist()

    return x, y, tvdkb


class Bayes2D():
    """
    Class for performing Bayesian classification using two attributes.
    """

    def __init__(self):
        pass

    def make_pdf_nonpar(self, datax, datay, nbins=25, smoothing=0.5):
        """
        Method to build a 2D non-parametric PDF from two vectors of data
        """

        # Build a 2D histogram from two 1D data vectors
        hist1, binsx, binsy = np.histogram2d(datax, datay, bins=[nbins, nbins])
        hist1 = hist1 / hist1.sum()
        hist1 = hist1.T # normalize by area under histogram

        # Smooth histogram using 2D gaussian filter
        hist2 = sp.ndimage.filters.gaussian_filter(hist1, smoothing,
                                                   mode='constant')
        hist2 = hist2/hist2.sum() # normalize by area under histogram

        # Calculate bin centers from bin edges
        dx = np.diff(binsx)
        dy = np.diff(binsy)
        binsxc = binsx[:-1]+dx/2
        binsyc = binsy[:-1]+dy/2

        f = interp2d(binsxc, binsyc, hist2.T, kind='cubic', bounds_error=False,
                     fill_value=0.0)
        X = np.linspace(-4, 4, nbins*4)
        Y = np.linspace(-4, 4, nbins*4)

        pdf2 = []
        for i in range(len(X)):
            for j in range(len(Y)):
                pdf2.append(f(X[i], Y[j]))
        pdf2 = np.array(pdf2).reshape([len(X), len(Y)])

        self.hist = hist1
        self.kernel = hist2
        self.datax = datax
        self.datay = datay
        self.binsx = binsx
        self.binsy = binsy
        self.binsxc = binsxc
        self.binsyc = binsyc
        self.pdf = pdf2
        self.pdf_func = f

    def calc_prob(self, x, y, prior=1.0):
        x = np.array(x)
        y = np.array(y)

        ns = len(x)
        prob = np.zeros(ns)
        for i in range(ns):
            prob[i] = self.pdf_func(x[i], y[i]) * prior

        return prob


def get_closest_point(px, py, x, y, tol=-1.0):
    """
    Function to return the index of closest sample in coordinate vectors (x, y)
    to point (px, py). 
    """
    
    r = np.sqrt((px-x)**2 + (py-y)**2)
    idx = np.argmin(r)
    
    if (tol >= 0.0) & (r[idx] > tol):
        idx = np.nan
    
    return idx


def generate_tfw(tif_path, world_ext=".tfw"):
    """
    Create a world file from the georeferencing in a GeoTIFF.

    Parameters
    ----------
    tif_path : str or Path
        Path to the input GeoTIFF.
    world_ext : str, optional
        Extension for the world file. Default is ".tfw".

    Returns
    -------
    Path
        Path to the created world file.
    """
    tif_path = Path(tif_path)

    with rasterio.open(tif_path) as src:
        transform = src.transform

    # Rasterio affine transform:
    #
    # x = a * col + b * row + c
    # y = d * col + e * row + f
    #
    # World file:
    #   Line 1: A  x pixel size
    #   Line 2: D  row rotation
    #   Line 3: B  column rotation
    #   Line 4: E  y pixel size (usually negative)
    #   Line 5: C  x coordinate of CENTER of upper-left pixel
    #   Line 6: F  y coordinate of CENTER of upper-left pixel

    A = transform.a
    D = transform.d
    B = transform.b
    E = transform.e

    # Rasterio's c/f refer to the upper-left CORNER.
    # World files require the CENTER of the upper-left pixel.
    C = transform.c + 0.5 * (transform.a + transform.b)
    F = transform.f + 0.5 * (transform.d + transform.e)

    world_path = tif_path.with_suffix(world_ext)

    with open(world_path, "w") as f:
        f.write(f"{A:.15g}\n")
        f.write(f"{D:.15g}\n")
        f.write(f"{B:.15g}\n")
        f.write(f"{E:.15g}\n")
        f.write(f"{C:.15g}\n")
        f.write(f"{F:.15g}\n")

    return world_path


def generate_tfw2(geotiff_file, gen_prj=False):
    """
    Python function to ouput world file and projection file from a geotiff. The
    world and projection files will share the same basefilename as the geotiff
    file but will be assigned extensions .tfw or .prj respectively. Modified
    from an original script from MerseyViking found at:
    
    https://gis.stackexchange.com/questions/9421/creating-tfw-and-prj-files-for-folder-of-geotiff-files
    
    Note from original author:
        Written by MerseyViking (mersey dot viking at gmail dot com), 2011.
        Released into the public domain - May 8, 2011
        I accept no responsibility for any errors or loss of data, revenue, or
        life this script may cause. Use at your own risk.

    Parameters
    ----------
    geotiff_file : string
        Path to geotiff file the world file will be created from.
    gen_prj : boolean, optional
        If True, also output the ESRI format projection file and output in 
        wellknown text format. The default is False.

    Returns
    -------
    None.

    """

    import osgeo.gdal as gdal
    import osgeo.osr as osr
    import os

    src = gdal.Open(geotiff_file)
    xform = src.GetGeoTransform()
    
    if gen_prj:
        src_srs = osr.SpatialReference()
        src_srs.ImportFromWkt(src.GetProjection())
        src_srs.MorphToESRI()
        src_wkt = src_srs.ExportToWkt()
        
        prj_file = os.path.splitext(geotiff_file)[0] + '.prj'
        with open(prj_file, 'w') as prj:
            prj.write(src_wkt)

    edit1=xform[0]+xform[1]/2
    edit2=xform[3]+xform[5]/2
    
    tfw_file = os.path.splitext(geotiff_file)[0] + '.tfw'
    with open(tfw_file, 'w') as tfw:
        tfw.write("%0.8f\n" % xform[1])
        tfw.write("%0.8f\n" % xform[2])
        tfw.write("%0.8f\n" % xform[4])
        tfw.write("%0.8f\n" % xform[5])
        tfw.write("%0.8f\n" % edit1)
        tfw.write("%0.8f\n" % edit2)
        

def generate_tfw_old(infile, outfile=None):
    """
    Create a tfw World file from a geotiff with coordinate information stored
    in the tiff headers. Requires the osgeo package.
    
    Note: There are limitations with this current implementation. If the 
          reference corner of the image is anything other than the upper left
          I think this will fail. Other scenarios need to be validated.
          Wes Hamlyn, 2025-03-19

    Parameters
    ----------
    infile : String
        Path to the geotiff file that the tfw will be created for.
    outfile : String, optional
        Path to the output tfw file. If no outfile path is specified, the tfw
        file will be written to the same directory as the input geotiff file
        with the same base filename and extension "tfw. The default is None.

    Returns
    -------
    None.

    """
        
    if outfile == None:
        outfile = os.path.splitext(infile)[0] + '.tfw'
        
    dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    
    gt = dataset.GetGeoTransform()
    
    nr = dataset.RasterYSize
    nc = dataset.RasterXSize
    
    print(gt)
    
    xref = gt[0]
    dxc = gt[1]
    dxr = gt[2]
    yref = gt[3]
    dyc = gt[4]
    dyr = gt[5]
    
    print('')
    print(xref)
    print(dxc)
    print(dxr)
    print(yref)
    print(dyc)
    print(dyr)
    print('')    
    
    if dxc < 0: # origin is on right
        if dyr < 0: # origin is on top
            origin = 'ur'
        else: # origin is on bottom
            origin = 'lr'
    else:  # origin is on left
        if dyr < 0: # origin is on top
            origin = 'ul'
        else: # origin is on bottom
            origin = 'll'
    
    # calculate the x and y of the upper left corner
    
    if origin == 'ur':
        #xref_ul = xref + nc*dxc + 0*dxr
        #yref_ul = yref + nc*dyc + 0*dyr
        xref_ul = xref
        yref_ul = yref
        dxr = dxr*-1
        dyc = dyc*-1
        
    elif origin == 'ul':
        xref_ul = xref
        yref_ul = yref
        
    elif origin == 'lr':
        pass
    
    elif origin == 'll':
        pass    
    
    xref0 = xref_ul + dxc/2 + dxr/2
    yref0 = yref_ul + dyc/2 + dyr/2
    
    print('')
    print(xref)
    print(dxc)
    print(dxr)
    print(yref)
    print(dyc)
    print(dyr)
    print('') 
    
    #print(xref0)
    #print(yref0)
    
    #x = xref0 + dxc*nc + dxr*nr
    #y = yref0 + dyc*nc + dyr*nr
    
    #print(x, y)
        
    # TFW output
    with open(outfile, 'w') as fd:
        fd.write('%f\n' % dxc)
        fd.write('%f\n' % dxr)
        fd.write('%f\n' % dyc)
        fd.write('%f\n' % dyr)
        fd.write('%f\n' % xref0)
        fd.write('%f\n' % yref0)
        
        
        
def calculate_grid_convergence(lat_deg, lon_deg, utm_zone, datum='NAD27',
                               northern_hemisphere=True):
    """
    Calculate the grid convergence angle for a UTM zone given the lat and long
    """
    
    # Define geographic CRS (NAD27)
    nad27 = CRS.from_epsg(4267)

    # Define UTM CRS based on NAD27
    utm_crs = CRS.from_proj4(f"+proj=utm +zone={utm_zone} +datum={datum} +units=m +{'north' if northern_hemisphere else 'south'}")

    # Create transformer from NAD27 geographic to NAD27 UTM
    transformer = Transformer.from_crs(nad27, utm_crs, always_xy=True)

    # Transform to UTM coordinates
    easting, northing = transformer.transform(lon_deg, lat_deg)

    # Get central meridian of the UTM zone
    central_meridian = (utm_zone - 1) * 6 - 180 + 3

    # Convert degrees to radians
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)
    cm_rad = math.radians(central_meridian)

    # Compute grid convergence
    gamma_rad = math.atan(math.tan(lon_rad - cm_rad) * math.sin(lat_rad))
    gamma_deg = math.degrees(gamma_rad)

    return gamma_deg


def grid_convergence_from_utm(easting, northing, utm_zone, datum='NAD27',
                              northern_hemisphere=True):
    """
    Calculate the grid convergence angle for a x/y location in a utm zone
    """    
    
    # Define UTM CRS based on NAD27
    utm_crs = CRS.from_proj4(f"+proj=utm +zone={utm_zone} +datum={datum} +units=m +{'north' if northern_hemisphere else 'south'}")

    # Define geographic CRS (NAD27)
    nad27 = CRS.from_epsg(4267)

    # Create transformer from UTM to geographic
    transformer = Transformer.from_crs(utm_crs, nad27, always_xy=True)

    # Transform UTM to geographic coordinates
    lon_deg, lat_deg = transformer.transform(easting, northing)

    # Get central meridian of the UTM zone
    central_meridian = (utm_zone - 1) * 6 - 180 + 3

    # Convert degrees to radians
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)
    cm_rad = math.radians(central_meridian)

    # Compute grid convergence
    gamma_rad = math.atan(math.tan(lon_rad - cm_rad) * math.sin(lat_rad))
    gamma_deg = math.degrees(gamma_rad)

    return gamma_deg


def nmo_curve(t0, v_nmo, x, eta=None):
    """
    Convenience function to calculate the normal moveout times for a 2nd or 4th
    order normal moveout curve
    """
    
    term2 = (x / v_nmo)**2
    t_nmo = (t0**2 + term2)**0.5
    
    if eta is not None:
        print('calculating 4th order term')
        term4 = -2*eta/t0**2 * (x / v_nmo)**4
        t_nmo = (t0**2 + term2 + term4)**0.5        
        #term4 = -2*eta*x**4 / ( v_nmo**2 * (t0**2*v_nmo**2 + (1+2*eta)*x**2) )
        #t_nmo = (t0**2 + term2 + term4)**0.5
    
    return t_nmo


def dix_vrms_to_vint(time, vrms):
    """
    Convert RMS velocities to interval velocities using the Dix equation.

    Parameters
    ----------
    time : array_like
        Two-way travel times. Must be monotonically increasing.
        Units can be seconds or milliseconds.
        
    vrms : array_like
        RMS velocities corresponding to each time sample.
        Typically in m/s.

    Returns
    -------
    vint : ndarray
        Dix interval velocities, with the same length as the input arrays.

    Notes
    -----
    The first interval is assumed to extend from t=0 to time[0],
    so its interval velocity is equal to vrms[0].
    """

    time = np.asarray(time, dtype=float)
    vrms = np.asarray(vrms, dtype=float)

    if time.ndim != 1 or vrms.ndim != 1:
        raise ValueError("time and vrms must be 1D arrays.")

    if len(time) != len(vrms):
        raise ValueError("time and vrms must have the same length.")

    if len(time) == 0:
        return np.array([])

    if np.any(np.diff(time) <= 0):
        raise ValueError("time must be monotonically increasing.")

    vint = np.empty_like(vrms)

    # First interval: from zero time to time[0]
    vint[0] = vrms[0]

    # Dix equation
    numerator = (
        vrms[1:]**2 * time[1:]
        - vrms[:-1]**2 * time[:-1]
    )

    denominator = time[1:] - time[:-1]

    vint[1:] = np.sqrt(numerator / denominator)

    return vint





## Below functions are related to SEGP1 file reading
from pathlib import Path
from typing import Any

import numpy as np


def _int_or_none(text: str) -> int | None:
    """
    Convert a fixed-width field to int.

    Returns None if the field is blank or cannot be converted.
    """
    text = text.strip()

    if not text:
        return None

    try:
        return int(text)
    except ValueError:
        return None


def _float_or_none(text: str) -> float | None:
    """
    Convert a fixed-width field to float.

    Returns None if the field is blank or cannot be converted.
    """
    text = text.strip()

    if not text:
        return None

    try:
        return float(text)
    except ValueError:
        return None


def _parse_latitude(text: str) -> dict[str, Any] | None:
    """
    Parse a standard SEG P1 latitude field.

    Expected format:

        DDMMSSssH

    where:

        DD   = degrees
        MM   = minutes
        SSss = seconds to two decimal places
        H    = N or S

    The decimal point in the seconds field is normally implicit.

    Example:

        17543354N

    represents:

        17 degrees
        54 minutes
        33.54 seconds North
    """
    raw = text
    text = text.strip()

    if not text:
        return None

    hemisphere = text[-1].upper()

    if hemisphere not in ("N", "S"):
        return {
            "raw": raw,
            "degrees": None,
            "minutes": None,
            "seconds": None,
            "hemisphere": None,
            "decimal_degrees": None,
        }

    numeric = text[:-1]

    try:
        degrees = int(numeric[0:2])
        minutes = int(numeric[2:4])

        seconds_string = numeric[4:]

        if "." in seconds_string:
            seconds = float(seconds_string)
        else:
            seconds = int(seconds_string) / 100.0

    except (ValueError, IndexError):
        return {
            "raw": raw,
            "degrees": None,
            "minutes": None,
            "seconds": None,
            "hemisphere": hemisphere,
            "decimal_degrees": None,
        }

    decimal_degrees = (
        degrees
        + minutes / 60.0
        + seconds / 3600.0
    )

    if hemisphere == "S":
        decimal_degrees *= -1.0

    return {
        "raw": raw,
        "degrees": degrees,
        "minutes": minutes,
        "seconds": seconds,
        "hemisphere": hemisphere,
        "decimal_degrees": decimal_degrees,
    }


def _parse_longitude(text: str) -> dict[str, Any] | None:
    """
    Parse a standard SEG P1 longitude field.

    Expected format:

        DDDMMSSssH

    where:

        DDD  = degrees
        MM   = minutes
        SSss = seconds to two decimal places
        H    = E or W

    The decimal point in the seconds field is normally implicit.

    Example:

        110445881E

    represents:

        110 degrees
        44 minutes
        58.81 seconds East
    """
    raw = text
    text = text.strip()

    if not text:
        return None

    hemisphere = text[-1].upper()

    if hemisphere not in ("E", "W"):
        return {
            "raw": raw,
            "degrees": None,
            "minutes": None,
            "seconds": None,
            "hemisphere": None,
            "decimal_degrees": None,
        }

    numeric = text[:-1]

    try:
        degrees = int(numeric[0:3])
        minutes = int(numeric[3:5])

        seconds_string = numeric[5:]

        if "." in seconds_string:
            seconds = float(seconds_string)
        else:
            seconds = int(seconds_string) / 100.0

    except (ValueError, IndexError):
        return {
            "raw": raw,
            "degrees": None,
            "minutes": None,
            "seconds": None,
            "hemisphere": hemisphere,
            "decimal_degrees": None,
        }

    decimal_degrees = (
        degrees
        + minutes / 60.0
        + seconds / 3600.0
    )

    if hemisphere == "W":
        decimal_degrees *= -1.0

    return {
        "raw": raw,
        "degrees": degrees,
        "minutes": minutes,
        "seconds": seconds,
        "hemisphere": hemisphere,
        "decimal_degrees": decimal_degrees,
    }


def _parse_data_record(
    line: str,
    line_number: int | None = None,
) -> dict[str, Any]:
    """
    Parse one SEG P1 80-character data record.

    SEG P1 column layout
    --------------------
    1       Data record identifier (blank)
    2-17    Line name
    18-25   Shotpoint/event parameter
    26      Reshoot code
    27-35   Latitude
    36-45   Longitude
    46-53   Map grid Easting
    54-61   Map grid Northing
    62-66   Water depth/elevation
    67-68   Year
    69-71   Day of year
    72-73   Hour GMT
    74-75   Minute GMT
    76-77   Second GMT
    78-80   Spare
    """

    record = line.ljust(80)[:80]

    return {
        "line_name": record[1:17].rstrip(),
        "shotpoint": _int_or_none(record[17:25]),
        "reshoot_code": record[25:26].strip() or None,

        "latitude": _parse_latitude(record[26:35]),
        "longitude": _parse_longitude(record[35:45]),

        "easting": _float_or_none(record[45:53]),
        "northing": _float_or_none(record[53:61]),
        "elevation": _float_or_none(record[61:66]),

        "year": _int_or_none(record[66:68]),
        "day_of_year": _int_or_none(record[68:71]),
        "hour": _int_or_none(record[71:73]),
        "minute": _int_or_none(record[73:75]),
        "second": _int_or_none(record[75:77]),

        "spare": record[77:80],

        "_line_number": line_number,
        "_raw": record,
    }


def read_segp1(filepath: str | Path) -> dict[str, Any]:
    """
    Read an ASCII SEG P1 (1983) positioning file.

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path to the ASCII SEG P1 file.

    Returns
    -------
    dict
        Dictionary containing:

        {
            "filepath": str,

            "header": {
                "records": [...],
                "text": str
            },

            "data": {
                "line_name": np.ndarray,
                "shotpoint": np.ndarray,
                "reshoot_code": np.ndarray,
                "latitude": np.ndarray,
                "longitude": np.ndarray,
                "easting": np.ndarray,
                "northing": np.ndarray,
                "elevation": np.ndarray,
                "year": np.ndarray,
                "day_of_year": np.ndarray,
                "hour": np.ndarray,
                "minute": np.ndarray,
                "second": np.ndarray,
                "spare": np.ndarray
            },

            "warnings": [...]
        }
    """

    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(
            f"SEG P1 file does not exist: {filepath}"
        )

    if not filepath.is_file():
        raise ValueError(
            f"Path is not a file: {filepath}"
        )

    header_records = []
    parsed_records = []
    warnings = []

    with filepath.open(
        "r",
        encoding="ascii",
        errors="replace",
        newline=None,
    ) as file:

        for line_number, raw_line in enumerate(file, start=1):

            # Remove newline characters only.
            # Do not strip spaces because SEG P1 is fixed-width.
            line = raw_line.rstrip("\r\n")

            if not line:
                continue

            if len(line) != 80:
                warnings.append(
                    {
                        "line_number": line_number,
                        "record_length": len(line),
                        "message": (
                            "Expected an 80-character SEG P1 record, "
                            f"but found {len(line)} characters."
                        ),
                    }
                )

            record = line.ljust(80)[:80]

            # Header record
            if record[0] == "H":
                header_records.append(record)
                continue

            # Data record
            if record[0] == " ":
                parsed_records.append(
                    _parse_data_record(
                        record,
                        line_number=line_number,
                    )
                )
                continue

            # Unknown record type
            warnings.append(
                {
                    "line_number": line_number,
                    "record_length": len(line),
                    "message": (
                        "Record does not begin with 'H' or a blank."
                    ),
                    "raw": record,
                }
            )

    # -------------------------------------------------------------
    # Convert list of record dictionaries into dictionary of arrays
    # -------------------------------------------------------------

    data = {
        "line_name": np.array(
            [r["line_name"] for r in parsed_records],
            dtype=str,
        ),

        "shotpoint": np.array(
            [
                np.nan if r["shotpoint"] is None
                else r["shotpoint"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "reshoot_code": np.array(
            [
                "" if r["reshoot_code"] is None
                else r["reshoot_code"]
                for r in parsed_records
            ],
            dtype=str,
        ),

        "latitude": np.array(
            [
                np.nan
                if r["latitude"] is None
                or r["latitude"]["decimal_degrees"] is None
                else r["latitude"]["decimal_degrees"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "longitude": np.array(
            [
                np.nan
                if r["longitude"] is None
                or r["longitude"]["decimal_degrees"] is None
                else r["longitude"]["decimal_degrees"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "easting": np.array(
            [
                np.nan if r["easting"] is None
                else r["easting"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "northing": np.array(
            [
                np.nan if r["northing"] is None
                else r["northing"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "elevation": np.array(
            [
                np.nan if r["elevation"] is None
                else r["elevation"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "year": np.array(
            [
                np.nan if r["year"] is None
                else r["year"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "day_of_year": np.array(
            [
                np.nan if r["day_of_year"] is None
                else r["day_of_year"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "hour": np.array(
            [
                np.nan if r["hour"] is None
                else r["hour"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "minute": np.array(
            [
                np.nan if r["minute"] is None
                else r["minute"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "second": np.array(
            [
                np.nan if r["second"] is None
                else r["second"]
                for r in parsed_records
            ],
            dtype=float,
        ),

        "spare": np.array(
            [r["spare"] for r in parsed_records],
            dtype=str,
        ),

        "_line_number": np.array(
            [r["_line_number"] for r in parsed_records],
            dtype=int,
        ),

        "_raw": np.array(
            [r["_raw"] for r in parsed_records],
            dtype=str,
        ),
    }

    return {
        "filepath": str(filepath),

        "header": {
            "records": header_records,
            "text": "\n".join(header_records),
        },

        "data": data,

        "warnings": warnings,
    }