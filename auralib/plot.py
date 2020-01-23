"""
auralib module containing plotting functions and related stuff...

Author:   Wes Hamlyn
Created:  20-Jan-2017
Last Mod: 20-Aug-2016

"""


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import MultiCursor
import numpy as np


def t2xy(d1, d2, d3, norm=True):
    
    c1 = np.array([0, 0])
    c2 = np.array([1, 0])
    c3 = np.array([0.5, 0.866])

    d1 = np.array(d1)
    d2 = np.array(d2)
    d3 = np.array(d3)
    
    # apply normalization
    if norm:
        datasum = np.sum(np.vstack([d1, d2, d3]), axis=0)
        d1 = d1/datasum
        d2 = d2/datasum
        d3 = d3/datasum

    px = c1[0]*d1 + c2[0]*d2 + c3[0]*d3
    py = c1[1]*d1 + c2[1]*d2 + c3[1]*d3
    
    return px, py


def tern(ax, lbls=['C1', 'C2', 'C3']):
    
    #  Corner points of triangular axes
    c1 = np.array([0, 0])
    c2 = np.array([1, 0])
    c3 = np.array([0.5, 0.866])
    
    #  Draw axes and add labels
    axbg_patch = mpatches.Polygon(np.vstack([c1, c2, c3]), closed=True, 
                                  fc='white', ec=None, zorder=1)
    ax.add_patch(axbg_patch)
    
    ax.plot([0, 1, 0.5, 0], [0, 0, 0.866, 0], 'k', lw=1.5, zorder=5)
    ax.text(c1[0], c1[1], lbls[0], ha='right', va='top', fontsize=14)
    ax.text(c2[0], c2[1], lbls[1], ha='left', va='top', fontsize=14)
    ax.text(c3[0], c3[1], lbls[2], ha='center', va='bottom', fontsize=14)
    
    #  Draw gridlines
    for i in np.arange(0.1, 1, 0.1):
        
        l1 = [i, i]
        l2 = 1.0 - i
        lx, ly = t2xy(l1, [0, l2], [l2, 0])
        ax.plot(lx, ly, ':', lw=1.0, color=u'0.4', zorder=2)
        ax.text(lx[-1]+0.01, ly[-1]-0.03, '%.2f' % i, ha='center', va='center', rotation=-60.0)

        l1 = [i, i]
        l2 = 1.0 - i
        lx, ly = t2xy([0, l2], l1,  [l2, 0])
        ax.plot(lx, ly, ':', lw=1.0, color=u'0.4', zorder=2)
        ax.text(lx[0]+0.005, ly[0]+0.03, '%.2f' % i, ha='left', va='center', rotation=60.0)

        l1 = [i, i]
        l2 = 1.0 - i
        lx, ly = t2xy([0, l2], [l2, 0], l1)
        ax.plot(lx, ly, ':', lw=1.0, color=u'0.4', zorder=2)
        ax.text(lx[-1]-0.01, ly[0], '%.2f' % i, ha='right', va='center')
    
    
    
    
    ax.set_xlim([-0.1, 1.1])
    ax.set_ylim([-0.1, 0.966])
    ax.set_aspect('equal')
    
    ax.set_axis_off()
    
    #ax.xaxis.set_visible(False)
    #ax.yaxis.set_visible(False)
    
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    

def tern_scatter(ax, d1, d2, d3, s=25, color=None, marker=None, cmap=None, 
                 lw=None, ec=None, alpha=1.0, label=None):
    
    #  Transform points from XY -> C1, C2, C3 coordinate system
    px, py = t2xy(d1, d2, d3)
    
    #  Plot points on
    pts = ax.scatter(px, py, s=s, c=color, marker=marker, cmap=cmap, lw=lw, 
                     edgecolor=ec, alpha=alpha, label=label,
                     zorder=10)
    
    return pts


def tern_line(ax, d1, d2, d3, c=None, lw=None, ls=None, label=None):
    
    #  Transform points from XY -> C1, C2, C3 coordinate system
    px, py = t2xy(d1, d2, d3)
    
    #  Plot points on
    hdl = ax.plot(px, py, c=c, lw=lw, ls=ls, label=label, zorder=10)
    
    return hdl


def plot_blocky(ax, data, zdata, linespec='b-', lw=1):
    """
    Convenience function for plotting a blocky log.
    
    Ensure that the zdata log has 1 more sample than the data log.
    """
    
    for i in range(0, len(data)):
        ax.plot([data[i], data[i]], [zdata[i], zdata[i+1]], linespec, lw=lw)
        
    for i in range(1, len(data)):
        ax.plot([data[i-1], data[i]], [zdata[i], zdata[i]], linespec, lw=lw)


def radar(ax, data, data_names, lw=1.0, color='k', ls='-', marker=None, label=None):
    """
    function to produce a radar plot
    
    ax = axis handle of a polar axis to do plotting in
    data = 1D array or list of data values
    data_names = 1D list of data names
    
    """
    
    # get number of values in data vector
    N = len(data)
        
    # append the first data value to the end of the list so we can make closed
    # polygonal regions for plotting and filling
    data = np.array(data)
    data = data.tolist()
    data.append(data[0])
    
    # What will be the angle of each axis in the plot?
    # (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]
     
    # Draw one axe per variable + add labels labels yet
    plt.sca(ax)
    plt.xticks(angles[:-1], data_names)
     
    # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks(color="grey", size=8)
     
    # Plot data
    ax.plot(angles, data, lw=lw, color=color, ls=ls, marker=marker, label=label)
    
    # Fill area
    ax.fill(angles, data, color=color, alpha=0.1)
    ax.grid(True, ls=':')



def plot_filled_logs(ax, logs, z, fill_colors, labels, lw=1.0, alpha=0.3):
    """
    Plot a series of logs using fill colors between each log. Designed to show
    a series of mineral fraction logs or fluid saturation logs
    """
    
    import matplotlib as mpl
    
    nlogs = len(logs)
    
    # take the log fill colors, and make them darker to be used for the log 
    # line colors. This will make them more prominent in the display
    log_colors = []
    for color in fill_colors:
        rgb = mpl.colors.to_rgb(color)
        rgb = np.array(rgb)
        rgb = rgb*0.3
        log_colors.append(rgb)
        
    # first plot the log fills
    logsum = 0 # cumulative log value
    for i in range(nlogs):
        curlog = logs[i] + logsum
        ax.fill_betweenx(z, curlog, logsum, where=curlog>=logsum, 
                         facecolor=fill_colors[i], alpha=alpha, label=labels[i])
        logsum = logsum + logs[i]
        
    # next plot the log curves to visually separate the filled areas
    logsum = 0 # cumulative log value
    for i in range(nlogs):
        curlog = logs[i] + logsum
        ax.plot(curlog, z, c=log_colors[i], lw=lw, alpha=alpha)
        logsum = logsum + logs[i]



def format_log_axes(axes, ylabel):
    """
    Function to format a series of axes displaying simple log curves.
    """
    
    for ax in axes:
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.grid(True, ls=':')
    
    for ax in axes[1:-1]:
        plt.setp(ax.get_yticklabels(), visible=False)
    
    axes[0].set_ylabel(ylabel)
    axes[-1].set_ylabel(ylabel)
    axes[-1].yaxis.set_ticks_position('right')
    axes[-1].yaxis.set_label_position('right')\
    
def make_wellview(ntrack=5, figsize=(5, 5)):
    """
    Function for creating a blank, multi-track well viewer with a well header
    area, a well log header area for each track, and a well log data area.
    """
    
    fig = plt.figure(num=1, figsize=figsize)
    fig.clf()

    nrow = 30
    ncol = ntrack
    ttl = 0
    hdr = 1
    dat = 4
    
    axd = [plt.subplot2grid((nrow, ncol), (dat, 0), rowspan=nrow-3)]
    axh = [plt.subplot2grid((nrow, ncol), (hdr, 0), rowspan=3, sharex=axd[0])]
    
    for i in range(1, ncol):
        axd.append(plt.subplot2grid((nrow, ncol), (dat, i), rowspan=nrow-2, sharey=axd[0]))
        axh.append(plt.subplot2grid((nrow, ncol), (hdr, i), rowspan=3, sharex=axd[i], sharey=axh[0]))
        
    axttl = plt.subplot2grid((nrow, ncol), (ttl, 0), colspan=ncol)
    
    for ax in fig.get_axes():
        ax.tick_params(which='both', direction='in')
        
    curs = MultiCursor(fig.canvas, axd, vertOn=False, horizOn=True,
                       lw=1, c=u'0.3')

    wview = {'fig': fig, 'axd': axd, 'axh': axh, 'axttl': axttl, 'curs': curs}
    
    return wview


def plot_log(wview, tracknum, aura_log, fmtstr='%.1f', numlogs=1, logpos=1, xscale='normal'):
    """
    Function to plot logs in a Well Viewer created using the 
    aura.plot.make_wellvew() function.
    
    Input logs are required to be aura.well.AuraLog() objects.
    
    """
    axd = wview['axd']
    axh = wview['axh']
    
    if xscale=='log':
        axd[tracknum].set_xscale('log')
        axh[tracknum].set_xscale('log')
        
    axd[tracknum].plot(aura_log.data, aura_log.zref, 
                       color=aura_log.c, lw=aura_log.lw)
    axd[tracknum].set_xlim(aura_log.plt_range)
    
    xlim = axd[tracknum].get_xlim()
    
    data_range = np.abs(xlim[1] - xlim[0])
    
    if xscale == 'log':
        data_lim_offset0 =  xlim[0]
        data_lim_offset1 =  xlim[1]
    else:
        data_lim_offset0 =  data_range*0.02 + xlim[0]
        data_lim_offset1 = -data_range*0.02 + xlim[1]
    
    axh[tracknum].plot([data_lim_offset0, data_lim_offset1], (logpos, logpos), 
                       c=aura_log.c, lw=aura_log.lw)
    
    bbox = dict(fc='white', ec='white', alpha=1.0, pad=0.5)
    
    axh[tracknum].text(data_lim_offset0, logpos, fmtstr % xlim[0], 
                       va='top', ha='left', color=aura_log.c, 
                       bbox=bbox, fontsize=aura_log.fs)
    
    axh[tracknum].text(data_lim_offset1, logpos, fmtstr % xlim[1], 
                       va='top', ha='right', color=aura_log.c, 
                       bbox=bbox, fontsize=aura_log.fs)
    
    if xscale=='log':
        xpos_logname = np.sqrt(np.cumsum(xlim))
    else:
        xpos_logname = np.mean(xlim)
        
    if len(aura_log.units)>0:
        axh[tracknum].text(xpos_logname, logpos, aura_log.name+' ('+aura_log.units+')', 
                           va='bottom', ha='center', color=aura_log.c, 
                           fontsize=aura_log.fs)
    else:
        axh[tracknum].text(xpos_logname, logpos, aura_log.name, 
                           va='bottom', ha='center', color=aura_log.c, 
                           fontsize=aura_log.fs)
    
    axh[tracknum].set_ylim([0, numlogs+1])
    axh[tracknum].set_xlim(xlim)
    
    


def format_wellview(wview, ylabel='Depth', title_text='Title', ylim='none'):
    """
    Once all well data are plotted in an aura wellviewer figure, call this
    function to make the well viewer figure look nice.
    """

    fig = wview['fig']
    axd = wview['axd']
    axh = wview['axh']
    axttl = wview['axttl']
    
    axd[0].invert_yaxis()
    if ylim != 'none':
        axd[0].set_ylim(ylim)
        
    ntrack = len(axd)
    count = 1
    for (axdi, axhi) in zip(axd, axh):
        
        axdi.grid(True, which='major', ls=':', lw=0.5)
        axdi.grid(True, which='minor', ls=':', lw=0.5)
        axdi.minorticks_on()
        axdi.xaxis.set_ticks_position('top')
        axdi.xaxis.set_label_position('top')
        [label.set_visible(False) for label in axdi.get_xticklabels()]
        
        axhi.set_facecolor('white')
        axhi.set_frame_on(True)
        axhi.grid(False)
        axhi.xaxis.set_visible(False)
        axhi.yaxis.set_visible(False)

        axdi.yaxis.set_ticks_position('both')
        
        if (count==1):
            axdi.set_ylabel(ylabel)
        
        elif count==ntrack:
            axdi.set_ylabel(ylabel)
            axdi.yaxis.set_ticks_position('right')
            axdi.yaxis.set_label_position('right')
            
        else:
            axdi.tick_params(labelright=True)
            [label.set_visible(False) for label in axdi.get_yticklabels()]
            
        count += 1
    
    axttl.set_facecolor('#fcfcc4')
    axttl.set_frame_on(True)
    axttl.grid(False)
    axttl.xaxis.set_visible(False)
    axttl.yaxis.set_visible(False)
    axttl.text(0.5, 0.5, title_text, ha='center', va='center', weight='normal')
    axttl.set_xlim([0, 1])
    axttl.set_ylim([0, 1])
    
    fig.tight_layout(w_pad=0.00, h_pad=0.00)


def add_van_krevelan_template(ax, lw=2, fs=14, c='k'):
    T1 = np.array([[  3.89993785,  97.72755495],
                   [  4.27284027, 140.94126413],
                   [  3.89993785, 182.6648454 ],
                   [  3.89993785, 227.36868248],
                   [  3.89993785, 260.15149634],
                   [  3.89993785, 291.44418229],
                   [  4.27284027, 324.22699615],
                   [  4.27284027, 376.38147274],
                   [  4.27284027, 419.59518192],
                   [  4.27284027, 456.84837949],
                   [  5.39154755, 489.63119334],
                   [  5.01864512, 514.96336769],
                   [  6.13735239, 559.66720477],
                   [  7.25605966, 601.39078604],
                   [  7.62896209, 629.70321619],
                   [  8.00186451, 659.50577425],
                   [ 10.61218148, 701.22935552],
                   [ 12.84959602, 735.50229728],
                   [ 15.08701057, 769.77523904],
                   [ 16.95152268, 793.61728548],
                   [ 20.68054692, 824.90997144],
                   [ 23.29086389, 856.20265739],
                   [ 25.90118086, 877.06444803],
                   [ 28.1385954 , 897.92623867],
                   [ 32.61342449, 921.76828511],
                   [ 37.08825357, 941.13994785],
                   [ 40.07147296, 962.00173848],
                   [ 43.8004972 , 976.90301751]])
    
    T2 = np.array([[  8.37476694,  99.21768285],
                   [  8.37476694, 123.05972929],
                   [  8.00186451, 160.31292686],
                   [  8.00186451, 202.03650813],
                   [  8.74766936, 242.26996151],
                   [  8.37476694, 286.97379858],
                   [  8.37476694, 321.24674035],
                   [  8.37476694, 355.51968211],
                   [  9.12057178, 382.34198435],
                   [ 10.61218148, 415.12479821],
                   [ 13.22249845, 450.88786788],
                   [ 16.95152268, 479.20029803],
                   [ 20.68054692, 498.57196076],
                   [ 24.78247359, 522.4140072 ],
                   [ 28.1385954 , 543.27579784],
                   [ 32.61342449, 567.11784428],
                   [ 36.34244873, 583.50925121],
                   [ 40.81727781, 601.39078604],
                   [ 49.02113114, 620.76244878],
                   [ 54.24176507, 635.6637278 ],
                   [ 60.58110628, 652.05513473],
                   [ 64.31013052, 660.99590215],
                   [ 69.53076445, 666.95641376],
                   [ 75.49720323, 675.89718117],
                   [ 80.71783717, 680.36756488]])
    
    T3 = np.array([[ 14.34120572,  39.61256675],
                   [ 18.07022996,  56.00397367],
                   [ 23.66376631,  67.9249969 ],
                   [ 34.47793661,  78.35589221],
                   [ 43.8004972 ,  85.80653173],
                   [ 51.25854568,  97.72755495],
                   [ 59.83530143, 108.15845027],
                   [ 67.29334991, 117.09921768],
                   [ 74.37849596, 124.5498572 ],
                   [ 83.32815413, 127.530113  ],
                   [ 96.00683654, 137.96100832],
                   [103.46488502, 142.43139203],
                   [115.39776259, 145.41164783],
                   [128.82224984, 152.86228735],
                   [138.89061529, 158.82279896],
                   [147.46737104, 163.29318266],
                   [159.77315103, 167.76356637],
                   [169.84151647, 172.23395008],
                   [179.90988191, 173.72407798],
                   [187.74083282, 178.19446169],
                   [193.33436917, 181.1747175 ]])
    
    L1 = np.array([[1.70913611e-01, 2.35936918e+00],
                   [1.99673710e+02, 9.99254936e+02]])
    
    L2 = np.array([[5.43816035e-01, 2.35936918e+00],
                   [1.99673710e+02, 5.98410530e+02]])
    
    ax.plot(T1[:, 0], T1[:, 1], color=c, lw=lw)
    ax.plot(T2[:, 0], T2[:, 1], color=c, lw=lw)
    ax.plot(T3[:, 0], T3[:, 1], color=c, lw=lw)
    ax.plot(L1[:, 0], L1[:, 1], color=c, ls='--', lw=lw)
    ax.plot(L2[:, 0], L2[:, 1], color=c, ls='--', lw=lw)
    
    
    ax.text(T1[-1, 0], T1[-1, 1], 'Type I', ha='left', va='center', fontsize=fs, color=c)
    ax.text(T2[-1, 0], T2[-1, 1], 'Type II', ha='left', va='center',  fontsize=fs, color=c)
    ax.text(T3[-1, 0], T3[-1, 1], 'Type III', ha='right', va='bottom',  fontsize=fs, color=c)
    ax.text(80, 35, 'Type IV', fontsize=fs, color=c)
    ax.text(110, 875, 'Oil', fontsize=fs, color=c)
    ax.text(157, 626, 'Mixed', fontsize=fs, color=c)
    ax.text(166, 382, 'Gas', fontsize=fs, color=c)
    
    
def plot_confusion_matrix(y_true, y_pred, classes, normalize=False, title=None,
                          cmap=plt.cm.Blues, fig=None, flip_text_color=False):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
	
    from sklearn.metrics import confusion_matrix
    from sklearn.utils.multiclass import unique_labels
    
    if not title:
        if normalize:
            title = 'Normalized Confusion Matrix'
        else:
            title = 'Confusion Matrix, Without Normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    
    # Only use the labels that appear in the data
    classes = classes[unique_labels(y_true, y_pred)]
    
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print('Normalized Confusion Matrix')
    else:
        print('Confusion Matrix, Without Normalization')

    print(cm)
    
    if not fig:
        fig = plt.figure(num=1)
    
    fig.clf()
    ax = fig.add_subplot(111)
    
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    plt.setp(ax.get_xticklabels(), rotation=-45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        
        for j in range(cm.shape[1]):
            
            if flip_text_color == True:
                ax.text(j, i, format(cm[i, j], fmt),
                        ha="center", va="center",
                        color="white" if cm[i, j] < thresh else "black")
                
            elif flip_text_color == False:
                ax.text(j, i, format(cm[i, j], fmt),
                        ha="center", va="center",
                        color="black" if cm[i, j] < thresh else "white")
    
    ax.grid(False)

    fig.tight_layout()
    
    return ax

