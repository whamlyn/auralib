"""
AuraQI module containing plotting functions and related stuff...

Author:   Wes Hamlyn
Created:  20-Jan-2017
Last Mod: 20-Aug-2016

Copyright 2017 Wes Hamlyn

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


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
        ax.plot(lx, ly, ':', color=u'0.4', zorder=5)
        ax.text(lx[-1]+0.01, ly[-1]-0.03, '%.2f' % i, ha='center', va='center', rotation=-60.0)

        l1 = [i, i]
        l2 = 1.0 - i
        lx, ly = t2xy([0, l2], l1,  [l2, 0])
        ax.plot(lx, ly, ':', color=u'0.4', zorder=5)
        ax.text(lx[0]+0.005, ly[0]+0.03, '%.2f' % i, ha='left', va='center', rotation=60.0)

        l1 = [i, i]
        l2 = 1.0 - i
        lx, ly = t2xy([0, l2], [l2, 0], l1)
        ax.plot(lx, ly, ':', color=u'0.4', zorder=5)
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
    

def tern_scatter(ax, d1, d2, d3, color='r', label=None):
    
    #  Transform points from XY -> C1, C2, C3 coordinate system
    px, py = t2xy(d1, d2, d3)
    
    #  Plot points on
    pts = ax.scatter(px, py, s=30, c=color, marker='o', lw=0, label=label,
                     zorder=4)
    
    return pts
