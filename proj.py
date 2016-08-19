"""
AuraQI module containing objects and functions to read data from an AuraQI
project structure.

Author:   Wes Hamlyn
Created:  25-Mar-2016
Last Mod: 17-Aug-2016
"""


import os
import numpy as np

def create(projpath):
    """
    Function to create a new empty project directory structure with the
    appropriate files initialized.
    """
    
    # Create the top-level project directory
    os.mkdir(projpath)
    
    # Create the project sub-directories
    os.mkdir(os.path.join(projpath, 'min'))
    os.mkdir(os.path.join(projpath, 'seis'))
    os.mkdir(os.path.join(projpath, 'seis', '2D'))
    os.mkdir(os.path.join(projpath, 'seis', '3D'))
    os.mkdir(os.path.join(projpath, 'well'))
    os.mkdir(os.path.join(projpath, 'wvlt'))
    os.mkdir(os.path.join(projpath, 'zone'))
    
    # Create the initial well_list.txt file
    head = "WELL,KB,GRD,X,Y"
    fd = open(os.path.join(projpath, 'well_list.txt'), 'w')
    fd.write(head)
    fd.close()


def get_well_heads(projpath):
    """
    Function to read the list of well headers into memory.
    """
    fpath = os.path.join(projpath, 'well_list.csv')
    fd = open(fpath, 'r')
    
    fd.close()


def write_blank_log(filename, nsamp, z0, dz, ltype='Misc'):
    """
    Function to write a blank log.
    """
        
    z1 = nsamp*dz + z0
    zref = np.arange(z0, z1, dz)
    data = zref * 0.0
    
    fd = open(filename, 'w')
    
    # write header
    fd.write('TYPE:%s\n' % ltype)
    fd.write('START:%f\n' % z0)
    fd.write('STEP:%f\n' % dz)
    
    # write log digits
    for i in data:
        buf = '%f\n' % i
        fd.write(buf)
    
    fd.close()
    


 
