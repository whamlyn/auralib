"""
auralib module with various functions that build on the segy module
functionality.
"""

import numpy as np
import matplotlib.pyplot as plt
import auralib as aura


class PseudoGather():
    """
    Class to generate pseudo gathers from partial stacks. Wraps much of the
    auralib.segy.Segy functionality for working with multiple volumes.
    """
    
    def __init__(self, infiles, angles, def_thead):
        
        self.infiles = infiles
        self.angles = angles
        self.num_trc = len(angles) # number of traces per CDP ensemble
        
        self.sgy_bufs = []
        for infile in infiles:
            self.sgy_bufs.append(aura.segy.Segy(infile, def_thead))
        
    def read_tdata(self, trcidx):
        
        tdata = []
        for buf in self.sgy_bufs:
            tdata_tmp = buf.read_tdata(trcidx)
            tdata.append(tdata_tmp)
        
        return tdata
    
    def calc_IG(self, trcidx):
        
        gdata = np.array(self.read_tdata(trcidx))
        
        na, ns = gdata.shape
        I = []
        G = []
        for i in range(ns):
            avo = gdata[:, i]
            Itmp, Gtmp = aura.avo.calc_IG(avo, self.angles)
            I.append(Itmp)
            G.append(Gtmp)
        
        I = np.array(I)
        G = np.array(G)
        
        return I, G
        