"""
auralib module defining well data classes for the auralib data model.

Author:   Wes Hamlyn
Created:  17-Aug-2016
Last Mod: 17-Aug-2016

"""

import numpy as np


class AuraWell():
    """
    Class for top-level well objects.
    """
    def __init__(self, wellid, type='Well'):
        self.type = type
        self.header = AuraWellHeader()
        self.header.wellid = wellid
        
        self.logs = AuraLogs()

        self.td = AuraTD()



class AuraWellHeader():
    """
    Class for well header data and operations
    """
    def __init__(self):
        pass



class AuraLogs():
    """
    Class for top-level well log data and operations
    """
    def __init__(self):
        self.type = 'Logs container'



class AuraLog():
    """
    Class to store log data.
    """
    
    def __init__(self, data, zref, name, ztype='md', units='', 
                 plt_range=[0, 0], c='k', lw=0.5, fs=8.0):
        """
        Constructor method...
        self.data = log digits
        self.zref = z-reference log (i.e. depth or time samples)
        self.name = textual name of well log (i.e. the log mnemonic)
        self.ztype = textual name of zref type (i.e. 'md', 'tvd', 'tvdss', 'twt')
        self.units = textual units string (i.e. 'ft/sec', 'API', 'g/cc', etc.)
        self.plt_range = min and max range of log in plot displays
        self.c = colour of log in plot displays
        self.lw = line width of log in plot displays
        self.fs = font size of log name in plot displays
        """
        
        self.data = np.array(data)
        self.zref = zref
        self.name = name
        self.ztype = ztype
        self.units = units
        self.plt_range = plt_range
        self.c = c
        self.lw = lw
        self.fs = fs



class AuraTD():
    """
    Class for top-level well time-depth data and operations
    """
    def __init__(self, tvd=[], twt=[]):
        
        # test to make sure there are the same number of time-depth pairs
        if len(tvd) == len(twt):
            self.tvd = tvd
            self.twt = twt
        else:
            print('Time and depth arrays have different numbers of samples!')



class AuraMarkers():
    """
    Class for storing well marker data.  Very simple for the moment, need to
    refine this for more rigorous usage.
    """
    def __init__(self):
        self.name = []
        self.well = []
        self.md = []
        self.tvdss = []
        self.utmx = []
        self.utmy = []


