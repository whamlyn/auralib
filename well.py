"""
AuraQI module defining well data classes for the auralib data model.

Author:   Wes Hamlyn
Created:  17-Aug-2016
Last Mod: 17-Aug-2016
"""

import numpy as np

class Well():
    """
    Class for top-level well objects.
    """
    def __init__(self, wellid, type='Well'):
        self.type = type
        self.header = WellHeader()
        self.header.wellid = wellid
        
        self.logs = Logs()

        self.td = TD()


class WellHeader():
    """
    Class for well header data and operations
    """
    def __init__(self):
        pass


class Logs():
    """
    Class for top-level well log data and operations
    """
    def __init__(self):
        self.type = 'Logs container'


class Log():
    def __init__(self, data, units='', ztype='md', 
                 start=np.nan, stop=np.nan, step=np.nan):
        self.data = data
        self.units = units
        self.ztype = ztype
        self.start = start
        self.stop = stop
        self.step = step


class TD():
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