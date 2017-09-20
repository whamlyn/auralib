"""
AuraQI module defining well data classes for the auralib data model.

Author:   Wes Hamlyn
Created:  17-Aug-2016
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
    """
    Class for regularly sampled (periodic) well log data.
    self.data = log sample values 
    self.units = unit string corresponding to log digits
    self.ztype = depth reference type (MD, TVDkb, TVDss, TVDml, TWT
                 OWT
    """
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

            
class Markers():
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
