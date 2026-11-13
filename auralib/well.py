import numpy as np
import matplotlib.pyplot as plt
import auralib as aura
import pandas as pd
from scipy.interpolate import interp1d, PchipInterpolator
from copy import copy


class Project():
    """
    Auralib well project class
    """
    
    def __init__(self):
        """
        Constructor method
        """
        self.wells = {}
        self.uwis = []
        
    def add_well(self, well):
        """
        Add an auralib.well.Well object to the project
        """
        self.wells[well.uwi] = well
        self.uwis = list(self.wells.keys())
        self.uwis.sort()
        
    def delete_well(self, uwi):
        """
        Delete an auralib.well.Well object from the project
        """
        del(self.wells[uwi])
        self.uwis = list(self.wells.keys())
        self.uwis.sort()
    
    def plot_map(self, s=30, c='k', labels=False, **kwargs):
        """
        Plot a very basic basemap to see where wells are located
        """
        fig = plt.figure(num=1)
        fig.clf()
        ax1 = fig.add_subplot(111)
        
        for uwi in self.uwis:
            well = self.wells[uwi]
            if well.pos != None:
                ax1.plot(well.pos.x, well.pos.y, c, lw=1, *kwargs)
                ax1.scatter(well.pos.x[-1], well.pos.y[-1], s, c, *kwargs)
                
                if labels:
                    ax1.annotate(uwi, [well.pos.x[-1], well.pos.y[-1]], *kwargs)
                    
            else:
                ax1.scatter(well.topx, well.topy, s, c, *kwargs)
                if labels:
                    ax1.annotate(uwi, [well.topx, well.topy], *kwargs)
            
        ax1.set_aspect('equal')
        
        plt.show()
        
        return ax1
    

class Log():
    """
    Auralib well log class
    """
    
    def __init__(self, name=None, data=None, zref=None, ztype=None, units=None,
                 kind=None):
        self.name = name
        self.data = data
        self.zref = zref
        self.ztype = ztype
        self.units = units
        self.kind = kind
    
    def get_subset(self, zmin, zmax):
        mask = (self.zref>=zmin) & (self.zref<=zmax)
        a = Log(self.name, self.data[mask], self.zref[mask], self.ztype, self.units,
                self.kind)
        return a
        


class Position():
    """
    Auralib position set class
    """
    
    def __init__(self, x, y, z):
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)


class Survey():
    """
    Auralib directional survey class
    """
    
    def __init__(self, md, inc, azm):
        self.md = md
        self.inc = inc
        self.azm = azm


class Well():
    """
    Auralib well class
    """
    
    def __init__(self, uwi=None, wellname=None, topx=None, topy=None,
                 kb=None, gl=None, logs=None, tops=None, dev=None):
        self.uwi = uwi
        self.wellname = wellname
        self.topx = topx
        self.topy = topy
        self.kb = kb
        self.gl = gl
        self.tops = {}
        self.logs = {}
        self.pos = None
        self.svy = None

    def add_top(self, name, zref, ztype="md"):
        """
        Add a top to a well
        """
        self.tops[name] = {'zref': zref, 'ztype': ztype}
    
    def del_top(self, name):
        """
        Delete a top from a well
        """
        if hasattr(self.tops, name):
            exec(f'del(self.tops.%s)' % (name))
        else:
            print('Top %s not found in well %s' % (name, self.uwi))
    
    def tops_dataframe(self):
        """
        Return a Pandas DataFrame containing top data
        """
        tops = pd.DataFrame(self.tops).T.sort_values('zref')
        return tops
    
    def add_log(self, name, data, zref, ztype=None, units=None, kind=None):
        """
        Add an Auralog object to a well
        """
        curlog = Log(name, data, zref, ztype, units, kind)
        self.logs[name] = curlog
    
    def del_log(self, name):
        """
        Delete an Auralog object from a well
        """
        
        if hasattr(self.logs, name):
            exec(f'del(self.logs.%s)' % (name))
        else:
            print('Log %s not found in well %s' % (name, self.uwi))
    
    def add_position_set(self, x, y, z):
        """
        Add a position set to the well
        """
        
        self.pos = Position(x, y, z)
    
    def md_to_tvd(self, md):
        """
        Convert an md value to a TVD value
        """
        
        
                