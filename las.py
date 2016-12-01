"""
AuraQI module for reading data stored in LAS format.

Author:   Wes Hamlyn
Created:  15-May-2015
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

class LASReader(object):
    """
    Python class for reading log data in LAS format.
    """
    
    def __init__(self, filename, null_subs=np.nan, max_lines=-1):
        """
        Constructor for LAS class
        """
        
        self.filename = filename
        self.block = 'none'
        self.vers_info = {}
        self.well_info = {}
        self.param_info = {}
        self.tops = {}
        self.other = {}
        self.curve_info = {}
        self.curves = {}
        self.curve_names = []
        
        # open LAS file for reading        
        self.fd = open(self.filename, 'r')
        
        # count number of lines in file and reset cursor to BOF
        self.num_lines = sum(1 for line in self.fd)
        self.fd.seek(0)
        
        # allow user to read a maximum number of lines in an LAS file in case
        # the files are massive and very slow to read
        if max_lines > 0:
            self.num_lines = max_lines
            
        # start reading data line by line
        for i in range(0, self.num_lines):
            
            self.buf = self.fd.readline()
            #print(self.buf[:-2])
            
            # determine which block in the LAS file the cursor is currently in
            if self.buf[0] == '~':  # Start of new section block in LAS file
                if self.buf[1] == 'V':
                    self.block = 'VERSION'
                    #print('*** In VERSION block ***')
                    continue
                
                elif self.buf[1] == 'W':
                    self.block = 'WELL'
                    #print('*** In WELL block ***')
                    continue
                
                elif self.buf[1] == 'P':
                    self.block = 'PARAMETER'
                    #print('*** In PARAMETER block ***')
                    continue
                
                elif self.buf[1].lower() == 't':
                    self.block = 'TOPS'
                    #print('*** In TOPS block ***')
                    continue
                
                elif self.buf[1] == 'O':
                    self.block = 'OTHER'
                    #print('*** In OTHER block ***')
                    continue
                
                elif self.buf[1] == 'C':
                    self.block = 'CURVE'
                    #print('*** In CURVE block ***')
                    continue
                
                elif self.buf[1] == 'A':
                    self.block = 'DATA'
                    #print('*** In DATA block ***')
                    #self.curve_names = self.buf.split()[0:]
                    
            
            elif self.buf.split(' ')[0] == '\n':  # skip blank lines
                #print('*** Found a blank line... skipping. ***')
                continue
            
            
            elif self.buf[0] == '#':  # skip comment lines
                #print('*** Found a comment line... skipping. ***')
                continue
            
            
            else:    # Update the various data dictionaries
                
                if self.block == 'VERSION':
                    #print('*** Updating VERSION dictionary ***')
                    self._get_version_info()
                        
                if self.block == 'WELL':
                    #print('*** Updating WELL dictionary ***')
                    self._get_well_info()
                    
                if self.block == 'PARAMETER':
                    #print('*** Updating WELL dictionary ***')
                    self._get_param_info()
                
                if self.block == 'TOPS':
                    #print('*** Updating TOPS dictionary ***')
                    self._get_tops()
                
                if self.block == 'OTHER':
                    #print('*** Updating OTHER dictionary ***')
                    self._get_other()
                
                if self.block == 'CURVE':
                    #print('*** Updating CURVE dictionary ***')
                    self._get_curve_info()
                
                if self.block == 'DATA':
                    #print('*** Reading LOG DATA ***')
                    self.log_data_start = i
                    break
        
        
        #  Read the log curves via the Numpy loadtxt function.  These values 
        #  will be stored in the self.curves dictionary
        
        #  first verify that the LAS file doesn't contain wrapped log values, this
        #  isn't supported by this module yet
        if self.vers_info['WRAP']['data'].lower() in set(['true', 'yes', 'y']):
            print('LAS file contains wrapped log values')
            print('This LAS reader does not support wrapped files')
            print('Aborting... ')
                    
        else:
            self._get_curve_data()
        
        #  Clean up some unnecessary object attributes
        del([self.buf, self.block, self.fd, self.log_data_start, self.num_lines])
        #del([self.curve_names])
    
    
    
    def _get_version_info(self):
        """
        Method to read data from the Version Information section of an LAS file
        """
        
        mnem = self.buf.split('.', 1)[0].strip()
        data = self.buf.split('.', 1)[1].split(':')[0].strip()
        desc = self.buf.split('.', 1)[1].split(':')[1].strip()
        
        self.vers_info[mnem] = {'data':data, 'desc':desc}
        
        
    
    def _get_well_info(self):
        """
        Method to read data from the Well Information section of an LAS file
        """
        mnem = self.buf.split('.', 1)[0].strip()
        unit = self.buf.split('.', 1)[1].split(' ', 1)[0].strip()
        data = self.buf.split('.', 1)[1].split(' ', 1)[1].split(':', 1)[0].strip()
        desc = self.buf.split('.', 1)[1].split(' ', 1)[1].split(':', 1)[1].strip()
        
        well_floats = set(['NULL', 'STRT', 'STOP', 'STEP'])
        if mnem in well_floats:
            data = float(data)
        
        self.well_info[mnem] = {'data':data, 'unit':unit, 'desc':desc}
        
        
    def _get_param_info(self):
        """
        Method to read data from the Parameter Information section of an LAS file
        """
        mnem = self.buf.split('.', 1)[0].strip()
        unit = self.buf.split('.', 1)[1].split(' ', 1)[0].strip()
        data = self.buf.split('.', 1)[1].split(' ', 1)[1].split(':', 1)[0].strip()
        desc = self.buf.split('.', 1)[1].split(' ', 1)[1].split(':', 1)[1].strip()
        
        self.param_info[mnem] = {'data':data, 'unit':unit, 'desc':desc}
        
        
    def _get_tops(self):
        """
        Method to read tops from the unofficial Tops section of an LAS file
        """
        top = self.buf.split()[0]
        depth = float(self.buf.split()[1])
        
        self.tops[top] = depth
        
    
    def _get_other(self):
        """
        Method to read data from the Other Information section of an LAS file
        """
        txt = 'LINE %i' % (len(self.other)+1)
        self.other[txt] = self.buf
    
    
    def _get_curve_info(self):
        """
        Method to read the Curve Information section of an LAS file
        """

        mnem = self.buf.split('.', 1)[0].strip()
        unit = self.buf.split('.', 1)[1].split(' ', 1)[0].strip()
        data = self.buf.split('.', 1)[1].split(' ', 1)[1].split(':', 1)[0].strip()
        desc = self.buf.split('.', 1)[1].split(' ', 1)[1].split(':', 1)[1].strip()
        
        self.curve_info[mnem] = {'api_code':data, 'unit':unit, 'desc':desc}
        self.curve_names.append(mnem)
    
    def _get_curve_data(self):
        """
        Method to read curve data from the Data section of an LAS file
        """
        tmp = np.loadtxt(self.filename, dtype='float', skiprows=self.log_data_start)
        idx = np.nonzero(tmp == self.well_info['NULL']['data'])
        tmp[idx] = np.nan
        for i in range(0, len(self.curve_names)):
            self.curves[self.curve_names[i]] = tmp[:,i]

