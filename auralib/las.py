"""
auralib module for reading data stored in LAS format.

Author:   Wes Hamlyn
Created:  15-May-2015
Last Mod: 17-Aug-2016

"""

import numpy as np

class LASReader(object):
    """
    Python class for reading log data in LAS format.
    """
    
    def __init__(self, filename, null_subs=np.nan, max_lines=-1, read_tops=True):
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
        self.read_tops = read_tops
        
        # open LAS file for reading
        with open(self.filename, 'r') as fd:
            self.dbuf = fd.readlines()
        
        # count number of lines in file
        self.num_lines = len(self.dbuf)
            
        # start parsing data line by line
        i = -1 # initialize a line counter
        for line in self.dbuf:
            
            i = i + 1    # increment counter
            
            self.buf = line
            self.buf = self.buf.strip() # remove leading/trailing whitespace                
            
            if len(self.buf)==0:
                # skip any blank lines
                continue
            
            elif self.buf[0] == '~':  # Start of new section block in LAS file
                # determine which block in the LAS file the cursor is currently in
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
                    continue
                
                elif self.buf[1] == 't':
                    self.block = 'TOPS'
                    #print('*** In Tops block ***')
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
                    self.log_data_start = i*1
                    #print('Log data starts on line %i' % self.log_data_start)
                    break
        
        
        #  Read the log curves via the Numpy loadtxt function.  These values 
        #  will be stored in the self.curves dictionary
        
        # Read the log digits using appropriate method for wrapped or unwrapped
        # LAS files
        if self.vers_info['WRAP']['data'].lower() in set(['true', 'yes', 'y']):
            self._get_curve_data_wrapped()
                    
        else:
            self._get_curve_data()
        
        #  Clean up some unnecessary object attributes
        #del([self.buf, self.block, self.fd, self.log_data_start, self.num_lines])
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
        Method to read tops from the unofficial Tops section of an LAS file.
        Spaces in top names are handled so long as they are single spaces.
        """
        if self.read_tops:
            
            line = self.buf.strip().split('  ')
            
            line_filt = []
            for each in line:
                if len(each)!=0:
                    line_filt.append(each)
                
            # reformat names/depth as appropriate
            cur_top_name = line_filt[0].strip()
            cur_top_depth = float(line_filt[1].strip())
            
            top = cur_top_name
            depth = cur_top_depth
            
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
        
        del(self.dbuf)
        
        
    def _get_curve_data_wrapped(self):
        """
        Method to read wrapped curve data from the Data section of an LAS file
        """
        ncurves = len(self.curve_names)
        data = self.dbuf[self.log_data_start:]
        dataflat = []
        for line in data:
            dataflat.extend(line.strip().split())
        
        dataflat = np.array(dataflat, dtype='float')
        dataflat = dataflat.reshape([-1, ncurves])
        
        idx = np.nonzero(dataflat == self.well_info['NULL']['data'])
        dataflat[idx] = np.nan
        
        for i in range(0, len(self.curve_names)):
            self.curves[self.curve_names[i]] = dataflat[:, i]
        
        del(self.dbuf)
       
        

