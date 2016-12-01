"""
AuraQI module for building and using indexing schemes for fast random access to 
data in SEG-Y format.

Author:   Wes Hamlyn
Created:   20-Nov-2016
Last Mod:  1-Dec-2016

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
import auralib.segy as segy
import os
from struct import pack, unpack

class SegyIndex(object):
    """
    Class for indexing and reading 3D CMP gathers in SEG-Y format.
    """
    
    def __init__(self, idx_file):
        """
        Initialize a SegyIndex object.
        """
        
        self.idx_file = idx_file
        self.vers_string = ''
        self.segy_file = ''
        self.il_min = -1
        self.xl_min = -1
        self.il_num = -1
        self.xl_num = -1
        
        if os.path.exists(self.idx_file):
            self.read_index_headers()
        
            
    def read_index_headers(self):
        
        with open(self.idx_file, 'rb') as fd:
            
            # read version string
            buf = fd.read(32)
            buf = unpack('32s', buf)[0]
            self.vers_string = buf.strip('\x00')
            
            # read path to segy file
            buf = fd.read(1024)
            buf = unpack('1024s', buf)[0]
            self.segyfile = buf.strip('\x00')
            
            # read min inline
            buf = fd.read(4)
            self.il_min = unpack('l', buf)[0]
            
            # read min xline
            buf = fd.read(4)
            self.xl_min = unpack('l', buf)[0]
            
            # read numer of inlines
            buf = fd.read(4)
            self.il_num = unpack('l', buf)[0]
            
            # read number of xlines
            buf = fd.read(4)
            self.xl_num = unpack('l', buf)[0]

    
    def build_segy_index(self, segy_file, def_bhead, def_thead):
        """
        Given a SEG-Y file and the auralib binary and trace header definitions
        build an index file.
        """
        
        self.segy_file = segy_file 
        
        #  Create the SEG-Y file object
        buf = segy.Segy(self.segy_file, def_bhead, def_thead)

        #  read defined headers from all SEG-Y traces
        trc_start = 0
        trc_stop = buf.num_traces
        thead = buf.read_thead2(trc_start, trc_stop, verbose=1000)
        
        #  build Python list of sequential trace numbers
        tnum = np.arange(0, buf.num_traces)
        tnum = tnum.tolist()
        
        #  assign trace header struct to lists for convenience
        il = np.array(thead['il'])
        xl = np.array(thead['xl'])
        
        #  get inline and crossline statistics
        il_min = min(il)
        il_max = max(il)
        il_num = il_max - il_min + 1
        
        xl_min = min(xl)
        xl_max = max(xl)
        xl_num = xl_max - xl_min + 1
        
        #  This section of code calculates the trace number in the SEG-Y file 
        #  that each cdp begins at.  The number of traces in that cdp (fold) is
        #  also calculated.
        
        fold_il = []
        fold_xl = []
        fold_count = []
        cdp_start_trace = []

        for i in set(il):
            for x in set(xl):
                
                idx = np.nonzero( (il==i) & (xl==x) )[0]
                fold_il.append(i)
                fold_xl.append(x)
                fold_count.append(len(idx))
                
                if len(idx) > 0:
                    cdp_start_trace.append(idx[0])
                else:
                    cdp_start_trace.append(-1)
        
        for i in xrange(len(fold_il)):
            print('IL: %d  XL: %d  SeqTr: %i  Fold: %i' % 
                  (fold_il[i], fold_xl[i], cdp_start_trace[i], fold_count[i]) )
            
        # This section of code writes the version 2 segy index file
        
        print('Writing segy index file...')
        with open(self.idx_file, 'wb') as fd:
            
            fd.seek(0)
        
            buf = pack('32s', 'AuraSegyIndex_v01               ')
            fd.write(buf)
            
            buf = pack('1024s', self.segy_file)
            fd.write(buf)
            
            buf = pack('l', il_min)
            fd.write(buf)
            
            buf = pack('l', xl_min)
            fd.write(buf)
            
            buf = pack('l', il_num)
            fd.write(buf)
            
            buf = pack('l', xl_num)
            fd.write(buf)
            
            counter = 0
            for i in xrange( len(cdp_start_trace) ):
                counter += 1
                if counter == 100:
                    print('Writing cdp %i of %i' %(i, len(cdp_start_trace)))
                    counter = 1
                    
                buf = pack('ll', cdp_start_trace[i], fold_count[i])
                fd.write(buf)
                
        print('... DONE!')
        
        
    def get_index_tnums(self, il, xl):
        
        with open(self.idx_file, 'rb') as fd:
            # make sure cursor is at beginning of file
            fd.seek(0)
            
            # seek to the start trace index in the index file
            offset = 32 + 1024 + 16 + (il-self.il_min)*self.xl_num*8 + (xl-self.xl_min)*8
            fd.seek(offset)
            
            buf = fd.read(8)
            buf = unpack('ll', buf)
            cdp_start_trace = buf[0]
            cdp_fold = buf[1]
            
        return cdp_start_trace, cdp_fold
    
    
    def get_trace_data(self, il, xl, def_bhead, def_thead):
        
        trace0, num_traces = self.get_index_tnums(il, xl)
        
        start_trace = trace0 - 1
        
        end_trace = start_trace + num_traces
        print(self.segy_file)
        buf = segy.Segy(self.segy_file, def_bhead, def_thead)
        tdata = buf.read_multi_trace_data(start_trace, end_trace)
        
        return tdata
        
    
    def get_trace_headers(self, il, xl, def_bhead, def_thead):
        
        with open(self.idx_file, 'rb') as fd:
            # make sure cursor is at beginning of file
            fd.seek(0)
            
            # read path to segy file
            buf = fd.read(1024)
            buf = unpack('1024s', buf)[0]
            self.segyfile = buf.strip('\x00')

        
        trace0, num_traces = self.get_index_tnums(il, xl)
        start_trace = trace0 - 1
        end_trace = start_trace + num_traces

        buf = segy.Segy(self.segy_file, def_bhead, def_thead)
        thead = buf.read_thead2(start_trace, end_trace)
        
        return thead
