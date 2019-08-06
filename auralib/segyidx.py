"""
AuraQI module for building and using indexing schemes for fast random access to 
data in SEG-Y format.

Author:   Wes Hamlyn
Created:   20-Nov-2016
Last Mod:  1-Dec-2016

"""

import numpy as np
import auralib.segy as segy
import os
from struct import pack, unpack
from time import time


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
            self._read_index_headers()
        
            
    def _read_index_headers(self):
        
        with open(self.idx_file, 'rb') as fd:
            
            # read version string
            buf = fd.read(32)
            buf = unpack('32s', buf)[0]
            self.vers_string = buf.strip('\x00'.encode())
            
            # read path to segy file
            buf = fd.read(1024)
            buf = unpack('1024s', buf)[0]
            self.segyfile = buf.strip('\x00'.encode())
            
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

    
    def build_segy_index(self, segy_file, def_thead):
        """
        Given a SEG-Y file and the auralib trace header definitions
        build an index file.
        
        trc, fold = build_segy_index(segy_file, def_thead)
        
        Inputs:
            segy_file - Path to segy file to be indexed
            def_thead - Dictionary of trace header kewords and format. Requires
                        keywords 'il' and 'xl' to perform indexing.
        
        Outputs:
            None
            
		
		Known Issues: 
		1) Apparently having an inline or crossline step other than 1 will cause
    	   this to break.  -Wes H. Sept 2017
        """
        
        self.segy_file = segy_file 
        
        #  Create the SEG-Y file object
        buf = segy.Segy(self.segy_file, def_thead)

        #  read defined headers from all SEG-Y traces
        trc_start = 0
        trc_stop = buf.num_traces
        thead = buf.read_thead_multi(trc_start, trc_stop, verbose=1000)
        print('Finished reading SEG-Y headers.')
        
        #  build Python list of sequential trace numbers
        tnum = np.arange(0, buf.num_traces)
        tnum = tnum.tolist()
        print('Finished building sequential trace numbers.')
        
        #  convert trace header dictionary to numpy arrays for convenience
        print('Converting trace header lists to numpy arrays...')
        tic = time()
        il = np.array(thead['il'])
        xl = np.array(thead['xl'])
        toc = time() - tic
        print('Done! (%e seconds)' % toc)
        
        #  get inline and crossline statistics
        print('Getting Inline and Crossline statistics...')
        self.il_min = min(il)
        self.il_max = max(il)
        self.il_num = self.il_max - self.il_min + 1
        
        self.xl_min = min(xl)
        self.xl_max = max(xl)
        self.xl_num = self.xl_max - self.xl_min + 1
        print('Done!')
        
        #  This section of code calculates the trace number in the SEG-Y file 
        #  that each cdp begins at.  The number of traces in that cdp (fold) is
        #  also calculated.
        
        print('Starting the loops for calculating the fold in each CDP')
        
        fold_il = []
        fold_xl = []
        fold_count = []
        cdp_start_trace = []

        for i in set(il):
            print('Working on inline: %i of %i' % (i, self.il_max))
            
            for x in set(xl):
                
                idx = np.nonzero( (il==i) & (xl==x) )[0]
                
                fold_il.append(i)
                fold_xl.append(x)
                fold_count.append(len(idx))
                
                if len(idx) > 0:
                    cdp_start_trace.append(idx[0])
                else:
                    cdp_start_trace.append(-1)
        
        
        # Set to True to have detailed info printed to command line
        if False:
            for i in range(len(fold_il)):
                print('IL: %d  XL: %d  SeqTr: %i  Fold: %i' % 
                      (fold_il[i], fold_xl[i], cdp_start_trace[i], fold_count[i]) )
        
        print('Done with CDP fold calculations!')
        
        
        # This section of code writes the version 2 segy index file
        
        print('Writing segy index file...')
        with open(self.idx_file, 'wb') as fd:
            
            fd.seek(0)
            
            print('Writing version string')
            buf = pack('32s', 'AuraSegyIndex_v01'.encode())
            fd.write(buf)
            
            print('Writing segy file path')
            buf = pack('1024s', self.segy_file.encode())
            fd.write(buf)
            
            buf = pack('l', self.il_min)
            fd.write(buf)
            
            buf = pack('l', self.xl_min)
            fd.write(buf)
            
            buf = pack('l', self.il_num)
            fd.write(buf)
            
            buf = pack('l', self.xl_num)
            fd.write(buf)
            
            counter = 0
            for i in range( len(cdp_start_trace) ):
                counter += 1
                if counter == 101:
                    print('Writing cdp %i of %i' %(i, len(cdp_start_trace)))
                    counter = 1
                   
                
                buf = pack('ll', cdp_start_trace[i], fold_count[i])
                fd.write(buf)
                
        print('... DONE!')

        
    def build_segy_index_v2(self, segy_file, def_thead):
        """
        Given a SEG-Y file and the auralib binary and trace header definitions
        build an index file.
        
        This method assumes data are sorted by inline and then crossline and 
        just does a straight running count through the headers rather than 
        using the np.nonzero() approach which becomes quite slow on large
        arrays.
        
        This method is still in development, needs to be finalized and tested
        -Wes
        """
        
        self.segy_file = segy_file 
        
        #  Create the SEG-Y file object
        buf = segy.Segy(self.segy_file, def_thead)

        #  read defined headers from all SEG-Y traces
        trc_start = 0
        trc_stop = buf.num_traces
        thead = buf.read_thead_multi(trc_start, trc_stop, verbose=1000)
        print('Finished reading SEG-Y headers.')
        
        #  build Python list of sequential trace numbers
        tnum = np.arange(0, buf.num_traces)
        tnum = tnum.tolist()
        print('Finished building sequential trace numbers.')
        
        #  convert trace header dictionary to numpy arrays for convenience
        print('Converting trace header lists to numpy arrays...')
        tic = time()
        il = np.array(thead['il'])
        xl = np.array(thead['xl'])
        toc = time() - tic
        print('Done! (%e seconds)' % toc)
        
        #  get inline and crossline statistics
        print('Getting Inline and Crossline statistics...')
        self.il_min = min(il)
        self.il_max = max(il)
        self.il_num = self.il_max - self.il_min + 1
        
        self.xl_min = min(xl)
        self.xl_max = max(xl)
        self.xl_num = self.xl_max - self.xl_min + 1
        print('Done!')
        
        #  This section of code calculates the trace number in the SEG-Y file 
        #  that each cdp begins at.  The number of traces in that cdp (fold) is
        #  also calculated.
        
        print('Starting the loops for calculating the fold in each CDP')
        
        fold_il = []
        fold_xl = []
        fold_count = []
        cdp_start_trace = []

        
        fold = 1
        start_trace = 0
        il_last = il[0]
        xl_last = xl[0]
        for i in range(1, buf.num_traces):
            
            il_cur = il[i]
            xl_cur = xl[i]
            
            if (il_cur==il_last) & (xl_cur==xl_last):
                fold += 1
                
            else:
                cdp_start_trace.append(start_trace)
                start_trace = i  # reset the starting trace to the current trace number
                
                fold_count.append(fold)
                fold_il.append(il_last)
                fold_xl.append(xl_last)
                fold = 1
                
            il_last = il_cur
            xl_last = xl_cur
        
        fold_count = np.array(fold_count)
        fold_il = np.array(fold_il)
        fold_xl = np.array(fold_xl)
        cdp_start_trace = np.array(cdp_start_trace)
        
        fold_il2 = []
        fold_xl2 = []
        fold_count_2 = []
        cdp_start_trace_2 = []
        for i in set(fold_il):
            for x in set(fold_xl):
                idx = np.nonzero( (fold_il==i) & (fold_xl==x))[0]
                if len(idx)==0:
                    fold_count_2.append(0)
                    cdp_start_trace_2.append(-1)
                else:
                    fold_count_2.append(fold_count[idx[0]])
                    cdp_start_trace_2.append(cdp_start_trace[idx[0]])
        
        cdp_start_trace = cdp_start_trace_2
        fold_count = fold_count_2
        fold_il = fold_il2
        fold_xl = fold_xl2
        
        print(fold_count)
        
        # Set to True to have detailed info printed to command line
        if True:
            for i in range(len(fold_il)):
                print('IL: %d  XL: %d  SeqTr: %i  Fold: %i' % 
                      (fold_il[i], fold_xl[i], cdp_start_trace[i], fold_count[i]) )
        
        print('Done with CDP fold calculations!')
        
        
        # This section of code writes the version 2 segy index file
        
        print('Writing segy index file...')
        with open(self.idx_file, 'wb') as fd:
            
            fd.seek(0)
            
            print('Writing version string')
            buf = pack('32s', 'AuraSegyIndex_v01'.encode())
            fd.write(buf)
            
            print('Writing segy file path')
            buf = pack('1024s', self.segy_file.encode())
            fd.write(buf)
            
            buf = pack('l', self.il_min)
            fd.write(buf)
            
            buf = pack('l', self.xl_min)
            fd.write(buf)
            
            buf = pack('l', self.il_num)
            fd.write(buf)
            
            buf = pack('l', self.xl_num)
            fd.write(buf)
            
            counter = 0
            for i in range( len(cdp_start_trace) ):
                counter += 1
                if counter == 101:
                    print('Writing cdp %i of %i' %(i, len(cdp_start_trace)))
                    counter = 1
                   
                
                buf = pack('ll', cdp_start_trace[i], fold_count[i])
                fd.write(buf)
                
        print('... DONE!')

        
        
    def get_index_tnums(self, il, xl):
        """
        Method to return the starting trace and fold of a particular
        inline/xline bin in 3D survey.  Operates both on post stack and
        prestack SEG-Y volumes.  Assumes that data are sorted first by INLINE, 
        next by XLINE, and lastly by OFFSET but may work in other cases
        
        trc, fold = get_idx_tnums(il, xl)
        
        Inputs:
            il - Inline number in segy file (integer value)
            xl - Xline number in segy file (integer value)
        
        Outputs:
            trc - Trace number (zero indexed) of the first occurrence of the 
                  specified Inline & Xline inputs
           fold - number of occrrences of the Inline & Xline inputs which
                  occur sequentially. For stack volumes, fold will either be
                  1 or -1 (i.e. unpopulated inline/xline bin).  For CMP gathers
                  fold will be either -1 (unpopulated inline/xline bin) or an 
                  integer value corresponding to number of traces in the CMP
                  ensemble.
        """
            
        il_max = self.il_min + self.il_num
        xl_max = self.il_min + self.il_num
        
        if (il>=self.il_min) & (il<=il_max) & (xl>=self.xl_min) & (xl<=xl_max):
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
            
        else:
            cdp_start_trace = -1
            cdp_fold = 0
        
        return cdp_start_trace, cdp_fold

