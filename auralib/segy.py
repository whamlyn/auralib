"""
AuraQI module for reading data stored in SEG-Y format.

Author:   Wes Hamlyn
Created:   1-Sep-2014
Last Mod:  9-Mar-2017

Copyright 2017 Wes Hamlyn

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

from struct import unpack, pack
from numpy import array, arange, flipud, round
import numpy as np
import os
import matplotlib.pyplot as plt
import sys



#  below are example of the Python dictionary structure used to define binary
#  and trace header keywords for SEG-Y format data files.  These can be 
#  copies and editied as necessary for individual segy files.  Note that the
#  binary header dictionary must have 'samp_rate', 'num_samp', and 'samp_fmt'
#  entries otherwise the Segy class will not intialize correctly

def_bhead = {'samp_rate':{'bpos':17, 'fmt':'h', 'nbyte':2},
             'num_samp':{'bpos':21, 'fmt':'h', 'nbyte':2},
             'samp_fmt':{'bpos':25, 'fmt':'h', 'nbyte':2}
             }

def_thead = {'il':{'bpos':189,  'fmt':'l', 'nbyte':4},
             'xl':{'bpos':193, 'fmt':'l', 'nbyte':4},
             'cmpx':{'bpos':181, 'fmt':'l', 'nbyte':4},
             'cmpy':{'bpos':185, 'fmt':'l', 'nbyte':4},
             'offset':{'bpos':37, 'fmt':'l', 'nbyte':4},
             'idcode':{'bpos':29, 'fmt':'h', 'nbyte':2}
             }


class Segy(object):
    """
    Python class for SEG-Y file manipulation
    """

    def __init__(self, filename, def_thead=def_thead, def_bhead=def_bhead):
        """
        Constructor for SEG-Y class.
        """

        self.filename = filename
        self.def_bhead = def_bhead
        self.def_thead = def_thead
        self.bhead = {}
        self.thead = {}
        self.ebcdic = ' '
        self.idx = {}

        for key in def_bhead.keys():
            self.bhead.update({key:[]})

        for key in def_thead.keys():
            self.thead.update({key:[]})

        self._calc_endian()
        self._get_ebcdic_head()
        self._get_bhead()
        self._calc_numbytes()

        self.filesize = os.path.getsize(filename)
        self.trace_size = self.bhead['num_samp']*self.numbytes + 240
        self.num_traces = (self.filesize-3600)//self.trace_size


    def _get_ebcdic_head(self):
        """
        Read EBCDIC header and convert to ASCII.
        """
        
        import codecs
        
        fd = open(self.filename, 'rb')
        
        buf = fd.read(3200)

        fd.close()

        ebhead = ''
        for i in range(0, 3120, 80):
            txt = '%s\n' % (codecs.decode(buf[i:(i+80)], 'cp500'))
            ebhead = ebhead + txt
            
        self.ebcdic = ebhead
        
        
    def _get_bhead(self):
        """
        Read binary header fields.
        """                
            
        fd = open(self.filename, 'rb')

        for key in self.bhead.keys():
            bpos = 3199 + self.def_bhead[key]['bpos']
            fd.seek(bpos, 0)
            buf = fd.read(self.def_bhead[key]['nbyte'])
            buf = unpack(self.def_bhead[key]['fmt'], buf)
            self.bhead[key] = buf[0]

        fd.close()


    def _calc_numbytes(self):
        """
        Calculate the number of bytes for each sample for a given sample format.
        """

        fmt_code = self.bhead['samp_fmt']

        if (fmt_code == 1) | (fmt_code == 6):
            self.fmt_str = 'ibm'
            self.numbytes = 4

        elif fmt_code == 2:
            self.fmt_str = 'l'
            self.numbytes = 4

        elif fmt_code == 3:
            self.fmt_str = 'h'
            self.numbytes = 2

        elif fmt_code == 5:
            self.fmt_str = 'f'
            self.numbytes = 4

        elif fmt_code == 8:
            self.fmt_str = 's'
            self.numbytes = 1

        if self.endian == 'big':
            self.fmt_str = '>' + self.fmt_str

        elif self.endian == 'little':
            self.fmt_str = '<' + self.fmt_str


    def _calc_endian(self):
        """
        Determine if SEG-Y file is big or little endian byte ordered. Adjust the
        appropriate format strings accordingly in the self.def_bhead and
        self.def_thead attributes.
        """

        fd = open(self.filename, 'rb')
        fd.seek(3224, 0)
        buf = fd.read(2)
        buf_big = unpack('>h', buf)[0]
        buf_little = unpack('<h', buf)[0]

        if (buf_big >= 1) & (buf_big <= 8):
            self.endian = 'big'

        elif (buf_little >= 1) & (buf_little <= 8):
            self.endian = 'little'
        
        # set endian format strings
        for key in self.def_bhead.keys():
            if (self.def_bhead[key]['fmt'][0] == '>') | (self.def_bhead[key]['fmt'][0] == '<'):
                next
            else:
                if self.endian == 'little':
                    self.def_bhead[key]['fmt'] = '<' + self.def_bhead[key]['fmt']

                if self.endian == 'big':
                    self.def_bhead[key]['fmt'] = '>' + self.def_bhead[key]['fmt']

        for key in self.def_thead.keys():
            if (self.def_thead[key]['fmt'][0] == '>') | (self.def_thead[key]['fmt'][0] == '<'):
                next
            else:
                if self.endian == 'little':
                    self.def_thead[key]['fmt'] = '<' + self.def_thead[key]['fmt']

                if self.endian == 'big':
                    self.def_thead[key]['fmt'] = '>' + self.def_thead[key]['fmt']        
        

    def get_ilxl(self, il0, xl0, multi=-1, verbose=False):
        """
        Find an inline and crossline number in a segy file using a simple
        binary search.
        """
        
        # create a composite inline-crossline number by scaling the inline 
        # number by a multiple of 10 and then adding the crossline number
        mult10 = 10000000
        ilxl0 = il0*mult10 + xl0
        
        # make an initial guess at a trace
        tmin = 0
        tmax = self.num_traces
        tg = int((tmax+tmin)/2)
        
        count = 0
        ilxlg = -1
        while (ilxl0 != ilxlg) & (count<100):
            count += 1
            #print('iteration %i' % count)
            
            self.read_thead(tg)
            ilg = self.thead['il'][0]
            xlg = self.thead['xl'][0]
            ilxlg = ilg*mult10 + xlg
            
            if verbose:
                print('ilxl0: %i\tilxlg: %i tg: %i' % (ilxl0, ilxlg, tg))
            
            if ilxlg == ilxl0:
                if verbose:
                    print('---Success: IL: %i XL: %i Trace: %i' % (il0, xl0, tg))
                
            elif ilxlg > ilxl0:
                tmax = tg*1
                tg = int(round((tmax + tmin)/2))
                
            elif ilxlg < ilxl0:
                tmin = tg*1
                tg = int(round((tmax + tmin)/2))
        
        if multi>0:
            # Search a number of traces around the "tg" trace to return
            # the trace numbers of the multi-traces
            if verbose:
                print('Using multi search...')
            
            tnum_start = tg-multi
            tnum_end = tg+multi
            tnums = np.arange(tnum_start, tnum_end, 1)
            
            self.read_thead_multi(tnum_start, tnum_end)
            ilg = np.array(self.thead['il'], dtype='int64')
            xlg = np.array(self.thead['xl'], dtype='int64')
            ilxlg = ilg*mult10 + xlg
            
            idx = np.nonzero(ilxl0==ilxlg)[0]
            tg = tnums[idx]
            tg = tg.tolist()
            
        return tg

    
    def read_tdata(self, tracenum):
        """
        Read a single trace from a SEG-Y file.
        Note:
            tracenum is the number of the trace in the SEG-Y file to be read.
            This starts at zero (e.g. first SEG-Y trace has tracenum = 0)
        """

        bpos = 3200 + 400 + 240 + self.trace_size * tracenum

        fd = open(self.filename, 'rb')
        fd.seek(bpos, 0)
    
        
        if self.fmt_str == '>ibm':
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                
                # move cursor to start of trace samples that you want to read
                bpos = 3600 + 240 + self.trace_size*tracenum
                fd.seek(bpos)
                
                # read all the samples for the trace
                ibm_float_trace = fd.read(self.bhead['num_samp'] * self.numbytes)
                tdata = self._ibm2ieee_b_new(ibm_float_trace)
                
                    
        elif self.fmt_str == '<ibm':
            # open file for binary read    
            with open(self.filename, 'rb') as fd:

                # move cursor to start of trace samples that you want to read
                bpos = 3600 + 240 + self.trace_size*tracenum
                fd.seek(bpos)
                
                # read all the samples for the trace
                ibm_float_trace = fd.read(self.bhead['num_samp'] * self.numbytes)
                tdata = self._ibm2ieee_l_new(ibm_float_trace)
                
                    
        else:
            #  build format string
            if self.endian == 'big':
                tr_fmt_str = '>'
            else:
                tr_fmt_str = '<'
            
            for i in range(0, self.bhead['num_samp']):
                tr_fmt_str = tr_fmt_str + self.fmt_str[1:]
            
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                # move cursor to start of trace samples that you want to read
                bpos = 3600 + 240 + self.trace_size*tracenum
                fd.seek(bpos)
                
                # read all the samples for the trace
                buf = fd.read(self.bhead['num_samp'] * self.numbytes)
                tdata = unpack(tr_fmt_str, buf)
    
        return tdata
    
    
    def read_tdata_multi(self, tr_start, tr_end, verbose=0):
        """
        Read multiple sequential traces from a SEG-Y file.  This method is faster
    		than read_trace_data() but is restricted to reading sequential traces.
        Note:
            - this is a development function that attempts to speed up the
              unpacking step
            - tr_start is the first trace to be read (zero-indexed)
            - tr_end is the final trace to be read (inclusive)
            - This method starts at zero (e.g. first SEG-Y trace has 
              tracenum = 0)
        """            
        
        tdata = []
        count = 0
        if self.fmt_str == '>ibm':
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                
                # create an empty python list to store trace amplitudes and begin
                # looping over the traces to be read.
                for tracenum in range(tr_start, tr_end):
                    count += 1
                    if (verbose > 0) & (count == verbose+1):
                        print('reading trace %i' % (tracenum))
                        count = 1
                    
                    # move cursor to start of trace samples that you want to read
                    bpos = 3600 + 240 + self.trace_size*tracenum
                    fd.seek(bpos)
                    
                    # read all the samples for the trace
                    ibm_float_trace = fd.read(self.bhead['num_samp'] * self.numbytes)
                    buf1 = self._ibm2ieee_b_new(ibm_float_trace)
                    
                    tdata.append(buf1)
                    
        elif self.fmt_str == '<ibm':
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                
                # create an empty python list to store trace amplitudes and begin
                # looping over the traces to be read.
                for tracenum in range(tr_start, tr_end):
                    count += 1
                    if (verbose > 0) & (count == verbose+1):
                        print('reading trace %i' % (tracenum))
                        count = 1
                    
                    # move cursor to start of trace samples that you want to read
                    bpos = 3600 + 240 + self.trace_size*tracenum
                    fd.seek(bpos)
                    
                    # read all the samples for the trace
                    ibm_float_trace = fd.read(self.bhead['num_samp'] * self.numbytes)
                    buf1 = self._ibm2ieee_l_new(ibm_float_trace)
                    
                    tdata.append(buf1)
                    
        else:
            #  build format string
            if self.endian == 'big':
                tr_fmt_str = '>'
            else:
                tr_fmt_str = '<'
            
            #for i in range(0, self.bhead['num_samp']):
            #    tr_fmt_str = tr_fmt_str + self.fmt_str[1:]
            tr_fmt_str = '%s%i%s' % (self.fmt_str[0], self.bhead['num_samp'], self.fmt_str[1])
            
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                
                # create an empty python list to store trace amplitudes and begin
                # looping over the traces to be read.
                for tracenum in range(tr_start, tr_end):
                    
                    count += 1
                    if (verbose > 0) & (count == verbose+1):
                        print('reading trace %i' % (tracenum))
                        count = 1                    
    
                    # move cursor to start of trace samples that you want to read
                    bpos = 3600 + 240 + self.trace_size*tracenum
                    fd.seek(bpos)
                    
                    # read all the samples for the trace
                    buf = fd.read(self.bhead['num_samp'] * self.numbytes)
                    buf1 = unpack(tr_fmt_str, buf)
                    
                    tdata.append(buf1)    
        
        return tdata
    

    def read_tdata_random(self, traces, tmin=-1, tmax=-1, verbose=0):
        """
        Read multiple random traces from a SEG-Y file.  This method is faster
    		than read_trace_data().
        Note:
            - tr_start is the first trace to be read (zero-indexed)
            - tr_end is the final trace to be read (inclusive)
            - This method starts at zero (e.g. first SEG-Y trace has 
              tracenum = 0)
        """
        
        if tmin != -1:
            samp_min = tmin / (self.bhead['samp_rate']*0.000001)
            samp_min = int(np.round(samp_min))
        else:
            samp_min = 0
        
        if tmax != -1:
            samp_max = tmax / (self.bhead['samp_rate']*0.000001)
            samp_max = int(np.round(samp_max))
        else:
            samp_max = self.bhead['num_samp']
        
        if (tmin!=-1) | (tmax!=-1):
            nsamp = samp_max - samp_min + 1
        else:
            nsamp = self.bhead['num_samp']
        
        # recast the array of trace numbers to int64, this seems to fix a
        # problem when calculating starting byte locations where the start byte
        # is greater than 2**32
        traces = np.array(traces, dtype='int64')
        
        tdata = []
        count = 0
        if self.fmt_str == '>ibm':
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                
                # create an empty python list to store trace amplitudes and begin
                # looping over the traces to be read.
                for tracenum in traces:
                    count += 1
                    if (verbose > 0) & (count == verbose+1):
                        print('reading trace %i' % (tracenum))
                        count = 1
                    
                    # move cursor to start of trace samples that you want to read
                    bpos = 3600 + 240 + self.trace_size*tracenum + samp_min*self.numbytes
                    fd.seek(bpos)
                    
                    # read all the samples for the trace
                    ibm_float_trace = fd.read(nsamp*self.numbytes)
                    buf1 = self._ibm2ieee_b_new(ibm_float_trace, nsamp)
                    
                    tdata.append(buf1)
                    
        elif self.fmt_str == '<ibm':
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                
                # create an empty python list to store trace amplitudes and begin
                # looping over the traces to be read.
                for tracenum in traces:
                    count += 1
                    if (verbose > 0) & (count == verbose+1):
                        print('reading trace %i' % (tracenum))
                        count = 1
                    
                    # move cursor to start of trace samples that you want to read
                    bpos = 3600 + 240 + self.trace_size*tracenum + samp_min*self.numbytes
                    fd.seek(bpos)
                    
                    # read all the samples for the trace
                    ibm_float_trace = fd.read(nsamp*self.numbytes)
                    buf1 = self._ibm2ieee_l_new(ibm_float_trace, nsamp)
                    
                    tdata.append(buf1)
                    
        else:
            #  build format string
            if self.endian == 'big':
                tr_fmt_str = '>'
            else:
                tr_fmt_str = '<'
            
            #for i in range(0, self.bhead['num_samp']):
            #    tr_fmt_str = tr_fmt_str + self.fmt_str[1:]
            tr_fmt_str = '%s%i%s' % (self.fmt_str[0], nsamp, self.fmt_str[1])
            
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                
                # create an empty python list to store trace amplitudes and begin
                # looping over the traces to be read.
                for tracenum in traces:
                    
                    count += 1
                    if (verbose > 0) & (count == verbose+1):
                        print('reading trace %i' % (tracenum))
                        count = 1                    
    
                    # move cursor to start of trace samples that you want to read
                    bpos = 3600 + 240 + self.trace_size*tracenum + samp_min*self.numbytes
                    fd.seek(bpos)
                    
                    # read all the samples for the trace
                    buf = fd.read(nsamp*self.numbytes)
                    buf1 = unpack(tr_fmt_str, buf)
                    
                    tdata.append(buf1)    
        
        return tdata
       
    
    def read_thead(self, tracenum):
        """
        Read a single trace header from a SEG-Y file.
        Note:
            tracenum is the number of the trace in the SEG-Y file to be read.
            Numbering starts at zero (e.g. first SEG-Y trace has tracenum = 0)
        """

        #   Make sure the thead attribute is empty
        self.thead = {}
        for key in self.def_thead.keys():
            self.thead.update({key:[]})

        fd = open(self.filename, 'rb')

        for key in self.def_thead.keys():

            if self.def_thead[key]['fmt'] == '>ibm':
                bpos = 3599 + tracenum*self.trace_size + self.def_thead[key]['bpos']
                fd.seek(bpos, 0)
                buf = fd.read(self.def_thead[key]['nbyte'])
                buf = self._ibm2ieee_b(buf)
                self.thead[key].append(buf)

            elif self.def_thead[key]['fmt'] == '<ibm':
                bpos = 3599 + tracenum*self.trace_size + self.def_thead[key]['bpos']
                fd.seek(bpos, 0)
                buf = fd.read(self.def_thead[key]['nbyte'])
                buf = self._ibm2ieee_l(buf)
                self.thead[key].append(buf)

            else:
                bpos = 3599 + tracenum*self.trace_size + self.def_thead[key]['bpos']
                fd.seek(bpos, 0)
                buf = unpack(self.def_thead[key]['fmt'], fd.read(self.def_thead[key]['nbyte']))
                self.thead[key].append(buf[0])


        fd.close()

        return self.thead


    def read_thead_multi(self,  tracenum1, tracenum2, verbose=0):
        """
        Read trace headers from sequential traces in a SEG-Y file.
        Note:
            tracenum1 is the number of the first trace to be read.
            tracenum2 is the number of the last trace to be read.
            Numbering starts at zero (e.g. first SEG-Y trace has tracenum = 0)
            verbose is the trace increment to write info to the command line
        """

        #   Make sure the thead attribute is empty
        self.thead = {}
        for key in self.def_thead.keys():
            self.thead.update({key:[]})

        fd = open(self.filename, 'rb')
        count = 0
        for i in range(tracenum1, tracenum2):
            
            count = count + 1
            if (verbose > 0) & (count == verbose+1):
                print('reading trace %i of %i' % (i, tracenum2))
                count = 1
            
            for key in self.def_thead.keys():

                if self.def_thead[key]['fmt'] == '>ibm':
                    bpos = 3599 + i*self.trace_size + self.def_thead[key]['bpos']
                    fd.seek(bpos, 0)
                    buf = fd.read(self.def_thead[key]['nbyte'])
                    buf = self._ibm2ieee_b(buf)
                    self.thead[key].append(buf)

                elif self.def_thead[key]['fmt'] == '<ibm':
                    bpos = 3599 + i*self.trace_size + self.def_thead[key]['bpos']
                    fd.seek(bpos, 0)
                    buf = fd.read(self.def_thead[key]['nbyte'])
                    buf = self._ibm2ieee_l(buf)
                    self.thead[key].append(buf)

                else:
                    bpos = 3599 + i*self.trace_size + self.def_thead[key]['bpos']
                    fd.seek(bpos, 0)
                    buf = unpack(self.def_thead[key]['fmt'], fd.read(self.def_thead[key]['nbyte']))
                    self.thead[key].append(buf[0])

        fd.close()

        return self.thead


    def read_thead_multi_devtest(self,  tracenum1, tracenum2, verbose=0):
        """
        Read trace headers from sequential traces in a SEG-Y file.
        Note:
            tracenum1 is the number of the first trace to be read.
            tracenum2 is the number of the last trace to be read.
            Numbering starts at zero (e.g. first SEG-Y trace has tracenum = 0)
            verbose is the trace increment to write info to the command line
        """

        #   Make sure the thead attribute is empty
        self.thead = {}
        for key in self.def_thead.keys():
            self.thead.update({key:[]})

        fd = open(self.filename, 'rb')
        
        count = 0
        for i in range(tracenum1, tracenum2):
            
            count = count + 1
            if (verbose > 0) & (count == verbose+1):
                print('reading trace %i of %i' % (i, tracenum2))
                count = 1
            
            for key in self.def_thead.keys():
                bpos = 3599 + i*self.trace_size + self.def_thead[key]['bpos']
                fd.seek(bpos, 0)
                buf = fd.read(self.def_thead[key]['nbyte'])
                self.thead[key].append(buf[0])

        fd.close()
                
        for key in self.def_thead.keys():
            if self.def_thead[key]['fmt'] == '>ibm':
                buf = self.thead[key]
                buf = self._ibm2ieee_b_new(buf)
                self.thead[key] = buf
            
            elif self.def_thead[key]['fmt'] == '<ibm':
                buf = self.thead[key]
                buf = self._ibm2ieee_l(buf)
                self.thead[key].append(buf)

            else:
                buf = self.thead[key]
                buf = b''.join(buf)
                ntrc = len(buf)
                fmt_str = '%s%i%s' % (self.def_thead[key]['fmt'][0], ntrc, self.def_thead[key]['fmt'][1:])
                #print(ntrc)
                #print(fmt_str)
                buf = unpack(fmt_str,buf)
                self.thead[key] = buf

        return self.thead

        
        
    def read_thead_random(self,  traces, verbose=0):
        """
        Read trace headers from random traces in a SEG-Y file.
        Note:
            traces is a list of the trace numbers to be read.
            Numbering starts at zero (e.g. first SEG-Y trace has tracenum = 0)
            verbose is the trace increment to write info to the command line
        """
        
        # recast the array of trace numbers to int64, this seems to fix a
        # problem when calculating starting byte locations where the start byte
        # is greater than 2**32
        traces = np.array(traces, dtype='int64')
        
        #   Make sure the thead attribute is empty
        self.thead = {}
        for key in self.def_thead.keys():
            self.thead.update({key:[]})

        fd = open(self.filename, 'rb')
        count = 0
        master_count = 0
        num_traces = len(traces)
        for i in traces:
            master_count += 1
            count += 1
            if (verbose > 0) & (count == verbose+1):
                print('reading trace %i of %i' % (i, num_traces))
                count = 1
            
            for key in self.def_thead.keys():

                if self.def_thead[key]['fmt'] == '>ibm':
                    bpos = 3599 + i*self.trace_size + self.def_thead[key]['bpos']
                    fd.seek(bpos, 0)
                    buf = fd.read(self.def_thead[key]['nbyte'])
                    buf = self._ibm2ieee_b(buf)
                    self.thead[key].append(buf)

                elif self.def_thead[key]['fmt'] == '<ibm':
                    bpos = 3599 + i*self.trace_size + self.def_thead[key]['bpos']
                    fd.seek(bpos, 0)
                    buf = fd.read(self.def_thead[key]['nbyte'])
                    buf = self._ibm2ieee_l(buf)
                    self.thead[key].append(buf)

                else:
                    bpos = 3599 + i*self.trace_size + self.def_thead[key]['bpos']
                    fd.seek(bpos, 0)
                    buf = unpack(self.def_thead[key]['fmt'], fd.read(self.def_thead[key]['nbyte']))
                    self.thead[key].append(buf[0])

        fd.close()

        return self.thead

        
    def write_ebcdic(self, ebcdic_text):
        """
        Writes a blank EBCDIC header
        """
        fd = open(self.filename, 'rb+')
        for i in range(0, 3200):
            buf = pack('c', ' ')
            fd.write(buf)
        fd.close()
    
    
    def write_bhead(self, def_bhead, bhead):
        """
        Writes a Binary header

        def_bhead = {'samp_rate':{'bpos':17, 'fmt':'>h', 'nbyte':2},
                'num_samp':{'bpos':21, 'fmt':'>h', 'nbyte':2},
                'samp_fmt':{'bpos':25, 'fmt':'>h', 'nbyte':2}
                }

        bhead = {'samp_rate':1000,
                'num_samp':1500,
                'samp_fmt':4
                }
        """
        fd = open(self.filename, 'rb+')

        for key in self.def_bhead.keys():

            bpos = self.def_bhead[key]['bpos'] +3199
            print(bpos)
            fd.seek(bpos, 0)
            buf = pack(self.def_bhead[key]['fmt'], self.bhead[key])
            fd.write(buf)

        fd.close()
    
    
    def write_thead(self, tracenum, bpos, fmt, data):
        """
        Writes a Trace header

        tracenum = trace number in file (zero indexed)
        bpos = byte position in trace header to be written
        fmt = f, l, h, c and > or <
        data = single value (int, float, double)
        """
        
        #  make sure the endian character is set in the format string
        if fmt[0] not in  ['>', '<']:
            if self.endian == 'big':
                fmt = '>' + fmt
            elif self.endian == 'little':
                fmt = '<' + 'fmt'
        
        with open(self.filename, 'rb+') as fd:
            abs_bpos = 3600 + bpos-1 + tracenum*self.trace_size
            fd.seek(abs_bpos, 0)
            buf = pack(fmt, data)
            fd.write(buf)
            

    def write_thead_multi(self, tracenums, bpos, fmt, data, verbose=-1):
        """
        Writes a Trace header

        tracenums = list or array of trace numbers (zero indexed)
        bpos = byte position in trace header to be written
        fmt = f, l, h, c and > or <
        data = list or array of header data corresponding to the traces in
               the tracenum list (int, float, double)
        """
        
        #  make sure the endian character is set in the format string
        if fmt[0] not in  ['>', '<']:
            if self.endian == 'big':
                fmt = '>' + fmt
            elif self.endian == 'little':
                fmt = '<' + 'fmt'
        
        #  pack the list/array of data into bytes
        fmt = fmt[0] + str(len(data)) + fmt[1:]
        bbuf = pack(fmt, *data)
        nbytes = int(len(bbuf)/len(data))
        
        with open(self.filename, 'rb+') as fd:
            
            count = 0
            ntraces = len(tracenums)
            for i in range(ntraces):
                
                count += 1
                if (verbose > 0) & (count == verbose+1):
                    print('Writing trace header %i of %i' % (tracenums[i], tracenums[-1]))
                    count = 1
                    
                abs_bpos = 3599 + tracenums[i]*self.trace_size + bpos
                fd.seek(abs_bpos, 0)
                b1 = nbytes*i
                b2 = b1+nbytes
                fd.write(bbuf[b1:b2])
        
    
    def write_trace_data(self, tracenum, tdata):
        """
        Writes trace data

        tracenum = trace number in file (zero indexed)
        data = single value (int, float, double)
        """

        with open(self.filename, 'rb+') as fd:
            
            bpos = 3840 + tracenum*self.trace_size
            fd.seek(bpos, 0)
            
            nsamp = self.bhead['num_samp']
            fmt_str = '%s%i%s' % (self.fmt_str[0], nsamp, self.fmt_str[1])
            
            buf = pack(fmt_str, *tdata)
            fd.write(buf)
    
    
    def write_trace_data_multi(self, tracenums, tdata):
        """
        Writes trace data

        tracenums = array of trace numbers in file corresponding to the tdata
                    trace array (zero indexed)
        tdata = 2d array (or list) containing traces to be written where axis 1
                is the trace number and axis 2 is the sample number.
        """
        
        nsamp = self.bhead['num_samp']
        fmt_str = '%s%i%s' % (self.fmt_str[0], nsamp, self.fmt_str[1])
        
        with open(self.filename, 'rb+') as fd:
            
            for tracenum, data in zip(tracenums, tdata):
                
                bpos = 3840 + tracenum*self.trace_size
                fd.seek(bpos, 0)
                
                buf = pack(fmt_str, *data)
                fd.write(buf)
    
    
    def _ibm2ieee_b(self, ibm_float):
        """
        Convert IBM Float (big endian byte order)
        Old method, works only on single floating point word.  New method will
        operate faster on lists of floating point words.
        """

        dividend = float(16**6)
        
        if ibm_float == 0:
            return 0.0
        
        istic, a, b, c = unpack('>BBBB', ibm_float)
        
        if istic >= 128:
            sign= -1.0
            istic = istic - 128
        else:
            sign = 1.0
            
        mant = float(a<<16) + float(b<<8) + float(c)

        return sign*16**(istic - 64)*(mant/dividend)
    
    
    def _ibm2ieee_l(self, ibm_float):
        """
        Convert IBM float (little endian byte order)
        Old method, works only on single floating point word.  New method will
        operate faster on lists of floating point words.
        """

        dividend = float(16**6)
        
        if ibm_float == 0:
            return 0.0
        c, b, a, istic = unpack('>BBBB', ibm_float)
        
        if istic >= 128:
            sign = -1.0
            istic = istic - 128
        else:
            sign = 1.0
        
        mant = float(a<<16) + float(b<<8) + float(c)
        
        return sign*16**(istic - 64)*(mant/dividend)
    
    
    def _ibm2ieee_b_new(self, ibm_float_trace, nsamp=-1):
        """
        Convert IBM Float (big endian byte order)
        New method.  More efficient when operating on lists of floating point
        words.
        """
        
        #  build trace format string
        
        #  The three commented lines below are an older way of building the 
        #  format string using a for loop.  Instead the format string will be
        #  built using the format code B preceeded by an integer equal to 4x
        #  the number of samples.
        
        #fmt_str = '>'
        #for i in range(0, self.bhead['num_samp']):
        #    fmt_str = fmt_str + 'BBBB'
        
        if nsamp == -1:
            nsamp = self.bhead['num_samp']
        
        fmt_str = '>%iB' % (nsamp*4)
        
        buf = unpack(fmt_str, ibm_float_trace)
        
        buf = np.array(buf, dtype='int32')
        
        istic = buf[0::4]
        a = buf[1::4]
        b = buf[2::4]
        c = buf[3::4]
        
        sign = np.ones(len(istic), dtype='int32')
        idx = np.nonzero(istic >= 128)
        sign[idx] = -1
        istic[idx] = istic[idx] - 128
        
        dividend = float(16**6)
        # the following two lines are equivalent results, one using
        # bit shifting, the other using multiplication
        #mant = a*2.0**16.0 + b*2.0**8.0 + c
        mant = np.left_shift(a, 16) + np.left_shift(b, 8) + c
        buf = sign * (16.0**(istic-64))*(mant/dividend)
        
        return buf.tolist()
    
    
    def _ibm2ieee_l_new(self, ibm_float_trace, nsamp=-1):
        """
        Convert IBM Float (little endian byte order)
        New method.  More efficient when operating on lists of floating point
        words.
        """
        
        #  build trace format string
        
        #  The three commented lines below are an older way of building the 
        #  format string using a for loop.  Instead the format string will be
        #  built using the format code B preceeded by an integer equal to 4x
        #  the number of samples.
        
        #fmt_str = '<'
        #for i in range(0, self.bhead['num_samp']):
        #    fmt_str = fmt_str + 'BBBB'
        
        if nsamp == -1:
            nsamp = self.bhead['num_samp']
        
        fmt_str = '<%iB' % (nsamp*4)
        buf = unpack(fmt_str, ibm_float_trace)
        
        buf = np.array(buf, dtype='int32')
        
        istic = buf[3::4]
        a = buf[2::4]
        b = buf[1::4]
        c = buf[0::4]
        
        sign = np.ones(len(istic), dtype='int32')
        idx = np.nonzero(istic >= 128)
        sign[idx] = -1
        istic[idx] = istic[idx] - 128
        
        dividend = float(16**6)
        # the following two lines are equivalent results, one using
        # bit shifting, the other using multiplication
        #mant = a*2.0**16.0 + b*2.0**8.0 + c
        mant = np.left_shift(a, 16) + np.left_shift(b, 8) + c
        buf = sign * (16.0**(istic-64))*(mant/dividend)
        
        return buf.tolist()



def write_blank_segy_file(filename, fmt_code, samp_rate, num_samp, num_trace, 
                          endian='big', verbose=500):
    """
    Writes a blank SEG-Y file with empty trace headers and zero sample 
    amplitudes.  Minimal binary header values (sample format, number of 
    samples, and sample rate) will be populated.
    
    filename = string containing full filename and path of output segy file
    fmt_code = sample encoding for trace amplitude samples (1=ibm float, 
               2=32-bit int, 3=16-bit int, 5=ieee float, 8=8-bit int)
    samp_rate = sample rate in microseconds (i.e. 2000 us = 2 ms)
    num_samp = number of samples per trace
    num_trace = number of traces in the segy file
    endian = endian format ('big' for UNIX byte order; 'little' for PC byte order)
    """
    
    #   Set appropriate format string for writing data samples
    if fmt_code == 1:
        fmt_str = 'ibm'
        sys.exit('fmt_code = 1: IBM floating point not yet supported for SEG-Y write')
        
    elif fmt_code == 2:
        fmt_str = 'l'

    elif fmt_code == 3:
        fmt_str = 'h'

    elif fmt_code == 5:
        fmt_str = 'f'

    elif fmt_code == 8:
        fmt_str = 's'

    #   Build null EBCDIC header
    ehead = []
    for i in range(0, 3200):
        ehead.append(' ')

    #   Build minimal Binary header
    bhead = []
    for i in range(0, 200):
        bhead.append(0)
    
    bhead[8] = samp_rate
    bhead[10] = num_samp
    bhead[12] = fmt_code

    #   Build null Trace header
    thead = []
    for i in range(0, 120):
        thead.append(0)
    thead[14] = 1  # set dead trace flag 0=dead, 1=live

    #   Build null Trace data
    tdata = []
    for i in range(0, num_samp):
        tdata.append(0)

    #   Start writing to disk
    with open(filename, 'w') as fd:
        for buf in ehead:
            fd.write(buf)
            
    with open(filename, 'ab') as fd:
        
        #   Do appropriate things for big endian
        if endian=='big':
            for buf in bhead:
                fd.write(pack('>h', buf))
    
            fmt_str = '>' + fmt_str
            
            count = 0
            for i in range(0, num_trace):
                count += 1
                if (verbose > 0) & (count == verbose+1):
                    print('writing trace %i of %i' % (i, num_trace))
                    count = 1
                    
                for buf in thead:
                    fd.write(pack('>h', buf))
                    
                for buf in tdata:
                    fd.write(pack(fmt_str, buf))
    
        #   Do appropriate things for little endian
        if endian=='little':
            for buf in bhead:
                fd.write(pack('<h', buf))
    
            for i in range(0, num_trace):
                count += 1
                if (verbose > 0) & (count == verbose+1):
                    print('writing trace %i of %i' % (i, num_trace))
                    count = 1
                    
                for buf in thead:
                    fd.write(pack('<h', buf))
    
                for buf in tdata:
                    fd.write(pack(fmt_str, buf))


def write_blank_segy_file_v2(filename, fmt_code, samp_rate, num_samp, num_trace, 
                          endian='big', verbose=100):
    """
    Writes a blank SEG-Y file with empty trace headers and zero sample 
    amplitudes.  Minimal binary header values (sample format, number of 
    samples, and sample rate) will be populated.
    
    filename = string containing full filename and path of output segy file
    fmt_code = sample encoding for trace amplitude samples (1=ibm float, 
               2=32-bit int, 3=16-bit int, 5=ieee float, 8=8-bit int)
    samp_rate = sample rate in microseconds (i.e. 2000 us = 2 ms)
    num_samp = number of samples per trace
    num_trace = number of traces in the segy file
    endian = endian format ('big' for UNIX byte order; 'little' for PC byte order)
    """
    
    #   Set appropriate format string for writing data samples
    if fmt_code == 1:
        fmt_str = 'ibm'
        sys.exit('fmt_code = 1: IBM floating point not yet supported for SEG-Y write')
        
    elif fmt_code == 2:
        fmt_str = '>l'

    elif fmt_code == 3:
        fmt_str = 'h'

    elif fmt_code == 5:
        fmt_str = 'f'

    elif fmt_code == 8:
        fmt_str = 's'
    
    #   Set appropriate endian format code
    if endian == 'big':
        fmt_str = '>' + fmt_str
    else:
        fmt_str = '<' + fmt_str

    #   Build null EBCDIC header
    ehead = []
    for i in range(0, 3200):
        ehead.append(' ')

    #   Build minimal Binary header
    bhead = []
    for i in range(0, 200):
        bhead.append(0)
    
    bhead[8] = samp_rate
    bhead[10] = num_samp
    bhead[12] = fmt_code

    #   Build null Trace header
    thead = []
    for i in range(0, 120):
        thead.append(0)
    thead[14] = 1  # set dead trace flag; 1=live 2=dead

    #   Build null Trace data
    tdata = []
    for i in range(0, num_samp):
        tdata.append(0)

    #   Start writing to disk
    with open(filename, 'w') as fd:
        for buf in ehead:
            fd.write(buf)
            
    with open(filename, 'ab') as fd:
        
        #   Do appropriate things for big endian
        if endian=='big':
            for buf in bhead:
                fd.write(pack('>h', buf))
            
            fmt_str = '%s%i%s' % (fmt_str[0], len(tdata), fmt_str[1])
            tdatap = pack(fmt_str, *tdata)
            
            count = 0
            for i in range(0, num_trace):
                count += 1
                if (verbose > 0) & (count == verbose+1):
                    print('writing trace %i of %i' % (i, num_trace))
                    count = 1
                    
                for buf in thead:
                    fd.write(pack('>h', buf))
                    
                fd.write(tdatap)
    
        #   Do appropriate things for little endian
        if endian=='little':
            for buf in bhead:
                fd.write(pack('<h', buf))
            
            fmt_str = '%s%i%s' % (fmt_str[0], len(tdata), fmt_str[1])
            tdatap = pack(fmt_str, *tdata)
            
            count = 0
            for i in range(0, num_trace):
                count += 1
                if (verbose > 0) & (count == verbose+1):
                    print('writing trace %i of %i' % (i, num_trace))
                    count = 1
                    
                for buf in thead:
                    fd.write(pack('<h', buf))
                
                fd.write(tdatap)



def plot_seis(ax, tdata, t_min=0, t_max=0, tr_min=0, tr_max=0, samp_rate=0.002,
              cmap=plt.cm.gray_r, amp_min=0, amp_max=0):
    """
    Convenience function for plotting seismic trace data to a maptlotlib axes.
    """
    
    from matplotlib.ticker import FormatStrFormatter
    majorFormatter = FormatStrFormatter('%i')
    
    tdata = array(tdata)
    num_trc, num_samp = tdata.shape
    
    if (amp_min==0) & (amp_max==0):
        amp_min = tdata.min()
        amp_max = tdata.max()
    
    if (tr_min == 0) & (tr_max==0):
        tr_min = 0
        tr_max = num_trc

    if (t_min==0) & (t_max==0):
        t_min = 0
        t_max = num_samp*samp_rate
        
    datalim = [tr_min, tr_max, t_min, t_max]
    
    im = ax.imshow(flipud(tdata.T), extent=datalim,
               cmap=cmap, vmin=amp_min, vmax=amp_max, interpolation='bilinear')
               
    ax.invert_yaxis()
    ax.set_aspect('auto')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_major_formatter(majorFormatter)
    
    return im


def plot_wigva(ax, tdata, t, trcstart=0, excursion=1, peak=False, trough=False, 
          line=True, lw=0.25, lcolor='k', pcolor=[0.2, 0.2, 1.0],
          tcolor=[1.0, 0.2, 0.2], alpha=0.5):
    """
    Plot wiggle traces with variable area fill.
    If using matplotlib earlier than version 2.1.0, use this function.
    """
    
    from scipy.interpolate import interp1d
    
    # estimate the sample rate; do it this way to avoid errors due to floating
    # point precision
    dt = int((t[1]-t[0])*1000)*0.001
    
    # set the new sample rate for the VA fill to 1/5 of the original sampling
    # This is done so that 
    dt2 = dt*0.25
    
    tdata = np.array(tdata)
    ntrc, nsamp = tdata.shape
    
    norm = np.max(np.abs([np.nanmax(tdata), np.nanmin(tdata)]))
    for i in range(0, ntrc):
        
        zeroval = i + trcstart
        trc = tdata[i, :]
        
        if np.isnan(trc).any():
            next
        
        else:
            #norm = max(abs([max(trc), min(trc)]))
            trc = trc/norm * excursion + zeroval
            
            if (peak == True) | (trough == True):
                # interpolate the new trace
                t2 = arange(t.min(), t.max(), dt2)
                f = interp1d(t, trc, kind='linear')
                trc2 = f(t2)
                
                # plot peak fill
                if peak==True:
                    ax.fill_betweenx(t2, zeroval, trc2, where=trc2>zeroval, 
                                     facecolor=pcolor, edgecolor=pcolor,
                                     alpha=alpha)
                
                # plot trough fill
                if trough==True:
                    ax.fill_betweenx(t2, zeroval, trc2, where=trc2<zeroval, 
                                     facecolor=tcolor,  edgecolor=tcolor,
                                     alpha=alpha)
            if line == True:
                ax.plot(trc, t, lcolor, lw=lw)


def plot_wigva2(ax, tdata, t, trcstart=0, excursion=1, peak=False, trough=False, 
          line=True, lw=0.25, lcolor='k', pcolor=[0.2, 0.2, 1.0],
          tcolor=[1.0, 0.2, 0.2], alpha=0.5):
    """
    Plot wiggle traces with variable area fill. 
    If using matplotlib version 2.1.0 and later, use this function.
    """
        
    tdata = np.array(tdata)
    ntrc, nsamp = tdata.shape
    
    norm = np.max(np.abs([np.nanmax(tdata), np.nanmin(tdata)]))
    for i in range(0, ntrc):
        
        zeroval = i + trcstart
        trc = tdata[i, :]
        
        if np.isnan(trc).any():
            next
        
        else:
            #norm = max(abs([max(trc), min(trc)]))
            trc = trc/norm * excursion + zeroval
            
            if (peak == True) | (trough == True):
                                
                # plot peak fill
                if peak==True:
                    ax.fill_betweenx(t, zeroval, trc, where=trc>zeroval, 
                                     facecolor=pcolor, edgecolor=pcolor,
                                     alpha=alpha, interpolate=True)
                
                # plot trough fill
                if trough==True:
                    ax.fill_betweenx(t, zeroval, trc, where=trc<zeroval, 
                                     facecolor=tcolor,  edgecolor=tcolor,
                                     alpha=alpha, interpolate=True)
            if line == True:
                ax.plot(trc, t, lcolor, lw=lw)
        
    
def time_to_samp(twt, dt):
    """
    Convert an array of TWT to an array of sample indicies
    
    samp = time_to_samp(twt, dt)
    """
    
    samp = round(twt / dt)
    samp = array(samp, dtype='int32')
    
    return samp
