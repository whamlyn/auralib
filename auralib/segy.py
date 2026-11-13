"""
auralib module for reading data stored in SEG-Y format.

Author:   Wes Hamlyn
Created:   1-Sep-2014
Last Mod:  9-Mar-2017

"""

from struct import unpack, pack
from numpy import array, arange, flipud, round
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import codecs
from numba import jit, njit


#  below are example of the Python dictionary structure used to define binary
#  and trace header keywords for SEG-Y format data files.  These can be 
#  copied and edited as necessary for individual segy files.  Note that the
#  binary header dictionary must have 'samp_rate', 'num_samp', and 'samp_fmt'
#  entries otherwise the Segy class will not initialize correctly

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

def print_fmtdef():
    """
    Convenience function to print a header format definition to the terminal
    window for easy copying into a script
    """
    
    print("def_thead = {'il':{'bpos':189,  'fmt':'l', 'nbyte':4},")
    print("             'xl':{'bpos':193, 'fmt':'l', 'nbyte':4},")
    print("             'cmpx':{'bpos':181, 'fmt':'l', 'nbyte':4},")
    print("             'cmpy':{'bpos':185, 'fmt':'l', 'nbyte':4},")
    print("             'offset':{'bpos':37, 'fmt':'l', 'nbyte':4},")
    print("             'idcode':{'bpos':29, 'fmt':'h', 'nbyte':2}}")
    

class Segy(object):
    """
    Python class for SEG-Y file manipulation
    """

    def __init__(self, filename, def_thead=def_thead, def_bhead=def_bhead, 
                 ebcdic_fmt='ebcdic'):
        """
        Constructor for SEG-Y class.
        """

        self.filename = filename
        self.def_bhead = def_bhead
        self.def_thead = def_thead
        self.bhead = {}
        self.thead = {}
        self.ebcdic_fmt = ebcdic_fmt
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
        
        if self.ebcdic_fmt == 'ebcdic':
            encoding = 'cp500'
        elif self.ebcdic_fmt == 'ascii':
            encoding = 'utf-8'
        else:
            encoding = 'cp500'
            
        with open(self.filename, 'rb') as fd:
            buf = fd.read(3200)

        ebhead = ''
        for i in range(0, 3181, 80):
            txt = '%s\n' % (codecs.decode(buf[i:(i+80)], encoding))
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
        ilxl0 = np.array(il0, dtype='int64')*mult10 + xl0
        
        # make an initial guess at a trace
        tmin = 0
        tmax = self.num_traces
        tg = int((tmax+tmin)/2)
        tg0 = tg*1
        
        max_iter = 100 # maximum number of iterations before deciding failure
        count = 0
        ilxlg = -1
        while (ilxl0 != ilxlg) & (count<max_iter):
            count += 1
            if verbose:
                print('iteration %i' % count)
            
            self.read_thead(tg)
            ilg = self.thead['il'][0]
            xlg = self.thead['xl'][0]
            ilxlg = ilg*mult10 + xlg
            
            if verbose:
                print('ilxl0: %i\tilxlg: %i tg: %i' % (ilxl0, ilxlg, tg))
            
            if ilxlg == ilxl0:
                if verbose:
                    print('---Success: IL: %i XL: %i Trace: %i' % (il0, xl0, tg))
                tg = [tg]
                
            elif (count == max_iter) | (np.abs(tg-tg0)==1):
                if verbose:
                    print('Failed to find desired inline. Terminating...')
                return []
                #break
                
            elif ilxlg > ilxl0:
                tmax = tg*1
                tg = int(round((tmax + tmin)/2))
                
            elif ilxlg < ilxl0:
                tmin = tg*1
                tg = int(round((tmax + tmin)/2))
            
            tg0 = tg*1
        
        if multi>0:
            # Search a number of traces around the "tg" trace to return
            # the trace numbers of the multi-traces
            if verbose:
                print('Using multi search...')
            
            tnum_start = tg[0] - multi
            tnum_end = tg[0] + multi
            
            if tnum_start < 0:
                tnum_start = 0
            
            if tnum_end > self.num_traces:
                tnum_end = self.num_traces
                
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
    
    
    def read_tdata_multi_v1(self, tr_start, tr_end, skip=1, verbose=0, as_array=True):
        """
        Read multiple sequential traces from a SEG-Y file.  This method is faster
    	than read_trace_data() but is restricted to reading sequential traces.
        
        Note:
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
                for tracenum in range(tr_start, tr_end, skip):
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
                for tracenum in range(tr_start, tr_end, skip):
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
                for tracenum in range(tr_start, tr_end, skip):
                    
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
    
    
    def _trace_numpy_dtype(self):
        """
        Return the NumPy dtype used for non-IBM SEG-Y trace samples.
        """
        fmt_code = self.bhead['samp_fmt']
    
        dtype_codes = {
            2: 'i4',   # four-byte two's-complement integer
            3: 'i2',   # two-byte two's-complement integer
            5: 'f4',   # four-byte IEEE floating point
            8: 'i1',   # one-byte two's-complement integer
        }
    
        if fmt_code not in dtype_codes:
            raise ValueError(
                'unsupported NumPy trace sample format code: %r'
                % fmt_code
            )
    
        if self.endian == 'big':
            endian = '>'
        else:
            endian = '<'
    
        return np.dtype(endian + dtype_codes[fmt_code])


    def read_tdata_multi(self, tr_start, tr_end, skip=1, verbose=0,
                         as_array=True):
        """
        Read sample data from multiple sequential SEG-Y traces.
    
        Trace samples are read through strided NumPy memory-map views. This
        avoids a seek/read/unpack operation for every trace.
    
        Parameters
        ----------
        tr_start : int
            First trace to read, zero indexed.
        tr_end : int
            Exclusive upper trace number, as for ``range``.
        skip : int, optional
            Trace increment.
        verbose : int, optional
            If greater than zero, process and report progress after this many
            output traces. For example, verbose=100000 reports progress every
            100,000 traces.
        as_array : bool, optional
            If True, return a two-dimensional NumPy array with shape
            ``(number_of_traces, number_of_samples)``. If False, return a
            nested Python list.
    
        Returns
        -------
        numpy.ndarray or list
            Trace data organized as trace by sample.
        """
        if skip <= 0:
            raise ValueError('skip must be greater than zero')
    
        if tr_start < 0 or tr_end < tr_start:
            raise ValueError('invalid trace range')
    
        if verbose < 0:
            raise ValueError('verbose must be zero or greater')
    
        # Clamp the requested range to the available traces.
        tr_end = min(tr_end, self.num_traces)
    
        ntrace = len(range(tr_start, tr_end, skip))
        nsamp = self.bhead['num_samp']
    
        if ntrace == 0:
            result = np.empty((0, nsamp), dtype=np.float64)
            return result if as_array else result.tolist()
    
        # With verbose=0, process all requested traces in one operation.
        # Otherwise, verbose is both the chunk size and reporting interval.
        chunk_size = int(verbose) if verbose > 0 else ntrace
    
        raw = np.memmap(
            self.filename,
            dtype=np.uint8,
            mode='r'
        )
    
        chunks = []
    
        try:
            for start in range(0, ntrace, chunk_size):
                stop = min(start + chunk_size, ntrace)
                chunk_ntrace = stop - start
    
                # Convert the output-array position into the corresponding
                # physical trace number in the SEG-Y file.
                chunk_first_trace = tr_start + start * skip
    
                sample_offset = (
                    3600
                    + chunk_first_trace * self.trace_size
                    + 240
                )
    
                trace_stride = self.trace_size * skip
    
                if self.fmt_str in ('>ibm', '<ibm'):
                    # IBM floating-point values have no native NumPy dtype.
                    # First gather the four bytes belonging to each sample.
                    byte_view = np.ndarray(
                        shape=(chunk_ntrace, nsamp, 4),
                        dtype=np.uint8,
                        buffer=raw,
                        offset=sample_offset,
                        strides=(trace_stride, 4, 1)
                    )
    
                    # Make the sample bytes contiguous before conversion.
                    packed = np.ascontiguousarray(byte_view).tobytes()
                    total_samples = chunk_ntrace * nsamp
    
                    if self.fmt_str == '>ibm':
                        values = self._ibm2ieee_b_new(
                            packed,
                            nsamp=total_samples
                        )
                    else:
                        values = self._ibm2ieee_l_new(
                            packed,
                            nsamp=total_samples
                        )
    
                    values = np.asarray(
                        values,
                        dtype=np.float64
                    ).reshape(chunk_ntrace, nsamp)
    
                else:
                    dtype = self._trace_numpy_dtype()
    
                    values = np.ndarray(
                        shape=(chunk_ntrace, nsamp),
                        dtype=dtype,
                        buffer=raw,
                        offset=sample_offset,
                        strides=(trace_stride, dtype.itemsize)
                    )
    
                    # Copy the data because the view is backed by the memmap.
                    values = np.array(values)
    
                chunks.append(values)
    
                if verbose > 0:
                    print(
                        'reading trace data: %i of %i'
                        % (stop, ntrace)
                    )
    
        finally:
            del raw
    
        if len(chunks) == 1:
            result = chunks[0]
        else:
            result = np.concatenate(chunks, axis=0)
    
        return result if as_array else result.tolist()

    
    def read_tdata_multi_dev(self, tr_start, tr_end, skip=1, verbose=0, as_array=True):
        """
        Read multiple sequential traces from a SEG-Y file.  This method is faster
    	than read_trace_data() but is restricted to reading sequential traces.
        
        Note:
            - tr_start is the first trace to be read (zero-indexed)
            - tr_end is the final trace to be read (inclusive)
            - This method starts at zero (e.g. first SEG-Y trace has 
              tracenum = 0)
        """            
        
        
        ntraces = tr_end - tr_start
        nsamp = self.bhead['num_samp']*1
        
        tdata = np.zeros((ntraces, nsamp))
        
        count = 0
        if self.fmt_str == '>ibm':
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                
                # create an empty python list to store trace amplitudes and begin
                # looping over the traces to be read.
                for tracenum in range(tr_start, tr_end, skip):
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
                    
                    tdata[count] = np.asarray(buf1)
                    
        elif self.fmt_str == '<ibm':
            # open file for binary read    
            with open(self.filename, 'rb') as fd:
                
                # create an empty python list to store trace amplitudes and begin
                # looping over the traces to be read.
                for tracenum in range(tr_start, tr_end, skip):
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
                    
                    tdata[count] = np.asarray(buf1)
                    
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
                for tracenum in range(tr_start, tr_end, skip):
                    
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
                    
                    tdata[count] = np.asarray(buf1)
        
        return tdata if as_array else tdata.tolist()

    

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


    def read_thead_multi_v1(self,  tracenum1, tracenum2, skip=1, verbose=0):
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
        for i in range(tracenum1, tracenum2, skip):
            
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


    def read_thead_multi(self, tracenum1, tracenum2, skip=1, verbose=0,
                         as_array=True):
        """
        Read trace headers from sequential traces in a SEG-Y file.
    
        Header values are read through strided NumPy memory-map views. This
        avoids a seek/read pair for every field on every trace.
    
        Parameters
        ----------
        tracenum1 : int
            First trace to read, zero indexed.
        tracenum2 : int
            Exclusive upper trace number, as for ``range``.
        skip : int, optional
            Trace increment.
        verbose : int, optional
            If greater than zero, process trace headers in chunks of this many
            traces and print a progress message after each chunk.
            For example, verbose=100000 reports progress every 100,000 traces.
        as_array : bool, optional
            If True, return header values as NumPy arrays. If False, return
            Python lists.
        """
        if skip <= 0:
            raise ValueError('skip must be greater than zero')
        if tracenum1 < 0 or tracenum2 < tracenum1:
            raise ValueError('invalid trace range')
    
        tracenum2 = min(tracenum2, self.num_traces)
        ntrace = len(range(tracenum1, tracenum2, skip))
    
        self.thead = {key: [] for key in self.def_thead}
    
        if ntrace == 0:
            return self.thead
    
        raw = np.memmap(self.filename, dtype=np.uint8, mode='r')
    
        try:
            # -------------------------------------------------------------
            # Determine chunk size.
            #
            # With verbose=0, read everything in one chunk so performance
            # remains equivalent to the current implementation.
            # -------------------------------------------------------------
            if verbose > 0:
                chunk_size = int(verbose)
            else:
                chunk_size = ntrace
    
            results = {key: [] for key in self.def_thead}
    
            for start in range(0, ntrace, chunk_size):
    
                stop = min(start + chunk_size, ntrace)
                chunk_ntrace = stop - start
    
                # Actual SEG-Y trace corresponding to the first element
                # of this chunk.
                chunk_first_trace = tracenum1 + start * skip
    
                trace_stride = self.trace_size * skip
    
                for key, definition in self.def_thead.items():
    
                    fmt = definition['fmt']
                    nbyte = definition['nbyte']
    
                    offset = (
                        3600
                        + chunk_first_trace * self.trace_size
                        + definition['bpos'] - 1
                    )
    
                    if fmt in ('>ibm', '<ibm'):
                        # IBM floating point has no NumPy dtype.
                        byte_view = np.ndarray(
                            shape=(chunk_ntrace, nbyte),
                            dtype=np.uint8,
                            buffer=raw,
                            offset=offset,
                            strides=(trace_stride, 1)
                        )
    
                        packed = np.ascontiguousarray(byte_view).tobytes()
    
                        if fmt == '>ibm':
                            values = self._ibm2ieee_b_new(
                                packed, nsamp=chunk_ntrace
                            )
                        else:
                            values = self._ibm2ieee_l_new(
                                packed, nsamp=chunk_ntrace
                            )
    
                        values = np.asarray(values)
    
                    else:
                        dtype = self._header_numpy_dtype(fmt, nbyte)
    
                        values = np.ndarray(
                            shape=(chunk_ntrace,),
                            dtype=dtype,
                            buffer=raw,
                            offset=offset,
                            strides=(trace_stride,)
                        )
    
                        # Copy because values is backed by the memmap.
                        values = np.array(values)
    
                    results[key].append(values)
    
                # ---------------------------------------------------------
                # Progress message
                # ---------------------------------------------------------
                if verbose > 0:
                    print(
                        'reading trace headers: %i of %i'
                        % (stop, ntrace)
                    )
    
            # Combine the individual chunks.
            for key in self.def_thead:
                values = np.concatenate(results[key])
    
                if as_array:
                    self.thead[key] = values
                else:
                    self.thead[key] = values.tolist()
    
        finally:
            del raw
    
        return self.thead


    @staticmethod
    def _header_numpy_dtype(fmt, nbyte):
        """Translate a scalar ``struct`` header format to a NumPy dtype."""
        endian = fmt[0] if fmt[0] in ('>', '<', '=', '@') else '='
        code = fmt[-1]
        dtype_codes = {
            'b': 'i1', 'B': 'u1', '?': '?',
            'h': 'i2', 'H': 'u2',
            'i': 'i4', 'I': 'u4',
            'l': 'i4', 'L': 'u4',
            'q': 'i8', 'Q': 'u8',
            'f': 'f4', 'd': 'f8',
        }

        if code == 'c' and nbyte == 1:
            return np.dtype('S1')
        if code not in dtype_codes:
            raise ValueError("unsupported trace-header format: %r" % fmt)

        dtype = np.dtype(endian + dtype_codes[code])
        if dtype.itemsize != nbyte:
            raise ValueError(
                "trace-header format %r requires %i bytes, definition says %i"
                % (fmt, dtype.itemsize, nbyte))
        return dtype

        
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
    
    
    def write_ebcdic_new(self, ebcdic_text_file):
        """
        Write or update an ebcdic header
        
        ebcdic_text = ASCII file containing ebcdic header text, must be 40
                      lines in length
        """
        
        with open(ebcdic_text_file, 'r') as fd:
            txt = fd.readlines()
        
        with open(self.filename, 'rb+') as fd:
            fd.seek(0, 0)
            for line in txt:
                # make sure no line is longer than 80 characters, truncate if
                # they are
                if len(line)>80:
                    line = line[0:80]
                
                # remove tabsm end of line characters, and that lines are all
                # 80 characters long
                line = '%-80s' % (line.replace('\t', '    ').replace('\n', ''))
                
                # make sure blank lines are written out as spaces
                if len(line)==0:
                    line = '%80s' % (80*' ')
                
                # convert from ascii to ebcdic
                line = codecs.encode(line, 'cp500')
                
                # write to segy
                fd.write(line)
        
        # re-read the updated ebcdic header to Segy object ebcdic attribute
        self._get_ebcdic_head()
    
    
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
        self.bhead = bhead
        fd = open(self.filename, 'rb+')
        
        for key in self.def_bhead.keys():

            bpos = self.def_bhead[key]['bpos'] + 3199
            #print(bpos)
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
        
        # recast tracenum to int64 to avoid overflowing long ints
        tracenum = np.array(tracenum, dtype='int64')

        #  make sure the endian character is set in the format string
        if fmt[0] not in  ['>', '<']:
            if self.endian == 'big':
                fmt = '>' + fmt
            elif self.endian == 'little':
                fmt = '<' + fmt
        
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

        # recast tracenums to int64 to avoid overflowing long ints
        tracenums = np.array(tracenums, dtype='int64')
        
        #  make sure the endian character is set in the format string
        if fmt[0] not in  ['>', '<']:
            if self.endian == 'big':
                fmt = '>' + fmt
            elif self.endian == 'little':
                fmt = '<' + fmt
        
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
        
    
    def write_tdata(self, tracenum, tdata):
        """
        Writes trace data

        tracenum = trace number in file (zero indexed)
        data = single value (int, float, double)
        """
        
        # recast tracenums to int64 to avoid overflowing long ints
        tracenum = np.array(tracenum, dtype='int64')
        
        with open(self.filename, 'rb+') as fd:
            
            bpos = 3840 + tracenum*self.trace_size
            fd.seek(bpos, 0)
            
            nsamp = self.bhead['num_samp']
            fmt_str = '%s%i%s' % (self.fmt_str[0], nsamp, self.fmt_str[1])
            
            buf = pack(fmt_str, *tdata)
            fd.write(buf)
    
    
    def write_tdata_multi(self, tracenums, tdata):
        """
        Writes trace data

        tracenums = array of trace numbers in file corresponding to the tdata
                    trace array (zero indexed)
        tdata = 2d array (or list) containing traces to be written where axis 1
                is the trace number and axis 2 is the sample number.
        """
        
        # recast tracenums to int64 to avoid overflowing long ints
        tracenums = np.array(tracenums, dtype='int64')
        
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
    

    def _ieee2ibm_b_new(ieee_float_trace, nsamp=-1):
        """
        Convert IEEE Float (big endian byte order) to IBM Float.
        """
        
        if nsamp == -1:
            # Assuming you have a way to determine nsamp from your data
            nsamp = len(ieee_float_trace) // 4 
    
        # We need to process each float individually
        ibm_data = []
    
        for i in range(nsamp):
            # Unpack the 4-byte IEEE float
            raw_bytes = ieee_float_trace[i*4:(i+1)*4]
            ieee_float = unpack('>f', raw_bytes)[0]
            
            # Handle zero
            if ieee_float == 0.0:
                ibm_data.extend([0, 0, 0, 0])
                continue
    
            # Get sign
            sign = 0
            if ieee_float < 0:
                sign = 128
                ieee_float = -ieee_float
    
            # Handle a subnormal number if necessary
            # Normalize the number for conversion to IBM's base-16
            if ieee_float < 1.0:
                 # Find the exponent and mantissa
                 exponent = 0
                 while ieee_float < 0.0625:  # 1/16
                     ieee_float *= 16.0
                     exponent -= 1
            else:
                 exponent = 0
                 while ieee_float >= 1.0:
                     ieee_float /= 16.0
                     exponent += 1
            
            # Adjust mantissa for normalization
            mantissa = ieee_float * (16**6)
            
            # Calculate the biased exponent
            ibm_exponent = exponent + 64
            
            # Unpack the mantissa into 3 bytes
            mant_bytes = pack('>I', int(mantissa))[1:]
            
            # Assemble the IBM float
            ibm_bytes = [sign + ibm_exponent] + list(mant_bytes)
            ibm_data.extend(ibm_bytes)
        
        return bytes(ibm_data)
    
    
    def _ieee2ibm_l_new(ieee_float_trace, nsamp=-1):
        """
        Convert IEEE Float (little endian byte order) to IBM Float.
        """
    
        if nsamp == -1:
            nsamp = len(ieee_float_trace) // 4
    
        ibm_data = []
    
        for i in range(nsamp):
            # Unpack the 4-byte IEEE float using little-endian byte order
            raw_bytes = ieee_float_trace[i*4:(i+1)*4]
            ieee_float = unpack('<f', raw_bytes)[0]
    
            # Handle zero
            if ieee_float == 0.0:
                ibm_data.extend([0, 0, 0, 0])
                continue
    
            # Get sign
            sign = 0
            if ieee_float < 0:
                sign = 128
                ieee_float = -ieee_float
    
            # Normalize the number for conversion to IBM's base-16
            if ieee_float < 1.0:
                 # Find the exponent
                 exponent = 0
                 while ieee_float < 0.0625:  # 1/16
                     ieee_float *= 16.0
                     exponent -= 1
            else:
                 exponent = 0
                 while ieee_float >= 1.0:
                     ieee_float /= 16.0
                     exponent += 1
            
            # Adjust mantissa for normalization
            mantissa = ieee_float * (16**6)
            
            # Calculate the biased exponent
            ibm_exponent = exponent + 64
            
            # Unpack the mantissa into 3 bytes
            mant_bytes = pack('<I', int(mantissa))[:3]
            
            # Assemble the IBM float in little-endian byte order
            ibm_bytes = [int(mant_bytes[0]), int(mant_bytes[1]), int(mant_bytes[2]), sign + ibm_exponent]
            ibm_data.extend(ibm_bytes)
    
        return bytes(ibm_data)



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
    
    bhead[8] = int(samp_rate)
    bhead[10] = int(num_samp)
    bhead[12] = int(fmt_code)
    
    print(bhead)
    
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
        for each in ehead:
            fd.write(each)
            
    with open(filename, 'a+b') as fd:
        
        # Seek past the EBCDIC header
        fd.seek(3199, 0)
        
        #   Do appropriate things for big endian
        if endian=='big':
            print('Working on big endian branch...')
            fmt_str = '>' + fmt_str
            
            for each in bhead:
                buf = pack('>h', each)
                fd.write(buf)
                
                fd.seek(-2, 1)
                buf = fd.read(2)
                eachin = unpack('>h', buf)[0]
                print(fd.tell(), each, eachin)
                fd.seek(2, 1)
            
            count = 0
            for i in range(0, num_trace):
                count += 1
                if (verbose > 0) & (count == verbose+1):
                    print('writing trace %i of %i' % (i, num_trace))
                    count = 1
                    
                for each in thead:
                    buf = pack('>h', each)
                    fd.write(buf)
                    
                for each in tdata:
                    buf = pack(fmt_str, each)
                    fd.write(buf)
    
        #   Do appropriate things for big endian
        if endian=='little':
            print('Working on little endian branch...')
            fmt_str = '<' + fmt_str
            
            for each in bhead:
                buf = pack('<h', each)
                fd.write(buf)

                fd.seek(-2, 1)
                buf = fd.read(2)
                eachin = unpack('<h', buf)[0]
                print(fd.tell(), each, eachin)
                fd.seek(2, 1)
            
            count = 0
            for i in range(0, num_trace):
                count += 1
                if (verbose > 0) & (count == verbose+1):
                    print('writing trace %i of %i' % (i, num_trace))
                    count = 1
                    
                for each in thead:
                    buf = pack('<h', each)
                    fd.write(buf)
                    
                for each in tdata:
                    buf = pack(fmt_str, each)
                    fd.write(buf)


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


def plot_wigva1(ax, tdata, t, trcstart=0, excursion=1, peak=False, trough=False, 
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
                                     facecolor=pcolor, edgecolor=None,
                                     alpha=alpha)
                
                # plot trough fill
                if trough==True:
                    ax.fill_betweenx(t2, zeroval, trc2, where=trc2<zeroval, 
                                     facecolor=tcolor,  edgecolor=None,
                                     alpha=alpha)
            if line == True:
                ax.plot(trc, t, lcolor, lw=lw)


def plot_wigva(ax, tdata, t, trcstart=0, excursion=1, peak=False, trough=False, 
          line=True, lw=0.25, lcolor='k', pcolor=[0.2, 0.2, 1.0],
          tcolor=[1.0, 0.2, 0.2], alpha=0.5, clipval=None):
    """
    Plot wiggle traces with variable area fill. 
    If using matplotlib version 2.1.0 and later, use this function.
    """
        
    tdata = np.array(tdata)
    ntrc, nsamp = tdata.shape
    
    norm = np.max(np.abs([np.nanmax(tdata), np.nanmin(tdata)]))
    if clipval != None:
        norm = clipval
        
    for i in range(0, ntrc):
        
        zeroval = i + trcstart
        trc = tdata[i, :]
        
        if np.isnan(trc).any():
            next
        
        else:
            #norm = max(abs([max(trc), min(trc)]))
            trc = trc/norm * excursion + zeroval
            
            if (peak==True) | (trough==True):
                                
                # plot peak fill
                if peak==True:
                    ax.fill_betweenx(t, zeroval, trc, where=trc>zeroval, 
                                     facecolor=pcolor, edgecolor=None,
                                     alpha=alpha, interpolate=True)
                
                # plot trough fill
                if trough==True:
                    ax.fill_betweenx(t, zeroval, trc, where=trc<zeroval, 
                                     facecolor=tcolor,  edgecolor=None,
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
