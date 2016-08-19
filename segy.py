"""
AuraQI module for reading data stored in SEG-Y format.

Author:   Wes Hamlyn
Created:   1-Sep-2014
Last Mod: 17-Aug-2016
"""

from struct import unpack, pack
from numpy import size, array, arange, abs, flipud, round, min, max
import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import sys



def_bhead = {'samp_rate':{'bpos':17, 'fmt':'h', 'nbyte':2},
             'num_samp':{'bpos':21, 'fmt':'h', 'nbyte':2},
             'samp_fmt':{'bpos':25, 'fmt':'h', 'nbyte':2}
             }

def_thead = {'il':{'bpos':189,  'fmt':'l', 'nbyte':4},
             'xl':{'bpos':193, 'fmt':'l', 'nbyte':4},
             'cmpx':{'bpos':181, 'fmt':'l', 'nbyte':4},
             'cmpy':{'bpos':185, 'fmt':'l', 'nbyte':4},
             'offset':{'bpos':37, 'fmt':'l', 'nbyte':4}
             }



class Segy(object):
    """
    Python class for SEG-Y file manipulation
    """

    def __init__(self, filename, def_bhead=def_bhead, def_thead=def_thead):
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

        if fmt_code == 1:
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
        


    def get_ilxl(self, il0, xl0):
        """
        Find an inline and crossline number in a segy file using a simple
        binary search.
        """
        
		# create a composite inline-crossline number by scaling the inline 
		# number by a multiple of 10 and then adding the crossline number
        ilxl0 = il0*10000 + xl0
        
        # make an initial guess at a trace
        tmin = 0
        tmax = self.num_traces
        
        count = 0
        ilxlg = -1
        while (ilxl0 != ilxlg) & (count<100):
            count += 1
            print('iteration %i' % count)
            tg = int((tmax+tmin)/2)
            self.read_thead1(tg)
            ilg = self.thead['il'][0]
            xlg = self.thead['xl'][0]
            ilxlg = ilg*10000 + xlg
            print('%i %i' % (ilg, xlg))
            
            if ilxlg == ilxl0:
                print('---Success: IL: %i XL: %i Trace: %i' % (il0, xl0, tg))
                
            elif ilxlg > ilxl0:
                tmax = tg
                
            else:
                tmin = tg
                
            trc_num = tg
            
            return trc_num
            
            
            
    def build_ilxl_index(self, verbose=1000):
        """
        Build indicies for 3D segy files that have IL and XL number in trace
        headers.  Output will be an object attribute (structure) containing 
        indices and inline/crossline ranges and counts.
        
        self.idx['ilmin'] = minimum inline number in segy file
        self.idx['ilmax'] = maximum inline number in segy file
        self.idx['numil'] = number of inlines in segy file
        
        self.idx['xlmin'] = minimum xline number in segy file
        self.idx['xlmax'] = maximum xline number in segy file
        self.idx['numxl'] = number of xlines in segy file
        
        self.idx['tidx'] = list of traces in each cdp
        
        Example:
        
            To find trace indicies (i.e. trace number - 1) at the 1st inline  
            and 5th xline: self.idx['tidx'][0][4]  (returns a 1D list)
        
        """
        
        from copy import deepcopy
        from numpy import min, max
        
        # read in all trace headers
        self.read_thead2(0, self.num_traces-1, verbose)
        
        # start building the index structure
        idx = {}
        idx['ilmin'] = min(self.thead['il'])
        idx['ilmax'] = max(self.thead['il'])
        idx['num_il'] = idx['ilmax'] - idx['ilmin'] + 1
        
        idx['xlmin'] = min(self.thead['xl'])
        idx['xlmax'] = max(self.thead['xl'])
        idx['num_xl'] = idx['xlmax'] - idx['xlmin'] + 1
        
        xlbuf = []
        for j in range(0, idx['num_xl']):
            xlbuf.append([])
        
        tidx = []
        for i in range(0, idx['num_il']):
            tidx.append(deepcopy(xlbuf))
        
        # now actually fill each element of the index lists
        count = 0
        for k in range(0, self.num_traces-1):

            count = count + 1
            if (verbose > 0) & (count == verbose+1):
                print("reading trace %i of %i" % (k+1, self.num_traces))
                count = 1
                
            ilidx = self.thead['il'][k] - idx['ilmin']
            xlidx = self.thead['xl'][k] - idx['xlmin']
            tidx[ilidx][xlidx].append(k)
            
        idx['tidx'] = tidx
        
        # save the index structure as an object attribute
        self.idx = idx
        
        #   Now clean up the thead attribute so it is back to empty
        self.thead = {}
        for key in self.def_thead.keys():
            self.thead.update({key:[]})
    
    
    
    def get_ilxl_idx(self, inline, xline):
        """
        Method to return the trace number(s) in a 3D segy file by using the idx
        object attribute if it has been intialized.
        """
        
        ilmin = self.idx['ilmin']
        ilmax = self.idx['ilmax']
        xlmin = self.idx['ilmin']
        xlmax = self.idx['xlmax']
        
        if ((inline >= ilmin) & (inline <= ilmax) & 
        (xline >= xlmin) & (xline <= xlmax)):
            
            il = inline
            xl = xline
            
            ili = il - self.idx['ilmin']
            xli = xl - self.idx['xlmin']
            
            tnum = self.idx['tidx'][ili][xli]
            
            if len(tnum) == 1:
                tnum = tnum[0]
                
            
        else:
            
            print('The inline or crossline value is out of range')
            print('Select inline values from %5i to %5i' % (ilmin, ilmax))
            print('Select xline values from  %5i to %5i' % (xlmin, xlmax))
        
        return tnum
    
    
    
    def read_inline(self, inline):
        """
        Method to read trace data from the specified inline.  Requires the idx
        segy object attribute to be initialized.
        """
        
        xlmin = self.idx['xlmin']
        xlmax = self.idx['xlmax']
        
        # get trace numbers in SEG-Y file
        buf_tnum = []
        for xline in range(xlmin, xlmax+1):
            tnum = self.get_ilxl_tnum(inline, xline)
            buf_tnum.append(tnum)

        tdata = []        
        # read trace data from SEG-Y file
        for tr in buf_tnum:
            tdata.append(self.read_trace_data(tr))
        
        return tdata



    def read_xline(self, xline):
        """
        Method to read trace data from the specified xline.  Requires the idx
        segy object attribute to be initialized.
        """
        
        ilmin = self.idx['ilmin']
        ilmax = self.idx['ilmax']
        
        # get trace numbers in SEG-Y file
        buf_tnum = []
        for inline in range(ilmin, ilmax+1):
            tnum = self.get_ilxl_tnum(inline, xline)
            buf_tnum.append(tnum)

        tdata = []        
        # read trace data from SEG-Y file
        for tr in buf_tnum:
            tdata.append(self.read_trace_data(tr))
        
        return tdata
    
    
    
    def plot_inline(self, inline, palette=plt.cm.bwr_r, amp_min=0, amp_max=0, clip=1.0,
                    tmin=0, tmax=0):
        """
        Method to quickly plot an inline
        """
        
        tdata = self.read_inline(inline)
        tdata = flipud(array(tdata).T)
        
        if (amp_min == 0) & (amp_max == 0):
            amp_min = tdata.min() * clip
            amp_max = tdata.max() * clip
        
        
        xmin = self.idx['xlmin']
        xmax = self.idx['xlmax']
        extent = [xmin, xmax,
                  0, self.bhead['num_samp'] * self.bhead['samp_rate']*0.001]
        
        if (tmin == 0) & (tmax == 0):
            tmin = extent[2]
            tmax = extent[3]
        
        fig = plt.figure(num=1)
        fig.clf()
        
        ax = fig.add_subplot(111)
        im = ax.imshow(tdata, cmap=palette, vmin=amp_min, vmax=amp_max, 
                       extent=extent, origin='upper')
        ax.set_aspect('auto')
        ax.set_ylim([tmax, tmin])
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_ticks_position('top')
        ax.set_xlabel('Crossline')
        ax.set_ylabel('TWT (ms)')
        
        fig.colorbar(im, ax=ax)
        
        plt.show()
    
    
    
    def read_trace_data(self, tracenum):
        """
        Read a single trace from a SEG-Y file.
        Note:
            tracenum is the number of the trace in the SEG-Y file to be read.
            This starts at zero (e.g. first SEG-Y trace has tracenum = 0)
        """

        bpos = 3200 + 400 + 240 + self.trace_size * tracenum

        fd = open(self.filename, 'rb')
        fd.seek(bpos, 0)

        tdata = []

        if self.fmt_str == '>ibm':
            for i in range(0, self.bhead['num_samp']):
                buf = fd.read(self.numbytes)
                buf = self._ibm2ieee_b(buf)
                tdata.append(buf)

        elif self.fmt_str == '<ibm':
            for i in range(0, self.bhead['num_samp']):
                buf = fd.read(self.numbytes)
                buf = self._ibm2ieee_l(buf)
                tdata.append(buf)

        else:
            for i in range(0, self.bhead['num_samp']):
                buf = fd.read(self.numbytes)
                buf = unpack(self.fmt_str, buf)
                tdata.append(buf[0])

        fd.close()

        return tdata


    def read_multi_trace_data(self, tr_start, tr_end, verbose=0):
        """
        Read multiple sequential traces from a SEG-Y file.  This method is faster
		than read_trace_data() but is restricted to reading sequential traces.
        Note:
            tr_start is the first trace to be read
            tr_end is the final trace to be read (inclusive)
            This function starts at zero (e.g. first SEG-Y trace has tracenum = 0)
        """            
        
        fd = open(self.filename, 'rb')

        tdata = []
        count = 0
        for tracenum in range(tr_start, tr_end):
            
            count += 1
            if (verbose > 0) & (count == verbose+1):
                print('reading trace %i' % (tracenum))
                count = 1
            
            bpos = 3200 + 400 + 240 + self.trace_size * tracenum 
            fd.seek(bpos, 0)
            
            tbuf = []
            
            if self.fmt_str == '>ibm':
                for i in range(0, self.bhead['num_samp']):
                    buf = fd.read(self.numbytes)
                    buf = self._ibm2ieee_b(buf)
                    tbuf.append(buf)
    
            elif self.fmt_str == '<ibm':
                for i in range(0, self.bhead['num_samp']):
                    buf = fd.read(self.numbytes)
                    buf = self._ibm2ieee_l(buf)
                    tbuf.append(buf)
    
            else:
                for i in range(0, self.bhead['num_samp']):
                    buf = fd.read(self.numbytes)
                    buf = unpack(self.fmt_str, buf)
                    tbuf.append(buf[0])
            
            tdata.append(tbuf)
            
        fd.close()

        return tdata



    def read_multi_trace_data_new(self, tr_start, tr_end, verbose=0):
        """
        Read multiple sequential traces from a SEG-Y file.  This method is faster
		than read_trace_data() but is restricted to reading sequential traces.
        Note:
            tr_start is the first trace to be read (zero-indexed)
            tr_end is the final trace to be read (inclusive)
            This method starts at zero (e.g. first SEG-Y trace has tracenum = 0)
        """            
        
        tdata = []
        count = 0
        if self.fmt_str == '>ibm':
            
            print('In >ibm block')
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
            print('In <ibm block')
            
        else:
            
            print('In Python supported sample format block')
            
            #  build format string
            if self.endian == 'big':
                tr_fmt_str = '>'
            else:
                tr_fmt_str = '<'
            
            for i in range(0, self.bhead['num_samp']):
                tr_fmt_str = tr_fmt_str + self.fmt_str[1:]
                
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
        

    def read_thead1(self, tracenum):
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


    def read_thead2(self,  tracenum1, tracenum2, verbose=0):
        """
        Read trace headers from sequential traces in a SEG-Y file.
        Note:
            tracenum2 is the number of the first trace to be read.
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
                print('reading trace %i of %i' % (i, tracenum2-tracenum1))
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
                

        fd = open(self.filename, 'rb+')

        abs_bpos = 3600 + bpos-1 + tracenum*self.trace_size
        fd.seek(abs_bpos, 0)
        buf = pack(fmt, data)
        fd.write(buf)

        fd.close()



    def write_trace(self, tracenum, data):
        """
        Writes trace data

        tracenum = trace number in file (zero indexed)
        bpos = byte position in trace header to be written
        fmt = f, l, h, c and > or <
        data = single value (int, float, double)
        """

        fd = open(self.filename, 'rb+')

        bpos = 3840 + tracenum*self.trace_size
        fd.seek(bpos, 0)

        for sample in data:
            buf = pack(self.fmt_str, sample)
            fd.write(buf)

        fd.close()


    def _ibm2ieee_b(self, ibm_float):
        """
        Convert IBM Float (big endian byte order)
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



    def _ibm2ieee_b_new(self, ibm_float_trace):
        """
        Convert IBM Float (big endian byte order)
        """
        
        #  build trace format string
        fmt_str = '>'
        for i in range(0, self.bhead['num_samp']):
            fmt_str = fmt_str + 'BBBB'
        
        buf = unpack(fmt_str, ibm_float_trace)
        
        buf = np.array(buf, dtype='int')
        
        istic = buf[0::4]
        a = buf[1::4]
        b = buf[2::4]
        c = buf[3::4]
        
        sign = np.ones(istic.shape)
        idx = np.nonzero(istic >= 128)
        sign[idx] = -1
        istic[idx] = istic[idx] - 128
        
        dividend = float(16**6)
        mant = a**(17.0) + b**(9.0) + c
        buf = sign* 16**(istic-64)*(mant/dividend)
        
        return buf.tolist()


    def _ibm2ieee_l_new(self, ibm_float_trace):
        """
        Convert IBM Float (big endian byte order)
        """
        
        #  build trace format string
        fmt_str = '>'
        for i in range(0, self.bhead['num_samp']):
            fmt_str = fmt_str + 'BBBB'
        
        buf = unpack(fmt_str, ibm_float_trace)
        
        dividend = float(16**6)
        
        buf = np.array(buf)
        
        istic = buf[0::4]
        a = buf[1::4]
        b = buf[2::4]
        c = buf[3::4]
        
        sign = np.ones(istic.shape)
        idx = np.nonzero(istic >= 128)
        sign[idx] = -1
        istic[idx] = istic[idx] - 128
            
        mant = a<<16 + b<<8 + c
        buf = sign*16**(istic - 64)*(mant/dividend)
        
        ieee_float_trace = buf.tolist()
        
        return ieee_float_trace        


def write_blank_segy_file(filename, fmt_code, samp_rate, num_samp, num_trace, endian='big'):
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
    with open(filename, 'wb') as fd:
        
        #   Do appropriate things for big endian
        if endian=='big':
            
            for buf in ehead:
                fd.write(buf)
    
            for buf in bhead:
                fd.write(pack('>h', buf))
    
            fmt_str = '>' + fmt_str
            
            for i in range(0, num_trace):
                for buf in thead:
                    fd.write(pack('>h', buf))
                    
                for buf in tdata:
                    fd.write(pack(fmt_str, buf))
    
        #   Do appropriate things for little endian
        if endian=='little':
            
            for buf in ehead:
                fd.write(buf)
    
            for buf in bhead:
                fd.write(pack('<h', buf))
    
            for i in range(0, num_trace):
                for buf in thead:
                    fd.write(pack('<h', buf))
    
                for buf in tdata:
                    fd.write(pack(fmt_str, buf))


#def segy_write_trace_data(filename, tracenum, tdata):




def plot_seis(ax, tdata, t_min=0, t_max=0, tr_min=0, tr_max=0, samp_rate=0.002,
              cmap=plt.cm.gray_r, amp_min=0, amp_max=0):
    
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



def plot_wigva(ax, tdata, t, trcstart=1, excursion=1, peak=False, trough=False, 
          lcolor='k', pcolor=[0.2, 0.2, 1.0], tcolor=[1.0, 0.2, 0.2]):
    """
    Python function to plot wiggle traces with variable area fill.
    """
    
    from scipy.interpolate import interp1d
    
    # estimate the sample rate; do it this way to avoid errors due to floating
    # point precision
    dt = int((t[1]-t[0])*1000)*0.001
    
    # set the new sample rate for the VA fill to 1/5 of the original sampling
    dt2 = dt*0.2
    
    tdata = np.array(tdata)
    ntrc, nsamp = tdata.shape
    
    norm = np.max(np.abs([np.max(tdata), np.min(tdata)]))
    for i in range(0, ntrc):
        
        zeroval = i + trcstart
        trc = tdata[i, :]
        #norm = max(abs([max(trc), min(trc)]))
        trc = trc/norm * excursion + zeroval
        
        if (peak == True) | (trough == True):
            # interpolate the new trace
            t2 = arange(t.min(), t.max(), dt2)
            f = interp1d(t, trc, kind='linear')
            trc2 = f(t2)
            
            # plot peak fill
            if peak==True:
                ax.fill_betweenx(t2, zeroval, trc2, where=trc2>=zeroval, 
                                 facecolor=pcolor, edgecolor=pcolor)
            
            # plot trough fill
            if trough==True:
                ax.fill_betweenx(t2, zeroval, trc2, where=trc2<=zeroval, 
                                 facecolor=tcolor,  edgecolor=tcolor)
        ax.plot(trc, t, lcolor)
        
    

def time_to_samp(twt, dt):
    """
    Convert an array of TWT to an array of sample indicies
    
    samp = time_to_samp(twt, dt)
    """
    
    samp = round(twt / dt)
    samp = array(samp, dtype='int32')
    
    return samp