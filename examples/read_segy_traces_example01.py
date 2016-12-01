"""
Example script to illustrate how auralib may be used to read trace data from
SEG-Y format data files.

Written by: Wes Hamlyn
Created:    1-Dec-2016
"""


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import auralib as aura


# Specify the path to an input SEG-Y file
infile = r'C:\Users\whamlyn\Dropbox\data\QSI\data\Project_data\cdps_line2.sgy'


# Create a Struct containing the important binary header fields, byte 
# positions, encoded formats, and number of bytes. Generally you should never 
# need to edit this.
def_bhead = {'samp_rate':{'bpos':17, 'fmt':'h', 'nbyte':2},
             'num_samp':{'bpos':21, 'fmt':'h', 'nbyte':2},
             'samp_fmt':{'bpos':25, 'fmt':'h', 'nbyte':2}
             }

# Create a Struct containing the trace header fields, byte positiosn, encoding 
# formats, and number of bytes. Generally you will always need to edit or 
# re-define this for each SEG-Y file.  Note that the more header fields you 
# define in this struct, the longer it will take to read from trace headers
def_thead = {'il':{'bpos':189,  'fmt':'l', 'nbyte':4},
             'xl':{'bpos':193, 'fmt':'l', 'nbyte':4},
             'cmpx':{'bpos':181, 'fmt':'l', 'nbyte':4},
             'cmpy':{'bpos':185, 'fmt':'l', 'nbyte':4},
             'offset':{'bpos':37, 'fmt':'l', 'nbyte':4}
             }



# Create an auralib Segy object.  Inputs required are the path to the input
# SEG-Y file, and the structs defining the binary and trace header definitions
segybuf = aura.segy.Segy(infile, def_bhead, def_thead)


# To view the ebcdic header, view the segybuf.ebcdic attribute
print(segybuf.ebcdic)



# To read a single trace from the SEG-Y file.  Returns a python list.
trc_number = 50
tdata = segybuf.read_trace_data_new(trc_number)


# To read multiple sequential traces (faster than reading one-by-one but
# traces have to be sequential in the SEG-Y file).  Returns a Python list
# with indexing as follows:  tdata[trace_number][sample_number]
trc_number_start = 50
trc_number_end = 200
tdata = segybuf.read_multi_trace_data_new(trc_number_start, trc_number_end)

# Note that the "_new" versions of read_trace_data() and 
# read_multi_trace_data() are still experimental. If they work, they will be 
# faster than the non "new" methods, but may not correctly handle all data 
# types correctly.  Use caution for SEG-Y files with sample amplitudes encoded 
# in IBM floating points format.



# To read the trace header there are two methods.  The read_thead1() method 
# reads the header fields defined in the def_thead Struct from a single trace.
# The read_thead2() method reads the header fields defined in the def_thead 
# Struct from a set of sequential traces.  Note that read_thead2() will be 
# faster to read multiple traces than read_thead1() but requires that the 
# traces be sequential. Trace header data are returned as a Struct with keys 
# corresponding to the def_thead key names.
thead = segybuf.read_thead1(trc_number)
thead = segybuf.read_thead2(trc_number_start, trc_number_end)



# Plotting example

# get some parameters for min/max time and traces and amplitude range
time_min = 0
time_max = segybuf.bhead['num_samp']*segybuf.bhead['samp_rate']*0.001
trc_min = trc_number_start
trc_max = trc_number_end
amp_max = np.max(np.abs(tdata));
amp_min = -amp_max


fig = plt.figure(num=1)
fig.clf()

ax1 = plt.subplot2grid((5, 9), (0, 0), colspan=8)
ax2 = plt.subplot2grid((5, 9), (1, 0), colspan=8, rowspan=4, sharex=ax1)
ax3 = plt.subplot2grid((5, 9), (1, 8), rowspan=4)

trc_nums = np.arange(trc_number_start, trc_number_end) + 1
ax1.plot(trc_nums, thead['offset'], 'k')
ax1.set_ylabel('Offset (m)')
ax1.grid(True)

img = ax2.imshow(np.transpose(tdata), cmap=cm.bwr_r, vmin=amp_min, vmax=amp_max,
                 extent=[trc_min, trc_max, time_max, time_min])
ax2.set_aspect('auto')
ax2.set_xlabel('Trace')
ax2.set_ylabel('TWT (ms)')

fig.colorbar(img, cax=ax3)
ax3.set_ylabel('amplitude')

fig.tight_layout()
plt.show()
