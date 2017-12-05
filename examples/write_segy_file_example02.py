# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 21:09:56 2017

@author: whamlyn
"""

#%%

# import some libraries
import numpy as np
import matplotlib.pyplot as plt
import auralib as aura


# First, lets create a seismic "trace" by using sine waves of different
# frequencies (10, 20, 50 Hz)

samp_rate = 0.002
num_samp = 251
t = np.arange(0, num_samp, 1)*samp_rate
a = np.sin(2*np.pi*30*t) + np.sin(2*np.pi*20*t) + np.sin(2*np.pi*50*t)


# Now, make 100 duplicate of that trace
num_traces = 100
tdata = []
for i in range(num_traces):
    tdata.append(a)


# Next create a SEG-Y file filled with zero-amplitude traces, a minimal  binary
# header and minimal trace headers.

outfile = r'C:\temp\output.sgy'
samp_fmt = 5 # writes out IEEE floating point traces
aura.segy.write_blank_segy_file_v2(outfile, samp_fmt, 2000, num_samp,
                                   num_traces, verbose=10)

# Now, create a segy object like we would use for reading segy files
buf = aura.segy.Segy(outfile)

# Create an array of trace numbers (0 to 99) and write the sine wave traces 
# to disk.
tracenums = np.arange(0, num_traces)
buf.write_trace_data_multi(tracenums, tdata)

# now read the traces from the SEG-Y file
tdata2 = buf.read_tdata_multi(0, buf.num_traces, 10)


# plot the data we just read in
fig = plt.figure(num=1)
fig.clf()

ax1 = fig.add_subplot(111)

aura.segy.plot_wigva(ax1, tdata2, t, trcstart=0,
                     excursion=1, peak=True, trough=True)
ax1.invert_yaxis()

fig.tight_layout()
plt.show()
