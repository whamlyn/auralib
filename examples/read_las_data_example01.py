"""
Example script to illustrate how auralib may be used to read data from
LAS format data files.

Written by: Wes Hamlyn
Created:    1-Dec-2016
"""


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import auralib as aura


# Specify the path to an input SEG-Y file
infile = r'C:\Users\whamlyn\Dropbox\data\Blackfoot\09-17.LAS'

# Create the LASReader object
lasbuf = aura.las.LASReader(infile)


# All fields in the LAS file can be accessed as attributes of the LASReader 
# object.  Each attribute will be a Python Struct where each key corresponds to
# the constants in the LAS header.  Each constant has its data, units, and 
# description stored as Structs.

print('\nPrint Well Info section of LAS file:')
print('-------------------------------------')
for key in lasbuf.well_info:
    print('%s: %s' % (key, lasbuf.well_info[key]))


print('\nPrint individual constants from the LAS file')
print('----------------------------------------------')
print('Well Name: %s' % lasbuf.well_info['WELL']['data'])
print('STRT: %s' % lasbuf.well_info['STRT']['data'])
print('STOP: %s' % lasbuf.well_info['STOP']['data'])
print('STEP: %s' % lasbuf.well_info['STEP']['data'])


# All constants read from the LAS file are stored as strings in the LASReader
# object.  If you want to do operations on these, they must first be converted
# to numeric types (e.g. floating point, integers, etc.).  Below is a simple 
# example of this.
print('\nTotal logged depth')
print('--------------------')
STRT = float(lasbuf.well_info['STRT']['data'])
STOP = float(lasbuf.well_info['STOP']['data'])
TLD = STOP - STRT
print('Total Logged Depth: %f' % TLD)


# List information about the LAS curves
print('\nPrint curve mnemonics in the LAS file:')
print('-------------------------------------')
for mnemonic in lasbuf.curve_names:
    print(mnemonic)

print('\nPrint Curve Info section of LAS file:')
print('-------------------------------------')
for key in lasbuf.curve_info:
    print('%s: %s' % (key, lasbuf.curve_info[key]))


# Well logs are just about the only type of data read from LAS files and 
# converted to numeric arrays in the LASReader object.  Below are some examples
# of how to access the log digits, do calculations on the logs, and plot them
# to a matplotlib figure.

# Set log digits equal to variables for convenience (note that the log 
# mnemonics may well be different in your LAS file).  Below we use mnemonics 
# for Depth log='DEPTH', sonic log='DT', and density log='RHOB'.
depth = lasbuf.curves['DEPTH']
psonic = lasbuf.curves['DT']
density = lasbuf.curves['RHOB'] * 0.001 # convert to g/cc

# convert p-wave sonic to p-wave velocity
vp = 1000.0/psonic # Vp will be in km/s

# compute p-wave impedance from velocity and density logs
pimpedance = vp * density


# All plotting code below here...
fig = plt.figure(num=1)
fig.clf()

ax1 = fig.add_subplot(141)
ax2 = fig.add_subplot(142, sharey=ax1)
ax3 = fig.add_subplot(143, sharey=ax1)
ax4 = fig.add_subplot(144, sharey=ax1)

ax1.plot(psonic, depth)
ax1.set_xlabel('P-Sonic')

ax2.plot(vp, depth)
ax2.set_xlabel('Vp')

ax3.plot(density, depth)
ax3.set_xlabel('Density')

ax4.plot(pimpedance, depth)
ax4.set_xlabel('P-Impednace')

for ax in fig.get_axes():
    ax.grid(True)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    
fig.tight_layout()
plt.show()

