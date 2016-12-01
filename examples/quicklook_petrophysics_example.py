"""
Example script to illustrate how auralib may be used to perform quick-look
petrophysical analysis.

Written by: Wes Hamlyn
Created:    1-Dec-2016
"""

import numpy as np
import matplotlib.pyplot as plt
import auralib as aura


# Read data from an LAS file

# Specify the path to an input LAS file
infile = r'C:\Users\whamlyn\Dropbox\data\Blackfoot\09-17.LAS'

# Create the LASReader object
lasbuf = aura.las.LASReader(infile)

# Assign log curves to variables for convenience
depth = lasbuf.curves['DEPTH']
gam = lasbuf.curves['GR']
rho = lasbuf.curves['RHOB'] * 0.001 # convert from kg/m3 to g/cc
resis = lasbuf.curves['ILD']

# calculate Vshale using Gamma Ray Index and Larionov (older rocks) methods
Vsh_gri = aura.pp.vsh_gr(gam, 20, 180)
Vsh_larold = aura.pp.vsh_larionov_older(gam, 20, 180)

# calculate density porosity using sandstone, limestone, and dolomite matricies
DPSS = aura.pp.dens_por(rho, Rmatrix=2.65, Rfluid=1.00)
DPLS = aura.pp.dens_por(rho, Rmatrix=2.71, Rfluid=1.00)
DPDL = aura.pp.dens_por(rho, Rmatrix=2.87, Rfluid=1.00)

# calculate water resistivity
Tc = 60.0       # Temperature (Celcius)
Wse = 60000.0   # Salinity (ppm)
Rw = aura.pp.calc_Rw2(Tc, Wse, verbose=True)

# calculate water saturation using Archie and Modified Simandoux method
Rsh = 11.0  # read resistivity from ILD log in a shaley zone
Sw_archie = aura.pp.sw_archie(resis, DPSS, Rw, a=1.0, m=2.0, n=2.0)
Sw_ms = aura.pp.sw_modsim(resis, DPSS, Vsh_larold, Rw, Rsh, a=1.0, m=2.0, n=2.0)


# Plotting code below here
fig = plt.figure(num=1)
fig.clf()

ax1 = fig.add_subplot(161)
ax2 = fig.add_subplot(162, sharey=ax1)
ax3 = fig.add_subplot(163, sharey=ax1)
ax4 = fig.add_subplot(164, sharey=ax1)
ax5 = fig.add_subplot(165, sharey=ax1)
ax6 = fig.add_subplot(166, sharey=ax1)

ax1.plot(gam, depth, label='GR')
ax1.set_xlabel('Gamma Ray')
ax1.legend(loc='upper right', fontsize=10)

ax2.plot(rho, depth, label='RHOB')
ax2.set_xlabel('Density')
ax2.legend(loc='upper right', fontsize=10)

ax3.semilogx(resis, depth, label='ILD')
ax3.set_xlabel('Resistivity')
ax3.grid(which='minor')
ax3.legend(loc='upper right', fontsize=10)

ax4.plot(Vsh_gri, depth, 'g', label='Gamma Ray Index')
ax4.plot(Vsh_larold, depth, 'c', label='Larionov (Older)')
ax4.set_xlabel('Vshale')
ax4.legend(loc='upper right', fontsize=10)

ax5.plot(DPSS, depth, 'b', label='DPSS')
ax5.plot(DPLS, depth, 'c', label='DPLS')
ax5.set_xlabel('Total Porosity')
ax5.legend(loc='upper right', fontsize=10)

ax6.plot(Sw_archie, depth, 'b', label='Archie')
ax6.plot(Sw_ms, depth, 'g', label='Modified Simandoux')
ax6.set_xlabel('Swater')
ax6.legend(loc='upper right', fontsize=10)
           
           
ax1.invert_yaxis()
ax1.set_ylabel('MD')
for ax in fig.get_axes():
    ax.grid(True)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    
fig.tight_layout()
plt.show()