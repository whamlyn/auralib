# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 17:09:07 2021

@author: wesha
"""

#%%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import MultiCursor
import auralib as aura
from scipy.interpolate import interp1d

# load vp-vs-rho logs
infile = r'D:\SERVICE\Tallman\04_RokDoc\exports\las\111041205125W300_MD.las'

buf = aura.las.LASReader(infile)
md = buf.curves['DEPT']
vp = buf.curves['Vp_SYN']
rho = buf.curves['Rho_SYN']

# 1) Strip leading and trailing nulls from Vp and Rho logs

leading_vp = np.argmin(np.isnan(vp))
leading_vp_idx = np.arange(leading_vp)
vp1 = np.delete(vp, leading_vp_idx)

trailing_vp = np.argmin(~np.isnan(vp1))
trailing_vp_idx = np.arange(trailing_vp, len(vp1))
vp1 = np.delete(vp1, trailing_vp_idx)
vp1z = np.delete(np.delete(md, leading_vp_idx), trailing_vp_idx)

leading_rho = np.argmin(np.isnan(rho))
leading_rho_idx = np.arange(leading_rho)
rho1 = np.delete(rho, leading_rho_idx)

trailing_rho = np.argmin(~np.isnan(rho1))
trailing_rho_idx = np.arange(trailing_rho, 0)
rho1 = np.delete(rho1, trailing_rho_idx)
rho1z = np.delete(np.delete(md, leading_rho_idx), trailing_rho_idx)


# 2) Resample logs such that they now cover the same Z-range

z_min = np.max([vp1z[0], rho1z[0]])
z_max = np.min([vp1z[-1], rho1z[-1]])
dz = 0.1
z = np.arange(z_min, z_max+dz, dz)

# note: The below interpolators ignore any remaining nulls meaning that any
#       gaps in the vp and rho logs will be interpolated. This may or may not
#       be your desired behaviour, particuarly if large gaps are present.
fvp = interp1d(vp1z[~np.isnan(vp1)], vp1[~np.isnan(vp1)], kind='linear',
               bounds_error=False, fill_value=(vp1[0], vp1[-1]))

frho= interp1d(rho1z[~np.isnan(rho1)], rho1[~np.isnan(rho1)], kind='linear',
               bounds_error=False, fill_value=(rho1[0], rho1[-1]))

vp_syn = fvp(z)
rho_syn = frho(z)

ai_syn = vp_syn*rho_syn
sm_len_samp = 51
ai_syn_sm = aura.utils.smooth_log(ai_syn, sm_len_samp)

# do some despiking of the AI log
if False:
    spike_tol = 1000
    idx = np.nonzero((ai_syn > ai_syn_sm+spike_tol) | (ai_syn < ai_syn_sm-spike_tol))
    ai_syn[idx] = np.nan

    f = interp1d(z[~np.isnan(ai_syn)], ai_syn[~np.isnan(ai_syn)], kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    ai_syn = f(z)
    ai_syn_sm = aura.utils.smooth_log(ai_syn, sm_len_samp)


# 3) Calculate the initial time-depth transform from the vp_syn log

twt_dz = np.cumsum(dz/(vp_syn))*2
dt = 0.002
twt_dt = np.arange(twt_dz.min(), twt_dz.max(), dt)


# 4) Convert log digits from depth sampling to time sampling. this makes use
#    of the Depth2Time() class

class Depth2Time():
    """
    Class for depth to time conversion
    """
    def __init__(self, twt_dz, twt_dt):
        """
        Establish relation beteen depth-sampled and time-sampled twt values
        """
        self.twt_dz = twt_dz
        self.twt_dt = twt_dt

    def conv(self, data_dz):
        """
        Resample depth logs to twt
        """
        f = interp1d(self.twt_dz, data_dz, kind='linear', bounds_error=True)
        data_dt = f(self.twt_dt)

        return data_dt

d2t = Depth2Time(twt_dz, twt_dt)

vp_syn_t  = d2t.conv(vp_syn)
rho_syn_t = d2t.conv(rho_syn)
ai_syn_t  = d2t.conv(ai_syn)
ai_syn_sm_t  = d2t.conv(ai_syn_sm)


# 5) Compute zero offset reflection coefficients and add a zero to the start of
#    the reflectivity arrays to account for no reflection at the top of the
#    first layer.

rpp_syn_sm_t = (ai_syn_sm_t[1:] - ai_syn_sm_t[:-1])/(ai_syn_sm_t[1:] + ai_syn_sm_t[:-1])
rpp_syn_sm_t = np.hstack([0, rpp_syn_sm_t])


# 6) Create a wavelet
f1 = 5.0
f2 = 12.0
f3 = 105.0
f4 = 125.0
phase = 0.0
wvlt_length_samp = 75
wvlt_length_sec = (wvlt_length_samp-1)*dt
wvlt_t, wvlt_a = aura.wvlt.wvlt_bpass(f1, f2, f3, f4, phase, dt, wvlt_length_sec)

wvlt_a = wvlt_a * aura.win.cosine(wvlt_length_samp, 10)


# 7) Convolve wavelet with reflectivity series

synth = np.convolve(wvlt_a, rpp_syn_sm_t, mode='same')


# 9)  Calibrate time-depth relationship

well_twt_shift = 0.242
twt_dt = twt_dt + well_twt_shift
twt_dz = twt_dz + well_twt_shift


#
# PLOTTING CODE BELOW...
#

fig = plt.figure(num=1)
fig.clf()

nr = 1; nc = 5
ax = [plt.subplot2grid((nr, nc), (0, 0))]
for i in range(1, nc):
    ax.append(plt.subplot2grid((nr, nc), (0, i), sharey=ax[0]))

# plot original AI logs in TWT but sampled in Depth
ax[0].plot(ai_syn, twt_dz, c='k', lw=0.75)
ax[0].plot(ai_syn_sm, twt_dz, c='r', lw=2)
#ax[0].set_xlim(3000, 10000)
ax[0].set_xlabel('AI\n(m/s*g/cc)')

# plot upscaled AI log sampled in TWT
ax[1].step(ai_syn_sm_t, twt_dt, where='pre', c='r')
#ax[1].set_xlim(3000, 10000)
ax[1].set_xlabel('AI\n(m/s*g/cc)')

# plot reflectivity series
ax[2].hlines(twt_dt, 0.0, rpp_syn_sm_t, 'r')
ax[2].plot(np.zeros_like(twt_dt), twt_dt, 'r')
ax[2].set_xlim(-0.5, 0.5)
ax[2].set_xlabel('Refl. Coeff.\n(unitless)')

# plot wavelet
wvlt_t2 = wvlt_t + 0.4
aura.syn.plot_wigva(ax[3], wvlt_a, wvlt_t2, repeat=1)
ax[3].set_xlabel('Wavelet\n(amplitude)')

# plot synthetic trace (replicated several times)
aura.syn.plot_wigva(ax[4], synth, twt_dt, repeat=5)
ax[4].set_xlabel('Synthetic\n(amplitude)')


for each in ax:
    each.invert_yaxis()

aura.plot.format_log_axes(ax, '')
#ax[0].set_ylim(1, 0.0)

curs = MultiCursor(fig.canvas, ax, horizOn=True, vertOn=False, c='k', ls='-', lw=1)
plt.show()