from kennett_cmp import model_kennett_cmp, ricker_wavelet, calc_min_n_slowness

import numpy as np
import matplotlib.pyplot as plt
import auralib as aura
import pandas as pd
import auralib as aura
from scipy.interpolate import PchipInterpolator, interp1d
from matplotlib.widgets import MultiCursor


def block_log(depth, values, block_boundaries):
    """
    Block a single well log using user-defined depth boundaries.

    Parameters
    ----------
    depth : array_like, shape (m,) or (1, m)
        Depth reference for the original log samples.

    values : array_like, shape (m,) or (1, m)
        Original log values. NaN values are ignored when calculating
        block averages.

    block_boundaries : array_like, shape (n,) or (1, n)
        Depths defining the block boundaries. The boundaries must be
        strictly increasing. There will be n - 1 output blocks.

    Returns
    -------
    blocked_depth : ndarray, shape (n - 1,)
        Midpoint depth of each block.

    blocked_values : ndarray, shape (n - 1,)
        Mean log value within each block. A block containing no valid
        log samples is assigned NaN.

    Notes
    -----
    Each block is defined as

        block_boundaries[i] <= depth < block_boundaries[i + 1]

    except for the final block, which includes its lower and upper
    boundaries.
    """
    depth = np.asarray(depth, dtype=float).reshape(-1)
    values = np.asarray(values, dtype=float).reshape(-1)
    boundaries = np.asarray(
        block_boundaries,
        dtype=float,
    ).reshape(-1)

    if depth.size != values.size:
        raise ValueError(
            "depth and values must contain the same number of samples."
        )

    if boundaries.size < 2:
        raise ValueError(
            "At least two block boundaries are required."
        )

    if not np.all(np.isfinite(depth)):
        raise ValueError("depth cannot contain NaN or infinite values.")

    if not np.all(np.isfinite(boundaries)):
        raise ValueError(
            "block_boundaries cannot contain NaN or infinite values."
        )

    if np.any(np.diff(depth) <= 0):
        raise ValueError("depth must be strictly increasing.")

    if np.any(np.diff(boundaries) <= 0):
        raise ValueError(
            "block_boundaries must be strictly increasing."
        )

    number_of_blocks = boundaries.size - 1

    # Use the midpoint as the representative depth of each block.
    #blocked_depth = 0.5 * ( boundaries[:-1] + boundaries[1:] )
    
    # Use the top of the interval as the representative depth of each block.
    blocked_depth = boundaries[:-1]

    blocked_values = np.full(
        number_of_blocks,
        np.nan,
        dtype=float,
    )

    for i in range(number_of_blocks):
        top = boundaries[i]
        base = boundaries[i + 1]

        # Make every block left-inclusive and right-exclusive.
        # Include the lower boundary of the final block as usual and
        # also include its upper boundary.
        if i < number_of_blocks - 1:
            in_block = (depth >= top) & (depth < base)
        else:
            in_block = (depth >= top) & (depth <= base)

        valid = in_block & np.isfinite(values)

        if np.any(valid):
            blocked_values[i] = np.mean(values[valid])

    return blocked_depth, blocked_values



# Depth is the top of each layer.
# The last layer extends downward as a half-space.


# Read in logs
lasfile = r'Q:\Wes\GEODATA\2026_IlleLake3D_Inversion\05_Scripts\MultipleModelling\103020904303W500\103020904303W500_MD.las'
buf = aura.las.LASReader(lasfile)

kb = 971.27
gl = 965.77

md = buf.curves['DEPT']
gr = buf.curves['HGR']
vp = buf.curves['Vp_Splice4Syn']
vs = buf.curves['Vs_Splice4Syn']
rho = buf.curves['Rho_Splice4Syn']


# Cut off some bad values at top of Vp/Vs/Rho logs
mask = md<=456.3
vp[mask] = np.nan
vs[mask] = np.nan
rho[mask] = np.nan

# Remove any nulls or values above ground
mask = (md>=(kb-gl)) & np.isfinite(gr) & np.isfinite(vp) & np.isfinite(vs) & np.isfinite(rho)
md = md[mask]
gr = gr[mask]
vp = vp[mask]
vs = vs[mask]
rho = rho[mask]


zstep = 0.1
z = np.arange(0, md[-1]+zstep, zstep)

fgr = PchipInterpolator(md, gr, extrapolate=True)
fgr = interp1d(md, gr, kind='linear', bounds_error=False,
               fill_value=(gr[0], gr[-1]))
gr = fgr(z)

fvp = PchipInterpolator(md, vp, extrapolate=True)
fvp = interp1d(md, vp, kind='linear', bounds_error=False,
               fill_value=(vp[0], vp[-1]))
vp = fvp(z)

fvs = PchipInterpolator(md, vs, extrapolate=True)
fvs = interp1d(md, vs, kind='linear', bounds_error=False,
               fill_value=(vs[0], vs[-1]))
vs = fvs(z)

frho = PchipInterpolator(md, rho, extrapolate=True)
frho = interp1d(md, rho, kind='linear', bounds_error=False,
                fill_value=(rho[0], rho[-1]))
rho = frho(z)

smooth_window_m = 2.0
smooth_window_samp = round(smooth_window_m / zstep)
gr = aura.utils.smooth_log(gr, smooth_window_samp)
vp = aura.utils.smooth_log(vp, smooth_window_samp)
vs = aura.utils.smooth_log(vs, smooth_window_samp)
rho = aura.utils.smooth_log(rho, smooth_window_samp)

block_step = 5.0
zblocks = np.arange(z[0], z[-1]+block_step, block_step)
zb, vpb = block_log(z, vp, zblocks)
zb, vsb = block_log(z, vs, zblocks)
zb, rhob = block_log(z, rho, zblocks)


dt = 0.002
wvlt_t, wvlt_a = aura.wvlt.wvlt_bpass(10.0, 12.0, 50.0, 60.0, 0.0, dt=dt, wvlt_length=0.100)
#wvlt_t, wvlt_a = aura.wvlt.ricker(50.0, 0.0, dt, 0.1)

# zb =  np.array([   0.0,  500.0,  750.0, 1200.0])
# vpb = np.array([2000.0, 3000.0, 2500.0, 3300.0])
# vsb = np.array([1000.0, 1350.0, 1150.0, 1800.0])
# rhob= np.array([2000.0, 2200.0, 2100.0, 2400.0])
offsets = np.arange(0.0, 4000.1, 25.0)


td_zbz = np.hstack((zb.copy(), 3000))
deltaZ = np.diff(td_zbz)
deltaT = np.hstack((0, deltaZ/vpb))
td_tbz = np.cumsum(deltaT)*2
td_tbt = np.arange(0, zb[-1]+dt, dt)
td_zbt = np.interp(td_tbt, td_tbz, td_zbz)

tz = np.interp(z, td_zbz, td_tbz)


def print_progress(fraction, completed, total, message):
    """Display modelling progress on one console line."""
    percent = 100.0 * fraction
    print(
        f"\r{percent:6.2f}%  "
        f"[{completed:,}/{total:,}]  "
        f"{message:<40}",
        end="",
        flush=True,
    )
    if fraction >= 1.0:
        print()

max_frequency = 30.0
n_slowness = calc_min_n_slowness(np.nanmax(vpb), max_frequency, np.max(offsets))
n_slowness = n_slowness + int(n_slowness*1.0)
n_slowness = 121

common_arguments = dict(
    depth=zb,
    vp=vpb,
    vs=vsb,
    rho=rhob,
    offsets=offsets,
    wavelet=wvlt_a,
    dt=dt,
    record_length=1.7,
    n_slowness=n_slowness,
    max_frequency=max_frequency,

    receiver_component="all",
    complex_frequency_damping="auto",
    slowness_taper_fraction=0.05,
    frequency_taper_fraction=0.05,
    nmo_corrected=True,
    stretch_mute=0.3,
    normalize=None,
    progress_callback=print_progress,
)


# Primaries only with no mode conversions (p-wavefield only)
result1 = model_kennett_cmp(
    **common_arguments,
    response_type="primary",  # no multiples
    free_surface=True,        # include free surface effect
    mode_coupling="p_only",   # no internal mode conversions (P-P only)
)

# Primaries only with mode conversions (p and sv wavefields)
result2 = model_kennett_cmp(
    **common_arguments,
    response_type="primary",  # no multiples
    free_surface=True,        # include free surface effect
    mode_coupling="full",     # include internal P and Sv mode conversions
)

# Primaries + surface multiples
result3 = model_kennett_cmp(
    **common_arguments,
    response_type="total",         # primaries and multiples
    free_surface=True,             # include free surface effect
    surface_multiple_order='all',  # include all orders of surface-related multiples
    internal_multiples=False,      # no internal multiples
    mode_coupling="p_only",        # no internal mode conversions (P-P only)
)

# Primaries + surface multiples + internal multiples
result4 = model_kennett_cmp(
    **common_arguments,
    response_type="total",         # primaries and multiples
    free_surface=True,             # include free surface effect
    surface_multiple_order='all',  # include all orders of surface-related multiples
    internal_multiples=True,       # include all internal multiples
    mode_coupling="p_only",        # no internal mode conversions (P-P only)
)

# Primaries + surface multiples + internal multiples + mode conversions
result5 = model_kennett_cmp(
    **common_arguments,
    response_type="total",         # primaries and multiples
    free_surface=True,             # include free surface effect
    surface_multiple_order='all',  # include all orders of surface-related multiples
    internal_multiples=True,       # include all internal multiples
    mode_coupling="full",          # include internal P and Sv mode conversions
)


#%%

#gather_type = "vertical_pp"
#gather_type = "vertical_ps"
gather_type = "vertical_velocity"

pp_primaries =         result1.gathers[gather_type] * -1
psv_primaries =        result2.gathers[gather_type] * -1
pp_multiples_fsonly =  result3.gathers[gather_type] * -1
pp_multiples_total =   result4.gathers[gather_type] * -1
psv_multiples_total =  result5.gathers[gather_type] * -1

if gather_type == 'vertical_pp':
    fig_title = 'Vertical Geophone Component | P Wavefield Only | Particle Velocity'
    
elif gather_type == 'vertical_ps':
    fig_title = 'Vertical Geophone Component | Sv Wavefield Only | Particle Velocity'

elif gather_type == 'vertical_velocity':
    fig_title = 'Vertical Geophone Component | P+Sv Wavefields | Particle Velocity'

# Converted wave events
psv_converted_waves = psv_primaries - pp_primaries

# Surface-related multiple events
surface_multiples = pp_multiples_fsonly - pp_primaries

# Internal multiple events
internal_multiples = pp_multiples_total - pp_multiples_fsonly

# Complete multiple and converted wave events
all_non_pp_primaries = psv_multiples_total - pp_primaries


clip = np.quantile(np.abs(psv_multiples_total), 0.99)
#clip = 0.002

fig = plt.figure(figsize=(10, 7), num=1)
fig.clf()

fig_title = fig_title + '\n%s' % buf.well_info['WELL']['data']
fig.suptitle(fig_title, fontsize=14, fontweight='bold')

nr = 1
nc = 17

ax1 = plt.subplot2grid((nr, nc), (0, 0))
ax2 = plt.subplot2grid((nr, nc), (0, 1), sharey=ax1)
ax3 = plt.subplot2grid((nr, nc), (0, 2), sharey=ax1)
ax4 = plt.subplot2grid((nr, nc), (0, 3), sharey=ax1, colspan=2)
ax5 = plt.subplot2grid((nr, nc), (0, 5), sharey=ax1, colspan=2)
ax6 = plt.subplot2grid((nr, nc), (0, 7), sharey=ax1, colspan=2)
ax7 = plt.subplot2grid((nr, nc), (0, 9), sharey=ax1, colspan=2)
ax8 = plt.subplot2grid((nr, nc), (0, 11), sharey=ax1, colspan=2)
ax9 = plt.subplot2grid((nr, nc), (0, 13), sharey=ax1, colspan=2)
ax10 = plt.subplot2grid((nr, nc), (0, 15), sharey=ax1, colspan=2)


ax1.plot(vp, tz, 'k', lw=1)
ax1.step(np.hstack((vpb, vpb[-1])), td_tbz, where='pre', c='r', lw=2)
ax1.set_title('P-wave Velocity')
ax1.set_xlabel('(m/s)')

ax2.plot(vs, tz, 'k', lw=1)
ax2.step(np.hstack((vsb, vsb[-1])), td_tbz, where='pre', c='r', lw=2)
ax2.set_title('S-wave Velocity')
ax2.set_xlabel('(m/s)')

ax3.plot(rho*1000, tz, 'k', lw=1)
ax3.step(np.hstack((rhob, rhob[-1]))*1000, td_tbz, where='pre', c='r', lw=2)
ax3.set_title('Bulk Density')
ax3.set_xlabel('(kg/m3)')

img4 = ax4.imshow(
    pp_primaries.T,
    cmap="bwr_r",
    aspect="auto",
    extent=[offsets[0], offsets[-1], result1.time[-1], result1.time[0]],
    vmin=-clip,
    vmax=clip,
    interpolation="bilinear",
)
ax4.set_xlabel("Offset (m)")
ax4.set_title("Primaries Only")

img5 = ax5.imshow(
    psv_primaries.T,
    cmap="bwr_r",
    aspect="auto",
    extent=[offsets[0], offsets[-1], result2.time[-1], result2.time[0]],
    vmin=-clip,
    vmax=clip,
    interpolation="bilinear",
)
ax5.set_xlabel("Offset (m)")
ax5.set_title("Primaries +\nMode Conversions")

img6 = ax6.imshow(
    pp_multiples_fsonly.T,
    cmap="bwr_r",
    aspect="auto",
    extent=[offsets[0], offsets[-1], result1.time[-1], result1.time[0]],
    vmin=-clip,
    vmax=clip,
    interpolation="bilinear",
)
ax6.set_xlabel("Offset (m)")
ax6.set_title("Primaries +\nFree Surf. Multiples")


img7 = ax7.imshow(
    pp_multiples_total.T,
    cmap="bwr_r",
    aspect="auto",
    extent=[offsets[0], offsets[-1], result1.time[-1], result1.time[0]],
    vmin=-clip,
    vmax=clip,
    interpolation="bilinear",
)
ax7.set_xlabel("Offset (m)")
ax7.set_title("Primaries +\nFree Surf. Multiples +\nInternal Multiples")


img8 = ax8.imshow(
    psv_multiples_total.T,
    cmap="bwr_r",
    aspect="auto",
    extent=[offsets[0], offsets[-1], result2.time[-1], result2.time[0]],
    vmin=-clip,
    vmax=clip,
    interpolation="bilinear",
)
ax8.set_xlabel("Offset (m)")
ax8.set_title("Primaries +\nMode Conversions +\nFree Surf. Multiples +\nInternal Multiples")

img9 = ax9.imshow(
    surface_multiples.T,
    cmap="bwr_r",
    aspect="auto",
    extent=[offsets[0], offsets[-1], result1.time[-1], result1.time[0]],
    vmin=-clip,
    vmax=clip,
    interpolation="bilinear",
)
ax9.set_xlabel("Offset (m)")
ax9.set_title("Free Surface Multiples\n(P-P Reflections Only)")

img10 = ax10.imshow(
    internal_multiples.T,
    cmap="bwr_r",
    aspect="auto",
    extent=[offsets[0], offsets[-1], result1.time[-1], result1.time[0]],
    vmin=-clip,
    vmax=clip,
    interpolation="bilinear",
)
ax10.set_xlabel("Offset (m)")
ax10.set_title("Internal Multiples\n(P-P Reflections Only)")

aura.plot.format_log_axes(fig.get_axes(), 'TWT (sec)')
fig.tight_layout()
fig.subplots_adjust(top=0.90)


curs = MultiCursor(fig.canvas, fig.get_axes(), horizOn=True, vertOn=False,
                   lw=1, c='k')
plt.show()

def plot_avo_curve(gather):
    p = plt.ginput(1)
    tidx = int(round(p[0][1]/dt))
    avo = gather[:, tidx]
    
    fig2 = plt.figure(num=2)
    fig2.clf()
    f2a1 = fig2.add_subplot(111)
    f2a1.plot(offsets, avo)
    f2a1.set_xlabel('Offset (m)')
    f2a1.set_ylabel('Amplitude')
    f2a1.set_title('AVO Response Curve')
    plt.show()
    
    
    