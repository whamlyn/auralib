"""
auralib module for creating synthetic seismograms.

Author:   Wes Hamlyn
Created:  16-Aug-2016
Last Mod: 17-Aug-2016

"""

import numpy as np
from scipy.interpolate import interp1d
from copy import deepcopy


def d2t(depth, data, td_depth, td_twt, target_twt, kind='linear', fill_value='extrapolate'):
    """
    Convenience function for converting a well log from depth to time given
    time-depth pairs.  This function both converts from depth to time and
    resamples the log data to a regular sampling rate in time.
    """

    f1 = interp1d(td_depth, td_twt, kind='linear', bounds_error=False, fill_value=fill_value)
    twt_logstep = f1(depth)

    f2 = interp1d(twt_logstep, data, kind=kind, bounds_error=False, fill_value=fill_value)
    data_twt = f2(target_twt)

    return data_twt



def t2d(twt, data, td_depth, td_twt, target_depth, kind='linear', fill_value='extrapolate'):
    """
    Convenience function for converting a well log from depth to time given
    time-depth pairs.  This function both converts from time to depth and
    resamples the log data to a regular sampling rate in depth.
    """

    f1 = interp1d(td_twt, td_depth, kind='linear', bounds_error=False, fill_value=fill_value)
    depth_timestep = f1(twt)

    f2 = interp1d(depth_timestep, data, kind=kind, bounds_error=False, fill_value=fill_value)
    data_depth = f2(target_depth)

    return data_depth



def backus_length(f_dom, Vs0, mode=1):
    """
    Returns the backus average length for scattering and transmission limits.
    
    L = backus_length(f_dom, Vs0, mode=1)
    
    f_dom = dominant frequency
    Vs0   = minimum shear velocity after upscaling
    mode  = 1 returns the scattering limit; 2 returns the transmission limit
    
    L = Backus averaging length

    """

    Vs0 = float(Vs0)
    f_dom = float(f_dom)

    Ls = Vs0/(3.0*f_dom)
    Lt = 2.0*Vs0/f_dom

    if mode == 1:
        L = Ls

    if mode == 2:
        L = Lt

    return L


def backus_average_v2(vp, vs, rho, nsamp, output='basic'):
    """
    Computes the backus average from vp, vs, and density logs.

    Vpv_b, Vph_b, Vsv_b, Vsh_b, rho_b = backus_average(vp, vs, rho, nsamp, output='basic')
    
    output: 'basic', 'thomsen', or 'all'
                  basic: output vertically propagating Vp, Vs, Rho logs
                thomsen: output basic logs + Thomsen epsilon, delta, and gamma
                    all: output vertical and horizontal propagation velocities plus Thompsen parameters
    """

    # make sure number of samples is odd
    nsamp = np.round(nsamp)
    if nsamp == 0:
        nsamp = 1
    if nsamp%2.0 == 0:
        nsamp = int(np.ceil(nsamp))

    print('Backus averaging over %i samples' % nsamp)

    filt = np.ones(nsamp) / nsamp

    # pad top and bottom of input logs to avoid filtering edge effects.  Do
    # this by duplicating the first and last sample in the vp, vs, and rho logs
    # by ~1/2 the Backus filtering length
    pad_length = int((nsamp)/2)

    vp = np.hstack([np.ones(pad_length)*vp[0], vp, np.ones(pad_length)*vp[-1]])
    vs = np.hstack([np.ones(pad_length)*vs[0], vs, np.ones(pad_length)*vs[-1]])
    rho = np.hstack([np.ones(pad_length)*rho[0], rho, np.ones(pad_length)*rho[-1]])

    mu = rho*(vs**2.0)
    K = rho*(vp**2.0) - 4.0/3.0*mu
    lame = K - 2.0/3.0*mu


    A1 = (4.0*mu*(lame+mu))/(lame+2.0*mu)
    A2 = 1.0/(lame+2.0*mu)
    A3 = lame/(lame+2.0*mu)

    A1 = np.convolve(filt, A1, mode='same')
    A2 = np.convolve(filt, A2, mode='same')
    A3 = np.convolve(filt, A3, mode='same')
    A = A1 + A2**-1 * A3**2

    B1 = (2.0*mu*lame)/(lame+2.0*mu)
    B2 = 1.0/(lame+2.0*mu)
    B3 = lame/(lame+2.0*mu)

    B1 = np.convolve(filt, B1, mode='same')
    B2 = np.convolve(filt, B2, mode='same')
    B3 = np.convolve(filt, B3, mode='same')
    B = B1 + B2**-1 * B3**2

    C = 1.0/(lame+2.0*mu)
    C = (np.convolve(filt, C, mode='same'))**-1

    F1 = 1.0/(lame+2.0*mu)
    F2 = lame/(lame+2.0*mu)

    F1 = np.convolve(filt, F1, mode='same')
    F2 = np.convolve(filt, F2, mode='same')
    F = F1**-1 * F2

    D = 1.0/mu
    D = (np.convolve(filt, D, mode='same'))**-1

    M = np.convolve(filt, mu, mode='same')

    rho_b = np.convolve(filt, rho, mode='same')

    Vpv_b = np.sqrt(C/rho_b) # vertically propagating p-wave
    Vph_b = np.sqrt(A/rho_b) # horizontally propagating p-wave
    Vshh_b = np.sqrt(M/rho_b) # horizontally polarized horizontally propagating s-wave
    Vsvh_b = np.sqrt(D/rho_b) # vertically polarized horizontally propagating s-wave
    Vshv_b = np.sqrt(D/rho_b) # any polarization vertically propagating s-wave
    Vsvv_b = np.sqrt(D/rho_b) # any polarization vertically propagating s-wave
    
    # calculate thompsen parameters
    c11 = A
    c33 = C
    c44 = D
    c66 = M
    c13 = F
    epsilon = (c11-c33)/(2*c33)
    gamma = (c66-c44)/(2*c44)
    delta = ((c13+c44)**2 - (c33-c44)**2)/(2*c33*(c33-c44))

    # remove padding added to avoid edge effects
    Vpv_b = Vpv_b[pad_length:-pad_length]
    Vph_b = Vph_b[pad_length:-pad_length]
    Vshh_b = Vshh_b[pad_length:-pad_length]
    Vsvh_b = Vsvh_b[pad_length:-pad_length]
    Vshv_b = Vshv_b[pad_length:-pad_length]
    Vsvv_b = Vsvv_b[pad_length:-pad_length]
    rho_b = rho_b[pad_length:-pad_length]    
    epsilon = epsilon[pad_length:-pad_length]
    gamma = gamma[pad_length:-pad_length]
    delta = delta[pad_length:-pad_length]
    
    if output == 'basic':
        return Vpv_b, Vsvv_b, rho_b
    
    elif output == 'thomsen':
        return Vpv_b, Vsvv_b, rho_b, epsilon, delta, gamma
    
    else:
        return Vpv_b, Vph_b, Vsvv_b, Vshh_b, rho_b, epsilon, delta, gamma


def plot_wigva(ax, tdata, twt, excursion=1.0, repeat=3, c='k', lw=0.5, peak=True,
               trough=True, pcolor='blue', tcolor='red', alpha=0.3):
    """
    Plot a single seismic trace in wiggle-variable area format.
    """

    tdata_norm = tdata / np.max(np.abs(tdata))
    for i in range(repeat):
        baseline = 1.0*i
        tdata_plt = tdata_norm*excursion + baseline

        if peak==True:
            ax.fill_betweenx(twt, baseline, tdata_plt, where=tdata_plt>=baseline,
                                fc=pcolor, alpha=alpha, interpolate=True)
        if trough==True:
            ax.fill_betweenx(twt, baseline, tdata_plt, where=tdata_plt<=baseline,
                                fc=tcolor, alpha=alpha, interpolate=True)
        ax.plot(tdata_plt, twt, c=c, lw=lw)


#  Multiple Modelling Functions

import numpy as np
from scipy.signal import fftconvolve


def next_power_of_two(n):
    """Return the smallest power of two greater than or equal to n."""
    return 1 << (int(n) - 1).bit_length()


def prepare_layer_model(depth, vp, rho):
    """
    Prepare a 1D acoustic layered model.

    Parameters
    ----------
    depth : array_like
        Depth to the top of each layer. Must increase monotonically.
        The last layer is treated as a half-space.
    vp : array_like
        P-wave velocity in each layer.
    rho : array_like
        Density in each layer.

    Returns
    -------
    dz : ndarray
        Thickness of each finite layer.
    impedance : ndarray
        Acoustic impedance of each layer.
    reflection_coefficients : ndarray
        Normal-incidence pressure reflection coefficient at each interface.
    transmission_product : ndarray
        Product t_down * t_up at each interface.
    """
    depth = np.asarray(depth, dtype=float)
    vp = np.asarray(vp, dtype=float)
    rho = np.asarray(rho, dtype=float)

    if not (depth.ndim == vp.ndim == rho.ndim == 1):
        raise ValueError("depth, vp and rho must be one-dimensional.")

    if not (depth.size == vp.size == rho.size):
        raise ValueError("depth, vp and rho must have the same length.")

    if depth.size < 2:
        raise ValueError("At least two layers are required.")

    if np.any(np.diff(depth) <= 0):
        raise ValueError("depth must be strictly increasing.")

    if np.any(vp <= 0):
        raise ValueError("All velocities must be positive.")

    if np.any(rho <= 0):
        raise ValueError("All densities must be positive.")

    # Finite-layer thicknesses. The final layer is a half-space.
    dz = np.diff(depth)

    impedance = vp * rho

    z1 = impedance[:-1]
    z2 = impedance[1:]

    # Reflection coefficient for incidence from above.
    reflection_coefficients = (z2 - z1) / (z2 + z1)

    # Pressure transmission coefficients:
    #
    # t_down = 2 * Z2 / (Z1 + Z2)
    # t_up   = 2 * Z1 / (Z1 + Z2)
    #
    # Their product equals 1 - r^2.
    transmission_product = 1.0 - reflection_coefficients**2

    return dz, impedance, reflection_coefficients, transmission_product


def primary_reflection_response(
    frequencies,
    depth,
    vp,
    rho,
    include_transmission=True,
):
    """
    Calculate the frequency-domain primary-only reflection response.

    This contains one reflection from each interface but no internal
    or surface-related multiples.

    Parameters
    ----------
    frequencies : ndarray
        Non-negative FFT frequencies in Hz.
    depth, vp, rho : array_like
        Layered acoustic model.
    include_transmission : bool, default True
        Include transmission losses through overlying interfaces.

    Returns
    -------
    response : complex ndarray
        Primary-only reflection response at each frequency.
    """
    frequencies = np.asarray(frequencies, dtype=float)

    dz, _, r, transmission_product = prepare_layer_model(
        depth, vp, rho
    )

    vp = np.asarray(vp, dtype=float)

    # One-way time through each finite layer.
    layer_time = dz / vp[:-1]

    # Two-way time from the top of the model to each interface.
    interface_twt = 2.0 * np.cumsum(layer_time)

    if include_transmission:
        amplitudes = np.empty_like(r)

        cumulative_transmission = 1.0

        for i in range(r.size):
            amplitudes[i] = cumulative_transmission * r[i]
            cumulative_transmission *= transmission_product[i]
    else:
        amplitudes = r.copy()

    omega = 2.0 * np.pi * frequencies

    # exp(-i * omega * time) represents a causal delay under NumPy's
    # inverse FFT convention.
    phase = np.exp(
        -1j * omega[:, None] * interface_twt[None, :]
    )

    response = phase @ amplitudes

    return response


def full_internal_multiple_response(
    frequencies,
    depth,
    vp,
    rho,
):
    """
    Calculate the exact 1D acoustic reflection response containing
    primaries and all orders of internal multiples.

    Surface-related multiples are not included.

    Parameters
    ----------
    frequencies : ndarray
        Non-negative FFT frequencies in Hz.
    depth, vp, rho : array_like
        Layered acoustic model.

    Returns
    -------
    response : complex ndarray
        Reflection response containing primaries and internal multiples.
    """
    frequencies = np.asarray(frequencies, dtype=float)

    dz, _, r, _ = prepare_layer_model(depth, vp, rho)

    vp = np.asarray(vp, dtype=float)
    omega = 2.0 * np.pi * frequencies

    # Begin with the deepest interface.
    effective_response = np.full(
        frequencies.shape,
        r[-1],
        dtype=complex,
    )

    # Work upward through the remaining interfaces.
    #
    # Interface i is separated from interface i+1 by layer i+1.
    for i in range(r.size - 2, -1, -1):
        layer_index = i + 1

        two_way_layer_time = (
            2.0 * dz[layer_index] / vp[layer_index]
        )

        propagation = np.exp(
            -1j * omega * two_way_layer_time
        )

        deeper_response = effective_response * propagation

        denominator = 1.0 + r[i] * deeper_response

        effective_response = (
            r[i] + deeper_response
        ) / denominator

    # Propagate from the top of layer 0 down to the first interface
    # and back.
    first_interface_twt = 2.0 * dz[0] / vp[0]

    top_propagation = np.exp(
        -1j * omega * first_interface_twt
    )

    return effective_response * top_propagation


def add_surface_multiples(
    subsurface_response,
    frequencies,
    free_surface_coefficient=-1.0,
    surface_delay=0.0,
):
    """
    Add all orders of surface-related multiples using a feedback series.

    Parameters
    ----------
    subsurface_response : complex ndarray
        Frequency-domain subsurface reflection response.
    frequencies : ndarray
        Non-negative FFT frequencies in Hz.
    free_surface_coefficient : float, default -1
        Pressure reflection coefficient at the free surface.
    surface_delay : float, default 0
        Additional two-way delay, in seconds, between the reference
        datum and the free surface. Normally zero for an idealized
        surface source and receiver.

    Returns
    -------
    response : complex ndarray
        Response containing the original subsurface events and all
        orders of surface-related multiples.
    """
    subsurface_response = np.asarray(
        subsurface_response,
        dtype=complex,
    )

    frequencies = np.asarray(frequencies, dtype=float)
    omega = 2.0 * np.pi * frequencies

    surface_operator = (
        free_surface_coefficient
        * np.exp(-1j * omega * surface_delay)
    )

    denominator = (
        1.0 - surface_operator * subsurface_response
    )

    return subsurface_response / denominator


def time_samples_to_depth(time, depth, vp):
    """
    Convert zero-offset two-way-time samples to depth using the
    layered interval-velocity model.

    Parameters
    ----------
    time : array_like
        Two-way-time samples in seconds.

    depth : array_like
        Depth to the top of each layer.

    vp : array_like
        P-wave interval velocity for each layer. The final layer is
        treated as a half-space.

    Returns
    -------
    sample_depth : ndarray
        Depth corresponding to each two-way-time sample.
    """
    time = np.asarray(time, dtype=float).reshape(-1)
    depth = np.asarray(depth, dtype=float).reshape(-1)
    vp = np.asarray(vp, dtype=float).reshape(-1)

    if depth.size != vp.size:
        raise ValueError(
            "depth and vp must contain the same number of values."
        )

    if np.any(np.diff(depth) <= 0):
        raise ValueError("depth must be strictly increasing.")

    if np.any(vp <= 0):
        raise ValueError("vp must contain only positive values.")

    # One-way time associated with each synthetic sample.
    one_way_time = time / 2.0

    # One-way time through each finite layer.
    layer_thickness = np.diff(depth)

    layer_time = (
        layer_thickness / vp[:-1]
    )

    # One-way time to the top of every layer.
    layer_top_time = np.concatenate((
        [0.0],
        np.cumsum(layer_time),
    ))

    # Identify the layer containing each time sample.
    layer_index = (
        np.searchsorted(
            layer_top_time,
            one_way_time,
            side="right",
        )
        - 1
    )

    # Any samples later than the finite model enter the final
    # half-space.
    layer_index = np.clip(
        layer_index,
        0,
        depth.size - 1,
    )

    time_inside_layer = (
        one_way_time
        - layer_top_time[layer_index]
    )

    sample_depth = (
        depth[layer_index]
        + vp[layer_index] * time_inside_layer
    )

    return sample_depth



def model_zero_offset_synthetic(
    depth,
    vp,
    rho,
    wavelet,
    dt=0.001,
    record_length=3.0,
    free_surface_coefficient=-1.0,
    include_transmission=True,
    nfft=None,
    primaries=True,
    surface_multiples=False,
    internal_multiples=False,
):
    """
    Model a zero-offset 1D acoustic synthetic seismic trace.

    Parameters
    ----------
    depth : array_like
        Depth to the top of each acoustic layer. The final layer is
        treated as a half-space.

    vp : array_like
        P-wave velocity for each layer.

    rho : array_like
        Density for each layer.

    wavelet : array_like
        Seismic wavelet. It should normally be an odd-length,
        zero-phase wavelet centred on its middle sample.

    dt : float, default 0.001
        Output sample interval in seconds.

    record_length : float, default 3.0
        Output record length in seconds.

    free_surface_coefficient : float, default -1.0
        Pressure reflection coefficient at the free surface.

    include_transmission : bool, default True
        Include transmission losses in the primary-only response.

    nfft : int or None, default None
        FFT size. If None, a zero-padded FFT size is selected
        automatically.

    primaries : bool, default True
        Include primary reflections in the returned trace.

    surface_multiples : bool, default False
        Include all orders of surface-related multiples.

    internal_multiples : bool, default False
        Include all orders of 1D internal multiples.

    Returns
    -------
    time : ndarray
        Output time axis in seconds.

    synthetic : ndarray
        Wavelet-convolved synthetic trace.

    impulse_response : ndarray
        Selected earth impulse response before wavelet convolution.

    frequency_response : complex ndarray
        Selected frequency-domain earth response.

    Notes
    -----
    The Boolean switches select individual response components:

    primaries=True
        Retains the primary reflection response.

    internal_multiples=True
        Adds the difference between the complete internal-multiple
        response and the primary-only response.

    surface_multiples=True
        Adds the surface-multiple component generated from the selected
        subsurface model.

    Consequently, setting primaries=False allows the multiple energy
    to be returned without the primary reflections.
    """
    depth = np.asarray(depth, dtype=float).reshape(-1)
    vp = np.asarray(vp, dtype=float).reshape(-1)
    rho = np.asarray(rho, dtype=float).reshape(-1)
    wavelet = np.asarray(wavelet, dtype=float).reshape(-1)

    if dt <= 0:
        raise ValueError("dt must be positive.")

    if record_length <= 0:
        raise ValueError("record_length must be positive.")

    if wavelet.size == 0:
        raise ValueError("wavelet cannot be empty.")

    nsamples = int(np.floor(record_length / dt)) + 1

    # Additional padding reduces circular wrap-around from late
    # multiple energy.
    if nfft is None:
        nfft = next_power_of_two(4 * nsamples)

    if nfft < nsamples:
        raise ValueError(
            "nfft must be greater than or equal to the number "
            "of output samples."
        )

    frequencies = np.fft.rfftfreq(nfft, d=dt)

    # Primary-only subsurface response.
    primary_response = primary_reflection_response(
        frequencies=frequencies,
        depth=depth,
        vp=vp,
        rho=rho,
        include_transmission=include_transmission,
    )

    # Calculate the complete subsurface response only when needed.
    if internal_multiples:
        subsurface_response = full_internal_multiple_response(
            frequencies=frequencies,
            depth=depth,
            vp=vp,
            rho=rho,
        )
    else:
        subsurface_response = primary_response

    # Start with an empty response and add the requested components.
    frequency_response = np.zeros_like(
        primary_response,
        dtype=complex,
    )

    if primaries:
        frequency_response += primary_response

    if internal_multiples:
        # Remove the primaries from the complete internal response,
        # leaving only internal-multiple energy.
        internal_component = (
            subsurface_response - primary_response
        )

        frequency_response += internal_component

    if surface_multiples:
        response_with_surface_multiples = add_surface_multiples(
            subsurface_response=subsurface_response,
            frequencies=frequencies,
            free_surface_coefficient=(
                free_surface_coefficient
            ),
        )

        # Retain only the additional energy generated by the
        # free-surface feedback.
        surface_component = (
            response_with_surface_multiples
            - subsurface_response
        )

        frequency_response += surface_component

    # Convert the selected frequency response into a causal earth
    # impulse response.
    full_impulse_response = np.fft.irfft(
        frequency_response,
        n=nfft,
    )

    impulse_response = full_impulse_response[:nsamples]

    # Convolve the earth response with the supplied wavelet.
    synthetic = fftconvolve(
        impulse_response,
        wavelet,
        mode="same",
    )

    # Two-way-time reference for every synthetic sample.
    time = np.arange(nsamples) * dt

    # Depth reference for every synthetic sample.
    sample_depth = time_samples_to_depth(
        time=time,
        depth=depth,
        vp=vp,
    )

    return (
        time,
        sample_depth,
        synthetic,
        impulse_response,
        frequency_response,
    )




