"""Kennett-style P-SV modelling for a 1-D isotropic layered earth.

This reference implementation builds stable reflection/transmission scattering
matrices, cascades them with a Redheffer star product, and synthesizes a 2-D
line-source CMP gather by numerical horizontal-slowness integration.

Conventions
-----------
* x is horizontal and z is positive downward.
* ``depth[i]`` is the top of layer i; the last layer is a half-space.
* Fourier synthesis uses exp(-i omega t) time dependence and NumPy irfft.
* Outputs may be modal P-P, modal P-SV, vertical particle velocity, or inline
  horizontal particle velocity.

The model is isotropic and lossless. Complex vertical slownesses retain
evanescent waves. SI units are recommended, although any internally consistent
depth/velocity/density system may be used.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Literal, Optional, Union

import numpy as np


Array = np.ndarray


@dataclass(frozen=True)
class KennettResult:
    """Results returned by :func:`model_kennett_cmp`."""

    time: Array
    offsets: Array
    gather: Union[Array, Dict[str, Array]]
    impulse_gather: Union[Array, Dict[str, Array]]
    frequency_gather: Union[Array, Dict[str, Array]]
    frequencies: Array
    slowness: Array
    uncorrected_gather: Union[Array, Dict[str, Array]]
    nmo_velocity: Optional[Array]
    receiver_component: str
    gathers: Dict[str, Array]
    earth_response: "KennettEarthResponse"


COMPONENTS = (
    "modal_pp",
    "modal_ps",
    "vertical_pp",
    "horizontal_pp",
    "vertical_ps",
    "horizontal_ps",
    "vertical_velocity",
    "horizontal_velocity",
)


@dataclass(frozen=True)
class KennettPrecomputed:
    """Frequency-independent quantities for one earth and p sampling."""

    depth: Array
    vp: Array
    vs: Array
    rho: Array
    slowness: Array
    thickness: Array
    qp: Array
    qs: Array
    down_modes: Array
    up_modes: Array
    interfaces: Array
    free_surface: Array


@dataclass(frozen=True)
class KennettEarthResponse:
    """Reusable p-omega response prior to offset and wavelet synthesis."""

    model: KennettPrecomputed
    frequencies: Array
    active_indices: Array
    responses: Dict[str, Array]
    internal_multiples: bool
    surface_multiples: bool
    surface_multiple_order: Optional[int]
    primaries: bool
    complex_frequency_damping: float


def rms_velocity_from_layers(time, depth, vp) -> Array:
    """Calculate vertical P-wave RMS velocity at zero-offset TWT samples.

    Times beyond the last finite layer are continued through the bottom
    half-space using its P-wave velocity.
    """
    time = np.asarray(time, dtype=float).reshape(-1)
    depth = np.asarray(depth, dtype=float).reshape(-1)
    vp = np.asarray(vp, dtype=float).reshape(-1)
    if depth.size != vp.size or depth.size < 2:
        raise ValueError("depth and vp must have equal lengths of at least two")
    if np.any(time < 0) or np.any(np.diff(time) < 0):
        raise ValueError("time must be nonnegative and increasing")

    finite_dt = 2.0 * np.diff(depth) / vp[:-1]
    layer_top_twt = np.concatenate(([0.0], np.cumsum(finite_dt)))
    vrms2_integral = np.zeros_like(time)

    # Integral of V_interval^2 dt from zero to each requested TWT.
    for i, velocity in enumerate(vp):
        top = layer_top_twt[i]
        if i < vp.size - 1:
            base = layer_top_twt[i + 1]
            duration = np.clip(time - top, 0.0, base - top)
        else:
            duration = np.maximum(time - top, 0.0)
        vrms2_integral += velocity**2 * duration

    vrms = np.empty_like(time)
    positive = time > 0.0
    vrms[positive] = np.sqrt(vrms2_integral[positive] / time[positive])
    vrms[~positive] = vp[0]
    return vrms


def apply_nmo_correction(
    gather,
    time,
    offsets,
    nmo_velocity,
    *,
    stretch_mute: Optional[float] = None,
    fill_value: float = 0.0,
) -> Array:
    """Apply conventional hyperbolic NMO correction to an offset gather.

    The output sample at zero-offset time ``t0`` is interpolated from

        t = sqrt(t0**2 + offset**2 / Vnmo(t0)**2).

    Parameters
    ----------
    gather : array_like, shape (n_offsets, n_time)
        Uncorrected offset-time gather.
    time : array_like, shape (n_time,)
        Two-way-time axis in seconds.
    offsets : array_like, shape (n_offsets,)
        Full source-receiver offsets, in units consistent with velocity.
    nmo_velocity : float or array_like, shape (n_time,)
        Scalar or time-varying NMO velocity.
    stretch_mute : float, optional
        Mute samples where ``(t - t0) / t0`` exceeds this fraction. For
        example, 0.30 mutes samples with more than 30 percent NMO stretch.
    fill_value : float
        Value assigned outside the input time range and within the mute.
    """
    gather = np.asarray(gather, dtype=float)
    time = np.asarray(time, dtype=float).reshape(-1)
    offsets = np.asarray(offsets, dtype=float).reshape(-1)
    velocity = np.asarray(nmo_velocity, dtype=float)
    if gather.shape != (offsets.size, time.size):
        raise ValueError("gather must have shape (len(offsets), len(time))")
    if time.size < 2 or np.any(np.diff(time) <= 0):
        raise ValueError("time must contain at least two increasing samples")
    if velocity.ndim == 0:
        velocity = np.full(time.size, float(velocity))
    else:
        velocity = velocity.reshape(-1)
    if velocity.size != time.size:
        raise ValueError("nmo_velocity must be scalar or have len(time) samples")
    if np.any(~np.isfinite(velocity)) or np.any(velocity <= 0):
        raise ValueError("nmo_velocity must be finite and positive")
    if stretch_mute is not None and stretch_mute < 0:
        raise ValueError("stretch_mute must be nonnegative")

    corrected = np.full_like(gather, fill_value, dtype=float)
    for ix, offset in enumerate(offsets):
        input_time = np.sqrt(time**2 + (offset / velocity) ** 2)
        corrected[ix] = np.interp(
            input_time, time, gather[ix], left=fill_value, right=fill_value
        )
        if stretch_mute is not None:
            stretch = np.full_like(time, np.inf)
            positive = time > 0.0
            stretch[positive] = (
                input_time[positive] - time[positive]
            ) / time[positive]
            corrected[ix, stretch > stretch_mute] = fill_value
    return corrected


def next_power_of_two(n: int) -> int:
    """Smallest power of two greater than or equal to *n*."""
    if n < 1:
        raise ValueError("n must be positive")
    return 1 << (int(n) - 1).bit_length()


def _validate_model(depth, vp, vs, rho):
    arrays = [np.asarray(x, dtype=float).reshape(-1) for x in (depth, vp, vs, rho)]
    depth, vp, vs, rho = arrays
    if not (depth.size == vp.size == vs.size == rho.size):
        raise ValueError("depth, vp, vs and rho must have equal lengths")
    if depth.size < 2:
        raise ValueError("at least two layers are required")
    if np.any(~np.isfinite(depth)) or np.any(np.diff(depth) <= 0):
        raise ValueError("depth must be finite and strictly increasing")
    if any(np.any(~np.isfinite(a)) or np.any(a <= 0) for a in (vp, vs, rho)):
        raise ValueError("vp, vs and rho must be finite and positive")
    if np.any(vs >= vp):
        raise ValueError("this isotropic implementation requires vs < vp")
    return depth, vp, vs, rho


def vertical_slownesses(p: float, vp: float, vs: float) -> tuple[complex, complex]:
    """Return causal P and SV vertical slownesses for horizontal slowness p."""
    # scimath.sqrt returns +i*sqrt(abs(x)) for negative real arguments. With
    # exp(i omega q h), this makes evanescent waves decay with depth.
    qp = np.lib.scimath.sqrt(1.0 / vp**2 - p**2)
    qs = np.lib.scimath.sqrt(1.0 / vs**2 - p**2)
    if np.imag(qp) < 0:
        qp = -qp
    if np.imag(qs) < 0:
        qs = -qs
    return complex(qp), complex(qs)


def modal_matrices(p: float, vp: float, vs: float, rho: float) -> tuple[Array, Array]:
    """State-vector columns for downgoing and upgoing P/SV modes.

    Rows are ``[u_x, u_z, sigma_xz/(i*omega), sigma_zz/(i*omega)]``.
    Columns are P and SV. Modal displacement vectors have unit length.
    """
    qp, qs = vertical_slownesses(p, vp, vs)
    mu = rho * vs**2
    lam = rho * vp**2 - 2.0 * mu

    def state(q: complex, ux: complex, uz: complex) -> Array:
        txz = mu * (q * ux + p * uz)
        tzz = lam * (p * ux + q * uz) + 2.0 * mu * q * uz
        return np.array([ux, uz, txz, tzz], dtype=complex)

    # SV polarization is chosen perpendicular to its propagation-slowness
    # vector. Its arbitrary overall sign is kept consistent at all interfaces.
    down = np.column_stack((
        state(qp, p * vp, qp * vp),
        state(qs, -qs * vs, p * vs),
    ))
    up = np.column_stack((
        state(-qp, p * vp, -qp * vp),
        state(-qs, qs * vs, p * vs),
    ))
    return down, up


def interface_scattering(
    p: float,
    vp_top: float,
    vs_top: float,
    rho_top: float,
    vp_bottom: float,
    vs_bottom: float,
    rho_bottom: float,
) -> Array:
    """Return a 4x4 elastic interface scattering matrix.

    It maps incoming ``[D_top, U_bottom]`` to outgoing
    ``[U_top, D_bottom]``; each entry is a two-component [P, SV] vector.
    """
    dt, ut = modal_matrices(p, vp_top, vs_top, rho_top)
    db, ub = modal_matrices(p, vp_bottom, vs_bottom, rho_bottom)
    lhs = np.column_stack((ut, -db))
    rhs = np.column_stack((-dt, ub))
    try:
        return np.linalg.solve(lhs, rhs)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(lhs) @ rhs


def propagation_scattering(p: float, omega: float, h: float, vp: float, vs: float) -> Array:
    """Scattering matrix for propagation through one uniform layer."""
    qp, qs = vertical_slownesses(p, vp, vs)
    # A positive travel-time delay has spectrum exp(-i*omega*t) under
    # NumPy's FFT convention.
    e = np.diag(np.exp(-1j * omega * np.array([qp, qs]) * h))
    z = np.zeros((2, 2), dtype=complex)
    return np.block([[z, e], [e, z]])


def star_product(a: Array, b: Array) -> Array:
    """Cascade scattering networks *a* (above) and *b* (below)."""
    a11, a12, a21, a22 = a[:2, :2], a[:2, 2:], a[2:, :2], a[2:, 2:]
    b11, b12, b21, b22 = b[:2, :2], b[:2, 2:], b[2:, :2], b[2:, 2:]
    eye = np.eye(2, dtype=complex)
    # Solves are used instead of explicit inverses.
    x = np.linalg.solve(eye - b11 @ a22, b11 @ a21)
    y = np.linalg.solve(eye - b11 @ a22, b12)
    s11 = a11 + a12 @ x
    s12 = a12 @ y
    s21 = b21 @ np.linalg.solve(eye - a22 @ b11, a21)
    s22 = b22 + b21 @ a22 @ y
    return np.block([[s11, s12], [s21, s22]])


def _cascade_no_internal(a: Array, b: Array) -> Array:
    """Cascade with inter-network reverberation terms suppressed."""
    a11, a12, a21, a22 = a[:2, :2], a[:2, 2:], a[2:, :2], a[2:, 2:]
    b11, b12, b21, b22 = b[:2, :2], b[:2, 2:], b[2:, :2], b[2:, 2:]
    s11 = a11 + a12 @ b11 @ a21
    s12 = a12 @ b12
    s21 = b21 @ a21
    s22 = b22 + b21 @ a22 @ b12
    return np.block([[s11, s12], [s21, s22]])


def stack_scattering(
    p: float,
    omega: float,
    depth,
    vp,
    vs,
    rho,
    internal_multiples: bool = True,
) -> Array:
    """Scattering matrix from the top of layer 0 into the bottom half-space."""
    depth, vp, vs, rho = _validate_model(depth, vp, vs, rho)
    cascade = star_product if internal_multiples else _cascade_no_internal
    total = propagation_scattering(p, omega, depth[1] - depth[0], vp[0], vs[0])
    for i in range(depth.size - 1):
        total = cascade(total, interface_scattering(
            p, vp[i], vs[i], rho[i], vp[i + 1], vs[i + 1], rho[i + 1]
        ))
        if i + 1 < depth.size - 1:
            total = cascade(total, propagation_scattering(
                p, omega, depth[i + 2] - depth[i + 1], vp[i + 1], vs[i + 1]
            ))
    return total


def free_surface_reflection(p: float, vp: float, vs: float, rho: float) -> Array:
    """Map an incident upgoing [P, SV] wave to reflected downgoing modes."""
    down, up = modal_matrices(p, vp, vs, rho)
    td, tu = down[2:, :], up[2:, :]
    try:
        return np.linalg.solve(td, -tu)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(td) @ (-tu)


def precompute_kennett_model(depth, vp, vs, rho, slowness) -> KennettPrecomputed:
    """Precompute every frequency-independent model quantity.

    This calculation is reusable for all frequencies, offsets, wavelets and
    receiver components that use the same layered model and p sampling.
    """
    depth, vp, vs, rho = _validate_model(depth, vp, vs, rho)
    p = np.asarray(slowness, dtype=float).reshape(-1)
    if p.size < 3 or np.any(~np.isfinite(p)) or np.any(np.diff(p) <= 0):
        raise ValueError("slowness must contain at least three increasing values")

    np_ = p.size
    nl = depth.size
    qp = np.empty((np_, nl), dtype=complex)
    qs = np.empty((np_, nl), dtype=complex)
    down = np.empty((np_, nl, 4, 2), dtype=complex)
    up = np.empty_like(down)

    # This loop runs once per earth model, rather than once per frequency.
    for jp, pj in enumerate(p):
        for il in range(nl):
            qp[jp, il], qs[jp, il] = vertical_slownesses(
                pj, vp[il], vs[il]
            )
            down[jp, il], up[jp, il] = modal_matrices(
                pj, vp[il], vs[il], rho[il]
            )

    interfaces = np.empty((np_, nl - 1, 4, 4), dtype=complex)
    for il in range(nl - 1):
        lhs = np.empty((np_, 4, 4), dtype=complex)
        rhs = np.empty_like(lhs)
        lhs[:, :, :2] = up[:, il]
        lhs[:, :, 2:] = -down[:, il + 1]
        rhs[:, :, :2] = -down[:, il]
        rhs[:, :, 2:] = up[:, il + 1]
        try:
            interfaces[:, il] = np.linalg.solve(lhs, rhs)
        except np.linalg.LinAlgError:
            interfaces[:, il] = np.array([
                np.linalg.pinv(lhs[j]) @ rhs[j] for j in range(np_)
            ])

    td = down[:, 0, 2:, :]
    tu = up[:, 0, 2:, :]
    try:
        free_surface = np.linalg.solve(td, -tu)
    except np.linalg.LinAlgError:
        free_surface = np.array([
            np.linalg.pinv(td[j]) @ (-tu[j]) for j in range(np_)
        ])

    return KennettPrecomputed(
        depth=depth,
        vp=vp,
        vs=vs,
        rho=rho,
        slowness=p,
        thickness=np.diff(depth),
        qp=qp,
        qs=qs,
        down_modes=down,
        up_modes=up,
        interfaces=interfaces,
        free_surface=free_surface,
    )


def _assemble_scattering(s11, s12, s21, s22) -> Array:
    """Assemble batched 2x2 blocks without the overhead of np.block."""
    shape = s11.shape[:-2] + (4, 4)
    result = np.empty(shape, dtype=complex)
    result[..., :2, :2] = s11
    result[..., :2, 2:] = s12
    result[..., 2:, :2] = s21
    result[..., 2:, 2:] = s22
    return result


def _batched_star_product(a: Array, b: Array) -> Array:
    """Redheffer cascade for all slownesses simultaneously."""
    a11, a12, a21, a22 = a[..., :2, :2], a[..., :2, 2:], a[..., 2:, :2], a[..., 2:, 2:]
    b11, b12, b21, b22 = b[..., :2, :2], b[..., :2, 2:], b[..., 2:, :2], b[..., 2:, 2:]
    eye = np.broadcast_to(np.eye(2, dtype=complex), a11.shape)
    left = eye - b11 @ a22
    x = np.linalg.solve(left, b11 @ a21)
    y = np.linalg.solve(left, b12)
    s11 = a11 + a12 @ x
    s12 = a12 @ y
    s21 = b21 @ np.linalg.solve(eye - a22 @ b11, a21)
    s22 = b22 + b21 @ a22 @ y
    return _assemble_scattering(s11, s12, s21, s22)


def _batched_cascade_no_internal(a: Array, b: Array) -> Array:
    a11, a12, a21, a22 = a[..., :2, :2], a[..., :2, 2:], a[..., 2:, :2], a[..., 2:, 2:]
    b11, b12, b21, b22 = b[..., :2, :2], b[..., :2, 2:], b[..., 2:, :2], b[..., 2:, 2:]
    return _assemble_scattering(
        a11 + a12 @ b11 @ a21,
        a12 @ b12,
        b21 @ a21,
        b22 + b21 @ a22 @ b12,
    )


def _batched_propagation(model: KennettPrecomputed, omega: float, layer: int) -> Array:
    phase = np.exp(
        -1j * omega
        * np.stack((model.qp[:, layer], model.qs[:, layer]), axis=-1)
        * model.thickness[layer]
    )
    n = model.slowness.size
    result = np.zeros((n, 4, 4), dtype=complex)
    result[:, 0, 2] = phase[:, 0]
    result[:, 1, 3] = phase[:, 1]
    result[:, 2, 0] = phase[:, 0]
    result[:, 3, 1] = phase[:, 1]
    return result


def _batched_stack_scattering(
    model: KennettPrecomputed,
    omega: float,
    internal_multiples: bool,
) -> Array:
    cascade = _batched_star_product if internal_multiples else _batched_cascade_no_internal
    total = _batched_propagation(model, omega, 0)
    for il in range(model.depth.size - 1):
        total = cascade(total, model.interfaces[:, il])
        if il + 1 < model.depth.size - 1:
            total = cascade(total, _batched_propagation(model, omega, il + 1))
    return total


def _batched_component_responses(
    model: KennettPrecomputed,
    omega: float,
    internal_multiples: bool,
    surface_multiples: bool,
    surface_multiple_order: Optional[int] = None,
) -> Dict[str, Array]:
    """Calculate all eight receiver outputs from one earth recursion."""
    earth = _batched_stack_scattering(model, omega, internal_multiples)[:, :2, :2]
    incident_p = earth[:, :, 0]
    if surface_multiples:
        if surface_multiple_order is None:
            eye = np.broadcast_to(np.eye(2, dtype=complex), earth.shape)
            up = np.linalg.solve(
                eye - earth @ model.free_surface,
                incident_p[..., None],
            )[..., 0]
        else:
            # Finite Neumann series. Order zero retains the initial subsurface
            # reflection; order one adds the first free-surface reverberation.
            term = incident_p.copy()
            up = term.copy()
            for _ in range(surface_multiple_order):
                down_term = np.einsum(
                    "pij,pj->pi", model.free_surface, term
                )
                term = np.einsum("pij,pj->pi", earth, down_term)
                up += term
        down_at_receiver = np.einsum("pij,pj->pi", model.free_surface, up)
    else:
        up = incident_p
        down_at_receiver = np.zeros_like(up)

    # Decompose receiver motion by the modal type present at the receiver.
    # A P-mode arrival may still contain earlier P-SV-P conversions in its
    # propagation history; this is an arriving-mode, not a full path, label.
    p_displacement = (
        model.up_modes[:, 0, :2, 0] * up[:, 0, None]
        + model.down_modes[:, 0, :2, 0] * down_at_receiver[:, 0, None]
    )
    sv_displacement = (
        model.up_modes[:, 0, :2, 1] * up[:, 1, None]
        + model.down_modes[:, 0, :2, 1] * down_at_receiver[:, 1, None]
    )
    p_velocity = -1j * omega * p_displacement
    sv_velocity = -1j * omega * sv_displacement
    return {
        "modal_pp": up[:, 0],
        "modal_ps": up[:, 1],
        "vertical_pp": p_velocity[:, 1],
        "horizontal_pp": p_velocity[:, 0],
        "vertical_ps": sv_velocity[:, 1],
        "horizontal_ps": sv_velocity[:, 0],
        "vertical_velocity": p_velocity[:, 1] + sv_velocity[:, 1],
        "horizontal_velocity": p_velocity[:, 0] + sv_velocity[:, 0],
    }


def build_kennett_earth_response(
    depth,
    vp,
    vs,
    rho,
    frequencies,
    slowness,
    *,
    max_frequency: Optional[float] = None,
    primaries: bool = True,
    internal_multiples: bool = True,
    surface_multiples: bool = False,
    surface_multiple_order: Optional[int] = None,
    complex_frequency_damping: Union[None, float, Literal["auto"]] = "auto",
    progress_callback: Optional[Callable[[float, int, int, str], None]] = None,
) -> KennettEarthResponse:
    """Build and cache all modal and motion responses in p-omega space."""
    frequencies = np.asarray(frequencies, dtype=float).reshape(-1)
    if frequencies.size < 2 or np.any(np.diff(frequencies) <= 0):
        raise ValueError("frequencies must be an increasing array")
    if max_frequency is None:
        max_frequency = frequencies[-1]
    if surface_multiple_order is not None:
        if not isinstance(surface_multiple_order, (int, np.integer)):
            raise TypeError("surface_multiple_order must be an integer or None")
        if surface_multiple_order < 0:
            raise ValueError("surface_multiple_order must be nonnegative")
    fft_duration = 1.0 / (frequencies[1] - frequencies[0])
    if complex_frequency_damping == "auto":
        # Suppress a periodic wraparound by 80 dB at the FFT boundary.
        damping = -np.log(1.0e-4) / fft_duration
    elif complex_frequency_damping is None:
        damping = 0.0
    else:
        damping = float(complex_frequency_damping)
        if not np.isfinite(damping) or damping < 0:
            raise ValueError("complex_frequency_damping must be nonnegative")
    active = np.flatnonzero((frequencies > 0.0) & (frequencies <= max_frequency))
    model = precompute_kennett_model(depth, vp, vs, rho, slowness)
    responses = {
        name: np.zeros((frequencies.size, model.slowness.size), dtype=complex)
        for name in COMPONENTS
    }
    total = active.size * model.slowness.size
    completed = 0

    def report(fraction, message):
        if progress_callback is not None:
            progress_callback(float(fraction), completed, total, message)

    report(0.0, "Precomputing Kennett model")
    for number, jf in enumerate(active, start=1):
        omega = 2.0 * np.pi * frequencies[jf] - 1j * damping
        values = _batched_component_responses(
            model, omega, internal_multiples, surface_multiples,
            surface_multiple_order,
        )
        if not primaries:
            baseline = _batched_component_responses(model, omega, False, False, None)
            values = {name: values[name] - baseline[name] for name in COMPONENTS}
        for name in COMPONENTS:
            responses[name][jf] = values[name]
        completed += model.slowness.size
        report(
            0.95 * completed / max(total, 1),
            f"Frequency {number} of {active.size}: {frequencies[jf]:.1f} Hz",
        )

    return KennettEarthResponse(
        model=model,
        frequencies=frequencies,
        active_indices=active,
        responses=responses,
        internal_multiples=internal_multiples,
        surface_multiples=surface_multiples,
        surface_multiple_order=surface_multiple_order,
        primaries=primaries,
        complex_frequency_damping=damping,
    )


def plane_wave_receiver_response(
    p: float,
    omega: float,
    depth,
    vp,
    vs,
    rho,
    *,
    internal_multiples: bool = True,
    surface_multiples: bool = False,
    primaries: bool = True,
    receiver_component: str = "modal_pp",
) -> complex:
    """Selected receiver response to a unit downgoing P plane wave.

    ``modal_pp`` and ``modal_ps`` return upgoing P and SV modal amplitudes.
    The PP and PS velocity outputs project each arriving mode separately;
    ``vertical_velocity`` and ``horizontal_velocity`` return their sums.

    ``primaries=False`` subtracts the response computed without internal or
    surface feedback, leaving the selected multiple component.
    """
    depth, vp, vs, rho = _validate_model(depth, vp, vs, rho)
    if receiver_component not in COMPONENTS:
        raise ValueError(
            f"receiver_component must be one of {COMPONENTS}"
        )

    def response(with_internal: bool, with_surface: bool) -> complex:
        earth = stack_scattering(p, omega, depth, vp, vs, rho, with_internal)[:2, :2]
        source = np.array([1.0, 0.0], dtype=complex)
        if with_surface:
            rfs = free_surface_reflection(p, vp[0], vs[0], rho[0])
            up = np.linalg.solve(np.eye(2) - earth @ rfs, earth @ source)
            down_at_receiver = rfs @ up
        else:
            up = earth @ source
            down_at_receiver = np.zeros(2, dtype=complex)

        if receiver_component == "modal_pp":
            return complex(up[0])
        if receiver_component == "modal_ps":
            return complex(up[1])

        down_modes, up_modes = modal_matrices(p, vp[0], vs[0], rho[0])
        p_displacement = (
            up_modes[:2, 0] * up[0]
            + down_modes[:2, 0] * down_at_receiver[0]
        )
        sv_displacement = (
            up_modes[:2, 1] * up[1]
            + down_modes[:2, 1] * down_at_receiver[1]
        )
        p_velocity = -1j * omega * p_displacement
        sv_velocity = -1j * omega * sv_displacement
        components = {
            "vertical_pp": p_velocity[1],
            "horizontal_pp": p_velocity[0],
            "vertical_ps": sv_velocity[1],
            "horizontal_ps": sv_velocity[0],
            "vertical_velocity": p_velocity[1] + sv_velocity[1],
            "horizontal_velocity": p_velocity[0] + sv_velocity[0],
        }
        return complex(components[receiver_component])

    selected = response(internal_multiples, surface_multiples)
    if not primaries:
        selected -= response(False, False)
    return selected


def plane_wave_pp_response(
    p: float,
    omega: float,
    depth,
    vp,
    vs,
    rho,
    *,
    internal_multiples: bool = True,
    surface_multiples: bool = False,
    primaries: bool = True,
) -> complex:
    """Backward-compatible alias for the modal P-P receiver response."""
    return plane_wave_receiver_response(
        p, omega, depth, vp, vs, rho,
        internal_multiples=internal_multiples,
        surface_multiples=surface_multiples,
        primaries=primaries,
        receiver_component="modal_pp",
    )


def _deposit_linear(trace: Array, sample_positions: Array, amplitudes: Array) -> None:
    """Deposit amplitudes at fractional samples using linear interpolation."""
    i0 = np.floor(sample_positions).astype(int)
    frac = sample_positions - i0
    valid0 = (i0 >= 0) & (i0 < trace.size)
    np.add.at(trace, i0[valid0], amplitudes[valid0] * (1.0 - frac[valid0]))
    i1 = i0 + 1
    valid1 = (i1 >= 0) & (i1 < trace.size)
    np.add.at(trace, i1[valid1], amplitudes[valid1] * frac[valid1])


def _cosine_edge_taper(size: int, fraction: float) -> Array:
    """Symmetric raised-cosine taper over a fraction of each array edge."""
    if not 0.0 <= fraction <= 0.5:
        raise ValueError("slowness_taper_fraction must lie in [0, 0.5]")
    taper = np.ones(size, dtype=float)
    if fraction == 0.0:
        return taper
    edge = min(int(np.ceil(fraction * size)), size // 2)
    if edge:
        ramp = 0.5 * (1.0 - np.cos(np.pi * np.arange(edge) / edge))
        taper[:edge] = ramp
        taper[-edge:] = ramp[::-1]
    return taper


def _upper_frequency_taper(frequencies, active_indices, fraction: float) -> Array:
    """Cosine taper over the upper fraction of the active frequency band."""
    if not 0.0 <= fraction < 1.0:
        raise ValueError("frequency_taper_fraction must lie in [0, 1)")
    taper = np.zeros_like(frequencies, dtype=float)
    taper[active_indices] = 1.0
    if fraction == 0.0 or active_indices.size == 0:
        return taper
    fmax = frequencies[active_indices[-1]]
    start = fmax * (1.0 - fraction)
    upper = active_indices[frequencies[active_indices] > start]
    if upper.size:
        phase = (frequencies[upper] - start) / max(fmax - start, np.finfo(float).eps)
        taper[upper] = 0.5 * (1.0 + np.cos(np.pi * phase))
    return taper


def synthesize_kennett_cmp(
    earth_response: KennettEarthResponse,
    offsets,
    wavelet,
    *,
    dt: float,
    record_length: float,
    receiver_component: str = "modal_pp",
    nmo_corrected: bool = False,
    nmo_velocity=None,
    stretch_mute: Optional[float] = None,
    normalize: Literal[None, "max"] = None,
    slowness_taper_fraction: float = 0.10,
    frequency_taper_fraction: float = 0.10,
    progress_callback: Optional[Callable[[float, int, int, str], None]] = None,
) -> KennettResult:
    """Synthesize offsets, wavelet and NMO from a cached earth response.

    Reusing ``earth_response`` avoids repeating the expensive Kennett recursion
    when changing offsets, wavelets, receiver components or NMO parameters.
    """
    offsets = np.asarray(offsets, dtype=float).reshape(-1)
    wavelet = np.asarray(wavelet, dtype=float).reshape(-1)
    if offsets.size == 0 or np.any(~np.isfinite(offsets)):
        raise ValueError("offsets must be a nonempty finite array")
    if wavelet.size == 0 or np.any(~np.isfinite(wavelet)):
        raise ValueError("wavelet must be a nonempty finite array")
    if receiver_component not in COMPONENTS + ("all",):
        raise ValueError("invalid receiver_component")
    if dt <= 0 or record_length <= 0:
        raise ValueError("dt and record_length must be positive")

    frequencies = earth_response.frequencies
    nfft = 2 * (frequencies.size - 1)
    expected = np.fft.rfftfreq(nfft, dt)
    if not np.allclose(frequencies, expected, rtol=1e-10, atol=1e-12):
        raise ValueError("dt is inconsistent with the cached frequency sampling")
    ns = int(np.floor(record_length / dt)) + 1
    if ns > nfft:
        raise ValueError("record_length exceeds the cached FFT duration")

    names = COMPONENTS if receiver_component == "all" else (receiver_component,)
    p = earth_response.model.slowness
    dp = np.diff(p)
    if not np.allclose(dp, dp[0]):
        raise ValueError("cached slowness samples must be uniformly spaced")
    weights = np.full(p.size, dp[0])
    weights[[0, -1]] *= 0.5
    weights *= _cosine_edge_taper(p.size, slowness_taper_fraction)
    frequency_taper = _upper_frequency_taper(
        frequencies, earth_response.active_indices, frequency_taper_fraction
    )

    frequency_gathers = {
        name: np.zeros((offsets.size, frequencies.size), dtype=complex)
        for name in names
    }
    for jf in earth_response.active_indices:
        omega = 2.0 * np.pi * frequencies[jf]
        kernel = np.exp(1j * omega * np.outer(offsets, p))
        weighted = np.stack([
            weights * earth_response.responses[name][jf] for name in names
        ])
        transformed = frequency_taper[jf] * (kernel @ weighted.T)
        for ic, name in enumerate(names):
            frequency_gathers[name][:, jf] = transformed[:, ic]

    if progress_callback is not None:
        total = earth_response.active_indices.size * p.size
        progress_callback(0.96, total, total, "Transforming frequency gathers to time")

    impulse_gathers = {
        name: np.fft.irfft(spectrum, n=nfft, axis=1)[:, :ns]
        for name, spectrum in frequency_gathers.items()
    }

    damping = earth_response.complex_frequency_damping
    if damping > 0.0:
        compensation = np.exp(damping * (np.arange(ns) * dt))
        impulse_gathers = {
            name: impulse * compensation[None, :]
            for name, impulse in impulse_gathers.items()
        }

    if progress_callback is not None:
        progress_callback(0.97, total, total, "Applying seismic wavelet")
    # Convolve the compensated causal impulse response linearly. This avoids
    # reintroducing circular wraparound during wavelet application.
    centre = wavelet.size // 2
    convolution_nfft = next_power_of_two(ns + wavelet.size - 1)
    wavelet_spectrum = np.fft.rfft(wavelet, n=convolution_nfft)
    uncorrected = {
        name: np.fft.irfft(
            np.fft.rfft(impulse, n=convolution_nfft, axis=1)
            * wavelet_spectrum[None, :],
            n=convolution_nfft,
            axis=1,
        )[:, centre:centre + ns]
        for name, impulse in impulse_gathers.items()
    }
    time = np.arange(ns) * dt

    used_nmo_velocity = None
    if nmo_corrected:
        if progress_callback is not None:
            progress_callback(0.98, total, total, "Applying NMO correction")
        if nmo_velocity is None:
            used_nmo_velocity = rms_velocity_from_layers(
                time, earth_response.model.depth, earth_response.model.vp
            )
        else:
            supplied = np.asarray(nmo_velocity, dtype=float)
            used_nmo_velocity = (
                np.full(ns, float(supplied))
                if supplied.ndim == 0 else supplied.reshape(-1)
            )
        gathers = {
            name: apply_nmo_correction(
                trace, time, offsets, used_nmo_velocity,
                stretch_mute=stretch_mute,
            )
            for name, trace in uncorrected.items()
        }
    else:
        gathers = {name: trace.copy() for name, trace in uncorrected.items()}

    if normalize == "max":
        if receiver_component == "all":
            # Preserve PP+PS summation identities and relative modal amplitudes.
            groups = (
                ("modal_pp", "modal_ps"),
                (
                    "vertical_pp", "horizontal_pp",
                    "vertical_ps", "horizontal_ps",
                    "vertical_velocity", "horizontal_velocity",
                ),
            )
        else:
            groups = (names,)
        for group in groups:
            scale = max(np.max(np.abs(gathers[name])) for name in group)
            if scale > 0:
                for name in group:
                    gathers[name] /= scale
                    uncorrected[name] /= scale
                    impulse_gathers[name] /= scale
    elif normalize is not None:
        raise ValueError("normalize must be None or 'max'")

    if progress_callback is not None:
        progress_callback(1.0, total, total, "Modelling complete")

    if receiver_component == "all":
        selected_gather = gathers
        selected_impulse = impulse_gathers
        selected_frequency = frequency_gathers
        selected_uncorrected = uncorrected
    else:
        selected_gather = gathers[receiver_component]
        selected_impulse = impulse_gathers[receiver_component]
        selected_frequency = frequency_gathers[receiver_component]
        selected_uncorrected = uncorrected[receiver_component]

    return KennettResult(
        time=time,
        offsets=offsets,
        gather=selected_gather,
        impulse_gather=selected_impulse,
        frequency_gather=selected_frequency,
        frequencies=frequencies,
        slowness=p,
        uncorrected_gather=selected_uncorrected,
        nmo_velocity=used_nmo_velocity,
        receiver_component=receiver_component,
        gathers=gathers,
        earth_response=earth_response,
    )


def model_kennett_cmp(
    depth,
    vp,
    vs,
    rho,
    offsets,
    wavelet,
    *,
    dt: float = 0.002,
    record_length: float = 3.0,
    p_max: Optional[float] = None,
    n_slowness: int = 121,
    nfft: Optional[int] = None,
    max_frequency: Optional[float] = None,
    primaries: bool = True,
    internal_multiples: bool = True,
    surface_multiples: bool = False,
    surface_multiple_order: Optional[int] = None,
    complex_frequency_damping: Union[None, float, Literal["auto"]] = "auto",
    receiver_component: str = "modal_pp",
    nmo_corrected: bool = False,
    nmo_velocity=None,
    stretch_mute: Optional[float] = None,
    normalize: Literal[None, "max"] = None,
    slowness_taper_fraction: float = 0.10,
    frequency_taper_fraction: float = 0.10,
    progress_callback: Optional[
        Callable[[float, int, int, str], None]
    ] = None,
    cached_earth_response: Optional[KennettEarthResponse] = None,
) -> KennettResult:
    """Synthesize a 2-D line-source P-P CMP gather.

    Parameters
    ----------
    depth, vp, vs, rho : array_like
        Tops and elastic properties for blocked layers. The final layer is a
        half-space.
    offsets : array_like
        Signed or unsigned source-receiver offsets.
    wavelet : array_like
        Centred, normally odd-length, zero-phase wavelet sampled at ``dt``.
    p_max : float, optional
        Maximum absolute horizontal slowness [s/distance]. The default is
        ``0.98 / max(vp)`` and therefore retains propagating P waves. Increase
        it cautiously to include wider-angle/evanescent content.
    n_slowness : int
        Number of symmetric p samples; an odd value is recommended.
    max_frequency : float, optional
        Frequencies above this value are skipped. Choose this just beyond the
        useful wavelet bandwidth for speed.
    primaries, internal_multiples, surface_multiples : bool
        Select response components. With ``primaries=False``, the practical
        no-feedback primary approximation is subtracted.
    surface_multiple_order : int or None
        Maximum number of free-surface reverberations. Zero retains no
        surface-related reverberation, one adds the first order, and so on.
        None uses the infinite matrix feedback series.
    complex_frequency_damping : {"auto", float, None}
        Exponential damping constant in 1/s used to suppress periodic FFT
        wraparound. ``"auto"`` gives 80 dB attenuation at the FFT boundary;
        a nonnegative float specifies the constant and None disables it.
    receiver_component : str
        Output representation: ``"modal_pp"`` for the arriving upgoing P
        mode; ``"modal_ps"`` for the arriving SV mode; ``"vertical_pp"``,
        ``"horizontal_pp"``, ``"vertical_ps"``, and ``"horizontal_ps"``
        for mode-separated particle velocities; and ``"vertical_velocity"``
        or ``"horizontal_velocity"`` for their PP+PS sums. Use ``"all"`` to
        synthesize all eight outputs from one earth recursion.
    nmo_corrected : bool
        If True, return a conventionally NMO-corrected ``gather``.
    nmo_velocity : float or array_like, optional
        Scalar or time-sampled NMO velocity. If omitted, a vertical RMS
        velocity function is calculated from the layered P velocities.
    stretch_mute : float, optional
        Fractional NMO-stretch mute, such as 0.30. Used only when
        ``nmo_corrected=True``.
    progress_callback : callable, optional
        Called as ``callback(fraction, completed, total, message)``. The
        fraction increases monotonically from 0.0 to 1.0. ``completed`` and
        ``total`` count the completed and required (frequency, slowness)
        response calculations. The callback may print progress, update a GUI,
        or emit a signal; modelling is silent when it is None.
    cached_earth_response : KennettEarthResponse, optional
        Reuse a response returned by an earlier result's ``earth_response``.
        This skips the expensive layered-earth recursion.
    slowness_taper_fraction : float
        Raised-cosine taper fraction applied at each p-aperture edge.
    frequency_taper_fraction : float
        Raised-cosine taper fraction applied below ``max_frequency``.

    Notes
    -----
    The offset transform is a direct quadrature of exp(i omega p x), appropriate
    to a 2-D line source. Absolute amplitudes depend on p sampling and source
    normalization; use the result primarily for relative amplitude studies.
    """
    depth, vp, vs, rho = _validate_model(depth, vp, vs, rho)
    offsets = np.asarray(offsets, dtype=float).reshape(-1)
    wavelet = np.asarray(wavelet, dtype=float).reshape(-1)
    if offsets.size == 0 or np.any(~np.isfinite(offsets)):
        raise ValueError("offsets must be a nonempty finite array")
    if wavelet.size == 0 or np.any(~np.isfinite(wavelet)):
        raise ValueError("wavelet must be a nonempty finite array")
    if dt <= 0 or record_length <= 0:
        raise ValueError("dt and record_length must be positive")
    if n_slowness < 3:
        raise ValueError("n_slowness must be at least 3")

    ns = int(np.floor(record_length / dt)) + 1
    if nfft is None:
        nfft = next_power_of_two(max(2 * ns, ns + wavelet.size))
    if nfft < ns:
        raise ValueError("nfft must be at least the output sample count")
    freqs = np.fft.rfftfreq(nfft, dt)
    if max_frequency is None:
        max_frequency = freqs[-1]
    if not (0 < max_frequency <= freqs[-1]):
        raise ValueError("max_frequency must lie in (0, Nyquist]")

    if p_max is None:
        p_max = 0.98 / np.max(vp)
    if p_max <= 0:
        raise ValueError("p_max must be positive")
    p = np.linspace(-p_max, p_max, n_slowness)
    if cached_earth_response is None:
        earth = build_kennett_earth_response(
            depth, vp, vs, rho, freqs, p,
            max_frequency=max_frequency,
            primaries=primaries,
            internal_multiples=internal_multiples,
            surface_multiples=surface_multiples,
            surface_multiple_order=surface_multiple_order,
            complex_frequency_damping=complex_frequency_damping,
            progress_callback=progress_callback,
        )
    else:
        earth = cached_earth_response
        cached_model = earth.model
        for supplied, cached, name in (
            (depth, cached_model.depth, "depth"),
            (vp, cached_model.vp, "vp"),
            (vs, cached_model.vs, "vs"),
            (rho, cached_model.rho, "rho"),
        ):
            if supplied.shape != cached.shape or not np.allclose(supplied, cached):
                raise ValueError(f"cached earth response has incompatible {name}")
        if not np.allclose(earth.frequencies, freqs):
            raise ValueError("cached earth response has incompatible FFT sampling")
        if not np.allclose(earth.model.slowness, p):
            raise ValueError("cached earth response has incompatible slowness sampling")
        expected_active = np.flatnonzero(
            (freqs > 0.0) & (freqs <= max_frequency)
        )
        if not np.array_equal(earth.active_indices, expected_active):
            raise ValueError("cached earth response has incompatible frequency band")
        if complex_frequency_damping == "auto":
            expected_damping = -np.log(1.0e-4) * (freqs[1] - freqs[0])
        elif complex_frequency_damping is None:
            expected_damping = 0.0
        else:
            expected_damping = float(complex_frequency_damping)
        if not np.isclose(
            earth.complex_frequency_damping, expected_damping,
            rtol=1e-12, atol=1e-14,
        ):
            raise ValueError("cached earth response uses different frequency damping")
        if not (
            earth.internal_multiples == internal_multiples
            and earth.surface_multiples == surface_multiples
            and earth.surface_multiple_order == surface_multiple_order
            and earth.primaries == primaries
        ):
            raise ValueError("cached earth response uses different multiple switches")
        if progress_callback is not None:
            total = earth.active_indices.size * p.size
            progress_callback(0.95, total, total, "Using cached Kennett earth response")

    return synthesize_kennett_cmp(
        earth,
        offsets,
        wavelet,
        dt=dt,
        record_length=record_length,
        receiver_component=receiver_component,
        nmo_corrected=nmo_corrected,
        nmo_velocity=nmo_velocity,
        stretch_mute=stretch_mute,
        normalize=normalize,
        slowness_taper_fraction=slowness_taper_fraction,
        frequency_taper_fraction=frequency_taper_fraction,
        progress_callback=progress_callback,
    )


def ricker_wavelet(frequency: float, dt: float, duration: float = 0.256) -> tuple[Array, Array]:
    """Return a centred zero-phase Ricker wavelet."""
    if frequency <= 0 or dt <= 0 or duration <= 0:
        raise ValueError("frequency, dt and duration must be positive")
    half = int(np.ceil(duration / (2.0 * dt)))
    t = np.arange(-half, half + 1) * dt
    a = np.pi * frequency * t
    return t, (1.0 - 2.0 * a**2) * np.exp(-a**2)


def calc_min_n_slowness(vp_max, freq_max, offset_max):
    """Estimate the minimum odd p-slowness sample count.

    Uses the phase-sampling criterion ``delta_p <= 1/(2*fmax*xmax)`` with
    the default aperture ``p_max = 0.98/vp_max``.
    """
    vp_max = float(vp_max)
    freq_max = float(freq_max)
    offset_max = float(offset_max)
    if not np.isfinite(vp_max) or vp_max <= 0.0:
        raise ValueError("vp_max must be finite and positive")
    if not np.isfinite(freq_max) or freq_max <= 0.0:
        raise ValueError("freq_max must be finite and positive")
    if not np.isfinite(offset_max) or offset_max < 0.0:
        raise ValueError("offset_max must be finite and nonnegative")
    n_slowness_min = int(np.ceil(1.0 + 4.0 * 0.98 / vp_max * freq_max * offset_max))
    if n_slowness_min % 2 == 0:
        n_slowness_min += 1
    return max(n_slowness_min, 3)
