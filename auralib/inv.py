"""One-dimensional sparse-spike inversion for seismic arrays.

The public :func:`sparse_spike_inversion` function accepts 1-D traces, 2-D
sections, 3-D volumes, and higher-dimensional arrays.  Inversion is performed
independently along ``sample_axis``; there is no lateral regularization.

Two solvers are provided:

``method="fista"``
    L1-regularized inversion using FISTA, with optional least-squares
    amplitude debiasing.

``method="ml"``
    A transparent greedy maximum-likelihood/MAP-style sparse-spike solver
    using Gaussian noise and a Bernoulli spike-occurrence prior.

Requires NumPy and SciPy.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.signal import convolve, correlate

__all__ = ["sparse_spike_inversion"]


def sparse_spike_inversion(
    seismic: Any,
    wavelet: Any,
    method: str = "fista",
    sample_axis: int = -1,
    lam: float = 0.01,
    niter: int = 1000,
    tol: float = 1e-7,
    noise_std: Any = None,
    spike_probability: float = 0.02,
    max_spikes: int | None = None,
    debias: bool = True,
    debias_threshold: float = 1e-3,
    debias_method: str = "cgls",
    debias_niter: int = 50,
    debias_tol: float = 1e-5,
    debias_damping: float = 0.0,
    fista_backend: str = "fft",
    chunk_size: int = 256,
    fft_workers: int | None = None,
    ml_candidates: int = 20,
    return_info: bool = True,
    progress: bool = False,
):
    """Estimate sparse reflectivity along one axis of a seismic array.

    Parameters
    ----------
    seismic : array_like
        Seismic data. Common shapes are ``(samples,)``, ``(traces, samples)``,
        and ``(inlines, crosslines, samples)``.
    wavelet : array_like, shape (n_wavelet,)
        One-dimensional seismic wavelet. Its zero-time sample is assumed to
        be at index ``(n_wavelet - 1) // 2``.
    method : {'fista', 'ml'}, default='fista'
        Inversion algorithm. ``'maximum_likelihood'`` is accepted as an alias
        for ``'ml'``.
    sample_axis : int, default=-1
        Axis containing time or depth samples.
    lam : float, default=0.01
        FISTA L1 penalty. Larger values generally produce fewer coefficients.
    niter : int, default=1000
        Maximum FISTA iterations per trace.
    tol : float, default=1e-7
        Relative FISTA convergence tolerance.
    noise_std : float, array_like, or None
        ML noise standard deviation. A scalar applies to all traces. An array
        must broadcast to the seismic shape with ``sample_axis`` removed. If
        None, it is estimated independently on each trace.
    spike_probability : float, default=0.02
        ML prior probability that a sample contains a reflector.
    max_spikes : int or None
        Maximum ML spikes per trace. None permits up to ``n_samples``.
    debias : bool, default=True
        Re-estimate FISTA amplitudes on the selected support by least squares.
    debias_threshold : float, default=1e-3
        Relative support threshold for FISTA debiasing. Coefficients smaller
        than this fraction of the maximum magnitude are excluded.
    debias_method : {'cgls', 'lstsq'}, default='cgls'
        Debiasing solver for the FFT backend. ``'cgls'`` solves all traces in
        a chunk with batched FFT operators and is normally much faster.
        ``'lstsq'`` retains the dense, trace-by-trace reference calculation.
    debias_niter : int, default=50
        Maximum CGLS iterations used for debiasing.
    debias_tol : float, default=1e-5
        Relative gradient tolerance for CGLS debiasing.
    debias_damping : float, default=0.0
        Optional non-negative L2 damping added to the restricted least-squares
        problem. Small values can stabilize strongly overlapping reflectors.
    fista_backend : {'fft', 'trace'}, default='fft'
        ``'fft'`` uses batched FFT convolution and is normally much faster for
        sections and volumes. ``'trace'`` retains the original reference
        implementation and is useful for comparison.
    chunk_size : int, default=256
        Number of traces processed simultaneously by batched FFT FISTA.
        Increase this when memory permits; decrease it for long traces or
        memory-constrained machines.
    fft_workers : int or None, default=None
        Number of workers supplied to SciPy's FFT routines. None uses SciPy's
        default. ``-1`` requests all available CPU workers.
    ml_candidates : int, default=20
        Number of strongest unused ML candidate locations tested per step.
    return_info : bool, default=True
        Return ``(reflectivity, info)`` when True, otherwise reflectivity only.
    progress : bool, default=False
        Print progress while processing multiple traces.

    Returns
    -------
    reflectivity : ndarray
        Reflectivity with the same shape and axis order as ``seismic``.
    info : dict, optional
        Diagnostics. ``synthetic`` and ``residual`` have the seismic shape.
        Per-trace diagnostic arrays have the spatial shape (all axes except
        ``sample_axis``). For a 1-D input, they are zero-dimensional arrays.

    Notes
    -----
    FISTA approximately minimizes

    ``0.5 * ||W r - s||_2**2 + lam * ||r||_1``.

    The ML option greedily minimizes

    ``||s - W r||_2**2 / (2 * noise_std**2) + K*log((1-p)/p)``.

    The ML implementation is educational and transparent; it is not a clone
    of a proprietary commercial sparse-spike inversion algorithm.
    """
    data = np.asarray(seismic, dtype=float)
    wav = np.asarray(wavelet, dtype=float).ravel()
    _validate_inputs(data, wav)

    if not isinstance(sample_axis, (int, np.integer)):
        raise TypeError("sample_axis must be an integer.")
    axis = np.core.numeric.normalize_axis_index(sample_axis, data.ndim)
    method_key = _normalize_method(method)

    if lam < 0 or not np.isfinite(lam):
        raise ValueError("lam must be finite and non-negative.")
    if niter < 1:
        raise ValueError("niter must be at least 1.")
    if tol < 0 or not np.isfinite(tol):
        raise ValueError("tol must be finite and non-negative.")
    if not 0 <= debias_threshold < 1:
        raise ValueError("debias_threshold must be in [0, 1).")
    if not isinstance(debias_method, str):
        raise TypeError("debias_method must be a string.")
    debias_method = debias_method.lower().strip()
    if debias_method not in {"cgls", "lstsq"}:
        raise ValueError("debias_method must be 'cgls' or 'lstsq'.")
    if debias_niter < 1:
        raise ValueError("debias_niter must be at least 1.")
    if debias_tol < 0 or not np.isfinite(debias_tol):
        raise ValueError("debias_tol must be finite and non-negative.")
    if debias_damping < 0 or not np.isfinite(debias_damping):
        raise ValueError("debias_damping must be finite and non-negative.")
    if ml_candidates < 1:
        raise ValueError("ml_candidates must be at least 1.")
    if not isinstance(fista_backend, str):
        raise TypeError("fista_backend must be a string.")
    fista_backend = fista_backend.lower().strip()
    if fista_backend not in {"fft", "trace"}:
        raise ValueError("fista_backend must be 'fft' or 'trace'.")
    if chunk_size < 1:
        raise ValueError("chunk_size must be at least 1.")

    work = np.moveaxis(data, axis, -1)
    spatial_shape = work.shape[:-1]
    n_samples = work.shape[-1]
    n_traces = int(np.prod(spatial_shape, dtype=int)) if spatial_shape else 1
    traces = work.reshape(n_traces, n_samples)
    noise_values = _prepare_noise_std(noise_std, spatial_shape, n_traces)

    reflectivity_flat = np.empty_like(traces)
    synthetic_flat = np.empty_like(traces)
    residual_flat = np.empty_like(traces)
    n_spikes = np.empty(n_traces, dtype=int)
    trace_info = np.empty(n_traces, dtype=object)

    if method_key == "fista":
        iterations = np.empty(n_traces, dtype=int)
        converged = np.empty(n_traces, dtype=bool)
        debias_iterations_all = np.zeros(n_traces, dtype=int)
        debias_converged_all = np.zeros(n_traces, dtype=bool)
    else:
        used_noise = np.empty(n_traces, dtype=float)
        supports = np.empty(n_traces, dtype=object)

    # Fast path: each chunk is inverted as one 2-D FFT batch. The ML solver
    # remains trace-wise because its greedy support search is sequential.
    if method_key == "fista" and fista_backend == "fft":
        fft_operator = _FFTConvolution(n_samples, wav, workers=fft_workers)
        for start in range(0, n_traces, chunk_size):
            stop = min(start + chunk_size, n_traces)
            result = _ssi_fista_fft_batch(
                traces[start:stop],
                fft_operator,
                lam,
                niter,
                tol,
                debias,
                debias_threshold,
                debias_method,
                debias_niter,
                debias_tol,
                debias_damping,
            )
            reflectivity_flat[start:stop] = result["reflectivity"]
            synthetic_flat[start:stop] = result["synthetic"]
            residual_flat[start:stop] = result["residual"]
            n_spikes[start:stop] = result["n_spikes"]
            iterations[start:stop] = result["iterations"]
            converged[start:stop] = result["converged"]
            debias_iterations_all[start:stop] = result["debias_iterations"]
            debias_converged_all[start:stop] = result["debias_converged"]
            for local, i in enumerate(range(start, stop)):
                trace_info[i] = {
                    "method": "fista",
                    "backend": "fft",
                    "iterations": int(result["iterations"][local]),
                    "converged": bool(result["converged"][local]),
                    "lambda": lam,
                    "step_size": result["step_size"],
                    "n_spikes": int(result["n_spikes"][local]),
                    "support": result["support"][local],
                    "synthetic": result["synthetic"][local],
                    "residual": result["residual"][local],
                    "objective": float(result["objective"][local]),
                    "debias_method": result["debias_method"],
                    "debias_iterations": int(result["debias_iterations"][local]),
                    "debias_converged": bool(result["debias_converged"][local]),
                }
            if progress:
                print(f"Inverted {stop:,} of {n_traces:,} traces")
    else:
        _invert_traces_reference(
            traces=traces,
            wavelet=wav,
            method=method_key,
            lam=lam,
            niter=niter,
            tol=tol,
            noise_values=noise_values,
            spike_probability=spike_probability,
            max_spikes=max_spikes,
            debias=debias,
            debias_threshold=debias_threshold,
            ml_candidates=ml_candidates,
            reflectivity=reflectivity_flat,
            synthetic=synthetic_flat,
            residual=residual_flat,
            n_spikes=n_spikes,
            trace_info=trace_info,
            iterations=iterations if method_key == "fista" else None,
            converged=converged if method_key == "fista" else None,
            used_noise=used_noise if method_key == "ml" else None,
            supports=supports if method_key == "ml" else None,
            progress=progress,
        )

    def restore_samples(array: np.ndarray) -> np.ndarray:
        shaped = array.reshape(work.shape)
        return np.moveaxis(shaped, -1, axis)

    reflectivity = restore_samples(reflectivity_flat)
    if not return_info:
        return reflectivity

    info = {
        "method": method_key,
        "input_shape": data.shape,
        "sample_axis": axis,
        "n_samples": n_samples,
        "n_traces": n_traces,
        "synthetic": restore_samples(synthetic_flat),
        "residual": restore_samples(residual_flat),
        "n_spikes": n_spikes.reshape(spatial_shape),
        "trace_info": trace_info.reshape(spatial_shape),
    }
    if method_key == "fista":
        info.update(
            backend=fista_backend,
            chunk_size=chunk_size if fista_backend == "fft" else None,
            debias_method=(debias_method if debias and fista_backend == "fft" else None),
            lambda_=lam,
            iterations=iterations.reshape(spatial_shape),
            converged=converged.reshape(spatial_shape),
            debias_iterations=debias_iterations_all.reshape(spatial_shape),
            debias_converged=debias_converged_all.reshape(spatial_shape),
        )
        info["lambda"] = info.pop("lambda_")
    else:
        info.update(
            spike_probability=spike_probability,
            noise_std=used_noise.reshape(spatial_shape),
            support=supports.reshape(spatial_shape),
        )
    return reflectivity, info


def _invert_traces_reference(
    *, traces, wavelet, method, lam, niter, tol, noise_values,
    spike_probability, max_spikes, debias, debias_threshold, ml_candidates,
    reflectivity, synthetic, residual, n_spikes, trace_info,
    iterations, converged, used_noise, supports, progress,
):
    """Original trace-at-a-time path, retained for ML and comparison."""
    n_traces = traces.shape[0]
    report_every = max(1, n_traces // 20)
    for i, trace in enumerate(traces):
        trace_noise = None if noise_values is None else float(noise_values[i])
        r, diag = _sparse_spike_inversion_1d(
            trace,
            wavelet,
            method=method,
            lam=lam,
            niter=niter,
            tol=tol,
            noise_std=trace_noise,
            spike_probability=spike_probability,
            max_spikes=max_spikes,
            debias=debias,
            debias_threshold=debias_threshold,
            ml_candidates=ml_candidates,
        )
        reflectivity[i] = r
        synthetic[i] = diag["synthetic"]
        residual[i] = diag["residual"]
        n_spikes[i] = diag["n_spikes"]
        trace_info[i] = diag

        if method == "fista":
            iterations[i] = diag["iterations"]
            converged[i] = diag["converged"]
        else:
            used_noise[i] = diag["noise_std"]
            supports[i] = diag["support"]

        if progress and ((i + 1) % report_every == 0 or i + 1 == n_traces):
            print(f"Inverted {i + 1:,} of {n_traces:,} traces")


def _sparse_spike_inversion_1d(
    seismic: np.ndarray,
    wavelet: np.ndarray,
    *,
    method: str,
    lam: float,
    niter: int,
    tol: float,
    noise_std: float | None,
    spike_probability: float,
    max_spikes: int | None,
    debias: bool,
    debias_threshold: float,
    ml_candidates: int,
):
    if method == "fista":
        return _ssi_fista(
            seismic, wavelet, lam, niter, tol, debias, debias_threshold
        )
    return _ssi_ml(
        seismic,
        wavelet,
        noise_std,
        spike_probability,
        max_spikes,
        ml_candidates,
    )


class _FFTConvolution:
    """Reusable batched FFT implementation of W and its exact adjoint."""

    def __init__(self, n_samples, wavelet, workers=None):
        self.n = int(n_samples)
        self.m = int(wavelet.size)
        self.origin = (self.m - 1) // 2
        self.full_length = self.n + self.m - 1
        self.nfft = next_fast_len(self.full_length)
        self.workers = workers
        self.wavelet = np.asarray(wavelet, dtype=float).copy()
        self.wavelet_fft = rfft(wavelet, self.nfft, workers=workers)
        self.reverse_wavelet_fft = rfft(wavelet[::-1], self.nfft, workers=workers)
        # The full-convolution spectral norm bounds the cropped operator norm.
        self.lipschitz = float(np.max(np.abs(self.wavelet_fft) ** 2))

    def forward(self, reflectivity):
        spectrum = rfft(reflectivity, self.nfft, axis=-1, workers=self.workers)
        full = irfft(
            spectrum * self.wavelet_fft,
            self.nfft,
            axis=-1,
            workers=self.workers,
        )[..., : self.full_length]
        return full[..., self.origin : self.origin + self.n]

    def adjoint(self, values):
        # C^T first embeds the cropped data into the full-convolution output;
        # convolution with the reversed wavelet then implements H^T.
        embedded = np.zeros(values.shape[:-1] + (self.nfft,), dtype=values.dtype)
        embedded[..., self.origin : self.origin + self.n] = values
        spectrum = rfft(embedded, self.nfft, axis=-1, workers=self.workers)
        full = irfft(
            spectrum * self.reverse_wavelet_fft,
            self.nfft,
            axis=-1,
            workers=self.workers,
        )
        return full[..., self.m - 1 : self.m - 1 + self.n]


def _ssi_fista_fft_batch(
    seismic, operator, lam, niter, tol, debias, support_threshold,
    debias_method, debias_niter, debias_tol, debias_damping,
):
    """Run FISTA simultaneously on a two-dimensional batch of traces."""
    n_traces, n_samples = seismic.shape
    if operator.lipschitz <= 0 or not np.isfinite(operator.lipschitz):
        raise ValueError("Could not determine a stable FISTA step size.")
    step = 1.0 / operator.lipschitz
    r = np.zeros_like(seismic)
    y = np.zeros_like(seismic)
    momentum = np.ones(n_traces)
    active = np.ones(n_traces, dtype=bool)
    converged = np.zeros(n_traces, dtype=bool)
    iterations = np.full(n_traces, niter, dtype=int)

    for iteration in range(niter):
        indices = np.flatnonzero(active)
        if indices.size == 0:
            break
        y_active = y[indices]
        residual = seismic[indices] - operator.forward(y_active)
        r_new = _soft_threshold(
            y_active + step * operator.adjoint(residual), step * lam
        )
        old_r = r[indices]
        old_momentum = momentum[indices]
        new_momentum = 0.5 * (
            1.0 + np.sqrt(1.0 + 4.0 * old_momentum**2)
        )
        y_new = r_new + (
            ((old_momentum - 1.0) / new_momentum)[:, None]
            * (r_new - old_r)
        )
        change = np.linalg.norm(r_new - old_r, axis=1) / np.maximum(
            np.linalg.norm(old_r, axis=1), 1e-12
        )
        r[indices] = r_new
        y[indices] = y_new
        momentum[indices] = new_momentum

        if iteration >= 5:
            finished_local = change < tol
            finished = indices[finished_local]
            converged[finished] = True
            iterations[finished] = iteration + 1
            active[finished] = False

    maximum = np.max(np.abs(r), axis=1)
    supports = [
        np.flatnonzero(np.abs(r[i]) > support_threshold * maximum[i])
        if maximum[i] > 0
        else np.empty(0, dtype=int)
        for i in range(n_traces)
    ]

    debias_iterations = np.zeros(n_traces, dtype=int)
    debias_converged = np.zeros(n_traces, dtype=bool)
    if debias:
        if debias_method == "cgls":
            support_mask = np.zeros_like(r, dtype=bool)
            for i, support in enumerate(supports):
                support_mask[i, support] = True
            r, debias_iterations, debias_converged = _debias_cgls_batch(
                seismic, r, support_mask, operator,
                debias_niter, debias_tol, debias_damping,
            )
        else:
            for i, support in enumerate(supports):
                if support.size:
                    dictionary = _build_wavelet_dictionary(
                        n_samples, operator.wavelet, support
                    )
                    amplitudes, *_ = np.linalg.lstsq(
                        dictionary, seismic[i], rcond=None
                    )
                    r[i] = 0.0
                    r[i, support] = amplitudes
                    debias_iterations[i] = 1
                    debias_converged[i] = True

    synthetic = operator.forward(r)
    residual = seismic - synthetic
    objective = (
        0.5 * np.sum(residual**2, axis=1)
        + lam * np.sum(np.abs(r), axis=1)
    )
    return {
        "reflectivity": r,
        "synthetic": synthetic,
        "residual": residual,
        "n_spikes": np.count_nonzero(r, axis=1),
        "support": supports,
        "iterations": iterations,
        "converged": converged,
        "objective": objective,
        "step_size": step,
        "debias_method": debias_method if debias else None,
        "debias_iterations": debias_iterations,
        "debias_converged": debias_converged,
    }


def _debias_cgls_batch(
    seismic, initial, support_mask, operator, niter, tol, damping
):
    """Solve support-restricted least squares with batched FFT CGLS."""
    x = np.where(support_mask, initial, 0.0).copy()
    data_residual = seismic - operator.forward(x)
    gradient = operator.adjoint(data_residual) - damping * x
    gradient *= support_mask
    direction = gradient.copy()
    gamma = np.sum(gradient * gradient, axis=1)
    initial_scale = np.maximum(np.sqrt(gamma), np.finfo(float).eps)

    active = gamma > np.finfo(float).eps
    converged = ~active
    iterations = np.zeros(seismic.shape[0], dtype=int)

    for iteration in range(niter):
        indices = np.flatnonzero(active)
        if indices.size == 0:
            break

        p = direction[indices]
        projected = operator.forward(p)
        denominator = np.sum(projected * projected, axis=1)
        if damping:
            denominator += damping * np.sum(p * p, axis=1)

        valid = denominator > np.finfo(float).eps
        alpha = np.zeros(indices.size)
        alpha[valid] = gamma[indices[valid]] / denominator[valid]
        x[indices] += alpha[:, None] * p
        x[indices] *= support_mask[indices]
        data_residual[indices] -= alpha[:, None] * projected

        gradient_new = operator.adjoint(data_residual[indices]) - damping * x[indices]
        gradient_new *= support_mask[indices]
        gamma_new = np.sum(gradient_new * gradient_new, axis=1)
        relative_gradient = np.sqrt(gamma_new) / initial_scale[indices]
        finished_local = (~valid) | (relative_gradient <= tol)

        beta = np.zeros(indices.size)
        can_update = (gamma[indices] > np.finfo(float).eps) & ~finished_local
        beta[can_update] = gamma_new[can_update] / gamma[indices[can_update]]
        direction[indices] = gradient_new + beta[:, None] * p
        gamma[indices] = gamma_new
        iterations[indices] = iteration + 1

        finished = indices[finished_local]
        converged[finished] = True
        active[finished] = False

    return x, iterations, converged


def _ssi_fista(seismic, wavelet, lam, niter, tol, debias, support_threshold):
    n = seismic.size
    # The full-convolution spectral norm bounds the cropped operator norm and
    # matches the step-size rule used by the batched FFT implementation.
    nfft = next_fast_len(n + wavelet.size - 1)
    lipschitz = float(
        np.max(np.abs(rfft(wavelet, nfft)) ** 2)
    )
    if lipschitz <= 0 or not np.isfinite(lipschitz):
        raise ValueError("Could not determine a stable FISTA step size.")
    step = 1.0 / lipschitz

    r = np.zeros(n)
    y = r.copy()
    momentum = 1.0
    objective = []
    converged = False

    for iteration in range(niter):
        residual_y = seismic - _forward(y, wavelet)
        r_new = _soft_threshold(
            y + step * _adjoint(residual_y, wavelet), step * lam
        )
        momentum_new = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * momentum**2))
        y = r_new + ((momentum - 1.0) / momentum_new) * (r_new - r)

        residual_new = seismic - _forward(r_new, wavelet)
        objective.append(
            0.5 * float(residual_new @ residual_new)
            + lam * float(np.sum(np.abs(r_new)))
        )
        change = np.linalg.norm(r_new - r) / max(np.linalg.norm(r), 1e-12)
        r = r_new
        momentum = momentum_new
        if iteration >= 5 and change < tol:
            converged = True
            break

    if np.any(r):
        cutoff = support_threshold * np.max(np.abs(r))
        support = np.flatnonzero(np.abs(r) > cutoff)
    else:
        support = np.empty(0, dtype=int)

    if debias and support.size:
        dictionary = _build_wavelet_dictionary(n, wavelet, support)
        amplitudes, *_ = np.linalg.lstsq(dictionary, seismic, rcond=None)
        r = np.zeros(n)
        r[support] = amplitudes

    synthetic = _forward(r, wavelet)
    residual = seismic - synthetic
    return r, {
        "method": "fista",
        "iterations": iteration + 1,
        "converged": converged,
        "lambda": lam,
        "step_size": step,
        "n_spikes": int(np.count_nonzero(r)),
        "support": np.flatnonzero(r),
        "synthetic": synthetic,
        "residual": residual,
        "objective": np.asarray(objective),
    }


def _ssi_ml(
    seismic,
    wavelet,
    noise_std,
    spike_probability,
    max_spikes,
    n_candidates,
):
    n = seismic.size
    if not 0.0 < spike_probability < 0.5:
        raise ValueError("spike_probability must be between 0 and 0.5.")
    if max_spikes is not None and max_spikes < 0:
        raise ValueError("max_spikes must be non-negative or None.")

    if noise_std is None:
        # First differences suppress slowly varying signal and make this more
        # useful than applying MAD directly to the seismic amplitudes.
        differences = np.diff(seismic)
        if differences.size:
            median = np.median(differences)
            noise_std = 1.4826 * np.median(np.abs(differences - median)) / np.sqrt(2.0)
        else:
            noise_std = 0.0
        if not np.isfinite(noise_std) or noise_std <= 0:
            noise_std = 0.01 * np.std(seismic) + np.finfo(float).eps
    elif not np.isfinite(noise_std) or noise_std <= 0:
        raise ValueError("noise_std must be positive and finite.")

    noise_var = float(noise_std) ** 2
    spike_penalty = float(np.log((1.0 - spike_probability) / spike_probability))
    limit = n if max_spikes is None else min(int(max_spikes), n)
    support: list[int] = []
    residual = seismic.copy()
    current_cost = float(residual @ residual) / (2.0 * noise_var)
    history = []

    for iteration in range(limit):
        scores = np.abs(_adjoint(residual, wavelet))
        if support:
            scores[np.asarray(support)] = -np.inf
        available = n - len(support)
        candidate_count = min(n_candidates, available)
        if candidate_count == 0:
            break
        candidates = np.argpartition(scores, -candidate_count)[-candidate_count:]
        candidates = candidates[np.argsort(scores[candidates])[::-1]]

        best = None
        best_cost = current_cost
        for candidate in candidates:
            trial_support = np.sort(
                np.append(np.asarray(support, dtype=int), int(candidate))
            )
            dictionary = _build_wavelet_dictionary(n, wavelet, trial_support)
            amplitudes, *_ = np.linalg.lstsq(dictionary, seismic, rcond=None)
            synthetic_trial = dictionary @ amplitudes
            residual_trial = seismic - synthetic_trial
            cost = (
                float(residual_trial @ residual_trial) / (2.0 * noise_var)
                + trial_support.size * spike_penalty
            )
            if cost < best_cost:
                best_cost = cost
                best = (trial_support.tolist(), amplitudes, residual_trial)

        if best is None:
            break
        previous_cost = current_cost
        support, amplitudes, residual = best
        current_cost = best_cost
        history.append(
            {
                "iteration": iteration + 1,
                "n_spikes": len(support),
                "cost": current_cost,
                "cost_reduction": previous_cost - current_cost,
            }
        )

    r = np.zeros(n)
    if support:
        support_array = np.asarray(support, dtype=int)
        dictionary = _build_wavelet_dictionary(n, wavelet, support_array)
        amplitudes, *_ = np.linalg.lstsq(dictionary, seismic, rcond=None)
        r[support_array] = amplitudes
    else:
        support_array = np.empty(0, dtype=int)

    synthetic = _forward(r, wavelet)
    residual = seismic - synthetic
    return r, {
        "method": "ml",
        "noise_std": float(noise_std),
        "spike_probability": spike_probability,
        "spike_penalty": spike_penalty,
        "n_spikes": len(support),
        "support": support_array,
        "synthetic": synthetic,
        "residual": residual,
        "cost_history": history,
    }


def _forward(reflectivity: np.ndarray, wavelet: np.ndarray) -> np.ndarray:
    """Cropped linear convolution with a documented wavelet origin."""
    n = reflectivity.size
    origin = (wavelet.size - 1) // 2
    full = convolve(reflectivity, wavelet, mode="full")
    return full[origin : origin + n]


def _adjoint(values: np.ndarray, wavelet: np.ndarray) -> np.ndarray:
    """Exact adjoint of :func:`_forward`, including its boundary cropping."""
    n = values.size
    origin = (wavelet.size - 1) // 2
    start = wavelet.size - 1 - origin
    full = correlate(values, wavelet, mode="full")
    return full[start : start + n]


def _soft_threshold(values: np.ndarray, threshold: float) -> np.ndarray:
    return np.sign(values) * np.maximum(np.abs(values) - threshold, 0.0)


def _build_wavelet_dictionary(n_samples, wavelet, locations):
    locations = np.asarray(locations, dtype=int)
    dictionary = np.empty((n_samples, locations.size), dtype=float)
    for column, location in enumerate(locations):
        spike = np.zeros(n_samples)
        spike[location] = 1.0
        dictionary[:, column] = _forward(spike, wavelet)
    return dictionary


def _prepare_noise_std(noise_std, spatial_shape, n_traces):
    if noise_std is None:
        return None
    values = np.asarray(noise_std, dtype=float)
    if values.ndim == 0:
        if not np.isfinite(values) or values <= 0:
            raise ValueError("noise_std must be positive and finite.")
        return np.full(n_traces, float(values))
    try:
        values = np.broadcast_to(values, spatial_shape).reshape(-1)
    except ValueError as exc:
        raise ValueError(
            f"noise_std must be scalar or broadcastable to {spatial_shape}."
        ) from exc
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("All noise_std values must be positive and finite.")
    return values


def _normalize_method(method: str) -> str:
    if not isinstance(method, str):
        raise TypeError("method must be a string.")
    key = method.lower().strip().replace("-", "_").replace(" ", "_")
    aliases = {"fista": "fista", "ml": "ml", "maximum_likelihood": "ml"}
    try:
        return aliases[key]
    except KeyError as exc:
        raise ValueError("method must be 'fista' or 'ml'.") from exc


def _validate_inputs(seismic: np.ndarray, wavelet: np.ndarray) -> None:
    if seismic.ndim == 0 or seismic.size == 0:
        raise ValueError("seismic must be a non-empty array with at least one axis.")
    if wavelet.size == 0:
        raise ValueError("wavelet must contain at least one sample.")
    if not np.all(np.isfinite(seismic)):
        raise ValueError("seismic contains NaN or infinite values.")
    if not np.all(np.isfinite(wavelet)):
        raise ValueError("wavelet contains NaN or infinite values.")
    if np.allclose(wavelet, 0.0):
        raise ValueError("wavelet cannot be identically zero.")


def reflectivity_to_impedance(
    reflectivity,
    initial_impedance=1.0,
    sample_axis=-1,
):
    """
    Convert reflection coefficients to impedance.

    Parameters
    ----------
    reflectivity : array_like
        Reflection-coefficient array.

    initial_impedance : float or array_like, default=1.0
        Impedance at the first sample.

        Use 1.0 to obtain relative impedance. Alternatively, provide an
        absolute starting impedance for every trace.

    sample_axis : int, default=-1
        Sample axis.

    Returns
    -------
    impedance : ndarray
        Relative or absolute impedance.
    """
    import numpy as np

    r = np.asarray(reflectivity, dtype=float)
    r_work = np.moveaxis(r, sample_axis, -1)

    #if np.any(np.abs(r_work) >= 1.0):
    #    raise ValueError(
    #        "All reflection coefficients must satisfy abs(r) < 1."
    #    )

    spatial_shape = r_work.shape[:-1]

    z0 = np.asarray(initial_impedance, dtype=float)
    z0 = np.broadcast_to(z0, spatial_shape)

    impedance = np.empty_like(r_work)
    impedance[..., 0] = z0

    ratios = (
        (1.0 + r_work[..., :-1])
        / (1.0 - r_work[..., :-1])
    )

    impedance[..., 1:] = (
        z0[..., None]
        * np.cumprod(ratios, axis=-1)
    )

    return np.moveaxis(
        impedance,
        -1,
        sample_axis,
    )