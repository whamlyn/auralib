"""Evaluate tensile and shear failure around the borehole wall."""

import numpy as np

from .failure import mohr_coulomb_margin, tensile_margin
from .wellbore import borehole_wall_stress, wall_principal_stresses


def mud_pressure_sweep(stress, theta, poisson, mud_pressure, pore_pressure,
                       ucs, friction_angle, tensile_strength):
    """Evaluate all combinations of mud pressure and wall angle.

    Supply one total borehole stress tensor (3, 3), scalar pore pressure,
    Poisson's ratio and strengths, and nonempty 1-D arrays of theta (degrees)
    and mud pressure. Use consistent stress and pressure units.

    Return (tensile, shear) margins, each shaped (n_pressures, n_angles).
    A negative margin exceeds the corresponding criterion and zero is onset.
    """
    stress = np.asarray(stress, dtype=float)
    theta = np.asarray(theta, dtype=float)
    mud_pressure = np.asarray(mud_pressure, dtype=float)
    if stress.shape != (3, 3):
        raise ValueError("stress must have shape (3, 3)")
    if theta.ndim != 1 or mud_pressure.ndim != 1 or not (theta.size and mud_pressure.size):
        raise ValueError("theta and mud_pressure must be nonempty 1-D arrays")
    if any(np.ndim(value) for value in
           (poisson, pore_pressure, ucs, friction_angle, tensile_strength)):
        raise ValueError("poisson, pore_pressure and strength parameters must be scalars")

    wall = borehole_wall_stress(
        stress, theta, poisson, mud_pressure[:, None], pore_pressure
    )
    principal = wall_principal_stresses(wall)
    sigma1, sigma3 = principal[..., 0], principal[..., 2]
    return (tensile_margin(sigma3, tensile_strength),
            mohr_coulomb_margin(sigma1, sigma3, ucs, friction_angle))


def mud_pressure_window(stress, theta, poisson, mud_pressure, pore_pressure,
                        ucs, friction_angle, tensile_strength):
    """Find the modeled pressure window over the supplied wall angles.

    Supply the same inputs as mud_pressure_sweep, with at least two
    strictly increasing mud pressures. The first and last pressures
    bound the search. Wall angles are evaluated at the supplied samples.

    Return (lower, upper), or None if no feasible pressure is found.
    Each bound is a dict with pressure, criterion, and theta (degrees).
    A bound at a search endpoint with a positive margin has criterion
    and theta set to None. For tied minima, return the first criterion
    (tensile before shear) and wall angle in input order.
    """
    from scipy.optimize import brentq, minimize_scalar

    pressures = np.asarray(mud_pressure, dtype=float)
    theta = np.asarray(theta, dtype=float)
    if (pressures.ndim != 1 or pressures.size < 2
            or not np.all(np.isfinite(pressures))
            or not np.all(np.diff(pressures) > 0)):
        raise ValueError("mud_pressure must contain increasing finite pressures")
    if not all(np.all(np.isfinite(value)) for value in (stress, theta, pore_pressure)):
        raise ValueError("stress, theta, and pore_pressure must be finite")

    def evaluate(pressure):
        return np.stack(mud_pressure_sweep(
            stress, theta, poisson, pressure, pore_pressure,
            ucs, friction_angle, tensile_strength,
        ))

    def margin(pressure):
        return evaluate([pressure]).min()

    margins = evaluate(pressures).min(axis=(0, 2))
    left, right = pressures[[0, -1]]
    best = pressures[np.argmax(margins)]
    if margins.max() < 0:
        # The minimum margin is concave in pressure. A narrow feasible
        # interval can lie entirely between two failed pressure samples.
        peak = minimize_scalar(
            lambda pressure: -margin(pressure), bounds=(left, right),
            method="bounded",
        )
        if not peak.success:
            raise RuntimeError("could not locate the peak failure margin")
        best = peak.x
        if margin(best) < 0:
            return None

    bounds = []
    for endpoint, endpoint_margin in ((left, margins[0]), (right, margins[-1])):
        if endpoint_margin > 0:
            bounds.append(dict(pressure=float(endpoint), criterion=None, theta=None))
            continue
        pressure = endpoint if endpoint_margin == 0 else brentq(
            margin, *sorted((endpoint, best))
        )
        wall_margins = evaluate([pressure])[:, 0, :]
        criterion, angle = np.unravel_index(np.argmin(wall_margins), wall_margins.shape)
        bounds.append(dict(pressure=float(pressure),
                           criterion=("tensile", "shear")[criterion],
                           theta=float(theta[angle])))
    return tuple(bounds)
