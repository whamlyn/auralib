"""Failure criteria for compression-positive effective stresses."""

import numpy as np


def tensile_margin(sigma3, tensile_strength):
    """Return sigma3 + tensile_strength, with zero at tensile failure onset.

    Tensile strength is a nonnegative magnitude in the same units as sigma3.
    A negative margin indicates that the tensile cutoff is exceeded.
    """
    tensile_strength = np.asarray(tensile_strength, dtype=float)
    if not np.all(np.isfinite(tensile_strength) & (tensile_strength >= 0)):
        raise ValueError("tensile_strength must be finite and nonnegative")
    return np.asarray(sigma3, dtype=float) + tensile_strength


def mohr_coulomb_margin(sigma1, sigma3, ucs, friction_angle):
    """Return the shear strength minus sigma1, with zero at failure onset.

    Supply UCS and effective principal stresses in consistent units.
    Friction angle is in degrees. A negative margin indicates that the
    Mohr-Coulomb criterion is exceeded.
    """
    ucs = np.asarray(ucs, dtype=float)
    friction_angle = np.asarray(friction_angle, dtype=float)
    if not np.all(np.isfinite(ucs) & (ucs >= 0)):
        raise ValueError("ucs must be finite and nonnegative")
    if not np.all((0 <= friction_angle) & (friction_angle < 90)):
        raise ValueError("friction_angle must be in [0, 90) degrees")

    sine = np.sin(np.deg2rad(friction_angle))
    confinement_factor = (1 + sine) / (1 - sine)
    return ucs + confinement_factor * np.asarray(sigma3, dtype=float) - sigma1
