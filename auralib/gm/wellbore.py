"""Elastic stresses at the wall of a circular borehole."""

import numpy as np


def borehole_wall_stress(stress, theta, poisson, mud_pressure, pore_pressure):
    """Return the effective wall stress tensor in (radial, hoop, axial) order.

    Supply the total far-field stress in borehole Cartesian coordinates
    as a symmetric array of shape (..., 3, 3). Stress and pressure use
    consistent units and compression is positive. Theta is in degrees
    from borehole +x toward +y. All leading dimensions broadcast with
    theta, Poisson's ratio, mud pressure, and formation pore pressure.

    The solution assumes a circular hole in an infinite, isotropic,
    linear elastic medium, with fixed formation pore pressure at the
    wall and a unit Biot coefficient. The result has shape (..., 3, 3).
    """
    stress = np.asarray(stress, dtype=float)
    if stress.shape[-2:] != (3, 3):
        raise ValueError("stress must have shape (..., 3, 3)")

    theta = np.deg2rad(np.asarray(theta, dtype=float))
    poisson = np.asarray(poisson, dtype=float)
    if not np.all((-1 < poisson) & (poisson < 0.5)):
        raise ValueError("poisson must be greater than -1 and less than 0.5")

    mud_pressure = np.asarray(mud_pressure, dtype=float)
    pore_pressure = np.asarray(pore_pressure, dtype=float)
    radial = mud_pressure - pore_pressure
    xx, yy, zz = stress[..., 0, 0], stress[..., 1, 1], stress[..., 2, 2]
    xy, xz, yz = stress[..., 0, 1], stress[..., 0, 2], stress[..., 1, 2]
    variation = 2 * (xx - yy) * np.cos(2 * theta) + 4 * xy * np.sin(2 * theta)

    hoop = xx + yy - 2 * pore_pressure - variation - radial
    axial = zz - pore_pressure - poisson * variation
    shear = 2 * (yz * np.cos(theta) - xz * np.sin(theta))

    radial, hoop, axial, shear = np.broadcast_arrays(radial, hoop, axial, shear)
    wall = np.zeros(radial.shape + (3, 3))
    wall[..., 0, 0] = radial
    wall[..., 1, 1] = hoop
    wall[..., 2, 2] = axial
    wall[..., 1, 2] = wall[..., 2, 1] = shear
    return wall


def wall_principal_stresses(wall_stress):
    """Return effective principal stresses in descending compressive order.

    Supply a symmetric wall tensor of shape (..., 3, 3). The final
    dimension of the result holds sigma1, sigma2, and sigma3.
    """
    wall_stress = np.asarray(wall_stress, dtype=float)
    if wall_stress.shape[-2:] != (3, 3):
        raise ValueError("wall_stress must have shape (..., 3, 3)")
    return np.linalg.eigvalsh(wall_stress)[..., ::-1]
