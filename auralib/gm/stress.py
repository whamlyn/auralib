"""Stress tensor calculations for geomechanics."""

import numpy as np

from .coordinates import borehole_rotation_matrix


def principal_stress_to_ned(shmax, shmin, sv, shmax_azimuth):
    """Build a North-East-Down stress tensor from principal stresses.

    Supply the maximum and minimum horizontal stresses (shmax and shmin)
    and vertical stress (sv) in consistent units and with a common sign
    convention. The shmax_azimuth is in degrees clockwise from north.
    Shmin is perpendicular to Shmax in the horizontal plane. Azimuths
    180 degrees apart describe the same stress orientation.

    The inputs can be scalars or broadcastable arrays. The returned
    tensor has shape (..., 3, 3) and preserves the stress units and signs.
    """
    shmax, shmin, sv, shmax_azimuth = np.broadcast_arrays(
        np.asarray(shmax, dtype=float),
        np.asarray(shmin, dtype=float),
        np.asarray(sv, dtype=float),
        np.asarray(shmax_azimuth, dtype=float),
    )
    if not np.all(np.isfinite(shmax_azimuth)):
        raise ValueError("shmax_azimuth must be finite")

    azimuth = np.deg2rad(np.remainder(shmax_azimuth, 180))
    cos_az, sin_az = np.cos(azimuth), np.sin(azimuth)

    stress = np.zeros(shmax.shape + (3, 3))
    stress[..., 0, 0] = shmax * cos_az**2 + shmin * sin_az**2
    stress[..., 1, 1] = shmax * sin_az**2 + shmin * cos_az**2
    stress[..., 0, 1] = stress[..., 1, 0] = (shmax - shmin) * cos_az * sin_az
    stress[..., 2, 2] = sv
    return stress


def rotate_stress_to_borehole(stress, azimuth, inclination):
    """Express a far-field stress tensor in borehole Cartesian coordinates.

    Supply the full tensor, including shear components, in North-East-Down
    order with shape (..., 3, 3). Use a diagonal tensor only when the
    principal stress axes align with North-East-Down. Use
    principal_stress_to_ned to construct the tensor from horizontal and
    vertical principal stresses and a Shmax azimuth. Azimuth and inclination
    follow the degree conventions in borehole_rotation_matrix. Leading stress
    dimensions and angle arrays broadcast against each other.

    Returns R @ stress @ R.T for each trajectory orientation, with axes
    ordered high side, transverse, and downhole. The rotation preserves the
    input stress units and sign convention.
    """
    stress = np.asarray(stress, dtype=float)
    if stress.shape[-2:] != (3, 3):
        raise ValueError("stress must have shape (..., 3, 3)")

    rotation = borehole_rotation_matrix(azimuth, inclination)
    return rotation @ stress @ rotation.swapaxes(-1, -2)
