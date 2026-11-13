"""Coordinate transforms for geomechanics."""

import numpy as np


def borehole_rotation_matrix(azimuth, inclination):
    """Return the rotation matrix from North-East-Down to borehole coordinates.

    Azimuth is clockwise from north, and inclination is measured from
    vertical down (0 degrees vertical, 90 degrees horizontal). Both angles
    are in degrees. Angles can be scalars or broadcastable arrays.

    The rows are the borehole axes expressed in North-East-Down coordinates.
    The z-axis points downhole, the x-axis points toward the high side, and
    the y-axis completes the right-handed frame. For vertical wells, where
    the high side is undefined, the supplied azimuth fixes the transverse
    axes by continuity. The result has shape (..., 3, 3).
    """
    azimuth, inclination = np.broadcast_arrays(
        np.asarray(azimuth, dtype=float), np.asarray(inclination, dtype=float)
    )
    if not np.all(np.isfinite(azimuth)):
        raise ValueError("azimuth must be finite")
    if not np.all((0 <= inclination) & (inclination <= 180)):
        raise ValueError("inclination must be between 0 and 180 degrees")

    azimuth = np.deg2rad(np.remainder(azimuth, 360))
    inclination = np.deg2rad(inclination)
    cos_az, sin_az = np.cos(azimuth), np.sin(azimuth)
    cos_inc, sin_inc = np.cos(inclination), np.sin(inclination)

    x = np.stack((cos_inc * cos_az, cos_inc * sin_az, -sin_inc), axis=-1)
    y = np.stack((-sin_az, cos_az, np.zeros_like(azimuth)), axis=-1)
    z = np.stack((sin_inc * cos_az, sin_inc * sin_az, cos_inc), axis=-1)
    return np.stack((x, y, z), axis=-2)
