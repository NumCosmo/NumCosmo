"""Geometry helpers for mock catalogs.

These are thin wrappers over NumCosmo's geometry primitives (:class:`Ncm.TriVec`)
plus the analytic solid angle of a sky rectangle. They keep a single source of
truth for the coordinate convention shared with the C library and are the natural
seam for a future move into C.
"""

import numpy as np

from numcosmo_py import Ncm


def spherical_to_cartesian(ra, dec, r):
    """Convert sky coordinates to Cartesian using NumCosmo's astro convention.

    :param ra: right ascension in degrees (scalar or array-like).
    :param dec: declination in degrees (scalar or array-like).
    :param r: radial distance (scalar or array-like).
    :return: tuple ``(x1, x2, x3)`` matching the inputs' broadcast shape (scalars
        for scalar inputs).

    Delegates to :func:`Ncm.TriVec.new_astro_ra_dec` per point so the convention
    is defined once, in C.
    """
    ra_a = np.asarray(ra, dtype=float)
    dec_a = np.asarray(dec, dtype=float)
    r_a = np.asarray(r, dtype=float)
    scalar = ra_a.ndim == 0 and dec_a.ndim == 0 and r_a.ndim == 0

    ra_b, dec_b, r_b = np.broadcast_arrays(ra_a, dec_a, r_a)
    x1 = np.empty(ra_b.shape, dtype=float)
    x2 = np.empty(ra_b.shape, dtype=float)
    x3 = np.empty(ra_b.shape, dtype=float)
    for idx in np.ndindex(ra_b.shape):
        coords = Ncm.TriVec.new_astro_ra_dec(
            float(r_b[idx]), float(ra_b[idx]), float(dec_b[idx])
        ).c
        x1[idx], x2[idx], x3[idx] = coords[0], coords[1], coords[2]

    if scalar:
        return float(x1), float(x2), float(x3)
    return x1, x2, x3


def rectangular_sky_area(
    ra_min: float, ra_max: float, dec_min: float, dec_max: float
) -> float:
    """Solid angle (square degrees) of an RA/DEC rectangle.

    The region is bounded by meridians and parallels (constant RA and DEC), which
    matches a catalog sampled uniformly in RA and ``sin(DEC)``. The exact solid
    angle is ``ΔRA · (sin δ_max − sin δ_min)``.
    """
    delta_ra = np.radians(ra_max - ra_min)
    delta_sin_dec = np.sin(np.radians(dec_max)) - np.sin(np.radians(dec_min))
    return delta_ra * delta_sin_dec * (180.0 / np.pi) ** 2
