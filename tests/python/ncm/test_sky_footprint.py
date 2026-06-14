#!/usr/bin/env python
#
# test_sky_footprint.py
#
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""Tests for the NcmSkyFootprint rectangular sky region."""

import math
import hashlib

import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


def _rect() -> Ncm.SkyFootprintRectangular:
    """A rectangular footprint over [10, 40] x [-5, 25] degrees."""
    return Ncm.SkyFootprintRectangular.new(10.0, 40.0, -5.0, 25.0)


def test_limits_roundtrip() -> None:
    """The RA/DEC limits are stored and read back."""
    rect = _rect()
    ra_min, ra_max = rect.get_ra_lim()
    dec_min, dec_max = rect.get_dec_lim()
    assert (ra_min, ra_max) == (10.0, 40.0)
    assert (dec_min, dec_max) == (-5.0, 25.0)


def test_area_matches_analytic() -> None:
    """The solid angle matches the spherical-rectangle area."""
    rect = _rect()
    expected = math.radians(40.0 - 10.0) * (
        math.sin(math.radians(25.0)) - math.sin(math.radians(-5.0))
    )
    assert rect.get_area() == expected
    assert isinstance(rect, Ncm.SkyFootprint)


def test_contains() -> None:
    """Membership respects the rectangle boundaries."""
    rect = _rect()
    assert rect.contains(25.0, 10.0)
    assert not rect.contains(100.0, 0.0)
    assert not rect.contains(25.0, 80.0)


def test_density_and_ln_density() -> None:
    """The density is cos(dec)-weighted inside and zero/-inf outside."""
    rect = _rect()
    inside = rect.density(25.0, 10.0)
    assert inside > 0.0
    assert math.isclose(rect.ln_density(25.0, 10.0), math.log(inside), rel_tol=1e-15)
    assert rect.density(100.0, 0.0) == 0.0
    assert rect.ln_density(100.0, 0.0) == -math.inf


def test_density_normalization() -> None:
    """The density integrates to one over the rectangle (grid quadrature)."""
    rect = _rect()
    ra = np.linspace(10.0, 40.0, 400)
    dec = np.linspace(-5.0, 25.0, 400)
    grid = np.array([[rect.density(r, d) for r in ra] for d in dec])
    integral = np.trapezoid(np.trapezoid(grid, ra, axis=1), dec)
    assert math.isclose(integral, 1.0, rel_tol=1e-4)


def test_gen_within_bounds_and_reproducible() -> None:
    """Sampling stays in the rectangle and is deterministic for a fixed seed."""
    rect = _rect()
    rng = Ncm.RNG.seeded_new(None, 123)
    ras = np.empty(5000)
    decs = np.empty(5000)
    for i in range(5000):
        ras[i], decs[i] = rect.gen_ra_dec(rng)

    assert ras.min() >= 10.0 and ras.max() <= 40.0
    assert decs.min() >= -5.0 and decs.max() <= 25.0

    digest = hashlib.sha256(np.concatenate([ras, decs]).tobytes()).hexdigest()[:16]
    assert digest == "419ea91d5c40d247"


def test_serialization_roundtrip() -> None:
    """The footprint survives a serialization round-trip."""
    rect = _rect()
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    dup = ser.from_string(ser.to_string(rect, True))

    assert isinstance(dup, Ncm.SkyFootprintRectangular)
    assert dup.get_ra_lim() == (10.0, 40.0)
    assert dup.get_dec_lim() == (-5.0, 25.0)
    assert dup.get_area() == rect.get_area()
