#!/usr/bin/env python
#
# test_halo_position.py
#
# Sun Nov 15 2026
# Copyright  2026  Caio Lima de Oliveira
# <caiolimadeoliveira@pm.me>
#
# test_halo_position.py
# Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

"""Tests on NcHaloPosition from Python."""

import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.mark.parametrize(
    "ra_h, dec_h, ra, dec",
    [(0.0, 0.0, 0.03, 0.02), (35.2, -4.1, 35.0, -4.3), (210.0, 60.0, 211.0, 60.5)],
)
def test_rotate_in_place(ra_h, dec_h, ra, dec):
    """NcHaloPosition.rotate() rotates the vector it is given in place: the halo goes
    to the z axis, and the spherical angles of a rotated position are those of
    polar_angles()."""
    hp = Nc.HaloPosition.new(Nc.Distance.new(2.0))
    hp.param_set_by_name("ra", ra_h)
    hp.param_set_by_name("dec", dec_h)

    h = Ncm.TriVec.new_astro_ra_dec(1.0, ra_h, dec_h)
    assert hp.rotate(h) is None
    assert_allclose(h.c, [0.0, 0.0, 1.0], atol=1.0e-14)

    v = Ncm.TriVec.new_astro_ra_dec(1.0, ra, dec)
    hp.rotate(v)
    _, theta_v, phi_v = v.get_spherical_coord()
    theta, phi = hp.polar_angles(ra, dec)

    assert theta_v == theta
    assert phi_v == phi
