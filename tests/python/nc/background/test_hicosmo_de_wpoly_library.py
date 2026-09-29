#!/usr/bin/env python
#
# test_hicosmo_de_wpoly_library.py
#
# Mon Sep 28 16:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_hicosmo_de_wpoly_library.py
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

"""Unit tests for NcHICosmoDEWPolyLibrary."""

import numpy as np
from numpy.testing import assert_allclose
import pytest

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

Z_F = 2.33
N_FIDUCIALS = 40


@pytest.fixture(name="lib", scope="module")
def fixture_lib() -> Nc.HICosmoDEWPolyLibrary:
    """A small library."""
    return Nc.HICosmoDEWPolyLibrary.new(
        seed=7, degree=6, z_f=Z_F, n_fiducials=N_FIDUCIALS
    )


def test_sizes(lib: Nc.HICosmoDEWPolyLibrary) -> None:
    """Every number of turning points has n_fiducials fiducials."""
    assert lib.get_max_turning_points() == 3
    for m in range(4):
        assert lib.len(m) == N_FIDUCIALS
    assert lib.get_n_draws() > 4 * N_FIDUCIALS


def test_count_and_positions(lib: Nc.HICosmoDEWPolyLibrary) -> None:
    """The stored count and positions of turning points match a recount on the fiducial."""
    for m in range(4):
        for i in range(lib.len(m)):
            cosmo = lib.get_fiducial(m, i)
            assert isinstance(cosmo, Nc.HICosmoDEWPoly)
            assert cosmo.get_ncoeffs() == 8
            assert cosmo.get_z_f() == Z_F
            tp = cosmo.turning_points()
            stored = lib.peek_turning_points(m, i)
            assert len(tp) == m
            assert len(stored) == m
            assert_allclose(tp, stored, rtol=1.0e-12, atol=1.0e-14)
            assert all(0.0 < z < Z_F for z in tp)
            assert np.all(np.diff(stored) > 0.0) if m > 1 else True


def test_amplitude_and_mean(lib: Nc.HICosmoDEWPolyLibrary) -> None:
    """The largest excursion equals the amplitude; the mean of w over the range is -1."""
    for m in range(4):
        for i in range(lib.len(m)):
            cosmo = lib.get_fiducial(m, i)
            amp = lib.get_amplitude(m, i)
            assert 0.0 <= amp <= 0.5
            stationary = [0.0] + list(lib.peek_turning_points(m, i)) + [Z_F]
            w = np.array([cosmo.w_de(z) for z in stationary])
            assert np.max(np.abs(np.diff(w))) == pytest.approx(amp, rel=1.0e-12)
            assert lib.peek_coeffs(m, i).get(0) == 0.0


def test_reproducible() -> None:
    """The same seed gives the same draw indices and coefficients."""
    lib1 = Nc.HICosmoDEWPolyLibrary.new(seed=11, degree=6, z_f=Z_F, n_fiducials=10)
    lib2 = Nc.HICosmoDEWPolyLibrary.new(seed=11, degree=6, z_f=Z_F, n_fiducials=10)
    lib3 = Nc.HICosmoDEWPolyLibrary.new(seed=12, degree=6, z_f=Z_F, n_fiducials=10)
    for m in range(4):
        for i in range(10):
            assert lib1.get_draw_index(m, i) == lib2.get_draw_index(m, i)
            assert_allclose(
                lib1.peek_coeffs(m, i).dup_array(), lib2.peek_coeffs(m, i).dup_array()
            )
    assert any(
        lib1.get_draw_index(m, i) != lib3.get_draw_index(m, i)
        for m in range(4)
        for i in range(10)
    )


def test_prefix_is_stable() -> None:
    """Draw indices are a property of the stream, not of n_fiducials."""
    small = Nc.HICosmoDEWPolyLibrary.new(seed=3, degree=6, z_f=Z_F, n_fiducials=10)
    large = Nc.HICosmoDEWPolyLibrary.new(seed=3, degree=6, z_f=Z_F, n_fiducials=30)
    small_draws = {(m, small.get_draw_index(m, i)) for m in range(4) for i in range(10)}
    large_draws = {(m, large.get_draw_index(m, i)) for m in range(4) for i in range(30)}
    # every monotonic draw is accepted, so the first 10 of them are shared
    small0 = sorted(d for m, d in small_draws if m == 0)
    large0 = sorted(d for m, d in large_draws if m == 0)
    assert small0 == large0[:10]


def test_properties() -> None:
    """Construct-only properties are kept and serialization round trips."""
    lib = Nc.HICosmoDEWPolyLibrary(
        seed=5, degree=4, zf=1.5, n_fiducials=5, max_turning_points=2, a_max=0.2
    )
    assert lib.get_max_turning_points() == 2
    for m in range(3):
        assert lib.len(m) == 5
        for i in range(5):
            assert lib.get_amplitude(m, i) <= 0.2
            assert lib.get_fiducial(m, i).get_ncoeffs() == 6
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    lib2 = ser.dup_obj(lib)
    assert isinstance(lib2, Nc.HICosmoDEWPolyLibrary)
    assert lib2.get_draw_index(2, 3) == lib.get_draw_index(2, 3)
