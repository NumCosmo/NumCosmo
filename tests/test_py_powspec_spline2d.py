#!/usr/bin/env python
#
# test_py_powspec_ml_spline.py
#
# Fri Mar 22 13:44:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_powspec_ml_spline.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Unit tests for NumCosmo powwer-spectra."""

import pytest

import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm

Ncm.cfg_init()


def _f(z, k):
    return 1.0e-9 * (k / 0.05) ** (0.96 - 1.0) * (1.0 + z)


@pytest.fixture(name="Pk2d")
def fixture_Pk2d() -> Ncm.Spline2d:
    """Fixture for k array."""
    ka = np.geomspace(1.0e-5, 1.0e1, 1000)
    za = np.linspace(0.0, 1.0, 50)

    za_v, ka_v = np.meshgrid(za, ka)

    Pk = np.log(_f(za_v, ka_v))

    Pk2d = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(za.tolist()),
        y_vector=Ncm.Vector.new_array(np.log(ka).tolist()),
        z_matrix=Ncm.Matrix.new_array(Pk.flatten().tolist(), len(za)),
    )

    return Pk2d


def test_constructor(Pk2d: Ncm.Spline2d) -> None:
    """Test the constructor of the power spectrum."""
    assert Pk2d is not None

    ps = Ncm.PowspecSpline2d.new(Pk2d)
    assert ps is not None
    assert isinstance(ps, Ncm.PowspecSpline2d)
    assert isinstance(ps, Ncm.Powspec)

    assert Pk2d == ps.peek_spline2d()


def test_eval(Pk2d: Ncm.Spline2d) -> None:
    """Test the evaluation of the power spectrum."""
    ps = Ncm.PowspecSpline2d.new(Pk2d)

    ps.prepare()

    assert_allclose(ps.eval(None, 0.0, 1.0), _f(0.0, 1.0))

    s2d = ps.peek_spline2d()

    za = s2d.peek_xv().dup_array()
    ka = np.exp(s2d.peek_yv().dup_array())

    ps_Pka = np.array([[ps.eval(None, z, k) for k in ka] for z in za])
    Pka = np.array([np.exp([s2d.eval(z, np.log(k)) for k in ka]) for z in za])
    fa = np.array([[_f(z, k) for k in ka] for z in za])

    assert_allclose(ps_Pka, Pka)
    assert_allclose(ps_Pka, fa)


def test_eval_vec(Pk2d: Ncm.Spline2d) -> None:
    """Test the evaluation of the power spectrum."""
    ps = Ncm.PowspecSpline2d.new(Pk2d)

    ps.prepare()

    assert_allclose(ps.eval(None, 0.0, 1.0), _f(0.0, 1.0))

    s2d = ps.peek_spline2d()

    za = s2d.peek_xv().dup_array()
    ka = np.exp(s2d.peek_yv().dup_array())

    kv = Ncm.Vector.new_array(ka.tolist())
    Pkv = kv.dup()

    def _evec(z):
        ps.eval_vec(None, z, kv, Pkv)
        return Pkv.dup_array()

    ps_Pka = np.array([_evec(z) for z in za])
    Pka = np.array([np.exp([s2d.eval(z, np.log(k)) for k in ka]) for z in za])
    fa = np.array([[_f(z, k) for k in ka] for z in za])

    assert_allclose(ps_Pka, Pka)
    assert_allclose(ps_Pka, fa)


def test_get_spline2d(Pk2d: Ncm.Spline2d) -> None:
    """Test the getter of the power spectrum."""
    ps = Ncm.PowspecSpline2d.new(Pk2d)

    assert ps is not None

    s2d = ps.peek_spline2d()

    assert s2d is not None
    assert s2d == Pk2d

    sd2 = ps.get_spline_2d()

    assert sd2 is not None
    assert sd2 == Pk2d


def test_serialization(Pk2d: Ncm.Spline2d) -> None:
    """Test the serialization of the power spectrum."""
    ps = Ncm.PowspecSpline2d.new(Pk2d)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    ps_dup = ser.dup_obj(ps)

    assert ps_dup is not None
    assert isinstance(ps_dup, Ncm.PowspecSpline2d)

    ps.prepare()
    ps_dup.prepare()

    zv = Pk2d.peek_xv()
    kv = Pk2d.peek_yv()
    za = zv.dup_array()
    Pkv = kv.dup()

    ps_Pka = []
    ps_dup_Pka = []
    for z in za:
        ps.eval_vec(None, z, kv, Pkv)
        ps_Pka.append(Pkv.dup_array())
        ps_dup.eval_vec(None, z, kv, Pkv)
        ps_dup_Pka.append(Pkv.dup_array())

    assert_allclose(ps_Pka, ps_dup_Pka)


def test_get_nknots(Pk2d: Ncm.Spline2d) -> None:
    """Test the getter of the number of knots."""
    ps = Ncm.PowspecSpline2d.new(Pk2d)

    assert ps is not None

    s2d = ps.peek_spline2d()

    assert s2d is not None
    assert s2d == Pk2d

    nknots = ps.get_nknots()

    assert nknots == (Pk2d.peek_xv().len(), Pk2d.peek_yv().len())


def test_extrapolation(Pk2d: Ncm.Spline2d) -> None:
    """Test the evaluation of the power spectrum."""
    ps = Ncm.PowspecSpline2d.new(Pk2d)

    ps.prepare()

    kmin = ps.get_kmin()
    kmax = ps.get_kmax()

    # Test extrapolation, decreasing k, must be smaller than the value at kmin
    for k_small in np.geomspace(kmin * 1.0e-5, kmin * 0.99, 100):
        assert ps.eval(None, 0.0, k_small) < ps.eval(None, 0.0, kmin)
    # Test extrapolation, increasing k, must be smaller than the value at kmax
    for k_large in np.geomspace(kmax * 1.01, kmax * 1.0e5, 100):
        assert ps.eval(None, 0.0, k_large) < ps.eval(None, 0.0, kmax)
