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
import numpy.typing as npt

from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import npa_to_seq

Ncm.cfg_init()


@pytest.fixture(name="Pk")
def fixture_Pk() -> Ncm.Spline:
    """Fixture for k array."""
    ka = np.geomspace(1.0e-5, 1.0e1, 1000, dtype=np.float64)
    Pk: npt.NDArray[np.float64] = np.array(
        1.0e-9 * (ka / 0.05) ** (0.96 - 1.0), dtype=np.float64
    )

    Pk_spline = Ncm.SplineCubicNotaknot.new()
    Pk_spline.set_array(npa_to_seq(ka), npa_to_seq(Pk), True)

    return Pk_spline


def test_constructor(Pk: Ncm.Spline) -> None:
    """Test the constructor of the power spectrum."""
    assert Pk is not None

    ps = Nc.PowspecMLSpline.new(Pk)
    assert ps is not None

    assert Pk == ps.peek_spline()


def test_eval_z0(Pk: Ncm.Spline) -> None:
    """Test the evaluation of the power spectrum."""
    ps = Nc.PowspecMLSpline.new(Pk)
    assert ps is not None

    cosmo = Nc.HICosmoDEXcdm()
    prim = Nc.HIPrimPowerLaw()
    cosmo.add_submodel(prim)

    ps.prepare(cosmo)

    ka = Pk.peek_xv().dup_array()

    ps_Pka = [ps.eval(cosmo, 0, k) for k in ka]
    Pka = [Pk.eval(k) for k in ka]

    assert_allclose(ps_Pka, Pka)


def test_eval_vec_z0(Pk: Ncm.Spline) -> None:
    """Test the evaluation of the power spectrum."""
    ps = Nc.PowspecMLSpline.new(Pk)
    assert ps is not None

    cosmo = Nc.HICosmoDEXcdm()
    prim = Nc.HIPrimPowerLaw()
    cosmo.add_submodel(prim)

    ps.prepare(cosmo)

    kv = Pk.peek_xv()
    ka = Pk.peek_xv().dup_array()
    Pkv = kv.dup()

    ps.eval_vec(cosmo, 0, kv, Pkv)
    ps_Pka = Pkv.dup_array()
    Pka = [Pk.eval(k) for k in ka]

    assert_allclose(ps_Pka, Pka)


def test_eval_za(Pk: Ncm.Spline) -> None:
    """Test the evaluation of the power spectrum."""
    ps = Nc.PowspecMLSpline.new(Pk)
    assert ps is not None

    cosmo = Nc.HICosmoDEXcdm()
    prim = Nc.HIPrimPowerLaw()
    cosmo.add_submodel(prim)

    gf = Nc.GrowthFunc()

    ps.prepare(cosmo)
    gf.prepare(cosmo)

    ka = Pk.peek_xv().dup_array()
    za = np.linspace(0, 1, 50)

    ps_Pka = [[ps.eval(cosmo, z, k) for k in ka] for z in za]
    Pka = [[Pk.eval(k) * gf.eval(cosmo, z) ** 2 for k in ka] for z in za]

    assert_allclose(ps_Pka, Pka)


def test_eval_vec_za(Pk: Ncm.Spline) -> None:
    """Test the evaluation of the power spectrum."""
    ps = Nc.PowspecMLSpline.new(Pk)
    assert ps is not None

    cosmo = Nc.HICosmoDEXcdm()
    prim = Nc.HIPrimPowerLaw()
    cosmo.add_submodel(prim)

    gf = Nc.GrowthFunc()

    ps.prepare(cosmo)
    gf.prepare(cosmo)

    kv = Pk.peek_xv()
    ka = Pk.peek_xv().dup_array()
    Pkv = kv.dup()

    za = np.linspace(0, 1, 50)

    ps_Pka = []
    for z in za:
        ps.eval_vec(cosmo, z, kv, Pkv)
        ps_Pka.append(Pkv.dup_array())

    Pka = [[Pk.eval(k) * gf.eval(cosmo, z) ** 2 for k in ka] for z in za]

    assert_allclose(ps_Pka, Pka)


def test_get_knots(Pk: Ncm.Spline) -> None:
    """Test the evaluation of the power spectrum."""
    ps = Nc.PowspecMLSpline.new(Pk)
    assert ps is not None

    zi = ps.get_zi()
    zf = ps.get_zf()

    assert zf > zi

    assert ps.get_nknots() == ((zf - zi) * 50, Pk.get_len())


def test_serialization(Pk: Ncm.Spline) -> None:
    """Test the serialization of the power spectrum."""
    ps = Nc.PowspecMLSpline.new(Pk)
    assert ps is not None

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    ps_dup = ser.dup_obj(ps)

    assert ps_dup is not None
    assert isinstance(ps_dup, Nc.PowspecMLSpline)

    cosmo = Nc.HICosmoDEXcdm()
    prim = Nc.HIPrimPowerLaw()
    cosmo.add_submodel(prim)

    ps.prepare(cosmo)
    ps_dup.prepare(cosmo)

    kv = Pk.peek_xv()
    Pkv = kv.dup()

    za = np.linspace(0, 1, 50)
    ps_Pka = []
    ps_dup_Pka = []
    for z in za:
        ps.eval_vec(cosmo, z, kv, Pkv)
        ps_Pka.append(Pkv.dup_array())
        ps_dup.eval_vec(cosmo, z, kv, Pkv)
        ps_dup_Pka.append(Pkv.dup_array())

    assert_allclose(ps_Pka, ps_dup_Pka)
