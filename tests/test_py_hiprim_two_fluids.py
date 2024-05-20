#!/usr/bin/env python
#
# test_py_hiprim_two_fluids.py
#
# Sun May 19 23:27:08 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_hiprim_two_fluids.py
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

"""Tests on NcHIPrimTwoFluids cosmological model."""

import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="two_fluids")
def fixture_two_fluids() -> Nc.HIPrimTwoFluids:
    """Fixture for NcHIPrimTwoFluids."""
    two_fluids = Nc.HIPrimTwoFluids.new()

    return two_fluids


def test_init(two_fluids):
    """Test NcHIPrimTwoFluids initialization."""
    assert two_fluids is not None
    assert isinstance(two_fluids, Nc.HIPrimTwoFluids)
    assert isinstance(two_fluids, Nc.HIPrim)
    assert isinstance(two_fluids, Ncm.Model)
    assert two_fluids.get_use_default_calib() is False


def test_set_use_default_calib(two_fluids):
    """Test NcHIPrimTwoFluids set_use_default_calib."""
    two_fluids.set_use_default_calib(True)
    assert two_fluids.get_use_default_calib() is True
    assert two_fluids.peek_lnk_lnw_spline() is None

    two_fluids.set_use_default_calib(False)
    assert two_fluids.get_use_default_calib() is False
    assert two_fluids.peek_lnk_lnw_spline() is None

    two_fluids.set_use_default_calib(True)
    assert two_fluids.get_use_default_calib() is True
    assert two_fluids.peek_lnk_lnw_spline() is None


def test_set_lnk_lnw_spline_init(two_fluids):
    """Test NcHIPrimTwoFluids set_lnk_lnw_spline."""
    lnk = np.log(np.geomspace(1.0e-5, 1.0e1, 1000))
    lnw = np.log(np.geomspace(1.0e-7, 1.0e-1, 330))

    lnk_v, lnw_v = np.meshgrid(lnk, lnw)

    lnPk = 2.0 * lnk_v + lnw_v

    lnPk2d = Ncm.Spline2dBicubic(
        init=True,
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(lnk.tolist()),
        y_vector=Ncm.Vector.new_array(lnw.tolist()),
        z_matrix=Ncm.Matrix.new_array(lnPk.flatten().tolist(), len(lnk)),
    )

    two_fluids.set_lnk_lnw_spline(lnPk2d)
    assert two_fluids.peek_lnk_lnw_spline() is not None
    assert two_fluids.peek_lnk_lnw_spline() is lnPk2d

    eval_lnSA = np.array([two_fluids.lnSA_powspec_lnk(lnki) for lnki in lnk])

    assert_allclose(
        eval_lnSA,
        two_fluids.props.ln10e10ASA
        - 10.0 * np.log(10.0)
        + 2.0 * (lnk - two_fluids.props.lnk0)
        + two_fluids.props.lnw,
        rtol=1.0e-8,
    )


def test_set_lnk_lnw_spline_no_init(two_fluids):
    """Test NcHIPrimTwoFluids set_lnk_lnw_spline."""
    lnk = np.log(np.geomspace(1.0e-5, 1.0e1, 1000))
    lnw = np.log(np.geomspace(1.0e-7, 1.0e-1, 330))

    lnk_v, lnw_v = np.meshgrid(lnk, lnw)

    lnPk = 2.0 * lnk_v + lnw_v

    lnPk2d = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(lnk.tolist()),
        y_vector=Ncm.Vector.new_array(lnw.tolist()),
        z_matrix=Ncm.Matrix.new_array(lnPk.flatten().tolist(), len(lnk)),
    )

    two_fluids.set_lnk_lnw_spline(lnPk2d)
    assert two_fluids.peek_lnk_lnw_spline() is not None
    assert two_fluids.peek_lnk_lnw_spline() is lnPk2d

    eval_lnSA = np.array([two_fluids.lnSA_powspec_lnk(lnki) for lnki in lnk])

    assert_allclose(
        eval_lnSA,
        two_fluids.props.ln10e10ASA
        - 10.0 * np.log(10.0)
        + 2.0 * (lnk - two_fluids.props.lnk0)
        + two_fluids.props.lnw,
        rtol=1.0e-8,
    )


def test_set_lnk_lnw_spline_serialize(two_fluids):
    """Test NcHIPrimTwoFluids set_lnk_lnw_spline."""
    lnk = np.log(np.geomspace(1.0e-5, 1.0e1, 1000))
    lnw = np.log(np.geomspace(1.0e-7, 1.0e-1, 330))

    lnk_v, lnw_v = np.meshgrid(lnk, lnw)

    lnPk = 2.0 * lnk_v + lnw_v

    lnPk2d = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(lnk.tolist()),
        y_vector=Ncm.Vector.new_array(lnw.tolist()),
        z_matrix=Ncm.Matrix.new_array(lnPk.flatten().tolist(), len(lnk)),
    )

    two_fluids.set_lnk_lnw_spline(lnPk2d)
    assert two_fluids.peek_lnk_lnw_spline() is not None
    assert two_fluids.peek_lnk_lnw_spline() is lnPk2d

    eval_lnSA = np.array([two_fluids.lnSA_powspec_lnk(lnki) for lnki in lnk])

    assert_allclose(
        eval_lnSA,
        two_fluids.props.ln10e10ASA
        - 10.0 * np.log(10.0)
        + 2.0 * (lnk - two_fluids.props.lnk0)
        + two_fluids.props.lnw,
        rtol=1.0e-8,
    )

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    two_fluids2 = ser.dup_obj(two_fluids)

    assert two_fluids2.peek_lnk_lnw_spline() is not None
    assert two_fluids2.peek_lnk_lnw_spline() is not lnPk2d

    eval_lnSA = np.array([two_fluids2.lnSA_powspec_lnk(lnki) for lnki in lnk])

    assert_allclose(
        eval_lnSA,
        two_fluids.props.ln10e10ASA
        - 10.0 * np.log(10.0)
        + 2.0 * (lnk - two_fluids.props.lnk0)
        + two_fluids.props.lnw,
        rtol=1.0e-8,
    )


def test_lnT_powespec_lnk(two_fluids):
    """Test the tensorial power spectrum."""
    lnk = np.log(np.geomspace(1.0e-5, 1.0e1, 1000))

    eval_lnT = [two_fluids.lnT_powspec_lnk(lnki) for lnki in lnk]

    assert_allclose(
        eval_lnT,
        two_fluids.props.ln10e10ASA
        - 10.0 * np.log(10.0)
        + np.log(two_fluids.props.T_SA_ratio)
        + two_fluids.props.n_T * (lnk - np.log(two_fluids.props.k_pivot)),
        rtol=1.0e-8,
    )
