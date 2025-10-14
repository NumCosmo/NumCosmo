#!/usr/bin/env python
#
# test_py_hicosmo_qgrw.py
#
# Tue Oct 14 16:18:15 2025
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_hicosmo_qgrw.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcHICosmoQGRW cosmological model."""

from itertools import product
import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="qgrw")
def fixture_qgrw() -> Nc.HICosmoQGRW:
    """Fixture for HICosmoQGRW."""
    cosmo = Nc.HICosmoQGRW()

    return cosmo


def test_init(qgrw: Nc.HICosmoQGRW) -> None:
    """Test HICosmoQGRW initialization."""
    assert qgrw is not None
    assert isinstance(qgrw, Nc.HICosmoQGRW)
    assert isinstance(qgrw, Nc.HICosmo)
    assert isinstance(qgrw, Ncm.Model)


@pytest.fixture(name="alpha_a")
def fixture_alpha_a() -> np.ndarray:
    """Compute an array of alpha values."""
    return np.linspace(-100.0, -1.0e-1, 1000)


def test_alpha_x(qgrw: Nc.HICosmoQGRW, alpha_a: np.ndarray) -> None:
    """Test HICosmoQGRW.tau_x."""
    x_a = np.array([qgrw.x_alpha(alpha) for alpha in alpha_a])
    assert all(np.isfinite(x_a))
    assert all(x_a > 0)
    assert all(np.diff(x_a) > 0)


def test_serialize(qgrw: Nc.HICosmoQGRW):
    """Test HICosmoQGRW serialization."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    qgrw2 = ser.dup_obj(qgrw)

    assert qgrw2 is not None
    assert isinstance(qgrw2, Nc.HICosmoQGRW)
    assert isinstance(qgrw2, Nc.HICosmo)
    assert isinstance(qgrw2, Ncm.Model)

    assert qgrw2 is not qgrw
    assert qgrw2.props.H0 == qgrw.props.H0
    assert qgrw2.props.Omegar == qgrw.props.Omegar


def test_iadiab_eval_at(qgrw: Nc.HICosmoQGRW, alpha_a: np.ndarray) -> None:
    """Evaluate HICosmoQGRW implementation of NcHIPertAdiab at a given time."""
    k_a = np.geomspace(1.0e-3, 1.0e3, 5)

    unit = Nc.HIPertIAdiab.eval_unit(qgrw)
    assert np.isfinite(unit)
    assert unit > 0.0
    for alpha, k in product(alpha_a, k_a):
        F1 = Nc.HIPertIAdiab.eval_F1(qgrw, alpha, k)
        m = Nc.HIPertIAdiab.eval_m(qgrw, alpha, k)
        nu = Nc.HIPertIAdiab.eval_nu(qgrw, alpha, k)
        xi = Nc.HIPertIAdiab.eval_xi(qgrw, alpha, k)
        x = Nc.HIPertIAdiab.eval_x(qgrw, alpha)
        p2Psi = Nc.HIPertIAdiab.eval_p2Psi(qgrw, alpha, k)
        p2drho = Nc.HIPertIAdiab.eval_p2drho(qgrw, alpha, k)
        lapse = Nc.HIPertIAdiab.eval_lapse(qgrw, alpha)
        for val in (F1, m, nu, xi, x, unit, p2Psi, p2drho, lapse):
            assert np.isfinite(val)

        assert_allclose(np.log(m * nu), xi, rtol=1.0e-11, atol=0.0)


def test_igw_eval_at(qgrw: Nc.HICosmoQGRW, alpha_a: np.ndarray) -> None:
    """Evaluate HICosmoQGRW implementation of NcHIPertIGW at a given time."""
    k_a = np.geomspace(1.0e-3, 1.0e3, 5)

    unit = Nc.HIPertIGW.eval_unit(qgrw)
    assert np.isfinite(unit)
    assert unit > 0.0
    for alpha, k in product(alpha_a, k_a):
        F1 = Nc.HIPertIGW.eval_F1(qgrw, alpha, k)
        m = Nc.HIPertIGW.eval_m(qgrw, alpha, k)
        nu = Nc.HIPertIGW.eval_nu(qgrw, alpha, k)
        xi = Nc.HIPertIGW.eval_xi(qgrw, alpha, k)
        x = Nc.HIPertIGW.eval_x(qgrw, alpha)

        for val in (F1, m, nu, xi, x, unit):
            assert np.isfinite(val)

        assert_allclose(np.log(m * nu), xi, rtol=1.0e-11, atol=0.0)


def _test_oem(oem: Nc.HIPertITwoFluidsEOM) -> None:
    assert isinstance(oem, Nc.HIPertITwoFluidsEOM)
    assert oem.gw1 > 0.0
    assert oem.gw2 > 0.0
    assert oem.gw1 != oem.gw2
    assert oem.Fnu >= 0.0
    assert oem.cs2 > 0.0
    assert oem.cm2 > 0.0
    assert oem.m_zeta > 0.0
    assert oem.m_s > 0.0
    assert oem.mnu2_zeta > 0.0
    assert oem.mnu2_s > 0.0
    assert_allclose(
        oem.m_zeta * oem.Fnu**2 * oem.cs2, oem.mnu2_zeta, rtol=1.0e-11, atol=0.0
    )
    assert_allclose(oem.m_s * oem.Fnu**2 * oem.cm2, oem.mnu2_s, rtol=1.0e-11, atol=0.0)


def _test_wkb(wkb: Nc.HIPertITwoFluidsWKB) -> None:
    assert isinstance(wkb, Nc.HIPertITwoFluidsWKB)
    for val in [
        wkb.mode1_zeta_scale,
        wkb.mode2_zeta_scale,
        wkb.mode1_Q_scale,
        wkb.mode2_Q_scale,
        wkb.mode1_Pzeta_scale,
        wkb.mode2_Pzeta_scale,
        wkb.mode1_PQ_scale,
        wkb.mode2_PQ_scale,
    ]:
        assert np.isfinite(val)

    state = wkb.peek_state()
    assert isinstance(state, Nc.HIPertITwoFluidsState)
    assert state.gw1 > 0.0
    assert state.gw2 > 0.0
    assert state.gw1 != state.gw2
    assert state.Fnu >= 0.0
    assert state.norma > 0.0


def test_itwo_fluids_wkb_eval(qgrw: Nc.HICosmoQGRW, alpha_a: np.ndarray) -> None:
    """Evaluate HICosmoQGRW implementation of NcHIPertITwoFluids at a given time."""
    k_a = np.geomspace(1.0e-3, 1.0e3, 5)

    unit = Nc.HIPertITwoFluids.eval_unit(qgrw)
    assert np.isfinite(unit)
    assert unit > 0.0

    for alpha, k in product(alpha_a, k_a):
        wkb = Nc.HIPertITwoFluids.wkb_eval(qgrw, alpha, k)
        _test_wkb(wkb)
    for k, alpha in product(k_a, alpha_a):
        wkb = Nc.HIPertITwoFluids.wkb_eval(qgrw, alpha, k)
        _test_wkb(wkb)


def test_itwo_fluids_oem_eval(qgrw: Nc.HICosmoQGRW, alpha_a: np.ndarray) -> None:
    """Evaluate HICosmoQGRW implementation of NcHIPertITwoFluids at a given time."""
    k_a = np.geomspace(1.0e-3, 1.0e3, 5)

    for alpha, k in product(alpha_a, k_a):
        oem = Nc.HIPertITwoFluids.eom_eval(qgrw, alpha, k)
        _test_oem(oem)
    for k, alpha in product(k_a, alpha_a):
        oem = Nc.HIPertITwoFluids.eom_eval(qgrw, alpha, k)
        _test_oem(oem)
