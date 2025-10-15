#!/usr/bin/env python
#
# test_py_hipert_two_fluids.py
#
# Mon May 20 00:24:30 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_hipert_two_fluids.py
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

"""Tests on NcHIPertTwoFluids perturbations module."""

import itertools as it
import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="two_fluids")
def fixture_two_fluids() -> Nc.HIPertTwoFluids:
    """Fixture for NcHIPertTwoFluids."""
    two_fluids = Nc.HIPertTwoFluids.new()
    two_fluids.set_initial_time(-160.0)
    two_fluids.set_final_time(1.0)
    two_fluids.set_wkb_reltol(1.0e-3)

    return two_fluids


@pytest.fixture(name="cosmo_qgrw")
def fixture_cosmo_qgrw() -> Nc.HICosmo:
    """Fixture for NcCosmoQGrw."""
    cosmo = Nc.HICosmoQGRW()
    cosmo["w"] = 1.0e-5
    cosmo["Omegar"] = 1.0 * (1.0e-5)
    cosmo["Omegaw"] = 1.0 * (1.0 - 1.0e-5)
    cosmo["xb"] = 1.0e30

    return cosmo


def test_init(two_fluids: Nc.HIPertTwoFluids) -> None:
    """Test NcHIPertTwoFluids initialization."""
    assert two_fluids is not None
    assert isinstance(two_fluids, Nc.HIPertTwoFluids)
    assert isinstance(two_fluids, Nc.HIPert)

    two_fluids.set_prepared(True)
    assert two_fluids.prepared() is True

    two_fluids.set_prepared(False)
    assert two_fluids.prepared() is False

    two_fluids.set_stiff_solver(True)
    two_fluids.set_stiff_solver(False)


def test_compute_full_spectrum(
    two_fluids: Nc.HIPertTwoFluids, cosmo_qgrw: Nc.HICosmo
) -> None:
    """Test NcHIPertTwoFluids compute_full_spectrum."""
    two_fluids.props.reltol = 1.0e-9

    def spec_params(Omega_rs=1.0e-5, w=1.0e-3, E0=1.0):
        cosmo_qgrw["w"] = w
        cosmo_qgrw["Omegar"] = E0 * Omega_rs
        cosmo_qgrw["Omegaw"] = E0 * (1.0 - Omega_rs)

        spec1 = two_fluids.compute_zeta_spectrum(
            cosmo_qgrw, 1, -cosmo_qgrw.abs_alpha(1.0e-14), -1.0, 1.0e-3, 1.0e8, 10
        )
        spec2 = two_fluids.compute_zeta_spectrum(
            cosmo_qgrw, 2, -cosmo_qgrw.abs_alpha(1.0e-14), -1.0, 1.0e-3, 1.0e8, 10
        )

        return spec1, spec2

    specs1 = []
    specs2 = []
    w_a = np.geomspace(1.0e-5, 1.0e-1, 6)
    for w in w_a:
        spec1, spec2 = spec_params(w=w)
        specs1.append(spec1)
        specs2.append(spec2)

    lnk_v = specs1[0].peek_xv()
    lnw_v = Ncm.Vector.new_array(np.log(w_a).tolist())
    zm = Ncm.Matrix.new(lnw_v.len(), lnk_v.len())

    for i, (spec1, spec2, w) in enumerate(zip(specs1, specs2, w_a)):
        lnPk_a1 = np.array(spec1.peek_yv().dup_array())
        lnPk_a2 = np.array(spec2.peek_yv().dup_array())

        lnPk0 = np.log(np.exp(lnPk_a1[0]) + np.exp(lnPk_a2[0]))
        lnPk = np.log(np.exp(lnPk_a1) + np.exp(lnPk_a2)) - lnPk0
        lnPk_v = Ncm.Vector.new_array(lnPk)

        zm.set_row(i, lnPk_v)

    s = Ncm.Spline2dBicubic.notaknot_new()
    s.set(lnk_v, lnw_v, zm, True)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    default_calib_file = Ncm.cfg_get_data_filename("hiprim_2f_spline.bin", True)
    s_calib = ser.from_binfile(default_calib_file)
    assert isinstance(s_calib, Ncm.Spline2dBicubic)

    for w, lnk in it.product(w_a, lnk_v.dup_array()):
        Pk0 = s.eval(lnk, np.log(w))
        Pk1 = s_calib.eval(lnk, np.log(w))
        assert_allclose(Pk0, Pk1, rtol=1.0e-3)


def test_evolve_array(two_fluids: Nc.HIPertTwoFluids, cosmo_qgrw: Nc.HICosmo):
    """Test NcHIPertTwoFluids evolve_array."""
    init_cond = Ncm.Vector.new_array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    two_fluids.props.reltol = 1.0e-9

    alpha_i = two_fluids.get_wkb_limit(cosmo_qgrw, 1, -90.0, 1.0e-8)
    assert alpha_i < 0.0
    two_fluids.get_init_cond_zetaS(cosmo_qgrw, alpha_i, 1, np.pi * 0.25, init_cond)
    two_fluids.set_init_cond(cosmo_qgrw, alpha_i, 1, False, init_cond)

    m = two_fluids.evolve_array(
        cosmo=cosmo_qgrw, alphaf=-1.0e-1, step_abstol=0.0, step_reltol=1.0e-5
    )

    m_a = np.array(m.dup_array())

    assert len(m_a) > 0
    assert all(np.isfinite(m_a))


def test_evol_mode_and_state_interp(
    two_fluids: Nc.HIPertTwoFluids, cosmo_qgrw: Nc.HICosmo
):
    """Test `evol_mode` and the returned NcHIPertTwoFluidsStateInterp."""
    two_fluids.props.reltol = 1.0e-9

    s_interp = two_fluids.evol_mode(cosmo_qgrw)
    assert s_interp is not None
    assert isinstance(s_interp, Nc.HIPertTwoFluidsStateInterp)

    ti = two_fluids.get_initial_time()
    tf = two_fluids.get_final_time()
    k = two_fluids.get_mode_k()

    assert ti < tf
    alpha_a = np.linspace(ti, tf, 100)

    s_interp2 = s_interp.dup()
    assert s_interp2 is not None
    assert isinstance(s_interp2, Nc.HIPertTwoFluidsStateInterp)
    assert s_interp2 is not s_interp

    for alpha in alpha_a:
        state = s_interp.eval(cosmo_qgrw, alpha)
        state2 = s_interp2.eval(cosmo_qgrw, alpha)
        assert isinstance(state, Nc.HIPertITwoFluidsState)
        assert state.alpha == alpha
        assert state.k == k
        assert state.gw1 > 0.0
        assert state.gw2 > 0.0
        assert state.gw1 != state.gw2
        assert state.Fnu >= 0.0
        assert state.norma > 0.0
        assert_allclose(state.gw1, state2.gw1, rtol=1.0e-11, atol=0.0)
        assert_allclose(state.gw2, state2.gw2, rtol=1.0e-11, atol=0.0)
        assert_allclose(state.Fnu, state2.Fnu, rtol=1.0e-11, atol=0.0)
        assert_allclose(state.norma, state2.norma, rtol=1.0e-11, atol=0.0)

        for mode, obs in it.product(
            [Nc.HIPertITwoFluidsObsMode.ONE, Nc.HIPertITwoFluidsObsMode.TWO],
            list(Nc.HIPertITwoFluidsObs),
        ):
            val = state.eval_mode(mode, obs)
            assert isinstance(val, Ncm.Complex)
            assert np.isfinite(val.Re())
            assert np.isfinite(val.Im())

        for obs1, obs2 in it.combinations_with_replacement(
            list(Nc.HIPertITwoFluidsObs), 2
        ):
            val1 = state.eval_obs(Nc.HIPertITwoFluidsObsMode.ONE, obs1, obs2)
            val2 = state.eval_obs(Nc.HIPertITwoFluidsObsMode.TWO, obs1, obs2)
            val = state.eval_obs(Nc.HIPertITwoFluidsObsMode.BOTH, obs1, obs2)
            assert isinstance(val1, float)
            assert isinstance(val2, float)
            assert isinstance(val, float)
            assert np.isfinite(val1)
            assert np.isfinite(val2)
            assert np.isfinite(val)
            assert_allclose(
                val1 + val2, val, rtol=1.0e-11, atol=1.0e-11 * (abs(val1) + abs(val2))
            )


def test_wkb_eval(cosmo_qgrw: Nc.HICosmo):
    """Test NcHIPertITwoFluidsWKB evaluation."""

    k_a = np.geomspace(1.0e-3, 1.0e3, 5)
    alpha_a = np.linspace(-100.0, 1.0, 5)

    for k, alpha in it.product(k_a, alpha_a):
        wkb = Nc.HIPertITwoFluids.wkb_eval(cosmo_qgrw, alpha, k)
        assert isinstance(wkb, Nc.HIPertITwoFluidsWKB)
        wkb2 = wkb.dup()
        assert isinstance(wkb2, Nc.HIPertITwoFluidsWKB)
        assert wkb2 is not wkb
        assert_allclose(wkb.mode1_zeta_scale, wkb2.mode1_zeta_scale, rtol=0.0, atol=0.0)
        assert_allclose(wkb.mode2_zeta_scale, wkb2.mode2_zeta_scale, rtol=0.0, atol=0.0)
        assert_allclose(wkb.mode1_Q_scale, wkb2.mode1_Q_scale, rtol=0.0, atol=0.0)
        assert_allclose(wkb.mode2_Q_scale, wkb2.mode2_Q_scale, rtol=0.0, atol=0.0)
        assert_allclose(
            wkb.mode1_Pzeta_scale, wkb2.mode1_Pzeta_scale, rtol=0.0, atol=0.0
        )
        assert_allclose(
            wkb.mode2_Pzeta_scale, wkb2.mode2_Pzeta_scale, rtol=0.0, atol=0.0
        )
        assert_allclose(wkb.mode1_PQ_scale, wkb2.mode1_PQ_scale, rtol=0.0, atol=0.0)
        assert_allclose(wkb.mode2_PQ_scale, wkb2.mode2_PQ_scale, rtol=0.0, atol=0.0)


def test_eom_eval(cosmo_qgrw: Nc.HICosmo):
    """Test NcHIPertITwoFluidsOEM evaluation."""
    k_a = np.geomspace(1.0e-3, 1.0e3, 5)
    alpha_a = np.linspace(-100.0, 1.0, 5)

    for k, alpha in it.product(k_a, alpha_a):
        eom = Nc.HIPertITwoFluids.eom_eval(cosmo_qgrw, alpha, k)
        assert isinstance(eom, Nc.HIPertITwoFluidsEOM)
        eom2 = eom.dup()
        assert isinstance(eom2, Nc.HIPertITwoFluidsEOM)
        assert eom2 is not eom
        assert_allclose(eom.gw1, eom2.gw1, rtol=0.0, atol=0.0)
        assert_allclose(eom.gw2, eom2.gw2, rtol=0.0, atol=0.0)
        assert_allclose(eom.Fnu, eom2.Fnu, rtol=0.0, atol=0.0)
        assert_allclose(eom.cs2, eom2.cs2, rtol=0.0, atol=0.0)
        assert_allclose(eom.cm2, eom2.cm2, rtol=0.0, atol=0.0)


def test_state(cosmo_qgrw: Nc.HICosmo):
    """Test NcHIPertITwoFluidsWKB evaluation."""

    k_a = np.geomspace(1.0e-3, 1.0e3, 5)
    alpha_a = np.linspace(-100.0, 1.0, 5)

    for k, alpha in it.product(k_a, alpha_a):
        wkb = Nc.HIPertITwoFluids.wkb_eval(cosmo_qgrw, alpha, k)
        state = wkb.peek_state()

        assert isinstance(state, Nc.HIPertITwoFluidsState)

        state2 = state.dup()
        assert isinstance(state2, Nc.HIPertITwoFluidsState)
        assert state2 is not state
        assert_allclose(state.gw1, state2.gw1, rtol=0.0, atol=0.0)
        assert_allclose(state.gw2, state2.gw2, rtol=0.0, atol=0.0)
        assert_allclose(state.Fnu, state2.Fnu, rtol=0.0, atol=0.0)
        assert_allclose(state.norma, state2.norma, rtol=0.0, atol=0.0)
