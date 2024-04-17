#!/usr/bin/env python
#
# test_py_hipert_adiab.py
#
# Tue Apr 16 10:24:29 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_hipert_adiab.py
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

"""Tests on NcHIPertAdiab class."""
import math
import pytest

from numpy.testing import assert_allclose
from scipy.special import hankel1e  # pylint: disable=no-name-in-module
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="adiab")
def fixture_adiab():
    """Fixture for NcHIPertAdiab."""
    adiab = Nc.HIPertAdiab.new()

    adiab.set_k(1.0)

    adiab.set_ti(-300.0)
    adiab.set_tf(-10.0)

    return adiab


@pytest.fixture(name="qgw")
def fixture_qgw() -> Nc.HICosmoQGW:
    """Fixture for NcHIPertAdiab model Nc.HICosmoQGW."""
    return Nc.HICosmoQGW()


def test_hipert_adiab_qgw(adiab):
    """Test basic functionality of NcHIPertAdiab."""
    assert_allclose(adiab.get_k(), 1.0)
    assert_allclose(adiab.get_ti(), -300.0)
    # We need to set the final time to a value away from the bounce since our
    # theoretical solution is not valid near the bounce.
    assert_allclose(adiab.get_tf(), -10.0)

    adiab.set_reltol(1.0e-6)
    adiab.set_abstol(1.0e-7)

    assert_allclose(adiab.get_reltol(), 1.0e-6)
    assert_allclose(adiab.get_abstol(), 1.0e-7)

    adiab.set_adiab_threshold(1.0e-3)
    adiab.set_prop_threshold(1.0e-3)

    assert_allclose(adiab.get_adiab_threshold(), 1.0e-3)
    assert_allclose(adiab.get_prop_threshold(), 1.0e-3)

    adiab.set_save_evol(True)
    assert adiab.get_save_evol()

    adiab.set_save_evol(False)
    assert not adiab.get_save_evol()

    adiab.set_initial_condition_type(Ncm.CSQ1DInitialStateType.AD_HOC)
    assert adiab.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.AD_HOC

    adiab.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC2)
    assert adiab.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC2

    adiab.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)
    assert adiab.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC4

    adiab.set_initial_condition_type(Ncm.CSQ1DInitialStateType.NONADIABATIC2)
    assert adiab.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.NONADIABATIC2


def test_initial_conditions_time_qgw(adiab, qgw):
    """Test initial conditions of NcHIPertAdiab."""
    limit_found, t_adiab = adiab.find_adiab_time_limit(qgw, -300.0, -1.0e-3, 1.0e-6)

    assert limit_found
    assert t_adiab >= adiab.get_ti()
    assert t_adiab <= adiab.get_tf()

    t_min, F1_min, t_lb, t_ub = adiab.find_adiab_max(qgw, -300.0, -1.0e1, 1.0e-1)

    assert_allclose(F1_min, adiab.eval_F1(qgw, t_min))
    assert math.fabs(F1_min - adiab.eval_F1(qgw, t_lb)) <= 1.0e-1
    assert math.fabs(F1_min - adiab.eval_F1(qgw, t_ub)) <= 1.0e-1


def _compute_analytical_solution_qgw(adiab, qgw, t_adiab):
    """Compute analytical solution for NcHIPertAdiab using NcHICosmoQGW."""
    w = qgw.props.w
    cs2 = w
    Omegaw = qgw.props.Omegaw
    sqrtOmegaw = np.sqrt(Omegaw)
    norma = np.sqrt(np.pi * cs2 / (6.0 * (1.0 + 3.0 * w) * (1.0 + w) * sqrtOmegaw))
    w = qgw.props.w
    alpha = 3.0 * (1.0 - w) / (2.0 * (1.0 + 3.0 * w))
    x = qgw.eval_x(t_adiab)
    E = np.abs(qgw.eval_hubble(t_adiab))
    eta = 2.0 * x / (E * (1.0 + 3.0 * w))
    cs = np.sqrt(w)
    csketa = cs * adiab.get_k() * eta
    hfnormm = norma * (2.0 / ((1.0 + 3.0 * w) * eta * sqrtOmegaw)) ** (alpha)

    # Analytical solution for phi and Pphi
    theo_phi = hfnormm * hankel1e(alpha, csketa)
    theo_Pphi = (csketa / hfnormm) * 0.25 * math.pi * hankel1e(1.0 + alpha, csketa)

    return theo_phi, theo_Pphi


def test_initial_conditions_adiabatic_qgw(adiab, qgw):
    """Test initial conditions of NcHIPertAdiab."""
    state = Ncm.CSQ1DState.new()

    for prec in np.geomspace(1.0e-14, 1.0e-6, 100):
        limit_found, t_adiab = adiab.find_adiab_time_limit(qgw, -300.0, -1.0e1, prec)

        assert limit_found

        # Getting the adiabatic solution
        state, _alpha_reltol, _dgamma_reltol = adiab.compute_adiab(qgw, t_adiab, state)
        adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ORIG)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        theo_phi, theo_Pphi = _compute_analytical_solution_qgw(adiab, qgw, t_adiab)

        # Compare with analytical solution
        assert_allclose(abs(phi), abs(theo_phi), rtol=1.0e-8)
        assert_allclose(abs(Pphi), abs(theo_Pphi), rtol=1.0e-8)


def test_evolution_qgw(adiab, qgw):
    """Test initial conditions of NcHIPertAdiab."""
    state = Ncm.CSQ1DState.new()

    adiab.set_save_evol(True)
    adiab.set_reltol(1.0e-10)
    adiab.set_abstol(0.0)
    adiab.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)
    adiab.set_vacuum_max_time(-1.0e1)
    adiab.set_vacuum_reltol(1.0e-8)
    adiab.prepare(qgw)

    t_a, _smaller_abst = adiab.get_time_array()

    for t in t_a:
        state = adiab.eval_at(qgw, t, state)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        # Analytical solution for phi and Pphi
        theo_phi, theo_Pphi = _compute_analytical_solution_qgw(adiab, qgw, t)

        # Compare with analytical solution
        assert_allclose(abs(phi), abs(theo_phi), rtol=1.0e-7)
        assert_allclose(abs(Pphi), abs(theo_Pphi), rtol=1.0e-7)

        J11, J12, J22 = state.get_J()

        assert_allclose(J11, 2.0 * abs(phi) ** 2, atol=1.0e-7)
        assert_allclose(J22, 2.0 * abs(Pphi) ** 2, atol=1.0e-7)
        assert_allclose(J12, (2.0 * phi * Pphi.conjugate()).real, atol=1.0e-7)


def test_evolution_adiabatic2_qgw(adiab, qgw):
    """Test initial conditions of NcHIPertAdiab."""
    adiab.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC2)
    adiab.set_vacuum_max_time(-1.0e0)
    adiab.set_vacuum_reltol(1.0e-3)
    state = Ncm.CSQ1DState.new()

    adiab.set_save_evol(True)
    adiab.set_reltol(1.0e-10)
    adiab.set_abstol(0.0)
    adiab.prepare(qgw)

    t_a, _smaller_abst = adiab.get_time_array()

    for t in t_a:
        state = adiab.eval_at(qgw, t, state)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        # Analytical solution for phi and Pphi
        theo_phi, theo_Pphi = _compute_analytical_solution_qgw(adiab, qgw, t)

        # Compare with analytical solution
        assert_allclose(abs(phi), abs(theo_phi), rtol=1.0e-7)
        assert_allclose(abs(Pphi), abs(theo_Pphi), rtol=1.0e-7)

        J11, J12, J22 = state.get_J()

        assert_allclose(J11, 2.0 * abs(phi) ** 2, atol=1.0e-7)
        assert_allclose(J22, 2.0 * abs(Pphi) ** 2, atol=1.0e-7)
        assert_allclose(J12, (2.0 * phi * Pphi.conjugate()).real, atol=1.0e-7)


@pytest.mark.parametrize(
    "frame",
    [Ncm.CSQ1DFrame.ORIG, Ncm.CSQ1DFrame.ADIAB1, Ncm.CSQ1DFrame.ADIAB2],
    ids=["orig", "adiab1", "adiab2"],
)
def test_evolution_frame_qgw(adiab, qgw, frame):
    """Test initial conditions of NcHIPertAdiab."""
    state = Ncm.CSQ1DState.new()
    ti = adiab.get_ti()
    adiab.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)

    limit_found, t_adiab_end = adiab.find_adiab_time_limit(qgw, ti, -1.0e0, 1.0e-1)
    assert limit_found
    adiab.set_tf(t_adiab_end)

    adiab.set_save_evol(True)
    adiab.set_reltol(1.0e-10)
    adiab.set_abstol(0.0)
    adiab.prepare(qgw)

    t_a, _smaller_abst = adiab.get_time_array()

    for t in t_a:
        state = adiab.eval_at_frame(qgw, frame, t, state)
        assert state.get_frame() == frame


def test_change_frame_orig_adiab1_qgw(adiab, qgw):
    """Test change_frame method of NcHIPertAdiab."""
    t = -1.4e2
    state = Ncm.CSQ1DState()
    state.set_ag(Ncm.CSQ1DFrame.ORIG, t, 0.1, 0.3)

    adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ORIG)

    assert state.get_frame() == Ncm.CSQ1DFrame.ORIG
    assert state.get_ag() == (0.1, 0.3)

    adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ADIAB1)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB1
    assert_allclose(state.get_ag(), (0.1, 0.3 - adiab.eval_xi(qgw, t)))

    adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ORIG)
    assert state.get_frame() == Ncm.CSQ1DFrame.ORIG
    assert_allclose(state.get_ag(), (0.1, 0.3))


def test_change_frame_orig_adiab2_qgw(adiab, qgw):
    """Test change_frame method of NcHIPertAdiab."""
    t = -1.4e2

    alpha0 = 0.1
    gamma0 = 0.3 + adiab.eval_xi(qgw, t)

    state = Ncm.CSQ1DState()
    state.set_ag(Ncm.CSQ1DFrame.ORIG, t, alpha0, gamma0)

    adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ORIG)

    assert state.get_frame() == Ncm.CSQ1DFrame.ORIG
    assert state.get_ag() == (alpha0, gamma0)

    adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ADIAB2)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB2
    assert not np.allclose(state.get_ag(), (alpha0, gamma0))

    adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ORIG)
    assert state.get_frame() == Ncm.CSQ1DFrame.ORIG
    assert_allclose(state.get_ag(), (alpha0, gamma0))


def test_change_frame_adiab1_adiab2_qgw(adiab, qgw):
    """Test change_frame method of NcHIPertAdiab."""
    t = -1.4e2
    state = Ncm.CSQ1DState()
    state.set_ag(Ncm.CSQ1DFrame.ADIAB1, t, 0.1, 0.3)

    adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ADIAB1)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB1
    assert state.get_ag() == (0.1, 0.3)

    adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ADIAB2)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB2
    assert not np.allclose(state.get_ag(), (0.1, 0.3))

    adiab.change_frame(qgw, state, Ncm.CSQ1DFrame.ADIAB1)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB1
    assert_allclose(state.get_ag(), (0.1, 0.3))
