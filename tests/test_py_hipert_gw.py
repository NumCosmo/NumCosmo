#!/usr/bin/env python
#
# test_py_hipert_gw.py
#
# Tue Apr 16 10:24:29 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_hipert_gw.py
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

"""Tests for the HIPertGW (tensor perturbation) interface and helpers.

This file exercises HIPertGW behavior (setters/getters, initial-condition
search, adiabatic solutions, evolution and power-spectrum evaluation) for
two example cosmologies used in the test-suite: `HICosmoVexp` and
`HICosmoQGRW`.
"""
import math
import pytest

from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="pgw_vexp")
def fixture_em() -> tuple[Nc.HIPertGW, Nc.HICosmoVexp]:
    """Provide a (HIPertGW, HICosmoVexp) pair configured for tests.

    The HIPertGW instance is configured with adiabatic initial-state
    parameters and time bounds derived from the `HICosmoVexp` model.
    """
    pgw = Nc.HIPertGW.new()
    vexp = Nc.HICosmoVexp()

    current_set = {
        "alphab": 7.4847e-3,  # Alpha (# of e-fold\s) at the bounce
        "sigmaphi": 100.0,  # Width of the Gaussian solution for the WdW equation
        "xb": 2.0e36,  # Inverse of the scale factor at the bounce (Initial condition)
        "dphi": -9.0e-4,  # Deviation of the Gaussian solution for the WdW equation
        "OmegaL": 1.0,  # H²(a when w=-1)/H²(a0). Basically gives the DE-dominated phase
        "Omegac": 1.0,  # Omega_d???
        "H0": 67.8,  # Hubble parameter today given by CMB observations
        "alphaem": 14.4,  # Amplitude of the EM gaussian coupling
        "betaem": 2.2,  # Width of the EM gaussian coupling
    }

    vexp.set_properties(**current_set)

    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)
    pgw.set_k(1.0)
    pgw.set_ti(vexp.tau_min())
    # Warning: this model has a singularity at the expanding phase, the final time
    # should be set prior to it.
    pgw.set_tf(vexp.tau_max())
    pgw.set_vacuum_max_time(-1.0e-1)
    pgw.set_vacuum_reltol(1.0e-8)

    return pgw, vexp


@pytest.fixture(name="pgw_qgrw")
def fixture_qgrw() -> tuple[Nc.HIPertGW, Nc.HICosmoQGRW]:
    """Provide a (HIPertGW, HICosmoQGRW) pair configured for tests.

    The HIPertGW instance is configured with adiabatic initial-state
    parameters and time bounds appropriate for the `HICosmoQGRW` model.
    """
    qgrw = Nc.HICosmoQGRW.new()
    qgrw["H0"] = 70.0
    qgrw["Omegar"] = 1.0e-5
    qgrw["Omegaw"] = 1.0 - 1.0e-5
    qgrw["w"] = 1.0e-5
    qgrw["xb"] = 1.0e30

    pgw = Nc.HIPertGW.new()
    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)
    pgw.set_k(1.0)
    pgw.set_ti(-qgrw.abs_alpha(1.0e-24))
    pgw.set_tf(1.0)

    pgw.set_vacuum_max_time(-1.0e1)
    pgw.set_vacuum_reltol(1.0e-8)

    return pgw, qgrw


def test_hipert_gw_vexp(pgw_vexp: tuple[Nc.HIPertGW, Nc.HICosmoVexp]) -> None:
    """Verify basic HIPertGW setters/getters and initial-condition modes.

    Exercises tolerances, thresholds, save flags and initial-condition
    type setters on a HIPertGW instance paired with a `HICosmoVexp`
    cosmology.
    """
    pgw, vexp = pgw_vexp

    assert_allclose(pgw.get_k(), 1.0)
    assert_allclose(pgw.get_ti(), vexp.tau_min())
    assert_allclose(pgw.get_tf(), vexp.tau_max())

    pgw.set_reltol(1.0e-6)
    pgw.set_abstol(1.0e-7)

    assert_allclose(pgw.get_reltol(), 1.0e-6)
    assert_allclose(pgw.get_abstol(), 1.0e-7)

    pgw.set_adiab_threshold(1.0e-3)
    pgw.set_prop_threshold(1.0e-3)

    assert_allclose(pgw.get_adiab_threshold(), 1.0e-3)
    assert_allclose(pgw.get_prop_threshold(), 1.0e-3)

    pgw.set_save_evol(True)
    assert pgw.get_save_evol()

    pgw.set_save_evol(False)
    assert not pgw.get_save_evol()

    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.AD_HOC)
    assert pgw.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.AD_HOC

    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC2)
    assert pgw.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC2

    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)
    assert pgw.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC4

    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.NONADIABATIC2)
    assert pgw.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.NONADIABATIC2


def test_initial_conditions_time_vexp(
    pgw_vexp: tuple[Nc.HIPertGW, Nc.HICosmoVexp],
) -> None:
    """Find and validate the adiabatic-time limit for the Vexp cosmology.

    Confirms the adiabatic time-finder returns a valid time inside the
    configured integration window and that the auxiliary `find_adiab_max`
    result matches evaluations of `eval_F1`.
    """
    pgw, vexp = pgw_vexp

    limit_found, t_adiab = pgw.find_adiab_time_limit(
        vexp, vexp.tau_min(), pgw.get_vacuum_max_time(), 1.0e-6
    )

    assert limit_found
    assert t_adiab >= pgw.get_ti()
    assert t_adiab <= pgw.get_tf()

    t_min, F1_min, t_lb, t_ub = pgw.find_adiab_max(
        vexp, vexp.tau_min(), pgw.get_vacuum_max_time(), 1.0e-1
    )

    assert_allclose(F1_min, pgw.eval_F1(vexp, t_min))
    assert math.fabs(F1_min - pgw.eval_F1(vexp, t_lb)) <= 1.0e-1
    assert math.fabs(F1_min - pgw.eval_F1(vexp, t_ub)) <= 1.0e-1


def test_initial_conditions_adiabatic_vexp(
    pgw_vexp: tuple[Nc.HIPertGW, Nc.HICosmoVexp],
) -> None:
    """Compute adiabatic initial state and verify the returned fields.

    For a range of target precisions, compute the adiabatic state and
    assert the returned complex field and its momentum are finite.
    """
    pgw, vexp = pgw_vexp

    state = Ncm.CSQ1DState.new()

    for prec in np.geomspace(1.0e-14, 1.0e-6, 10):
        limit_found, t_adiab = pgw.find_adiab_time_limit(
            vexp, vexp.tau_min(), pgw.get_vacuum_max_time(), prec
        )

        assert limit_found

        # Getting the adiabatic solution
        state, _alpha_reltol, _dgamma_reltol = pgw.compute_adiab(vexp, t_adiab, state)
        pgw.change_frame(vexp, state, Ncm.CSQ1DFrame.ORIG)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        # Compare with analytical solution
        assert np.isfinite(phi)
        assert np.isfinite(Pphi)


def test_evolution_vexp(pgw_vexp: tuple[Nc.HIPertGW, Nc.HICosmoVexp]) -> None:
    """Prepare and evolve HIPertGW for the Vexp cosmology and check outputs.

    Runs the evolution at configured times, verifies that the state
    invariants (J matrix) match expectations and that the instantaneous
    power spectrum computation returns consistent finite values.
    """
    pgw, vexp = pgw_vexp

    state = Ncm.CSQ1DState.new()

    pgw.set_tf(1.0)  # We do not want to evolve through the singularity
    pgw.prepare(vexp)

    t_a, _smaller_abst = pgw.get_time_array()
    unit = Nc.HIPertIGW.eval_unit(vexp)
    assert np.isfinite(unit)

    for t in t_a:
        state = pgw.eval_at(vexp, t, state)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        # Compare with analytical solution
        assert np.isfinite(abs(phi))
        assert np.isfinite(abs(Pphi))

        J11, J12, J22 = state.get_J()

        assert_allclose(J11, 2.0 * abs(phi) ** 2)
        assert_allclose(J22, 2.0 * abs(Pphi) ** 2)
        assert_allclose(J12, 2.0 * np.real(np.conj(phi) * Pphi))

        pk = pgw.eval_powspec_at(vexp, t)

        assert_allclose(pk, abs(unit * phi) ** 2 / (2.0 * math.pi**2))


def test_powspec_vexp(pgw_vexp: tuple[Nc.HIPertGW, Nc.HICosmoVexp]) -> None:
    """Evaluate the power spectrum spline for a range of k and assert finiteness.

    Ensures that `eval_powspec` returns an `Ncm.Spline` and that its
    evaluations produce finite numbers across a sampling grid.
    """
    pgw, vexp = pgw_vexp

    pgw.set_tf(1.0)  # We do not want to evolve through the singularity
    pgw.prepare(vexp)

    k_a = np.geomspace(1.0e-4, 1.0e4, 10)
    k_eval_a = np.geomspace(1.0e-4, 1.0e4, 50)

    pw_s = pgw.eval_powspec(vexp, 1.0, k_a)
    assert isinstance(pw_s, Ncm.Spline)
    pw_k_a = np.array([pw_s.eval(k) for k in k_eval_a])
    assert np.all(np.isfinite(pw_k_a))


def test_evolution_vexp_duplicate(pgw_vexp: tuple[Nc.HIPertGW, Nc.HICosmoVexp]) -> None:
    """Serialize/duplicate objects and run the evolution test on the copies.

    Verifies that serialization duplication preserves object types and
    that duplicated objects behave the same under evolution as the
    originals.
    """
    pgw, vexp = pgw_vexp

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    pgw_dup = ser.dup_obj(pgw)
    vexp_dup = ser.dup_obj(vexp)
    assert isinstance(pgw_dup, Nc.HIPertGW)
    assert isinstance(vexp_dup, Nc.HICosmoVexp)

    test_evolution_vexp((pgw_dup, vexp_dup))


def test_interface_eval_vexp(pgw_vexp: tuple[Nc.HIPertGW, Nc.HICosmoVexp]) -> None:
    """Directly exercise the HIPertIGW interface evaluation helpers.

    Calls static interface helpers (eval_unit, eval_F1, eval_m, eval_xi,
    eval_nu, eval_x) across a tau grid to ensure they return finite
    numerical values for the `HICosmoVexp` model.
    """
    _, vexp = pgw_vexp

    tau_a = np.linspace(vexp.tau_min() + 1.0, vexp.tau_max(), 1000)

    assert np.isfinite(Nc.HIPertIGW.eval_unit(vexp))
    for tau in tau_a:
        assert np.isfinite(Nc.HIPertIGW.eval_F1(vexp, tau, 1.0))
        assert np.isfinite(Nc.HIPertIGW.eval_m(vexp, tau, 1.0))
        assert np.isfinite(Nc.HIPertIGW.eval_xi(vexp, tau, 1.0))
        assert np.isfinite(Nc.HIPertIGW.eval_nu(vexp, tau, 1.0))
        assert np.isfinite(Nc.HIPertIGW.eval_x(vexp, tau))


def test_hipert_gw_qgrw(pgw_qgrw: tuple[Nc.HIPertGW, Nc.HICosmoQGRW]) -> None:
    """Verify HIPertGW setters/getters and state flags for QGRW cosmology.

    Mirrors the basic checks performed for the Vexp fixture but using
    `HICosmoQGRW` to ensure behavior is consistent across cosmology
    implementations.
    """
    pgw, qgrw = pgw_qgrw

    assert_allclose(pgw.get_k(), 1.0)
    assert_allclose(pgw.get_ti(), -qgrw.abs_alpha(1.0e-24))
    assert_allclose(pgw.get_tf(), 1.0)

    pgw.set_reltol(1.0e-6)
    pgw.set_abstol(1.0e-7)

    assert_allclose(pgw.get_reltol(), 1.0e-6)
    assert_allclose(pgw.get_abstol(), 1.0e-7)

    pgw.set_adiab_threshold(1.0e-3)
    pgw.set_prop_threshold(1.0e-3)

    assert_allclose(pgw.get_adiab_threshold(), 1.0e-3)
    assert_allclose(pgw.get_prop_threshold(), 1.0e-3)

    pgw.set_save_evol(True)
    assert pgw.get_save_evol()

    pgw.set_save_evol(False)
    assert not pgw.get_save_evol()

    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.AD_HOC)
    assert pgw.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.AD_HOC

    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC2)
    assert pgw.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC2

    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)
    assert pgw.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC4

    pgw.set_initial_condition_type(Ncm.CSQ1DInitialStateType.NONADIABATIC2)
    assert pgw.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.NONADIABATIC2


def test_initial_conditions_time_qgrw(
    pgw_qgrw: tuple[Nc.HIPertGW, Nc.HICosmoQGRW],
) -> None:
    """Find and validate the adiabatic-time limit for the QGRW cosmology.

    Checks that the time finder returns a valid time and that `find_adiab_max`
    results are consistent with `eval_F1` evaluations.
    """
    pgw, qgrw = pgw_qgrw

    limit_found, t_adiab = pgw.find_adiab_time_limit(
        qgrw, -qgrw.abs_alpha(1.0e-24), pgw.get_vacuum_max_time(), 1.0e-6
    )

    assert limit_found
    assert t_adiab >= pgw.get_ti()
    assert t_adiab <= pgw.get_tf()

    t_min, F1_min, t_lb, t_ub = pgw.find_adiab_max(
        qgrw, -qgrw.abs_alpha(1.0e-24), pgw.get_vacuum_max_time(), 1.0e-1
    )

    assert_allclose(F1_min, pgw.eval_F1(qgrw, t_min))
    assert math.fabs(F1_min - pgw.eval_F1(qgrw, t_lb)) <= 1.0e-1
    assert math.fabs(F1_min - pgw.eval_F1(qgrw, t_ub)) <= 1.0e-1


def test_initial_conditions_adiabatic_qgrw(
    pgw_qgrw: tuple[Nc.HIPertGW, Nc.HICosmoQGRW],
) -> None:
    """Compute adiabatic initial states for QGRW and assert finiteness.

    Iterates over a range of requested precisions, computes the adiabatic
    solution, and asserts the complex field and its momentum are finite.
    """
    pgw, qgrw = pgw_qgrw

    state = Ncm.CSQ1DState.new()

    for prec in np.geomspace(1.0e-14, 1.0e-6, 10):
        limit_found, t_adiab = pgw.find_adiab_time_limit(
            qgrw, -qgrw.abs_alpha(1.0e-24), pgw.get_vacuum_max_time(), prec
        )

        assert limit_found

        # Getting the adiabatic solution
        state, _alpha_reltol, _dgamma_reltol = pgw.compute_adiab(qgrw, t_adiab, state)
        pgw.change_frame(qgrw, state, Ncm.CSQ1DFrame.ORIG)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        # Compare with analytical solution
        assert np.isfinite(phi)
        assert np.isfinite(Pphi)


def test_evolution_qgrw(pgw_qgrw: tuple[Nc.HIPertGW, Nc.HICosmoQGRW]) -> None:
    """Prepare and evolve HIPertGW for the QGRW cosmology and validate outputs.

    Evolves the HIPertGW state, checks the J invariants and compares the
    instantaneous power spectrum with the analytical combination of the
    unit and mode field.
    """
    pgw, qgrw = pgw_qgrw

    state = Ncm.CSQ1DState.new()

    pgw.set_tf(1.0)
    pgw.prepare(qgrw)

    t_a, _smaller_abst = pgw.get_time_array()
    unit = Nc.HIPertIGW.eval_unit(qgrw)
    assert np.isfinite(unit)

    for t in t_a:
        state = pgw.eval_at(qgrw, t, state)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        # Compare with analytical solution
        assert np.isfinite(abs(phi))
        assert np.isfinite(abs(Pphi))

        J11, J12, J22 = state.get_J()

        assert_allclose(J11, 2.0 * abs(phi) ** 2)
        assert_allclose(J22, 2.0 * abs(Pphi) ** 2)
        assert_allclose(J12, 2.0 * np.real(np.conj(phi) * Pphi))

        pk = pgw.eval_powspec_at(qgrw, t)
        assert_allclose(pk, abs(unit * phi) ** 2 / (2.0 * math.pi**2))


def test_powspec_qgrw(pgw_qgrw: tuple[Nc.HIPertGW, Nc.HICosmoQGRW]) -> None:
    """Evaluate the power spectrum spline for QGRW and assert finiteness.

    Ensures `eval_powspec` returns an `Ncm.Spline` and that its evaluations
    over a sampling grid are finite numbers.
    """
    pgw, qgrw = pgw_qgrw

    pgw.set_tf(1.0)
    pgw.prepare(qgrw)

    k_a = np.geomspace(1.0e-4, 1.0e4, 10)
    k_eval_a = np.geomspace(1.0e-4, 1.0e4, 50)

    pw_s = pgw.eval_powspec(qgrw, 1.0, k_a)
    assert isinstance(pw_s, Ncm.Spline)
    pw_k_a = np.array([pw_s.eval(k) for k in k_eval_a])
    assert np.all(np.isfinite(pw_k_a))


def test_evolution_qgrw_duplicate(pgw_qgrw: tuple[Nc.HIPertGW, Nc.HICosmoQGRW]) -> None:
    """Serialize/duplicate QGRW objects and run the evolution on the copies.

    Verifies duplicated objects preserve types and that evolution behaves
    correctly on serialized copies.
    """
    pgw, qgrw = pgw_qgrw

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    pgw_dup = ser.dup_obj(pgw)
    qgrw_dup = ser.dup_obj(qgrw)
    assert isinstance(pgw_dup, Nc.HIPertGW)
    assert isinstance(qgrw_dup, Nc.HICosmoQGRW)

    test_evolution_qgrw((pgw_dup, qgrw_dup))


def test_interface_eval_qgrw(pgw_qgrw: tuple[Nc.HIPertGW, Nc.HICosmoQGRW]) -> None:
    """Exercise HIPertIGW interface evaluation helpers for QGRW cosmology.

    Calls the static interface helpers (eval_unit, eval_F1, eval_m, eval_xi,
    eval_nu, eval_x) across a tau grid to assert they return finite values.
    """
    _, qgrw = pgw_qgrw

    tau_a = np.linspace(-qgrw.abs_alpha(1.0e-24), 1.0, 1000)

    assert np.isfinite(Nc.HIPertIGW.eval_unit(qgrw))
    for tau in tau_a:
        assert np.isfinite(Nc.HIPertIGW.eval_F1(qgrw, tau, 1.0))
        assert np.isfinite(Nc.HIPertIGW.eval_m(qgrw, tau, 1.0))
        assert np.isfinite(Nc.HIPertIGW.eval_xi(qgrw, tau, 1.0))
        assert np.isfinite(Nc.HIPertIGW.eval_nu(qgrw, tau, 1.0))
        assert np.isfinite(Nc.HIPertIGW.eval_x(qgrw, tau))
