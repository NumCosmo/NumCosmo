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

"""Tests on NcHIPertEM class."""
import math
import pytest

from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="pgw_vexp")
def fixture_em():
    """Fixture for NcHIPertEM."""
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
        "Bem": -1.0 / 4.0 + 1.0e-12,  # Amplitude of the EM gaussian coupling
        "betaem": 3.6e-1,  # Width of the EM gaussian coupling
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


def test_hipert_gw_vexp(pgw_vexp):
    """Test basic functionality of NcHIPertAdiab."""
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


def test_initial_conditions_time_vexp(pgw_vexp):
    """Test initial conditions of NcHIPertAdiab."""
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


def test_initial_conditions_adiabatic_vexp(pgw_vexp):
    """Test initial conditions of NcHIPertAdiab."""
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


def test_evolution_vexp(pgw_vexp):
    """Test initial conditions of NcHIPertAdiab."""
    pgw, vexp = pgw_vexp

    state = Ncm.CSQ1DState.new()

    pgw.set_tf(1.0)  # We do not want to evolve through the singularity
    pgw.prepare(vexp)

    t_a, _smaller_abst = pgw.get_time_array()

    for t in t_a:
        state = pgw.eval_at(vexp, t, state)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        # Compare with analytical solution
        assert np.isfinite(abs(phi))
        assert np.isfinite(abs(Pphi))

        J11, J12, J22 = state.get_J()

        assert np.isfinite(J11)
        assert np.isfinite(J22)
        assert np.isfinite(J12)


def test_evolution_vexp_duplicate(pgw_vexp):
    """Test initial conditions of NcHIPertAdiab."""
    pgw, vexp = pgw_vexp

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    pgw_dup = ser.dup_obj(pgw)
    vexp_dup = ser.dup_obj(vexp)

    test_evolution_vexp((pgw_dup, vexp_dup))
