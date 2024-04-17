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

"""Tests on NcHIPertEM class."""
import math
import pytest

from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(
    name="pem_vexp",
    params=[
        Nc.HICosmoVexpEMCoupling.NONE,
        Nc.HICosmoVexpEMCoupling.CAUCHY,
        Nc.HICosmoVexpEMCoupling.GAUSS,
    ],
    ids=["none", "cauchy", "gauss"],
)
def fixture_em(request):
    """Fixture for NcHIPertEM."""
    pem = Nc.HIPertEM.new()
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
    vexp.set_em_coupling(request.param)

    # No coupling has F1 exactly zero, so computing F2 numerically can be problematic
    # using ADIABATIC2 avoids this problem.
    # Gaussians have a very small F1 (suppressed by the Gaussian term), so it is also
    # problematic to compute F2 numerically.
    if request.param == Nc.HICosmoVexpEMCoupling.NONE:
        pem.set_abstol(1.0e-200)
        pem.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC2)
    elif request.param == Nc.HICosmoVexpEMCoupling.GAUSS:
        pem.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC2)
    else:
        pem.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)

    pem.set_k(1.0)

    pem.set_ti(vexp.tau_min())
    # Warning: this model has a singularity at the expanding phase, the final time
    # should be set prior to it.
    pem.set_tf(vexp.tau_max())
    pem.set_vacuum_max_time(-1.0e-1)
    pem.set_vacuum_reltol(1.0e-8)

    return pem, vexp


def test_hipert_adiab_vexp(pem_vexp):
    """Test basic functionality of NcHIPertAdiab."""
    pem, vexp = pem_vexp

    assert_allclose(pem.get_k(), 1.0)
    assert_allclose(pem.get_ti(), vexp.tau_min())
    assert_allclose(pem.get_tf(), vexp.tau_max())

    pem.set_reltol(1.0e-6)
    pem.set_abstol(1.0e-7)

    assert_allclose(pem.get_reltol(), 1.0e-6)
    assert_allclose(pem.get_abstol(), 1.0e-7)

    pem.set_adiab_threshold(1.0e-3)
    pem.set_prop_threshold(1.0e-3)

    assert_allclose(pem.get_adiab_threshold(), 1.0e-3)
    assert_allclose(pem.get_prop_threshold(), 1.0e-3)

    pem.set_save_evol(True)
    assert pem.get_save_evol()

    pem.set_save_evol(False)
    assert not pem.get_save_evol()

    pem.set_initial_condition_type(Ncm.CSQ1DInitialStateType.AD_HOC)
    assert pem.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.AD_HOC

    pem.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC2)
    assert pem.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC2

    pem.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)
    assert pem.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC4

    pem.set_initial_condition_type(Ncm.CSQ1DInitialStateType.NONADIABATIC2)
    assert pem.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.NONADIABATIC2


def test_initial_conditions_time_vexp(pem_vexp):
    """Test initial conditions of NcHIPertAdiab."""
    pem, vexp = pem_vexp

    limit_found, t_adiab = pem.find_adiab_time_limit(
        vexp, vexp.tau_min(), pem.get_vacuum_max_time(), 1.0e-6
    )

    assert limit_found
    assert t_adiab >= pem.get_ti()
    assert t_adiab <= pem.get_tf()

    t_min, F1_min, t_lb, t_ub = pem.find_adiab_max(
        vexp, vexp.tau_min(), pem.get_vacuum_max_time(), 1.0e-1
    )

    assert_allclose(F1_min, pem.eval_F1(vexp, t_min))
    assert math.fabs(F1_min - pem.eval_F1(vexp, t_lb)) <= 1.0e-1
    assert math.fabs(F1_min - pem.eval_F1(vexp, t_ub)) <= 1.0e-1


def test_initial_conditions_adiabatic_vexp(pem_vexp):
    """Test initial conditions of NcHIPertAdiab."""
    pem, vexp = pem_vexp

    state = Ncm.CSQ1DState.new()

    for prec in np.geomspace(1.0e-14, 1.0e-6, 10):
        limit_found, t_adiab = pem.find_adiab_time_limit(
            vexp, vexp.tau_min(), pem.get_vacuum_max_time(), prec
        )

        assert limit_found

        # Getting the adiabatic solution
        state, _alpha_reltol, _dgamma_reltol = pem.compute_adiab(vexp, t_adiab, state)
        pem.change_frame(vexp, state, Ncm.CSQ1DFrame.ORIG)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        # Compare with analytical solution
        assert np.isfinite(phi)
        assert np.isfinite(Pphi)


def test_evolution_vexp(pem_vexp):
    """Test initial conditions of NcHIPertAdiab."""
    pem, vexp = pem_vexp

    state = Ncm.CSQ1DState.new()

    pem.set_tf(1.0)  # We do not want to evolve through the singularity
    pem.prepare(vexp)

    t_a, _smaller_abst = pem.get_time_array()

    for t in t_a:
        state = pem.eval_at(vexp, t, state)
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
