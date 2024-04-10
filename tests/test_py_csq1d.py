#!/usr/bin/env python
#
# test_csq1d.py
#
# thu Dec 28 14:47:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_csq1d.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcmCSQ1D class."""
import math
import pytest

from numpy.testing import assert_allclose
from scipy.special import hankel1e, jv, yv  # pylint: disable=no-name-in-module
import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


class BesselTest(Ncm.CSQ1D):
    """Test class for NcmCSQ1D."""

    def __init__(self, alpha=2.0, k=1.0, adiab=True):
        """Constructor for BesselTest class."""
        Ncm.CSQ1D.__init__(self)
        self.alpha = alpha
        self.k = k
        if adiab:
            self.t_sign = -1.0
        else:
            self.t_sign = 1.0

    def set_k(self, k):
        """Set k parameter."""
        self.k = k

    def get_k(self):
        """Get k parameter."""
        return self.k

    def do_eval_m(self, _model, x):  # pylint: disable=arguments-differ
        """Evaluate m function, m = (-x)**(1+2alpha)."""
        return (self.t_sign * x) ** (1.0 + 2.0 * self.alpha)

    def do_eval_nu(self, _model, _x):  # pylint: disable=arguments-differ
        """Evaluate nu2 function, nu = k."""
        return self.k

    def do_eval_nu2(self, _model, _x):  # pylint: disable=arguments-differ
        """Evaluate nu2 function, nu2 = k**2."""
        return self.k**2

    def do_eval_xi(self, _model, x):  # pylint: disable=arguments-differ
        """Evaluate xi function, xi = ln(m*nu)."""
        return math.log(self.k) + (1.0 + 2.0 * self.alpha) * math.log(self.t_sign * x)

    def do_eval_F1(self, _model, x):  # pylint: disable=arguments-differ
        """Evaluate F1 function, F1 = xi'/(2nu)."""
        return 0.5 * (1.0 + 2.0 * self.alpha) / (x * self.k)

    def do_eval_F2(self, _model, x):  # pylint: disable=arguments-differ
        """Evaluate F2 function, F2 = F1'/(2nu)."""
        return -0.25 * (1.0 + 2.0 * self.alpha) / (x * self.k) ** 2

    def do_eval_int_1_m(self, _model, t):  # pylint: disable=arguments-differ
        """Evaluate int_1_m function, int_1_m."""
        return (
            -self.t_sign * (t * self.t_sign) ** (-2.0 * self.alpha) / (2.0 * self.alpha)
        )

    def do_eval_int_mnu2(self, _model, t):  # pylint: disable=arguments-differ
        """Evaluate int_mnu2 function, int_mnu2."""
        return (self.k**2 * t * (t * self.t_sign) ** (2.0 * self.alpha + 1.0)) / (
            2.0 * (self.alpha + 1.0)
        )

    def do_eval_int_qmnu2(self, _model, t):  # pylint: disable=arguments-differ
        """Evaluate int_qmnu2 function, int_qmnu2."""
        return -((self.k * t) ** 2 / (4.0 * self.alpha))

    def do_eval_int_q2mnu2(self, _model, t):  # pylint: disable=arguments-differ
        """Evaluate int_q2mnu2 function, int_q2mnu2."""
        return ((self.k * t) ** 2 * (t * self.t_sign) ** (-2.0 * self.alpha)) / (
            8.0 * self.t_sign * self.alpha**2 * (1.0 - self.alpha)
        )

    def do_prepare(self, _model):  # pylint: disable=arguments-differ
        pass


def test_csq1d():
    """Test basic functionality of NcmCSQ1D."""

    bs = BesselTest(alpha=2.0)
    bs.set_k(1.0)

    bs.set_ti(-100.0)
    bs.set_tf(-1.0e-3)

    assert_allclose(bs.get_k(), 1.0)
    assert_allclose(bs.get_ti(), -100.0)
    assert_allclose(bs.get_tf(), -1.0e-3)

    bs.set_reltol(1.0e-6)
    bs.set_abstol(1.0e-7)

    assert_allclose(bs.get_reltol(), 1.0e-6)
    assert_allclose(bs.get_abstol(), 1.0e-7)

    bs.set_adiab_threshold(1.0e-3)
    bs.set_prop_threshold(1.0e-3)

    assert_allclose(bs.get_adiab_threshold(), 1.0e-3)
    assert_allclose(bs.get_prop_threshold(), 1.0e-3)

    bs.set_save_evol(True)
    assert bs.get_save_evol()

    bs.set_save_evol(False)
    assert not bs.get_save_evol()

    bs.set_initial_condition_type(Ncm.CSQ1DInitialStateType.AD_HOC)
    assert bs.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.AD_HOC

    bs.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC2)
    assert bs.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC2

    bs.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)
    assert bs.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.ADIABATIC4

    bs.set_initial_condition_type(Ncm.CSQ1DInitialStateType.NONADIABATIC2)
    assert bs.get_initial_condition_type() == Ncm.CSQ1DInitialStateType.NONADIABATIC2


def test_initial_conditions_time():
    """Test initial conditions of NcmCSQ1D."""

    bs = BesselTest(alpha=2.0)
    bs.set_k(1.0)
    bs.set_ti(-100.0)
    bs.set_tf(-1.0e-3)

    limit_found, t_adiab = bs.find_adiab_time_limit(None, -1.0e3, -1.0e1, 1.0e-6)

    assert limit_found
    assert t_adiab >= bs.get_ti()
    assert t_adiab <= bs.get_tf()

    t_min, F1_min, t_lb, t_ub = bs.find_adiab_max(None, -1.0e3, -1.0e1, 1.0e-1)

    assert_allclose(F1_min, bs.eval_F1(None, t_min))
    assert math.fabs(F1_min - bs.eval_F1(None, t_lb)) <= 1.0e-1
    assert math.fabs(F1_min - bs.eval_F1(None, t_ub)) <= 1.0e-1


def test_initial_conditions_adiabatic():
    """Test initial conditions of NcmCSQ1D."""

    bs = BesselTest(alpha=2.0)
    bs.set_k(1.0)
    bs.set_ti(-100.0)
    bs.set_tf(-1.0e-3)
    state = Ncm.CSQ1DState.new()

    for prec in np.geomspace(1.0e-14, 1.0e-6, 100):
        limit_found, t_adiab = bs.find_adiab_time_limit(None, -1.0e4, -1.0e1, prec)

        assert limit_found

        # Getting the adiabatic solution
        state, _alpha_reltol, _dgamma_reltol = bs.compute_adiab(None, t_adiab, state)
        bs.change_frame(None, state, Ncm.CSQ1DFrame.ORIG)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        kt = bs.get_k() * t_adiab
        hfnormm = 0.5 * math.sqrt(math.pi) * (-t_adiab) ** (-bs.alpha)
        hfnormp = 0.5 * math.sqrt(math.pi) * (-t_adiab) ** (+bs.alpha)

        # Analytical solution for phi and Pphi
        theo_phi = hfnormm * hankel1e(bs.alpha, kt)
        theo_Pphi = kt * hfnormp * hankel1e(1.0 + bs.alpha, kt)

        # Compare with analytical solution
        assert_allclose(abs(phi), abs(theo_phi), rtol=1.0e-6)
        assert_allclose(abs(Pphi), abs(theo_Pphi), rtol=1.0e-6)


def test_evolution():
    """Test initial conditions of NcmCSQ1D."""

    bs = BesselTest(alpha=2.0)
    bs.set_k(1.0)
    bs.set_ti(-100.0)
    bs.set_tf(-1.0e-3)
    state = Ncm.CSQ1DState.new()

    limit_found, t_adiab = bs.find_adiab_time_limit(None, -1.0e4, -1.0e1, 1.0e-8)

    assert limit_found

    bs.set_save_evol(True)
    bs.set_reltol(1.0e-10)
    bs.set_abstol(0.0)
    bs.set_init_cond_adiab(None, t_adiab)
    bs.prepare(None)

    t_a, _smaller_abst = bs.get_time_array()

    for t in t_a:
        state = bs.eval_at(None, t, state)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        kt = bs.get_k() * t
        hfnormm = 0.5 * math.sqrt(math.pi) * (-t) ** (-bs.alpha)
        hfnormp = 0.5 * math.sqrt(math.pi) * (-t) ** (+bs.alpha)

        # Analytical solution for phi and Pphi
        theo_phi = hfnormm * hankel1e(bs.alpha, kt)
        theo_Pphi = kt * hfnormp * hankel1e(1.0 + bs.alpha, kt)

        # Compare with analytical solution
        assert_allclose(abs(phi), abs(theo_phi), rtol=1.0e-7)
        assert_allclose(abs(Pphi), abs(theo_Pphi), rtol=1.0e-7)

        J11, J12, J22 = state.get_J()

        assert_allclose(J11, 2.0 * abs(phi) ** 2, atol=1.0e-7)
        assert_allclose(J22, 2.0 * abs(Pphi) ** 2, atol=1.0e-7)
        assert_allclose(J12, (2.0 * phi * Pphi.conjugate()).real, atol=1.0e-7)


def test_evolution_adiabatic2():
    """Test initial conditions of NcmCSQ1D."""

    bs = BesselTest(alpha=2.0)
    bs.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC2)
    bs.set_k(1.0)
    bs.set_ti(-1.0e5)
    bs.set_tf(-1.0e-3)
    bs.set_vacuum_max_time(-1.0e0)
    bs.set_vacuum_reltol(1.0e-3)
    state = Ncm.CSQ1DState.new()

    limit_found, t_adiab = bs.find_adiab_time_limit(None, -1.0e5, -1.0e0, 1.0e-3)

    assert limit_found

    bs.set_save_evol(True)
    bs.set_reltol(1.0e-10)
    bs.set_abstol(0.0)
    bs.set_init_cond_adiab(None, t_adiab)
    bs.prepare(None)

    t_a, _smaller_abst = bs.get_time_array()

    for t in t_a:
        state = bs.eval_at(None, t, state)
        phi_vec, Pphi_vec = state.get_phi_Pphi()

        phi = phi_vec[0] + 1.0j * phi_vec[1]
        Pphi = Pphi_vec[0] + 1.0j * Pphi_vec[1]

        kt = bs.get_k() * t
        hfnormm = 0.5 * math.sqrt(math.pi) * (-t) ** (-bs.alpha)
        hfnormp = 0.5 * math.sqrt(math.pi) * (-t) ** (+bs.alpha)

        # Analytical solution for phi and Pphi
        theo_phi = hfnormm * hankel1e(bs.alpha, kt)
        theo_Pphi = kt * hfnormp * hankel1e(1.0 + bs.alpha, kt)

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
def test_evolution_frame(frame):
    """Test initial conditions of NcmCSQ1D."""

    bs = BesselTest(alpha=2.0)
    bs.set_k(1.0)
    bs.set_ti(-100.0)
    state = Ncm.CSQ1DState.new()

    limit_found, t_adiab = bs.find_adiab_time_limit(None, -1.0e4, -1.0e1, 1.0e-8)
    assert limit_found
    limit_found, t_adiab_end = bs.find_adiab_time_limit(None, -1.0e4, -1.0e0, 1.0e-1)
    assert limit_found
    bs.set_tf(t_adiab_end)

    bs.set_save_evol(True)
    bs.set_reltol(1.0e-10)
    bs.set_abstol(0.0)
    bs.set_init_cond_adiab(None, t_adiab)
    bs.prepare(None)

    t_a, _smaller_abst = bs.get_time_array()

    for t in t_a:
        state = bs.eval_at_frame(None, frame, t, state)
        assert state.get_frame() == frame


def test_change_frame_orig_adiab1():
    """Test change_frame method of NcmCSQ1D."""

    bs = BesselTest(alpha=2.0)
    bs.set_k(1.0)
    bs.set_ti(-100.0)
    bs.set_tf(-1.0e-3)

    t = -2.0e1
    state = Ncm.CSQ1DState()
    state.set_ag(Ncm.CSQ1DFrame.ORIG, t, 0.1, 0.3)

    bs.change_frame(None, state, Ncm.CSQ1DFrame.ORIG)

    assert state.get_frame() == Ncm.CSQ1DFrame.ORIG
    assert state.get_ag() == (0.1, 0.3)

    bs.change_frame(None, state, Ncm.CSQ1DFrame.ADIAB1)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB1
    assert_allclose(state.get_ag(), (0.1, 0.3 - bs.eval_xi(None, t)))

    bs.change_frame(None, state, Ncm.CSQ1DFrame.ORIG)
    assert state.get_frame() == Ncm.CSQ1DFrame.ORIG
    assert_allclose(state.get_ag(), (0.1, 0.3))


def test_change_frame_orig_adiab2():
    """Test change_frame method of NcmCSQ1D."""

    bs = BesselTest(alpha=2.0)
    bs.set_k(1.0)
    bs.set_ti(-100.0)
    bs.set_tf(-1.0e-3)

    t = -2.0e1
    state = Ncm.CSQ1DState()
    state.set_ag(Ncm.CSQ1DFrame.ORIG, t, 0.1, 0.3)

    bs.change_frame(None, state, Ncm.CSQ1DFrame.ORIG)

    assert state.get_frame() == Ncm.CSQ1DFrame.ORIG
    assert state.get_ag() == (0.1, 0.3)

    bs.change_frame(None, state, Ncm.CSQ1DFrame.ADIAB2)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB2
    assert not np.allclose(state.get_ag(), (0.1, 0.3 - bs.eval_xi(None, t)))

    bs.change_frame(None, state, Ncm.CSQ1DFrame.ORIG)
    assert state.get_frame() == Ncm.CSQ1DFrame.ORIG
    assert_allclose(state.get_ag(), (0.1, 0.3))


def test_change_frame_adiab1_adiab2():
    """Test change_frame method of NcmCSQ1D."""

    bs = BesselTest(alpha=2.0)
    bs.set_k(1.0)
    bs.set_ti(-100.0)
    bs.set_tf(-1.0e-3)

    t = -2.0e1
    state = Ncm.CSQ1DState()
    state.set_ag(Ncm.CSQ1DFrame.ADIAB1, t, 0.1, 0.3)

    bs.change_frame(None, state, Ncm.CSQ1DFrame.ADIAB1)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB1
    assert state.get_ag() == (0.1, 0.3)

    bs.change_frame(None, state, Ncm.CSQ1DFrame.ADIAB2)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB2
    assert not np.allclose(state.get_ag(), (0.1, 0.3))

    bs.change_frame(None, state, Ncm.CSQ1DFrame.ADIAB1)
    assert state.get_frame() == Ncm.CSQ1DFrame.ADIAB1
    assert_allclose(state.get_ag(), (0.1, 0.3))


def test_nonadiab_prop():
    """Test basic functionality of NcmCSQ1D."""

    k = 8.0
    alpha = 0.5
    bs = BesselTest(alpha=alpha, adiab=False)
    bs.set_k(k)

    ti = 0.0
    tii = 1.0e-7
    tf = 2.0

    bs.prepare_prop(None, ti, tii, tf)
    state0 = Ncm.CSQ1DState.new()
    state1 = Ncm.CSQ1DState.new()
    state0.set_up(Ncm.CSQ1DFrame.NONADIAB1, ti, 0.0, 0.0)

    def theo_phi(t):
        """Theoretical phi."""
        prefactor = 0.5 * (1.0 - 1.0j) * math.sqrt(math.pi / 2.0) * (k * t) ** (-alpha)
        return prefactor * (
            jv(alpha, k * t) - 1.0j * k ** (2.0 * alpha) * yv(alpha, k * t)
        )

    test_t = np.geomspace(tii, 1.0e-2, 10)

    prop_J11 = [
        bs.evolve_prop_vector(None, state0, Ncm.CSQ1DFrame.ORIG, t, state1).get_J()[0]
        for t in test_t
    ]

    analytic_J11 = [2.0 * np.abs(theo_phi(t)) ** 2 for t in test_t]

    assert_allclose(prop_J11, analytic_J11, rtol=1.0e-9)


def test_nonadiab_evol():
    """Test basic functionality of NcmCSQ1D."""

    k = 8.0
    alpha = 0.5
    bs = BesselTest(alpha=alpha, adiab=False)
    bs.set_k(k)

    ti = 0.0
    tii = 1.0e-7
    tf = 2.0

    bs.prepare_prop(None, ti, tii, tf)
    state0 = Ncm.CSQ1DState.new()
    state1 = Ncm.CSQ1DState.new()
    state = Ncm.CSQ1DState.new()
    state0.set_up(Ncm.CSQ1DFrame.NONADIAB1, ti, 0.0, 0.0)

    def theo_phi(t):
        """Theoretical phi."""
        prefactor = 0.5 * (1.0 - 1.0j) * math.sqrt(math.pi / 2.0) * (k * t) ** (-alpha)
        return prefactor * (
            jv(alpha, k * t) - 1.0j * k ** (2.0 * alpha) * yv(alpha, k * t)
        )

    bs.evolve_prop_vector(None, state0, Ncm.CSQ1DFrame.ORIG, 1.0e-4, state1)
    bs.set_init_cond(None, Ncm.CSQ1DEvolState.UP, state1)
    bs.set_tf(10.0)
    bs.set_reltol(1.0e-14)
    bs.prepare(None)
    t_a, _smaller_abst = bs.get_time_array()

    evol_J11 = [bs.eval_at(None, t, state).get_J()[0] for t in t_a]
    analytic_J11 = [2.0 * np.abs(theo_phi(t)) ** 2 for t in t_a]

    assert_allclose(evol_J11, analytic_J11, rtol=1.0e-7)


def test_state_circle():
    """Test circle functionality of NcmCSQ1DState."""

    state = Ncm.CSQ1DState.new()
    state.set_up(Ncm.CSQ1DFrame.ORIG, 1.0, 0.1, 0.2)

    for r in np.geomspace(1.0e-4, 1.0e4, 100):
        for theta in np.linspace(-4.0 * math.pi, 4.0 * math.pi, 120):
            state1 = state.get_circle(r, theta)
            assert_allclose(state.compute_distance(state1), r)


def test_state_phi_Pphi():
    """Test phi and Pphi functionality of NcmCSQ1DState."""

    state = Ncm.CSQ1DState.new()
    alpha_a = np.linspace(-2.0, 2.0, 100)
    gamma_a = np.linspace(-2.0, 2.0, 100)

    for alpha, gamma in zip(alpha_a, gamma_a):
        state.set_ag(Ncm.CSQ1DFrame.ORIG, 0.0, alpha, gamma)

        phi_v, Pphi_v = state.get_phi_Pphi()
        J11, J12, J22 = state.get_J()

        phi = phi_v[0] + 1.0j * phi_v[1]
        Pphi = Pphi_v[0] + 1.0j * Pphi_v[1]
        phic = phi_v[0] - 1.0j * phi_v[1]
        Pphic = Pphi_v[0] - 1.0j * Pphi_v[1]

        assert_allclose(2.0 * phi * phic, J11)
        assert_allclose(2.0 * Pphi * Pphic, J22)
        assert_allclose(phi * Pphic + phic * Pphi, J12)
        assert_allclose(phi * Pphic - phic * Pphi, 1.0j)


def test_state_um():
    """Test um functionality of NcmCSQ1DState."""

    state = Ncm.CSQ1DState.new()
    alpha_a = np.linspace(-2.0, 2.0, 100)
    gamma_a = np.linspace(-2.0, 2.0, 100)

    for alpha, gamma in zip(alpha_a, gamma_a):
        state.set_ag(Ncm.CSQ1DFrame.ORIG, 0.0, alpha, gamma)

        chi, Um = state.get_um()
        alpha1, gamma1 = state.get_ag()

        state.set_um(Ncm.CSQ1DFrame.ORIG, 0.0, chi, Um)
        alpha2, gamma2 = state.get_ag()

        assert_allclose(alpha, alpha1)
        assert_allclose(alpha, alpha2)
        assert_allclose(gamma, gamma1)
        assert_allclose(gamma, gamma2)
        assert_allclose(chi, np.sinh(alpha))
        assert_allclose(Um, -gamma + np.log(np.cosh(alpha)))


def test_state_up():
    """Test up functionality of NcmCSQ1DState."""

    state = Ncm.CSQ1DState.new()
    alpha_a = np.linspace(-2.0, 2.0, 100)
    gamma_a = np.linspace(-2.0, 2.0, 100)

    for alpha, gamma in zip(alpha_a, gamma_a):
        state.set_ag(Ncm.CSQ1DFrame.ORIG, 0.0, alpha, gamma)

        chi, Up = state.get_up()
        alpha1, gamma1 = state.get_ag()

        state.set_up(Ncm.CSQ1DFrame.ORIG, 0.0, chi, Up)
        alpha2, gamma2 = state.get_ag()

        assert_allclose(alpha, alpha1)
        assert_allclose(alpha, alpha2)
        assert_allclose(gamma, gamma1)
        assert_allclose(gamma, gamma2)
        assert_allclose(chi, np.sinh(alpha))
        assert_allclose(Up, gamma + np.log(np.cosh(alpha)))


def test_state_poincare_half_plane():
    """Test poincare_half_plane functionality of NcmCSQ1DState."""

    state = Ncm.CSQ1DState.new()
    alpha_a = np.linspace(-2.0, 2.0, 100)
    gamma_a = np.linspace(-2.0, 2.0, 100)

    for alpha, gamma in zip(alpha_a, gamma_a):
        state.set_ag(Ncm.CSQ1DFrame.ORIG, 0.0, alpha, gamma)

        x, lny = state.get_poincare_half_plane()
        chi, Up = state.get_up()

        assert_allclose(x, chi * np.exp(-Up))
        assert_allclose(lny, -Up)


def test_state_poincare_disc():
    """Test poincare_disc functionality of NcmCSQ1DState."""

    state = Ncm.CSQ1DState.new()
    alpha_a = np.linspace(-2.0, 2.0, 100)
    gamma_a = np.linspace(-2.0, 2.0, 100)

    for alpha, gamma in zip(alpha_a, gamma_a):
        state.set_ag(Ncm.CSQ1DFrame.ORIG, 0.0, alpha, gamma)

        x, y = state.get_poincare_disc()
        chi, Up = state.get_up()
        chi, Um = state.get_um()

        assert_allclose(x, chi / (1.0 + 0.5 * (np.exp(Up) + np.exp(Um))))
        assert_allclose(
            y,
            -0.5 * (np.exp(Up) - np.exp(Um)) / (1.0 + 0.5 * (np.exp(Up) + np.exp(Um))),
        )


if __name__ == "__main__":
    test_csq1d()
    test_initial_conditions_time()
    test_initial_conditions_adiabatic()
    test_evolution()
    test_evolution_adiabatic2()
    for frame0 in [Ncm.CSQ1DFrame.ORIG, Ncm.CSQ1DFrame.ADIAB1, Ncm.CSQ1DFrame.ADIAB2]:
        test_evolution_frame(frame0)
    test_change_frame_orig_adiab1()
    test_change_frame_orig_adiab2()
    test_change_frame_adiab1_adiab2()
    test_nonadiab_prop()
    test_nonadiab_evol()
    test_state_circle()
    test_state_phi_Pphi()
    test_state_um()
    test_state_up()
