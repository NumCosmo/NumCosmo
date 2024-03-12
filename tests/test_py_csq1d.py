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
from numpy.testing import assert_allclose
from scipy.special import hankel1e  # pylint: disable=no-name-in-module
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

    bs.set_max_order_2(True)
    assert bs.get_max_order_2()

    bs.set_max_order_2(False)
    assert not bs.get_max_order_2()


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
        state = bs.eval_at(t, state)
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


if __name__ == "__main__":
    test_csq1d()
    test_initial_conditions_time()
    test_initial_conditions_adiabatic()
    test_evolution()
    test_change_frame_orig_adiab1()
    test_change_frame_orig_adiab2()
    test_change_frame_adiab1_adiab2()
