#!/usr/bin/env python
#
# test_py_csq1d.py
#
# thu Dec 28 14:47:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_csq1d.py
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
import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


class BesselTest(Ncm.CSQ1D):
    """Test class for NcmCSQ1D."""

    def __init__(self, alpha=2.0, k=1.0, adiab=True):
        """Initialize BesselTest."""
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
        """Prepare method, nothing to do."""


class BesselTestWithIntNu(BesselTest):
    """BesselTest subclass with analytic eval_int_nu override."""

    def do_eval_int_nu(self, _model, t):  # pylint: disable=arguments-differ
        """Evaluate the integral of nu analytically: k * (t - ti)."""
        return self.k * (t - self.get_ti())


def test_eval_int_nu_override():
    """Test that an overridden eval_int_nu is used in place of the ODE spline."""
    bs = BesselTestWithIntNu(alpha=2.0)
    k = 1.0
    ti = -1.0e4

    bs.set_k(k)
    bs.set_ti(ti)
    bs.set_tf(-1.0e-3)

    # The override should work without calling prepare()
    test_times = np.linspace(ti, -1.0e-3, 20)[1:-1]
    for t in test_times:
        assert_allclose(bs.eval_int_nu(None, t), k * (t - ti))

    # After prepare(), the override should still be used (not the ODE spline)
    bs.set_save_evol(True)
    bs.set_reltol(1.0e-10)
    bs.set_abstol(0.0)
    bs.set_initial_condition_type(Ncm.CSQ1DInitialStateType.ADIABATIC4)
    bs.set_vacuum_max_time(-1.0e1)
    bs.set_vacuum_reltol(1.0e-8)
    bs.prepare(None)

    for t in test_times:
        assert_allclose(bs.eval_int_nu(None, t), k * (t - ti))


if __name__ == "__main__":
    test_eval_int_nu_override()
