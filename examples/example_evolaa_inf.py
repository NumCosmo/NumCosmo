#!/usr/bin/env python
#
# example_evolaa_inf.py
#
# Sun May 21 21:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_evolaa_inf.py
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

"""Example of using the HOAA ODE solver a scalar field model.
"""

import math
import numpy as np
from numcosmo_py import Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


class PyHOAATest(Ncm.HOAA):
    """Python implementation of the HOAA ODE solver for a scalar field model."""

    def __init__(self):
        Ncm.HOAA.__init__(self, opt=Ncm.HOAAOpt.DLNMNU_ONLY)

    def do_eval_nu(self, _model, _t, k):
        """Evaluate the nu function."""
        return k

    def do_eval_m(self, _model, t, _k):
        """Evaluate the m function."""
        return 1.0 / (t * t)

    def do_eval_mnu(self, _model, t, k):
        """Evaluate the mnu function."""
        return k / (t * t)

    def do_eval_dlnmnu(self, _model, t, _k):
        """Evaluate the dlnmnu function."""
        return -2.0 / t

    def do_eval_system(self, _model, t, k):
        """Evaluate the system function."""
        return k, -2.0 / t, 0.0

    def do_nsing(self, _model, _k):
        """Evaluate the nsing function."""
        return 1

    def do_get_sing_info(self, _model, _k, _sing):
        """Evaluate the get_sing_info function."""

        return 0.0, -1.0, +1.0, Ncm.HOAASingType.INF

    def do_eval_sing_mnu(self, _model, t, k, _sing):
        """Evaluate the sing_mnu function."""
        return k / (t * t)

    def do_eval_sing_dlnmnu(self, _model, t, _k, _sing):
        """Evaluate the sing_dlnmnu function."""
        return -2.0 / t

    def do_eval_sing_system(self, _model, t, k, _sing):
        """Evaluate the sing_system function."""
        return k, -2.0 / t, 0.0


def sol_q(k, t):
    """Analytical solution for the q function."""
    a = k * t
    return (-t / k**2) * (math.cos(a) - math.sin(a) / a)


def sol_p(k, t):
    """Analytical solution for the p function."""
    a = k * t
    return math.sin(a) / a


def test():
    """Test the HOAA ODE solver for a scalar field model."""
    hoaa = PyHOAATest()

    k = 1.0
    hoaa.set_ti(-1.0e10)
    hoaa.set_tf(+1.0e10)
    hoaa.set_k(k)
    hoaa.set_reltol(1.0e-14)

    ti = -10.0

    S1 = sol_q(k, ti)
    PS1 = sol_p(k, ti)

    hoaa.prepare()

    (t0, t1) = hoaa.get_t0_t1()

    (Aq, Av) = hoaa.eval_solution(None, ti, S1, PS1)

    print("# ", t0, t1)
    print("# ", Aq, Av)

    ta = np.linspace(-42.0, -39.00, 100000)

    for t in ta:
        (q, v, Pq, Pv) = hoaa.eval_QV(None, t)
        (epsilon, _gamma, _sin_thetab, _cos_thetab) = hoaa.eval_AA(None, t)

        S = Aq * q + Av * v
        PS = Aq * Pq + Av * Pv

        mnu = hoaa.eval_mnu(None, t, k)
        _nu = hoaa.eval_nu(None, t, k)

        I = 0.5 * (mnu * q**2 + Pq**2 / mnu)
        J = 0.5 * (mnu * v**2 + Pv**2 / mnu)

        print(
            t,
            S,
            sol_q(k, t),
            PS,
            sol_p(k, t),
            I,
            J,
            math.sqrt(I * J),
            math.cosh(epsilon),
        )


if __name__ == "__main__":
    test()
