#!/usr/bin/env python
#
# example_ode_spline.py
#
# Fri May 19 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_ode_spline.py
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

"""Example of using the OdeSpline object to create a spline from an ODE."""

import math
import numpy as np
from numcosmo_py import Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


class RhsF:
    """Test function for the ODE, it describes the right-hand
    side of the equation $y_x = f(y,x)$."""

    def __init__(self, alpha):
        self.alpha = alpha

    # C closure _data is not supported in Python, so we use a class
    # to store the data.
    def rhs(self, y, _x, _data):
        r"""Test function for the ODE, it describes the right-hand
        side of the equation $y_x = f(y,x)$, here we use $f(y,x) = \alpha y$."""
        return self.alpha * y


def test_ode_spline() -> None:
    """Test function for the ODE spline."""

    f = RhsF(2.5)

    s = Ncm.SplineCubicNotaknot.new()

    os = Ncm.OdeSpline.new(s, f.rhs)
    os.set_reltol(1.0e-3)

    x_0 = 0.0
    y_0 = 1.0
    x_1 = 5.0

    os.props.xi = x_0
    os.props.xf = x_1
    os.props.yi = y_0

    # Prepare the ODE solver to calculate the spline.
    # Since C closures are not supported in Python, we need to pass
    # None as the data argument.
    os.prepare(None)

    # Gets the results from the ODE solver.
    ss = os.peek_spline()
    assert isinstance(ss, Ncm.SplineCubic)

    # Print the spline coefficients.
    for i in range(ss.len):
        print(
            f"{i} {ss.xv.get(i): 22.15g} {ss.yv.get(i): 22.15g} "
            f"{ss.b.get(i): 22.15g} {ss.c.get(i): 22.15g} {ss.d.get(i): 22.15g}"
        )

    # Evaluate the spline at some points.
    # Compare with the exponential function.
    x_a = np.linspace(x_0, x_1, 100)
    for x in x_a:
        expx = math.exp(x * f.alpha)
        odex = ss.eval(x)
        print(
            f"x = {x: 8.5}, exp(x) = {expx: 22.15g}, ode(x) = {odex: 22.15g}, "
            f"relative difference {math.fabs((expx - odex) / expx): 22.15e}"
        )


if __name__ == "__main__":
    test_ode_spline()
