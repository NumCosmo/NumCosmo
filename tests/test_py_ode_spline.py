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

"""Tests for OdeSpline object to create a spline from an ODE."""

from numpy.testing import assert_allclose
import numpy as np
from numcosmo_py import Ncm

Ncm.cfg_init()


class RhsF:
    """Test function for the ODE, it describes the right-hand
    side of the equation $y_x = f(y,x)$."""

    def __init__(self, alpha):
        self.alpha = alpha

    def rhs(self, y, _x, _data):
        r"""Test function for the ODE, it describes the right-hand
        side of the equation $y_x = f(y,x)$, here we use $f(y,x) = \alpha y$."""
        return self.alpha * y


def test_ode_spline() -> None:
    """Test function for the ODE spline."""

    f = RhsF(2.5)

    s = Ncm.SplineCubicNotaknot.new()

    os = Ncm.OdeSpline.new(s, f.rhs)
    os.set_reltol(1.0e-10)

    x_0 = 0.0
    y_0 = 1.0
    x_1 = 5.0

    os.props.xi = x_0
    os.props.xf = x_1
    os.props.yi = y_0

    os.prepare(None)

    ss = os.peek_spline()
    assert isinstance(ss, Ncm.SplineCubic)

    x_a = np.linspace(x_0, x_1, 100)

    assert_allclose([ss.eval(x) for x in x_a], np.exp(x_a * f.alpha), rtol=1e-7)


def test_ode_spline_copy_props() -> None:
    """Test function for the ODE spline."""

    f = RhsF(2.5)

    s = Ncm.SplineCubicNotaknot.new()

    os0 = Ncm.OdeSpline.new(s, f.rhs)
    os0.set_reltol(1.0e-10)

    x_0 = 0.0
    y_0 = 1.0
    x_1 = 5.0

    os0.props.xi = x_0
    os0.props.xf = x_1
    os0.props.yi = y_0

    os = Ncm.OdeSpline.new(s, f.rhs)
    for prop in os0.list_properties():
        if prop.name == "dydx" or prop.name == "yf":
            continue
        print(prop.name)
        os.set_property(prop.name, os0.get_property(prop.name))

    assert isinstance(os, Ncm.OdeSpline)

    os.prepare(None)

    ss = os.peek_spline()
    assert isinstance(ss, Ncm.SplineCubic)

    x_a = np.linspace(x_0, x_1, 100)

    assert_allclose([ss.eval(x) for x in x_a], np.exp(x_a * f.alpha), rtol=1e-7)


if __name__ == "__main__":
    test_ode_spline()
    test_ode_spline_copy_props()
