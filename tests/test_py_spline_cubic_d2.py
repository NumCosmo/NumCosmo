#!/usr/bin/env python
#
# test_py_spline_cubic_d2.py
#
# Fri May 19 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_spline_cubic_d2.py
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

"""Tests on NcmSplineCubicD2 class."""

import math
import pytest
import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


def f_x2(x: float) -> float:
    """Function x^2."""
    return x**2


def df_x2(x: float) -> float:
    """Derivative of function x^2."""
    return 2.0 * x


def d2f_x2(_: float) -> float:
    """Second derivative of function x^2."""
    return 2.0


def f_x3(x: float) -> float:
    """Function x^3."""
    return x**3


def df_x3(x: float) -> float:
    """Derivative of function x^3."""
    return 3.0 * x**2


def d2f_x3(x: float) -> float:
    """Second derivative of function x^3."""
    return 6.0 * x


def f_cos(x: float) -> float:
    """Function cos(x)."""
    return math.cos(x)


def df_cos(x: float) -> float:
    """Derivative of function cos(x)."""
    return -math.sin(x)


def d2f_cos(x: float) -> float:
    """Second derivative of function cos(x)."""
    return -math.cos(x)


@pytest.mark.parametrize(
    "f,df,d2f,interp_nknots,tol",
    [
        (f_x2, df_x2, d2f_x2, 10, (1.0e-10, 1.0e-10)),
        (f_x3, df_x3, d2f_x3, 12, (1.0e-10, 1.0e-10)),
        (f_cos, df_cos, d2f_cos, 1000, (1.0e-5, 1.0e-5)),
    ],
)
def test_func_x2(f, df, d2f, interp_nknots, tol) -> None:
    """Test the numcosmo library to calculate derivatives of functions."""

    interp_knots = np.linspace(0.0, 4.0, interp_nknots)

    s = Ncm.SplineCubicD2.new(
        Ncm.Vector.new_array(interp_knots),
        Ncm.Vector.new_array([f(x) for x in interp_knots]),
        Ncm.Vector.new_array([d2f(x) for x in interp_knots]),
        True,
    )

    knots = np.linspace(0.0, 4.0, interp_nknots * 10)
    f_interp = [s.eval(x) for x in knots]
    f_true = [f(x) for x in knots]

    df_interp = [s.eval_deriv(x) for x in knots]
    df_true = [df(x) for x in knots]

    d2f_interp = [s.eval_deriv2(x) for x in knots]
    d2f_true = [d2f(x) for x in knots]

    np.testing.assert_allclose(f_interp, f_true, rtol=tol[0], atol=tol[1])
    np.testing.assert_allclose(df_interp, df_true, rtol=tol[0], atol=tol[1])

    d2f_interp = [s.eval_deriv2(x) for x in knots]
    np.testing.assert_allclose(d2f_interp, d2f_true, rtol=tol[0], atol=tol[1])
