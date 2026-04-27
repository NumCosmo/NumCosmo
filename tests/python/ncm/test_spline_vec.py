#!/usr/bin/env python
#
# test_py_spline_vec.py
#
# Sat Mar 15 19:53:22 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_spline_vec.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcmSplineVec class."""

import math
import pytest
import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


def f0(x: float) -> float:
    """Compute x^2."""
    return x**2


def f1(x: float) -> float:
    """Compute x^3."""
    return x**3


def f2(x: float) -> float:
    """Compute cos(x)."""
    return math.cos(x)


def df0(x: float) -> float:
    """Compute the derivative of x^2."""
    return 2.0 * x


def df1(x: float) -> float:
    """Compute the derivative of x^3."""
    return 3.0 * x**2


def df2(x: float) -> float:
    """Compute the derivative of cos(x)."""
    return -math.sin(x)


def int_f0(x0: float, x1: float) -> float:
    """Integrate x^2 from x0 to x1."""
    return (x1**3 - x0**3) / 3.0


def int_f1(x0: float, x1: float) -> float:
    """Integrate x^3 from x0 to x1."""
    return (x1**4 - x0**4) / 4.0


def int_f2(x0: float, x1: float) -> float:
    """Integrate cos(x) from x0 to x1."""
    return math.sin(x1) - math.sin(x0)


@pytest.fixture(name="spline_vec")
def fixture_spline_vec() -> Ncm.SplineVec:
    """Fixture for NcmSplineVec with three components."""
    nknots = 100
    xi = 0.0
    xf = 10.0

    # Create x vector
    x_arr = np.linspace(xi, xf, nknots)
    xv = Ncm.Vector.new_array(x_arr.tolist())

    # Create y matrix with 3 components: x^2, x^3, cos(x)
    ym = Ncm.Matrix.new(3, nknots)
    for i, x in enumerate(x_arr):
        ym.set(0, i, f0(x))
        ym.set(1, i, f1(x))
        ym.set(2, i, f2(x))

    # Create spline vector
    s = Ncm.SplineCubicNotaknot.new()
    sv = Ncm.SplineVec.new(s, xv, ym, True)

    return sv


@pytest.fixture(name="spline_vec_gpa")
def fixture_spline_vec_gpa() -> Ncm.SplineVec:
    """Fixture for NcmSplineVec using GPtrArray constructor."""
    nknots = 100
    xi = 0.0
    xf = 10.0

    # Create x vector
    x_arr = np.linspace(xi, xf, nknots)
    xv = Ncm.Vector.new_array(x_arr.tolist())

    # Create y vectors in a list
    yv_list = []
    y0 = Ncm.Vector.new(nknots)
    y1 = Ncm.Vector.new(nknots)
    y2 = Ncm.Vector.new(nknots)

    for i, x in enumerate(x_arr):
        y0.set(i, f0(x))
        y1.set(i, f1(x))
        y2.set(i, f2(x))

    yv_list.append(y0)
    yv_list.append(y1)
    yv_list.append(y2)

    # Create spline vector
    s = Ncm.SplineCubicNotaknot.new()
    sv = Ncm.SplineVec.new_gpa(s, xv, yv_list, True)

    return sv


def test_spline_vec_new(spline_vec: Ncm.SplineVec) -> None:
    """Test creation of NcmSplineVec."""
    assert isinstance(spline_vec, Ncm.SplineVec)
    assert spline_vec.get_len() == 3
    assert spline_vec.is_init()


def test_spline_vec_new_gpa(spline_vec_gpa: Ncm.SplineVec) -> None:
    """Test creation of NcmSplineVec using GPtrArray."""
    assert isinstance(spline_vec_gpa, Ncm.SplineVec)
    assert spline_vec_gpa.get_len() == 3
    assert spline_vec_gpa.is_init()


def test_spline_vec_get_nknots(spline_vec: Ncm.SplineVec) -> None:
    """Test get_nknots method."""
    # The fixture creates a spline with 100 knots
    assert spline_vec.get_nknots() == 100

    # Verify it matches the underlying spline length
    s0 = spline_vec.peek_spline(0)
    assert spline_vec.get_nknots() == s0.get_len()


def test_spline_vec_get_nknots_gpa(spline_vec_gpa: Ncm.SplineVec) -> None:
    """Test get_nknots method with GPtrArray constructor."""
    # The fixture creates a spline with 100 knots
    assert spline_vec_gpa.get_nknots() == 100

    # All component splines should have the same number of knots
    for i in range(spline_vec_gpa.get_len()):
        s_i = spline_vec_gpa.peek_spline(i)
        assert spline_vec_gpa.get_nknots() == s_i.get_len()


def test_spline_vec_eval(spline_vec: Ncm.SplineVec) -> None:
    """Test evaluation of NcmSplineVec."""
    x = 5.0

    res = spline_vec.eval_array(x)

    # Check that results match individual spline evaluations
    for i in range(3):
        s_i = spline_vec.peek_spline(i)
        expected = s_i.eval(x)
        computed = res[i]
        assert abs(computed - expected) < 1e-15


def test_spline_vec_eval_values(spline_vec: Ncm.SplineVec) -> None:
    """Test that evaluated values are close to expected function values."""
    x = 5.0

    res = spline_vec.eval_array(x)

    # For well-sampled cubic splines, interpolation should be very accurate
    assert abs(res[0] - f0(x)) < 1e-10  # x^2 is exactly interpolated
    assert abs(res[1] - f1(x)) < 1e-10  # x^3 is exactly interpolated
    assert abs(res[2] - f2(x)) < 1e-4  # cos(x) has interpolation error


def test_spline_vec_deriv(spline_vec: Ncm.SplineVec) -> None:
    """Test derivative of NcmSplineVec."""
    x = 5.0

    res = spline_vec.deriv_array(x)

    # Check that results match individual spline derivatives
    for i in range(3):
        s_i = spline_vec.peek_spline(i)
        expected = s_i.eval_deriv(x)
        computed = res[i]
        assert abs(computed - expected) < 1e-15


def test_spline_vec_deriv_values(spline_vec: Ncm.SplineVec) -> None:
    """Test that derivative values are close to expected derivative values."""
    x = 5.0

    res = spline_vec.deriv_array(x)

    # Check against analytical derivatives
    assert abs(res[0] - df0(x)) < 1e-10
    assert abs(res[1] - df1(x)) < 1e-9
    assert abs(res[2] - df2(x)) < 1e-3


def test_spline_vec_integ(spline_vec: Ncm.SplineVec) -> None:
    """Test integration of NcmSplineVec."""
    xi = 2.0
    xf = 7.0

    res = spline_vec.integ_array(xi, xf)

    # Check that results match individual spline integrations
    for i in range(3):
        s_i = spline_vec.peek_spline(i)
        expected = s_i.eval_integ(xi, xf)
        computed = res[i]
        assert abs(computed - expected) < 1e-15


def test_spline_vec_integ_values(spline_vec: Ncm.SplineVec) -> None:
    """Test that integral values are close to expected integral values."""
    xi = 2.0
    xf = 7.0

    res = spline_vec.integ_array(xi, xf)

    # Check against analytical integrals
    assert abs(res[0] - int_f0(xi, xf)) < 1e-10
    assert abs(res[1] - int_f1(xi, xf)) < 1e-9
    assert abs(res[2] - int_f2(xi, xf)) < 1e-3


def test_spline_vec_peek_spline(spline_vec: Ncm.SplineVec) -> None:
    """Test accessing individual splines."""
    for i in range(3):
        s_i = spline_vec.peek_spline(i)
        assert isinstance(s_i, Ncm.Spline)
        assert s_i.is_init()


def test_spline_vec_both_constructors_equivalent() -> None:
    """Test that both constructors produce equivalent results."""
    nknots = 50
    xi = 0.0
    xf = 5.0

    # Create x vector
    x_arr = np.linspace(xi, xf, nknots)
    xv = Ncm.Vector.new_array(x_arr.tolist())

    # Matrix constructor
    ym = Ncm.Matrix.new(3, nknots)
    for i, x in enumerate(x_arr):
        ym.set(0, i, f0(x))
        ym.set(1, i, f1(x))
        ym.set(2, i, f2(x))

    s1 = Ncm.SplineCubicNotaknot.new()
    sv1 = Ncm.SplineVec.new(s1, xv, ym, True)

    # GPtrArray constructor
    yv_list = []
    y0 = Ncm.Vector.new(nknots)
    y1 = Ncm.Vector.new(nknots)
    y2 = Ncm.Vector.new(nknots)

    for i, x in enumerate(x_arr):
        y0.set(i, f0(x))
        y1.set(i, f1(x))
        y2.set(i, f2(x))

    yv_list.append(y0)
    yv_list.append(y1)
    yv_list.append(y2)

    s2 = Ncm.SplineCubicNotaknot.new()
    sv2 = Ncm.SplineVec.new_gpa(s2, xv, yv_list, True)

    # Compare evaluations
    x = 2.5

    res1 = sv1.eval_array(x)
    res2 = sv2.eval_array(x)

    for i in range(3):
        assert abs(res1[i] - res2[i]) < 1e-15


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
