#!/usr/bin/env python
#
# test_py_diff.py
#
# Fri May 19 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_diff.py
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

"""Tests on NcmDiff class."""

import math
import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.helper import npa_to_seq

Ncm.cfg_init()


#
# Functions to be differentiated
#
def ftest(x_v, y_v, *args):
    """Compute function to be differentiated."""
    p0 = args[0]
    roffpad = args[1]
    x_a = x_v.dup_array()

    y_v.set(0, math.sin(math.pi * x_a[0] * x_a[1]))
    y_v.set(1, roffpad + math.cos(math.pi * x_a[0] / x_a[1]))
    y_v.set(2, math.exp(x_a[0] * p0))
    y_v.set(3, (x_a[0] - 10.0) ** (3.0))


def ftest_py(x_a, *args):
    """Compute function to be differentiated."""
    p0 = args[0]
    roffpad = args[1]
    return np.array(
        [
            math.sin(math.pi * x_a[0] * x_a[1]),
            roffpad + math.cos(math.pi * x_a[0] / x_a[1]),
            math.exp(x_a[0] * p0),
            (x_a[0] - 10.0) ** (3.0),
        ]
    )


def ftest2(x_v, *_):
    """Compute function to be differentiated."""
    x_a = x_v.dup_array()
    return math.sin(math.pi * (x_a[0] * x_a[1] * x_a[2] - 1.0))


#
# Comparison + log function
#
def cmp_array(a, b, err):
    """Compare two arrays and print the results."""
    a = np.array(a)
    b = np.array(b)

    for a_i, b_i, err_i in zip(a, b, err):
        eerr = 0.0
        aerr = 0.0
        if b_i == 0.0:
            aerr = math.fabs(a_i)
            eerr = err_i
        else:
            aerr = math.fabs(a_i / b_i - 1.0)
            eerr = math.fabs(err_i / b_i)

        assert eerr >= aerr


def test_diff() -> None:
    """Test the numcosmo library to calculate derivatives of functions."""
    diff = Ncm.Diff.new()
    _run_diff_checks(diff)


def test_diff_dual_series() -> None:
    """Test derivatives with the dual-series scheme enabled."""
    diff = Ncm.Diff.new()
    diff.set_dual_series(True)
    assert diff.get_dual_series()
    _run_diff_checks(diff)


def test_diff_dual_series_property() -> None:
    """Test the dual-series property default and round trip."""
    diff = Ncm.Diff.new()
    assert not diff.get_dual_series()

    diff_dual = Ncm.Diff(**{"dual-series": True})
    assert diff_dual.get_dual_series()
    assert diff_dual.props.dual_series

    diff_dual.set_dual_series(False)
    assert not diff_dual.get_dual_series()


def test_diff_dual_series_1_to_1() -> None:
    """Test the scalar entry points with the dual-series scheme enabled."""
    diff = Ncm.Diff(**{"dual-series": True})
    x0 = 1.3

    for method, exact in (
        (diff.rc_d1_1_to_1, math.cos(x0)),
        (diff.rf_d1_1_to_1, math.cos(x0)),
        (diff.rc_d2_1_to_1, -math.sin(x0)),
    ):
        val, err = method(x0, lambda x, *_: math.sin(x), None)
        assert math.fabs(val - exact) <= err


def test_diff_spectral_window_property() -> None:
    """Test the spectral-window property default and round trip."""
    diff = Ncm.Diff.new()
    assert diff.get_spectral_window() == 1.0

    diff.set_spectral_window(0.5)
    assert diff.get_spectral_window() == 0.5

    diff2 = Ncm.Diff(**{"spectral-window": 2.0})
    assert diff2.props.spectral_window == 2.0


def test_diff_spectral_1_to_1() -> None:
    """Test the spectral scalar derivatives against analytical results."""
    diff = Ncm.Diff.new()

    for x0, f, d1, d2 in (
        (1.3, math.sin, math.cos(1.3), -math.sin(1.3)),
        (1.0, math.exp, math.e, math.e),
        (0.0, math.sin, 1.0, 0.0),
        (2.0, lambda x: x**3, 12.0, 12.0),
    ):
        val, err = diff.sc_d1_1_to_1(x0, lambda x, *_: f(x), None)
        assert math.fabs(val - d1) <= err

        val, err = diff.sc_d2_1_to_1(x0, lambda x, *_: f(x), None)
        assert math.fabs(val - d2) <= err


def test_diff_spectral_scale_detection() -> None:
    """Test that the spectral window search finds narrow features."""
    diff = Ncm.Diff.new()

    # Varies on scale 1e-3 around x0 = 1; the initial window is 1.
    val, err = diff.sc_d1_1_to_1(1.0, lambda x, *_: math.atan(1000.0 * (x - 1.0)), None)
    assert math.fabs(val - 1000.0) <= err
    assert math.fabs(val / 1000.0 - 1.0) < 1.0e-6


def test_diff_spectral_offset_cancellation() -> None:
    """Test the spectral derivative of a function with a large constant offset."""
    diff = Ncm.Diff.new()

    val, err = diff.sc_d1_1_to_1(1.3, lambda x, *_: 1.0e10 + math.sin(x), None)
    assert math.fabs(val - math.cos(1.3)) <= err
    assert math.fabs(val / math.cos(1.3) - 1.0) < 1.0e-3


def test_diff_spectral_N_to_M() -> None:
    """Test the spectral vector derivatives against analytical results."""
    diff = Ncm.Diff.new()
    x0_a = np.array([3.0, 1.01234])

    def ftest_sc(x_v, y_v, *_):
        x_a = x_v.dup_array()
        y_v.set(0, math.sin(math.pi * x_a[0] * x_a[1]))
        y_v.set(1, math.exp(x_a[0] * 1.0e-2))
        y_v.set(2, (x_a[0] - 10.0) ** 3)

    df_a, err_a = diff.sc_d1_N_to_M(npa_to_seq(x0_a), 3, ftest_sc, None)
    dfE_a = [
        math.pi * x0_a[1] * math.cos(math.pi * x0_a[0] * x0_a[1]),
        1.0e-2 * math.exp(x0_a[0] * 1.0e-2),
        3.0 * (x0_a[0] - 10.0) ** 2,
        math.pi * x0_a[0] * math.cos(math.pi * x0_a[0] * x0_a[1]),
        0.0,
        0.0,
    ]
    cmp_array(df_a, dfE_a, err_a)

    df_a, err_a = diff.sc_d2_N_to_M(npa_to_seq(x0_a), 3, ftest_sc, None)
    dfE_a = [
        -math.pi**2 * x0_a[1] ** 2 * math.sin(math.pi * x0_a[0] * x0_a[1]),
        1.0e-4 * math.exp(x0_a[0] * 1.0e-2),
        6.0 * (x0_a[0] - 10.0),
        -math.pi**2 * x0_a[0] ** 2 * math.sin(math.pi * x0_a[0] * x0_a[1]),
        0.0,
        0.0,
    ]
    cmp_array(df_a, dfE_a, err_a)


def test_diff_spectral_N_to_1_and_1_to_M() -> None:
    """Test the spectral N-to-1 and 1-to-M wrappers."""
    diff = Ncm.Diff.new()
    x0_a = np.array([1.5, 2.0])

    def fscal(x_v, *_):
        x_a = x_v.dup_array()
        return math.sin(x_a[0]) * math.exp(0.1 * x_a[1])

    df_a, err_a = diff.sc_d1_N_to_1(npa_to_seq(x0_a), fscal, None)
    dfE_a = [
        math.cos(x0_a[0]) * math.exp(0.1 * x0_a[1]),
        0.1 * math.sin(x0_a[0]) * math.exp(0.1 * x0_a[1]),
    ]
    cmp_array(df_a, dfE_a, err_a)

    df_a, err_a = diff.sc_d2_N_to_1(npa_to_seq(x0_a), fscal, None)
    dfE_a = [
        -math.sin(x0_a[0]) * math.exp(0.1 * x0_a[1]),
        0.01 * math.sin(x0_a[0]) * math.exp(0.1 * x0_a[1]),
    ]
    cmp_array(df_a, dfE_a, err_a)

    def fvec(x, y_v, *_):
        y_v.set(0, math.sin(x))
        y_v.set(1, x**2)

    df_a, err_a = diff.sc_d1_1_to_M(1.3, 2, fvec, None)
    cmp_array(df_a, [math.cos(1.3), 2.6], err_a)

    df_a, err_a = diff.sc_d2_1_to_M(1.3, 2, fvec, None)
    cmp_array(df_a, [-math.sin(1.3), 2.0], err_a)


def _run_diff_checks(diff: Ncm.Diff) -> None:
    """Check all derivative methods of diff against analytical results."""

    roffpad = 0.0
    x0_a = np.array([3.0, 1.01234])

    #
    # Point where to calculate the derivative + function parameters
    #
    # x0_a = [3.12345, 1.012345]
    L = 1.0e-2

    #
    # First derivative: Forward method + Richardson extrapolation
    #
    df_a, err_a = diff.rf_d1_N_to_M(npa_to_seq(x0_a), 4, ftest, L, roffpad)
    # Analytical derivative
    dfE_a = [
        math.pi * x0_a[1] * math.cos(math.pi * x0_a[0] * x0_a[1]),
        -math.pi / x0_a[1] * math.sin(math.pi * x0_a[0] / x0_a[1]),
        L * math.exp(x0_a[0] * L),
        3.0 * (x0_a[0] - 10.0) ** 2,
        math.pi * x0_a[0] * math.cos(math.pi * x0_a[0] * x0_a[1]),
        +math.pi * x0_a[0] / x0_a[1] ** 2 * math.sin(math.pi * x0_a[0] / x0_a[1]),
        0.0,
        0.0,
    ]
    cmp_array(df_a, dfE_a, err_a)

    #
    # First derivative: Central method + Richardson extrapolation
    #
    df_a, err_a = diff.rc_d1_N_to_M(npa_to_seq(x0_a), 4, ftest, L, roffpad)
    # Analytical derivative
    dfE_a = [
        math.pi * x0_a[1] * math.cos(math.pi * x0_a[0] * x0_a[1]),
        -math.pi / x0_a[1] * math.sin(math.pi * x0_a[0] / x0_a[1]),
        L * math.exp(x0_a[0] * L),
        3.0 * (x0_a[0] - 10.0) ** 2,
        math.pi * x0_a[0] * math.cos(math.pi * x0_a[0] * x0_a[1]),
        +math.pi * x0_a[0] / x0_a[1] ** 2 * math.sin(math.pi * x0_a[0] / x0_a[1]),
        0.0,
        0.0,
    ]
    cmp_array(df_a, dfE_a, err_a)

    #
    # Second derivative: Central method + Richardson extrapolation
    #
    df_a, err_a = diff.rc_d2_N_to_M(npa_to_seq(x0_a), 4, ftest, L, roffpad)
    # Analytical derivative
    dfE_a = [
        -math.pi**2 * x0_a[1] ** 2 * math.sin(math.pi * x0_a[0] * x0_a[1]),
        -math.pi**2 / x0_a[1] ** 2 * math.cos(math.pi * x0_a[0] / x0_a[1]),
        (L) ** 2 * math.exp(x0_a[0] * L),
        6.0 * (x0_a[0] - 10.0),
        -math.pi**2 * x0_a[0] ** 2 * math.sin(math.pi * x0_a[0] * x0_a[1]),
        -2.0 * math.pi * x0_a[0] / x0_a[1] ** 3 * math.sin(math.pi * x0_a[0] / x0_a[1])
        - (math.pi * x0_a[0] / x0_a[1] ** 2) ** 2
        * math.cos(math.pi * x0_a[0] / x0_a[1]),
        0.0,
        0.0,
    ]
    cmp_array(df_a, dfE_a, err_a)

    #
    # Hessian matrix: Forward method + Richardson extrapolation
    #
    x0_a = np.array([1.5, 2.0, 3.0])
    H_a, Herr_a = diff.rf_Hessian_N_to_1(npa_to_seq(x0_a), ftest2, None, None)
    # Analytical derivative
    HE_a = [
        -((math.pi * x0_a[1] * x0_a[2]) ** 2)
        * math.sin(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0)),
        math.pi * x0_a[2] * math.cos(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0))
        - math.pi**2
        * x0_a[0]
        * x0_a[1]
        * x0_a[2] ** 2
        * math.sin(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0)),
        math.pi * x0_a[1] * math.cos(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0))
        - math.pi**2
        * x0_a[0]
        * x0_a[1] ** 2
        * x0_a[2]
        * math.sin(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0)),
        math.pi * x0_a[2] * math.cos(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0))
        - math.pi**2
        * x0_a[0]
        * x0_a[1]
        * x0_a[2] ** 2
        * math.sin(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0)),
        -((math.pi * x0_a[0] * x0_a[2]) ** 2)
        * math.sin(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0)),
        math.pi * x0_a[0] * math.cos(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0))
        - math.pi**2
        * x0_a[0] ** 2
        * x0_a[1]
        * x0_a[2]
        * math.sin(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0)),
        math.pi * x0_a[1] * math.cos(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0))
        - math.pi**2
        * x0_a[0]
        * x0_a[1] ** 2
        * x0_a[2]
        * math.sin(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0)),
        math.pi * x0_a[0] * math.cos(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0))
        - math.pi**2
        * x0_a[0] ** 2
        * x0_a[1]
        * x0_a[2]
        * math.sin(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0)),
        -((math.pi * x0_a[0] * x0_a[1]) ** 2)
        * math.sin(math.pi * (x0_a[0] * x0_a[1] * x0_a[2] - 1.0)),
    ]
    cmp_array(H_a, HE_a, Herr_a)
