#!/usr/bin/env python
#
# example_diff.py
#
# Fri May 19 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_diff.py
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

"""Example of how to use the numcosmo library to calculate derivatives
of functions. Also compares the results with the analytical derivatives
and with the derivatives calculated by the python library numdifftools."""

import math

import numpy as np
import numdifftools as nd

from numcosmo_py import Ncm


#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


#
# Functions to be differentiated
#
def ftest(x_v, y_v, *args):
    """Function to be differentiated."""
    p0 = args[0]
    roffpad = args[1]
    x_a = x_v.dup_array()

    y_v.set(0, math.sin(math.pi * x_a[0] * x_a[1]))
    y_v.set(1, roffpad + math.cos(math.pi * x_a[0] / x_a[1]))
    y_v.set(2, math.exp(x_a[0] * p0))
    y_v.set(3, (x_a[0] - 10.0) ** (3.0))


def ftest_py(x_a, *args):
    """Function to be differentiated."""
    p0 = args[0]
    roffpad = args[1]
    # print ((x_a - x0_a)/x0_a)
    return np.array(
        [
            math.sin(math.pi * x_a[0] * x_a[1]),
            roffpad + math.cos(math.pi * x_a[0] / x_a[1]),
            math.exp(x_a[0] * p0),
            (x_a[0] - 10.0) ** (3.0),
        ]
    )


def ftest2(x_v, *_):
    """Function to be differentiated."""
    x_a = x_v.dup_array()
    return math.sin(math.pi * (x_a[0] * x_a[1] * x_a[2] - 1.0))


#
# Comparison + log function
#
def cmp_array(a, b, err, title=None):
    """Compare two arrays and print the results."""
    a = np.array(a)
    b = np.array(b)

    if title:
        print(f"#\n# {title}\n#")
    print(
        " Numerical derivative   Analytical derivative  Error"
        "                  Estimate error"
    )
    for a_i, b_i, err_i in zip(a, b, err):
        eerr = 0.0
        aerr = 0.0
        if b_i == 0.0:
            aerr = math.fabs(a_i)
            eerr = err_i
        else:
            aerr = math.fabs(a_i / b_i - 1.0)
            eerr = math.fabs(err_i / b_i)
        print(
            f"{a_i: 22.15e} {b_i: 22.15e} {aerr: 22.15e} "
            f"{eerr: 22.15e} [{'PASS' if eerr >= aerr else 'FAIL'}]"
        )
    print("")


def test_diff() -> None:
    """Test the numcosmo library to calculate derivatives of functions."""

    #
    # New diff object
    #
    diff = Ncm.Diff.new()
    # diff.set_ini_h (2.0)
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
    df_a, err_a = diff.rf_d1_N_to_M(x0_a.tolist(), 4, ftest, L, roffpad)

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

    print(np.log1p(x0_a).clip(1.0) * 2.0)
    cmp_array(
        df_a,
        dfE_a,
        err_a,
        "First derivative: Forward method + Richardson extrapolation",
    )

    def flatten(l):
        """Flatten a list of lists."""
        return [item for sublist in l for item in sublist]

    jac = nd.Jacobian(ftest_py)
    cmp_array(flatten(np.transpose(jac(x0_a, L, roffpad))), dfE_a, err_a)

    #
    # First derivative: Central method + Richardson extrapolation
    #
    df_a, err_a = diff.rc_d1_N_to_M(x0_a.tolist(), 4, ftest, L, roffpad)

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

    cmp_array(
        df_a,
        dfE_a,
        err_a,
        "First derivative: Central method + Richardson extrapolation",
    )

    #
    # Second derivative: Central method + Richardson extrapolation
    #
    df_a, err_a = diff.rc_d2_N_to_M(x0_a.tolist(), 4, ftest, L, roffpad)

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

    cmp_array(
        df_a,
        dfE_a,
        err_a,
        "Second derivative: Central method + Richardson extrapolation",
    )

    #
    # Hessian matrix: Forward method + Richardson extrapolation
    #
    x0_a = np.array([1.5, 2.0, 3.0])
    H_a, Herr_a = diff.rf_Hessian_N_to_1(x0_a.tolist(), ftest2, None, None)

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

    cmp_array(
        H_a, HE_a, Herr_a, "Hessian matrix: Forward method + Richardson extrapolation"
    )


if __name__ == "__main__":
    test_diff()
