#!/usr/bin/env python
#
# example_simple.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_simple.py
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

"""Simple example computing cosmological distances."""

from typing import Tuple

import numpy as np
import timeit

from numcosmo_py import Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


class TestIntegralND(Ncm.Integralnd):
    """Test class for IntegralND."""

    def do_get_dimensions(self) -> Tuple[int, int]:
        """Get number of dimensions."""
        return 3, 3

    def do_integrand(
        self,
        x_vec: Ncm.Vector,
        dim: int,
        npoints: int,
        fdim: int,
        fval_vec: Ncm.Vector,
    ) -> None:
        """Integrand function."""

        x = np.array(x_vec.dup_array())
        fval = (3.0 * x**2 + 2.0 * x) * np.sin(1.0 * x)

        fval_vec.set_array(fval.tolist())


def test_integralnd() -> None:
    """Example computing cosmological distances."""
    test_f = TestIntegralND(method=Ncm.IntegralndMethod.H_V)

    res = Ncm.Vector.new(3)
    err = Ncm.Vector.new(3)

    def test_int():
        test_f.eval(
            Ncm.Vector.new_array([0.0, 0.0, 0.0]),
            Ncm.Vector.new_array([1.0, 1.0, 1.0]),
            res,
            err,
        )

    test_f.set_method(Ncm.IntegralndMethod.P_V)
    print(timeit.timeit(test_int, number=10))

    test_f.set_method(Ncm.IntegralndMethod.H_V)
    print(timeit.timeit(test_int, number=10))

    test_f.set_method(Ncm.IntegralndMethod.P)
    print(timeit.timeit(test_int, number=10))

    test_f.set_method(Ncm.IntegralndMethod.H)
    print(timeit.timeit(test_int, number=10))

    print(res.dup_array(), err.dup_array())


if __name__ == "__main__":
    test_integralnd()
