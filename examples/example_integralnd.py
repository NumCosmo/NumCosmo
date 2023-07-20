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

from numcosmo_py import Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


class TestIntegralND(Ncm.Integralnd):
    """Test class for IntegralND."""

    def do_integrand(self, x: list[float], dim, npoints, fdim) -> list[float]:
        """Integrand function."""
        return [1.0]


def test_integralnd() -> None:
    """Example computing cosmological distances."""
    test_f = TestIntegralND()

    test_f.eval([0.0, 0.0], [1.0, 1.0])


if __name__ == "__main__":
    test_integralnd()
