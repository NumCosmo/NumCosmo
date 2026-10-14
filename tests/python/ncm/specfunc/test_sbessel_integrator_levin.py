#!/usr/bin/env python
#
# test_py_sbessel_integrator_levin.py
#
# Sat Jan 25 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sbessel_integrator_levin.py
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

"""Unit tests for Levin spherical Bessel integrator."""

import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy.special import spherical_jn
from scipy.integrate import quad

from numcosmo_py import Ncm


def _jl_second_deriv(ell: int, x: np.ndarray) -> np.ndarray:
    """j_l'' from the homogeneous ODE and scipy's j_l, j_l'."""
    j = spherical_jn(ell, x)
    jp = spherical_jn(ell, x, derivative=True)
    return -2.0 / x * jp + (ell * (ell + 1.0) / x**2 - 1.0) * j


class TestSBesselIntegratorLevinDeriv:
    """Tests for the derivative-weighted integrals of NcmSBesselIntegratorLevin.

    integrate_deriv computes int_a^b K(x, k) j_l^{(d)}(k x) dx with the
    derivative taken with respect to the Bessel argument x = k chi.
    """

    @pytest.mark.parametrize("l_val", [0, 2, 10, 30])
    @pytest.mark.parametrize("deriv", [1, 2])
    def test_gaussian_vs_quadrature(self, l_val: int, deriv: int) -> None:
        """Gaussian K against scipy adaptive quadrature, including the turning point."""
        a, b = 1.0, 8.0
        center, sigma = 4.0, 1.2
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        result = Ncm.Vector.new(1)

        def f_gauss(x: float, _k: float) -> float:
            return np.exp(-0.5 * ((x - center) / sigma) ** 2)

        # k values placing the turning point x ~ l inside, below and above the window
        k_list = [0.5, 2.0]
        if l_val > 0:
            k_list += [l_val / center, l_val / a]

        for k in k_list:
            integrator.integrate_deriv(f_gauss, a, b, k, deriv, result)

            if deriv == 1:

                def integrand(x: float) -> float:
                    return f_gauss(x, k) * spherical_jn(l_val, k * x, derivative=True)

            else:

                def integrand(x: float) -> float:
                    return f_gauss(x, k) * _jl_second_deriv(l_val, np.asarray(k * x))

            expected, _ = quad(integrand, a, b, limit=800, epsabs=1e-14, epsrel=1e-13)
            assert_allclose(
                result.get(0),
                expected,
                rtol=1.0e-8,
                atol=1.0e-13,
                err_msg=f"l={l_val}, deriv={deriv}, k={k}",
            )
