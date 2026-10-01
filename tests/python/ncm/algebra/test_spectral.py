#
# test_py_spectral.py
#
# Tue Feb 04 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_spectral.py
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

"""Tests of the NcmSpectral bindings: Python callbacks and returned arrays.

The numerics are tested in tests/c/ncm/algebra/test_ncm_spectral.c.
"""

import math

from numcosmo_py import Ncm

Ncm.cfg_init()


def _exp(_user_data, x: float) -> float:
    return math.exp(x)


def test_coefficients_from_python_callbacks() -> None:
    """Each expansion takes a Python callback and returns a new array."""
    spectral = Ncm.Spectral.new()

    coeffs = spectral.compute_chebyshev_coeffs(_exp, -1.0, 1.0, 24, None)
    assert len(coeffs) == 24

    k, coeffs = spectral.compute_chebyshev_coeffs_adaptive(
        _exp, -1.0, 1.0, 2, 1.0e-12, None
    )
    assert len(coeffs) == (1 << k) + 1

    k, coeffs = spectral.compute_chebyshev_coeffs_adaptive_full(
        _exp, -1.0, 1.0, 2, 1.0e-12, 0.0, None
    )
    assert len(coeffs) == (1 << k) + 1

    k, coeffs, converged = spectral.compute_chebyshev_coeffs_adaptive_try(
        _exp, -1.0, 1.0, 2, 10, 1.0e-12, 0.0, None
    )
    assert converged
    assert len(coeffs) == (1 << k) + 1

    integral = Ncm.Spectral.chebyshev_integrate(coeffs, -1.0, 1.0)
    assert math.isclose(integral, math.e - 1.0 / math.e, rel_tol=1.0e-15)


def test_batch_from_python_callback() -> None:
    """The batch callback fills the NcmVector it is given."""
    spectral = Ncm.Spectral.new()

    def batch(_user_data, x: float, y: Ncm.Vector) -> None:
        y.set(0, math.exp(x))
        y.set(1, math.cos(x))

    k, coeffs = spectral.compute_chebyshev_coeffs_batch_adaptive(
        batch, 2, -1.0, 1.0, 2, 1.0e-12, 0.0, None
    )
    assert coeffs.nrows() == 2
    assert coeffs.ncols() == (1 << k) + 1

    k, coeffs = spectral.compute_chebyshev_coeffs_batch_adaptive_cap(
        batch, 2, -1.0, 1.0, 2, 10, 1.0e-12, 0.0, False, None
    )
    assert k > 0
    assert coeffs.ncols() == (1 << k) + 1


def test_arrays_returned_by_conversions() -> None:
    """Conversions, rebasing and matrices return new objects."""
    spectral = Ncm.Spectral.new()
    c = [1.0, 0.5, 0.25, 0.125]

    assert len(Ncm.Spectral.chebT_to_gegenbauer_alpha1(c)) == 4
    assert len(Ncm.Spectral.chebT_to_gegenbauer_alpha2(c)) == 4
    assert len(Ncm.Spectral.chebT_deriv_to_gegenbauer_alpha2(c)) == 3
    assert len(Ncm.Spectral.chebT_deriv2_to_gegenbauer_alpha2(c)) == 2
    assert len(Ncm.Spectral.gegenbauer_alpha2_xmul(c, 1.0, 0.0)) == 5

    norm, rebased = spectral.chebyshev_rebase(c, 0, -1.0, 1.0, 0.0, 1.0)
    assert len(rebased) == 4
    assert norm >= abs(Ncm.Spectral.chebyshev_eval(c, 0.7))

    assert Ncm.Spectral.get_d2_matrix(8).nrows() == 8
