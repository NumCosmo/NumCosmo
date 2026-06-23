#!/usr/bin/env python
#
# test_spline_curvature.py
#
# Sat Jun 21 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_spline_curvature.py
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

"""Tests on the generic NcmSpline curvature L_p functionals."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm

Ncm.cfg_init()

XI, XF = -1.0, 1.0


@pytest.fixture(name="spline_x2")
def fixture_spline_x2() -> Ncm.Spline:
    """A not-a-knot cubic spline of f(x) = x^2 (reproduced exactly)."""
    x = np.linspace(XI, XF, 41)
    y = x**2
    xv = Ncm.Vector.new_array(x.tolist())
    yv = Ncm.Vector.new_array(y.tolist())
    return Ncm.SplineCubicNotaknot.new_full(xv, yv, True)


def _lp_reference(c, p: float) -> float:
    """Domain-normalized L_p norm of curvature density on a fine grid."""
    x = np.linspace(XI, XF, 200001)
    return (np.trapezoid(np.abs(c(x)) ** p, x) / (XF - XI)) ** (1.0 / p)


@pytest.mark.parametrize("p", [2.0, 4.0, 8.0, 16.0])
def test_d2_constant(spline_x2: Ncm.Spline, p: float) -> None:
    """For f = x^2, f'' = 2 everywhere, so the D2 L_p norm is 2 for all p."""
    n_p = spline_x2.curvature_lp_norm(Ncm.SplineCurvatureType.D2, p, XI, XF)
    assert_allclose(n_p, 2.0, rtol=1e-7)


def test_d2_max(spline_x2: Ncm.Spline) -> None:
    """For f = x^2, the maximum |f''| is 2."""
    cmax = spline_x2.curvature_max(Ncm.SplineCurvatureType.D2, XI, XF)
    assert_allclose(cmax, 2.0, rtol=1e-7)


@pytest.mark.parametrize("p", [2.0, 4.0, 8.0, 16.0])
def test_geometric_matches_reference(spline_x2: Ncm.Spline, p: float) -> None:
    """Geometric L_p norm matches a fine-grid numerical reference."""

    def kappa(x):
        return 2.0 / (1.0 + 4.0 * x**2) ** 1.5

    n_p = spline_x2.curvature_lp_norm(Ncm.SplineCurvatureType.GEOMETRIC, p, XI, XF)
    assert_allclose(n_p, _lp_reference(kappa, p), rtol=1e-4)


def test_geometric_max(spline_x2: Ncm.Spline) -> None:
    """Geometric curvature of x^2 peaks at x=0 with value 2."""
    cmax = spline_x2.curvature_max(Ncm.SplineCurvatureType.GEOMETRIC, XI, XF)
    assert_allclose(cmax, 2.0, rtol=1e-4)


def test_lp_converges_to_max(spline_x2: Ncm.Spline) -> None:
    """L_p geometric norm increases monotonically toward the max as p grows.

    Convergence is slow because the geometric curvature of x^2 is sharply peaked
    at x=0, so we only require a monotone increase that stays a lower bound of
    the max and is within ~10% by p=32.
    """
    cmax = spline_x2.curvature_max(Ncm.SplineCurvatureType.GEOMETRIC, XI, XF)
    norms = [
        spline_x2.curvature_lp_norm(Ncm.SplineCurvatureType.GEOMETRIC, p, XI, XF)
        for p in (2.0, 8.0, 32.0)
    ]
    assert norms[0] < norms[1] < norms[2] <= cmax * (1.0 + 1e-9)
    assert_allclose(norms[2], cmax, rtol=1e-1)


def test_density_pointwise(spline_x2: Ncm.Spline) -> None:
    """Pointwise curvature density matches the analytic expressions."""
    for x in (-0.7, -0.2, 0.0, 0.3, 0.9):
        d2 = spline_x2.curvature_density(Ncm.SplineCurvatureType.D2, x)
        geo = spline_x2.curvature_density(Ncm.SplineCurvatureType.GEOMETRIC, x)
        assert_allclose(d2, 2.0, rtol=1e-7)
        assert_allclose(geo, 2.0 / (1.0 + 4.0 * x**2) ** 1.5, rtol=1e-6)
