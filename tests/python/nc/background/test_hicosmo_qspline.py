#!/usr/bin/env python
#
# test_hicosmo_qspline.py
#
# Sat Jun 21 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_hicosmo_qspline.py
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

"""Unit tests for NcHICosmoQSpline curvature functionals."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()


@pytest.fixture(name="cosmo_q")
def fixture_cosmo_q() -> Nc.HICosmoQSpline:
    """A QSpline model with a non-trivial q(z)."""
    n_knots = 10
    s = Ncm.SplineCubicNotaknot.new()
    cosmo = Nc.HICosmoQSpline.new(s, n_knots, 2.0)
    # Seeded, genuinely curved q(z) profile (non-zero second derivative).
    rng = np.random.default_rng(42)
    q_array = rng.uniform(-1.0, 0.5, size=n_knots)
    cosmo.props.qparam = Ncm.Vector.new_array(q_array.tolist())
    return cosmo


@pytest.mark.parametrize("name", ["lp_kappa", "lp_q2"])
def test_qspline_lp_mset_func(cosmo_q: Nc.HICosmoQSpline, name: str) -> None:
    """The lp_* MSetFuncList entries take p = x[0] and match the C accessor."""
    ctype = (
        Ncm.SplineCurvatureType.GEOMETRIC
        if name == "lp_kappa"
        else Ncm.SplineCurvatureType.D2
    )
    mset = Ncm.MSet.new_array([cosmo_q])
    func = Ncm.MSetFuncList.new(f"NcHICosmoQSpline:{name}", None)

    assert func.get_nvar() == 1
    for p in (2.0, 8.0, 16.0):
        n_p = func.eval1(mset, p)
        assert np.isfinite(n_p)
        assert n_p >= 0.0
        assert_allclose(n_p, cosmo_q.lp_norm(ctype, p))


def test_qspline_lp_norm_nondecreasing(cosmo_q: Nc.HICosmoQSpline) -> None:
    """L_p norm is non-decreasing in p for both curvature types."""
    for ctype in (Ncm.SplineCurvatureType.D2, Ncm.SplineCurvatureType.GEOMETRIC):
        norms = [cosmo_q.lp_norm(ctype, p) for p in (2.0, 8.0, 32.0)]
        assert norms[0] <= norms[1] * (1.0 + 1e-9)
        assert norms[1] <= norms[2] * (1.0 + 1e-9)


def test_qspline_mean_kappa(cosmo_q: Nc.HICosmoQSpline) -> None:
    """mean_kappa is the p=2 geometric case of lp_norm, and matches its MSetFunc."""
    mean_kappa = cosmo_q.mean_kappa()
    assert_allclose(mean_kappa, cosmo_q.lp_norm(Ncm.SplineCurvatureType.GEOMETRIC, 2.0))

    func = Ncm.MSetFuncList.new("NcHICosmoQSpline:mean_kappa", None)
    assert func.get_nvar() == 0
    mset = Ncm.MSet.new_array([cosmo_q])
    assert_allclose(mean_kappa, func.eval0(mset))


def test_qspline_lp_prior(cosmo_q: Nc.HICosmoQSpline) -> None:
    """A curvature L_p norm can drive a Gaussian prior via the var = p slot."""
    func = Ncm.MSetFuncList.new("NcHICosmoQSpline:lp_q2", None)
    prior = Ncm.PriorGaussFunc.new(func, 0.0, 1.0, 16.0)
    assert prior is not None
