#!/usr/bin/env python
#
# test_py_hicosmo_de_wspline.py
#
# Mon Sep 22 18:14:50 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_hicosmo_de_wspline.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Unit tests for the NcmMSet class."""

import pytest
import numpy as np
from numpy.testing import assert_allclose
from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="cosmo_w", params=[3, 4, 6, 12])
def fixture_cosmo_w(request) -> Nc.HICosmoDEWSpline:
    """Returns a DEWSpline model with preset knots and values."""
    return Nc.HICosmoDEWSpline.new(nknots=request.param, z_f=2.0)


def test_wspline_knots(cosmo_w: Nc.HICosmoDEWSpline) -> None:
    """Test that knots are set and retrieved correctly."""

    alpha = cosmo_w.get_alpha()

    assert alpha.len() > 0

    alpha_array = alpha.dup_array()
    # Now we check if it is monotonically increasing
    assert np.all(np.diff(alpha_array) > 0)


def test_wspline_eval_w(cosmo_w: Nc.HICosmoDEWSpline) -> None:
    """Test evaluation of w(z) at knot points and in between."""
    alpha = cosmo_w.get_alpha()
    w_len = alpha.len()
    w_array = np.random.uniform(-1.5, -0.5, size=w_len).astype(float)
    z_test = np.expm1(alpha.dup_array())

    for z in z_test:
        w_test = cosmo_w.w_de(z)
        assert_allclose(w_test, -1.0)  # Default value is -1.0

    w_vec = Ncm.Vector.new_array(w_array)
    cosmo_w.props.w = w_vec

    for z, w in zip(z_test, w_array):
        w_test = cosmo_w.w_de(z)
        assert_allclose(w_test, w)


def test_wspline_mean_kappa(cosmo_w: Nc.HICosmoDEWSpline) -> None:
    """Test evaluation of mean kappa at knot points and in between."""
    alpha = cosmo_w.get_alpha()
    w_len = alpha.len()
    w_array = np.random.uniform(-1.5, -0.5, size=w_len).astype(float)

    kappa_mean = cosmo_w.mean_kappa()
    assert_allclose(kappa_mean, 0.0)  # Default value is 0.0

    w_vec = Ncm.Vector.new_array(w_array)
    cosmo_w.props.w = w_vec

    kappa_mean = cosmo_w.mean_kappa()
    assert isinstance(kappa_mean, float)
    assert not np.isnan(kappa_mean)
    assert not np.isinf(kappa_mean)
    assert kappa_mean != 0.0  # Should be different from default now


def test_wspline_mean_kappa_mset_func(cosmo_w: Nc.HICosmoDEWSpline) -> None:
    """Test evaluation of mean kappa at knot points and in between."""
    alpha = cosmo_w.get_alpha()
    w_len = alpha.len()
    w_array = np.random.uniform(-1.5, -0.5, size=w_len).astype(float)
    w_vec = Ncm.Vector.new_array(w_array)
    cosmo_w.props.w = w_vec

    kappa_mean = cosmo_w.mean_kappa()
    kappa_mean_func = Ncm.MSetFuncList.new("NcHICosmoDEWSpline:mean_kappa", None)

    mset = Ncm.MSet.new_array([cosmo_w])
    assert_allclose(kappa_mean, kappa_mean_func.eval0(mset))


def test_wspline_lp_norm_equals_mean_kappa(cosmo_w: Nc.HICosmoDEWSpline) -> None:
    """mean_kappa is the p=2 geometric case of lp_norm."""
    w_len = cosmo_w.get_alpha().len()
    cosmo_w.props.w = Ncm.Vector.new_array(
        np.random.uniform(-1.5, -0.5, size=w_len).tolist()
    )
    assert_allclose(
        cosmo_w.mean_kappa(),
        cosmo_w.lp_norm(Ncm.SplineCurvatureType.GEOMETRIC, 2.0),
    )


@pytest.mark.parametrize("name", ["lp_kappa", "lp_w2"])
def test_wspline_lp_mset_func(cosmo_w: Nc.HICosmoDEWSpline, name: str) -> None:
    """The lp_* MSetFuncList entries take p = x[0] and match the C accessor."""
    rng = np.random.default_rng(42)
    w_len = cosmo_w.get_alpha().len()
    cosmo_w.props.w = Ncm.Vector.new_array(
        rng.uniform(-1.5, -0.5, size=w_len).tolist()
    )
    ctype = (
        Ncm.SplineCurvatureType.GEOMETRIC
        if name == "lp_kappa"
        else Ncm.SplineCurvatureType.D2
    )
    mset = Ncm.MSet.new_array([cosmo_w])
    func = Ncm.MSetFuncList.new(f"NcHICosmoDEWSpline:{name}", None)

    assert func.get_nvar() == 1
    for p in (2.0, 8.0, 16.0):
        assert_allclose(func.eval1(mset, p), cosmo_w.lp_norm(ctype, p))

    # L_p norm is non-decreasing in p. The check stays at moderate p where the
    # |c|^p quadrature is reliable -- at very large p the integrand becomes too
    # peaked for the adaptive quadrature to resolve and can underestimate the
    # norm. The slack covers the near-constant curvature degenerate case (e.g.
    # few knots) where successive norms coincide exactly.
    norms = [func.eval1(mset, p) for p in (2.0, 8.0, 16.0)]
    assert norms[0] <= norms[1] * (1.0 + 1e-9)
    assert norms[1] <= norms[2] * (1.0 + 1e-9)


def test_wspline_knots_default_is_chebyshev() -> None:
    """The default knot placement is Chebyshev (clustered toward the endpoints)."""
    cosmo = Nc.HICosmoDEWSpline.new(nknots=8, z_f=2.0)
    assert cosmo.props.knots == Nc.HICosmoSplineKnots.CHEBYSHEV

    alpha = np.array(cosmo.get_alpha().dup_array())
    diffs = np.diff(alpha)
    # Chebyshev clusters at the ends: the central gap exceeds the boundary gaps.
    assert diffs[len(diffs) // 2] > diffs[0]
    assert diffs[len(diffs) // 2] > diffs[-1]


def test_wspline_knots_uniform_is_evenly_spaced() -> None:
    """The uniform option spaces knots evenly in alpha and spreads the high-z end."""
    n = 8
    cheb = Nc.HICosmoDEWSpline.new(nknots=n, z_f=2.0)
    uni = Nc.HICosmoDEWSpline(
        zf=2.0, w_length=n, knots=Nc.HICosmoSplineKnots.UNIFORM
    )
    assert uni.props.knots.value_nick == "uniform"

    a_cheb = np.array(cheb.get_alpha().dup_array())
    a_uni = np.array(uni.get_alpha().dup_array())

    # Uniform: constant spacing in alpha.
    assert_allclose(np.diff(a_uni), np.diff(a_uni)[0])
    # Both span the same range; uniform's last gap is wider than Chebyshev's.
    assert_allclose(a_uni[0], a_cheb[0])
    assert_allclose(a_uni[-1], a_cheb[-1])
    assert (a_uni[-1] - a_uni[-2]) > (a_cheb[-1] - a_cheb[-2])


def test_wspline_knots_survive_serialization() -> None:
    """The knot-placement choice round-trips through serialization."""
    uni = Nc.HICosmoDEWSpline(
        zf=2.0, w_length=8, knots=Nc.HICosmoSplineKnots.UNIFORM
    )
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dup = ser.dup_obj(uni)
    assert dup.props.knots.value_nick == "uniform"
    assert_allclose(dup.get_alpha().dup_array(), uni.get_alpha().dup_array())
