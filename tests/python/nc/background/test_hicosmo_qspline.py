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

from pathlib import Path

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
    """L_p norm is non-decreasing in p for both curvature types.

    Kept at moderate p where the |c|^p quadrature is reliable; at very large p
    the integrand is too peaked for the adaptive quadrature to resolve.
    """
    for ctype in (Ncm.SplineCurvatureType.D2, Ncm.SplineCurvatureType.GEOMETRIC):
        norms = [cosmo_q.lp_norm(ctype, p) for p in (2.0, 8.0, 16.0)]
        assert norms[0] <= norms[1] * (1.0 + 1e-9)
        assert norms[1] <= norms[2] * (1.0 + 1e-9)


def _z_spline(z_max: float, values) -> Ncm.Spline:
    """A not-a-knot weight spline W(z) on a uniform grid over [0, z_max]."""
    z = np.linspace(0.0, z_max, len(values))
    return Ncm.Spline.new(
        Ncm.SplineCubicNotaknot.new(),
        Ncm.Vector.new_array(z.tolist()),
        Ncm.Vector.new_array(np.asarray(values, dtype=float).tolist()),
        True,
    )


@pytest.mark.parametrize("ctype_name", ["GEOMETRIC", "D2"])
def test_qspline_weighted_lp_constant_weight(
    cosmo_q: Nc.HICosmoQSpline, ctype_name: str
) -> None:
    """A constant weight reduces the weighted L_p norm to the plain one."""
    ctype = getattr(Ncm.SplineCurvatureType, ctype_name)
    weight = _z_spline(2.0, [1.0] * 16)
    for p in (2.0, 4.0, 8.0):
        assert_allclose(
            cosmo_q.weighted_lp_norm(ctype, p, weight),
            cosmo_q.lp_norm(ctype, p),
            rtol=1e-6,
        )


def test_qspline_wlp_mset_func_serializes(cosmo_q: Nc.HICosmoQSpline) -> None:
    """wlp_* carries its weight spline through serialization and stays evaluable."""
    rng = np.random.default_rng(13)
    weight = _z_spline(2.0, rng.uniform(0.2, 2.0, 16))
    func = Ncm.MSetFuncList.new("NcHICosmoQSpline:wlp_q2", weight)
    assert func.get_nvar() == 1
    mset = Ncm.MSet.new_array([cosmo_q])
    direct = func.eval1(mset, 4.0)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    reloaded = ser.from_yaml(ser.to_yaml(func))
    assert_allclose(reloaded.eval1(mset, 4.0), direct)


def test_qspline_wlp_kappa_mset_func_evaluates(cosmo_q: Nc.HICosmoQSpline) -> None:
    """The geometric-curvature wlp_kappa accessor evaluates against a weight spline."""
    weight = _z_spline(2.0, [1.0] * 16)
    func = Ncm.MSetFuncList.new("NcHICosmoQSpline:wlp_kappa", weight)
    assert func.get_nvar() == 1
    mset = Ncm.MSet.new_array([cosmo_q])
    value = func.eval1(mset, 2.0)
    assert_allclose(
        value, cosmo_q.weighted_lp_norm(Ncm.SplineCurvatureType.GEOMETRIC, 2.0, weight)
    )


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


def test_qspline_q_transition() -> None:
    """The q=0 crossing matches the analytic root of a linear q(z)."""
    n_knots, z_max = 8, 2.0
    cosmo = Nc.HICosmoQSpline.new(Ncm.SplineCubicNotaknot.new(), n_knots, z_max)
    z_knots = np.array(cosmo.props.spline.peek_xv().dup_array())
    # q(z) = 0.5 - z crosses zero at z = 0.5 (exactly, linear -> cubic spline).
    cosmo.props.qparam = Ncm.Vector.new_array((0.5 - z_knots).tolist())

    assert_allclose(cosmo.q_transition(), 0.5, atol=1e-6)
    # And it matches the registered MSetFunc.
    func = Ncm.MSetFuncList.new("NcHICosmoQSpline:q_transition", None)
    mset = Ncm.MSet.new_array([cosmo])
    assert_allclose(cosmo.q_transition(), func.eval0(mset))


def test_qspline_q_transition_no_crossing() -> None:
    """With q(z) > 0 everywhere there is no transition, returning NaN."""
    cosmo = Nc.HICosmoQSpline.new(Ncm.SplineCubicNotaknot.new(), 8, 2.0)
    cosmo.props.qparam = Ncm.Vector.new_array([0.5] * 8)
    assert np.isnan(cosmo.q_transition())


def test_qspline_band_function_grid(tmp_path: Path) -> None:
    """A bound-z q(z) function evaluates via eval0 and round-trips through YAML."""
    cosmo = Nc.HICosmoQSpline.new(Ncm.SplineCubicNotaknot.new(), 8, 2.0)
    cosmo.props.qparam = Ncm.Vector.new_array(list(np.linspace(0.5, -0.6, 8)))
    mset = Ncm.MSet.new_array([cosmo])

    oa = Ncm.ObjArray.new()
    z_nodes = [0.0, 0.5, 1.0, 1.5]
    for z in z_nodes:
        func = Ncm.MSetFuncList.new("NcHICosmo:q", None)
        func.set_eval_x([z])
        oa.add(func)

    direct = [cosmo.q(z) for z in z_nodes]
    assert_allclose([oa.get(i).eval0(mset) for i in range(oa.len())], direct)

    path = tmp_path / "funcs.yaml"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    ser.array_to_yaml_file(oa, path.absolute().as_posix())
    reloaded = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).array_from_yaml_file(
        path.absolute().as_posix()
    )
    assert_allclose(
        [reloaded.get(i).eval0(mset) for i in range(reloaded.len())], direct
    )


def test_qspline_knots_default_is_uniform() -> None:
    """The default q(z) knot placement is uniform in z."""
    cosmo = Nc.HICosmoQSpline.new(Ncm.SplineCubicNotaknot.new(), 8, 2.0)
    assert cosmo.props.knots == Nc.HICosmoSplineKnots.UNIFORM

    z = np.array(cosmo.props.spline.peek_xv().dup_array())
    assert_allclose(np.diff(z), np.diff(z)[0])  # evenly spaced in z


def test_qspline_knots_chebyshev_clusters() -> None:
    """The Chebyshev option clusters q(z) knots toward the endpoints in z."""
    n, z_f = 8, 2.0
    uni = Nc.HICosmoQSpline.new(Ncm.SplineCubicNotaknot.new(), n, z_f)
    cheb = Nc.HICosmoQSpline(
        spline=Ncm.SplineCubicNotaknot.new(),
        qparam_length=n,
        zf=z_f,
        knots=Nc.HICosmoSplineKnots.CHEBYSHEV,
    )
    assert cheb.props.knots == Nc.HICosmoSplineKnots.CHEBYSHEV

    z_uni = np.array(uni.props.spline.peek_xv().dup_array())
    z_cheb = np.array(cheb.props.spline.peek_xv().dup_array())
    assert_allclose(z_cheb[0], z_uni[0])
    assert_allclose(z_cheb[-1], z_uni[-1])
    # Chebyshev's high-z gap is tighter than uniform's.
    assert (z_cheb[-1] - z_cheb[-2]) < (z_uni[-1] - z_uni[-2])
    # and the central gap exceeds the boundary gaps.
    d = np.diff(z_cheb)
    assert d[len(d) // 2] > d[0]
