#!/usr/bin/env python
#
# test_galaxy_shape_factor_moments_tilt.py
#
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
# Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Tests for NcGalaxyShapeFactorMomentsTilt.

The class solves the tilt's moment equations exactly and interpolates the
solution on a dyadic Chebyshev mesh in ghat = min(|g|, 1/|g|). These tests
target what that construction promises:

* the closed-form target moments, in both ellipticity conventions, including
  the facts that do not depend on the population (E[x] = ghat exactly under
  TRACE_DET; the collapse to (1, 1 + sn^2, sn^2) at ghat = 1);
* g = 0 returning the exact zero-shear marginal bitwise, which also makes
  NcGalaxyShapeFactorFixedQuad a valid pointwise reference there, and only
  there -- at g > 0 the tilt family and the full convolution are different
  models;
* invariance under g -> 1/g, at observed ellipticities with a cross
  component, in both conventions;
* agreement with the frozen values of the removed perturbative
  NcGalaxyShapeFactorTiltedSeries at small shear
  (data/truth_tables/moments/golden.json, written once by
  tests/tools/make_moments_compat_fixtures.py), to that class's truncation
  error;
* the reachable-range restriction and its guard, and that serialized tables
  are adopted without rebuilding and rejected when stale.
"""

import json
import os
import subprocess
import sys

import numpy as np
import pytest
from scipy import integrate, special

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

_CONVS = {
    "trace": Nc.GalaxyWLObsEllipConv.TRACE,
    "trace-det": Nc.GalaxyWLObsEllipConv.TRACE_DET,
}


def _build_mset(sigma_pop):
    """Same recipe as tests/tools/make_moments_compat_fixtures.py."""
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma_pop
    hms.param_set_by_name("log10MDelta", 14.0)
    hp.param_set_by_name("z", 0.2)
    hp.prepare(cosmo)
    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop):
        mset.set(model)
    mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())
    mset.set(Nc.GalaxyRedshiftObsGauss.new())
    return mset, pop


def _make_data(gsf, mset, std_noise, ra=0.0, dec=0.0):
    posf = Nc.GalaxyPositionFactorFlat.new(-1.0, 1.0, -1.0, 1.0)
    pos_data = Nc.GalaxyPositionFactorData.new(posf, mset)
    pos_data.ra, pos_data.dec = ra, dec
    zf = Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0)
    z_data = Nc.GalaxyRedshiftFactorData.new(zf, mset)
    data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)
    gsf.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsf.prepare_data_array(mset, [data], True, True)
    return data, (pos_data, z_data)


@pytest.fixture(name="mset_pop")
def fixture_mset_pop():
    return _build_mset(0.3)


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_construct(conv):
    gsf = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    assert isinstance(gsf, Nc.GalaxyShapeFactor)
    assert gsf.props.ellip_conv == conv
    assert "gate=" in gsf.get_desc()


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_boundary_moments(conv, mset_pop):
    """At ghat = 1 both maps collapse every source to 1."""
    mset, pop = mset_pop
    gsf = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    sn = 0.1
    data, _keep = _make_data(gsf, mset, sn)
    mu, ex2, ey2 = gsf.exact_moments(pop, data, 1.0)
    assert mu == pytest.approx(1.0, abs=1.0e-12)
    assert ex2 == pytest.approx(1.0 + sn * sn, abs=1.0e-12)
    assert ey2 == pytest.approx(sn * sn, abs=1.0e-12)


def test_trace_det_mean_is_ghat(mset_pop):
    """epsilon is an unbiased shear estimator: E[x] = ghat exactly."""
    mset, pop = mset_pop
    gsf = Nc.GalaxyShapeFactorMomentsTilt.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _keep = _make_data(gsf, mset, 0.05)
    for gh in (0.0, 0.1, 0.37, 0.8, 0.999):
        assert gsf.exact_moments(pop, data, gh)[0] == gh


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_zero_shear_is_exact(conv, mset_pop):
    """At g = 0 the tilt vanishes bitwise and the marginal is P_0 itself."""
    mset, pop = mset_pop
    gsf = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    data, _keep = _make_data(gsf, mset, 0.1)
    assert gsf.eval_tilt(pop, data, 0.0, 0.0) == (0.0, 0.0, 0.0, 0.0)

    ref = Nc.GalaxyShapeFactorFixedQuad.new(conv)
    rdata, _keep2 = _make_data(ref, mset, 0.1)
    for r in (0.05, 0.3, 0.6):
        a = gsf.eval_ln_marginal(pop, data, 0.0, 0.0, r, 0.0)
        b = ref.eval_ln_marginal(pop, rdata, 0.0, 0.0, r, 0.0)
        assert a == pytest.approx(b, rel=1.0e-6)


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_duality(conv, mset_pop):
    """P(g) = P(1/g), including a cross component, which the fold's
    conjugation of the observed ellipticity flips."""
    mset, pop = mset_pop
    gsf = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    data, _keep = _make_data(gsf, mset, 0.1)
    for g in (0.1, 0.5, 0.9):
        for e1, e2 in ((0.3, 0.25), (-0.5, -0.4), (0.1, -0.7)):
            a = gsf.eval_ln_marginal(pop, data, g, 0.0, e1, e2)
            b = gsf.eval_ln_marginal(pop, data, 1.0 / g, 0.0, e1, e2)
            assert a == pytest.approx(b, abs=1.0e-11)


def _golden():
    filename = Ncm.cfg_get_data_filename("truth_tables/moments/golden.json", True)
    with open(filename, encoding="utf-8") as f:
        return json.load(f)


@pytest.mark.parametrize("conv_name", _CONVS.keys())
def test_matches_frozen_tilted_series(conv_name):
    """Same object as the perturbative class it replaced, to its truncation
    error -- the check that keeps the TRACE_DET derivation honest."""
    golden = _golden()
    tol = golden["tol"]["tilted_series"][conv_name]
    conv = _CONVS[conv_name]
    cache = {}
    for case in golden["cases"]:
        if case["conv"] != conv_name:
            continue
        key = (case["sigma_pop"], case["std_noise"])
        if key not in cache:
            mset, pop = _build_mset(case["sigma_pop"])
            gsf = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
            data, keep = _make_data(gsf, mset, case["std_noise"])
            cache[key] = (gsf, pop, data, keep, mset)
        gsf, pop, data, _keep, _mset = cache[key]
        v = gsf.eval_ln_marginal(pop, data, case["g"], 0.0, case["e1"], case["e2"])
        assert abs(v - case["tilted_series"]) <= tol[str(case["g"])], case


def test_reachable_range_and_guard():
    """A galaxy that cannot reach the critical curve gets a truncated mesh,
    equal to the full one where both are defined, and refuses to
    extrapolate past it."""
    mset, pop = _build_mset(0.4658)
    hms = Nc.HaloDensityProfile.peek_mass_summary(mset.peek(Nc.HaloDensityProfile.id()))
    hms.param_set_desc("log10MDelta", {"upper-bound": 14.5})

    rest = Nc.GalaxyShapeFactorMomentsTilt(ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE)
    full = Nc.GalaxyShapeFactorMomentsTilt(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE, restrict_range=False
    )
    d_r, _a = _make_data(rest, mset, 0.1, ra=0.8)
    d_f, _b = _make_data(full, mset, 0.1, ra=0.8)
    rest.prepare(mset)
    rest.data_prepare(mset, d_r, 3.0)

    n_r, top, _deg = rest.peek_layout(pop, d_r)
    n_f, top_f, _deg_f = full.peek_layout(pop, d_f)
    assert top < 1.0 and top_f == 1.0 and n_r < n_f

    for gh in np.linspace(0.0, top, 40):
        a = rest.eval_ln_marginal(pop, d_r, gh, 0.0, 0.4, 0.1)
        b = full.eval_ln_marginal(pop, d_f, gh, 0.0, 0.4, 0.1)
        assert a == b

    rest.reset_range_error_count()
    assert rest.eval_marginal(pop, d_r, top + 0.05, 0.0, 0.4, 0.1) == 0.0
    assert rest.get_range_error_count() == 1


def test_range_covers_shear_calibration():
    """The table is read at the calibrated shear (1 + m) g_t e^{2i phi} + c,
    so the reachable range must bound that, not the bare reduced shear. At
    this radius the bare bound sits just under the first panel top, 0.5, and
    a 10% multiplicative bias carries the shear past it."""
    mset, pop = _build_mset(0.4658)
    hms = Nc.HaloDensityProfile.peek_mass_summary(mset.peek(Nc.HaloDensityProfile.id()))
    hms.param_set_desc("log10MDelta", {"upper-bound": 15.5})
    conv = Nc.GalaxyWLObsEllipConv.TRACE

    def layout(c1, c2, m):
        gsf = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
        d, _keep = _make_data(gsf, mset, 0.1, ra=0.042)
        gsf.data_set(d, 0.0, 0.0, 0.1, c1, c2, m, Nc.WLEllipticityFrame.CELESTIAL)
        gsf.prepare(mset)
        gsf.data_prepare(mset, d, 3.0)
        n, top, _deg = gsf.peek_layout(pop, d)
        return n, top

    assert layout(0.0, 0.0, 0.0) == (1, 0.5)
    assert layout(0.0, 0.0, 0.1) == (2, 0.75)
    # |c| = 1 alone can take the shear to the critical curve.
    assert layout(0.6, 0.8, 0.0)[1] == 1.0
    # A negative m only tightens the bound.
    assert layout(0.0, 0.0, -0.5) == (1, 0.5)


def test_serialized_tables():
    """Stored tables are adopted without rebuilding; a different population
    rejects them."""
    mset, pop = _build_mset(0.3)
    src = Nc.GalaxyShapeFactorMomentsTilt.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    src.prepare(mset)
    data, _keep = _make_data(src, mset, 0.1)
    ref = src.eval_ln_marginal(pop, data, 0.4, 0.0, 0.3, -0.2)
    blob = Ncm.Serialize.new(Ncm.SerializeOpt.NONE).to_string(src, True)

    dst = Ncm.Serialize.new(Ncm.SerializeOpt.NONE).from_string(blob)
    dst.prepare(mset)
    d2, _keep2 = _make_data(dst, mset, 0.1)
    assert dst.eval_ln_marginal(pop, d2, 0.4, 0.0, 0.3, -0.2) == ref
    assert dst.get_table_build_count() == 0

    mset2, pop2 = _build_mset(0.2)
    stale = Ncm.Serialize.new(Ncm.SerializeOpt.NONE).from_string(blob)
    stale.prepare(mset2)
    d3, _keep3 = _make_data(stale, mset2, 0.1)
    stale.eval_ln_marginal(pop2, d3, 0.4, 0.0, 0.3, -0.2)
    assert stale.get_table_build_count() == 1


@pytest.mark.parametrize(
    "cls",
    [Nc.GalaxyShapeFactorMomentsTilt, Nc.GalaxyShapeFactorCGF],
    ids=["moments-tilt", "cgf-no-hook"],
)
def test_data_prefetch_is_inert(cls, mset_pop):
    """Every prefetch stage is safe before the galaxy has a table, after it
    has one, and for a class without the hook, and never changes a value."""
    mset, pop = mset_pop
    gsf = cls.new(Nc.GalaxyWLObsEllipConv.TRACE)
    data, _keep = _make_data(gsf, mset, 0.1)
    points = [(0.05, 0.0, 0.3, 0.25), (0.4, 0.1, -0.5, 0.2)]

    for stage in (0, 1, 2):
        gsf.data_prefetch(data, stage)
    before = [gsf.eval_ln_marginal(pop, data, *p) for p in points]
    for stage in (0, 1, 2):
        gsf.data_prefetch(data, stage)
    after = [gsf.eval_ln_marginal(pop, data, *p) for p in points]

    assert before == after


def _polar_rule(sn, n_gl=40, n_phi=96):
    """Polar product rule on the observed-ellipticity plane, fine across the
    |chi| = 1 layer the noise smears, and out to 12 sigma_nu beyond it."""
    edges = np.concatenate(
        [
            np.linspace(0.0, 0.8, 5),
            np.linspace(0.8, 1.2, 9)[1:],
            np.linspace(1.2, 1.0 + 12.0 * sn, 3)[1:],
        ]
    )
    x, w = np.polynomial.legendre.leggauss(n_gl)
    r = np.concatenate(
        [0.5 * (a + b) + 0.5 * (b - a) * x for a, b in zip(edges[:-1], edges[1:])]
    )
    wr = np.concatenate([0.5 * (b - a) * w for a, b in zip(edges[:-1], edges[1:])])
    phi = 2.0 * np.pi * np.arange(n_phi) / n_phi
    return r, wr * r * (2.0 * np.pi / n_phi), phi


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_tilt_is_exponential_family(conv, mset_pop):
    """ln P(g) - ln P(0) = lambda.T - W, with T = (x, x^2, y^2) in the frame
    of g: eval_tilt() returns exactly what the marginal uses."""
    mset, pop = mset_pop
    gsf = Nc.GalaxyShapeFactorMomentsTilt(ellip_conv=conv, restrict_range=False)
    data, _keep = _make_data(gsf, mset, 0.1)
    for g in (0.2, 0.6, 0.9):
        l1, l2, l3, W = gsf.eval_tilt(pop, data, g, 0.0)
        for e1, e2 in ((0.3, 0.2), (-0.4, 0.5), (0.05, -0.1)):
            a = gsf.eval_ln_marginal(pop, data, g, 0.0, e1, e2)
            a -= gsf.eval_ln_marginal(pop, data, 0.0, 0.0, e1, e2)
            b = l1 * e1 + l2 * e1 * e1 + l3 * e2 * e2 - W
            assert a == pytest.approx(b, abs=1.0e-12)
        # g and 1/g are the same point of the fold.
        assert gsf.eval_tilt(pop, data, 1.0 / g, 0.0) == pytest.approx(
            (l1, l2, l3, W), rel=1.0e-12, abs=1.0e-14
        )


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_normalized_and_moment_matched(conv, mset_pop):
    """The model integrates to one and reproduces the exact moments it was
    matched to, both to far better than the accuracy gate; including the
    part of the plane beyond |chi| = 1 that only the noise reaches."""
    mset, pop = mset_pop
    sn = 0.1
    gsf = Nc.GalaxyShapeFactorMomentsTilt(ellip_conv=conv, restrict_range=False)
    data, _keep = _make_data(gsf, mset, sn)
    r, w, phi = _polar_rule(sn)
    gate = gsf.props.accuracy_gate
    for g in (0.2, 0.6, 0.9):
        Z, M = 0.0, np.zeros(3)
        for ri, wi in zip(r, w):
            for p in phi:
                e1, e2 = ri * np.cos(p), ri * np.sin(p)
                P = gsf.eval_marginal(pop, data, g, 0.0, e1, e2)
                Z += wi * P
                M += wi * P * np.array([e1, e1 * e1, e2 * e2])
        assert abs(np.log(Z)) < 0.1 * gate
        np.testing.assert_allclose(
            M / Z, gsf.exact_moments(pop, data, g), rtol=0.0, atol=1.0e-5
        )


def test_matches_gauss_moments(mset_pop):
    """Both classes match the same closed-form moments."""
    mset, pop = mset_pop
    for conv in _CONVS.values():
        tilt = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
        gauss = Nc.GalaxyShapeFactorMomentsGauss.new(conv)
        d_t, _a = _make_data(tilt, mset, 0.15)
        d_g, _b = _make_data(gauss, mset, 0.15)
        for gh in (0.0, 0.3, 0.7, 1.0):
            mu, ex2, ey2 = tilt.exact_moments(pop, d_t, gh)
            mu_g, var_g, ey2_g = gauss.exact_moments(pop, d_g, gh)
            assert mu_g == pytest.approx(mu, abs=1.0e-14)
            assert var_g + mu_g * mu_g == pytest.approx(ex2, abs=1.0e-14)
            assert ey2_g == pytest.approx(ey2, abs=1.0e-14)


def _ln_P0_reference(pop, pop_data, R, sn):
    """ln P_0(R): the population's r-marginal convolved with the isotropic
    noise, the angle integrated out analytically into I_0."""

    def f(r):
        return (
            pop.eval_p(pop_data, r)
            * np.exp(-((R - r) ** 2) / (2.0 * sn * sn))
            * special.ive(0, R * r / (sn * sn))
            / (2.0 * np.pi * sn * sn)
        )

    brk = [max(0.0, R - 5.0 * sn), min(1.0, R + 5.0 * sn)]
    val, _err = integrate.quad(
        f, 0.0, 1.0, points=brk, limit=400, epsabs=0.0, epsrel=1.0e-12
    )
    return np.log(val)


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_zero_shear_marginal_up_to_and_past_unit_circle(conv, mset_pop):
    """At g = 0 the marginal is the exact noisy P_0 at every radius,
    including right at |chi| = 1 and past it, where only the noise reaches
    and where a fixed rule in the source ellipticity is least accurate."""
    mset, pop = mset_pop
    sn = 0.1
    gsf = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    data, _keep = _make_data(gsf, mset, sn)
    for R in (0.3, 0.9, 0.999, 1.001, 1.05, 1.3):
        a = gsf.eval_ln_marginal(pop, data, 0.0, 0.0, R, 0.0)
        assert a == pytest.approx(
            _ln_P0_reference(pop, data.pop_data, R, sn), abs=1.0e-7
        )


def test_tighter_gate_refines_and_agrees(mset_pop):
    """A tighter accuracy gate raises the local degrees and moves ln P by no
    more than the looser gate allows."""
    mset, pop = mset_pop
    conv = Nc.GalaxyWLObsEllipConv.TRACE
    loose = Nc.GalaxyShapeFactorMomentsTilt(ellip_conv=conv, restrict_range=False)
    tight = Nc.GalaxyShapeFactorMomentsTilt(
        ellip_conv=conv, restrict_range=False, accuracy_gate=1.0e-9, max_degree=64
    )
    d_l, _a = _make_data(loose, mset, 0.05)
    d_t, _b = _make_data(tight, mset, 0.05)
    n_l, _top_l, deg_l = loose.peek_layout(pop, d_l)
    n_t, _top_t, deg_t = tight.peek_layout(pop, d_t)
    assert n_l == n_t
    assert all(t >= l for t, l in zip(deg_t, deg_l))
    assert sum(deg_t) > sum(deg_l)
    for g in (0.1, 0.5, 0.8, 0.97):
        a = loose.eval_ln_marginal(pop, d_l, g, 0.0, 0.3, 0.1)
        b = tight.eval_ln_marginal(pop, d_t, g, 0.0, 0.3, 0.1)
        assert a == pytest.approx(b, abs=loose.props.accuracy_gate)


def test_counters():
    """The build and solve counters count and reset."""
    mset, pop = _build_mset(0.3)
    gsf = Nc.GalaxyShapeFactorMomentsTilt.new(Nc.GalaxyWLObsEllipConv.TRACE)
    data, _keep = _make_data(gsf, mset, 0.1)
    gsf.eval_ln_marginal(pop, data, 0.4, 0.0, 0.3, 0.1)
    assert gsf.get_table_build_count() == 1
    assert gsf.get_solve_error_count() == 0
    gsf.reset_table_build_count()
    gsf.reset_solve_error_count()
    assert gsf.get_table_build_count() == 0
    gsf.eval_ln_marginal(pop, data, 0.5, 0.0, 0.3, 0.1)
    assert gsf.get_table_build_count() == 0


# A narrow population with little noise: the moment equations become too
# stiff for the Newton solve at some nodes. Found by scanning (sigma_pop,
# sigma_nu); TRACE_DET only.
_STIFF = (0.01, 0.01)


def test_solve_failure_is_counted():
    """Failed nodes are counted, the table is still built and evaluates to
    finite numbers, and the count resets."""
    mset, pop = _build_mset(_STIFF[0])
    gsf = Nc.GalaxyShapeFactorMomentsTilt(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET, restrict_range=False
    )
    data, _keep = _make_data(gsf, mset, _STIFF[1])
    v = gsf.eval_ln_marginal(pop, data, 0.99, 0.0, 0.3, 0.1)
    assert np.isfinite(v)
    n_failed = gsf.get_solve_error_count()
    assert n_failed > 0

    gsf.reset_solve_error_count()
    assert gsf.get_solve_error_count() == 0


def test_strict_solve_aborts():
    """With strict-solve a failed node is fatal."""
    code = f"""
import sys
sys.path.insert(0, {repr(os.path.dirname(__file__))})
from numcosmo_py import Nc
import test_galaxy_shape_factor_moments_tilt as T
mset, pop = T._build_mset({_STIFF[0]})
gsf = Nc.GalaxyShapeFactorMomentsTilt(
    ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET, restrict_range=False, strict_solve=True
)
data, keep = T._make_data(gsf, mset, {_STIFF[1]})
gsf.eval_ln_marginal(pop, data, 0.99, 0.0, 0.3, 0.1)
"""
    res = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    assert res.returncode != 0
    assert "Newton" in res.stderr


def test_range_edges():
    """The reachable range at its edges: a foreground galaxy gets the
    smallest mesh, a halo centre free enough to reach the galaxy the full
    one, and without prepare() there is no bound at all."""
    mset, pop = _build_mset(0.3)
    hp = mset.peek(Nc.HaloPosition.id())
    conv = Nc.GalaxyWLObsEllipConv.TRACE

    gsf = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    d, _a = _make_data(gsf, mset, 0.1, ra=0.3)
    gsf.prepare(mset)
    gsf.data_prepare(mset, d, 0.1)  # z_max below z_cl = 0.2
    n, top, _deg = gsf.peek_layout(pop, d)
    assert (n, top) == (1, 0.5)
    assert gsf.eval_ln_marginal(pop, d, 0.0, 0.0, 0.3, 0.1) == pytest.approx(
        Nc.GalaxyShapeFactorMomentsTilt(
            ellip_conv=conv, restrict_range=False
        ).eval_ln_marginal(pop, _make_data(gsf, mset, 0.1)[0], 0.0, 0.0, 0.3, 0.1),
        rel=1.0e-12,
    )

    # A support with no upper end: the shear there is not finite, so it
    # bounds nothing and the galaxy keeps the full mesh.
    unbounded = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    d_inf, _i = _make_data(unbounded, mset, 0.1, ra=0.3)
    unbounded.prepare(mset)
    unbounded.data_prepare(mset, d_inf, np.inf)
    assert unbounded.peek_layout(pop, d_inf)[1] == 1.0

    fresh = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    pos = Nc.GalaxyPositionFactorData.new(
        Nc.GalaxyPositionFactorFlat.new(-1.0, 1.0, -1.0, 1.0), mset
    )
    pos.ra, pos.dec = 0.3, 0.0
    zd = Nc.GalaxyRedshiftFactorData.new(
        Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0), mset
    )
    d2 = Nc.GalaxyShapeFactorData.new(fresh, mset, pos, zd)
    fresh.data_set(d2, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    fresh.data_prepare(mset, d2, 3.0)  # never prepared: no bound
    assert fresh.peek_layout(pop, d2)[1] == 1.0

    for name in ("ra", "dec"):
        hp.param_set_desc(name, {"fit": True, "lower-bound": -0.5, "upper-bound": 0.5})
    free = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    d3, _c = _make_data(free, mset, 0.1, ra=0.3)
    free.prepare(mset)
    free.data_prepare(mset, d3, 3.0)
    assert free.peek_layout(pop, d3)[1] == 1.0


def test_range_without_free_concentration():
    """A concentration-mass relation has no cDelta to put at its bounds;
    the bound is still taken, and still agrees with the full mesh."""
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMDuffy08.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hms.param_set_desc("log10MDelta", {"value": 14.0, "upper-bound": 14.5})
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    hp.param_set_by_name("z", 0.2)
    hp.prepare(cosmo)
    pop = Nc.GalaxyShapePopGauss.new()
    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, Nc.WLSurfaceMassDensity.new(dist), pop):
        mset.set(model)
    mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())
    mset.set(Nc.GalaxyRedshiftObsGauss.new())

    conv = Nc.GalaxyWLObsEllipConv.TRACE
    rest = Nc.GalaxyShapeFactorMomentsTilt.new(conv)
    full = Nc.GalaxyShapeFactorMomentsTilt(ellip_conv=conv, restrict_range=False)
    d_r, _a = _make_data(rest, mset, 0.1, ra=0.8)
    d_f, _b = _make_data(full, mset, 0.1, ra=0.8)
    rest.prepare(mset)
    rest.data_prepare(mset, d_r, 3.0)
    _n, top, _deg = rest.peek_layout(pop, d_r)
    assert top < 1.0
    for gh in np.linspace(0.0, top, 7):
        assert rest.eval_ln_marginal(
            pop, d_r, gh, 0.0, 0.4, 0.1
        ) == full.eval_ln_marginal(pop, d_f, gh, 0.0, 0.4, 0.1)


def _stored_tables(conv):
    """A built instance's tables, its stamp and a reference value."""
    mset, pop = _build_mset(0.3)
    src = Nc.GalaxyShapeFactorMomentsTilt(ellip_conv=conv, restrict_range=False)
    src.prepare(mset)
    data, _keep = _make_data(src, mset, 0.1)
    ref = src.eval_ln_marginal(pop, data, 0.4, 0.0, 0.3, -0.2)
    tables = src.props.tables
    assert tables.len() == 1
    return mset, pop, tables.peek(0).dup_array(), src.props.tables_stamp, ref


def _corruptions(v):
    """One malformed copy of a stored table per check the loader makes."""
    h, n_panels = 6, int(v[3])
    yield "short", v[:h]
    for k, bad in (
        (0, 0.0),
        (0, 5.0),
        (3, 0.0),
        (3, 9.0),
        (4, 0.0),
        (4, 1.5),
        (h, 0.0),
        (h, 256.0),
    ):
        c = list(v)
        c[k] = bad
        yield f"v[{k}]={bad}", c
    yield "long", v + [0.0]
    assert n_panels >= 1


def test_corrupted_tables_are_rebuilt():
    """A stored table whose shape does not add up is rebuilt, not trusted;
    a non-vector entry is skipped; a sound one is adopted."""
    conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    mset, pop, v, stamp, ref = _stored_tables(conv)

    def load(entries):
        arr = Ncm.ObjArray.new()
        for e in entries:
            arr.add(e)
        gsf = Nc.GalaxyShapeFactorMomentsTilt(ellip_conv=conv, restrict_range=False)
        gsf.props.tables_stamp = stamp
        gsf.props.tables = arr
        gsf.prepare(mset)
        data, keep = _make_data(gsf, mset, 0.1)
        val = gsf.eval_ln_marginal(pop, data, 0.4, 0.0, 0.3, -0.2)
        return gsf.get_table_build_count(), val

    assert load([Ncm.Matrix.new(1, 1), Ncm.Vector.new_array(v)]) == (0, ref)
    for name, bad in _corruptions(v):
        assert load([Ncm.Vector.new_array(bad)]) == (1, ref), name


def test_narrow_population_small_noise():
    """A narrow population seen with little noise: the moment solve needs its
    finest angular rule, yet every node converges, the tilt identity and the
    g -> 1/g symmetry hold, and, given the degree it asks for (the default
    cap of 16 is sized for real shape noise, see
    test_small_noise_realistic_population), the table is within its gate
    of a much tighter one."""
    mset, pop = _build_mset(0.05)
    conv = Nc.GalaxyWLObsEllipConv.TRACE
    gsf = Nc.GalaxyShapeFactorMomentsTilt(
        ellip_conv=conv, restrict_range=False, max_degree=64
    )
    ref = Nc.GalaxyShapeFactorMomentsTilt(
        ellip_conv=conv, restrict_range=False, accuracy_gate=1.0e-6, max_degree=64
    )
    data, _keep = _make_data(gsf, mset, 0.01)
    d_ref, _keep2 = _make_data(ref, mset, 0.01)
    e1, e2 = 0.06, 0.03
    for g in (0.02, 0.1, 0.5):
        l1, l2, l3, W = gsf.eval_tilt(pop, data, g, 0.0)
        a = gsf.eval_ln_marginal(pop, data, g, 0.0, e1, e2)
        b = gsf.eval_ln_marginal(pop, data, 0.0, 0.0, e1, e2)
        b += l1 * e1 + l2 * e1 * e1 + l3 * e2 * e2 - W
        assert a == pytest.approx(b, abs=1.0e-9)
        assert a == pytest.approx(
            gsf.eval_ln_marginal(pop, data, 1.0 / g, 0.0, e1, e2), abs=1.0e-9
        )
        assert a == pytest.approx(
            ref.eval_ln_marginal(pop, d_ref, g, 0.0, e1, e2),
            abs=gsf.props.accuracy_gate,
        )
    assert gsf.get_solve_error_count() == 0
    assert ref.get_solve_error_count() == 0


def _fine_polar_rule(sn, n_rad=40, n_phi=128):
    """As _polar_rule, resolving the sigma_nu-wide ridge at |chi| = 1 in both
    radius and angle: at small noise and large shear the density crowds
    against the unit circle around the shear axis."""
    edges = np.unique(
        np.concatenate(
            [
                np.linspace(0.0, 1.0 - 10.0 * sn, 16),
                np.linspace(1.0 - 10.0 * sn, 1.0 + 12.0 * sn, n_rad),
            ]
        )
    )
    x, w = np.polynomial.legendre.leggauss(12)
    r = np.concatenate(
        [0.5 * (a + b) + 0.5 * (b - a) * x for a, b in zip(edges[:-1], edges[1:])]
    )
    wr = np.concatenate([0.5 * (b - a) * w for a, b in zip(edges[:-1], edges[1:])])
    p_edges = np.unique(
        np.concatenate(
            [
                np.linspace(-np.pi, np.pi, n_phi // 8 + 1),
                np.linspace(-0.2, 0.2, n_phi // 8 + 1),
            ]
        )
    )
    xp, wp = np.polynomial.legendre.leggauss(8)
    phi = np.concatenate(
        [0.5 * (a + b) + 0.5 * (b - a) * xp for a, b in zip(p_edges[:-1], p_edges[1:])]
    )
    wphi = np.concatenate(
        [0.5 * (b - a) * wp for a, b in zip(p_edges[:-1], p_edges[1:])]
    )
    return r, wr * r, phi, wphi


def test_small_noise_realistic_population():
    """Real catalogues reach sigma_nu ~ 0.01 with e_rms ~ 0.4: there the
    default table is normalized and within its gate of a much tighter one,
    around the lensed mean where the likelihood is read."""
    mset, pop = _build_mset(0.3)
    sn = 0.01
    conv = Nc.GalaxyWLObsEllipConv.TRACE
    gsf = Nc.GalaxyShapeFactorMomentsTilt(ellip_conv=conv, restrict_range=False)
    ref = Nc.GalaxyShapeFactorMomentsTilt(
        ellip_conv=conv, restrict_range=False, accuracy_gate=1.0e-7, max_degree=64
    )
    data, _keep = _make_data(gsf, mset, sn)
    d_ref, _keep2 = _make_data(ref, mset, sn)
    gate = gsf.props.accuracy_gate
    r, wr, phi, wphi = _fine_polar_rule(sn)
    cphi, sphi = np.cos(phi), np.sin(phi)
    for g in (0.3, 0.97):
        Z = 0.0
        for ri, wi in zip(r, wr):
            Z += wi * sum(
                wj * gsf.eval_marginal(pop, data, g, 0.0, ri * c, ri * s)
                for c, s, wj in zip(cphi, sphi, wphi)
            )
        assert abs(np.log(Z)) < 0.1 * gate

        mu, ex2, ey2 = gsf.exact_moments(pop, data, g)
        sx, sy = np.sqrt(ex2 - mu * mu), np.sqrt(ey2)
        for dx in (-3.0, -1.0, 0.0, 1.0, 3.0):
            for dy in (-3.0, 0.0, 3.0):
                x, y = mu + dx * sx, dy * sy
                assert gsf.eval_ln_marginal(pop, data, g, 0.0, x, y) == pytest.approx(
                    ref.eval_ln_marginal(pop, d_ref, g, 0.0, x, y), abs=gate
                )


@pytest.mark.parametrize(
    "sn,n_panels", [(1.0, 2), (0.4, 2), (0.1, 3), (0.01, 7), (1.0e-4, 8)]
)
def test_panel_count_rule(sn, n_panels, mset_pop):
    """K = max(2, round(log2(1/sigma_nu))) dyadic panels, capped at the
    layout's eight, over all of [0, 1] when the range is not restricted."""
    mset, pop = mset_pop
    gsf = Nc.GalaxyShapeFactorMomentsTilt(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE, restrict_range=False
    )
    data, _keep = _make_data(gsf, mset, sn)
    n, top, deg = gsf.peek_layout(pop, data)
    assert (n, top, len(deg)) == (n_panels, 1.0, n_panels)
