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

import numpy as np
import pytest

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
