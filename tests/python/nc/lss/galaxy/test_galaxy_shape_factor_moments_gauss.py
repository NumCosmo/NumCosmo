#!/usr/bin/env python
#
# test_galaxy_shape_factor_moments_gauss.py
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

"""Tests for NcGalaxyShapeFactorMomentsGauss.

A Gaussian matched to the exact moments of the lensed, noisy marginal. These
tests check that the tabulated moments are the exact ones to the table's
tolerance, right up to the critical curve; that the tables do not depend on
the noise, so one serves a whole population; the TRACE_DET structure (mean
ghat, isotropic covariance); invariance under g -> 1/g; and agreement with
the frozen values of the removed truncated-series Gaussian,
NcGalaxyShapeFactorMomentSeries, to its truncation error
(data/truth_tables/moments/golden.json).
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


def _golden():
    filename = Ncm.cfg_get_data_filename("truth_tables/moments/golden.json", True)
    with open(filename, encoding="utf-8") as f:
        return json.load(f)


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_table_is_exact(conv):
    """Interpolated moments against direct evaluation, up to ghat = 1."""
    mset, pop = _build_mset(0.4658)
    gsf = Nc.GalaxyShapeFactorMomentsGauss.new(conv)
    data, _keep = _make_data(gsf, mset, 0.02)
    gh = np.concatenate(
        [np.linspace(0.0, 0.99, 120), 1.0 - np.logspace(-2, -9, 40), [1.0]]
    )
    for g in gh:
        a = np.array(gsf.eval_moments(pop, data, float(g), 0.0))
        b = np.array(gsf.exact_moments(pop, data, float(g)))
        assert np.max(np.abs(a - b)) < 1.0e-9, g


def test_one_table_per_population():
    mset, pop = _build_mset(0.3)
    gsf = Nc.GalaxyShapeFactorMomentsGauss.new(Nc.GalaxyWLObsEllipConv.TRACE)
    keep = []
    for sn in np.linspace(0.008, 0.4, 20):
        data, k = _make_data(gsf, mset, float(sn))
        keep.append((data, k))
        gsf.eval_ln_marginal(pop, data, 0.3, 0.0, 0.2, 0.1)
    assert gsf.get_table_build_count() == 1


def test_trace_det_structure():
    """Under TRACE_DET the matched Gaussian is centred on ghat and isotropic."""
    mset, pop = _build_mset(0.3)
    gsf = Nc.GalaxyShapeFactorMomentsGauss.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _keep = _make_data(gsf, mset, 0.05)
    for g in (0.0, 0.2, 0.6, 0.95):
        mu, var_x, e_y2 = gsf.exact_moments(pop, data, g)
        assert mu == g
        assert var_x == pytest.approx(e_y2, abs=1.0e-14)


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_duality(conv):
    mset, pop = _build_mset(0.3)
    gsf = Nc.GalaxyShapeFactorMomentsGauss.new(conv)
    data, _keep = _make_data(gsf, mset, 0.1)
    for g in (0.1, 0.5, 0.9):
        for e1, e2 in ((0.3, 0.25), (-0.5, -0.4), (0.1, -0.7)):
            a = gsf.eval_ln_marginal(pop, data, g, 0.0, e1, e2)
            b = gsf.eval_ln_marginal(pop, data, 1.0 / g, 0.0, e1, e2)
            assert a == pytest.approx(b, abs=1.0e-11)


@pytest.mark.parametrize("conv_name", _CONVS.keys())
def test_matches_frozen_moment_series(conv_name):
    golden = _golden()
    tol = golden["tol"]["moment_series"][conv_name]
    conv = _CONVS[conv_name]
    cache = {}
    for case in golden["cases"]:
        if case["conv"] != conv_name:
            continue
        key = (case["sigma_pop"], case["std_noise"])
        if key not in cache:
            mset, pop = _build_mset(case["sigma_pop"])
            gsf = Nc.GalaxyShapeFactorMomentsGauss.new(conv)
            data, keep = _make_data(gsf, mset, case["std_noise"])
            cache[key] = (gsf, pop, data, keep, mset)
        gsf, pop, data, _keep, _mset = cache[key]
        v = gsf.eval_ln_marginal(pop, data, case["g"], 0.0, case["e1"], case["e2"])
        assert abs(v - case["moment_series"]) <= tol[str(case["g"])], case


def test_data_prefetch_is_inert():
    """Every prefetch stage is safe before the galaxy has a table and after
    it has one, and never changes a value."""
    mset, pop = _build_mset(0.3)
    gsf = Nc.GalaxyShapeFactorMomentsGauss.new(Nc.GalaxyWLObsEllipConv.TRACE)
    data, _keep = _make_data(gsf, mset, 0.1)
    points = [(0.05, 0.0, 0.3, 0.25), (0.4, 0.1, -0.5, 0.2), (0.9, 0.0, 0.1, 0.0)]

    for stage in (0, 1, 2):
        gsf.data_prefetch(data, stage)
    before = [gsf.eval_ln_marginal(pop, data, *p) for p in points]
    for stage in (0, 1, 2):
        gsf.data_prefetch(data, stage)
    after = [gsf.eval_ln_marginal(pop, data, *p) for p in points]

    assert before == after


@pytest.mark.parametrize("conv", _CONVS.values(), ids=_CONVS.keys())
def test_is_the_matched_gaussian(conv):
    """ln P is the normalized Gaussian of the tabulated moments, in the
    frame rotated by -arg g, for a shear with both components."""
    mset, pop = _build_mset(0.3)
    gsf = Nc.GalaxyShapeFactorMomentsGauss.new(conv)
    data, _keep = _make_data(gsf, mset, 0.1)
    for gmod, ang in ((0.3, np.pi / 3.0), (0.8, -2.0), (0.05, 0.4)):
        g1, g2 = gmod * np.cos(ang), gmod * np.sin(ang)
        mu, c_t, c_x = gsf.eval_moments(pop, data, g1, g2)
        for e1, e2 in ((0.3, 0.2), (-0.4, 0.5), (0.0, 0.0)):
            x = np.cos(ang) * e1 + np.sin(ang) * e2
            y = -np.sin(ang) * e1 + np.cos(ang) * e2
            ref = (
                -((x - mu) ** 2) / (2.0 * c_t)
                - y * y / (2.0 * c_x)
                - 0.5 * np.log(c_t * c_x)
                - np.log(2.0 * np.pi)
            )
            assert gsf.eval_ln_marginal(pop, data, g1, g2, e1, e2) == pytest.approx(
                ref, abs=1.0e-12
            )
            assert gsf.eval_marginal(pop, data, g1, g2, e1, e2) == pytest.approx(
                np.exp(ref), rel=1.0e-12
            )


def test_noise_enters_additively():
    """The table holds the source moments; each galaxy adds its own
    sigma_nu^2 to both second moments and leaves the mean alone."""
    mset, pop = _build_mset(0.3)
    gsf = Nc.GalaxyShapeFactorMomentsGauss.new(Nc.GalaxyWLObsEllipConv.TRACE)
    d1, _a = _make_data(gsf, mset, 0.05)
    d2, _b = _make_data(gsf, mset, 0.2)
    for g in (0.0, 0.4, 0.9):
        m1 = gsf.eval_moments(pop, d1, g, 0.0)
        m2 = gsf.eval_moments(pop, d2, g, 0.0)
        assert m2[0] == m1[0]
        assert m2[1] - m1[1] == pytest.approx(0.2**2 - 0.05**2, abs=1.0e-14)
        assert m2[2] - m1[2] == pytest.approx(0.2**2 - 0.05**2, abs=1.0e-14)


def test_properties_layout_and_counter():
    """The build knobs are properties that survive a save and load; the
    mesh covers all of [0, 1] within the degree cap, and a tighter
    tolerance needs no lower degree; the build counter resets."""
    gsf = Nc.GalaxyShapeFactorMomentsGauss(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE, moment_tol=1.0e-8, max_degree=32
    )
    assert gsf.props.moment_tol == 1.0e-8
    assert gsf.props.max_degree == 32
    dup = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).dup_obj(gsf)
    assert (dup.props.moment_tol, dup.props.max_degree) == (1.0e-8, 32)

    mset, pop = _build_mset(0.3)
    data, _keep = _make_data(gsf, mset, 0.1)
    n, top, deg = gsf.peek_layout(pop, data)
    assert top == 1.0
    assert n == len(deg) >= 2
    assert all(1 <= d <= 32 for d in deg)

    tight = Nc.GalaxyShapeFactorMomentsGauss(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE, moment_tol=1.0e-13, max_degree=32
    )
    d_t, _keep2 = _make_data(tight, mset, 0.1)
    n_t, _top_t, deg_t = tight.peek_layout(pop, d_t)
    assert n_t == n
    assert all(a >= b for a, b in zip(deg_t, deg))

    assert gsf.get_table_build_count() == 1
    gsf.reset_table_build_count()
    assert gsf.get_table_build_count() == 0
