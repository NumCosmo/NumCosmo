#!/usr/bin/env python
#
# test_galaxy_shape_factor_series_lensed.py
#
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

"""Tests for the lensed-frame series shape factor calculator, both
ellipticity conventions.

``NcGalaxyShapeFactorSeriesLensed`` is an alternative to the noise-side
scheme once shipped as ``NcGalaxyShapeFactorSeries`` (removed once this class
superseded it): instead of expanding the noise kernel in `g`
(whose Taylor coefficients scale like ``1/sigma_noise^2`` and combinatorially
blow up at small per-galaxy ``std_noise`` -- a real, non-roundoff failure
mode), this class works in the LENSED frame (same
substitution ``NcGalaxyShapeFactorQuad`` already uses) and expands the
POPULATION in `g` instead, whose coefficients scale with
``1/sigma_pop^2`` -- a population/prior parameter this project already
hard-constrains away from being pathologically small. See
docs/theory/wl_shape_marginalization_series.qmd and
dev-notes/wl_shape_series_marginalization_derivation.py (sections 9-11) for
the derivation and symbolic/numeric verification this class mirrors.

Gaussian population only in this version (matches
``NcGalaxyShapeFactorVarAdd``'s own existing precedent).
"""

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

_CONVS = [Nc.GalaxyWLObsEllipConv.TRACE_DET, Nc.GalaxyWLObsEllipConv.TRACE]


def _build_mset(pop):
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)

    hms.param_set_by_name("log10MDelta", 14.0)
    hp.param_set_by_name("z", 0.2)
    hp.prepare(cosmo)

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop):
        mset.set(model)

    mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())
    mset.set(Nc.GalaxyRedshiftObsGauss.new())

    return mset


def _build_factor_data(gsf, mset):
    posf = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    pos_data = Nc.GalaxyPositionFactorData.new(posf, mset)
    zf = Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0)
    z_data = Nc.GalaxyRedshiftFactorData.new(zf, mset)
    data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)
    return data, pos_data, z_data


def _make_series_lensed(ellip_conv, order, pop, mset, std_noise):
    gsfsl = Nc.GalaxyShapeFactorSeriesLensed.new(ellip_conv, order)
    data, _, _ = _build_factor_data(gsfsl, mset)
    gsfsl.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfsl.prepare_data_array(mset, [data], True, True)
    return gsfsl, data


def _make_quad(ellip_conv, pop, mset, std_noise):
    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfq.prepare_data_array(mset, [data], True, True)
    return gsfq, data


def _eval_both(pop, ellip_conv, order, g, eps_obs, std_noise):
    mset = _build_mset(pop)
    gsfsl, data_sl = _make_series_lensed(ellip_conv, order, pop, mset, std_noise)
    lensed_val = gsfsl.eval_marginal(
        pop, data_sl, g.real, g.imag, eps_obs.real, eps_obs.imag
    )

    gsfq, data_q = _make_quad(ellip_conv, pop, mset, std_noise)
    quad_val = gsfq.eval_marginal(
        pop, data_q, g.real, g.imag, eps_obs.real, eps_obs.imag
    )

    return lensed_val, quad_val


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_quad_moderate_g_large_noise(ellip_conv):
    """Moderate noise (std_noise=0.3) and moderate real g -- both
    conventions, default order (4). TRACE converges more slowly than
    TRACE_DET at the same order (steeper own O(g) response), hence the
    looser tolerance here than a TRACE_DET-only test would need."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    g = 0.09 + 0.0j
    eps_obs = 0.4 * np.exp(1j * 0.7)

    lensed_val, quad_val = _eval_both(pop, ellip_conv, 4, g, eps_obs, 0.3)

    assert_allclose(lensed_val, quad_val, rtol=2.0e-3)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_quad_general_complex_g(ellip_conv):
    """A genuinely complex g (not aligned with any axis) exercises the
    gauge-fixing rotation (rotate (g, eps_obs) together by -arg(g) before
    applying the real-g-only closed forms), for both conventions."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.2
    g = 0.06 - 0.04j
    eps_obs = -0.15 + 0.25j

    lensed_val, quad_val = _eval_both(pop, ellip_conv, 4, g, eps_obs, 0.25)

    assert_allclose(lensed_val, quad_val, rtol=5.0e-3)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_quad_small_std_noise_stress_case(ellip_conv):
    """THE HEADLINE CASE this class exists for: small std_noise (0.005)
    combined with moderate g -- the exact regime that made the noise-side
    scheme's expansion (once shipped as NcGalaxyShapeFactorSeries, removed
    once this class superseded it) combinatorially blow up (up to ~850%+
    error at higher truncation order) since its coefficients scale like
    1/std_noise^2. This class' coefficients scale with 1/sigma_pop^2
    instead, so accuracy here should be unaffected by how small std_noise
    is -- confirmed at the default order (4) matching Quad tightly."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.28
    std_noise = 0.005
    g = 0.09 + 0.0j
    eps_obs = 0.35 + 0.1j

    lensed_val, quad_val = _eval_both(pop, ellip_conv, 4, g, eps_obs, std_noise)

    assert_allclose(lensed_val, quad_val, rtol=1.0e-2)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_higher_order_improves_accuracy_at_moderate_g(ellip_conv):
    """Increasing trunc-order should reduce the error against Quad
    monotonically-ish (checked as order=8 clearly beating order=2, not
    requiring strict monotonicity at every intermediate step) -- for BOTH
    conventions."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    g = 0.3 + 0.0j
    eps_obs = 0.4 * np.exp(1j * 0.7)
    std_noise = 0.3

    mset = _build_mset(pop)
    gsfq, data_q = _make_quad(ellip_conv, pop, mset, std_noise)
    quad_val = gsfq.eval_marginal(
        pop, data_q, g.real, g.imag, eps_obs.real, eps_obs.imag
    )

    errs = {}

    for order in (2, 8):
        gsfsl, data_sl = _make_series_lensed(ellip_conv, order, pop, mset, std_noise)
        val = gsfsl.eval_marginal(
            pop, data_sl, g.real, g.imag, eps_obs.real, eps_obs.imag
        )
        errs[order] = abs(val - quad_val) / abs(quad_val)

    assert (
        errs[8] < errs[2] / 5.0
    ), f"order=8 ({errs[8]:.2e}) should clearly beat order=2 ({errs[2]:.2e})"


def test_caching_consistent_across_repeated_calls():
    """Simulate multiple z-nodes evaluated on the SAME per-galaxy data object
    (same R, phi, std_noise; varying only g, as NcDataClusterWLFactor's fixed
    z-nodes do) -- the per-order radial integral is cached and reused across
    calls, so this exercises that cache directly: each call must still match
    a fresh single-shot computation, not a stale or partially-updated one."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    eps_obs = 0.4 * np.exp(1j * 0.7)
    std_noise = 0.3

    mset = _build_mset(pop)
    gsfsl, data_sl = _make_series_lensed(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, 4, pop, mset, std_noise
    )

    gs = [0.01, 0.05, 0.02, 0.09, 0.001]  # deliberately non-monotonic, same phase (0)
    cached_vals = [
        gsfsl.eval_marginal(pop, data_sl, g, 0.0, eps_obs.real, eps_obs.imag)
        for g in gs
    ]

    for g, cached_val in zip(gs, cached_vals):
        gsfsl_fresh, data_fresh = _make_series_lensed(
            Nc.GalaxyWLObsEllipConv.TRACE_DET, 4, pop, mset, std_noise
        )
        fresh_val = gsfsl_fresh.eval_marginal(
            pop, data_fresh, g, 0.0, eps_obs.real, eps_obs.imag
        )
        assert_allclose(cached_val, fresh_val, rtol=1.0e-12)

    # a call with a different g-phase must still be correct (cache miss + recompute)
    gsfq, data_q = _make_quad(Nc.GalaxyWLObsEllipConv.TRACE_DET, pop, mset, std_noise)
    g_rot = 0.05 + 0.03j
    cached_after_rot = gsfsl.eval_marginal(
        pop, data_sl, g_rot.real, g_rot.imag, eps_obs.real, eps_obs.imag
    )
    quad_val = gsfq.eval_marginal(
        pop, data_q, g_rot.real, g_rot.imag, eps_obs.real, eps_obs.imag
    )
    assert_allclose(cached_after_rot, quad_val, rtol=5.0e-3)


def test_cache_invalidates_on_population_parameter_change():
    """Regression test for the cache-staleness gap found (and left unfixed,
    as a separate follow-up, in the now-removed NcGalaxyShapeFactorSeries)
    during this class' own development: NcGalaxyShapeFactorSeriesLensed's
    cache is keyed
    to nc_galaxy_shape_factor_get_pop_hash() in addition to (R,phi,
    std_noise), so changing a population parameter (sigma here) BETWEEN two
    evaluations on the SAME data object must be picked up, not silently
    served from a stale cache."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    mset = _build_mset(pop)
    eps_obs = 0.3 - 0.1j
    std_noise = 0.3
    g = 0.05
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    order = 4

    gsfsl, data_sl = _make_series_lensed(ellip_conv, order, pop, mset, std_noise)
    val_before = gsfsl.eval_marginal(pop, data_sl, g, 0.0, eps_obs.real, eps_obs.imag)

    pop["sigma"] = 0.35
    gsfsl.prepare_data_array(mset, [data_sl], True, True)
    val_after_same_data = gsfsl.eval_marginal(
        pop, data_sl, g, 0.0, eps_obs.real, eps_obs.imag
    )

    gsfsl_fresh, data_fresh = _make_series_lensed(
        ellip_conv, order, pop, mset, std_noise
    )
    val_fresh = gsfsl_fresh.eval_marginal(
        pop, data_fresh, g, 0.0, eps_obs.real, eps_obs.imag
    )

    assert not np.isclose(val_after_same_data, val_before, rtol=1.0e-6), (
        "cache did not invalidate on population-parameter change -- "
        "returned the OLD (stale) value"
    )
    assert_allclose(val_after_same_data, val_fresh, rtol=1.0e-10)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_zero_shear_matches_quad(ellip_conv):
    """g=0: no series correction at all, should reduce to the exact base
    Gaussian-kernel/population convolution, matching Quad tightly."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.28
    eps_obs = 0.1 + 0.05j

    lensed_val, quad_val = _eval_both(pop, ellip_conv, 4, 0.0 + 0.0j, eps_obs, 0.2)

    assert_allclose(lensed_val, quad_val, rtol=1.0e-4)


def test_ln_marginal_consistency():
    """The ln marginal is the log of the linear marginal."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    mset = _build_mset(pop)
    gsfsl, data = _make_series_lensed(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, 4, pop, mset, 0.3
    )

    p = gsfsl.eval_marginal(pop, data, 0.08, 0.03, 0.2, -0.1)
    lnp = gsfsl.eval_ln_marginal(pop, data, 0.08, 0.03, 0.2, -0.1)

    assert p > 0.0
    assert_allclose(lnp, math.log(p), rtol=1.0e-10)


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("order", [1, 2, 4, 8])
def test_constructs_and_evaluates_at_various_orders(ellip_conv, order):
    """trunc-order is a plain GObject uint property (minimum 1, enforced by
    the property system itself) -- just check construction and evaluation
    work cleanly across the range this module's other tests rely on."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    mset = _build_mset(pop)
    gsfsl, data = _make_series_lensed(ellip_conv, order, pop, mset, 0.3)

    val = gsfsl.eval_marginal(pop, data, 0.05, 0.0, 0.2, 0.1)

    assert np.isfinite(val)
    assert val > 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
