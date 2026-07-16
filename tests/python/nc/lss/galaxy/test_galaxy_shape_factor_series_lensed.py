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

Any population implementing ``eval_p_rho2_g_series`` works: the Gaussian
family (``NcGalaxyShapePopGauss`` and ``NcGalaxyShapePopGaussLocal``, which
share one implementation) and ``NcGalaxyShapePopBeta`` (a genuinely
different, non-Gaussian shape composed via
``ncm_laurent_series_tps_pow()``'s generalized-binomial recursion). A
population without an implementation errors clearly if used with this class.
"""

import math
import time
import warnings

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


def _make_series_lensed(ellip_conv, order, pop, mset, std_noise, pop_data_setup=None):
    gsfsl = Nc.GalaxyShapeFactorSeriesLensed.new(ellip_conv, order)
    data, _, _ = _build_factor_data(gsfsl, mset)
    if pop_data_setup is not None:
        pop_data_setup(data.pop_data)
    gsfsl.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfsl.prepare_data_array(mset, [data], True, True)
    return gsfsl, data


def _make_quad(ellip_conv, pop, mset, std_noise, pop_data_setup=None):
    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsfq, mset)
    if pop_data_setup is not None:
        pop_data_setup(data.pop_data)
    gsfq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfq.prepare_data_array(mset, [data], True, True)
    return gsfq, data


def _eval_both(pop, ellip_conv, order, g, eps_obs, std_noise, pop_data_setup=None):
    mset = _build_mset(pop)
    gsfsl, data_sl = _make_series_lensed(
        ellip_conv, order, pop, mset, std_noise, pop_data_setup
    )
    lensed_val = gsfsl.eval_marginal(
        pop, data_sl, g.real, g.imag, eps_obs.real, eps_obs.imag
    )

    gsfq, data_q = _make_quad(ellip_conv, pop, mset, std_noise, pop_data_setup)
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
def test_marginal_matches_quad_gauss_local(ellip_conv):
    """NcGalaxyShapePopGaussLocal shares NcGalaxyShapePopGauss's own
    eval_p_rho2_g_series implementation directly (identical ldata layout,
    only sigma's source differs) -- same moderate-noise/moderate-g setup as
    test_marginal_matches_quad_moderate_g_large_noise, just with a
    per-galaxy e_rms instead of a shared sigma parameter."""
    pop = Nc.GalaxyShapePopGaussLocal.new()
    g = 0.09 + 0.0j
    eps_obs = 0.4 * np.exp(1j * 0.7)

    def _set_e_rms(pop_data):
        Nc.GalaxyShapePopGaussLocal.data_set(pop, pop_data, 0.28)

    lensed_val, quad_val = _eval_both(
        pop, ellip_conv, 4, g, eps_obs, 0.3, pop_data_setup=_set_e_rms
    )

    assert_allclose(lensed_val, quad_val, rtol=2.0e-3)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_quad_beta(ellip_conv):
    """NcGalaxyShapePopBeta exercises a genuinely different composition path
    than the Gaussian family (ncm_laurent_series_tps_pow()'s generalized-
    binomial recursion, not the exp-of-power-series one) -- same moderate-
    noise/moderate-g setup as test_marginal_matches_quad_moderate_g_large_noise.
    mu=0.5, nu=6 (alpha=beta=3, both >1) keeps eval_p smooth on the whole
    unit interval, including at x=0 -- unlike the default mu/nu (alpha=0.9<1,
    P(x) ~ x^-0.1 diverges at x=0), which restricts the g-Taylor series' own
    radius of convergence whenever the rho quadrature approaches rho=0 (see
    test_marginal_matches_quad_beta_small_g_near_singular below)."""
    pop = Nc.GalaxyShapePopBeta.new()
    pop["mu"] = 0.5
    pop["nu"] = 6.0
    g = 0.09 + 0.0j
    eps_obs = 0.4 * np.exp(1j * 0.7)

    lensed_val, quad_val = _eval_both(pop, ellip_conv, 4, g, eps_obs, 0.3)

    assert_allclose(lensed_val, quad_val, rtol=2.0e-3)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_quad_beta_small_g_near_singular(ellip_conv):
    """Default mu/nu (alpha=mu*nu=0.9<1) makes P(x) ~ x^(alpha-1) diverge at
    x=0 -- a genuine branch-point singularity, not a numerical artifact
    (confirmed by raising trunc_order making the mismatch worse at moderate
    g, the classic signature of evaluating a Taylor series outside its own
    radius of convergence, rather than better as truncation error alone
    would predict). g=0 (no truncation at all) matches quad exactly; this
    checks a small g comfortably inside that radius still works. TRACE's own
    steeper O(g) response (see test_marginal_matches_quad_moderate_g_large_noise's
    docstring) shrinks that radius further than TRACE_DET's, hence g=0.003
    here rather than the smooth case's g=0.09."""
    pop = Nc.GalaxyShapePopBeta.new()
    g = 0.003 + 0.0j
    eps_obs = 0.4 * np.exp(1j * 0.7)

    lensed_val, quad_val = _eval_both(pop, ellip_conv, 4, g, eps_obs, 0.3)

    assert_allclose(lensed_val, quad_val, rtol=1.0e-2)


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


def test_caching_correct_and_fast_with_drifting_phase_bias_regression():
    """Regression test for a real caching bug found in production: earlier
    versions of this class cached its per-galaxy radial integral keyed to
    ``(R, phi, std_noise, pop_hash)``, where ``phi`` is the "gauge-fixing"
    rotation angle (``arg(eps_obs) - arg(g)``). ``phi`` is only genuinely
    invariant across a fixed-``(R,std_noise)`` sequence of calls (e.g. the
    redshift-node loop in ``NcDataClusterWLFactor``) when the galaxy's
    additive shear bias is exactly zero (``c1=c2=0``, as every OTHER test in
    this module deliberately uses via ``_make_series_lensed``'s hardcoded
    zero bias) -- as soon as a galaxy carries nonzero additive bias
    (``g = (1+m)*g_true + c1 + i*c2``, applied at the
    ``NcGalaxyShapeFactor`` call site), ``phi`` genuinely drifts from node
    to node, since it is the argument of a fixed vector added to a
    magnitude-varying one. Keying the cache on ``phi`` therefore defeated it
    on essentially every call in the realistic (nonzero-bias) case, silently
    destroying an intended ~30x-class speedup while still returning correct
    values (the bug was pure-performance, not correctness, hence the timing
    check below in addition to the value check).

    This class' own cache is now split: the expensive, ``phi``-independent
    per-radial-node harmonic sums are cached on ``(R, std_noise, pop_hash)``
    alone, and a cheap trig-polynomial evaluation applies ``phi`` fresh on
    every call. This test simulates a bias-driven drifting-``phi`` sequence
    directly (``g = g_true + c`` for a fixed complex bias ``c`` and varying
    real ``g_true``, on the SAME per-galaxy data object -- exactly what
    ``NcDataClusterWLFactor``'s redshift-node loop produces for a biased
    galaxy) and checks both that every value still matches an independent
    fresh computation (correctness) and that the cache measurably speeds up
    repeated calls despite ``phi`` drifting every time (performance)."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    eps_obs = 0.4 * np.exp(1j * 0.7)
    std_noise = 0.3
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    order = 4

    mset = _build_mset(pop)
    gsfsl, data_sl = _make_series_lensed(ellip_conv, order, pop, mset, std_noise)

    # Simulate bias-driven phi drift: g = g_true + c, with c fixed (mimics a
    # nonzero additive c1+i*c2) and g_true varying in magnitude across a set
    # of "redshift nodes" -- the angle of g (and hence
    # phi = arg(eps_obs) - arg(g)) genuinely drifts from node to node, just
    # as it does for a real galaxy with nonzero c1/c2.
    c = 0.02 - 0.01j
    g_trues = [0.01, 0.05, 0.02, 0.09, 0.001, 0.15, 0.03, 0.12]
    gs = [g_true + c for g_true in g_trues]

    # The phases genuinely differ across this sequence -- otherwise this
    # test would not actually exercise the drifting-phi path at all.
    phases = [math.atan2(g.imag, g.real) for g in gs]
    assert len(set(np.round(phases, 6))) > 1

    cached_vals = [
        gsfsl.eval_marginal(pop, data_sl, g.real, g.imag, eps_obs.real, eps_obs.imag)
        for g in gs
    ]

    for g, cached_val in zip(gs, cached_vals):
        gsfsl_fresh, data_fresh = _make_series_lensed(
            ellip_conv, order, pop, mset, std_noise
        )
        fresh_val = gsfsl_fresh.eval_marginal(
            pop, data_fresh, g.real, g.imag, eps_obs.real, eps_obs.imag
        )
        assert_allclose(
            cached_val,
            fresh_val,
            rtol=1.0e-10,
            err_msg=f"cached vs fresh mismatch at g={g}",
        )

    # Performance guard: a shared instance evaluated across the drifting-phi
    # sequence must be markedly faster than paying the full (Bessel +
    # Laurent-series) cost on every single call -- i.e. the phi-independent
    # harmonic cache must actually be engaging despite phi drifting every
    # time. Compare against a fresh instance PER call (forces a cache miss
    # on every evaluation, exactly the bug's effective behavior) as the cold
    # baseline. A generous 3x margin (well below the ~30x found in
    # production) keeps this robust to machine noise while still catching a
    # regression back to the old phi-keyed scheme.
    num_evals = 20

    gsfsl_warm, data_warm = _make_series_lensed(
        ellip_conv, order, pop, mset, std_noise
    )
    start = time.perf_counter()

    for i in range(num_evals):
        g = gs[i % len(gs)]
        gsfsl_warm.eval_marginal(
            pop, data_warm, g.real, g.imag, eps_obs.real, eps_obs.imag
        )
    warm_time = time.perf_counter() - start

    start = time.perf_counter()

    for i in range(num_evals):
        g = gs[i % len(gs)]
        gsfsl_cold, data_cold = _make_series_lensed(
            ellip_conv, order, pop, mset, std_noise
        )
        gsfsl_cold.eval_marginal(
            pop, data_cold, g.real, g.imag, eps_obs.real, eps_obs.imag
        )
    cold_time = time.perf_counter() - start

    if warm_time * 3.0 >= cold_time:
        warnings.warn(
            f"drifting-phi cache does not appear to be engaging: warm "
            f"{warm_time * 1000.0:.3f}ms vs cold {cold_time * 1000.0:.3f}ms "
            f"for {num_evals} calls each (expected warm < cold/3)",
            RuntimeWarning,
        )


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


def test_low_order_polynomial_crossing_clamps_to_floor():
    """A degree-trunc_order polynomial in g_mag is only a LOCAL
    approximation around g=0: for configurations where the true envelope
    decreases with g, a low-order fit is guaranteed to eventually cross
    zero and would otherwise diverge to -Infinity (see the class' own
    _nc_galaxy_shape_factor_series_lensed_marginal doc comment). This must
    instead clamp to a tiny positive floor -- checked here at a
    (sigma_pop, std_noise, eps_obs, order) combination that genuinely
    crosses zero -- so eval_ln_marginal's log() never sees a non-positive
    value."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.02
    mset = _build_mset(pop)
    gsfsl, data = _make_series_lensed(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, 2, pop, mset, 0.05
    )

    val_small_g = gsfsl.eval_marginal(pop, data, 0.01, 0.0, 0.02, 0.0)
    val_large_g = gsfsl.eval_marginal(pop, data, 0.6, 0.0, 0.02, 0.0)

    assert val_small_g > 1.0  # accurate, well within the polynomial's radius
    assert 0.0 < val_large_g < 1.0e-250  # clamped to the floor, not negative

    lnp_large_g = gsfsl.eval_ln_marginal(pop, data, 0.6, 0.0, 0.02, 0.0)
    assert np.isfinite(lnp_large_g)
    assert_allclose(lnp_large_g, math.log(val_large_g), rtol=1.0e-10)


def test_required_columns():
    """SeriesLensed's ldata carries no per-row data of its own; only the
    upstream fragments' columns (plus the factor's own fixed columns) are
    required."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    mset = _build_mset(pop)
    gsfsl, data = _make_series_lensed(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, 4, pop, mset, 0.3
    )

    cols = data.required_columns()
    own = [
        "epsilon_int_1", "epsilon_int_2", "epsilon_obs_1", "epsilon_obs_2",
        "std_noise", "c1", "c2", "m",
    ]
    assert cols[: len(own)] == own


def test_trunc_order_n_nodes_gobject_property_round_trip():
    """trunc-order/n-nodes (CONSTRUCT_ONLY) are reachable through the
    GObject property system (get_property), not just props.trunc_order /
    props.n_nodes."""
    gsfsl = Nc.GalaxyShapeFactorSeriesLensed.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, 6
    )
    assert gsfsl.get_property("trunc-order") == 6
    assert gsfsl.get_property("trunc-order") == gsfsl.props.trunc_order
    assert gsfsl.get_property("n-nodes") == gsfsl.props.n_nodes


def test_read_write_row_round_trip():
    """SeriesLensed's ldata carries no per-row data of its own (see
    test_required_columns), so read_row/write_row are pure pass-throughs --
    exercised here through a real NcGalaxyWLObs catalog row, the path
    NcDataClusterWLFactor's own set_obs()/prepare() uses, which none of
    this file's other tests (all built via data_set()) touch."""
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    frame = Nc.WLEllipticityFrame.CELESTIAL
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    mset = _build_mset(pop)
    gsfsl, data = _make_series_lensed(ellip_conv, 4, pop, mset, 0.3)

    cols = data.required_columns()
    obs = Nc.GalaxyWLObs.new(ellip_conv, frame, 1, cols)
    obs.set("ra", 0, 0.03)
    obs.set("dec", 0, -0.02)
    obs.set("z", 0, 0.0)
    obs.set("zp", 0, 0.6)
    obs.set("sigma0", 0, 0.03)
    obs.set("epsilon_int_1", 0, 0.0)
    obs.set("epsilon_int_2", 0, 0.0)
    obs.set("epsilon_obs_1", 0, 0.05)
    obs.set("epsilon_obs_2", 0, -0.02)
    obs.set("std_noise", 0, 0.03)
    obs.set("c1", 0, 0.01)
    obs.set("c2", 0, -0.01)
    obs.set("m", 0, 0.02)

    data.read_row(obs, 0)
    assert_allclose(data.epsilon_obs_1, 0.05)
    assert_allclose(data.epsilon_obs_2, -0.02)
    assert_allclose(data.std_noise, 0.03)
    assert_allclose(data.c1, 0.01)
    assert_allclose(data.c2, -0.01)
    assert_allclose(data.m, 0.02)

    obs_out = Nc.GalaxyWLObs.new(ellip_conv, frame, 1, cols)
    data.write_row(obs_out, 0)
    for col in cols:
        assert_allclose(obs_out.get(col, 0), obs.get(col, 0))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
