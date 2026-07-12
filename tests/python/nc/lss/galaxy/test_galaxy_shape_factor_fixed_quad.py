#!/usr/bin/env python
#
# test_galaxy_shape_factor_fixed_quad.py
#
# Thu Jul 9 2026
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Tests for the fixed-node lens-domain quadrature shape factor.

``NcGalaxyShapeFactorFixedQuad`` computes the same exact integral as
``NcGalaxyShapeFactorQuad`` (no series truncation in ``g``, no
map-linearization), but via a FIXED node count over the noise-disk /
unit-disc lens intersection in the lensed frame, instead of Quad's adaptive
Divonne cubature over a big generic box -- see that class' own docs for the
domain construction, the four branches, and the one known limitation
(narrow populations, same blind spot Quad's own predecessor had).
"""

import math

import pytest
import numpy as np
from scipy import integrate
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

_CONVS = [Nc.GalaxyWLObsEllipConv.TRACE, Nc.GalaxyWLObsEllipConv.TRACE_DET]


def _build_mset(pop):
    """One mset serving both engines: lens models shared by construction, the
    population model read only by the new calculator."""
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


def _shear_map(ellip_conv, g, chi):
    """Weak-branch forward shear map (see docs/theory/wl_ellipticity.qmd)."""
    if ellip_conv == Nc.GalaxyWLObsEllipConv.TRACE:
        denom = 1 + abs(g) ** 2 + 2 * (g * np.conj(chi)).real
        return (chi + g * (g * np.conj(chi) + 2)) / denom
    return (chi + g) / (1 + np.conj(g) * chi)


def _scipy_exact_marginal(pop, pop_data, ellip_conv, g, eps_obs, std_noise):
    """Independent truth-table oracle: direct polar integration over the
    physical disc -- no relation to FixedQuad's own lens-domain machinery,
    same oracle test_galaxy_shape_factor_quad.py uses for Quad."""

    def integrand(r, theta):
        chi = r * np.exp(1j * theta)
        p_pop = pop.eval_p(pop_data, chi.real**2 + chi.imag**2) / np.pi
        delta = eps_obs - _shear_map(ellip_conv, g, chi)
        noise = np.exp(-(delta.real**2 + delta.imag**2) / (2 * std_noise**2)) / (
            2 * np.pi * std_noise**2
        )
        return p_pop * noise * r

    result, _ = integrate.dblquad(
        integrand, 0, 2 * np.pi, 0, 1.0 - 1.0e-10, epsabs=1.0e-13, epsrel=1.0e-10
    )
    return result


def _make(ellip_conv, sigma_pop, std_noise):
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma_pop
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsffq, mset)
    gsffq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq.prepare_data_array(mset, [data], True, True)
    return gsffq, pop, data


_CASES = [
    # (g_1, g_2, eps_obs_1, eps_obs_2, sigma, std_noise)
    (0.15, 0.05, 0.42, -0.18, 0.30, 0.25),
    (0.02, 0.0, 0.1, -0.05, 0.30, 0.25),
    (0.10, -0.05, -0.2, 0.15, 0.15, 0.20),
    (0.05, 0.02, 0.05, 0.05, 0.25, 0.25),
    # production-regime crash configuration this class exists for
    (0.3, 0.0, -0.368837, 0.101348, 0.28, 0.03),
]


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("case", _CASES)
def test_marginal_matches_scipy_truth_table(ellip_conv, case):
    """FixedQuad's marginal matches an independent scipy disc integral."""
    g_1, g_2, eps_obs_1, eps_obs_2, sigma, std_noise = case
    gsffq, pop, data = _make(ellip_conv, sigma, std_noise)

    val = gsffq.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)
    exact = _scipy_exact_marginal(
        pop, data.pop_data, ellip_conv, g_1 + 1j * g_2, eps_obs_1 + 1j * eps_obs_2, std_noise
    )

    assert_allclose(val, exact, rtol=2.0e-4)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_quad(ellip_conv):
    """Cross-check against the real Quad class (not just scipy) at the
    production crash-regime configuration."""
    sigma, std_noise = 0.28, 0.03
    eps_obs_1, eps_obs_2 = -0.368837, 0.101348

    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(ellip_conv)
    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data_fq, _, _ = _build_factor_data(gsffq, mset)
    data_q, _, _ = _build_factor_data(gsfq, mset)

    for gsf, data in ((gsffq, data_fq), (gsfq, data_q)):
        gsf.data_set(
            data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
        )
        gsf.prepare_data_array(mset, [data], True, True)

    for g_mag in (0.05, 0.15, 0.3, 0.5):
        fq_val = gsffq.eval_marginal(pop, data_fq, g_mag, 0.0, eps_obs_1, eps_obs_2)
        q_val = gsfq.eval_marginal(pop, data_q, g_mag, 0.0, eps_obs_1, eps_obs_2)
        assert_allclose(fq_val, q_val, rtol=5.0e-4)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_extreme_g_stays_accurate(ellip_conv):
    """Headline regression: unlike Series/SeriesLensed, this class has no
    truncated polynomial to cross zero -- stays accurate through g=0.99,
    real and complex, at the actual crash-regime configuration. Cross-checked
    against the real Quad class rather than the scipy oracle here: scipy's
    dblquad fails to converge (hits its subdivision limit) for the most
    extreme, sharply-peaked-near-the-disc-boundary cases, while Quad's own
    adaptive Divonne cubature is already extensively validated there."""
    sigma, std_noise = 0.28, 0.03
    eps_obs_1, eps_obs_2 = -0.368837, 0.101348

    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(ellip_conv)
    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data_fq, _, _ = _build_factor_data(gsffq, mset)
    data_q, _, _ = _build_factor_data(gsfq, mset)

    for gsf, data in ((gsffq, data_fq), (gsfq, data_q)):
        gsf.data_set(
            data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
        )
        gsf.prepare_data_array(mset, [data], True, True)

    for g_1, g_2 in [
        (0.6, 0.0), (0.8, 0.0), (0.95, 0.0), (0.99, 0.0),
        (0.3, 0.2), (0.5, -0.3), (0.1, 0.6), (-0.4, 0.4),
    ]:
        fq_val = gsffq.eval_marginal(pop, data_fq, g_1, g_2, eps_obs_1, eps_obs_2)
        q_val = gsfq.eval_marginal(pop, data_q, g_1, g_2, eps_obs_1, eps_obs_2)
        assert fq_val > 0.0
        assert_allclose(fq_val, q_val, rtol=5.0e-3)


def test_deep_tail_underflows_to_floor():
    """eps_obs far enough outside the unit disc (~50 sigma at this
    std_noise) that the marginal genuinely underflows double precision --
    this floor is a defensive underflow guard now (see the class docs),
    not a branch-selection cutoff: there is no longer any domain branch
    that returns a fixed floor regardless of how deep in the tail the
    observation is (see test_no_artificial_jump_at_old_no_overlap_boundary
    for the regression this replaced)."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.01)

    val = gsffq.eval_marginal(pop, data, 0.1, 0.0, 1.5, 0.0)
    assert val <= 1.0e-290


def test_no_artificial_jump_at_old_no_overlap_boundary():
    """Regression test: this class used to have a separate "no overlap"
    domain branch (d >= R1+R2) that returned a fixed floor value regardless
    of how close d was to that threshold, causing an artificial jump in the
    marginal (and hence -2lnL) of over 1000 for an infinitesimal change in
    a single galaxy's observed ellipticity -- found via direct inspection
    at std_noise=0.03 (R=1.239999 -> R=1.240000). Fixed by growing the
    EFFECTIVE noise-disk radius the genuine-lens branch integrates over
    (not the fixed n_sigma window used for branch selection) just enough
    to guarantee real overlap depth into the unit disc, so that threshold
    is never actually reached -- see NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_
    NSIGMA_TAIL's docs in the .c file for why a full-disc quadrature
    fallback was tried first and rejected (fails silently at smaller
    std_noise, see test_adaptive_window_accurate_at_small_std_noise below).
    Checked here as: the change in -2lnL between two R values 1e-4 apart,
    straddling the old boundary, must be a small, ordinary amount
    (comparable to the change over the same step size elsewhere), not a
    four-digit jump."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.03)

    std_noise = 0.03
    old_boundary = 1.0 + 8.0 * std_noise

    def m2lnp_at(R):
        val = gsffq.eval_marginal(pop, data, 0.1, 0.0, R, 0.0)
        return -2.0 * math.log(val)

    step = 1.0e-4
    for R in (old_boundary - 2 * step, old_boundary - step, old_boundary, old_boundary + step):
        delta = abs(m2lnp_at(R + step) - m2lnp_at(R))
        assert delta < 5.0, f"jump too large ({delta}) straddling R={R}"


def test_adaptive_window_accurate_at_small_std_noise():
    """Regression test for a real flaw found in an earlier fix attempt: a
    full-disc quadrature fallback near the old no-overlap boundary looked
    fine at std_noise=0.03 (a few percent error) but failed silently at
    smaller std_noise -- off by ~2350x at std_noise=0.01, ~1e14 at
    std_noise=0.005 -- because the surviving probability there is a deep
    exp(-x^2) tail concentrated in a sliver near the disc edge, narrower
    the smaller std_noise is, which a modest uniformly-spaced full-disc
    grid has no reason to resolve. The shipped fix (growing the lens
    branch's own effective noise-disk radius instead of switching
    quadrature schemes) stays accurate here because it's the same
    already-validated Coons-patch machinery, just reaching further."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.01)

    std_noise = 0.01
    boundary = 1.0 + 8.0 * std_noise
    exact = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.1 + 0.0j, boundary + 0.0j, std_noise
    )
    val = gsffq.eval_marginal(pop, data, 0.1, 0.0, boundary, 0.0)
    assert_allclose(val, exact, rtol=0.2)


def test_branch_noise_contained_matches_scipy():
    """Small std_noise, eps_obs well inside the unit disc: the noise-disk-
    contained-in-unit-disc branch (the production regime)."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.03)

    val = gsffq.eval_marginal(pop, data, 0.2, 0.1, 0.3, 0.0)
    exact = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.2 + 0.1j, 0.3 + 0.0j, 0.03
    )
    assert_allclose(val, exact, rtol=2.0e-4)


def test_branch_unit_contained_matches_scipy():
    """Large std_noise, small eps_obs: the unit-disc-contained-in-noise-disk
    branch."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.2)

    val = gsffq.eval_marginal(pop, data, 0.1, 0.0, 0.3, 0.0)
    exact = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.1 + 0.0j, 0.3 + 0.0j, 0.2
    )
    assert_allclose(val, exact, rtol=2.0e-4)


def test_branch_genuine_lens_matches_scipy():
    """Moderate std_noise, moderate eps_obs: the genuine two-circle
    partial-overlap ("lens") branch."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.1)

    val = gsffq.eval_marginal(pop, data, 0.2, 0.0, 0.5, 0.0)
    exact = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.2 + 0.0j, 0.5 + 0.0j, 0.1
    )
    assert_allclose(val, exact, rtol=2.0e-4)


def test_cache_consistent_across_repeated_g_calls():
    """Repeated eval_marginal calls at different g on the same data object
    each match a fresh single-shot computation -- the cached domain nodes
    (g-independent) must not go stale or get reused incorrectly."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.03)

    for g_mag in (0.05, 0.3, 0.1, 0.5, 0.02):
        val = gsffq.eval_marginal(pop, data, g_mag, 0.0, 0.3, 0.0)
        exact = _scipy_exact_marginal(
            pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, g_mag + 0.0j, 0.3 + 0.0j, 0.03
        )
        assert_allclose(val, exact, rtol=2.0e-4)


def test_no_special_handling_needed_when_sigma_pop_changes():
    """Unlike SeriesLensed (whose cache must additionally track pop_hash),
    FixedQuad's domain cache depends only on (R, phi, std_noise) -- changing
    sigma_pop between two calls on the SAME data object needs no cache
    invalidation and must immediately reflect the new value. (pop_data
    itself must still be refreshed via prepare_data_array after mutating the
    model parameter -- that's the population model's own contract, not
    something specific to this class' domain cache, which is the point
    being tested here.)"""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.28
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsffq, mset)
    gsffq.data_set(
        data, 0.0, 0.0, 0.03, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq.prepare_data_array(mset, [data], True, True)

    val_1 = gsffq.eval_marginal(pop, data, 0.2, 0.0, 0.3, 0.0)
    exact_1 = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.2 + 0.0j, 0.3 + 0.0j, 0.03
    )
    assert_allclose(val_1, exact_1, rtol=2.0e-4)

    pop["sigma"] = 0.35
    gsffq.prepare_data_array(mset, [data], True, True)
    val_2 = gsffq.eval_marginal(pop, data, 0.2, 0.0, 0.3, 0.0)
    exact_2 = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.2 + 0.0j, 0.3 + 0.0j, 0.03
    )
    assert_allclose(val_2, exact_2, rtol=2.0e-4)
    assert abs(val_2 - val_1) / val_1 > 1.0e-3


def test_narrow_population_is_a_documented_limitation():
    """Documented (not silently claimed away) limitation: a fixed grid
    cannot resolve a population much narrower than its node spacing. This
    asserts the KNOWN large error exists, so a future accidental fix to this
    doesn't go unnoticed, and so this class never silently claims parity
    with Quad in a regime it was never validated for."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.02, 0.03)

    val = gsffq.eval_marginal(pop, data, 0.3, 0.0, -0.368837, 0.101348)
    exact = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.3 + 0.0j,
        -0.368837 + 0.101348j, 0.03,
    )
    rel_err = abs(val - exact) / exact
    assert rel_err > 0.5, "expected the known narrow-population blind spot"


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_ln_marginal_consistency(ellip_conv):
    """The ln marginal is the log of the linear marginal."""
    gsffq, pop, data = _make(ellip_conv, 0.3, 0.25)

    p = gsffq.eval_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)
    lnp = gsffq.eval_ln_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)

    assert p > 0.0
    assert_allclose(lnp, math.log(p), rtol=1.0e-10)


def test_zero_shear_matches_scipy():
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.1)

    val = gsffq.eval_marginal(pop, data, 0.0, 0.0, 0.1, -0.05)
    exact = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.0 + 0.0j, 0.1 - 0.05j, 0.1
    )
    assert_allclose(val, exact, rtol=2.0e-4)


def test_n_lens_forced_odd():
    """n-lens is always rounded up to the next odd number (guarantees a
    node exactly on the u=0.5 symmetry line, see the class docs)."""
    gsffq_default = Nc.GalaxyShapeFactorFixedQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    assert gsffq_default.props.n_lens % 2 == 1

    gsffq_even = Nc.GalaxyShapeFactorFixedQuad(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET, n_lens=10
    )
    assert gsffq_even.props.n_lens == 11


def _make_auto(ellip_conv, sigma_pop, std_noise, auto_lens_nodes=True):
    """Same as _make(), but constructs with auto-lens-nodes (CONSTRUCT_ONLY,
    so it must be passed to the constructor, not set afterward)."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma_pop
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad(
        ellip_conv=ellip_conv, auto_lens_nodes=auto_lens_nodes
    )
    data, _, _ = _build_factor_data(gsffq, mset)
    gsffq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq.prepare_data_array(mset, [data], True, True)
    return gsffq, pop, data


# (std_noise, R) pairs landing in the genuine-lens branch's "expensive
# middle" (R1=1, R2=8*std_noise, both discs partially overlapping) --
# same regime dev-notes/wl_fixed_quad_lens_nlens_calibration.py validated.
_EXPENSIVE_MIDDLE = [
    (std_noise, R)
    for std_noise in (0.06, 0.09, 0.12, 0.15, 0.2)
    for R in (0.1, 0.3, 0.5, 0.7)
    if abs(1.0 - 8.0 * std_noise) < R < 1.0 + 8.0 * std_noise
]


@pytest.mark.parametrize("std_noise,R", _EXPENSIVE_MIDDLE)
def test_auto_lens_nodes_matches_scipy_truth_in_expensive_middle(std_noise, R):
    """auto-lens-nodes' calibrated node count must still match the scipy
    oracle, within a tolerance consistent with the default lens-node-reltol
    (looser than the fixed default's own rtol=2e-4 -- this honestly
    reflects the calibrated tradeoff, not a tightened requirement)."""
    gsffq, pop, data = _make_auto(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, std_noise)

    val = gsffq.eval_marginal(pop, data, 0.15, 0.0, R, 0.0)
    exact = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.15 + 0.0j, R + 0.0j, std_noise
    )
    assert_allclose(val, exact, rtol=5.0e-3)


@pytest.mark.parametrize("std_noise,R", _EXPENSIVE_MIDDLE[:6])
def test_auto_lens_nodes_generalizes_across_g(std_noise, R):
    """The node count is calibrated once at a fixed g_calib (see
    _calibrate_n_lens()'s docs) but reused for every g a fit tries --
    confirm it stays accurate at g values other than the calibration
    point, not just at g_calib itself."""
    gsffq, pop, data = _make_auto(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, std_noise)

    for g_test in (0.05, 0.25, 0.4):
        val = gsffq.eval_marginal(pop, data, g_test, 0.0, R, 0.0)
        exact = _scipy_exact_marginal(
            pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET,
            g_test + 0.0j, R + 0.0j, std_noise,
        )
        assert_allclose(val, exact, rtol=5.0e-3)


def test_auto_lens_nodes_matches_fixed_default_closely():
    """auto-lens-nodes' result must be numerically close to the fixed
    n-lens=41 default on the same data -- both are approximations of the
    same exact integral, calibrated to a tolerance stricter than the gap
    between them should ever be in practice."""
    std_noise, R = 0.12, 0.3
    gsffq_fixed, pop_fixed, data_fixed = _make(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, std_noise
    )
    gsffq_auto, pop_auto, data_auto = _make_auto(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, std_noise
    )

    for g_test in (0.05, 0.15, 0.35):
        val_fixed = gsffq_fixed.eval_marginal(pop_fixed, data_fixed, g_test, 0.0, R, 0.0)
        val_auto = gsffq_auto.eval_marginal(pop_auto, data_auto, g_test, 0.0, R, 0.0)
        assert_allclose(val_auto, val_fixed, rtol=2.0e-3)


def test_auto_lens_nodes_off_by_default():
    """auto-lens-nodes defaults False -- a freshly-constructed instance
    behaves exactly like before this feature existed."""
    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    assert gsffq.props.auto_lens_nodes is False
    assert_allclose(gsffq.props.lens_node_reltol, 1.0e-4)


def test_required_columns():
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)
    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsffq, mset)

    cols = Nc.GalaxyShapeFactorData.required_columns(data)
    own = [
        "epsilon_int_1", "epsilon_int_2", "epsilon_obs_1", "epsilon_obs_2",
        "std_noise", "c1", "c2", "m",
    ]
    assert cols[: len(own)] == own


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
