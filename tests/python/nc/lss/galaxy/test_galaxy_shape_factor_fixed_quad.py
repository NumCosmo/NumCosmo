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
        p_pop = pop.eval_p(pop_data, abs(chi)) / (2.0 * np.pi * abs(chi))
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


def _make_with_pop(pop, ellip_conv, std_noise, **fixed_quad_kwargs):
    """Same as _make(), but takes an already-constructed pop (so the caller
    controls its parameters/type) and forwards extra kwargs (e.g.
    use_marginal_spline, spline_g_max) to the FixedQuad constructor."""
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad(ellip_conv=ellip_conv, **fixed_quad_kwargs)
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
        pop,
        data.pop_data,
        ellip_conv,
        g_1 + 1j * g_2,
        eps_obs_1 + 1j * eps_obs_2,
        std_noise,
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
def test_marginal_matches_scipy_truth_table_beta_alpha_below_one(ellip_conv):
    """NcGalaxyShapePopBeta with alpha<2 (P(x) diverges at x=0, a genuine
    singularity -- see the class doc: the r=sqrt(x)~Beta(alpha,beta)
    reparametrization's Jacobian makes this the relevant threshold, not
    alpha<1): unlike SeriesLensed, whose g-Taylor series loses its radius of
    convergence here since the population stops being analytic at x=0,
    FixedQuad's direct lens-domain quadrature has no series expansion to
    break and stays accurate through the singularity. This is a real
    fitting regime -- alpha=1.2 sits comfortably above the class's own
    alpha>=1 floor (on r) while still being alpha<2 (divergent in x)."""
    alpha, beta, std_noise = 1.2, 4.0, 0.03
    g_1, g_2 = 0.1, 0.05
    eps_obs_1, eps_obs_2 = 0.15, -0.1

    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = alpha
    pop["beta"] = beta
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsffq, mset)
    gsffq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq.prepare_data_array(mset, [data], True, True)

    val = gsffq.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)
    exact = _scipy_exact_marginal(
        pop,
        data.pop_data,
        ellip_conv,
        g_1 + 1j * g_2,
        eps_obs_1 + 1j * eps_obs_2,
        std_noise,
    )

    assert math.isfinite(val)
    assert_allclose(val, exact, rtol=2.0e-4)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_scipy_truth_table_beta_alpha_below_one_g_scan(ellip_conv):
    """Same alpha<2 Beta population as above, scanned over several shears
    against the scipy oracle directly.

    Not cross-checked against Quad: known accuracy bug for alpha<2 Beta
    populations (P(x) divergent at x=0) near g~0.18 (see
    test_galaxy_shape_factor_quad.py).
    """
    alpha, beta, std_noise = 1.2, 4.0, 0.03
    eps_obs_1, eps_obs_2 = 0.15, -0.1

    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = alpha
    pop["beta"] = beta
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsffq, mset)
    gsffq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq.prepare_data_array(mset, [data], True, True)

    for g_mag in (0.05, 0.14, 0.16, 0.18, 0.3, 0.5):
        fq_val = gsffq.eval_marginal(pop, data, g_mag, 0.0, eps_obs_1, eps_obs_2)
        exact = _scipy_exact_marginal(
            pop,
            data.pop_data,
            ellip_conv,
            g_mag + 0j,
            eps_obs_1 + 1j * eps_obs_2,
            std_noise,
        )
        assert math.isfinite(fq_val)
        assert_allclose(fq_val, exact, rtol=8.0e-3)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_is_rotation_covariant_for_narrow_beta(ellip_conv):
    """P(eps_obs*e^{ia} | g*e^{ia}) = P(eps_obs | g) exactly, for any angle
    a: the shear map is equivariant (f_{g*e^{ia}}(chi*e^{ia}) = f_g(chi)*e^{ia}
    for both ellip_conv), P_pop depends only on r=|chi_I| (never on its
    angle), and the noise kernel is isotropic -- see
    docs/theory/wl_shape_factor_history.md's derivation.

    _regen_noise_contained()/_regen_unit_contained() used to place their
    theta nodes on an ABSOLUTE grid (no dependence on phi=arg(eps_obs)),
    unlike _regen_lens()'s own local-frame construction: for a smooth
    (e.g. Gaussian) population this cost nothing (equally-spaced/
    Gauss-Legendre sampling of a smooth periodic integrand is offset-
    insensitive), but for a population near-singular at chi_I=0 it broke
    this exact identity at the class' default node count -- verified
    directly (dev session 2026-07-29): swings of up to 96% of the mean
    across rotations at n_angular=15 for alpha=1.55. Fixed by offsetting
    both branches' theta grid by phi, matching _regen_lens()."""
    alpha, beta, std_noise = 1.55, 1.62, 0.3
    gt, eps_mag, eps_arg0 = 0.15, 0.3, 0.7

    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = alpha
    pop["beta"] = beta
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsffq, mset)
    gsffq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq.prepare_data_array(mset, [data], True, True)

    baseline = None

    for a in (0.0, 0.3, 0.7, 1.2, 2.1, 3.0, -1.5):
        g = gt * np.exp(1j * a)
        eps = eps_mag * np.exp(1j * (eps_arg0 + a))
        val = gsffq.eval_marginal(pop, data, g.real, g.imag, eps.real, eps.imag)

        assert math.isfinite(val)

        if baseline is None:
            baseline = val
        else:
            assert_allclose(val, baseline, rtol=1.0e-9)


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
        (0.6, 0.0),
        (0.8, 0.0),
        (0.95, 0.0),
        (0.99, 0.0),
        (0.3, 0.2),
        (0.5, -0.3),
        (0.1, 0.6),
        (-0.4, 0.4),
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
    for R in (
        old_boundary - 2 * step,
        old_boundary - step,
        old_boundary,
        old_boundary + step,
    ):
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
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.1 + 0.0j,
        boundary + 0.0j,
        std_noise,
    )
    val = gsffq.eval_marginal(pop, data, 0.1, 0.0, boundary, 0.0)
    assert_allclose(val, exact, rtol=0.2)


def test_branch_noise_contained_matches_scipy():
    """Small std_noise, eps_obs well inside the unit disc: the noise-disk-
    contained-in-unit-disc branch (the production regime)."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.03)

    val = gsffq.eval_marginal(pop, data, 0.2, 0.1, 0.3, 0.0)
    exact = _scipy_exact_marginal(
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.2 + 0.1j,
        0.3 + 0.0j,
        0.03,
    )
    assert_allclose(val, exact, rtol=2.0e-4)


def test_branch_unit_contained_matches_scipy():
    """Large std_noise, small eps_obs: the unit-disc-contained-in-noise-disk
    branch."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.2)

    val = gsffq.eval_marginal(pop, data, 0.1, 0.0, 0.3, 0.0)
    exact = _scipy_exact_marginal(
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.1 + 0.0j,
        0.3 + 0.0j,
        0.2,
    )
    assert_allclose(val, exact, rtol=2.0e-4)


def test_branch_genuine_lens_matches_scipy():
    """Moderate std_noise, moderate eps_obs: the genuine two-circle
    partial-overlap ("lens") branch."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.1)

    val = gsffq.eval_marginal(pop, data, 0.2, 0.0, 0.5, 0.0)
    exact = _scipy_exact_marginal(
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.2 + 0.0j,
        0.5 + 0.0j,
        0.1,
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
            pop,
            data.pop_data,
            Nc.GalaxyWLObsEllipConv.TRACE_DET,
            g_mag + 0.0j,
            0.3 + 0.0j,
            0.03,
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
    gsffq.data_set(data, 0.0, 0.0, 0.03, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsffq.prepare_data_array(mset, [data], True, True)

    val_1 = gsffq.eval_marginal(pop, data, 0.2, 0.0, 0.3, 0.0)
    exact_1 = _scipy_exact_marginal(
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.2 + 0.0j,
        0.3 + 0.0j,
        0.03,
    )
    assert_allclose(val_1, exact_1, rtol=2.0e-4)

    pop["sigma"] = 0.35
    gsffq.prepare_data_array(mset, [data], True, True)
    val_2 = gsffq.eval_marginal(pop, data, 0.2, 0.0, 0.3, 0.0)
    exact_2 = _scipy_exact_marginal(
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.2 + 0.0j,
        0.3 + 0.0j,
        0.03,
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
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.3 + 0.0j,
        -0.368837 + 0.101348j,
        0.03,
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
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.0 + 0.0j,
        0.1 - 0.05j,
        0.1,
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
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.15 + 0.0j,
        R + 0.0j,
        std_noise,
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
            pop,
            data.pop_data,
            Nc.GalaxyWLObsEllipConv.TRACE_DET,
            g_test + 0.0j,
            R + 0.0j,
            std_noise,
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
        val_fixed = gsffq_fixed.eval_marginal(
            pop_fixed, data_fixed, g_test, 0.0, R, 0.0
        )
        val_auto = gsffq_auto.eval_marginal(pop_auto, data_auto, g_test, 0.0, R, 0.0)
        assert_allclose(val_auto, val_fixed, rtol=2.0e-3)


def test_auto_lens_nodes_off_by_default():
    """auto-lens-nodes defaults False -- a freshly-constructed instance
    behaves exactly like before this feature existed."""
    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    assert gsffq.props.auto_lens_nodes is False
    assert_allclose(gsffq.props.lens_node_reltol, 1.0e-4)


def test_auto_lens_nodes_skips_calibration_below_n_lens_5():
    """_calibrate_n_lens() short-circuits (returns n-lens unchanged) when
    n-lens<=5 -- too few candidates for the geometric-bracket/bisection
    search to be worth running at all. Just needs to not crash and to match
    the (necessarily coarse, n-lens=5) fixed-n_lens result exactly, since
    both paths use the same n_lens in this regime."""
    std_noise, R = 0.12, 0.3
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.28
    mset = _build_mset(pop)

    gsffq_auto = Nc.GalaxyShapeFactorFixedQuad(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET, auto_lens_nodes=True, n_lens=5
    )
    data_auto, _, _ = _build_factor_data(gsffq_auto, mset)
    gsffq_auto.data_set(
        data_auto, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq_auto.prepare_data_array(mset, [data_auto], True, True)

    gsffq_fixed = Nc.GalaxyShapeFactorFixedQuad(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET, auto_lens_nodes=False, n_lens=5
    )
    data_fixed, _, _ = _build_factor_data(gsffq_fixed, mset)
    gsffq_fixed.data_set(
        data_fixed, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq_fixed.prepare_data_array(mset, [data_fixed], True, True)

    val_auto = gsffq_auto.eval_marginal(pop, data_auto, 0.15, 0.0, R, 0.0)
    val_fixed = gsffq_fixed.eval_marginal(pop, data_fixed, 0.15, 0.0, R, 0.0)
    assert_allclose(val_auto, val_fixed, rtol=1.0e-12)


def test_auto_lens_nodes_falls_back_to_cap_when_reltol_unreachable():
    """When lens-node-reltol is tighter than any n below the n-lens cap can
    achieve (relative to the stricter 2*n_lens+1 reference), the geometric
    bracket search reaches the ceiling without finding a passing candidate
    and _calibrate_n_lens() falls back to n-lens itself -- must still
    produce a finite, correct-ish result, not crash or hang."""
    std_noise, R = 0.12, 0.3
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.28
    mset = _build_mset(pop)

    gsffq = Nc.GalaxyShapeFactorFixedQuad(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET,
        auto_lens_nodes=True,
        lens_node_reltol=1.0e-14,
    )
    data, _, _ = _build_factor_data(gsffq, mset)
    gsffq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq.prepare_data_array(mset, [data], True, True)

    val = gsffq.eval_marginal(pop, data, 0.15, 0.0, R, 0.0)
    exact = _scipy_exact_marginal(
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        0.15 + 0.0j,
        R + 0.0j,
        std_noise,
    )
    assert np.isfinite(val)
    assert_allclose(val, exact, rtol=5.0e-3)


def test_required_columns():
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)
    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsffq, mset)

    cols = Nc.GalaxyShapeFactorData.required_columns(data)
    own = [
        "epsilon_int_1",
        "epsilon_int_2",
        "epsilon_obs_1",
        "epsilon_obs_2",
        "std_noise",
        "c1",
        "c2",
        "m",
    ]
    assert cols[: len(own)] == own


def test_n_radial_n_angular_gobject_property_round_trip():
    """n-radial/n-angular (CONSTRUCT_ONLY) are reachable through the
    GObject property system (get_property), not just props.n_radial /
    props.n_angular."""
    gsffq = Nc.GalaxyShapeFactorFixedQuad(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET, n_radial=17, n_angular=9
    )
    assert gsffq.get_property("n-radial") == 17
    assert gsffq.get_property("n-angular") == 9
    assert gsffq.get_property("n-radial") == gsffq.props.n_radial
    assert gsffq.get_property("n-angular") == gsffq.props.n_angular


def test_read_write_row_round_trip():
    """FixedQuad's ldata carries no per-row data of its own (see
    test_required_columns), so read_row/write_row are pure pass-throughs --
    exercised here through a real NcGalaxyWLObs catalog row, the path
    NcDataClusterWLFactor's own set_obs()/prepare() uses, which none of this
    file's other tests (all built via data_set()) touch."""
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    frame = Nc.WLEllipticityFrame.CELESTIAL
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    mset = _build_mset(pop)
    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsffq, mset)

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


def test_use_marginal_spline_off_by_default():
    """use-marginal-spline defaults False -- a freshly-constructed instance
    behaves exactly like before this feature existed."""
    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    assert gsffq.props.use_marginal_spline is False
    assert_allclose(gsffq.props.spline_g_max, 0.3)
    assert_allclose(gsffq.props.spline_rel_err, 1.0e-4)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_use_marginal_spline_matches_direct_in_box(ellip_conv):
    """use-marginal-spline's cached ln(marginal) surface must closely match
    the direct (uncached) path everywhere inside its own box, for
    populations without a genuine density singularity (Gauss, and a
    Beta with alpha>=2 -- see the alpha<2 caveat test below for the one
    known-hard exception)."""
    rng = np.random.default_rng(11)
    std_noise, eps_obs_1, eps_obs_2, g_max = 0.25, 0.42, -0.18, 0.15

    for pop_direct, pop_spline in (
        (Nc.GalaxyShapePopGauss.new(), Nc.GalaxyShapePopGauss.new()),
        (Nc.GalaxyShapePopBeta.new(), Nc.GalaxyShapePopBeta.new()),
    ):
        if isinstance(pop_direct, Nc.GalaxyShapePopGauss):
            pop_direct["sigma"] = pop_spline["sigma"] = 0.3
        else:
            pop_direct["alpha"] = pop_spline["alpha"] = 3.0
            pop_direct["beta"] = pop_spline["beta"] = 6.0

        gsffq_direct, pop_d, data_direct = _make_with_pop(
            pop_direct, ellip_conv, std_noise
        )
        gsffq_spline, pop_s, data_spline = _make_with_pop(
            pop_spline,
            ellip_conv,
            std_noise,
            use_marginal_spline=True,
            spline_g_max=g_max,
        )

        for _ in range(30):
            g_1 = rng.uniform(-g_max * 0.95, g_max * 0.95)
            g_2 = rng.uniform(-g_max * 0.95, g_max * 0.95)
            val_direct = gsffq_direct.eval_marginal(
                pop_d, data_direct, g_1, g_2, eps_obs_1, eps_obs_2
            )
            val_spline = gsffq_spline.eval_marginal(
                pop_s, data_spline, g_1, g_2, eps_obs_1, eps_obs_2
            )
            assert_allclose(val_spline, val_direct, rtol=5.0e-3)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_use_marginal_spline_falls_back_exactly_outside_box(ellip_conv):
    """g outside [-spline-g-max,spline-g-max]^2 always uses the exact direct
    computation (never extrapolated from the cached surface), so it must
    match the non-spline path to machine precision, not just to the
    cache's own tolerance."""
    std_noise, eps_obs_1, eps_obs_2, g_max = 0.25, 0.42, -0.18, 0.15
    pop_direct = Nc.GalaxyShapePopGauss.new()
    pop_spline = Nc.GalaxyShapePopGauss.new()
    pop_direct["sigma"] = pop_spline["sigma"] = 0.3

    gsffq_direct, pop_d, data_direct = _make_with_pop(pop_direct, ellip_conv, std_noise)
    gsffq_spline, pop_s, data_spline = _make_with_pop(
        pop_spline, ellip_conv, std_noise, use_marginal_spline=True, spline_g_max=g_max
    )

    for g_1, g_2 in ((0.2, 0.0), (0.0, -0.25), (0.29, 0.29)):
        val_direct = gsffq_direct.eval_marginal(
            pop_d, data_direct, g_1, g_2, eps_obs_1, eps_obs_2
        )
        val_spline = gsffq_spline.eval_marginal(
            pop_s, data_spline, g_1, g_2, eps_obs_1, eps_obs_2
        )
        assert_allclose(val_spline, val_direct, rtol=1.0e-12)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_use_marginal_spline_invalidates_on_pop_pkey_change(ellip_conv):
    """Changing pop's own parameters (and re-preparing, as any caller must
    do whenever pop's parameters change -- see this class' own prepare
    protocol) after the g-spline was already built must be picked up on
    the next call, not silently served from a stale cache built under the
    old parameters."""
    std_noise, eps_obs_1, eps_obs_2, g_max = 0.25, 0.42, -0.18, 0.15
    g_1, g_2 = 0.1, -0.05

    pop_spline = Nc.GalaxyShapePopGauss.new()
    pop_spline["sigma"] = 0.3
    mset_spline = _build_mset(pop_spline)
    gsffq_spline = Nc.GalaxyShapeFactorFixedQuad(
        ellip_conv=ellip_conv, use_marginal_spline=True, spline_g_max=g_max
    )
    data_spline, _, _ = _build_factor_data(gsffq_spline, mset_spline)
    gsffq_spline.data_set(
        data_spline, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsffq_spline.prepare_data_array(mset_spline, [data_spline], True, True)

    # First call builds the g-spline under sigma=0.3.
    gsffq_spline.eval_marginal(pop_spline, data_spline, g_1, g_2, eps_obs_1, eps_obs_2)

    # Changing sigma bumps pop's own pkey; re-preparing must invalidate the
    # already-built g-spline rather than silently reuse it.
    pop_spline["sigma"] = 0.15
    gsffq_spline.prepare_data_array(mset_spline, [data_spline], True, True)
    val_spline = gsffq_spline.eval_marginal(
        pop_spline, data_spline, g_1, g_2, eps_obs_1, eps_obs_2
    )

    pop_direct = Nc.GalaxyShapePopGauss.new()
    pop_direct["sigma"] = 0.15
    gsffq_direct, pop_d, data_direct = _make_with_pop(pop_direct, ellip_conv, std_noise)
    val_direct = gsffq_direct.eval_marginal(
        pop_d, data_direct, g_1, g_2, eps_obs_1, eps_obs_2
    )

    assert_allclose(val_spline, val_direct, rtol=5.0e-3)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_use_marginal_spline_beta_alpha_below_one_caveat(ellip_conv):
    """Known caveat, documented in the class docs: a Beta population with
    alpha<2 has a genuinely divergent density at x=0 (the same regime
    test_marginal_matches_scipy_truth_table_beta_alpha_below_one_g_scan
    documents as a known accuracy issue for NcGalaxyShapeFactorQuad itself,
    near g~0.18). use-marginal-spline's cache does not hit its usual tight
    tolerance here -- this test only asserts the result stays finite,
    positive, and within a coarse factor of the direct path, not that it
    matches closely."""
    std_noise, eps_obs_1, eps_obs_2, g_max = 0.03, 0.15, -0.1, 0.15
    pop_direct = Nc.GalaxyShapePopBeta.new()
    pop_spline = Nc.GalaxyShapePopBeta.new()
    pop_direct["alpha"] = pop_spline["alpha"] = 1.2
    pop_direct["beta"] = pop_spline["beta"] = 4.0

    gsffq_direct, pop_d, data_direct = _make_with_pop(pop_direct, ellip_conv, std_noise)
    gsffq_spline, pop_s, data_spline = _make_with_pop(
        pop_spline, ellip_conv, std_noise, use_marginal_spline=True, spline_g_max=g_max
    )

    g_1, g_2 = 0.09, -0.04
    val_direct = gsffq_direct.eval_marginal(
        pop_d, data_direct, g_1, g_2, eps_obs_1, eps_obs_2
    )
    val_spline = gsffq_spline.eval_marginal(
        pop_s, data_spline, g_1, g_2, eps_obs_1, eps_obs_2
    )
    assert math.isfinite(val_spline)
    assert val_spline > 0.0
    assert_allclose(val_spline, val_direct, rtol=1.0)  # coarse: within 2x


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_use_marginal_spline_beta_alpha_below_two_never_crashes_and_caches(ellip_conv):
    """Regression for the production crash this fix addresses: a Beta
    population with 1<alpha<2 has a genuinely divergent AREA DENSITY at
    r=0 (see nc_galaxy_shape_pop_beta.c's own docs), and EVERY domain node
    has some g mapping to chi_I=0 (see _build_g_spline()'s docs) -- so a
    box of any reasonable size is guaranteed to contain at least one point
    where ln(marginal) gets arbitrarily sharp. The adaptive autoknots
    build cannot resolve that and aborts via g_error; _build_g_spline()
    detects this population is unsafe for the adaptive path (via the
    log-log slope probe of eval_p(r) near r=0) and builds a plain FIXED
    knot grid instead (_build_g_spline_fixed_knots()), which cannot abort
    no matter how sharp the sampled function gets. This test exercises a
    dense sweep of the whole cached box (no crash anywhere) plus a
    same-point repeat call (must reproduce bit-identically, confirming the
    built spline -- not a fresh direct evaluation -- is what's actually
    being reused)."""
    alpha, beta, std_noise, g_max = 1.55, 1.62, 0.3, 0.3
    eps_obs_1, eps_obs_2 = 0.05, -0.02

    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = alpha
    pop["beta"] = beta
    gsffq, pop, data = _make_with_pop(
        pop, ellip_conv, std_noise, use_marginal_spline=True, spline_g_max=g_max
    )

    n = 25
    grid = np.linspace(-g_max, g_max, n)
    for g_1 in grid:
        for g_2 in grid:
            val = gsffq.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)
            assert math.isfinite(val)
            assert val > 0.0

    val_a = gsffq.eval_marginal(pop, data, 0.1, 0.05, eps_obs_1, eps_obs_2)
    val_b = gsffq.eval_marginal(pop, data, 0.1, 0.05, eps_obs_1, eps_obs_2)
    assert val_a == val_b


def test_spline_g_max_rel_err_gobject_property_round_trip():
    """spline-g-max/spline-rel-err (CONSTRUCT_ONLY) are reachable through
    the GObject property system, not just props.spline_g_max /
    props.spline_rel_err."""
    gsffq = Nc.GalaxyShapeFactorFixedQuad(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET,
        spline_g_max=0.2,
        spline_rel_err=1.0e-3,
    )
    assert_allclose(gsffq.get_property("spline-g-max"), 0.2)
    assert_allclose(gsffq.get_property("spline-rel-err"), 1.0e-3)
    assert_allclose(gsffq.get_property("spline-g-max"), gsffq.props.spline_g_max)
    assert_allclose(gsffq.get_property("spline-rel-err"), gsffq.props.spline_rel_err)


def test_pop_correction_off_by_default():
    """use-pop-correction defaults False -- a freshly-constructed instance
    behaves exactly like before this feature existed."""
    gsffq = Nc.GalaxyShapeFactorFixedQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    assert gsffq.props.use_pop_correction is False
    assert_allclose(gsffq.props.pop_correction_eps1, 0.15)
    assert_allclose(gsffq.props.pop_correction_eps2, 0.35)
    assert gsffq.props.pop_correction_n_radial == 8
    assert gsffq.props.pop_correction_n_angular == 16


def test_pop_correction_gobject_property_round_trip():
    """pop-correction-* (CONSTRUCT_ONLY) are reachable through the GObject
    property system, not just props.pop_correction_*."""
    gsffq = Nc.GalaxyShapeFactorFixedQuad(
        ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET,
        use_pop_correction=True,
        pop_correction_eps1=0.1,
        pop_correction_eps2=0.4,
        pop_correction_n_radial=6,
        pop_correction_n_angular=12,
    )
    assert gsffq.get_property("use-pop-correction") is True
    assert_allclose(gsffq.get_property("pop-correction-eps1"), 0.1)
    assert_allclose(gsffq.get_property("pop-correction-eps2"), 0.4)
    assert gsffq.get_property("pop-correction-n-radial") == 6
    assert gsffq.get_property("pop-correction-n-angular") == 12


def test_pop_correction_reduces_error_for_real_hard_case():
    """Real hard case found via a full-catalog scan (dev session 2026-07-29,
    galaxy #6897 of a production exp_007_v5-like catalog): a Beta
    alpha=1.55 population where eps_obs happens to sit close to f_g(0)
    (the chi_I=0 image) inside the noise-contained branch's densely-sampled
    region. The pointwise P_pop(r)/(2*pi*r) divergence there -- not a
    general narrow-population problem, see
    test_narrow_population_is_a_documented_limitation for that different
    case -- is what a coarse fixed grid (n_radial=n_angular=15, the class
    default) cannot resolve: 18.5% off the independent scipy oracle here,
    vs. a median error of ~1e-9 across the other 10705 galaxies in that same
    real catalog at the same resolution. use-pop-correction brings this one
    hard case back down near the class' own usual accuracy without adding
    more nodes to the (still cached, still cheap) main grid."""
    alpha, beta = 1.55, 1.62
    std_noise = 0.09130239486694336
    g_1, g_2 = 0.15, 0.05
    eps_obs_1, eps_obs_2 = 0.2183132767677307, -0.04060542210936546
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE

    pop_plain = Nc.GalaxyShapePopBeta.new()
    pop_plain["alpha"] = alpha
    pop_plain["beta"] = beta
    pop_corr = Nc.GalaxyShapePopBeta.new()
    pop_corr["alpha"] = alpha
    pop_corr["beta"] = beta

    gsffq_plain, pop_p, data_plain = _make_with_pop(pop_plain, ellip_conv, std_noise)
    gsffq_corr, pop_c, data_corr = _make_with_pop(
        pop_corr, ellip_conv, std_noise, use_pop_correction=True
    )

    exact = _scipy_exact_marginal(
        pop_p,
        data_plain.pop_data,
        ellip_conv,
        g_1 + 1j * g_2,
        eps_obs_1 + 1j * eps_obs_2,
        std_noise,
    )
    val_plain = gsffq_plain.eval_marginal(
        pop_p, data_plain, g_1, g_2, eps_obs_1, eps_obs_2
    )
    val_corr = gsffq_corr.eval_marginal(
        pop_c, data_corr, g_1, g_2, eps_obs_1, eps_obs_2
    )

    rel_err_plain = abs(val_plain - exact) / exact
    rel_err_corr = abs(val_corr - exact) / exact

    assert rel_err_plain > 0.1, "expected the known hard-case error without correction"
    assert_allclose(val_corr, exact, rtol=5.0e-3)
    assert rel_err_corr < rel_err_plain / 10.0


def test_pop_correction_improves_g_spline_accuracy_for_real_hard_case():
    """use-pop-correction was motivated by exactly this: use-marginal-spline's
    fixed-knot fallback for alpha<2 Beta populations (see
    NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_UNSAFE_SPLINE_N_KNOTS's own docs, which
    document a bounded but real worst-case interpolation error) samples
    _direct_marginal_at_g() at each knot -- the SAME shared function
    use-pop-correction fixes, so enabling both together must interpolate a
    meaningfully more accurate underlying surface, with no extra wiring.
    Checked directly against the independent scipy oracle at the same real
    hard-case configuration as test_pop_correction_reduces_error_for_real_
    hard_case."""
    alpha, beta = 1.55, 1.62
    std_noise = 0.09130239486694336
    g_1, g_2 = 0.15, 0.05
    eps_obs_1, eps_obs_2 = 0.2183132767677307, -0.04060542210936546
    g_max = 0.3
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE

    pop_old = Nc.GalaxyShapePopBeta.new()
    pop_old["alpha"] = alpha
    pop_old["beta"] = beta
    gsffq_old, pop_o, data_old = _make_with_pop(
        pop_old, ellip_conv, std_noise, use_marginal_spline=True, spline_g_max=g_max
    )

    pop_new = Nc.GalaxyShapePopBeta.new()
    pop_new["alpha"] = alpha
    pop_new["beta"] = beta
    gsffq_new, pop_n, data_new = _make_with_pop(
        pop_new,
        ellip_conv,
        std_noise,
        use_marginal_spline=True,
        spline_g_max=g_max,
        use_pop_correction=True,
    )

    exact = _scipy_exact_marginal(
        pop_o,
        data_old.pop_data,
        ellip_conv,
        g_1 + 1j * g_2,
        eps_obs_1 + 1j * eps_obs_2,
        std_noise,
    )
    val_old = gsffq_old.eval_marginal(pop_o, data_old, g_1, g_2, eps_obs_1, eps_obs_2)
    val_new = gsffq_new.eval_marginal(pop_n, data_new, g_1, g_2, eps_obs_1, eps_obs_2)

    ln_err_old = abs(math.log(val_old) - math.log(exact))
    ln_err_new = abs(math.log(val_new) - math.log(exact))

    assert math.isfinite(val_new)
    assert val_new > 0.0
    assert ln_err_new < ln_err_old / 5.0


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_pop_correction_is_rotation_covariant(ellip_conv):
    """Same exact identity as test_marginal_is_rotation_covariant_for_narrow_beta
    above, now with use-pop-correction enabled: the outer taper depends only
    on r_i=|chi_I(chi_L_i,g)|, already rotation-covariant via the domain's
    own phi offset, and the inner correction integrates in native chi_I
    polar coordinates with no absolute angular reference of its own (unlike
    the bug that test documents, the population term here has no
    theta-dependence at all before the shear map is applied) -- verified
    numerically in a dev session prototype before shipping this, checked
    here against the real C implementation."""
    alpha, beta, std_noise = 1.55, 1.62, 0.09130239486694336
    gt, eps_mag, eps_arg0 = 0.15, 0.2220573959987452, -0.18398692

    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = alpha
    pop["beta"] = beta

    gsffq, pop, data = _make_with_pop(
        pop,
        ellip_conv,
        std_noise,
        use_pop_correction=True,
        pop_correction_eps1=0.15,
        pop_correction_eps2=0.35,
    )

    baseline = None

    for a in (0.0, 0.3, 0.7, 1.2, 2.1, 3.0, -1.5):
        g = gt * np.exp(1j * a)
        eps = eps_mag * np.exp(1j * (eps_arg0 + a))
        val = gsffq.eval_marginal(pop, data, g.real, g.imag, eps.real, eps.imag)

        assert math.isfinite(val)

        if baseline is None:
            baseline = val
        else:
            assert_allclose(val, baseline, rtol=1.0e-7)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_pop_correction_no_regression_for_smooth_population(ellip_conv):
    """A Gaussian population has no adjustable parameter that makes
    P_pop(r)/(2*pi*r) diverge at r=0 (unlike Beta's alpha) -- there is
    nothing for use-pop-correction to fix here, and it must not make an
    already-accurate result meaningfully worse."""
    sigma, std_noise = 0.28, 0.03
    eps_obs_1, eps_obs_2 = -0.368837, 0.101348
    g_1, g_2 = 0.2, 0.1

    pop_plain = Nc.GalaxyShapePopGauss.new()
    pop_plain["sigma"] = sigma
    pop_corr = Nc.GalaxyShapePopGauss.new()
    pop_corr["sigma"] = sigma

    gsffq_plain, pop_p, data_plain = _make_with_pop(pop_plain, ellip_conv, std_noise)
    gsffq_corr, pop_c, data_corr = _make_with_pop(
        pop_corr, ellip_conv, std_noise, use_pop_correction=True
    )

    val_plain = gsffq_plain.eval_marginal(
        pop_p, data_plain, g_1, g_2, eps_obs_1, eps_obs_2
    )
    val_corr = gsffq_corr.eval_marginal(
        pop_c, data_corr, g_1, g_2, eps_obs_1, eps_obs_2
    )

    assert_allclose(val_corr, val_plain, rtol=5.0e-3)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
