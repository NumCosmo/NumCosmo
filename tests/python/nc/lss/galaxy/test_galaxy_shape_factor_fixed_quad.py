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

"""Tests for the fixed-node quadrature shape factor.

``NcGalaxyShapeFactorFixedQuad`` computes the same exact integral as
``NcGalaxyShapeFactorQuad`` (no series truncation in ``g``, no
map-linearization), via one of two FIXED quadrature schemes -- a two-panel
psi-reparametrization or a native chi_I-polar grid -- chosen per galaxy from
(eps_obs, std_noise) alone. See the class' own docs for the two schemes and
the switch criterion, and
nc_galaxy_shape_factor_fixed_quad_eval_two_panel()/_eval_chi_i_native()'s
docs for the test-only accessors used below to exercise each branch in
isolation.
"""

import math
import subprocess
import sys
from pathlib import Path

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
    physical disc -- no relation to FixedQuad's own machinery, same oracle
    test_galaxy_shape_factor_quad.py uses for Quad."""

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


def _quad_exact_marginal(
    pop, mset, ellip_conv, g_1, g_2, eps_obs_1, eps_obs_2, std_noise
):
    """Independent oracle using the real adaptive NcGalaxyShapeFactorQuad
    cubature instead of scipy's dblquad -- orders of magnitude faster (no
    Python-level integrand callback per node), while already cross-validated
    against the scipy oracle above by test_marginal_matches_quad and
    test_extreme_g_stays_accurate."""
    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfq.prepare_data_array(mset, [data], True, True)
    return gsfq.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)


def _use_chi_i_native(eps_obs, std_noise, nsigma=8.0):
    """Python mirror of the C class' own _use_chi_i_native() switch
    criterion -- used only to document/assert which branch a given test
    case lands in, never to reimplement any accuracy-relevant computation.
    Also TRUE for |eps_obs|>=1: the two-panel branch's shear_at_origin()
    is undefined there (see the class docs)."""
    if abs(eps_obs) >= 1.0:
        return True

    return (1.0 + abs(eps_obs)) <= nsigma * std_noise


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
    # production-regime configuration this class exists for
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
    production configuration."""
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


# Pinned production debug cases from devel/compare_shape_factor_methods.py,
# itself built to reproduce every FixedQuad case that failed in a real WL
# catalog run. (g_1, g_2, eps_obs_1, eps_obs_2, std_noise, alpha, beta,
# ellip_conv). The script's Gaussian truth-table cases are already covered
# by _CASES above and are not repeated here.
_REAL_DATA_DEBUG_CASES = [
    pytest.param(
        0.15,
        0.05,
        0.2183132767677307,
        -0.04060542210936546,
        0.09130239486694336,
        1.55,
        1.62,
        Nc.GalaxyWLObsEllipConv.TRACE,
        id="beta_alpha_1p55_divergent_hard_case",  # galaxy #6897, exp_007_v5-like catalog
    ),
    pytest.param(
        0.1,
        0.05,
        0.12,
        -0.08,
        0.15,
        3.0,
        2.0,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        id="beta_alpha_3p0_safe_side_of_boundary",
    ),
    pytest.param(
        -0.000282866646198319,
        -0.000353319013837193,
        0.00872796028852463,
        0.00156980438623577,
        0.010328422300517538,
        2.677998528261235,
        1.620880072781353,
        Nc.GalaxyWLObsEllipConv.TRACE,
        id="beta_alpha_2p68_tiny_shear_and_noise",
    ),
    pytest.param(
        -0.000484624097136588,
        -0.00310900893031676,
        -0.00650996621698141,
        -0.00432653771713376,
        0.015406767837703242,
        1.861177381982181,
        1.485017027354532,
        Nc.GalaxyWLObsEllipConv.TRACE,
        id="beta_alpha_1p86_tiny_shear_and_noise",
    ),
]


@pytest.mark.parametrize(
    "g_1,g_2,eps_obs_1,eps_obs_2,std_noise,alpha,beta,ellip_conv",
    _REAL_DATA_DEBUG_CASES,
)
def test_marginal_matches_quad_real_data_debug_cases(
    g_1, g_2, eps_obs_1, eps_obs_2, std_noise, alpha, beta, ellip_conv
):
    """Regression for real Beta populations (both sides of the alpha<2
    divergent-area-density boundary) at their actual production (g, eps_obs,
    std_noise), pinned from devel/compare_shape_factor_methods.py. Uses the
    Quad oracle (see _quad_exact_marginal) rather than scipy: much faster,
    with no accuracy tradeoff since it is already cross-validated against
    scipy elsewhere in this file."""
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
    exact = _quad_exact_marginal(
        pop, mset, ellip_conv, g_1, g_2, eps_obs_1, eps_obs_2, std_noise
    )

    assert math.isfinite(val)
    assert_allclose(val, exact, rtol=5.0e-3)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_scipy_truth_table_beta_alpha_below_one(ellip_conv):
    """NcGalaxyShapePopBeta with alpha<2 (the area density P_pop(r)/(2*pi*r)
    diverges at r=0): unlike SeriesLensed, whose g-Taylor series loses its
    radius of convergence here since the population stops being analytic at
    r=0, FixedQuad's direct quadrature has no series expansion to break and
    stays accurate through the singularity, via whichever branch this case
    lands in (see the class docs). alpha=1.2 sits comfortably above the
    class's own alpha>=1 floor (on r) while still being alpha<2."""
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
    against the scipy oracle directly."""
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
    angle), and the noise kernel is isotropic. Both branches' theta grids
    are offset by phi=arg(eps_obs) specifically so a coarse grid stays
    exactly rotation-covariant even for a population near-singular at
    chi_I=0 (verified directly, dev session 2026-07-29: swings of up to 96%
    of the mean across rotations at n_angular=15 without that offset, for a
    predecessor of this design)."""
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
    real and complex, at the actual production configuration. Cross-checked
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


def test_no_discontinuity_across_the_branch_switch():
    """std_noise=(1+|eps_obs|)/8 is exactly where the class' own auto-switch
    (_use_chi_i_native(), see the class docs) flips from two-panel psi to
    native chi_I at a fixed, disc-interior eps_obs: the two branches are
    genuinely different quadrature schemes, so nothing guarantees a priori
    that they agree closely right at the boundary. Checked as: the change
    in -2lnL between two std_noise values 1e-5 apart, straddling the
    switch, must be a small, ordinary amount (comparable to the change over
    the same step size elsewhere), not a large jump -- regression test for
    exactly the kind of artificial likelihood discontinuity an earlier
    (now-removed) branch design had at its own boundary."""
    eps_obs = 0.3
    switch_std_noise = (1.0 + eps_obs) / 8.0

    def m2lnp_at(std_noise):
        gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, std_noise)
        val = gsffq.eval_marginal(pop, data, 0.1, 0.0, eps_obs, 0.0)
        return -2.0 * math.log(val)

    step = 1.0e-5
    for std_noise in (
        switch_std_noise - 2 * step,
        switch_std_noise - step,
        switch_std_noise,
        switch_std_noise + step,
    ):
        delta = abs(m2lnp_at(std_noise + step) - m2lnp_at(std_noise))
        assert delta < 5.0, f"jump too large ({delta}) straddling std_noise={std_noise}"


def test_accurate_near_disc_boundary_at_small_std_noise():
    """eps_obs close to (but strictly inside) the unit disc boundary, at a
    small std_noise: the switch keeps this in the two-panel branch (the
    noise kernel is far too localized for the whole-disc condition to
    trigger), whose radial split must still resolve the thin sliver of
    surviving probability near the disc edge -- an earlier fix attempt
    (full-disc quadrature near this regime) looked fine at larger std_noise
    but failed silently here (off by orders of magnitude), because a modest
    uniformly-spaced full-disc grid has no reason to resolve that sliver."""
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    eps_obs, std_noise = 0.99, 0.01
    gsffq, pop, data = _make(ellip_conv, 0.28, std_noise)

    assert not _use_chi_i_native(eps_obs, std_noise)

    exact = _scipy_exact_marginal(
        pop,
        data.pop_data,
        ellip_conv,
        0.1 + 0.0j,
        eps_obs + 0.0j,
        std_noise,
    )
    val = gsffq.eval_marginal(pop, data, 0.1, 0.0, eps_obs, 0.0)
    assert_allclose(val, exact, rtol=5.0e-3)


def test_localized_noise_uses_two_panel_and_matches_scipy():
    """Small std_noise, eps_obs well inside the unit disc: the switch picks
    the two-panel psi branch (verified via _use_chi_i_native() directly),
    and eval_marginal() must match eval_two_panel() bit-identically (same
    code path) and scipy closely."""
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    gsffq, pop, data = _make(ellip_conv, 0.28, 0.03)
    eps_obs_1, eps_obs_2, std_noise = 0.3, 0.0, 0.03

    assert not _use_chi_i_native(eps_obs_1 + 1j * eps_obs_2, std_noise)

    val = gsffq.eval_marginal(pop, data, 0.2, 0.1, eps_obs_1, eps_obs_2)
    val_two_panel = gsffq.eval_two_panel(pop, data, 0.2, 0.1, eps_obs_1, eps_obs_2)
    exact = _scipy_exact_marginal(
        pop,
        data.pop_data,
        ellip_conv,
        0.2 + 0.1j,
        eps_obs_1 + 1j * eps_obs_2,
        std_noise,
    )
    assert val == val_two_panel
    assert_allclose(val, exact, rtol=2.0e-4)


def test_broad_noise_uses_chi_i_native_and_matches_scipy():
    """Large std_noise, small eps_obs: the switch picks the native chi_I
    branch instead, and eval_marginal() must match eval_chi_i_native()
    bit-identically and scipy closely."""
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    gsffq, pop, data = _make(ellip_conv, 0.28, 0.2)
    eps_obs_1, eps_obs_2, std_noise = 0.3, 0.0, 0.2

    assert _use_chi_i_native(eps_obs_1 + 1j * eps_obs_2, std_noise)

    val = gsffq.eval_marginal(pop, data, 0.1, 0.0, eps_obs_1, eps_obs_2)
    val_chi_i_native = gsffq.eval_chi_i_native(
        pop, data, 0.1, 0.0, eps_obs_1, eps_obs_2
    )
    exact = _scipy_exact_marginal(
        pop,
        data.pop_data,
        ellip_conv,
        0.1 + 0.0j,
        eps_obs_1 + 1j * eps_obs_2,
        std_noise,
    )
    assert val == val_chi_i_native
    assert_allclose(val, exact, rtol=2.0e-4)


def test_hybrid_switch_avoids_each_branchs_own_hard_regime():
    """Headline finding of the two-branch design (840-case sweep, see the
    class docs): each branch is excellent in its OWN regime but fails badly
    when forced into the OTHER's -- the switch exists specifically to avoid
    ever landing in either failure mode. Two cases, one per regime:

    - "hard corner" (|eps_obs| close to 1 AND std_noise broad enough that
      the whole unit disc sits within a few sigma of eps_obs): the psi
      reparametrization's own Jacobian develops a huge dynamic range near
      the disc boundary there. Forcing two-panel gives ~9% error; the
      switch (via native chi_I) gives ~0.01%.
    - "far from pop pull" (small std_noise, eps_obs well inside the disc):
      the native chi_I branch has no mechanism at all for resolving a
      narrow noise kernel wherever it happens to sit relative to chi_I=0.
      Forcing chi_i_native gives ~99% error (it returns close to zero,
      since none of its fixed nodes land anywhere near the noise peak);
      the switch (via two-panel) matches scipy to ~1e-4.
    """
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET

    pop_hard = Nc.GalaxyShapePopBeta.new()
    pop_hard["alpha"] = 1.3
    pop_hard["beta"] = 1.7
    std_noise_hard = 0.35
    eps_hard = 0.95 * np.exp(1j * 0.7)
    gsffq_hard, _, data_hard = _make_with_pop(pop_hard, ellip_conv, std_noise_hard)
    g_1, g_2 = 0.1, 0.05

    exact_hard = _scipy_exact_marginal(
        pop_hard,
        data_hard.pop_data,
        ellip_conv,
        g_1 + 1j * g_2,
        eps_hard,
        std_noise_hard,
    )
    forced_two_panel = gsffq_hard.eval_two_panel(
        pop_hard, data_hard, g_1, g_2, eps_hard.real, eps_hard.imag
    )
    switch_result = gsffq_hard.eval_marginal(
        pop_hard, data_hard, g_1, g_2, eps_hard.real, eps_hard.imag
    )
    assert abs(forced_two_panel - exact_hard) / exact_hard > 0.05
    assert abs(switch_result - exact_hard) / exact_hard < 1.0e-3

    pop_far = Nc.GalaxyShapePopGauss.new()
    pop_far["sigma"] = 0.28
    std_noise_far = 0.03
    eps_far = -0.368837 + 0.101348j
    gsffq_far, _, data_far = _make_with_pop(pop_far, ellip_conv, std_noise_far)

    exact_far = _scipy_exact_marginal(
        pop_far, data_far.pop_data, ellip_conv, 0.3 + 0.0j, eps_far, std_noise_far
    )
    forced_chi_i_native = gsffq_far.eval_chi_i_native(
        pop_far, data_far, 0.3, 0.0, eps_far.real, eps_far.imag
    )
    switch_result_far = gsffq_far.eval_marginal(
        pop_far, data_far, 0.3, 0.0, eps_far.real, eps_far.imag
    )
    assert abs(forced_chi_i_native - exact_far) / exact_far > 0.5
    assert abs(switch_result_far - exact_far) / exact_far < 1.0e-3


def test_eps_obs_outside_disc_uses_chi_i_native_and_converges():
    """Regression for a real production crash: measurement noise added to
    an unbounded TRACE distortion can legitimately push the observed
    epsilon_obs_1/epsilon_obs_2 outside the unit disc (|eps_obs|>=1), where
    the two-panel branch's shear_at_origin() is undefined (NaN for TRACE,
    no longer a genuine disc automorphism for TRACE_DET) -- found via a
    real catalog case with |eps_obs|=1.135, which silently floored to
    MIN_MARGINAL for every shear tried, corrupting -2lnL by orders of
    magnitude for that galaxy. The switch now always routes such an
    eps_obs to the native chi_I branch (see _use_chi_i_native()), which
    needs no shear_at_origin at all and converges cleanly with node count
    to the same independent scipy oracle used everywhere else in this
    file."""
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE
    eps_obs_1, eps_obs_2 = 0.392602682113647, -1.06484353542328
    std_noise = 0.15
    g_1, g_2 = 0.1, 0.05

    assert _use_chi_i_native(eps_obs_1 + 1j * eps_obs_2, std_noise)

    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.28
    mset = _build_mset(pop)
    exact = None

    for n in (15, 60):
        gsffq = Nc.GalaxyShapeFactorFixedQuad(
            ellip_conv=ellip_conv, n_radial=n, n_angular=n
        )
        data, _, _ = _build_factor_data(gsffq, mset)
        gsffq.data_set(
            data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
        )
        gsffq.prepare_data_array(mset, [data], True, True)

        if exact is None:
            exact = _scipy_exact_marginal(
                pop,
                data.pop_data,
                ellip_conv,
                g_1 + 1j * g_2,
                eps_obs_1 + 1j * eps_obs_2,
                std_noise,
            )

        val = gsffq.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)

        assert math.isfinite(val)
        assert val > 1.0e-10  # not the MIN_MARGINAL underflow floor

        if n == 60:
            assert_allclose(val, exact, rtol=1.0e-3)


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


def test_cache_switches_branch_correctly_across_calls():
    """A single `data` object reused across calls whose eps_obs/std_noise
    cross the branch switch must rebuild the domain (and pick the new
    branch) each time -- not reuse whichever branch happened to be cached
    from the previous call. Alternates two configurations landing in
    opposite branches, on the SAME data object, and checks both forced
    accessors AND eval_marginal() agree throughout."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.28, 0.03)

    localized = (0.3, 0.0, 0.03)  # eps_obs_1, eps_obs_2, std_noise -> two-panel
    broad = (0.3, 0.0, 0.2)  # -> chi_i_native

    mset = _build_mset(pop)

    for eps_obs_1, eps_obs_2, std_noise in (localized, broad, localized, broad):
        gsffq.data_set(
            data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
        )
        gsffq.prepare_data_array(mset, [data], True, True)

        eps_obs = eps_obs_1 + 1j * eps_obs_2
        expect_chi_i_native = _use_chi_i_native(eps_obs, std_noise)
        val_auto = gsffq.eval_marginal(pop, data, 0.1, 0.0, eps_obs_1, eps_obs_2)
        val_forced = (
            gsffq.eval_chi_i_native(pop, data, 0.1, 0.0, eps_obs_1, eps_obs_2)
            if expect_chi_i_native
            else gsffq.eval_two_panel(pop, data, 0.1, 0.0, eps_obs_1, eps_obs_2)
        )
        assert val_auto == val_forced


def test_no_special_handling_needed_when_sigma_pop_changes():
    """Unlike SeriesLensed (whose cache must additionally track pop_hash),
    FixedQuad's domain cache depends only on (eps_obs, std_noise) --
    changing sigma_pop between two calls on the SAME data object needs no
    cache invalidation and must immediately reflect the new value. (pop_data
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
    """Documented (not silently claimed away) limitation, shared by both
    branches: a fixed grid cannot resolve a population much narrower than
    its node spacing. This asserts the KNOWN large error exists, so a
    future accidental fix to this doesn't go unnoticed, and so this class
    never silently claims parity with Quad in a regime it was never
    validated for. sigma=0.02 (which triggered the blind spot at the old
    n-radial/n-angular default of 15) no longer does at the current default
    of 21 (rel_err drops to ~12%); 0.015 is narrow enough to still trigger
    it here."""
    gsffq, pop, data = _make(Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.015, 0.03)

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
    props.n_angular. n-radial is PER PANEL for the two-panel branch, or the
    single grid's own count for the native chi_I branch (see the class
    docs) -- this test only checks the round trip, not either branch's own
    resulting node count."""
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


def test_peek_domain_reflects_the_branch_actually_used():
    """peek_domain() exposes the SAME node array eval_marginal() builds and
    caches internally (see its docs) -- not the full per-node contribution
    to the marginal (that also needs the population density and, for
    two-panel, the inverse-map Jacobian and puncture-correction baseline,
    so summing the returned weights is explicitly NOT expected to reproduce
    eval_marginal()'s result). What IS a genuine, branch-specific invariant
    worth checking: the two-panel branch's node positions are g-INDEPENDENT
    (built once from eps_obs/std_noise alone, see the class docs), while
    the native chi_I branch's positions are the FORWARD map of a g-
    independent chi_I grid, so they must change with g."""
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.28

    for eps_obs_1, eps_obs_2, std_noise, expect_g_independent in (
        (0.3, 0.0, 0.03, True),  # two-panel regime
        (0.3, 0.0, 0.2, False),  # native chi_I regime
    ):
        gsffq, pop_i, data = _make_with_pop(pop, ellip_conv, std_noise)

        chi_l_a, w_a = gsffq.peek_domain(pop_i, data, 0.1, 0.0)
        chi_l_b, w_b = gsffq.peek_domain(pop_i, data, 0.4, -0.2)

        assert chi_l_a.to_numpy().shape[1] == 2
        assert np.all(np.isfinite(w_a.to_numpy())) and np.all(w_a.to_numpy() >= 0.0)
        assert np.all(np.isfinite(w_b.to_numpy())) and np.all(w_b.to_numpy() >= 0.0)

        positions_equal = np.allclose(chi_l_a.to_numpy(), chi_l_b.to_numpy())
        assert positions_equal == expect_g_independent


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
    alpha<2 has a genuinely divergent area density at r=0. use-marginal-
    spline's cache does not hit its usual tight tolerance here -- this test
    only asserts the result stays finite, positive, and within a coarse
    factor of the direct path, not that it matches closely."""
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
    r=0 (see nc_galaxy_shape_pop_beta.c's own docs), and (for the two-panel
    branch specifically) EVERY domain node has some g mapping to chi_I=0,
    so a box of any reasonable size is guaranteed to contain at least one
    point where ln(marginal) gets arbitrarily sharp. The adaptive autoknots
    build cannot resolve that and aborts via g_error; _build_g_spline()
    detects this population is unsafe for the adaptive path (via a
    log-log slope probe of eval_p(r) near r=0, independent of which branch
    is actually in use for this galaxy) and builds a plain FIXED knot grid
    instead, which cannot abort no matter how sharp the sampled function
    gets. This test exercises a dense sweep of the whole cached box (no
    crash anywhere) plus a same-point repeat call (must reproduce
    bit-identically, confirming the built spline -- not a fresh direct
    evaluation -- is what's actually being reused)."""
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


def test_marginal_spline_falls_back_to_direct_when_pop_density_underflows_near_origin():
    """_build_g_spline()'s own log-log slope probe (r1=1e-6, r2=1e-3) can
    itself underflow to exactly 0 at r1 for a Beta population concentrated
    tightly enough near r=1 (alpha=60 here): unlike the alpha<2 case (whose
    AREA DENSITY genuinely diverges and is handled by the fixed-knots
    grid), this is p1==0.0 itself, so the probe is not "well defined" and
    the spline is never built at all (ldata->g_spline_built stays FALSE) --
    every marginal call for this population instead falls back to
    _direct_marginal_at_g(), same as use_marginal_spline=False. Picking
    eps_obs near the population's own peak ring (~0.983) keeps the direct
    2D integral itself safely away from underflow, isolating this fallback
    branch from the separate non-finite-marginal hard error."""
    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = 60.0
    pop["beta"] = 2.0
    std_noise = 0.05
    eps_obs_1, eps_obs_2 = 0.9, 0.0

    gsffq_spline, pop_s, data_spline = _make_with_pop(
        pop, Nc.GalaxyWLObsEllipConv.TRACE_DET, std_noise, use_marginal_spline=True
    )
    gsffq_direct, pop_d, data_direct = _make_with_pop(
        pop, Nc.GalaxyWLObsEllipConv.TRACE_DET, std_noise, use_marginal_spline=False
    )
    gsffq_spline.data_set(
        data_spline,
        0.0,
        0.0,
        std_noise,
        eps_obs_1,
        eps_obs_2,
        0.0,
        Nc.WLEllipticityFrame.CELESTIAL,
    )
    gsffq_direct.data_set(
        data_direct,
        0.0,
        0.0,
        std_noise,
        eps_obs_1,
        eps_obs_2,
        0.0,
        Nc.WLEllipticityFrame.CELESTIAL,
    )

    g_1, g_2 = 0.05, 0.02
    val_spline = gsffq_spline.eval_marginal(
        pop_s, data_spline, g_1, g_2, eps_obs_1, eps_obs_2
    )
    val_direct = gsffq_direct.eval_marginal(
        pop_d, data_direct, g_1, g_2, eps_obs_1, eps_obs_2
    )

    assert math.isfinite(val_spline)
    assert val_spline > 0.0
    assert val_spline == val_direct

    # A second call at the same point must reuse the same fallback path
    # (g_spline_valid latches TRUE even though g_spline_built stays FALSE)
    # and reproduce bit-identically.
    val_spline_again = gsffq_spline.eval_marginal(
        pop_s, data_spline, g_1, g_2, eps_obs_1, eps_obs_2
    )
    assert val_spline_again == val_spline


def test_eval_marginal_floors_deep_tail_underflow():
    """A zero (underflowed) marginal is silently floored to MIN_MARGINAL
    (1e-300), not a hard error: a Beta population this concentrated
    (alpha=400) combined with an eps_obs far from its peak ring and a
    std_noise narrow enough (0.001) that even the CLOSEST quadrature node
    sits ~460 sigma away makes every node's noise value underflow to
    exactly 0 in double precision (a wider std_noise=0.01, sufficient
    before the n-radial/n-angular default was raised 15->21, no longer
    underflows: the extra nodes' closest distance stays inside the
    ~466-sigma point where exp() still returns a tiny but nonzero double).
    _direct_marginal_at_g()'s result is finite here (a sum of exact 0.0
    terms), so it takes the floor branch, not the !isfinite() g_error
    branch covered by test_eval_marginal_aborts_on_non_finite_marginal()
    below."""
    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = 400.0
    pop["beta"] = 2.0
    gsffq, pop, data = _make_with_pop(
        pop, Nc.GalaxyWLObsEllipConv.TRACE_DET, 0.001, use_marginal_spline=False
    )

    val = gsffq.eval_marginal(pop, data, -0.9, 0.9, 0.9, 0.9)

    assert math.isfinite(val)
    assert val == 1.0e-300


def test_eval_marginal_aborts_on_non_finite_marginal():
    """A genuinely non-finite (NaN) marginal is still a hard error: with
    std_noise=1e-170, sig2=std_noise**2 underflows to exactly 0.0 in
    double precision, so _noise_val()'s exp(-d2/(2*sig2))/(2*pi*sig2)
    divides by zero -- 0/0 (NaN) at d2==0, inf/0 (still 0, but 0/(2*pi*0)
    is 0/0 = NaN) elsewhere -- unlike the deep-tail case above, where
    every term is a clean, finite 0.0. Forces the chi_i_native branch via
    eval_chi_i_native() for a deterministic method. Fatal via g_error, a
    real process abort, so it is checked in a subprocess (see
    test_galaxy_shape_factor_cgf.py's test_non_gaussian_population_gate
    for the same pattern)."""
    script = (
        "import sys\n"
        "sys.path.insert(0, 'tests/python/nc/lss/galaxy')\n"
        "from numcosmo_py import Nc, Ncm\n"
        "Ncm.cfg_init()\n"
        "import test_galaxy_shape_factor_fixed_quad as m\n"
        "pop = Nc.GalaxyShapePopBeta.new()\n"
        "pop['alpha'] = 2.0\n"
        "pop['beta'] = 2.0\n"
        "gsffq, pop, data = m._make_with_pop(\n"
        "    pop, Nc.GalaxyWLObsEllipConv.TRACE_DET, 1e-170,\n"
        "    use_marginal_spline=False,\n"
        ")\n"
        "gsffq.eval_chi_i_native(pop, data, 0.0, 0.0, 0.5, 0.5)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=str(Path(__file__).resolve().parents[5]),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "non-finite" in result.stderr
    assert "_marginal_chi_i_native_debug" in result.stderr


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_eval_two_panel_aborts_on_non_finite_marginal(ellip_conv):
    """Same NaN scenario as test_eval_marginal_aborts_on_non_finite_marginal(),
    but forcing the two-panel branch directly via eval_two_panel() instead of
    the chi_i_native branch: this is the only way to reach
    _marginal_two_panel_debug() (_direct_marginal_at_g()'s TWO_PANEL
    diagnostic dump, the two-panel counterpart of _marginal_chi_i_native_debug()
    exercised above). Parametrized over both ellipticity conventions since
    the debug dump's own kernel choice (TRACE vs TRACE_DET) is a separate
    branch inside it. Fatal via g_error, checked in a subprocess as above."""
    script = (
        "import sys\n"
        "sys.path.insert(0, 'tests/python/nc/lss/galaxy')\n"
        "from numcosmo_py import Nc, Ncm\n"
        "Ncm.cfg_init()\n"
        "import test_galaxy_shape_factor_fixed_quad as m\n"
        "pop = Nc.GalaxyShapePopBeta.new()\n"
        "pop['alpha'] = 2.0\n"
        "pop['beta'] = 2.0\n"
        "gsffq, pop, data = m._make_with_pop(\n"
        f"    pop, Nc.GalaxyWLObsEllipConv.{ellip_conv.value_nick.upper().replace('-', '_')}, 1e-170,\n"
        "    use_marginal_spline=False,\n"
        ")\n"
        "gsffq.eval_two_panel(pop, data, 0.0, 0.0, 0.5, 0.5)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=str(Path(__file__).resolve().parents[5]),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "non-finite" in result.stderr
    assert "_marginal_two_panel_debug" in result.stderr
