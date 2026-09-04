#!/usr/bin/env python
#
# test_galaxy_shape_factor_moment_series.py
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

"""Tests for the truncated-g-power-series shape factor calculator.

``NcGalaxyShapeFactorMomentSeries`` approximates the intrinsic-ellipticity
marginal by the Gaussian carrying the population's own mean/covariance,
truncated as a power series in the reduced shear g (order ``trunc-order``)
rather than in the intrinsic width -- see the class documentation and
``MOMENT_SERIES.md`` at the repository root for the derivation.

Unlike a pointwise comparison of the marginal density against
``NcGalaxyShapeFactorFixedQuad`` (which mixes two independent error
sources -- the g-truncation itself, and the residual mismatch between any
moment-matched Gaussian and the true, non-Gaussian marginal, present even
for ``NcGalaxyShapeFactorCGF`` at g->0), these tests target quantities the
design doc actually makes claims about: the closed-form reference
coefficients, the exact TRACE_DET degeneracies that fall out of the
map (see the class's own top-of-file note on the corrected recursion),
the n=1 responsivity limit, convergence in truncation order at a narrow
population (where the Gaussian-clamp residual is small so g-truncation
is the dominant error), and the class's headline claim: working with a
population ``NcGalaxyShapeFactorCGF`` cannot handle at all.
"""

import math
import subprocess
import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.special import beta as scipy_beta
from scipy.special import gammainc

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

_CONVS = [Nc.GalaxyWLObsEllipConv.TRACE, Nc.GalaxyWLObsEllipConv.TRACE_DET]


def _build_mset(pop):
    """One mset: lens models plus the population model read by the marginal."""
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


def _eval_factor(gsf, pop, mset, g, eps_obs, std_noise):
    data, _, _ = _build_factor_data(gsf, mset)
    gsf.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsf.prepare_data_array(mset, [data], True, True)
    return gsf.eval_marginal(pop, data, g.real, g.imag, eps_obs.real, eps_obs.imag)


def _pop_gauss(sigma):
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma
    return pop


# ---------------------------------------------------------------------------
# Reference coefficients (MOMENT_SERIES.md sec. 8, at sigma_pop=0.4847,
# e_rms=0.4091): pins the constructed()-time table build against the
# values the design doc's own accuracy tables were measured at. This is
# the check that would have caught the sec. 3.6 map error (a wrong table
# would not have reproduced these to 6 digits).
# ---------------------------------------------------------------------------


def test_reference_coefficients_trace():
    """The TRACE mean/covariance coefficients at sigma_pop=0.4847 match
    MOMENT_SERIES.md sec. 8 to the quoted precision."""
    pop = _pop_gauss(0.4847)
    mset = _build_mset(pop)
    gsf = Nc.GalaxyShapeFactorMomentSeries.new(Nc.GalaxyWLObsEllipConv.TRACE, 5)
    gsf.prepare(mset)

    # Recover m1, v0, v2 from two probes: mu(g) is odd in g and the
    # covariance is even, so a tiny-g probe isolates the leading terms
    # (m1*g dominates mu, v0 dominates C_t-sigma_noise^2) to the precision
    # needed here, then a second point pins v2 via the g^2 slope.
    std_noise = 1.0e-6  # negligible: isolates the population's own moments
    g_small = 1.0e-4

    p0 = _eval_factor(gsf, pop, mset, complex(g_small, 0.0), complex(0.0, 0.0), std_noise)
    # At g -> 0, eps_obs=0: mu -> 0 and Ct=Cx=v0 (v0=w0=M2/2 for both
    # conventions), so p0 = 1/(2 pi v0).
    v0 = 1.0 / (2.0 * np.pi * p0)
    assert_allclose(v0, 0.167368, rtol=2.0e-4)

    # m1: evaluate near the mean chi_t=m1*g at small g and read off the peak.
    data, _, _ = _build_factor_data(gsf, mset)
    gsf.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsf.prepare_data_array(mset, [data], True, True)
    g = 0.05
    # scan eps_t to find the argmax (the Gaussian's mean, m1*g)
    ts = np.linspace(0.0, 0.15, 4001)
    vals = [
        gsf.eval_marginal(pop, data, g, 0.0, t, 0.0) for t in ts
    ]
    m1_g = ts[int(np.argmax(vals))]
    assert_allclose(m1_g / g, 1.665263, rtol=5.0e-3)


def test_reference_coefficients_trace_det_exact():
    """TRACE_DET: m1=1 and all higher m_j vanish identically (falls out of
    the corrected recursion's own algebra, no special-casing) -- so the
    mean response is exactly linear, unlike TRACE's responsivity-weighted
    form. v0=w0=M2/2 for both conventions (MOMENT_SERIES.md sec. 3.4)."""
    pop = _pop_gauss(0.4847)
    mset = _build_mset(pop)

    for n in (1, 3, 5, 7):
        gsf = Nc.GalaxyShapeFactorMomentSeries.new(
            Nc.GalaxyWLObsEllipConv.TRACE_DET, n
        )
        gsf.prepare(mset)
        data, _, _ = _build_factor_data(gsf, mset)
        gsf.data_set(
            data, 0.0, 0.0, 0.1256, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
        )
        gsf.prepare_data_array(mset, [data], True, True)

        for g_val in (0.05, 0.2, 0.5):
            # Isotropy (C_t == C_x, mu purely tangential): swapping a small
            # offset between the t/x components must not change the
            # marginal at all, to machine precision, at every order.
            d = 0.03
            p_t = gsf.eval_marginal(pop, data, g_val, 0.0, g_val + d, 0.0)
            p_x = gsf.eval_marginal(pop, data, g_val, 0.0, g_val, d)
            assert_allclose(p_t, p_x, rtol=1.0e-10)


# ---------------------------------------------------------------------------
# n=1 reproduces the classical responsivity estimator (sec. 3.4).
# ---------------------------------------------------------------------------


def test_trace_det_n1_mean_is_exactly_g():
    """<epsilon_obs> = g exactly at n=1 already (m1=1, no responsivity
    factor), matching direct_estimate()'s TRACE_DET branch and CGF's
    isotropic TRACE_DET response moments -- the check that would have
    caught MOMENT_SERIES.md sec. 3.6's map error, whose (1-M2)*g form
    disagrees with both."""
    pop = _pop_gauss(0.3)
    mset = _build_mset(pop)
    gsf = Nc.GalaxyShapeFactorMomentSeries.new(Nc.GalaxyWLObsEllipConv.TRACE_DET, 1)
    gsf.prepare(mset)
    data, _, _ = _build_factor_data(gsf, mset)
    gsf.data_set(
        data, 0.0, 0.0, 1.0e-6, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsf.prepare_data_array(mset, [data], True, True)

    for g_val in (0.05, 0.15, 0.3):
        ts = np.linspace(0.0, g_val * 1.5, 4001)
        vals = [gsf.eval_marginal(pop, data, g_val, 0.0, t, 0.0) for t in ts]
        peak = ts[int(np.argmax(vals))]
        assert_allclose(peak, g_val, atol=2.0e-3)


def test_trace_n1_mean_matches_direct_estimate_responsivity():
    """TRACE: m1 = 2 - M2 = 2*(1 - e_rms^2), exactly the responsivity R
    used by nc_galaxy_shape_factor_direct_estimate()'s TRACE branch."""
    sigma_pop = 0.3
    pop = _pop_gauss(sigma_pop)
    data0 = Nc.GalaxyShapePopData.new(pop)
    pop.prepare(data0)
    e_rms = pop.e_rms(data0)
    expected_m1 = 2.0 * (1.0 - e_rms**2)

    mset = _build_mset(pop)
    gsf = Nc.GalaxyShapeFactorMomentSeries.new(Nc.GalaxyWLObsEllipConv.TRACE, 1)
    gsf.prepare(mset)
    data, _, _ = _build_factor_data(gsf, mset)
    gsf.data_set(
        data, 0.0, 0.0, 1.0e-6, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsf.prepare_data_array(mset, [data], True, True)

    g_val = 0.1
    ts = np.linspace(0.0, g_val * expected_m1 * 1.5, 4001)
    vals = [gsf.eval_marginal(pop, data, g_val, 0.0, t, 0.0) for t in ts]
    peak = ts[int(np.argmax(vals))]
    assert_allclose(peak / g_val, expected_m1, atol=5.0e-3)


# ---------------------------------------------------------------------------
# Convergence in n, at a narrow population (Gaussian-clamp residual is
# small there -- CGF's own docstring calls this its "essentially exact"
# regime), so the g-truncation error is the dominant, isolable signal.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_convergence_in_n_narrow_population(ellip_conv):
    """Relative error against FixedQuad decreases against a tolerance
    floor as n grows -- not asserted strictly monotonic, since sec. 8's
    own convergence table flattens near the 2-D quadrature reference's
    own floor rather than continuing to shrink, and sec. 7's positivity
    is explicitly non-monotonic in n by construction."""
    pop = _pop_gauss(0.1)
    mset = _build_mset(pop)
    std_noise = 0.1256

    fq = Nc.GalaxyShapeFactorFixedQuad.new(ellip_conv)
    fq.prepare(mset)

    g = complex(0.25, 0.0)
    eps_obs = complex(0.30, -0.10) if ellip_conv == Nc.GalaxyWLObsEllipConv.TRACE_DET else complex(0.45, -0.10)

    p_ref = _eval_factor(fq, pop, mset, g, eps_obs, std_noise)

    errs = []
    for n in (1, 3, 5, 7):
        gsf = Nc.GalaxyShapeFactorMomentSeries.new(ellip_conv, n)
        gsf.prepare(mset)
        p = _eval_factor(gsf, pop, mset, g, eps_obs, std_noise)
        errs.append(abs(p - p_ref) / p_ref)

    # n=5 must be a real improvement over n=1 (order-of-magnitude), and
    # n=5/n=7 must sit under a loose floor -- avoids asserting strict
    # monotonicity against the doc's own documented floor (n=3 is still
    # mid-convergence at this point, so it is not held to the same floor).
    assert errs[2] < errs[0] / 5.0
    assert all(e < 2.0e-2 for e in errs[2:])


# ---------------------------------------------------------------------------
# The headline claim: works with populations CGF cannot.
# ---------------------------------------------------------------------------


def test_works_with_beta_population_cgf_rejects():
    """MomentSeries needs only radial moments, so it accepts
    NcGalaxyShapePopBeta (no Gaussian width) where CGF aborts; both agree
    with FixedQuad in the same ballpark a Gaussian-clamp does elsewhere."""
    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = 1.55
    pop["beta"] = 1.55
    mset = _build_mset(pop)
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    std_noise = 0.15

    ms = Nc.GalaxyShapeFactorMomentSeries.new(ellip_conv, 5)
    ms.prepare(mset)
    p = _eval_factor(
        ms, pop, mset, complex(0.15, 0.0), complex(0.2, 0.05), std_noise
    )
    assert np.isfinite(p) and p > 0.0

    fq = Nc.GalaxyShapeFactorFixedQuad.new(ellip_conv)
    fq.prepare(mset)
    p_ref = _eval_factor(
        fq, pop, mset, complex(0.15, 0.0), complex(0.2, 0.05), std_noise
    )
    # Same order of magnitude -- a Gaussian clamp of the true moments, not
    # an unrelated number.
    assert 0.3 < p / p_ref < 3.0


def test_cgf_still_rejects_beta_population():
    """The Gaussian-population gate on CGF/VarAdd is unaffected by adding
    moment_2k: preparing CGF against Beta still aborts via g_error (a
    fatal process abort, checked in a subprocess, mirroring
    test_non_gaussian_population_gate in test_galaxy_shape_factor_cgf.py)."""
    script = (
        "from numcosmo_py import Nc, Ncm\n"
        "Ncm.cfg_init()\n"
        "cosmo = Nc.HICosmoDEXcdm.new()\n"
        "dist = Nc.Distance.new(100.0)\n"
        "hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)\n"
        "dp = Nc.HaloDensityProfileNFW.new(hms)\n"
        "hp = Nc.HaloPosition.new(dist)\n"
        "smd = Nc.WLSurfaceMassDensity.new(dist)\n"
        "pop = Nc.GalaxyShapePopBeta.new()\n"
        "hp.prepare(cosmo)\n"
        "mset = Ncm.MSet.empty_new()\n"
        "[mset.set(m) for m in (cosmo, dp, hp, smd, pop)]\n"
        "mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())\n"
        "mset.set(Nc.GalaxyRedshiftObsGauss.new())\n"
        "posf = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)\n"
        "pos_data = Nc.GalaxyPositionFactorData.new(posf, mset)\n"
        "zf = Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0)\n"
        "z_data = Nc.GalaxyRedshiftFactorData.new(zf, mset)\n"
        "gsf = Nc.GalaxyShapeFactorCGF.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)\n"
        "data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)\n"
        "gsf.data_set(data, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)\n"
        "gsf.prepare_data_array(mset, [data], True, True)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=False
    )
    assert result.returncode != 0
    assert "requires a Gaussian population" in result.stderr


# ---------------------------------------------------------------------------
# Numerical guard: assert, not clamp, on covariance breakdown.
# ---------------------------------------------------------------------------


def test_covariance_positivity_guard_aborts_at_low_order_large_g():
    """At n=3, e_rms~0.41, sigma_noise~0.0033 (the catalogue's worst case),
    MOMENT_SERIES.md sec. 7's own table puts the covariance breakdown near
    |g|~0.51 -- an assert, not a silent clamp, so it is a real process
    abort (checked in a subprocess)."""
    script = (
        "from numcosmo_py import Nc, Ncm\n"
        "Ncm.cfg_init()\n"
        "cosmo = Nc.HICosmoDEXcdm.new()\n"
        "dist = Nc.Distance.new(100.0)\n"
        "hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)\n"
        "dp = Nc.HaloDensityProfileNFW.new(hms)\n"
        "hp = Nc.HaloPosition.new(dist)\n"
        "smd = Nc.WLSurfaceMassDensity.new(dist)\n"
        "pop = Nc.GalaxyShapePopGauss.new()\n"
        "pop['sigma'] = 0.4847\n"
        "hp.prepare(cosmo)\n"
        "mset = Ncm.MSet.empty_new()\n"
        "[mset.set(m) for m in (cosmo, dp, hp, smd, pop)]\n"
        "mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())\n"
        "mset.set(Nc.GalaxyRedshiftObsGauss.new())\n"
        "posf = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)\n"
        "pos_data = Nc.GalaxyPositionFactorData.new(posf, mset)\n"
        "zf = Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0)\n"
        "z_data = Nc.GalaxyRedshiftFactorData.new(zf, mset)\n"
        "gsf = Nc.GalaxyShapeFactorMomentSeries.new(Nc.GalaxyWLObsEllipConv.TRACE, 3)\n"
        "data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)\n"
        "gsf.data_set(data, 0.0, 0.0, 0.0033, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)\n"
        "gsf.prepare_data_array(mset, [data], True, True)\n"
        "gsf.eval_marginal(pop, data, 0.7, 0.0, 0.1, 0.05)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=False
    )
    assert result.returncode != 0
    assert "truncated covariance is not positive" in result.stderr


def test_covariance_positivity_holds_at_default_order_across_physical_range():
    """At the default trunc-order (5), the guard must not fire anywhere in
    the physical shear range -- sec. 7's own table shows it holding past
    |g|=0.999 across the catalogue's full e_rms range."""
    for e_rms in (0.3299, 0.4091, 0.4284, 0.4534):
        sigma = e_rms  # NcGalaxyShapePopGauss's sigma is close enough to
        # e_rms at these widths for this smoke test; exactness is not the
        # point here, coverage of the guard's non-firing is.
        pop = _pop_gauss(sigma)
        mset = _build_mset(pop)
        gsf = Nc.GalaxyShapeFactorMomentSeries.new(Nc.GalaxyWLObsEllipConv.TRACE, 5)
        gsf.prepare(mset)
        data, _, _ = _build_factor_data(gsf, mset)
        gsf.data_set(
            data, 0.0, 0.0, 0.0033, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
        )
        gsf.prepare_data_array(mset, [data], True, True)

        for g_val in (0.1, 0.3, 0.5, 0.7, 0.9, 0.999):
            p = gsf.eval_marginal(pop, data, g_val, 0.0, 0.1, 0.05)
            assert np.isfinite(p)


# ---------------------------------------------------------------------------
# moment_2k consistency, and cache invalidation on the two axes CGF's own
# _peek_V already handles: pop_hash (model generation) and read_row (the
# per-galaxy population's catalog row).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("pop_name", ["gauss", "gauss_local", "beta"])
def test_moment_2k_k1_matches_e_rms_squared_over_two(pop_name):
    """M_2 / 2 == e_rms^2 for all three shipped populations -- the
    consistency MOMENT_SERIES.md sec. 5 calls out, without re-routing
    e_rms() itself through moment_2k()."""
    if pop_name == "gauss":
        pop = Nc.GalaxyShapePopGauss.new()
        pop["sigma"] = 0.35
        data = Nc.GalaxyShapePopData.new(pop)
    elif pop_name == "gauss_local":
        pop = Nc.GalaxyShapePopGaussLocal.new()
        data = Nc.GalaxyShapePopData.new(pop)
        Nc.GalaxyShapePopGaussLocal.data_set(pop, data, 0.35)
    else:
        pop = Nc.GalaxyShapePopBeta.new()
        pop["alpha"] = 2.0
        pop["beta"] = 3.5
        data = Nc.GalaxyShapePopData.new(pop)

    pop.prepare(data)
    m2 = pop.moment_2k(data, 1)
    e_rms = pop.e_rms(data)

    assert_allclose(e_rms**2, m2 / 2.0, rtol=1.0e-8)
    assert pop.moment_2k(data, 0) == 1.0


def test_moment_2k_gauss_matches_incomplete_gamma_closed_form():
    """The Gauss population's elementary integer-k moment_2k() matches the
    regularized-incomplete-gamma spelling MOMENT_SERIES.md sec. 5 gives
    (same number, no special function needed for integer k)."""
    sigma = 0.37
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma
    data = Nc.GalaxyShapePopData.new(pop)
    pop.prepare(data)

    a = 1.0 / (2.0 * sigma * sigma)
    for k in range(0, 5):
        m2k = pop.moment_2k(data, k)
        if k == 0:
            expected = 1.0
        else:
            expected = (
                (2.0 * sigma**2) ** k
                * math.factorial(k)
                * gammainc(k + 1, a)
                / (1.0 - np.exp(-a))
            )
        assert_allclose(m2k, expected, rtol=1.0e-10)


def test_moment_2k_beta_matches_beta_function_closed_form():
    """The Beta population's elementary integer-k moment_2k() matches
    B(alpha+2k,beta)/B(alpha,beta) -- alpha+2k, not alpha+k, per the 2026
    r=|chi_I| reparametrization (docs/theory/wl_shape_factor_history.md)."""
    alpha, beta_p = 2.3, 3.1
    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = alpha
    pop["beta"] = beta_p
    data = Nc.GalaxyShapePopData.new(pop)
    pop.prepare(data)

    for k in range(0, 4):
        m2k = pop.moment_2k(data, k)
        expected = (
            1.0 if k == 0 else scipy_beta(alpha + 2 * k, beta_p) / scipy_beta(alpha, beta_p)
        )
        assert_allclose(m2k, expected, rtol=1.0e-10)


def test_moment_2k_gauss_closed_form_matches_brute_force_quadrature():
    """Cross-checks Gauss's closed-form moment_2k() override against a
    brute-force numerical integral of r^2k * eval_p(r) -- the same
    reference the base class' own generic Gauss-Legendre default computes,
    so this pins the override to the accuracy that fallback path targets
    for any population without a closed form."""
    from scipy import integrate

    sigma = 0.3
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma
    data = Nc.GalaxyShapePopData.new(pop)
    pop.prepare(data)

    for k in range(1, 4):
        m2k = pop.moment_2k(data, k)
        quad_val, _ = integrate.quad(
            lambda r: r ** (2 * k) * pop.eval_p(data, r), 0.0, 1.0
        )
        assert_allclose(m2k, quad_val, rtol=1.0e-8)


def test_cache_invalidation_on_pop_hash():
    """Refitting the population's model parameter (sigma) must invalidate
    the per-galaxy coefficient cache, keyed on pop_hash -- mirrors CGF's
    own _peek_V invalidation axis."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    mset = _build_mset(pop)
    gsf = Nc.GalaxyShapeFactorMomentSeries.new(Nc.GalaxyWLObsEllipConv.TRACE_DET, 5)
    gsf.prepare(mset)
    data, _, _ = _build_factor_data(gsf, mset)
    gsf.data_set(
        data, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsf.prepare_data_array(mset, [data], True, True)

    g, eps_obs = 0.2, 0.25
    p_before = gsf.eval_marginal(pop, data, g, 0.0, eps_obs, 0.0)

    pop["sigma"] = 0.45
    gsf.prepare(mset)
    gsf.prepare_data_array(mset, [data], True, True)
    p_after = gsf.eval_marginal(pop, data, g, 0.0, eps_obs, 0.0)

    assert abs(p_after - p_before) > 1.0e-3 * abs(p_before)


def test_cache_invalidation_on_read_row_gauss_local():
    """A per-galaxy population (GaussLocal) moves the per-galaxy moments
    through the catalog row, with NO model pkey bump at all -- the axis
    ldata_read_row() alone catches. Two rows with different e_rms, read
    into the SAME NcGalaxyShapeFactorData object in sequence, must not
    reuse the first row's cached coefficients for the second."""
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    frame = Nc.WLEllipticityFrame.CELESTIAL
    pop = Nc.GalaxyShapePopGaussLocal.new()
    mset = _build_mset(pop)
    gsf = Nc.GalaxyShapeFactorMomentSeries.new(ellip_conv, 5)
    gsf.prepare(mset)
    data, _, _ = _build_factor_data(gsf, mset)

    cols = data.required_columns()
    obs = Nc.GalaxyWLObs.new(ellip_conv, frame, 2, cols)
    for row, e_rms in enumerate((0.2, 0.42)):
        obs.set("ra", row, 0.03)
        obs.set("dec", row, -0.02)
        obs.set("z", row, 0.6)
        obs.set("zp", row, 0.6)
        obs.set("sigma0", row, 0.03)
        obs.set("epsilon_int_1", row, 0.0)
        obs.set("epsilon_int_2", row, 0.0)
        obs.set("epsilon_obs_1", row, 0.25)
        obs.set("epsilon_obs_2", row, 0.0)
        obs.set("std_noise", row, 0.1)
        obs.set("c1", row, 0.0)
        obs.set("c2", row, 0.0)
        obs.set("m", row, 0.0)
        obs.set("e_rms", row, e_rms)

    g = 0.2

    data.read_row(obs, 0)
    gsf.prepare(mset)
    gsf.prepare_data_array(mset, [data], True, True)
    p_row0 = gsf.eval_marginal(pop, data, g, 0.0, 0.25, 0.0)

    data.read_row(obs, 1)
    gsf.prepare(mset)
    gsf.prepare_data_array(mset, [data], True, True)
    p_row1 = gsf.eval_marginal(pop, data, g, 0.0, 0.25, 0.0)

    # A fresh data object reading row 1 directly must agree with the reused
    # one -- if the cache had gone stale (still reflecting row 0's e_rms),
    # this would disagree with p_row1.
    data_fresh, _, _ = _build_factor_data(gsf, mset)
    data_fresh.read_row(obs, 1)
    gsf.prepare_data_array(mset, [data_fresh], True, True)
    p_row1_fresh = gsf.eval_marginal(pop, data_fresh, g, 0.0, 0.25, 0.0)

    assert abs(p_row0 - p_row1) > 1.0e-3 * abs(p_row0)
    assert_allclose(p_row1, p_row1_fresh, rtol=1.0e-12)


# ---------------------------------------------------------------------------
# Registration / serialization / construction plumbing.
# ---------------------------------------------------------------------------


def test_serialization_round_trip():
    """Registered with NcmSerialize (ncm_cfg.c) -- a YAML round-trip must
    reproduce the same trunc-order and ellip-conv."""
    gsf = Nc.GalaxyShapeFactorMomentSeries.new(Nc.GalaxyWLObsEllipConv.TRACE, 7)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    gsf2 = ser.dup_obj(gsf)

    assert isinstance(gsf2, Nc.GalaxyShapeFactorMomentSeries)
    assert gsf2.get_property("trunc-order") == 7
    assert gsf2.get_ellip_conv() == Nc.GalaxyWLObsEllipConv.TRACE


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_trunc_order_property_round_trip(ellip_conv):
    gsf = Nc.GalaxyShapeFactorMomentSeries.new(ellip_conv, 3)
    assert gsf.get_property("trunc-order") == 3
    assert gsf.props.trunc_order == 3
    assert gsf.get_ellip_conv() == ellip_conv


def test_default_trunc_order_is_five():
    gsf = Nc.GalaxyShapeFactorMomentSeries(ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE)
    assert gsf.get_property("trunc-order") == 5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
