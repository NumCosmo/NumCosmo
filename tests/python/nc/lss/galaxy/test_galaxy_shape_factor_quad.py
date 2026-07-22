#!/usr/bin/env python
#
# test_galaxy_shape_factor_quad.py
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Tests for the shape factor calculator with the exact quadrature method.

``NcGalaxyShapeFactorQuad`` computes the intrinsic-ellipticity marginal
exactly (no shear-map linearization, no untruncated-Gaussian approximation),
by integrating the population density over the physical unit disc via a
smooth substitution onto the whole plane. There is no legacy oracle for this
exact method (the legacy classes only implement the variance-add
approximation), so correctness is validated against an independent
scipy-based truth table (direct polar integration over the disc) rather than
bit parity, and cross-checked for the expected convergence to
``NcGalaxyShapeFactorVarAdd`` as the reduced shear shrinks (the regime where
VarAdd's map-linearization approximation becomes exact).
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

    # Redshift models needed only to build the composed redshift fragment.
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
    physical disc, using the population's own (generically dispatched)
    eval_p and the exact shear map -- no relation to Quad's own (u,v)
    substitution or cubature machinery."""

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


_CASES = [
    # (g_1, g_2, eps_obs_1, eps_obs_2, sigma, std_noise)
    (0.15, 0.05, 0.42, -0.18, 0.30, 0.25),
    (0.02, 0.0, 0.1, -0.05, 0.30, 0.25),
    (0.10, -0.05, -0.2, 0.15, 0.15, 0.20),
    (0.05, 0.02, 0.05, 0.05, 0.05, 0.25),
]


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("case", _CASES)
def test_marginal_matches_scipy_truth_table_gauss(ellip_conv, case):
    """Quad's marginal matches an independent scipy disc integral, for the
    Gauss population."""
    g_1, g_2, eps_obs_1, eps_obs_2, sigma, std_noise = case

    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma
    mset = _build_mset(pop)

    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfq.prepare_data_array(mset, [data], True, True)

    quad_val = gsfq.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)
    exact_val = _scipy_exact_marginal(
        pop,
        data.pop_data,
        ellip_conv,
        g_1 + 1j * g_2,
        eps_obs_1 + 1j * eps_obs_2,
        std_noise,
    )

    assert_allclose(quad_val, exact_val, rtol=2.0e-4)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_scipy_truth_table_beta(ellip_conv):
    """Quad works generically for a non-Gaussian population (Beta), unlike
    VarAdd which requires nc_galaxy_shape_pop_get_sigma()."""
    pop = Nc.GalaxyShapePopBeta.new()
    mset = _build_mset(pop)

    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(data, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsfq.prepare_data_array(mset, [data], True, True)

    g = 0.1 - 0.05j
    eps_obs = 0.15 + 0.1j
    quad_val = gsfq.eval_marginal(pop, data, g.real, g.imag, eps_obs.real, eps_obs.imag)
    exact_val = _scipy_exact_marginal(pop, data.pop_data, ellip_conv, g, eps_obs, 0.2)

    assert_allclose(quad_val, exact_val, rtol=2.0e-4)


# (g_1, g_2, eps_obs_1, eps_obs_2, sigma, std_noise): populations narrower
# than std_noise, combined with an off-center g/eps_obs. An earlier version
# of this class (h-adaptive, then noise-scaled Cuhre) silently returned
# results wrong by orders of magnitude here -- its fixed-degree base rule had
# no way to notice an isolated narrow feature it never sampled near. The
# current design (Divonne seeded with explicit peak hints) resolves this
# directly rather than by tuning box size against the noise scale.
_NARROW_CASES = [
    (0.1, 0.05, 0.2, -0.1, 0.02, 0.25),
    (0.1, 0.05, 0.2, -0.1, 0.005, 0.25),
    (0.3, 0.1, -0.2, 0.15, 0.002, 0.10),
    (0.3, 0.1, -0.2, 0.15, 0.001, 0.10),
]


@pytest.mark.parametrize("case", _NARROW_CASES)
def test_marginal_matches_scipy_truth_table_narrow_gauss(case):
    """Regression test: Quad must not silently fail for a population much
    narrower than the noise, combined with an off-center g/eps_obs (see the
    class doc for the failure mode this guards against)."""
    g_1, g_2, eps_obs_1, eps_obs_2, sigma, std_noise = case

    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma
    mset = _build_mset(pop)

    gsfq = Nc.GalaxyShapeFactorQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfq.prepare_data_array(mset, [data], True, True)

    quad_val = gsfq.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)
    exact_val = _scipy_exact_marginal(
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE_DET,
        g_1 + 1j * g_2,
        eps_obs_1 + 1j * eps_obs_2,
        std_noise,
    )

    assert_allclose(quad_val, exact_val, rtol=1.0e-4)


def test_marginal_matches_scipy_truth_table_narrow_trace_convention():
    """Regression test: the population-peak hint is f_g(0), not literally g
    -- the two coincide exactly in the TRACE_DET (ellipticity) convention,
    but f_g(0) = 2g/(1+|g|^2) in the TRACE (distortion) convention. An
    earlier version of this hint hardcoded g regardless of convention, which
    silently missed the peak for a narrow, off-center TRACE population (this
    exact case returned 0 instead of the true ~1.4e-11)."""
    g_1, g_2, eps_obs_1, eps_obs_2, sigma, std_noise = 0.3, 0.1, -0.2, 0.15, 0.007, 0.10

    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma
    mset = _build_mset(pop)

    gsfq = Nc.GalaxyShapeFactorQuad.new(Nc.GalaxyWLObsEllipConv.TRACE)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfq.prepare_data_array(mset, [data], True, True)

    quad_val = gsfq.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)
    exact_val = _scipy_exact_marginal(
        pop,
        data.pop_data,
        Nc.GalaxyWLObsEllipConv.TRACE,
        g_1 + 1j * g_2,
        eps_obs_1 + 1j * eps_obs_2,
        std_noise,
    )

    assert_allclose(quad_val, exact_val, rtol=1.0e-4)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_scipy_truth_table_beta_peaked_off_center(ellip_conv):
    """Places eps_obs at the population's peak ring (mode_x != 0) to
    sanity-check the mode_x hint; not itself proof the hint is load-bearing
    (see class doc)."""
    alpha, beta, std_noise = 700.0, 300.0, 0.02
    g = 0.1 + 0.05j
    theta = 0.3

    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = alpha
    pop["beta"] = beta
    mset = _build_mset(pop)

    mode_x = (alpha - 1.0) / (alpha + beta - 2.0)
    rho_mode = np.sqrt(mode_x)
    chi_i_peak = rho_mode * np.exp(1j * theta)
    eps_obs = _shear_map(ellip_conv, g, chi_i_peak)

    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfq.prepare_data_array(mset, [data], True, True)

    quad_val = gsfq.eval_marginal(pop, data, g.real, g.imag, eps_obs.real, eps_obs.imag)
    exact_val = _scipy_exact_marginal(
        pop, data.pop_data, ellip_conv, g, eps_obs, std_noise
    )

    assert_allclose(quad_val, exact_val, rtol=1.0e-4)


def test_marginal_matches_scipy_truth_table_concentrated_beta():
    """Regression test: a concentrated (narrow-ring) Beta population, well
    beyond where the old per-factor pow()/norm evaluation would overflow to
    NaN (see nc_galaxy_shape_pop_beta.c), integrates correctly through Quad."""
    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = 300.0
    pop["beta"] = 700.0
    mset = _build_mset(pop)

    gsfq = Nc.GalaxyShapeFactorQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(data, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsfq.prepare_data_array(mset, [data], True, True)

    g = 0.1 + 0.05j
    eps_obs = 0.2 - 0.1j
    quad_val = gsfq.eval_marginal(pop, data, g.real, g.imag, eps_obs.real, eps_obs.imag)
    exact_val = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, g, eps_obs, 0.25
    )

    assert_allclose(quad_val, exact_val, rtol=1.0e-4)


def test_marginal_alpha_below_one_known_accuracy_bug():
    """Known accuracy bug: Quad loses ~11% vs scipy near g~0.18 for alpha<1
    Beta populations (FixedQuad stays ~0.2%, see
    docs/theory/wl_shape_factor_history.md); pins current behavior against
    regression."""
    pop = Nc.GalaxyShapePopBeta.new()
    pop["alpha"] = 0.6
    pop["beta"] = 4.0
    mset = _build_mset(pop)

    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    std_noise = 0.03
    eps_obs = 0.15 - 0.1j

    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfq.prepare_data_array(mset, [data], True, True)

    g = 0.18 + 0.0j
    quad_val = gsfq.eval_marginal(pop, data, g.real, g.imag, eps_obs.real, eps_obs.imag)
    exact_val = _scipy_exact_marginal(
        pop, data.pop_data, ellip_conv, g, eps_obs, std_noise
    )

    rel_err = abs(quad_val - exact_val) / exact_val
    assert rel_err > 0.05, (
        "Quad's alpha<1 accuracy bug appears fixed -- update/remove this "
        "test and re-enable Quad as a trusted cross-check for alpha<1 Beta "
        "populations."
    )
    assert rel_err < 0.20, "Quad's alpha<1 accuracy bug got worse, investigate."


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_ln_marginal_consistency(ellip_conv):
    """The ln marginal is the log of the linear marginal."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)

    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(data, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsfq.prepare_data_array(mset, [data], True, True)

    p = gsfq.eval_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)
    lnp = gsfq.eval_ln_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)

    assert p > 0.0
    assert_allclose(lnp, math.log(p), rtol=1.0e-10)


def test_required_columns():
    """Quad's ldata carries no per-row data of its own; only the upstream
    fragments' columns (plus the factor's own fixed columns) are required."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)
    gsfq = Nc.GalaxyShapeFactorQuad.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsfq, mset)

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
    for col in ("ra", "dec", "zp"):
        assert col in cols


def test_bound_reltol_properties():
    """bound/reltol getters and setters round-trip."""
    gsfq = Nc.GalaxyShapeFactorQuad.new(Nc.GalaxyWLObsEllipConv.TRACE)

    assert gsfq.get_bound() == 8.0
    assert gsfq.get_reltol() == pytest.approx(1.0e-7)

    gsfq.set_bound(12.0)
    gsfq.set_reltol(1.0e-9)
    assert gsfq.get_bound() == 12.0
    assert gsfq.get_reltol() == pytest.approx(1.0e-9)


def test_bound_reltol_gobject_property_round_trip():
    """bound/reltol are also reachable through the GObject property system
    (get_property/set_property), not just the plain getter/setter wrappers
    already checked by test_bound_reltol_properties."""
    gsfq = Nc.GalaxyShapeFactorQuad.new(Nc.GalaxyWLObsEllipConv.TRACE)

    assert gsfq.get_property("bound") == gsfq.get_bound()
    assert gsfq.get_property("reltol") == pytest.approx(gsfq.get_reltol())

    gsfq.set_property("bound", 15.0)
    gsfq.set_property("reltol", 1.0e-8)
    assert gsfq.get_bound() == 15.0
    assert gsfq.get_reltol() == pytest.approx(1.0e-8)
    assert gsfq.get_property("bound") == 15.0
    assert gsfq.get_property("reltol") == pytest.approx(1.0e-8)


def test_read_write_row_round_trip():
    """Quad's ldata carries no per-row data of its own (see
    test_required_columns), so read_row/write_row are pure pass-throughs to
    the upstream position/redshift/population fragments and the factor's own
    fixed columns -- exercised here through a real NcGalaxyWLObs catalog row,
    the path NcDataClusterWLFactor's own set_obs()/prepare() uses, which
    none of this file's other tests (all built via data_set()) touch."""
    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    frame = Nc.WLEllipticityFrame.CELESTIAL
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.25
    mset = _build_mset(pop)
    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data, _, _ = _build_factor_data(gsfq, mset)

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


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_converges_to_var_add_as_shear_shrinks(ellip_conv):
    """VarAdd's map-linearization approximation becomes more accurate as the
    reduced shear shrinks; Quad (exact) should track that trend, i.e. the
    relative disagreement should shrink monotonically with |g|."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)

    gsfva = Nc.GalaxyShapeFactorVarAdd.new(ellip_conv)
    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data_va, _, _ = _build_factor_data(gsfva, mset)
    data_q, _, _ = _build_factor_data(gsfq, mset)

    for gsf, data in ((gsfva, data_va), (gsfq, data_q)):
        gsf.data_set(
            data, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
        )
        gsf.prepare_data_array(mset, [data], True, True)

    eps_obs_1, eps_obs_2 = 0.1, -0.05
    rel_diffs = []

    for g_mag in (0.20, 0.10, 0.05, 0.01):
        va = gsfva.eval_marginal(pop, data_va, g_mag, 0.0, eps_obs_1, eps_obs_2)
        q = gsfq.eval_marginal(pop, data_q, g_mag, 0.0, eps_obs_1, eps_obs_2)
        rel_diffs.append(abs(va - q) / q)

    # Not required to be strictly monotonic step-by-step (the linearization
    # error can have a sign crossing at a particular g for a particular
    # eps_obs), but the largest shear must disagree markedly more than the
    # smallest, and the weak-shear limit must be tight.
    assert rel_diffs[0] > 5.0 * rel_diffs[-1]
    assert rel_diffs[-1] < 1.0e-2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
