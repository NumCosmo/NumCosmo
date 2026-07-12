#!/usr/bin/env python
#
# test_galaxy_shape_factor_knots.py
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

"""Tests for the fixed-knot shape factor calculator playground.

``NcGalaxyShapeFactorKnots`` is not a validated likelihood backend (see its
class documentation): it is a scaffold for trying different fixed-node
schemes against ``NcGalaxyShapeFactorQuad``'s Cuhre reference. These tests
check the scaffold's plumbing (property round-trips, node generation sanity)
and that the plain fixed scheme (CARTESIAN) is in the right ballpark against
the same scipy truth table used for Quad, at a generous tolerance -- not
meant to certify accuracy there. The adaptive scheme (CARTESIAN_ADAPTIVE),
which re-centers and rescales the same node count onto the joint mode of the
integrand, is checked more tightly: it is specifically meant to resolve
narrow/off-center features the fixed scheme cannot at any practical node
count (see the class documentation).
"""

import pytest
import numpy as np
from scipy import integrate
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def _build_mset(pop):
    """One mset serving the calculator; lens models are fixed, only the
    population model is read by the marginal (see NcGalaxyShapeFactorQuad's
    test file, which this mirrors)."""
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
    physical disc (same oracle used by NcGalaxyShapeFactorQuad's tests)."""

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


def test_bound_n_method_properties():
    """bound/n/method getters and setters round-trip."""
    gsfk = Nc.GalaxyShapeFactorKnots.new(Nc.GalaxyWLObsEllipConv.TRACE)

    assert gsfk.get_bound() == 8.0
    assert gsfk.get_n() == 32
    assert gsfk.get_method() == Nc.GalaxyShapeFactorKnotsMethod.CARTESIAN

    gsfk.set_bound(6.0)
    gsfk.set_n(16)
    assert gsfk.get_bound() == 6.0
    assert gsfk.get_n() == 16


def test_nodes_regenerate_on_property_change():
    """peek_nodes/peek_weights reflect bound/n after they change, and the 1D
    rule integrates a constant exactly to 2*bound (Gauss-Legendre exactness)."""
    gsfk = Nc.GalaxyShapeFactorKnots.new(Nc.GalaxyWLObsEllipConv.TRACE)

    for n, bound in ((32, 8.0), (10, 3.0), (64, 5.0)):
        gsfk.set_n(n)
        gsfk.set_bound(bound)
        nodes = gsfk.peek_nodes()
        weights = gsfk.peek_weights()

        assert nodes.len() == n
        assert weights.len() == n
        assert_allclose(sum(weights.get(i) for i in range(n)), 2.0 * bound, rtol=1.0e-12)
        assert nodes.get(0) == pytest.approx(-nodes.get(n - 1), abs=1.0e-10)


def test_adaptive_reference_nodes_are_canonical():
    """The CARTESIAN_ADAPTIVE scheme's peek_nodes/peek_weights are the
    REFERENCE rule on [-1,1] (rescaled per evaluation, see the class doc),
    unlike CARTESIAN's fixed [-bound,bound] layout."""
    gsfk = Nc.GalaxyShapeFactorKnots.new(Nc.GalaxyWLObsEllipConv.TRACE)
    gsfk.set_method(Nc.GalaxyShapeFactorKnotsMethod.CARTESIAN_ADAPTIVE)
    gsfk.set_n(16)
    gsfk.set_bound(6.0)

    nodes = gsfk.peek_nodes()
    weights = gsfk.peek_weights()

    assert_allclose(sum(weights.get(i) for i in range(16)), 2.0, rtol=1.0e-12)
    assert nodes.get(0) == pytest.approx(-nodes.get(15), abs=1.0e-10)


@pytest.mark.parametrize("ellip_conv", [Nc.GalaxyWLObsEllipConv.TRACE, Nc.GalaxyWLObsEllipConv.TRACE_DET])
def test_marginal_matches_scipy_truth_table_narrow_adaptive(ellip_conv):
    """Regression test for the CARTESIAN_ADAPTIVE scheme: a population much
    narrower than CARTESIAN's fixed [-bound,bound] box, combined with an
    off-center g/eps_obs. CARTESIAN itself returns ~0 here even at n=64 (its
    uniform grid never samples near the isolated peak, the same failure mode
    documented for the pre-hint version of NcGalaxyShapeFactorQuad);
    CARTESIAN_ADAPTIVE, re-centered and rescaled on the joint mode, matches
    the scipy truth table tightly already at a modest node count."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.02
    mset = _build_mset(pop)

    g_1, g_2, eps_obs_1, eps_obs_2, std_noise = 0.1, 0.05, 0.2, -0.1, 0.25

    gsfk = Nc.GalaxyShapeFactorKnots.new(ellip_conv)
    gsfk.set_method(Nc.GalaxyShapeFactorKnotsMethod.CARTESIAN_ADAPTIVE)
    gsfk.set_n(32)
    data, _, _ = _build_factor_data(gsfk, mset)
    gsfk.data_set(data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsfk.prepare_data_array(mset, [data], True, True)

    knots_val = gsfk.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)
    exact_val = _scipy_exact_marginal(
        pop, data.pop_data, ellip_conv, g_1 + 1j * g_2, eps_obs_1 + 1j * eps_obs_2, std_noise,
    )

    assert_allclose(knots_val, exact_val, rtol=1.0e-6)


def test_marginal_matches_scipy_truth_table_concentrated_beta_adaptive():
    """Same regression as above, for a concentrated ring-peaked Beta
    population (mode away from chi_I=0): CARTESIAN_ADAPTIVE's joint-mode
    re-centering (nc_galaxy_shape_intrinsic_mode_find()) has to find the
    population's peak ring, not just the naive noiseless inverse map."""
    mu, nu, std_noise = 0.7, 1.0e3, 0.02
    alpha, beta = mu * nu, (1.0 - mu) * nu
    mode_x = (alpha - 1.0) / (alpha + beta - 2.0)
    rho_mode = np.sqrt(mode_x)
    theta = 0.3
    g = 0.1 + 0.05j

    pop = Nc.GalaxyShapePopBeta.new()
    pop["mu"] = mu
    pop["nu"] = nu
    mset = _build_mset(pop)

    chi_i_peak = rho_mode * np.exp(1j * theta)
    eps_obs = _shear_map(Nc.GalaxyWLObsEllipConv.TRACE_DET, g, chi_i_peak)

    gsfk = Nc.GalaxyShapeFactorKnots.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    gsfk.set_method(Nc.GalaxyShapeFactorKnotsMethod.CARTESIAN_ADAPTIVE)
    gsfk.set_n(32)
    data, _, _ = _build_factor_data(gsfk, mset)
    gsfk.data_set(data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsfk.prepare_data_array(mset, [data], True, True)

    knots_val = gsfk.eval_marginal(pop, data, g.real, g.imag, eps_obs.real, eps_obs.imag)
    exact_val = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET, g, eps_obs, std_noise
    )

    assert_allclose(knots_val, exact_val, rtol=1.0e-6)


def test_marginal_matches_scipy_truth_table_gauss_loose():
    """The Cartesian GL scheme lands in the right ballpark of the exact
    disc integral for a moderate-width Gaussian population. This is a loose
    sanity check on the integrand/re-centering plumbing, not an accuracy
    claim (see NcGalaxyShapeFactorQuad's tight version of this same test)."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)

    gsfk = Nc.GalaxyShapeFactorKnots.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    gsfk.set_n(64)
    data, _, _ = _build_factor_data(gsfk, mset)
    gsfk.data_set(data, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsfk.prepare_data_array(mset, [data], True, True)

    g_1, g_2, eps_obs_1, eps_obs_2, std_noise = 0.1, 0.05, 0.2, -0.1, 0.25
    knots_val = gsfk.eval_marginal(pop, data, g_1, g_2, eps_obs_1, eps_obs_2)
    exact_val = _scipy_exact_marginal(
        pop, data.pop_data, Nc.GalaxyWLObsEllipConv.TRACE_DET,
        g_1 + 1j * g_2, eps_obs_1 + 1j * eps_obs_2, std_noise,
    )

    assert_allclose(knots_val, exact_val, rtol=5.0e-2)


def test_ln_marginal_consistency():
    """The ln marginal is the log of the linear marginal."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)

    gsfk = Nc.GalaxyShapeFactorKnots.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsfk, mset)
    gsfk.data_set(data, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsfk.prepare_data_array(mset, [data], True, True)

    p = gsfk.eval_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)
    lnp = gsfk.eval_ln_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)

    assert p > 0.0
    assert_allclose(lnp, np.log(p), rtol=1.0e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
