#!/usr/bin/env python
#
# test_galaxy_shape_factor_laplace.py
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

"""Tests for the Laplace-approximate shape factor calculator.

``NcGalaxyShapeFactorLaplace`` approximates the intrinsic-ellipticity
marginal with a single Gaussian expansion around the joint mode of the
integrand (see the class documentation): it is exact only in the limit of an
infinitely sharp, single-peaked population/noise combination, so these tests
compare it against ``NcGalaxyShapeFactorQuad`` (already validated against an
independent scipy truth table) rather than the scipy oracle directly, at
tolerances that reflect how peaked each case is -- tight for concentrated
populations, loose for deliberately broad ones.
"""

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

_CONVS = [Nc.GalaxyWLObsEllipConv.TRACE, Nc.GalaxyWLObsEllipConv.TRACE_DET]


def _shear_map(ellip_conv, g, chi):
    """Weak-branch forward shear map (see docs/theory/wl_ellipticity.qmd)."""
    if ellip_conv == Nc.GalaxyWLObsEllipConv.TRACE:
        denom = 1 + abs(g) ** 2 + 2 * (g * np.conj(chi)).real
        return (chi + g * (g * np.conj(chi) + 2)) / denom
    return (chi + g) / (1 + np.conj(g) * chi)


def _build_mset(pop):
    """One mset serving both engines: lens models shared by construction, the
    population model read only by the marginal (mirrors the Quad test file)."""
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


def _eval_both(pop, ellip_conv, g, eps_obs, std_noise):
    mset = _build_mset(pop)

    gsfl = Nc.GalaxyShapeFactorLaplace.new(ellip_conv)
    data_l, _, _ = _build_factor_data(gsfl, mset)
    gsfl.data_set(
        data_l, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfl.prepare_data_array(mset, [data_l], True, True)
    laplace_val = gsfl.eval_marginal(
        pop, data_l, g.real, g.imag, eps_obs.real, eps_obs.imag
    )

    gsfq = Nc.GalaxyShapeFactorQuad.new(ellip_conv)
    data_q, _, _ = _build_factor_data(gsfq, mset)
    gsfq.data_set(
        data_q, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsfq.prepare_data_array(mset, [data_q], True, True)
    quad_val = gsfq.eval_marginal(
        pop, data_q, g.real, g.imag, eps_obs.real, eps_obs.imag
    )

    return laplace_val, quad_val


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_matches_quad_narrow_gauss(ellip_conv):
    """Narrow, off-center Gauss population: the Laplace approximation should
    already be excellent here (a sharp, single Gaussian-like peak)."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.02
    g = 0.1 + 0.05j
    eps_obs = 0.2 - 0.1j

    laplace_val, quad_val = _eval_both(pop, ellip_conv, g, eps_obs, 0.25)

    assert_allclose(laplace_val, quad_val, rtol=1.0e-3)


def test_marginal_matches_quad_concentrated_beta_ring():
    """Concentrated, ring-peaked Beta population (mode away from chi_I=0):
    the joint-mode search has to resolve both radius and angle correctly."""
    mu, nu, std_noise = 0.7, 1.0e3, 0.02
    alpha, beta = mu * nu, (1.0 - mu) * nu
    mode_x = (alpha - 1.0) / (alpha + beta - 2.0)
    rho_mode = np.sqrt(mode_x)
    theta = 0.3
    g = 0.1 + 0.05j

    pop = Nc.GalaxyShapePopBeta.new()
    pop["mu"] = mu
    pop["nu"] = nu

    chi_i_peak = rho_mode * np.exp(1j * theta)
    eps_obs = _shear_map(Nc.GalaxyWLObsEllipConv.TRACE_DET, g, chi_i_peak)

    laplace_val, quad_val = _eval_both(
        pop, Nc.GalaxyWLObsEllipConv.TRACE_DET, g, eps_obs, std_noise
    )

    assert_allclose(laplace_val, quad_val, rtol=5.0e-3)


def test_marginal_broad_beta_within_documented_tolerance():
    """A broad, barely-peaked Beta population: the class doc is explicit that
    accuracy degrades here (a single Gaussian is a poor description of a
    near-flat density) -- this checks it degrades gracefully (percent level),
    not catastrophically."""
    pop = Nc.GalaxyShapePopBeta.new()
    pop["mu"] = 0.5
    pop["nu"] = 10.0
    g = 0.1 - 0.05j
    eps_obs = 0.15 + 0.1j

    laplace_val, quad_val = _eval_both(
        pop, Nc.GalaxyWLObsEllipConv.TRACE_DET, g, eps_obs, 0.2
    )

    assert_allclose(laplace_val, quad_val, rtol=0.1)


def test_ln_marginal_consistency():
    """The ln marginal is the log of the linear marginal."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)

    gsfl = Nc.GalaxyShapeFactorLaplace.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsfl, mset)
    gsfl.data_set(data, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsfl.prepare_data_array(mset, [data], True, True)

    p = gsfl.eval_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)
    lnp = gsfl.eval_ln_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)

    assert p > 0.0
    assert_allclose(lnp, math.log(p), rtol=1.0e-10)


def test_marginal_nan_for_non_peaked_population():
    """A U-shaped Beta density (alpha<1 and beta<1, diverging at BOTH disc
    boundary points x=0 and x=1) has no interior maximum along rho for a
    weak enough shear/noise combination -- the class doc says
    eval_marginal() returns NAN rather than a silently wrong number when the
    found point is not a genuine local maximum (non-positive-definite
    Hessian)."""
    pop = Nc.GalaxyShapePopBeta.new()
    pop["mu"] = 0.5
    pop["nu"] = 0.2  # alpha=beta=0.1, both < 1: diverges at x=0 and x=1.
    mset = _build_mset(pop)

    gsfl = Nc.GalaxyShapeFactorLaplace.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsfl, mset)
    gsfl.data_set(data, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsfl.prepare_data_array(mset, [data], True, True)

    val = gsfl.eval_marginal(pop, data, 0.01, 0.0, 0.02, 0.0)

    assert math.isnan(val)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
