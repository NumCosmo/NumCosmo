#!/usr/bin/env python
#
# test_galaxy_shape_factor_cgf.py
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

"""Tests for the CGF-expansion shape factor calculator.

``NcGalaxyShapeFactorCGF`` approximates the intrinsic-ellipticity marginal by a
single Gaussian in the observed ellipticity, built from the analytic response
moments of the forward shear map plus the truncated per-component intrinsic
variance ``V = e_rms**2`` (see the class documentation). It is a second-moment
Gaussian pushforward, so it is exact in the ``g -> 0`` / ``V -> 0`` limit and
degrades as the reduced shear or the intrinsic width grows.

These tests validate it against ``NcGalaxyShapeFactorQuad`` (the exact truncated
marginal, itself validated against an independent scipy truth table), at
tolerances that reflect the measured accuracy envelope:

* noise-dominated / narrow population (``sigma_pop=0.1``, ``std_noise=0.3``):
  sub-0.1% -- essentially exact;
* production regime (``sigma_pop in {0.2, 0.3}``, ``std_noise=0.3``,
  ``|g| <= 0.35``): ~0.2%-1%;
* degradation edge (``|g|=0.35``, ``std_noise=0.03``, wide population): a few
  percent up to ~11%, driven by the map's nonlinearity over the intrinsic
  spread with too little noise to mask it.

The envelope was measured directly (both conventions) when the tolerances below
were set; they carry headroom over the measured maxima to stay non-flaky.
"""

import math
import subprocess
import sys

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


def _eval_factor(cls, pop, mset, ellip_conv, g, eps_obs, std_noise):
    gsf = cls.new(ellip_conv)
    data, _, _ = _build_factor_data(gsf, mset)
    gsf.data_set(
        data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsf.prepare_data_array(mset, [data], True, True)
    return gsf.eval_marginal(pop, data, g.real, g.imag, eps_obs.real, eps_obs.imag)


def _eval_cgf_and_quad(sigma_pop, ellip_conv, g, eps_obs, std_noise):
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma_pop
    mset = _build_mset(pop)
    cgf = _eval_factor(
        Nc.GalaxyShapeFactorCGF, pop, mset, ellip_conv, g, eps_obs, std_noise
    )
    quad = _eval_factor(
        Nc.GalaxyShapeFactorQuad, pop, mset, ellip_conv, g, eps_obs, std_noise
    )
    return cgf, quad


# Observation placed near the pushed-forward population peak so the density is
# non-tiny (mirrors the intrinsic peak at chi_I = 0.05 e^{i 0.3}).
def _eps_obs_near_peak(ellip_conv, g):
    return _shear_map(ellip_conv, g, 0.05 * np.exp(1j * 0.3))


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("sigma_pop", [0.2, 0.3])
@pytest.mark.parametrize("gmag", [0.05, 0.1, 0.2, 0.35])
def test_marginal_matches_quad_production_regime(ellip_conv, sigma_pop, gmag):
    """Production regime (std_noise=0.3): CGF tracks the exact marginal to ~1%
    across the full sigma_pop/g grid this project actually uses."""
    g = gmag * np.exp(1j * 0.6)
    eps_obs = _eps_obs_near_peak(ellip_conv, g)

    cgf, quad = _eval_cgf_and_quad(sigma_pop, ellip_conv, g, eps_obs, 0.3)

    assert_allclose(cgf, quad, rtol=1.5e-2)


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("gmag", [0.05, 0.1, 0.2, 0.35])
def test_marginal_matches_quad_noise_dominated(ellip_conv, gmag):
    """Narrow population under dominant noise (sigma_pop=0.1, std_noise=0.3):
    the pushforward's curvature terms are negligible, so CGF is near-exact."""
    g = gmag * np.exp(1j * 0.6)
    eps_obs = _eps_obs_near_peak(ellip_conv, g)

    cgf, quad = _eval_cgf_and_quad(0.1, ellip_conv, g, eps_obs, 0.3)

    assert_allclose(cgf, quad, rtol=1.0e-3)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_marginal_degradation_edge_within_documented_envelope(ellip_conv):
    """Worst documented case: large shear (|g|=0.35), tiny noise
    (std_noise=0.03), wide population (sigma_pop=0.3). A second-moment Gaussian
    pushforward genuinely degrades here (measured ~5-11%); this pins that it
    degrades gracefully within the envelope, not catastrophically."""
    g = 0.35 * np.exp(1j * 0.6)
    eps_obs = _eps_obs_near_peak(ellip_conv, g)

    cgf, quad = _eval_cgf_and_quad(0.3, ellip_conv, g, eps_obs, 0.03)

    assert_allclose(cgf, quad, rtol=0.15)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_ln_marginal_consistency(ellip_conv):
    """The ln marginal is the log of the linear marginal."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)

    gsf = Nc.GalaxyShapeFactorCGF.new(ellip_conv)
    data, _, _ = _build_factor_data(gsf, mset)
    gsf.data_set(data, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsf.prepare_data_array(mset, [data], True, True)

    p = gsf.eval_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)
    lnp = gsf.eval_ln_marginal(pop, data, 0.1, 0.05, 0.2, -0.1)

    assert p > 0.0
    assert_allclose(lnp, math.log(p), rtol=1.0e-12)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_eval_at_nodes_matches_integrand(ellip_conv):
    """The fixed-node likelihood path agrees with the integrand path (both
    routed generically by the base engine; the CGF marginal is the only
    subclass-specific piece)."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)

    gsf = Nc.GalaxyShapeFactorCGF.new(ellip_conv)
    data, pos_data, _ = _build_factor_data(gsf, mset)
    # Offset from the cluster centre: at radius 0 the reduced shear diverges.
    pos_data.ra = 0.03
    pos_data.dec = 0.02
    gsf.data_set(
        data, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05, Nc.WLEllipticityFrame.CELESTIAL
    )

    z_nodes = Ncm.Vector.new_array(np.linspace(0.05, 1.5, 25).tolist())
    out = Ncm.Vector.new(z_nodes.len())

    gsf.prepare_data_array_at_nodes(mset, [data], [z_nodes], True, True, True)
    gsf.eval_at_nodes(mset, data, z_nodes, out)

    gsf.prepare_data_array(mset, [data], True, True)
    integ = gsf.integ(mset, False)
    expected = [integ.eval(z, data) for z in z_nodes.dup_array()]

    assert_allclose(out.dup_array(), expected, rtol=1.0e-9)


def test_required_columns():
    """The factor requires its own columns plus the upstream fragments'."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)
    gsf = Nc.GalaxyShapeFactorCGF.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)
    data, _, _ = _build_factor_data(gsf, mset)

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


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_eval_marginal_dispatch(ellip_conv):
    """The marginal hooks are callable directly and mutually consistent."""
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = 0.3
    mset = _build_mset(pop)
    gsf = Nc.GalaxyShapeFactorCGF.new(ellip_conv)
    data, _, _ = _build_factor_data(gsf, mset)

    gsf.data_set(
        data, 0.04, -0.03, 0.03, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )
    gsf.prepare_data_array(mset, [data], True, True)

    p = gsf.eval_marginal(pop, data, 0.02, 0.01, 0.04, -0.03)
    lnp = gsf.eval_ln_marginal(pop, data, 0.02, 0.01, 0.04, -0.03)
    assert p > 0.0
    assert_allclose(lnp, math.log(p), rtol=1.0e-12)


def test_non_gaussian_population_gate():
    """CGF is a Gaussian-population method: preparing it against a population
    with no sigma support (Beta) aborts with a clear g_error. That is a fatal
    GLib error (process abort), so it is checked in a subprocess."""
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
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "NcGalaxyShapeFactorCGF" in result.stderr


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
