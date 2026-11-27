#!/usr/bin/env python
#
# test_galaxy_shape_factor_wtg.py
#
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

"""Tests for the WtG shear-estimator shape factor and the shear-scatter
populations.

``NcGalaxyShapeFactorWtG`` evaluates the Applegate et al. (2014) per-galaxy
likelihood P(g_hat_t | mu_t) = K(g_hat_t - mu_t), with
g_hat_t = Re(eps e^{-2i phi}) / (s R) and mu_t the tangential component of the
calibrated shear (1 + m) g_t e^{2i phi} + c. The kernels K come from
``NcGalaxyShapePopShearScatterGauss`` and ``NcGalaxyShapePopShearScatterVoigt``.
Every expected value here is recomputed independently from those formulas.
"""

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.special import voigt_profile

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

Z_CL = 0.2

# (ra, dec, z, epsilon_obs_1, epsilon_obs_2, std_noise, c1, c2, m, responsivity, e_rms)
_GALAXIES = [
    (0.03, 0.02, 0.90, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05, 0.84, 0.36),
    (-0.10, 0.15, 0.45, -0.04, 0.01, 0.05, -0.002, 0.004, -0.10, 0.80, 0.40),
    (0.05, -0.08, 0.15, 0.02, 0.03, 0.04, 0.0, 0.0, 0.0, 0.90, 0.30),  # z below z_cl
]

_CONVS = [Nc.GalaxyWLObsEllipConv.TRACE, Nc.GalaxyWLObsEllipConv.TRACE_DET]
_S = {Nc.GalaxyWLObsEllipConv.TRACE: 2.0, Nc.GalaxyWLObsEllipConv.TRACE_DET: 1.0}


def _gauss_ln_p(delta, var):
    return -0.5 * delta**2 / var - 0.5 * math.log(2.0 * math.pi * var)


def _voigt_ln_p(delta, sigma, gamma):
    return math.log(voigt_profile(delta, sigma, gamma))


def _make_pop(kind) -> Nc.GalaxyShapePopShearScatter:
    pop: Nc.GalaxyShapePopShearScatter
    if kind == "gauss":
        pop = Nc.GalaxyShapePopShearScatterGauss.new(False)
        pop.param_set_by_name("sigma", 0.27)
    elif kind == "gauss_catalog":
        pop = Nc.GalaxyShapePopShearScatterGauss.new(True)
    else:
        pop = Nc.GalaxyShapePopShearScatterVoigt.new()
        pop.param_set_by_name("sigma", 0.22)
        pop.param_set_by_name("gamma", 0.03)
    return pop


def _expected_ln_p(pop, kind, delta, var_noise_g, eps_to_g, e_rms):
    if not pop.get_include_noise():
        var_noise_g = 0.0
    if kind == "gauss":
        return _gauss_ln_p(delta, pop.param_get_by_name("sigma") ** 2 + var_noise_g)
    if kind == "gauss_catalog":
        return _gauss_ln_p(delta, (e_rms * eps_to_g) ** 2 + var_noise_g)
    sigma_t = math.sqrt(pop.param_get_by_name("sigma") ** 2 + var_noise_g)
    return _voigt_ln_p(delta, sigma_t, pop.param_get_by_name("gamma"))


def _build_mset(pop):
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)

    hms.param_set_by_name("log10MDelta", 14.5)
    hp.param_set_by_name("z", Z_CL)
    hp.prepare(cosmo)

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop):
        mset.set(model)
    mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())
    mset.set(Nc.GalaxyRedshiftObsGauss.new())

    return mset, cosmo, dp, smd


def _build_galaxy(mset, ellip_conv, galaxy, frame=Nc.WLEllipticityFrame.CELESTIAL):
    ra, dec, z, e1, e2, std_noise, c1, c2, m, resp, e_rms = galaxy
    posf = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    pos_data = Nc.GalaxyPositionFactorData.new(posf, mset)
    zf = Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0)
    z_data = Nc.GalaxyRedshiftFactorData.new(zf, mset)
    gsf = Nc.GalaxyShapeFactorWtG.new(ellip_conv)
    data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)

    pos_data.ra = ra
    pos_data.dec = dec
    z_data.z = z
    gsf.data_set(data, e1, e2, std_noise, c1, c2, m, frame)
    gsf.data_set_responsivity(data, resp)
    data.pop_data.e_rms = e_rms

    gsf.prepare_data_array(mset, [data], True, True)
    return gsf, data


_KINDS = ["gauss", "gauss_catalog", "voigt"]


@pytest.mark.parametrize("include_noise", [True, False])
@pytest.mark.parametrize("kind", _KINDS)
def test_kernel_matches_formula(kind, include_noise):
    """The scatter kernels are the normalized Gaussian and Voigt densities,
    with the noise added to the core only when include-noise is set."""
    pop = _make_pop(kind)
    pop.set_include_noise(include_noise)
    pop_data = Nc.GalaxyShapePopData.new(pop)
    pop_data.e_rms = 0.37
    pop.prepare(pop_data)

    for delta in np.linspace(-2.5, 2.5, 41):
        for var_noise_g, eps_to_g in [(0.0, 0.6), (0.01, 0.6), (0.04, 1.0)]:
            assert_allclose(
                pop.eval_ln_p(pop_data, delta, var_noise_g, eps_to_g),
                _expected_ln_p(pop, kind, delta, var_noise_g, eps_to_g, 0.37),
                rtol=1.0e-12,
                atol=1.0e-13,
            )


@pytest.mark.parametrize("kind", _KINDS)
def test_kernel_normalized(kind):
    """The kernel integrates to one over the real line."""
    from scipy.integrate import quad

    pop = _make_pop(kind)
    pop_data = Nc.GalaxyShapePopData.new(pop)
    pop_data.e_rms = 0.37
    pop.prepare(pop_data)

    total, _ = quad(
        lambda d: pop.eval_p(pop_data, d, 0.01, 0.6), -np.inf, np.inf, limit=400
    )
    assert_allclose(total, 1.0, rtol=1.0e-7)


def test_voigt_matches_legacy():
    """Without noise the Voigt kernel equals the legacy
    NcReducedShearClusterMass P(g_obs | g_th)."""
    pop = _make_pop("voigt")
    pop.set_include_noise(False)
    pop_data = Nc.GalaxyShapePopData.new(pop)
    pop.prepare(pop_data)

    rscm = Nc.ReducedShearClusterMass.new()
    rscm.param_set_by_name("sigma", 0.22)
    rscm.param_set_by_name("Gamma", 0.03)
    cosmo = Nc.HICosmoDEXcdm.new()

    for g_obs in np.linspace(-1.0, 1.0, 21):
        g_th = 0.05
        assert_allclose(
            pop.eval_p(pop_data, g_obs - g_th, 0.123, 1.0),
            rscm.P_z_gth_gobs(cosmo, 0.5, g_th, g_obs),
            rtol=1.0e-6,
        )


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
@pytest.mark.parametrize("kind", _KINDS)
def test_integrand_matches_formula(ellip_conv, galaxy, kind):
    """The z-integrand is K(g_hat_t - mu_t(z)), with g_hat_t from the
    responsivity and the convention, and mu_t the tangential calibrated
    shear."""
    pop = _make_pop(kind)
    mset, cosmo, dp, smd = _build_mset(pop)
    gsf, data = _build_galaxy(mset, ellip_conv, galaxy)
    _, _, z, e1, e2, std_noise, c1, c2, m, resp, e_rms = galaxy

    eps_to_g = 1.0 / (_S[ellip_conv] * resp)
    rot_re, rot_im = data.get_rotation()
    radius = data.get_radius()
    hat_g_t = (e1 * rot_re + e2 * rot_im) * eps_to_g
    assert_allclose(gsf.data_get_shear_estimate(data), hat_g_t, rtol=1.0e-14)

    integ_ln = gsf.integ(mset, True)
    integ = gsf.integ(mset, False)
    for zs in [Z_CL + 0.05, 0.5, 1.0, 2.0]:
        g_t = smd.reduced_shear(dp, cosmo, radius, zs, Z_CL, Z_CL)
        mu_t = (1.0 + m) * g_t + (c1 * rot_re + c2 * rot_im)
        expected = _expected_ln_p(
            pop, kind, hat_g_t - mu_t, (std_noise * eps_to_g) ** 2, eps_to_g, e_rms
        )
        assert_allclose(integ_ln.eval(zs, data), expected, rtol=1.0e-9)
        assert_allclose(integ.eval(zs, data), math.exp(expected), rtol=1.0e-9)


def test_accepts_pop():
    """WtG takes only shear-scatter populations; every other scheme refuses
    them."""
    conv = Nc.GalaxyWLObsEllipConv.TRACE_DET
    wtg = Nc.GalaxyShapeFactorWtG.new(conv)
    var_add = Nc.GalaxyShapeFactorVarAdd.new(conv)
    disc_pop = Nc.GalaxyShapePopGauss.new()

    for kind in _KINDS:
        scatter_pop = _make_pop(kind)
        assert wtg.accepts_pop(scatter_pop)
        assert not var_add.accepts_pop(scatter_pop)

    assert not wtg.accepts_pop(disc_pop)
    assert var_add.accepts_pop(disc_pop)


def test_responsivity_column_roundtrip():
    """The responsivity is a required column, written and read back."""
    pop = _make_pop("voigt")
    mset, _, _, _ = _build_mset(pop)
    gsf, data = _build_galaxy(mset, Nc.GalaxyWLObsEllipConv.TRACE, _GALAXIES[0])
    cols = Nc.GalaxyShapeFactorData.required_columns(data)
    assert "responsivity" in cols

    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE,
        Nc.WLEllipticityFrame.CELESTIAL,
        1,
        cols,
        [Ncm.CatalogColType.DOUBLE] * len(cols),
    )
    gsf.data_set_responsivity(data, 0.77)
    data.write_row(obs, 0)
    assert obs.get("responsivity", 0) == 0.77

    gsf.data_set_responsivity(data, 0.5)
    data.read_row(obs, 0)
    assert gsf.data_get_responsivity(data) == 0.77


@pytest.mark.parametrize(
    "frame", [Nc.WLEllipticityFrame.CELESTIAL, Nc.WLEllipticityFrame.CARTESIAN]
)
@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_gen_distribution(ellip_conv, frame):
    """gen() draws g_hat_t about mu_t with the kernel's scatter, in either
    ellipticity frame: the mean and variance of the estimates match."""
    pop = _make_pop("gauss")
    mset, cosmo, dp, smd = _build_mset(pop)
    galaxy = _GALAXIES[0]
    gsf, data = _build_galaxy(mset, ellip_conv, galaxy, frame)
    _, _, z, _, _, std_noise, c1, c2, m, resp, _ = galaxy

    rng = Ncm.RNG.seeded_new(None, 1234)
    eps_to_g = 1.0 / (_S[ellip_conv] * resp)
    n = 20000
    est = np.empty(n)
    for i in range(n):
        gsf.gen(mset, data, rng)
        est[i] = gsf.data_get_shear_estimate(data)

    rot_re, rot_im = data.get_rotation()
    g_t = smd.reduced_shear(dp, cosmo, data.get_radius(), z, Z_CL, Z_CL)
    mu_t = (1.0 + m) * g_t + (c1 * rot_re + c2 * rot_im)
    sd = math.sqrt(0.27**2 + (std_noise * eps_to_g) ** 2)

    assert abs(est.mean() - mu_t) < 4.0 * sd / math.sqrt(n)
    assert_allclose(est.std(), sd, rtol=0.03)


def _build_dataset(kind, ellip_conv, n_gal, seed, log10m):
    pop = _make_pop(kind)
    mset, _, _, _ = _build_mset(pop)
    hms = mset.peek_by_name("NcHaloDensityProfile").peek_mass_summary()
    hms.param_set_by_name("log10MDelta", log10m)

    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.15, 0.15, -0.15, 0.15)
    redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(0.0, 5.0)
    shape_factor = Nc.GalaxyShapeFactorWtG.new(ellip_conv)

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)
    obs = Nc.GalaxyWLObs.new(
        ellip_conv,
        Nc.WLEllipticityFrame.CELESTIAL,
        n_gal,
        cols,
        [Ncm.CatalogColType.DOUBLE] * len(cols),
    )
    for i in range(n_gal):
        for col, val in [
            ("std_noise", 0.05),
            ("c1", 0.0),
            ("c2", 0.0),
            ("m", 0.0),
            ("sigma0", 0.03),
            ("responsivity", 0.85),
            ("e_rms", 0.36),
        ]:
            if col in cols:
                obs.set(col, i, val)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)
    dcwlf.set_cut(0.3, 3.0)
    rng = Ncm.RNG.seeded_new(None, seed)
    dcwlf.resample(mset, rng)
    return dcwlf, mset, hms


@pytest.mark.parametrize("kind", _KINDS)
def test_dataset_recovers_mass(kind):
    """On a catalog generated with the WtG model, -2lnL is finite and its
    minimum over the mass lies within 3 sigma of the input mass, sigma from
    the curvature at the input mass. Over 40 seeds the score at the input mass
    averages to zero and its variance matches the mean curvature (ratio 1.00
    Gauss, 0.90 Voigt), so the likelihood is that of the generator."""
    log10m = 14.8
    dcwlf, mset, hms = _build_dataset(
        kind, Nc.GalaxyWLObsEllipConv.TRACE, 3000, 42, log10m
    )
    grid = np.linspace(14.2, 15.4, 25)
    values = []
    for lm in grid:
        hms.param_set_by_name("log10MDelta", lm)
        values.append(dcwlf.m2lnL_val(mset))
    m2lnL = np.array(values)
    assert np.all(np.isfinite(m2lnL))

    h = 0.02
    f = []
    for d in (-h, 0.0, h):
        hms.param_set_by_name("log10MDelta", log10m + d)
        f.append(dcwlf.m2lnL_val(mset))
    curvature = 0.5 * (f[0] - 2.0 * f[1] + f[2]) / h**2
    assert curvature > 0.0

    best = grid[np.argmin(m2lnL)]
    assert abs(best - log10m) < 3.0 / math.sqrt(curvature)
