#!/usr/bin/env python
#
# test_data_cluster_wl_binned.py
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

"""Tests of the binned cluster weak-lensing likelihood ``NcDataClusterWLBinned``.

The estimators are checked against direct numpy implementations of their equations
(Umetsu et al. 2020 for Delta Sigma, the weighted radial profile and the Schrabback et
al. 2018 model for g_t, also checked against CLMM), and the exact
model against noiseless catalogs, for which the measured profile equals the expected
value of the estimator.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.cluster_wl import set_pz_from_redshift_factor
from numcosmo_py.helper import duplicate_via_serialization

Ncm.cfg_init()

Z_CL = 0.3
LOG10M = 14.5
EDGES = np.geomspace(0.3, 3.0, 6)
SIGMA_POP = 0.26


def _build_mset():
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hms.param_set_by_name("log10MDelta", LOG10M)
    hms.param_set_by_name("cDelta", 4.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    hp.param_set_by_name("z", Z_CL)
    hp.prepare(cosmo)
    smd.prepare(cosmo)
    scatter = Nc.GalaxyWLBinnedScatterGauss.new()

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, scatter):
        mset.set(model)
    mset.prepare_fparam_map()

    return mset


def _models(mset):
    return (
        mset.peek(Nc.HICosmo.id()),
        mset.peek(Nc.HaloPosition.id()),
        mset.peek(Nc.WLSurfaceMassDensity.id()),
        mset.peek(Nc.HaloDensityProfile.id()),
    )


def _galaxy_geometry(mset, ra, dec):
    cosmo, hp, _, _ = _models(mset)
    radius = np.empty(len(ra))
    phi = np.empty(len(ra))
    for i, (r, d) in enumerate(zip(ra, dec)):
        theta, phi[i] = hp.polar_angles(r, d)
        radius[i] = hp.projected_radius(cosmo, theta)
    return radius, phi


def _reduced_shear(mset, radius, z):
    """Reduced tangential shear of each galaxy (zero in front of the cluster)."""
    cosmo, _, smd, dp = _models(mset)
    g = np.zeros(len(radius))
    for i, (r, zs) in enumerate(zip(radius, z)):
        if zs > Z_CL:
            s = 1.0 / smd.sigma_critical(cosmo, zs, Z_CL, Z_CL)
            sig = smd.sigma(dp, cosmo, r, Z_CL)
            dsig = smd.sigma_excess(dp, cosmo, r, Z_CL)
            g[i] = s * dsig / (1.0 - s * sig)
    return g


def _build_catalog(
    mset,
    n=600,
    seed=1,
    ellip_conv=Nc.GalaxyWLObsEllipConv.TRACE_DET,
    noise=0.0,
    with_bias=True,
):
    """A catalog with point redshifts and epsilon_obs = (1 + m) g + c (+ noise)."""
    rng = np.random.default_rng(seed)
    ra = rng.uniform(-0.25, 0.25, n)
    dec = rng.uniform(-0.25, 0.25, n)
    z = rng.uniform(0.1, 1.5, n)
    std_noise = rng.uniform(0.05, 0.15, n)
    e_rms = rng.uniform(0.3, 0.4, n)
    m = rng.uniform(-0.05, 0.05, n) if with_bias else np.zeros(n)
    c1 = rng.uniform(-1e-3, 1e-3, n) if with_bias else np.zeros(n)
    c2 = rng.uniform(-1e-3, 1e-3, n) if with_bias else np.zeros(n)

    radius, phi = _galaxy_geometry(mset, ra, dec)
    g_t = _reduced_shear(mset, radius, z)
    g = g_t * np.exp(2.0j * phi)
    e = (1.0 + m) * g + (c1 + 1.0j * c2)
    e += noise * (rng.normal(size=n) + 1.0j * rng.normal(size=n))

    cols = [
        "ra",
        "dec",
        "z",
        "epsilon_obs_1",
        "epsilon_obs_2",
        "std_noise",
        "e_rms",
        "m",
        "c1",
        "c2",
    ]
    obs = Nc.GalaxyWLObs.new(ellip_conv, Nc.WLEllipticityFrame.CELESTIAL, n, cols)
    for i in range(n):
        obs.set("ra", i, ra[i])
        obs.set("dec", i, dec[i])
        obs.set("z", i, z[i])
        obs.set("epsilon_obs_1", i, e.real[i])
        obs.set("epsilon_obs_2", i, e.imag[i])
        obs.set("std_noise", i, std_noise[i])
        obs.set("e_rms", i, e_rms[i])
        obs.set("m", i, m[i])
        obs.set("c1", i, c1[i])
        obs.set("c2", i, c2[i])

    gal = {
        "ra": ra,
        "dec": dec,
        "z": z,
        "e1": e.real,
        "e2": e.imag,
        "std_noise": std_noise,
        "e_rms": e_rms,
        "m": m,
        "c1": c1,
        "c2": c2,
        "radius": radius,
        "phi": phi,
    }
    return obs, gal


def _tangential(gal):
    rot = np.exp(-2.0j * gal["phi"])
    e = (gal["e1"] + 1.0j * gal["e2"]) * rot
    c = (gal["c1"] + 1.0j * gal["c2"]) * rot
    return e.real, c.real


def _inv_sigma_crit(mset, z):
    cosmo, _, smd, _ = _models(mset)
    return np.array(
        [
            1.0 / smd.sigma_critical(cosmo, zs, Z_CL, Z_CL) if zs > Z_CL else 0.0
            for zs in z
        ]
    )


def _bin_index(radius):
    sel = (radius >= EDGES[0]) & (radius < EDGES[-1])
    b = np.clip(np.searchsorted(EDGES, radius, side="right") - 1, 0, len(EDGES) - 2)
    return sel, b


def _col(binned, name):
    return np.asarray(binned.get_column(name).dup_array())


def _cov(binned):
    cov = binned.peek_cov()
    n = binned.len()
    return np.array([[cov.get(i, j) for j in range(n)] for i in range(n)])


def _ref_delta_sigma(mset, gal, shear_scale, use_resp, subtract=True):
    """Umetsu et al. 2020, Eqs. 4, 6-9, 11, 14, 23."""
    e_t, c_t = _tangential(gal)
    s = _inv_sigma_crit(mset, gal["z"])
    w = np.where(s > 0, s**2 / (gal["std_noise"] ** 2 + gal["e_rms"] ** 2), 0.0)
    sel, b = _bin_index(gal["radius"])
    out = {k: [] for k in ("value", "radius", "var", "inv_sigma_crit", "offset")}
    for k in range(len(EDGES) - 1):
        idx = sel & (b == k) & (w > 0)
        if not idx.any():
            continue
        wk = w[idx]
        W = wk.sum()
        mult = (wk * (1.0 + gal["m"][idx])).sum() / W
        resp = 1.0 - (wk * gal["e_rms"][idx] ** 2).sum() / W if use_resp else 1.0
        scale = shear_scale * resp
        ds = (wk * e_t[idx] / s[idx]).sum() / (scale * W)
        ct = (wk * c_t[idx] / s[idx]).sum() / W
        out["offset"].append(0.0 if subtract else ct / mult)
        if subtract:
            ds -= ct
        out["value"].append(ds / mult)
        out["radius"].append(W / (wk / gal["radius"][idx]).sum())
        out["var"].append(1.0 / (scale**2 * mult**2 * W))
        out["inv_sigma_crit"].append((wk * s[idx]).sum() / W)
    return {k: np.array(v) for k, v in out.items()}


def _ref_gt(mset, gal, shear_scale, use_resp, error_type, beta_weighted=False):
    """Weighted radial profile with the Mandelbaum et al. 2018a calibration."""
    cosmo, _, smd, _ = _models(mset)
    e_t, _ = _tangential(gal)
    sel, b = _bin_index(gal["radius"])
    w = 1.0 / (SIGMA_POP**2 + 2.0 * gal["std_noise"] ** 2)
    ws = w[sel]
    m_hat = (ws * gal["m"][sel]).sum() / ws.sum()
    c1_hat = (ws * gal["c1"][sel]).sum() / ws.sum()
    c2_hat = (ws * gal["c2"][sel]).sum() / ws.sum()
    resp = 1.0 - (ws * gal["e_rms"][sel] ** 2).sum() / ws.sum() if use_resp else 1.0
    scale = shear_scale * resp
    c_hat_t = ((c1_hat + 1.0j * c2_hat) * np.exp(-2.0j * gal["phi"])).real
    g_t = (e_t / scale - c_hat_t) / (1.0 + m_hat)
    sci = smd.sigma_critical_infinity(cosmo, Z_CL, Z_CL)
    beta = sci * _inv_sigma_crit(mset, gal["z"][sel])
    out = {k: [] for k in ("value", "radius", "var")}
    for k in range(len(EDGES) - 1):
        idx = sel & (b == k)
        if not idx.any():
            continue
        wk = w[idx]
        W = wk.sum()
        gm = (wk * g_t[idx]).sum() / W
        if error_type == Nc.GalaxyWLBinnedEstimatorGtError.STD_ERR:
            # Standard error of the weighted mean (CLMM error_model='ste')
            var = ((wk * g_t[idx] ** 2).sum() / W - gm**2) * (wk**2).sum() / W**2
        else:
            v = (gal["e_rms"][idx] ** 2 + gal["std_noise"][idx] ** 2) / (
                scale * (1.0 + m_hat)
            ) ** 2
            var = (wk**2 * v).sum() / W**2
        out["value"].append(gm)
        out["radius"].append((wk * gal["radius"][idx]).sum() / W)
        out["var"].append(var)
    res = {k: np.array(v) for k, v in out.items()}
    wb = ws if beta_weighted else np.ones_like(ws)
    res["beta_s_mean"] = (wb * beta).sum() / wb.sum()
    res["beta_s2_mean"] = (wb * beta**2).sum() / wb.sum()
    return res


@pytest.mark.parametrize(
    "ellip_conv, shear_scale, use_resp",
    [
        (Nc.GalaxyWLObsEllipConv.TRACE_DET, 1.0, False),
        (Nc.GalaxyWLObsEllipConv.TRACE, 2.0, True),
    ],
)
@pytest.mark.parametrize("subtract", [True, False])
def test_delta_sigma_build_matches_reference(
    ellip_conv, shear_scale, use_resp, subtract
):
    mset = _build_mset()
    obs, gal = _build_catalog(mset, ellip_conv=ellip_conv, noise=0.2)
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()
    est.set_subtract_additive_bias(subtract)
    binned = est.build(mset, obs, EDGES.tolist())
    ref = _ref_delta_sigma(mset, gal, shear_scale, use_resp, subtract)

    assert binned.peek_estimator() == "delta_sigma"
    assert_allclose(_col(binned, "value"), ref["value"], rtol=1e-11)
    assert_allclose(_col(binned, "radius"), ref["radius"], rtol=1e-12)
    assert_allclose(_col(binned, "inv_sigma_crit"), ref["inv_sigma_crit"], rtol=1e-12)
    assert_allclose(_col(binned, "offset"), ref["offset"], rtol=1e-11, atol=1e-30)
    assert_allclose(np.diag(_cov(binned)), ref["var"], rtol=1e-12)
    assert_allclose(_cov(binned) - np.diag(np.diag(_cov(binned))), 0.0)


@pytest.mark.parametrize(
    "ellip_conv, shear_scale, use_resp",
    [
        (Nc.GalaxyWLObsEllipConv.TRACE_DET, 1.0, False),
        (Nc.GalaxyWLObsEllipConv.TRACE, 2.0, True),
    ],
)
@pytest.mark.parametrize(
    "error_type",
    [
        Nc.GalaxyWLBinnedEstimatorGtError.STD_ERR,
        Nc.GalaxyWLBinnedEstimatorGtError.SHAPE_NOISE,
    ],
)
@pytest.mark.parametrize("beta_weighted", [False, True])
def test_gt_build_matches_reference(
    ellip_conv, shear_scale, use_resp, error_type, beta_weighted
):
    mset = _build_mset()
    obs, gal = _build_catalog(mset, ellip_conv=ellip_conv, noise=0.2)
    est = Nc.GalaxyWLBinnedEstimatorGt.new()
    est.set_error_type(error_type)
    est.set_beta_weighted(beta_weighted)
    assert est.get_sigma_pop() == SIGMA_POP
    binned = est.build(mset, obs, EDGES.tolist())
    ref = _ref_gt(mset, gal, shear_scale, use_resp, error_type, beta_weighted)

    assert binned.peek_estimator() == "gt"
    assert_allclose(_col(binned, "value"), ref["value"], rtol=1e-11)
    assert_allclose(_col(binned, "radius"), ref["radius"], rtol=1e-12)
    assert_allclose(np.diag(_cov(binned)), ref["var"], rtol=1e-10)
    assert_allclose(_col(binned, "beta_s_mean"), ref["beta_s_mean"], rtol=1e-12)
    assert_allclose(_col(binned, "beta_s2_mean"), ref["beta_s2_mean"], rtol=1e-12)


def test_delta_sigma_model_and_m2lnL():
    """Umetsu et al. 2020, Eqs. 21-22, and the Gaussian -2 ln L."""
    mset = _build_mset()
    cosmo, _, smd, dp = _models(mset)
    obs, _ = _build_catalog(mset, noise=0.2)
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    binned = data.peek_binned()

    radius = _col(binned, "radius")
    s = _col(binned, "inv_sigma_crit")
    sig = np.array([smd.sigma(dp, cosmo, r, Z_CL) for r in radius])
    dsig = np.array([smd.sigma_excess(dp, cosmo, r, Z_CL) for r in radius])
    ref_model = dsig / (1.0 - s * sig)

    model = Ncm.Vector.new(binned.len())
    data.eval_model(mset, model)
    assert_allclose(model.dup_array(), ref_model, rtol=1e-12)

    diff = _col(binned, "value") - ref_model
    chi2 = diff @ np.linalg.solve(_cov(binned), diff)
    assert_allclose(data.m2lnL_val(mset), chi2, rtol=1e-10)


def test_gt_model():
    """Schrabback et al. 2018, Eq. 12 (CLMM order2)."""
    mset = _build_mset()
    cosmo, _, smd, dp = _models(mset)
    obs, _ = _build_catalog(mset, noise=0.2)
    est = Nc.GalaxyWLBinnedEstimatorGt.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    binned = data.peek_binned()

    radius = _col(binned, "radius")
    b1 = _col(binned, "beta_s_mean")
    b2 = _col(binned, "beta_s2_mean")
    sci = smd.sigma_critical_infinity(cosmo, Z_CL, Z_CL)
    kappa = np.array([smd.sigma(dp, cosmo, r, Z_CL) for r in radius]) / sci
    gamma = np.array([smd.sigma_excess(dp, cosmo, r, Z_CL) for r in radius]) / sci
    ref_model = (
        b1 * gamma / (1.0 - b1 * kappa) * (1.0 + (b2 / b1**2 - 1.0) * b1 * kappa)
    )

    model = Ncm.Vector.new(binned.len())
    data.eval_model(mset, model)
    assert_allclose(model.dup_array(), ref_model, rtol=1e-12)


@pytest.mark.parametrize(
    "est_type", [Nc.GalaxyWLBinnedEstimatorDeltaSigma, Nc.GalaxyWLBinnedEstimatorGt]
)
@pytest.mark.parametrize("with_bias", [False, True])
def test_exact_model_is_expected_value_of_noiseless_data(est_type, with_bias):
    """For a noiseless catalog the profile is the expected value of the estimator."""
    mset = _build_mset()
    obs, _ = _build_catalog(mset, with_bias=with_bias)
    est = est_type.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    binned = data.peek_binned()
    value = _col(binned, "value")

    data.set_exact(True)
    model = Ncm.Vector.new(binned.len())
    data.eval_model(mset, model)
    assert_allclose(model.dup_array(), value, rtol=1e-10)

    # The classic approximations are close but not exact; the g_t model uses the
    # lensing efficiency of the whole sample in every bin, which differs from the one
    # of the few galaxies of a bin
    data.set_exact(False)
    data.eval_model(mset, model)
    assert_allclose(
        model.dup_array(),
        value,
        rtol=0.05 if est_type == Nc.GalaxyWLBinnedEstimatorDeltaSigma else 0.25,
    )


def _pz_table(zp, sigma):
    z = np.linspace(0.0, 3.0, 120)
    p = np.exp(-0.5 * ((z - zp) / sigma) ** 2)
    return Nc.GalaxyPzBSplineSqrt.new_from_table(
        Ncm.Vector.new_array(z), Ncm.Vector.new_array(p)
    )


def test_redshift_density_average():
    """<Sigma_cr^-1> over a p(z) against a fine trapezoid rule (Eq. 6)."""
    mset = _build_mset()
    cosmo, _, smd, _ = _models(mset)
    obs, _ = _build_catalog(mset, n=4)
    pz = _pz_table(0.7, 0.2)
    obs.set_pz(1, pz)
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()

    zmin, zmax = pz.get_range()
    zz = np.linspace(Z_CL, zmax, 20001)[1:]
    pp = np.array([pz.eval_p(z) for z in zz])
    ss = np.array([1.0 / smd.sigma_critical(cosmo, z, Z_CL, Z_CL) for z in zz])
    ref1 = np.trapezoid(pp * ss, zz) / pz.get_norm()
    ref2 = np.trapezoid(pp * ss**2, zz) / pz.get_norm()

    mean1, mean2 = est.mean_inv_sigma_crit(mset, obs, 1)
    assert_allclose(mean1, ref1, rtol=1e-6)
    assert_allclose(mean2, ref2, rtol=1e-6)


def test_redshift_density_average_two_narrow_peaks():
    """The adaptive integrals find narrow, separated peaks of a p(z)."""
    mset = _build_mset()
    cosmo, _, smd, _ = _models(mset)
    obs, _ = _build_catalog(mset, n=4)
    z = np.linspace(0.0, 6.0, 1201)
    p = np.exp(-0.5 * ((z - 0.6) / 0.02) ** 2)
    p += 0.2 * np.exp(-0.5 * ((z - 4.3) / 0.02) ** 2)
    pz = Nc.GalaxyPzBSplineSqrt.new_from_table(
        Ncm.Vector.new_array(z), Ncm.Vector.new_array(p)
    )
    obs.set_pz(1, pz)
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()

    zz = np.linspace(Z_CL, 6.0, 200001)[1:]
    pp = np.array([pz.eval_p(z) for z in zz])
    ss = np.array([1.0 / smd.sigma_critical(cosmo, z, Z_CL, Z_CL) for z in zz])
    mean1, mean2 = est.mean_inv_sigma_crit(mset, obs, 1)
    assert_allclose(mean1, np.trapezoid(pp * ss, zz) / pz.get_norm(), rtol=1e-5)
    assert_allclose(mean2, np.trapezoid(pp * ss**2, zz) / pz.get_norm(), rtol=1e-5)


def test_exact_model_with_redshift_densities():
    """With p(z) for every galaxy, the exact model converges in exact-intervals and the
    classic model stays close to it."""
    mset = _build_mset()
    obs, _ = _build_catalog(mset, n=200, noise=0.2)
    rng = np.random.default_rng(7)
    for i in range(obs.len()):
        obs.set_pz(i, _pz_table(rng.uniform(0.4, 1.2), 0.15))
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    n = data.peek_binned().len()

    data.set_exact(True)
    m_20 = Ncm.Vector.new(n)
    m_80 = Ncm.Vector.new(n)
    data.eval_model(mset, m_20)
    est.set_exact_intervals(80)
    data.eval_model(mset, m_80)
    assert_allclose(m_20.dup_array(), m_80.dup_array(), rtol=1e-6)

    data.set_exact(False)
    m_classic = Ncm.Vector.new(n)
    data.eval_model(mset, m_classic)
    assert_allclose(m_classic.dup_array(), m_80.dup_array(), rtol=0.02)


def test_empty_and_foreground_bins_are_dropped():
    mset = _build_mset()
    obs, gal = _build_catalog(mset, n=300)
    edges = np.array([0.3, 0.5, 0.50001, 1.0, 2.0, 50.0, 60.0])
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()
    binned = est.build(mset, obs, edges.tolist())
    # [0.5, 0.50001) and [50, 60) hold no galaxy
    assert_allclose(_col(binned, "r_lo"), [0.3, 0.50001, 1.0, 2.0])
    n_gal = _col(binned, "n_gal")
    for k, (lo, hi) in enumerate([(0.3, 0.5), (0.50001, 1.0), (1.0, 2.0), (2.0, 50.0)]):
        assert n_gal[k] == ((gal["radius"] >= lo) & (gal["radius"] < hi)).sum()

    # A bin of galaxies in front of the cluster has zero weight for Delta Sigma
    z = np.array([obs.get("z", i) for i in range(obs.len())])
    for i in np.flatnonzero((gal["radius"] >= 1.0) & (gal["radius"] < 2.0)):
        obs.set("z", i, 0.1)
    binned = est.build(mset, obs, edges.tolist())
    assert_allclose(_col(binned, "r_lo"), [0.3, 0.50001, 2.0])
    members = binned.peek_members()
    rows = [members.get_int("row", j) for j in range(members.len())]
    assert all(z[r] > Z_CL and not (1.0 <= gal["radius"][r] < 2.0) for r in rows)


def test_external_profile():
    mset = _build_mset()
    radius = [0.5, 1.0, 2.0]
    value = [0.03, 0.02, 0.01]
    cov = Ncm.Matrix.new(3, 3)
    cov.set_identity()
    cov.scale(1.0e-4)
    binned = Nc.GalaxyWLObsBinned.new_from_arrays(
        radius, value, cov, ["beta_s_mean", "beta_s2_mean"]
    )
    binned.set_column("beta_s_mean", [0.5, 0.5, 0.5])
    binned.set_column("beta_s2_mean", (0.3, 0.3, 0.3))
    assert not binned.has_members()

    est = Nc.GalaxyWLBinnedEstimatorGt.new()
    data = Nc.DataClusterWLBinned.new(est, binned)
    m2lnL = data.m2lnL_val(mset)
    assert np.isfinite(m2lnL)
    assert data.get_length() == 3


def test_serialization_round_trip():
    mset = _build_mset()
    obs, _ = _build_catalog(mset, n=200, noise=0.2)
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    data.set_exact(True)
    data_dup = duplicate_via_serialization(data)

    assert isinstance(data_dup, Nc.DataClusterWLBinned)
    assert data_dup.get_exact()
    assert data_dup.peek_binned().has_members()
    assert_allclose(data_dup.m2lnL_val(mset), data.m2lnL_val(mset), rtol=1e-14)


def test_scatter_gauss_lnnorm():
    mset = _build_mset()
    obs, _ = _build_catalog(mset, noise=0.2)
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    cov = _cov(data.peek_binned())
    m2lnL = data.m2lnL_val(mset)

    scatter = mset.peek(Nc.GalaxyWLBinnedScatter.id())
    scatter.set_include_lnnorm(True)
    lnnorm = len(cov) * np.log(2.0 * np.pi) + np.linalg.slogdet(cov)[1]
    assert_allclose(data.m2lnL_val(mset), m2lnL + lnnorm, rtol=1e-12)


def test_resample_is_unbiased():
    mset = _build_mset()
    obs, _ = _build_catalog(mset, noise=0.2)
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    binned = data.peek_binned()
    n = binned.len()
    model = Ncm.Vector.new(n)
    data.eval_model(mset, model)
    sigma = np.sqrt(np.diag(_cov(binned)))

    rng = Ncm.RNG.seeded_new(None, 1234)
    n_real = 4000
    draws = np.empty((n_real, n))
    for k in range(n_real):
        data.resample(mset, rng)
        draws[k] = _col(binned, "value")

    pulls = (draws - model.dup_array()) / sigma
    assert_allclose(pulls.mean(axis=0), 0.0, atol=5.0 / np.sqrt(n_real))
    assert_allclose(pulls.std(axis=0), 1.0, atol=0.05)


@pytest.mark.parametrize(
    "est_type", [Nc.GalaxyWLBinnedEstimatorDeltaSigma, Nc.GalaxyWLBinnedEstimatorGt]
)
def test_fit_recovers_mass_from_noiseless_exact(est_type):
    mset = _build_mset()
    obs, _ = _build_catalog(mset)
    est = est_type.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    data.set_exact(True)
    if est_type == Nc.GalaxyWLBinnedEstimatorGt:
        # A noiseless catalog has no scatter for the standard error
        binned = est.build(mset, obs, EDGES.tolist())
        est.set_error_type(Nc.GalaxyWLBinnedEstimatorGtError.SHAPE_NOISE)
        binned = est.build(mset, obs, EDGES.tolist())
        data.set_binned(binned)

    hms = mset.peek(Nc.HaloMassSummary.id())
    hms.param_set_desc(
        "log10MDelta", {"fit": True, "lower-bound": 13.0, "upper-bound": 16.0}
    )
    hms.param_set_by_name("log10MDelta", 14.0)
    mset.prepare_fparam_map()

    dset = Ncm.Dataset.new_array([data])
    lh = Ncm.Likelihood.new(dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )
    fit.run(Ncm.FitRunMsgs.NONE)
    assert_allclose(hms.param_get_by_name("log10MDelta"), LOG10M, atol=1e-4)


@pytest.mark.parametrize(
    "est_type", [Nc.GalaxyWLBinnedEstimatorDeltaSigma, Nc.GalaxyWLBinnedEstimatorGt]
)
def test_cartesian_frame_matches_celestial(est_type):
    """The same galaxies in the Cartesian frame (conjugated ellipticities and additive
    biases) give the same profile and the same exact model."""
    mset = _build_mset()
    obs_cel, gal = _build_catalog(mset, noise=0.2)
    n = obs_cel.len()
    cols = list(obs_cel.peek_columns())
    obs_car = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, Nc.WLEllipticityFrame.CARTESIAN, n, cols
    )
    sign = {"epsilon_obs_2": -1.0, "c2": -1.0}
    for col in cols:
        for i in range(n):
            obs_car.set(col, i, sign.get(col, 1.0) * obs_cel.get(col, i))

    edges = EDGES.tolist()
    est = est_type.new()
    d_cel = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs_cel, edges)
    d_car = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs_car, edges)
    for col in ("value", "radius", "offset"):
        assert_allclose(
            _col(d_car.peek_binned(), col),
            _col(d_cel.peek_binned(), col),
            rtol=1e-12,
            atol=1e-300,
        )

    m_cel = Ncm.Vector.new(d_cel.get_length())
    m_car = Ncm.Vector.new(d_car.get_length())
    d_cel.set_exact(True)
    d_car.set_exact(True)
    d_cel.eval_model(mset, m_cel)
    d_car.eval_model(mset, m_car)
    assert_allclose(m_car.dup_array(), m_cel.dup_array(), rtol=1e-12)


def test_register_shared_shares_source():
    """Copies made through a serializer after register_shared() reference the source
    catalog of the profile instead of copying it, and keep their own profile."""
    mset = _build_mset()
    obs, _ = _build_catalog(mset, n=200, noise=0.2)
    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    data.register_shared(ser)
    dup = ser.dup_obj(data)

    assert dup.peek_binned() is not data.peek_binned()
    assert dup.peek_binned().peek_source() is obs

    rng = Ncm.RNG.seeded_new(None, 5)
    before = _col(data.peek_binned(), "value")
    dup.resample(mset, rng)
    assert_allclose(_col(data.peek_binned(), "value"), before)


def test_gt_matches_clmm():
    """The g_t profile and model against CLMM's make_radial_profile and its order2
    reduced tangential shear, given the same per-galaxy shears, weights and lensing
    efficiencies."""
    clmm_stat = pytest.importorskip("clmm.utils.statistic")
    mset = _build_mset()
    cosmo, _, smd, dp = _models(mset)
    obs, gal = _build_catalog(mset, noise=0.2)
    est = Nc.GalaxyWLBinnedEstimatorGt.new()
    data = Nc.DataClusterWLBinned.new_from_obs(est, mset, obs, EDGES.tolist())
    binned = data.peek_binned()

    # Per-galaxy calibrated shears, as fed to CLMM
    e_t, _ = _tangential(gal)
    sel, _ = _bin_index(gal["radius"])
    w = 1.0 / (SIGMA_POP**2 + 2.0 * gal["std_noise"] ** 2)
    ws = w[sel]
    m_hat = (ws * gal["m"][sel]).sum() / ws.sum()
    c_hat = (
        (ws * gal["c1"][sel]).sum() + 1.0j * (ws * gal["c2"][sel]).sum()
    ) / ws.sum()
    c_hat_t = (c_hat * np.exp(-2.0j * gal["phi"])).real
    g_t = (e_t - c_hat_t) / (1.0 + m_hat)

    mean_r, mean_g, err_g, _, _, _ = clmm_stat.compute_radial_averages(
        gal["radius"][sel], g_t[sel], EDGES, error_model="ste", weights=ws
    )
    assert_allclose(_col(binned, "radius"), mean_r, rtol=1e-12)
    assert_allclose(_col(binned, "value"), mean_g, rtol=1e-11)
    assert_allclose(np.sqrt(np.diag(_cov(binned))), err_g, rtol=1e-10)

    # CLMM order2 (parent_class._eval_reduced_tangential_shear), with the same
    # gamma_inf, kappa_inf and <beta_s>, <beta_s^2>
    sci = smd.sigma_critical_infinity(cosmo, Z_CL, Z_CL)
    b1 = _col(binned, "beta_s_mean")
    b2 = _col(binned, "beta_s2_mean")
    radius = _col(binned, "radius")
    gamma_inf = np.array([smd.sigma_excess(dp, cosmo, r, Z_CL) for r in radius]) / sci
    kappa_inf = np.array([smd.sigma(dp, cosmo, r, Z_CL) for r in radius]) / sci
    gt_clmm = b1 * gamma_inf / (1.0 - b1 * kappa_inf)
    gt_clmm *= 1.0 + (b2 / (b1 * b1) - 1.0) * b1 * kappa_inf
    model = Ncm.Vector.new(binned.len())
    data.eval_model(mset, model)
    assert_allclose(model.dup_array(), gt_clmm, rtol=1e-12)


def test_pz_from_redshift_factor_matches_direct_integral():
    """The mocks of the CLI get a NcGalaxyPz tabulated from their redshift factor: its
    <Sigma_cr^-1> is the average over p(z) = P_pop(z) P_obs(zp | z) / W(z)."""
    zp_min, zp_max = 0.0, 5.0
    mset = _build_mset()
    cosmo, _, smd, _ = _models(mset)
    pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    obs_z = Nc.GalaxyRedshiftObsGauss.new()
    mset.set(pop_z)
    mset.set(obs_z)
    mset.prepare_fparam_map()

    obs, _ = _build_catalog(mset, n=4)
    cols = list(obs.peek_columns()) + ["zp", "sigma0"]
    obs_c = Nc.GalaxyWLObs.new(obs.get_ellip_conv(), obs.get_coord(), obs.len(), cols)
    for i in range(obs.len()):
        for col in obs.peek_columns():
            obs_c.set(col, i, obs.get(col, i))
        obs_c.set("zp", i, [0.5, 0.8, 1.1, 0.35][i])
        obs_c.set("sigma0", i, 0.05)

    est = Nc.GalaxyWLBinnedEstimatorDeltaSigma.new()
    set_pz_from_redshift_factor(
        mset, obs_c, Nc.GalaxyRedshiftFactorComposed.new(zp_min, zp_max), 1.0e-6
    )

    zz = np.linspace(1.0e-4, 6.0, 60001)
    ss = np.array(
        [
            1.0 / smd.sigma_critical(cosmo, z, Z_CL, Z_CL) if z > Z_CL else 0.0
            for z in zz
        ]
    )
    data_z = Nc.GalaxyRedshiftObsData.new(obs_z)
    for i in range(obs_c.len()):
        data_z.read_row(obs_c, i)
        pp = np.array(
            [
                pop_z.eval(z)
                * obs_z.eval(data_z, z)
                / obs_z.window_mass(data_z, z, zp_min, zp_max)
                for z in zz
            ]
        )
        ref = np.trapezoid(pp * ss, zz) / np.trapezoid(pp, zz)
        mean1, _ = est.mean_inv_sigma_crit(mset, obs_c, i)
        assert_allclose(mean1, ref, rtol=1e-5)
