#!/usr/bin/env python
#
# test_data_cluster_wl_factor_parity_combinations.py
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

"""Cross-check of the remaining ``NcDataClusterWLFactor`` vs legacy
``NcDataClusterWL`` combinations at the full orchestrator level.

``test_data_cluster_wl_factor_parity.py`` only ever exercises ONE
combination: shape=VarAdd+GaussGlobal x redshift=Composed+Gauss x
position=Flat. Legacy only ever provides an oracle for the *variance-add*
shape approximation (``VarAdd``), in either its Global (single shared
sigma, ``NcGalaxySDShapeHSMGaussGlobal``) or per-galaxy (``e_rms``/
``std_shape`` catalog column, ``NcGalaxySDShapeHSMGauss``) flavour, and for
the redshift scheme in either its Composed (population x observable,
``NcGalaxySDObsRedshiftGauss`` with ``use_true_z``) or Spline (per-galaxy
pre-tabulated p(z), ``NcGalaxySDObsRedshiftPz``) flavour -- position only
ever had one flavour (Flat). This file covers the three legacy-comparable
combinations left untested at the orchestrator level:

* shape=GaussGlobal x redshift=Spline   (redshift axis, new)
* shape=GaussLocal  x redshift=Spline   (both axes, new)
* shape=GaussLocal  x redshift=Composed (shape axis, new)

(GaussGlobal x Composed is the existing file's combination.) The no-legacy-
oracle schemes (Quad/SeriesLensed/FixedQuad/Laplace shape, Beta pop) are out
of scope here by design -- see the cross-check scoping discussion.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

Z_CL = 0.2
ZP_MIN, ZP_MAX = 0.0, 5.0
ELLIP_CONV = Nc.GalaxyWLObsEllipConv.TRACE_DET
FRAME = Nc.WLEllipticityFrame.CELESTIAL
SIGMA_INT_GLOBAL = 0.3
PZ_HALF_WIDTH_SIGMAS = 6.0

# (ra, dec, zp, sigma0, e_rms, e1, e2, std_noise, c1, c2, m) -- includes
# nonzero multiplicative (m) and additive (c1, c2) calibration bias on every
# row, since those terms are exactly what a bug in the shear/bias wiring
# would show up in.
_GALAXIES = [
    (0.03, 0.02, 0.60, 0.030, 0.28, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05),
    (-0.10, 0.15, 0.90, 0.040, 0.35, -0.04, 0.01, 0.05, -0.002, 0.004, -0.10),
    (0.05, -0.08, 0.15, 0.020, 0.22, 0.02, 0.03, 0.04, 0.008, -0.006, 0.08),
    (0.08, 0.05, 1.50, 0.050, 0.30, 0.03, -0.01, 0.04, 0.001, -0.001, 0.02),
    (-0.05, -0.03, 0.40, 0.030, 0.18, -0.01, 0.02, 0.035, -0.004, 0.007, -0.03),
]

_LOG10_MDELTA_GRID = [13.5, 14.0, 14.5]


def _build_lens_mset():
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)

    hp.param_set_by_name("z", Z_CL)
    hp.prepare(cosmo)

    return cosmo, dp, hp, smd, hms


def _pz_spline(zp, sigma0, n=200):
    """A per-galaxy p(z) spline for the Spline/Pz combos; domain scaled to
    sigma0 to avoid a near-delta-function integrand (see
    adaptive-spline-degenerate-integrand-oom)."""
    z_min = max(1.0e-3, zp - PZ_HALF_WIDTH_SIGMAS * sigma0)
    z_max = zp + PZ_HALF_WIDTH_SIGMAS * sigma0
    zs = np.linspace(z_min, z_max, n)
    pz_vals = np.exp(-0.5 * ((zs - zp) / sigma0) ** 2) + 1.0e-6

    xv = Ncm.Vector.new_array(zs.tolist())
    yv = Ncm.Vector.new_array(pz_vals.tolist())
    spline = Ncm.SplineCubicNotaknot.new()
    spline.set(xv, yv, True)
    return spline


def _build_mset(shape_kind, z_kind):
    """One shared mset carrying both engines' models, keyed by disjoint
    MAIN ids -- exactly the strangler-fig pattern used throughout this
    refactor (see calculators-must-not-hold-models)."""
    cosmo, dp, hp, smd, hms = _build_lens_mset()
    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd):
        mset.set(model)

    # Legacy NcDataClusterWL peeks NcGalaxySDPosition from mset by ID even
    # though Flat carries no fitted params -- required regardless of combo.
    sd_pos = Nc.GalaxySDPositionFlat.new(-0.2, 0.2, -0.2, 0.2)
    mset.set(sd_pos)

    if shape_kind == "global":
        pop_shape = Nc.GalaxyShapePopGauss.new()
        pop_shape.param_set_by_name("sigma", SIGMA_INT_GLOBAL)
        sd_shape = Nc.GalaxySDShapeHSMGaussGlobal.new(ELLIP_CONV)
        sd_shape.param_set_by_name("sigma", SIGMA_INT_GLOBAL)
    elif shape_kind == "local":
        pop_shape = Nc.GalaxyShapePopGaussLocal.new()
        sd_shape = Nc.GalaxySDShapeHSMGauss.new(ELLIP_CONV)
    else:
        raise ValueError(shape_kind)

    mset.set(pop_shape)
    mset.set(sd_shape)

    if z_kind == "composed":
        pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
        obs_z = Nc.GalaxyRedshiftObsGauss.new()
        sd_true_z = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
        sd_obs_z = Nc.GalaxySDObsRedshiftGauss.new(sd_true_z, ZP_MIN, ZP_MAX)
        mset.set(pop_z)
        mset.set(obs_z)
        mset.set(sd_obs_z)  # sd_true_z auto-registers as its submodel
    elif z_kind == "spline":
        sd_obs_z = Nc.GalaxySDObsRedshiftPz.new()
        mset.set(sd_obs_z)
    else:
        raise ValueError(z_kind)

    mset.prepare_fparam_map()
    return mset, hms


def _set_common_columns(obs, i, galaxy):
    ra, dec, _zp, _sigma0, _e_rms, e1, e2, std_noise, c1, c2, m = galaxy
    obs.set("ra", i, ra)
    obs.set("dec", i, dec)
    obs.set("z", i, 0.0)
    obs.set("epsilon_int_1", i, 0.0)
    obs.set("epsilon_int_2", i, 0.0)
    obs.set("epsilon_obs_1", i, e1)
    obs.set("epsilon_obs_2", i, e2)
    obs.set("std_noise", i, std_noise)
    obs.set("c1", i, c1)
    obs.set("c2", i, c2)
    obs.set("m", i, m)


def _build_new_obs(mset, shape_kind, z_kind):
    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)

    if z_kind == "composed":
        redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(ZP_MIN, ZP_MAX)
    elif z_kind == "spline":
        redshift_factor = Nc.GalaxyRedshiftFactorSpline.new()
    else:
        raise ValueError(z_kind)

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, len(_GALAXIES), cols)

    for i, galaxy in enumerate(_GALAXIES):
        _, _, zp, sigma0, e_rms, *_ = galaxy
        _set_common_columns(obs, i, galaxy)

        if z_kind == "composed":
            obs.set("zp", i, zp)
            obs.set("sigma0", i, sigma0)
        else:
            obs.set_pz(i, _pz_spline(zp, sigma0))

        if shape_kind == "local":
            obs.set("e_rms", i, e_rms)

    return position_factor, redshift_factor, shape_factor, obs


def _build_legacy_obs(mset, shape_kind, z_kind):
    if z_kind == "composed":
        sd_true_z = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
        sd_obs_z = Nc.GalaxySDObsRedshiftGauss.new(sd_true_z, ZP_MIN, ZP_MAX)
    elif z_kind == "spline":
        sd_obs_z = Nc.GalaxySDObsRedshiftPz.new()
    else:
        raise ValueError(z_kind)

    sd_pos = Nc.GalaxySDPositionFlat.new(-0.2, 0.2, -0.2, 0.2)

    if shape_kind == "global":
        sd_shape = Nc.GalaxySDShapeHSMGaussGlobal.new(ELLIP_CONV)
    elif shape_kind == "local":
        sd_shape = Nc.GalaxySDShapeHSMGauss.new(ELLIP_CONV)
    else:
        raise ValueError(shape_kind)

    z_data = Nc.GalaxySDObsRedshiftData.new(sd_obs_z)
    pos_data = Nc.GalaxySDPositionData.new(sd_pos, z_data)
    s_data = Nc.GalaxySDShapeData.new(sd_shape, pos_data)
    cols = Nc.GalaxySDShapeData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, len(_GALAXIES), cols)

    for i, galaxy in enumerate(_GALAXIES):
        _, _, zp, sigma0, e_rms, *_ = galaxy
        _set_common_columns(obs, i, galaxy)

        if z_kind == "composed":
            obs.set("zp", i, zp)
            obs.set("sigma_0", i, sigma0)
            obs.set("sigma_z", i, sigma0 * (1.0 + zp))
        else:
            obs.set_pz(i, _pz_spline(zp, sigma0))

        if shape_kind == "local":
            obs.set("std_shape", i, e_rms)

    return obs


_COMBOS = [
    ("global", "spline"),
    ("local", "spline"),
    ("local", "composed"),
]
_COMBO_IDS = ["shape=global,z=spline", "shape=local,z=spline", "shape=local,z=composed"]


@pytest.mark.parametrize("shape_kind,z_kind", _COMBOS, ids=_COMBO_IDS)
@pytest.mark.parametrize("log10_mdelta", _LOG10_MDELTA_GRID)
def test_m2lnL_parity(shape_kind, z_kind, log10_mdelta):
    """NcDataClusterWLFactor's total -2lnL matches legacy's (LNINT) across
    masses, for each of the three previously-uncovered combinations."""
    mset, hms = _build_mset(shape_kind, z_kind)
    hms.param_set_by_name("log10MDelta", log10_mdelta)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(
        mset, shape_kind, z_kind
    )
    legacy_obs = _build_legacy_obs(mset, shape_kind, z_kind)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    dcwl = Nc.DataClusterWL.new()
    dcwl.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    new_m2lnL = dcwlf.m2lnL_val(mset)
    old_m2lnL = dcwl.m2lnL_val(mset)

    assert_allclose(new_m2lnL, old_m2lnL, rtol=1.0e-5)


@pytest.mark.parametrize("shape_kind,z_kind", _COMBOS, ids=_COMBO_IDS)
def test_m2lnL_gal_parity(shape_kind, z_kind):
    """Per-galaxy breakdown localizes any mismatch to a single galaxy."""
    mset, hms = _build_mset(shape_kind, z_kind)
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(
        mset, shape_kind, z_kind
    )
    legacy_obs = _build_legacy_obs(mset, shape_kind, z_kind)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    dcwl = Nc.DataClusterWL.new()
    dcwl.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    n = len(_GALAXIES)
    new_gal = Ncm.Vector.new(n)
    old_gal = Ncm.Vector.new(n)

    dcwlf.eval_m2lnP_gal(mset, new_gal)
    dcwl.eval_m2lnP_gal(mset, old_gal)

    assert_allclose(new_gal.dup_array(), old_gal.dup_array(), rtol=1.0e-5)


@pytest.mark.parametrize("shape_kind,z_kind", _COMBOS, ids=_COMBO_IDS)
def test_m2lnL_parity_fixed_nodes(shape_kind, z_kind):
    """FIXED_NODES (both engines' default) agrees with legacy's own default
    exactly, for each combination -- checks the orchestrator's fixed-node
    wiring is correct for schemes other than Composed+GaussGlobal too."""
    mset, hms = _build_mset(shape_kind, z_kind)
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(
        mset, shape_kind, z_kind
    )
    legacy_obs = _build_legacy_obs(mset, shape_kind, z_kind)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    dcwl = Nc.DataClusterWL.new()  # legacy defaults to FIXED_NODES already
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    new_m2lnL = dcwlf.m2lnL_val(mset)
    old_m2lnL = dcwl.m2lnL_val(mset)

    assert_allclose(new_m2lnL, old_m2lnL, rtol=1.0e-6)


@pytest.mark.parametrize("shape_kind,z_kind", _COMBOS, ids=_COMBO_IDS)
def test_resample_matches_legacy(shape_kind, z_kind):
    """``resample()`` -- redshift draw, then position draw+radius rejection,
    then shape draw, all sharing one RNG stream -- must reproduce legacy
    seed-for-seed under the same seed for every combination, not just the
    already-covered GaussGlobal+Composed one. Compares the raw regenerated
    columns directly (ra/dec/epsilon_obs_1/2, identically named on both
    sides for every combo; zp too when z=composed) and the resulting -2lnL.

    shape=local is not bit-identical (rtol=0): both sides draw the same
    underlying uniforms, but GaussLocal resolves its per-galaxy sigma from
    e_rms via an independent bisection (self-consistent with its own
    forward map), while legacy's HSMGauss resolves sigma from std_shape via
    a closed-form inversion of an algebraically equivalent but differently
    arranged formula -- the same rtol=1e-8 caveat already established at
    the single-galaxy level (test_galaxy_shape_pop_gauss_local.py). ra/dec
    (position, shape-independent) and shape=global remain exactly bit
    identical (rtol=0).
    """
    mset, hms = _build_mset(shape_kind, z_kind)
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(
        mset, shape_kind, z_kind
    )
    legacy_obs = _build_legacy_obs(mset, shape_kind, z_kind)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    dcwl = Nc.DataClusterWL.new()
    dcwl.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    n = len(_GALAXIES)
    bit_exact_cols = ["ra", "dec"]
    shape_cols = ["epsilon_obs_1", "epsilon_obs_2"]
    if z_kind == "composed":
        bit_exact_cols.append("zp")
    shape_rtol = 0.0 if shape_kind == "global" else 1.0e-8

    for seed in (200, 201, 202):
        rng_new = Ncm.RNG.seeded_new(None, seed)
        rng_old = Ncm.RNG.seeded_new(None, seed)

        dcwlf.resample(mset, rng_new)
        dcwl.resample(mset, rng_old)

        for col in bit_exact_cols:
            new_vals = [new_obs.get(col, i) for i in range(n)]
            old_vals = [legacy_obs.get(col, i) for i in range(n)]
            assert_allclose(new_vals, old_vals, rtol=0.0, atol=0.0)

        for col in shape_cols:
            new_vals = [new_obs.get(col, i) for i in range(n)]
            old_vals = [legacy_obs.get(col, i) for i in range(n)]
            assert_allclose(new_vals, old_vals, rtol=shape_rtol, atol=1.0e-12)

        new_m2lnL = dcwlf.m2lnL_val(mset)
        old_m2lnL = dcwl.m2lnL_val(mset)
        assert_allclose(new_m2lnL, old_m2lnL, rtol=1.0e-5)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
