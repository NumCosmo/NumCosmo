#!/usr/bin/env python
#
# test_data_cluster_wl_factor_parity.py
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

"""Tier 1 parity test: ``NcDataClusterWLFactor`` vs the legacy ``NcDataClusterWL``.

``NcDataClusterWLFactor`` is the new per-galaxy-Factor orchestrator
(``NcGalaxyPositionFactorFlat`` / ``NcGalaxyRedshiftFactorComposed`` /
``NcGalaxyShapeFactorVarAdd``); it must reproduce legacy's ``-2lnL`` exactly
(within adaptive-quadrature tolerance) when configured with ``VarAdd``, which
is itself already proven bit-identical to the legacy
``NcGalaxySDShapeHSMGaussGlobal`` shape model at the single-galaxy level (see
``test_galaxy_shape_factor_var_add.py``), and ``Composed`` +
``GalaxyRedshiftPopLSSTSRD`` + ``GalaxyRedshiftObsGauss`` already proven
bit-identical to ``NcGalaxySDObsRedshiftGauss`` wrapping
``NcGalaxySDTrueRedshiftLSSTSRD`` (see
``test_galaxy_redshift_factor_composed.py``). This test is what actually
exercises the *new* code -- ``NcDataClusterWLFactor``'s own ``prepare()``
rebuild/hash-refresh cascade and its ``m2lnL_val`` (the z-integral, the z_cl
split, the log-sum-exp panel recombination) -- against legacy configured with
``LNINT`` (not its default ``FIXED_NODES``), so both sides are independently
coded adaptive log-domain 1D integrators, not merely re-running the same
code.

Two independently-built ``NcGalaxyWLObs`` catalogs feed the same physical
galaxy values into each engine, one per engine's own column-name schema
(legacy stores ``sigma_0``/``sigma_z`` explicitly, spelled with an
underscore; the new ``Composed`` scheme stores only ``sigma0`` and derives
``sigma_z = sigma0 * (1 + zp)`` internally) -- deliberately not one shared
catalog, to sidestep that naming mismatch entirely.
"""

import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

SIGMA_INT = 0.3
Z_CL = 0.2
ZP_MIN, ZP_MAX = 0.0, 5.0
ELLIP_CONV = Nc.GalaxyWLObsEllipConv.TRACE_DET
FRAME = Nc.WLEllipticityFrame.CELESTIAL

# (ra, dec, zp, sigma0, e1, e2, std_noise, c1, c2, m)
_GALAXIES = [
    (0.03, 0.02, 0.60, 0.030, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05),
    (-0.10, 0.15, 0.90, 0.040, -0.04, 0.01, 0.05, -0.002, 0.004, -0.10),
    (0.05, -0.08, 0.15, 0.020, 0.02, 0.03, 0.04, 0.0, 0.0, 0.0),  # zp below z_cl
    (0.08, 0.05, 1.50, 0.050, 0.03, -0.01, 0.04, 0.001, -0.001, 0.02),
    (-0.05, -0.03, 0.35, 0.025, -0.01, 0.02, 0.035, 0.0, 0.0, 0.0),  # straddles z_cl
]

_LOG10_MDELTA_GRID = [13.5, 14.0, 14.5, 15.0]


def _build_mset():
    """One shared mset: lens-side models common to both engines, plus each
    engine's own galaxy-population-side models (peeked by model ID, exactly
    like NcDataClusterWL peeks NcGalaxySDShape/NcGalaxySDObsRedshift/
    NcGalaxySDPosition; see nc_data_cluster_wl.c:1017-1019)."""
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)

    hp.param_set_by_name("z", Z_CL)
    hp.prepare(cosmo)

    # New-side population/observable models.
    pop_shape = Nc.GalaxyShapePopGauss.new()
    pop_shape.param_set_by_name("sigma", SIGMA_INT)
    pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    obs_z = Nc.GalaxyRedshiftObsGauss.new()

    # Legacy-side models, peeked from mset by NcDataClusterWL itself.
    sd_true_z = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
    sd_obs_z = Nc.GalaxySDObsRedshiftGauss.new(sd_true_z, ZP_MIN, ZP_MAX)
    sd_pos = Nc.GalaxySDPositionFlat.new(-0.2, 0.2, -0.2, 0.2)
    sd_shape = Nc.GalaxySDShapeHSMGaussGlobal.new(ELLIP_CONV)
    sd_shape.param_set_by_name("sigma", SIGMA_INT)

    mset = Ncm.MSet.empty_new()
    for model in (
        cosmo,
        dp,
        hp,
        smd,
        pop_shape,
        pop_z,
        obs_z,
        sd_obs_z,
        sd_pos,
        sd_shape,  # sd_true_z auto-registers as sd_obs_z's submodel
    ):
        mset.set(model)
    mset.prepare_fparam_map()

    return mset, hms


def _build_new_obs(mset):
    """New-engine Factor objects + one NcGalaxyWLObs built via required_columns."""
    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(ZP_MIN, ZP_MAX)
    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, len(_GALAXIES), cols)

    for i, (ra, dec, zp, sigma0, e1, e2, std_noise, c1, c2, m) in enumerate(_GALAXIES):
        obs.set("ra", i, ra)
        obs.set("dec", i, dec)
        obs.set("z", i, 0.0)  # true z: latent, unread by the likelihood
        obs.set("zp", i, zp)
        obs.set("sigma0", i, sigma0)
        obs.set("epsilon_int_1", i, 0.0)  # intrinsic draw: unread by the likelihood
        obs.set("epsilon_int_2", i, 0.0)
        obs.set("epsilon_obs_1", i, e1)
        obs.set("epsilon_obs_2", i, e2)
        obs.set("std_noise", i, std_noise)
        obs.set("c1", i, c1)
        obs.set("c2", i, c2)
        obs.set("m", i, m)

    return position_factor, redshift_factor, shape_factor, obs


def _build_legacy_obs():
    """Legacy SD objects + one NcGalaxyWLObs built via its own required_columns."""
    sd_true_z = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
    sd_obs_z = Nc.GalaxySDObsRedshiftGauss.new(sd_true_z, ZP_MIN, ZP_MAX)
    sd_pos = Nc.GalaxySDPositionFlat.new(-0.2, 0.2, -0.2, 0.2)
    sd_shape = Nc.GalaxySDShapeHSMGaussGlobal.new(ELLIP_CONV)

    z_data = Nc.GalaxySDObsRedshiftData.new(sd_obs_z)
    pos_data = Nc.GalaxySDPositionData.new(sd_pos, z_data)
    s_data = Nc.GalaxySDShapeData.new(sd_shape, pos_data)
    cols = Nc.GalaxySDShapeData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, len(_GALAXIES), cols)

    for i, (ra, dec, zp, sigma0, e1, e2, std_noise, c1, c2, m) in enumerate(_GALAXIES):
        obs.set("ra", i, ra)
        obs.set("dec", i, dec)
        obs.set("z", i, 0.0)
        obs.set("zp", i, zp)
        obs.set("sigma_0", i, sigma0)
        obs.set("sigma_z", i, sigma0 * (1.0 + zp))
        obs.set("epsilon_int_1", i, 0.0)
        obs.set("epsilon_int_2", i, 0.0)
        obs.set("epsilon_obs_1", i, e1)
        obs.set("epsilon_obs_2", i, e2)
        obs.set("std_noise", i, std_noise)
        obs.set("c1", i, c1)
        obs.set("c2", i, c2)
        obs.set("m", i, m)

    return obs


@pytest.mark.parametrize("log10_mdelta", _LOG10_MDELTA_GRID)
def test_m2lnL_parity(log10_mdelta):
    """NcDataClusterWLFactor's total -2lnL matches legacy's (LNINT) across masses."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", log10_mdelta)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)
    legacy_obs = _build_legacy_obs()

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)

    dcwl = Nc.DataClusterWL.new()
    dcwl.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    new_m2lnL = dcwlf.m2lnL_val(mset)
    old_m2lnL = dcwl.m2lnL_val(mset)

    assert_allclose(new_m2lnL, old_m2lnL, rtol=1.0e-5)


@pytest.mark.parametrize("log10_mdelta", _LOG10_MDELTA_GRID)
def test_m2lnL_parity_fixed_nodes(log10_mdelta):
    """NcDataClusterWLFactor's FIXED_NODES matches legacy's own (default)
    FIXED_NODES exactly -- both sides call the same underlying VarAdd-style
    closed-form eval_at_nodes/make_fixed_nodes machinery, so this checks the
    orchestrator-level wiring (grid construction, control-variate combine),
    not a numerically-independent scheme like the LNINT comparison above."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", log10_mdelta)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)
    legacy_obs = _build_legacy_obs()

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    dcwl = Nc.DataClusterWL.new()  # legacy defaults to FIXED_NODES already
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    new_m2lnL = dcwlf.m2lnL_val(mset)
    old_m2lnL = dcwl.m2lnL_val(mset)

    assert_allclose(new_m2lnL, old_m2lnL, rtol=1.0e-8)


def test_m2lnL_fixed_nodes_matches_lnint():
    """Cross-check within the new class alone: FIXED_NODES and LNINT are two
    independently-coded, mathematically exact schemes for the same integral
    -- they should agree with each other within adaptive-quadrature
    tolerance, independent of legacy."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf_fixed = Nc.DataClusterWLFactor.new(
        position_factor, redshift_factor, shape_factor
    )
    dcwlf_fixed.set_obs(new_obs)
    dcwlf_fixed.set_prec(1.0e-8)
    dcwlf_fixed.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    dcwlf_lnint = Nc.DataClusterWLFactor.new(
        position_factor, redshift_factor, shape_factor
    )
    dcwlf_lnint.set_obs(new_obs)
    dcwlf_lnint.set_prec(1.0e-8)
    dcwlf_lnint.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    fixed_m2lnL = dcwlf_fixed.m2lnL_val(mset)
    lnint_m2lnL = dcwlf_lnint.m2lnL_val(mset)

    assert_allclose(fixed_m2lnL, lnint_m2lnL, rtol=1.0e-6)


def test_fixed_nodes_cache_consistency_across_revisits():
    """FIXED_NODES matches legacy across a sequence that revisits mass and
    lens-redshift values -- the scenario that catches stale per-galaxy
    fixed-node-grid caches (see submodel-pkey-does-not-propagate lessons)."""
    mset, hms = _build_mset()
    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)
    legacy_obs = _build_legacy_obs()

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    dcwl = Nc.DataClusterWL.new()
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    for log10_mdelta in (13.5, 14.5, 13.5, 15.0, 14.5):
        hms.param_set_by_name("log10MDelta", log10_mdelta)
        assert_allclose(dcwlf.m2lnL_val(mset), dcwl.m2lnL_val(mset), rtol=1.0e-8)

    hms.param_set_by_name("log10MDelta", 14.0)
    hp = mset.peek(Nc.HaloPosition.id())
    cosmo = mset.peek(Nc.HICosmo.id())

    for z_cl in (0.2, 0.35, 0.2, 0.5):
        hp.param_set_by_name("z", z_cl)
        hp.prepare(cosmo)
        assert_allclose(dcwlf.m2lnL_val(mset), dcwl.m2lnL_val(mset), rtol=1.0e-8)


def test_fixed_nodes_correct_under_angular_only_changes():
    """Changing halo_position's RA/Dec (not z) must not corrupt the
    fixed-node grid: the grid depends only on z_cl, but halo_position's
    whole-model hash (used elsewhere for radius_changed) bumps on *any* of
    its parameters, including RA/Dec -- e.g. in a miscentering model. This
    guards against conflating "halo_position changed" with "z_cl changed"."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)
    legacy_obs = _build_legacy_obs()

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    dcwl = Nc.DataClusterWL.new()
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    hp = mset.peek(Nc.HaloPosition.id())
    cosmo = mset.peek(Nc.HICosmo.id())

    for ra in (0.0, 0.05, -0.03, 0.05):
        hp.param_set_by_name("ra", ra)
        hp.prepare(cosmo)
        assert_allclose(dcwlf.m2lnL_val(mset), dcwl.m2lnL_val(mset), rtol=1.0e-8)


def test_fixed_nodes_correct_after_switching_integ_method_mid_run():
    """Switching integ-method from LNINT to FIXED_NODES on the *same*
    instance, with no other mset change in between, must still populate
    every FIXED_NODES-only per-galaxy cache (grid, crit, sigma) from
    scratch. Found as a real bug during review: `sigma`'s gate
    (`radius_changed || optzs_changed`) had no way to detect "FIXED_NODES
    has never run on this object before" the way `fixed_grid_changed`'s
    GSL_NAN-seeded sentinel does -- radius_hash/optzs_hash were already kept
    current by the preceding LNINT cycles' own (integ-method-agnostic)
    steps, so neither flag fired on the switch, leaving `sigma_cache`
    uninitialized. Fixed by also gating sigma on `fixed_grid_changed`."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)
    legacy_obs = _build_legacy_obs()

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    dcwl = Nc.DataClusterWL.new()
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    # Prime LNINT state -- radius_hash/optzs_hash become current here, with
    # FIXED_NODES-only caches (grid/crit/sigma) never touched.
    dcwlf.m2lnL_val(mset)
    dcwlf.m2lnL_val(mset)  # idempotent second call, nothing changed

    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)
    new_m2lnL = dcwlf.m2lnL_val(mset)
    old_m2lnL = dcwl.m2lnL_val(mset)

    assert_allclose(new_m2lnL, old_m2lnL, rtol=1.0e-8)

    # And switching back to LNINT must still match the original LNINT value.
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    assert_allclose(dcwlf.m2lnL_val(mset), new_m2lnL, rtol=1.0e-5)


def test_fixed_nodes_resample_reuse_matches_lnint():
    """Regression test for a real, severe bug found during a mass-recovery
    timing sweep: `ncm_data_resample()` calls `prepare()` *before* invoking
    the `resample` virtual (see `ncm_data_resample()` in `ncm_data.c`), so
    FIXED_NODES's per-galaxy z-grid could get built from whatever data was
    current *before* resample() overwrote it (e.g. set_obs()'s placeholder
    zeros on a fresh instance). Because the grid is deliberately gated only
    by z_cl/n-nodes/rule-n (not by any per-galaxy data hash, so that a
    mass-only fit doesn't rebuild it every iteration), that staleness was
    never naturally corrected afterwards -- reusing one instance across
    repeated resample()+fit cycles (exactly what every MC/timing harness in
    this project does) silently corrupted every likelihood evaluation,
    driving fits to the mass parameter's bound. This used to be documented
    here as a "known shared limitation" inherited from legacy
    `NcDataClusterWL` (which has the identical bug); it is now fixed by
    having `resample()` mark `obs_changed` so the very next prepare() cycle
    redoes every per-galaxy cache from the fresh data, regardless of which
    integ-method is active or whether any model parameter has moved."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf_reused = Nc.DataClusterWLFactor.new(
        position_factor, redshift_factor, shape_factor
    )
    dcwlf_reused.set_obs(new_obs)
    dcwlf_reused.set_prec(1.0e-8)
    # integ-method set to FIXED_NODES *before* the first resample() -- the
    # exact ordering that used to poison the grid with pre-resample data.
    dcwlf_reused.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    dcwlf_lnint = Nc.DataClusterWLFactor.new(
        position_factor, redshift_factor, shape_factor
    )
    dcwlf_lnint.set_obs(new_obs)
    dcwlf_lnint.set_prec(1.0e-8)
    dcwlf_lnint.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    # Several resample() cycles on the *same* reused instances, mirroring a
    # real MC harness -- must agree with LNINT on every single one, not just
    # the first.
    for seed in (42, 43, 44):
        rng_reused = Ncm.RNG.seeded_new(None, seed)
        dcwlf_reused.resample(mset, rng_reused)
        reused_val = dcwlf_reused.m2lnL_val(mset)

        rng_lnint = Ncm.RNG.seeded_new(None, seed)
        dcwlf_lnint.resample(mset, rng_lnint)
        lnint_val = dcwlf_lnint.m2lnL_val(mset)

        assert_allclose(reused_val, lnint_val, rtol=1.0e-5)


def test_resample_matches_legacy():
    """``resample()`` itself -- not just static ``m2lnL_val`` on a fixed
    catalog -- must generate the *same* realizations as legacy, seed for
    seed. Both engines call gen() in the identical order (redshift,
    position, shape -- see ``_nc_data_cluster_wl_factor_resample`` and
    ``_nc_data_cluster_wl_resample``), and each individual gen() step was
    already proven bit-parity in isolation (shape:
    ``test_galaxy_shape_factor_var_add.py::test_gen_parity_legacy``); this
    is the first check that the *whole* per-galaxy resample pipeline
    (redshift z draw -> position ra/dec draw+radius rejection -> shape
    draw, all sharing one RNG stream) reproduces legacy bit-for-bit once
    composed together in the orchestrator. Compares the raw regenerated
    columns directly (ra/dec/zp/epsilon_obs_1/2, all identically named on
    both sides) as well as the resulting -2lnL."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)
    legacy_obs = _build_legacy_obs()

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)

    dcwl = Nc.DataClusterWL.new()
    dcwl.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    dcwl.set_obs(legacy_obs)
    dcwl.set_prec(1.0e-8)

    n = len(_GALAXIES)

    for seed in (100, 101, 102):
        rng_new = Ncm.RNG.seeded_new(None, seed)
        rng_old = Ncm.RNG.seeded_new(None, seed)

        dcwlf.resample(mset, rng_new)
        dcwl.resample(mset, rng_old)

        for col in ("ra", "dec", "zp", "epsilon_obs_1", "epsilon_obs_2"):
            new_vals = [new_obs.get(col, i) for i in range(n)]
            old_vals = [legacy_obs.get(col, i) for i in range(n)]
            assert_allclose(new_vals, old_vals, rtol=0.0, atol=0.0)

        new_m2lnL = dcwlf.m2lnL_val(mset)
        old_m2lnL = dcwl.m2lnL_val(mset)
        assert_allclose(new_m2lnL, old_m2lnL, rtol=1.0e-5)

        new_gal = Ncm.Vector.new(n)
        old_gal = Ncm.Vector.new(n)
        dcwlf.eval_m2lnP_gal(mset, new_gal)
        dcwl.eval_m2lnP_gal(mset, old_gal)
        assert_allclose(new_gal.dup_array(), old_gal.dup_array(), rtol=1.0e-5)


def test_fixed_nodes_correct_after_changing_n_nodes_rule_n_mid_run():
    """`n-nodes`/`rule-n` are orchestrator properties, not NcmModel
    parameters -- changing them bumps no pkey at all, so
    `fixed_grid_changed` cannot rely on any hash to detect a resize. Found
    as a real bug during review: without dedicated tracking, changing
    n-nodes/rule-n mid-run left z_nodes_per_galaxy/fixed_bg_nodes sized to
    the *old* panel count while `_eval_m2lnP_fixed` computed n_total_nodes
    fresh from the *new* n-nodes/rule-n every call -- a size mismatch that
    corrupted the heap via ncm_integral_fixed_integ_vec_mult() on a
    stale-sized subvector. Fixed by tracking last-seen n-nodes/rule-n
    (0 forces the first rebuild, mirroring fixed_nodes_zcl's GSL_NAN) and
    folding them into fixed_grid_changed."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)
    legacy_obs = _build_legacy_obs()

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)
    dcwlf.set_n_nodes(10)
    dcwlf.set_rule_n(5)
    dcwlf.m2lnL_val(mset)

    for n_nodes, rule_n in ((6, 3), (12, 7), (2, 1), (10, 5)):
        dcwlf.set_n_nodes(n_nodes)
        dcwlf.set_rule_n(rule_n)
        new_m2lnL = dcwlf.m2lnL_val(mset)

        dcwl = Nc.DataClusterWL.new()
        dcwl.set_obs(legacy_obs)
        dcwl.set_prec(1.0e-8)
        dcwl.set_property("n-nodes", n_nodes)
        dcwl.set_property("rule-n", rule_n)
        old_m2lnL = dcwl.m2lnL_val(mset)

        assert_allclose(new_m2lnL, old_m2lnL, rtol=1.0e-8)


def test_fixed_nodes_correct_after_swapping_obs_to_different_sized_catalog():
    """set_obs() to a catalog with a *different galaxy count*, with no other
    mset change, must fully rebuild every FIXED_NODES-only cache. Found as a
    real, pre-existing Phase-1 bug exposed by this review, not introduced by
    FIXED_NODES itself: `self->obs_changed` was read (in the `pos_changed`/
    `z_changed`/`radius_changed`/`optzs_changed`/`pop_changed` computations)
    *after* prepare() had already reset it to FALSE a few lines above --
    always evaluating the `self->obs_changed ||` term as a no-op. This was
    masked for those five flags because their own hash fields start at 0,
    which any real hash is virtually certain to differ from, so "first ever"
    was still detected by the hash comparison alone. It was not masked for
    `fixed_grid_changed`: a pure obs-swap changes no NcmModel at all, so
    every one of its own comparisons (z_cl, n_nodes, rule_n) still matched
    their last-seen values, and the (silently dead) `self->obs_changed` term
    was the only thing that could have signaled "rebuild anyway" -- leaving
    the grid/crit/sigma caches sized and built for the *previous* catalog
    while shape_data itself had already been correctly resized, corrupting
    every result via a stale-grid/fresh-data mismatch. Fixed by capturing
    obs_changed's value in a local *before* the reset and using that local
    throughout the flag computations."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    original_val = dcwlf.m2lnL_val(mset)

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    small_galaxies = _GALAXIES[:3]
    small_obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, len(small_galaxies), cols)

    for i, (ra, dec, zp, sigma0, e1, e2, std_noise, c1, c2, m) in enumerate(
        small_galaxies
    ):
        small_obs.set("ra", i, ra)
        small_obs.set("dec", i, dec)
        small_obs.set("z", i, 0.0)
        small_obs.set("zp", i, zp)
        small_obs.set("sigma0", i, sigma0)
        small_obs.set("epsilon_int_1", i, 0.0)
        small_obs.set("epsilon_int_2", i, 0.0)
        small_obs.set("epsilon_obs_1", i, e1)
        small_obs.set("epsilon_obs_2", i, e2)
        small_obs.set("std_noise", i, std_noise)
        small_obs.set("c1", i, c1)
        small_obs.set("c2", i, c2)
        small_obs.set("m", i, m)

    dcwlf.set_obs(small_obs)
    swapped_val = dcwlf.m2lnL_val(mset)

    dcwlf_fresh = Nc.DataClusterWLFactor.new(
        position_factor, redshift_factor, shape_factor
    )
    dcwlf_fresh.set_obs(small_obs)
    dcwlf_fresh.set_prec(1.0e-8)
    dcwlf_fresh.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)
    fresh_val = dcwlf_fresh.m2lnL_val(mset)

    assert_allclose(swapped_val, fresh_val, rtol=1.0e-8)

    # And swapping back to the original catalog must match the original value.
    dcwlf.set_obs(new_obs)
    assert_allclose(dcwlf.m2lnL_val(mset), original_val, rtol=1.0e-8)


def test_m2lnL_gal_parity():
    """Per-galaxy breakdown localizes any mismatch to a single galaxy."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)
    legacy_obs = _build_legacy_obs()

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)

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


def test_prepare_rebuild_is_idempotent():
    """A second prepare() with nothing changed does not alter the result."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)

    first = dcwlf.m2lnL_val(mset)
    second = dcwlf.m2lnL_val(mset)
    assert first == second


def test_prepare_reacts_to_mass_change():
    """Changing log10MDelta between prepare() cycles changes the likelihood."""
    mset, hms = _build_mset()
    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)

    hms.param_set_by_name("log10MDelta", 14.0)
    m2lnL_a = dcwlf.m2lnL_val(mset)

    hms.param_set_by_name("log10MDelta", 14.5)
    m2lnL_b = dcwlf.m2lnL_val(mset)

    assert m2lnL_a != m2lnL_b


def test_auto_nodes_matches_fixed_and_lnint():
    """auto-nodes calibrates a per-galaxy fixed-node configuration instead of
    using the global (n-nodes, rule-n) for every galaxy -- it must still
    agree with the plain FIXED_NODES path and with LNINT to within a
    tolerance consistent with the default node-reltol."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf_auto = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf_auto.set_obs(new_obs)
    dcwlf_auto.set_prec(1.0e-8)
    dcwlf_auto.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)
    dcwlf_auto.set_auto_nodes(True)

    dcwlf_fixed = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf_fixed.set_obs(new_obs)
    dcwlf_fixed.set_prec(1.0e-8)
    dcwlf_fixed.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    dcwlf_lnint = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf_lnint.set_obs(new_obs)
    dcwlf_lnint.set_prec(1.0e-8)
    dcwlf_lnint.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    auto_m2lnL = dcwlf_auto.m2lnL_val(mset)
    fixed_m2lnL = dcwlf_fixed.m2lnL_val(mset)
    lnint_m2lnL = dcwlf_lnint.m2lnL_val(mset)

    assert_allclose(auto_m2lnL, fixed_m2lnL, rtol=1.0e-3)
    assert_allclose(auto_m2lnL, lnint_m2lnL, rtol=1.0e-3)


def test_auto_nodes_mid_run_property_changes_do_not_corrupt_state():
    """auto-nodes/node-reltol/max-total-nodes are orchestrator properties,
    not NcmModel parameters -- changing them bumps no pkey, exactly like
    n-nodes/rule-n (see test_fixed_nodes_correct_after_changing_n_nodes_rule_n_mid_run).
    Toggling them mid-run must not corrupt the per-galaxy fixed-node grid,
    and each configuration must still agree with LNINT."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    dcwlf_lnint = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf_lnint.set_obs(new_obs)
    dcwlf_lnint.set_prec(1.0e-8)
    dcwlf_lnint.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    lnint_m2lnL = dcwlf_lnint.m2lnL_val(mset)

    dcwlf.m2lnL_val(mset)  # first cycle, auto-nodes off

    for auto_nodes, node_reltol, max_total_nodes in (
        (True, 1.0e-3, 2000),
        (True, 1.0e-5, 500),
        (False, 1.0e-4, 2000),
        (True, 1.0e-4, 2000),
    ):
        dcwlf.set_auto_nodes(auto_nodes)
        dcwlf.set_node_reltol(node_reltol)
        dcwlf.set_max_total_nodes(max_total_nodes)

        new_m2lnL = dcwlf.m2lnL_val(mset)

        assert_allclose(new_m2lnL, lnint_m2lnL, rtol=1.0e-2)


def test_cubature_matches_lnint_and_fixed():
    """CUBATURE (an NcmIntegralND cubature over the linear-domain product,
    reusing the same linear-domain integrand twins auto-nodes calibration
    needs) is a third, independently-coded exact scheme -- it must agree
    with both LNINT and FIXED_NODES to a tight tolerance."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf_cub = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf_cub.set_obs(new_obs)
    dcwlf_cub.set_prec(1.0e-8)
    dcwlf_cub.set_integ_method(Nc.DataClusterWLIntegMethod.CUBATURE)

    dcwlf_lnint = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf_lnint.set_obs(new_obs)
    dcwlf_lnint.set_prec(1.0e-8)
    dcwlf_lnint.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    dcwlf_fixed = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf_fixed.set_obs(new_obs)
    dcwlf_fixed.set_prec(1.0e-8)
    dcwlf_fixed.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    cub_m2lnL = dcwlf_cub.m2lnL_val(mset)
    lnint_m2lnL = dcwlf_lnint.m2lnL_val(mset)
    fixed_m2lnL = dcwlf_fixed.m2lnL_val(mset)

    assert_allclose(cub_m2lnL, lnint_m2lnL, rtol=1.0e-6)
    assert_allclose(cub_m2lnL, fixed_m2lnL, rtol=1.0e-6)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
