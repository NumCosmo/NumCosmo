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
split, the log-sum-exp panel recombination) -- originally against legacy
configured with ``LNINT`` (not its default ``FIXED_NODES``), so both sides
were independently coded adaptive log-domain 1D integrators, not merely
re-running the same code.

FROZEN REFERENCE VALUES: the dual-engine agreement documented above was
proven by running both engines live and is captured, not re-derived, in the
functions below. Values were captured from an actual passing run of this
file's original legacy-comparison code, at git rev ``77313f22``
(2026-07-16), then legacy (``NcDataClusterWL``) construction was removed so
these tests no longer depend on legacy at runtime -- legacy is slated for
deletion in a follow-up PR. Each frozen assertion keeps the tolerance
(``rtol``/``atol``) that the original live comparison used, so the same
numerical claim ("new engine matches legacy to within adaptive-quadrature
tolerance") is still being checked, just against a constant instead of a
second live engine.
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
    """One shared mset: lens-side models plus the new engine's own
    galaxy-population-side models (peeked by model ID, exactly like
    NcDataClusterWL peeks NcGalaxySDShape/NcGalaxySDObsRedshift/
    NcGalaxySDPosition; see nc_data_cluster_wl.c:1017-1019, historically --
    legacy is no longer built in this file, see module docstring)."""
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

    mset = Ncm.MSet.empty_new()
    for model in (
        cosmo,
        dp,
        hp,
        smd,
        pop_shape,
        pop_z,
        obs_z,
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


# Frozen legacy (LNINT) -2lnL values, keyed by log10_mdelta -- see module
# docstring for provenance.
_M2LNL_LNINT_FROZEN = {
    13.5: -18.165656396467302,
    14.0: -18.160408926654355,
    14.5: -18.144406481325497,
    15.0: -18.094430558394198,
}


@pytest.mark.parametrize("log10_mdelta", _LOG10_MDELTA_GRID)
def test_m2lnL_parity(log10_mdelta):
    """NcDataClusterWLFactor's total -2lnL matches legacy's (LNINT) across
    masses -- checked here against frozen legacy values (see module
    docstring)."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", log10_mdelta)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)

    new_m2lnL = dcwlf.m2lnL_val(mset)

    assert_allclose(new_m2lnL, _M2LNL_LNINT_FROZEN[log10_mdelta], rtol=1.0e-5)


# Frozen legacy (default FIXED_NODES) -2lnL values, keyed by log10_mdelta.
_M2LNL_FIXED_NODES_FROZEN = {
    13.5: -18.165656396162543,
    14.0: -18.160408925652295,
    14.5: -18.144406477942447,
    15.0: -18.094430546800766,
}


@pytest.mark.parametrize("log10_mdelta", _LOG10_MDELTA_GRID)
def test_m2lnL_parity_fixed_nodes(log10_mdelta):
    """NcDataClusterWLFactor's FIXED_NODES matches legacy's own (default)
    FIXED_NODES exactly -- both sides call the same underlying VarAdd-style
    closed-form eval_at_nodes/make_fixed_nodes machinery, so this checked the
    orchestrator-level wiring (grid construction, control-variate combine),
    not a numerically-independent scheme like the LNINT comparison above.
    Now checked against frozen legacy values (see module docstring)."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", log10_mdelta)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    new_m2lnL = dcwlf.m2lnL_val(mset)

    assert_allclose(new_m2lnL, _M2LNL_FIXED_NODES_FROZEN[log10_mdelta], rtol=1.0e-8)


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


# Frozen legacy FIXED_NODES -2lnL values for the revisit sequences below (see
# module docstring for provenance). Lists align positionally with the
# (log10_mdelta,)/(z_cl,) sequences iterated in the test.
_CACHE_REVISIT_MASS_SEQ = (13.5, 14.5, 13.5, 15.0, 14.5)
_CACHE_REVISIT_MASS_FROZEN = [
    -18.165656396162543,
    -18.144406477942447,
    -18.165656396162543,
    -18.094430546800766,
    -18.144406477942447,
]
_CACHE_REVISIT_ZCL_SEQ = (0.2, 0.35, 0.2, 0.5)
_CACHE_REVISIT_ZCL_FROZEN = [
    -18.160408925652295,
    -18.16466857347916,
    -18.160408925652295,
    -18.16670389900757,
]


def test_fixed_nodes_cache_consistency_across_revisits():
    """FIXED_NODES matches legacy across a sequence that revisits mass and
    lens-redshift values -- the scenario that catches stale per-galaxy
    fixed-node-grid caches (see submodel-pkey-does-not-propagate lessons).
    Now checked against frozen legacy values (see module docstring)."""
    mset, hms = _build_mset()
    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    for log10_mdelta, frozen in zip(_CACHE_REVISIT_MASS_SEQ, _CACHE_REVISIT_MASS_FROZEN):
        hms.param_set_by_name("log10MDelta", log10_mdelta)
        assert_allclose(dcwlf.m2lnL_val(mset), frozen, rtol=1.0e-8)

    hms.param_set_by_name("log10MDelta", 14.0)
    hp = mset.peek(Nc.HaloPosition.id())
    cosmo = mset.peek(Nc.HICosmo.id())

    for z_cl, frozen in zip(_CACHE_REVISIT_ZCL_SEQ, _CACHE_REVISIT_ZCL_FROZEN):
        hp.param_set_by_name("z", z_cl)
        hp.prepare(cosmo)
        assert_allclose(dcwlf.m2lnL_val(mset), frozen, rtol=1.0e-8)


# Frozen legacy FIXED_NODES -2lnL values for the RA-only-change sequence
# below (see module docstring for provenance); aligned positionally with the
# `ra` sequence iterated in the test.
_ANGULAR_ONLY_RA_SEQ = (0.0, 0.05, -0.03, 0.05)
_ANGULAR_ONLY_FROZEN = [
    -18.160408925652295,
    -18.16212591226027,
    -18.150993979018196,
    -18.16212591226027,
]


def test_fixed_nodes_correct_under_angular_only_changes():
    """Changing halo_position's RA/Dec (not z) must not corrupt the
    fixed-node grid: the grid depends only on z_cl, but halo_position's
    whole-model hash (used elsewhere for radius_changed) bumps on *any* of
    its parameters, including RA/Dec -- e.g. in a miscentering model. This
    guards against conflating "halo_position changed" with "z_cl changed".
    Now checked against frozen legacy values (see module docstring)."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    hp = mset.peek(Nc.HaloPosition.id())
    cosmo = mset.peek(Nc.HICosmo.id())

    for ra, frozen in zip(_ANGULAR_ONLY_RA_SEQ, _ANGULAR_ONLY_FROZEN):
        hp.param_set_by_name("ra", ra)
        hp.prepare(cosmo)
        assert_allclose(dcwlf.m2lnL_val(mset), frozen, rtol=1.0e-8)


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

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    # Prime LNINT state -- radius_hash/optzs_hash become current here, with
    # FIXED_NODES-only caches (grid/crit/sigma) never touched.
    dcwlf.m2lnL_val(mset)
    dcwlf.m2lnL_val(mset)  # idempotent second call, nothing changed

    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)
    new_m2lnL = dcwlf.m2lnL_val(mset)

    # Frozen legacy FIXED_NODES value at log10MDelta=14.0, z_cl=0.2 (matches
    # _M2LNL_FIXED_NODES_FROZEN[14.0]; see module docstring for provenance).
    frozen_fixed_nodes = -18.160408925652295
    assert_allclose(new_m2lnL, frozen_fixed_nodes, rtol=1.0e-8)

    # And switching back to LNINT must still match the value just obtained
    # under FIXED_NODES (loose rtol -- the two schemes agree closely but are
    # not identical; this mirrors the original live comparison, which
    # checked the post-switch LNINT value against the just-computed
    # FIXED_NODES value, not an independently frozen LNINT constant).
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    assert_allclose(dcwlf.m2lnL_val(mset), frozen_fixed_nodes, rtol=1.0e-5)


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


# Frozen legacy resample() output, keyed by seed (see module docstring for
# provenance). For seeds 100/101/102, legacy's regenerated ra/dec/zp/
# epsilon_obs_1/epsilon_obs_2 columns were bit-identical (rtol=0, atol=0) to
# the new engine's, and the resulting -2lnL and per-galaxy breakdown agreed
# to rtol=1e-5 -- both are frozen here at the values legacy actually
# produced.
_RESAMPLE_FROZEN = {
    100: {
        "cols": {
            "ra": [
                0.01055298587307335,
                0.0497247557155788,
                0.05295343436300755,
                0.05867658639326691,
                -0.023098387103527795,
            ],
            "dec": [
                -0.03019290320912554,
                0.1565286570695748,
                -0.09037035664249453,
                -0.09902930641926685,
                -0.18540952509591038,
            ],
            "zp": [
                0.4777536975753619,
                0.7638973337747148,
                0.6600705768089372,
                0.9701997732398961,
                0.49017092510689864,
            ],
            "epsilon_obs_1": [
                -0.3172410411898537,
                -0.20498262157952554,
                -0.30412239584524453,
                -0.16957200003042172,
                0.2914492480836165,
            ],
            "epsilon_obs_2": [
                -0.055600021566330905,
                -0.06113458768895565,
                0.09135402976695699,
                -0.47695467362272886,
                -0.34268955891820446,
            ],
        },
        "m2lnL": -14.093228702890663,
        "gal": [
            -3.557916867979317,
            -3.7399762475326623,
            -3.431647383974465,
            -0.7374323566877973,
            -2.6262558467164214,
        ],
    },
    101: {
        "cols": {
            "ra": [
                0.193779047485441,
                -0.16972724087536334,
                0.044639803748577844,
                0.04105259515345097,
                0.07187970345839857,
            ],
            "dec": [
                -0.18861026707202516,
                0.02169099315840281,
                0.09159695642327917,
                0.19772715238506186,
                -0.05711015623362085,
            ],
            "zp": [
                0.3177101960527326,
                0.2887145181110548,
                0.6068695443468395,
                0.8006118456521737,
                0.08148659910639647,
            ],
            "epsilon_obs_1": [
                0.22406866894888353,
                0.06120785610161111,
                0.44372863619670483,
                0.08717573236579651,
                -0.02348215563242228,
            ],
            "epsilon_obs_2": [
                -0.10874634715113454,
                0.1585961894846795,
                -0.048557488147818476,
                0.45047214219564263,
                -0.09709542061085731,
            ],
        },
        "m2lnL": -14.562186081501187,
        "gal": [
            -4.065729532204247,
            -4.271628487741558,
            -2.5175692245590118,
            -1.8309841961002156,
            -1.8762746408961548,
        ],
    },
    102: {
        "cols": {
            "ra": [
                0.09229610580950975,
                -0.003920704964548355,
                -0.18824651734903458,
                -0.03952054223045706,
                -0.1716628957539797,
            ],
            "dec": [
                -0.07591722087144458,
                -0.10741032766193949,
                -0.04850003632055123,
                -0.051722427234190994,
                0.015806831370924382,
            ],
            "zp": [
                0.4930734468193961,
                0.5439188689953817,
                0.24344587800819395,
                0.36697771392115247,
                0.5692915550930165,
            ],
            "epsilon_obs_1": [
                -0.2968287095919376,
                -0.16670040740001915,
                -0.1304271847516241,
                -0.10740386699634384,
                0.2338141017198752,
            ],
            "epsilon_obs_2": [
                -0.21837640908748684,
                0.088507079912425,
                0.3235980182479213,
                0.08241525800245833,
                -0.36964458619803375,
            ],
        },
        "m2lnL": -17.858798273307382,
        "gal": [
            -3.2956756702660766,
            -4.342752966909957,
            -3.061495558486798,
            -4.545358328268286,
            -2.613515749376264,
        ],
    },
}


def test_resample_matches_legacy():
    """``resample()`` itself -- not just static ``m2lnL_val`` on a fixed
    catalog -- must generate the *same* realizations as legacy, seed for
    seed. Both engines called gen() in the identical order (redshift,
    position, shape -- see ``_nc_data_cluster_wl_factor_resample`` and
    ``_nc_data_cluster_wl_resample``), and each individual gen() step was
    already proven bit-parity in isolation (shape:
    ``test_galaxy_shape_factor_var_add.py::test_gen_parity_legacy``); this
    was the first check that the *whole* per-galaxy resample pipeline
    (redshift z draw -> position ra/dec draw+radius rejection -> shape
    draw, all sharing one RNG stream) reproduced legacy bit-for-bit once
    composed together in the orchestrator. Now checks the raw regenerated
    columns (ra/dec/zp/epsilon_obs_1/2) and the resulting -2lnL against
    frozen legacy values (see module docstring)."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)

    n = len(_GALAXIES)

    for seed in (100, 101, 102):
        rng_new = Ncm.RNG.seeded_new(None, seed)

        dcwlf.resample(mset, rng_new)

        frozen = _RESAMPLE_FROZEN[seed]

        for col in ("ra", "dec", "zp", "epsilon_obs_1", "epsilon_obs_2"):
            new_vals = [new_obs.get(col, i) for i in range(n)]
            # rtol=1e-12 (not bit-exact): these columns are all RNG draws
            # that route through platform/libm-sensitive transcendentals
            # (asin() in the sky-footprint sampler, an ODE-spline inverse-CDF
            # for redshift, ncm_rng_gaussian_gen() for shape noise) -- see
            # the same reasoning in test_galaxy_redshift_factor_spline_legacy_parity.py.
            assert_allclose(new_vals, frozen["cols"][col], rtol=1.0e-12, atol=0.0)

        new_m2lnL = dcwlf.m2lnL_val(mset)
        assert_allclose(new_m2lnL, frozen["m2lnL"], rtol=1.0e-5)

        new_gal = Ncm.Vector.new(n)
        dcwlf.eval_m2lnP_gal(mset, new_gal)
        assert_allclose(new_gal.dup_array(), frozen["gal"], rtol=1.0e-5)


# Frozen legacy FIXED_NODES -2lnL values, keyed by (n_nodes, rule_n) (see
# module docstring for provenance).
_N_NODES_RULE_N_FROZEN = {
    (6, 3): -18.16036639525873,
    (12, 7): -18.160408926654412,
    (2, 1): -18.130207439694868,
    (10, 5): -18.160408925652295,
}


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
    folding them into fixed_grid_changed. Now checked against frozen legacy
    values (see module docstring)."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

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

        assert_allclose(
            new_m2lnL, _N_NODES_RULE_N_FROZEN[(n_nodes, rule_n)], rtol=1.0e-8
        )


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


# Frozen legacy per-galaxy -2lnP breakdown at log10MDelta=14.0, z_cl=0.2,
# LNINT (see module docstring for provenance).
_GAL_PARITY_FROZEN = [
    -4.652123769387357,
    -3.798160022033027,
    -3.473934061131708,
    -1.4650061634643785,
    -4.771184910637885,
]


def test_m2lnL_gal_parity():
    """Per-galaxy breakdown localizes any mismatch to a single galaxy. Now
    checked against a frozen legacy breakdown (see module docstring)."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(mset)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)

    n = len(_GALAXIES)
    new_gal = Ncm.Vector.new(n)

    dcwlf.eval_m2lnP_gal(mset, new_gal)

    assert_allclose(new_gal.dup_array(), _GAL_PARITY_FROZEN, rtol=1.0e-5)


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
