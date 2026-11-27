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
configured with ``LNINT``, so both sides
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

import numpy as np
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


def _build_new_obs(mset, scale=1.0):
    """New-engine Factor objects + one NcGalaxyWLObs built via required_columns.

    @scale multiplies every galaxy's angular offset from the lens centre, so a
    test can move the whole sample to smaller radii (and hence larger reduced
    shear) without inventing a second catalog. scale=1.0 reproduces _GALAXIES
    exactly, which is what every frozen value in this file was captured at.
    """
    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(ZP_MIN, ZP_MAX)
    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, len(_GALAXIES), cols)

    for i, (ra, dec, zp, sigma0, e1, e2, std_noise, c1, c2, m) in enumerate(_GALAXIES):
        obs.set("ra", i, ra * scale)
        obs.set("dec", i, dec * scale)
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


CHEB = Nc.DataClusterWLIntegMethod.CHEB_WEIGHTS
LNINT = Nc.DataClusterWLIntegMethod.LNINT


def _dataset(factors, obs, method):
    """A dataset on @obs, with the position, redshift and shape factors @factors, a
    precision of 1e-8 and the integration method @method."""
    dcwlf = Nc.DataClusterWLFactor.new(*factors)
    dcwlf.set_obs(obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(method)
    return dcwlf


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


def test_m2lnL_cheb_matches_lnint():
    """CHEB_WEIGHTS and LNINT, two independently coded schemes for the same integral,
    agree within the adaptive tolerance."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)
    *factors, obs = _build_new_obs(mset)

    assert_allclose(
        _dataset(factors, obs, CHEB).m2lnL_val(mset),
        _dataset(factors, obs, LNINT).m2lnL_val(mset),
        rtol=1.0e-6,
    )


def _assert_fresh(dcwlf, factors, obs, mset):
    """The value of @dcwlf equals that of a dataset built at the current point."""
    fresh = _dataset(factors, obs, dcwlf.get_integ_method())
    assert_allclose(dcwlf.m2lnL_val(mset), fresh.m2lnL_val(mset), rtol=1.0e-12)


def test_cheb_cache_consistency_across_revisits():
    """A sequence that revisits mass and lens-redshift values gives at each point the
    value of a dataset built there: the grid and the node caches follow z_cl, the
    mass-dependent caches follow the mass."""
    mset, hms = _build_mset()
    *factors, obs = _build_new_obs(mset)
    dcwlf = _dataset(factors, obs, CHEB)
    hp = mset.peek(Nc.HaloPosition.id())
    cosmo = mset.peek(Nc.HICosmo.id())

    for log10_mdelta in (13.5, 14.5, 13.5, 15.0, 14.5):
        hms.param_set_by_name("log10MDelta", log10_mdelta)
        _assert_fresh(dcwlf, factors, obs, mset)

    hms.param_set_by_name("log10MDelta", 14.0)
    for z_cl in (0.2, 0.35, 0.2, 0.5):
        hp.param_set_by_name("z", z_cl)
        hp.prepare(cosmo)
        _assert_fresh(dcwlf, factors, obs, mset)


def test_cheb_correct_under_angular_only_changes():
    """Moving the halo center (RA, not z) keeps the grid, which depends on z_cl only,
    and gives at each point the value of a dataset built there."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)
    *factors, obs = _build_new_obs(mset)
    dcwlf = _dataset(factors, obs, CHEB)
    hp = mset.peek(Nc.HaloPosition.id())
    cosmo = mset.peek(Nc.HICosmo.id())

    for ra in (0.0, 0.05, -0.03, 0.05):
        hp.param_set_by_name("ra", ra)
        hp.prepare(cosmo)
        _assert_fresh(dcwlf, factors, obs, mset)


def test_cheb_correct_after_switching_integ_method_mid_run():
    """LNINT, then CHEB_WEIGHTS on the same dataset with no model change between them:
    the switch takes the grid and the node caches again, so the value is that of a
    dataset built under CHEB_WEIGHTS; switching back gives the first LNINT value."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)
    *factors, obs = _build_new_obs(mset)
    dcwlf = _dataset(factors, obs, LNINT)

    lnint = dcwlf.m2lnL_val(mset)
    dcwlf.m2lnL_val(mset)

    dcwlf.set_integ_method(CHEB)
    _assert_fresh(dcwlf, factors, obs, mset)

    dcwlf.set_integ_method(LNINT)
    assert_allclose(dcwlf.m2lnL_val(mset), lnint, rtol=1.0e-12)


def test_cheb_resample_reuse_matches_lnint():
    """One dataset under CHEB_WEIGHTS, set before the first resample(), agrees with LNINT
    on every one of several resample() cycles. ncm_data_resample() calls prepare()
    before the resample virtual; resample() marks the catalog as changed, so the next
    prepare() takes the grid of the new rows."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)
    *factors, obs = _build_new_obs(mset)
    dcwlf_reused = _dataset(factors, obs, CHEB)
    dcwlf_lnint = _dataset(factors, obs, LNINT)

    for seed in (42, 43, 44):
        dcwlf_reused.resample(mset, Ncm.RNG.seeded_new(None, seed))
        dcwlf_lnint.resample(mset, Ncm.RNG.seeded_new(None, seed))

        assert_allclose(
            dcwlf_reused.m2lnL_val(mset), dcwlf_lnint.m2lnL_val(mset), rtol=1.0e-5
        )


# Frozen legacy resample() output, keyed by seed (see module docstring for
# provenance). For seeds 100/101/102, legacy's regenerated ra/dec/zp/
# epsilon_obs_1/epsilon_obs_2 columns were bit-identical (rtol=0, atol=0) to
# the new engine's, and the resulting -2lnL and per-galaxy breakdown agreed
# to rtol=1e-5 -- both are frozen here at the values legacy actually
# produced.
#
# The epsilon_obs_1/2 columns are re-frozen from the current engine
# (2026-09-25): ncm_trivec_get_spherical_coord() now takes the polar angle
# from atan2 instead of acos, which lost up to 8e-10 relative accuracy for
# galaxies near the halo axis, and ncm_quaternion_set_to_rotate_to_z() was
# rewritten with atan2 half-angles. The ellipticities moved by at most 4e-11
# relative; ra/dec/zp and the -2lnL checks still hold the legacy values.
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
                -0.3172410411890312,
                -0.2049826215795677,
                -0.3041223958451971,
                -0.16957200003045658,
                0.29144924808360806,
            ],
            "epsilon_obs_2": [
                -0.05560002156554589,
                -0.06113458768892292,
                0.09135402976705671,
                -0.4769546736227781,
                -0.3426895589182039,
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
                0.22406866894888364,
                0.061207856101618996,
                0.44372863619665204,
                0.08717573236579043,
                -0.02348215563242228,
            ],
            "epsilon_obs_2": [
                -0.10874634715112702,
                0.15859618948467746,
                -0.048557488147717584,
                0.4504721421956441,
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
                -0.29682870959192836,
                -0.166700407400143,
                -0.13042718475162582,
                -0.10740386699619837,
                0.23381410171989778,
            ],
            "epsilon_obs_2": [
                -0.21837640908751518,
                0.08850707991243387,
                0.3235980182479202,
                0.08241525800189359,
                -0.36964458619803503,
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
            # for redshift, ncm_rng_gaussian_gen() for shape noise).
            assert_allclose(new_vals, frozen["cols"][col], rtol=1.0e-12, atol=0.0)

        new_m2lnL = dcwlf.m2lnL_val(mset)
        assert_allclose(new_m2lnL, frozen["m2lnL"], rtol=1.0e-5)

        new_gal = Ncm.Vector.new(n)
        dcwlf.eval_m2lnP_gal(mset, new_gal)
        assert_allclose(new_gal.dup_array(), frozen["gal"], rtol=1.0e-5)


def test_cheb_correct_after_swapping_obs_to_different_sized_catalog():
    """set_obs() to a catalog of another galaxy count, with no model change, takes the
    grid of the new rows: the value is that of a dataset built on that catalog, and
    swapping back gives the first value."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)
    *factors, new_obs = _build_new_obs(mset)
    position_factor, redshift_factor, shape_factor = factors
    dcwlf = _dataset(factors, new_obs, CHEB)

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
    _assert_fresh(dcwlf, factors, small_obs, mset)

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


# Radial squeeze applied to the sample below. Pulling every galaxy 33x closer to
# the lens centre is what produces reduced shears large enough to expose the
# pathology this test guards; at the unscaled radii the old code was fine.
_WALL_SCALE = 0.03
_WALL_LOG10M_GRID = (15.0, 15.2, 15.4, 15.6, 15.8, 16.0, 16.2)


def test_no_low_prob_wall_at_high_mass():
    """At reduced shears near the critical curve (every galaxy pulled 33 times closer to
    the lens), CHEB_WEIGHTS gives a finite -2lnL with no NC_GALAXY_LOW_PROB substitution
    and no jump across masses. A substitution adds about 1e6 to -2lnL, a flat wall that
    confines an optimizer."""
    mset, hms = _build_mset()
    *factors, new_obs = _build_new_obs(mset, scale=_WALL_SCALE)
    dcwlf = _dataset(factors, new_obs, CHEB)

    values = []
    for log10_mdelta in _WALL_LOG10M_GRID:
        hms.param_set_by_name("log10MDelta", log10_mdelta)
        m2lnL = dcwlf.m2lnL_val(mset)

        assert (
            dcwlf.get_low_prob_count() == 0
        ), f"NC_GALAXY_LOW_PROB substituted at log10MDelta={log10_mdelta}"
        assert np.isfinite(m2lnL)
        # Nowhere near the 1e6 flat fallback -- one galaxy hitting it would
        # dominate this 5-galaxy total by five orders of magnitude.
        assert abs(m2lnL) < 1.0e4
        values.append(m2lnL)

    # No isolated spike: the wall showed up as a single mass step jumping by
    # ~1e6 and coming straight back down.
    steps = np.abs(np.diff(values))
    assert steps.max() < 1.0e3, f"discontinuous -2lnL across the mass grid: {values}"


def test_cubature_matches_lnint_and_cheb():
    """CUBATURE, an NcmIntegralND cubature over the linear-domain product, agrees with
    LNINT and CHEB_WEIGHTS to a tight tolerance."""
    mset, hms = _build_mset()
    hms.param_set_by_name("log10MDelta", 14.0)
    *factors, obs = _build_new_obs(mset)

    cub_m2lnL = _dataset(factors, obs, Nc.DataClusterWLIntegMethod.CUBATURE).m2lnL_val(
        mset
    )

    assert_allclose(
        cub_m2lnL, _dataset(factors, obs, LNINT).m2lnL_val(mset), rtol=1.0e-6
    )
    assert_allclose(
        cub_m2lnL, _dataset(factors, obs, CHEB).m2lnL_val(mset), rtol=1.0e-6
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
