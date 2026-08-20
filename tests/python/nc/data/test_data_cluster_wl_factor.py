#!/usr/bin/env python
#
# test_data_cluster_wl_factor.py
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

"""Unit tests for ``NcDataClusterWLFactor``'s own accessors.

The parity (``test_data_cluster_wl_factor_parity.py``) and mc-bias tests
already exercise the class's numerical core (``prepare``/``m2lnL_val``)
end-to-end; this file covers what those never touch: the GObject property
get/set round trip and the plain getter/setter/``ref`` wrappers.
"""

import pathlib
import subprocess
import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import duplicate_via_serialization

Ncm.cfg_init()

ZP_MIN, ZP_MAX = 0.0, 5.0
ELLIP_CONV = Nc.GalaxyWLObsEllipConv.TRACE_DET
FRAME = Nc.WLEllipticityFrame.CELESTIAL


def _build_factors_and_obs():
    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(ZP_MIN, ZP_MAX)
    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)

    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hms.param_set_by_name("log10MDelta", 14.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    hp.param_set_by_name("z", 0.2)
    hp.prepare(cosmo)

    pop_shape = Nc.GalaxyShapePopGauss.new()
    pop_shape.param_set_by_name("sigma", 0.3)
    pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    obs_z = Nc.GalaxyRedshiftObsGauss.new()

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop_shape, pop_z, obs_z):
        mset.set(model)
    mset.prepare_fparam_map()

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, 1, cols)
    obs.set("ra", 0, 0.03)
    obs.set("dec", 0, 0.02)
    obs.set("z", 0, 0.0)
    obs.set("zp", 0, 0.6)
    obs.set("sigma0", 0, 0.03)
    obs.set("epsilon_int_1", 0, 0.0)
    obs.set("epsilon_int_2", 0, 0.0)
    obs.set("epsilon_obs_1", 0, 0.05)
    obs.set("epsilon_obs_2", 0, -0.02)
    obs.set("std_noise", 0, 0.03)
    obs.set("c1", 0, 0.0)
    obs.set("c2", 0, 0.0)
    obs.set("m", 0, 0.0)

    return position_factor, redshift_factor, shape_factor, obs


def _build_multi_galaxy_setup(n):
    """N distinct galaxies (varying zp/ellipticity), for bootstrap tests
    where a single, uniform galaxy would make every resample indistinguishable."""
    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(ZP_MIN, ZP_MAX)
    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)

    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hms.param_set_by_name("log10MDelta", 14.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    hp.param_set_by_name("z", 0.2)
    hp.prepare(cosmo)

    pop_shape = Nc.GalaxyShapePopGauss.new()
    pop_shape.param_set_by_name("sigma", 0.3)
    pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    obs_z = Nc.GalaxyRedshiftObsGauss.new()

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop_shape, pop_z, obs_z):
        mset.set(model)
    mset.prepare_fparam_map()

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, n, cols)

    # Distinct galaxies (varying ra/dec/zp/ellipticity), matching
    # test_data_cluster_wl_factor_parity.py's own known-good _GALAXIES
    # values -- a uniform/identical-per-galaxy catalog was found to trip a
    # separate, unrelated LNINT numerical edge case unrelated to bootstrap.
    _GALAXIES = [
        (0.03, 0.02, 0.60, 0.030, 0.05, -0.02, 0.03),
        (-0.10, 0.15, 0.90, 0.040, -0.04, 0.01, 0.05),
        (0.05, -0.08, 0.15, 0.020, 0.02, 0.03, 0.04),
        (0.08, 0.05, 1.50, 0.050, 0.03, -0.01, 0.04),
        (-0.05, -0.03, 0.35, 0.025, -0.01, 0.02, 0.035),
        (0.02, 0.09, 0.75, 0.035, 0.01, -0.03, 0.045),
    ]
    assert n <= len(_GALAXIES)

    for i in range(n):
        ra, dec, zp, sigma0, e1, e2, std_noise = _GALAXIES[i]
        obs.set("ra", i, ra)
        obs.set("dec", i, dec)
        obs.set("z", i, 0.0)
        obs.set("zp", i, zp)
        obs.set("sigma0", i, sigma0)
        obs.set("epsilon_int_1", i, 0.0)
        obs.set("epsilon_int_2", i, 0.0)
        obs.set("epsilon_obs_1", i, e1)
        obs.set("epsilon_obs_2", i, e2)
        obs.set("std_noise", i, std_noise)
        obs.set("c1", i, 0.0)
        obs.set("c2", i, 0.0)
        obs.set("m", i, 0.0)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)

    return dcwlf, mset


@pytest.mark.parametrize(
    "integ_method",
    [
        Nc.DataClusterWLIntegMethod.LNINT,
        Nc.DataClusterWLIntegMethod.FIXED_NODES,
        Nc.DataClusterWLIntegMethod.CUBATURE,
    ],
)
def test_bootstrap_matches_manual_resum(integ_method):
    """Under bootstrap, m2lnL_val() must equal the sum of the per-galaxy
    m2lnP values (from the un-bootstrapped eval_m2lnP_gal) at the
    bootstrap-drawn indices, each counted with its own multiplicity --
    exercised for all three integration methods."""
    n = 6
    dcwlf, mset = _build_multi_galaxy_setup(n)
    dcwlf.set_integ_method(integ_method)

    m2lnP_gal = Ncm.Vector.new(n)
    dcwlf.eval_m2lnP_gal(mset, m2lnP_gal)
    per_galaxy = [m2lnP_gal.get(i) for i in range(n)]

    dcwlf.bootstrap_create()
    assert dcwlf.bootstrap_enabled() is True

    rng = Ncm.RNG.seeded_new(None, 4321)
    dcwlf.bootstrap_resample(rng)
    bstrap = dcwlf.peek_bootstrap()

    expected = sum(per_galaxy[bstrap.get(i)] for i in range(bstrap.get_bsize()))
    got = dcwlf.m2lnL_val(mset)

    assert_allclose(got, expected, rtol=1.0e-10)

    dcwlf.bootstrap_remove()
    assert dcwlf.bootstrap_enabled() is False
    assert_allclose(dcwlf.m2lnL_val(mset), sum(per_galaxy), rtol=1.0e-10)


def test_serialize_deserialize():
    """A round trip through NcmSerialize preserves the construct-only
    factors and every scalar property."""
    position_factor, redshift_factor, shape_factor, obs = _build_factors_and_obs()
    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)
    dcwlf.set_cut(0.1, 3.0)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.CUBATURE)
    dcwlf.set_n_nodes(5)
    dcwlf.set_rule_n(3)
    dcwlf.set_auto_nodes(True)
    dcwlf.set_node_reltol(1.0e-4)
    dcwlf.set_max_total_nodes(100)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dcwlf2 = duplicate_via_serialization(dcwlf, ser)

    assert isinstance(dcwlf2, Nc.DataClusterWLFactor)
    assert dcwlf2 is not dcwlf
    assert dcwlf2.props.len == dcwlf.props.len
    assert dcwlf2.props.r_min == dcwlf.props.r_min
    assert dcwlf2.props.r_max == dcwlf.props.r_max
    assert dcwlf2.get_integ_method() == dcwlf.get_integ_method()
    assert dcwlf2.get_n_nodes() == dcwlf.get_n_nodes()
    assert dcwlf2.get_rule_n() == dcwlf.get_rule_n()
    assert dcwlf2.get_auto_nodes() == dcwlf.get_auto_nodes()
    assert dcwlf2.get_node_reltol() == dcwlf.get_node_reltol()
    assert dcwlf2.get_max_total_nodes() == dcwlf.get_max_total_nodes()


def test_register_shared_anchors_obs_across_dups():
    """register_shared() anchors the obs catalog in the given NcmSerialize
    under a custom name, so it survives ser.reset(True) (autosave-only
    cleanup) and every later dup_obj() through that same ser still
    resolves it to the identical obs object instead of deep-copying it --
    exactly the sequence NcmFitMC's per-worker loop runs (dup, then
    reset(True), for every worker) to share a large read-only per-galaxy
    catalog across parallel workers rather than duplicating it once per
    worker."""
    position_factor, redshift_factor, shape_factor, obs = _build_factors_and_obs()
    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)

    dset = Ncm.Dataset.new()
    dset.append_data(dcwlf)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dset.register_shared(ser)

    dup1 = duplicate_via_serialization(dcwlf, ser)
    ser.reset(True)
    dup2 = duplicate_via_serialization(dcwlf, ser)

    assert dup1 is not dup2
    assert dup1.peek_obs() is dup2.peek_obs()


def test_without_register_shared_obs_is_deep_copied_across_dups():
    """Negative control for the test above: skipping register_shared()
    leaves obs unanchored, so it goes back to being an independent deep
    copy after each dup_obj()+reset(True) pass -- proving the sharing above
    is caused by register_shared(), not some other dup_obj()/reset()
    side effect."""
    position_factor, redshift_factor, shape_factor, obs = _build_factors_and_obs()
    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    dup1 = duplicate_via_serialization(dcwlf, ser)
    ser.reset(True)
    dup2 = duplicate_via_serialization(dcwlf, ser)

    assert dup1 is not dup2
    assert dup1.peek_obs() is not dup2.peek_obs()


def test_construct_only_properties_readable():
    """The construct-only factor properties round-trip through the property
    system (peek_obs/set_obs cover the "obs" property; these three only
    ever go through get_property, since they cannot be re-set)."""
    position_factor, redshift_factor, shape_factor, obs = _build_factors_and_obs()
    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)

    assert dcwlf.props.position_factor == position_factor
    assert dcwlf.props.redshift_factor == redshift_factor
    assert dcwlf.props.shape_factor == shape_factor
    assert dcwlf.props.obs == obs
    assert dcwlf.props.len == 1
    assert dcwlf.peek_obs() == obs
    assert dcwlf.peek_data_array() is not None


def test_obs_property_setter():
    """Setting "obs" through the GObject property system (not set_obs()
    directly) exercises the PROP_OBS case in set_property."""
    position_factor, redshift_factor, shape_factor, obs = _build_factors_and_obs()
    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)

    dcwlf.props.obs = obs

    assert dcwlf.peek_obs() == obs
    assert dcwlf.props.len == 1


def test_scalar_properties_round_trip():
    """Every read/write scalar property round-trips through both the
    property system and its dedicated getter/setter methods."""
    position_factor, redshift_factor, shape_factor, obs = _build_factors_and_obs()
    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)

    dcwlf.props.r_min = 0.5
    dcwlf.props.r_max = 1.5
    dcwlf.props.prec = 1.0e-7
    assert dcwlf.props.r_min == 0.5
    assert dcwlf.props.r_max == 1.5
    assert dcwlf.props.prec == 1.0e-7

    dcwlf.set_cut(0.2, 2.0)
    assert dcwlf.props.r_min == 0.2
    assert dcwlf.props.r_max == 2.0

    dcwlf.props.resample_flag = Nc.DataClusterWLResampleFlag.POSITION
    assert dcwlf.props.resample_flag == Nc.DataClusterWLResampleFlag.POSITION
    assert dcwlf.get_resample_flag() == Nc.DataClusterWLResampleFlag.POSITION
    dcwlf.set_resample_flag(Nc.DataClusterWLResampleFlag.ALL)
    assert dcwlf.props.resample_flag == Nc.DataClusterWLResampleFlag.ALL
    assert dcwlf.get_resample_flag() == Nc.DataClusterWLResampleFlag.ALL

    dcwlf.props.integ_method = Nc.DataClusterWLIntegMethod.CUBATURE
    assert dcwlf.props.integ_method == Nc.DataClusterWLIntegMethod.CUBATURE
    assert dcwlf.get_integ_method() == Nc.DataClusterWLIntegMethod.CUBATURE
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    assert dcwlf.props.integ_method == Nc.DataClusterWLIntegMethod.LNINT
    assert dcwlf.get_integ_method() == Nc.DataClusterWLIntegMethod.LNINT

    dcwlf.props.n_nodes = 5
    assert dcwlf.props.n_nodes == 5
    assert dcwlf.get_n_nodes() == 5
    dcwlf.set_n_nodes(7)
    assert dcwlf.props.n_nodes == 7
    assert dcwlf.get_n_nodes() == 7

    dcwlf.props.rule_n = 3
    assert dcwlf.props.rule_n == 3
    assert dcwlf.get_rule_n() == 3
    dcwlf.set_rule_n(4)
    assert dcwlf.props.rule_n == 4
    assert dcwlf.get_rule_n() == 4

    dcwlf.props.auto_nodes = True
    assert dcwlf.props.auto_nodes is True
    assert dcwlf.get_auto_nodes() is True
    dcwlf.set_auto_nodes(False)
    assert dcwlf.props.auto_nodes is False
    assert dcwlf.get_auto_nodes() is False

    dcwlf.props.node_reltol = 1.0e-3
    assert dcwlf.props.node_reltol == 1.0e-3
    assert dcwlf.get_node_reltol() == 1.0e-3
    dcwlf.set_node_reltol(1.0e-4)
    assert dcwlf.props.node_reltol == 1.0e-4
    assert dcwlf.get_node_reltol() == 1.0e-4

    dcwlf.props.max_total_nodes = 100
    assert dcwlf.props.max_total_nodes == 100
    assert dcwlf.get_max_total_nodes() == 100
    dcwlf.set_max_total_nodes(200)
    assert dcwlf.props.max_total_nodes == 200
    assert dcwlf.get_max_total_nodes() == 200


def test_eval_m2lnP_gal_matches_m2lnL_val():
    """eval_m2lnP_gal fills the per-galaxy vector consistently with the
    aggregate m2lnL_val across every integration method."""
    position_factor, redshift_factor, shape_factor, obs = _build_factors_and_obs()
    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)
    dcwlf.set_prec(1.0e-6)

    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hms.param_set_by_name("log10MDelta", 14.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    hp.param_set_by_name("z", 0.2)
    hp.prepare(cosmo)

    pop_shape = Nc.GalaxyShapePopGauss.new()
    pop_shape.param_set_by_name("sigma", 0.3)
    pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    obs_z = Nc.GalaxyRedshiftObsGauss.new()

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop_shape, pop_z, obs_z):
        mset.set(model)
    mset.prepare_fparam_map()

    for integ_method in (
        Nc.DataClusterWLIntegMethod.LNINT,
        Nc.DataClusterWLIntegMethod.FIXED_NODES,
        Nc.DataClusterWLIntegMethod.CUBATURE,
    ):
        dcwlf.set_integ_method(integ_method)

        m2lnP_gal = Ncm.Vector.new(1)
        dcwlf.eval_m2lnP_gal(mset, m2lnP_gal)
        total = dcwlf.m2lnL_val(mset)

        assert m2lnP_gal.get(0) == total


# ---------------------------------------------------------------------------
# NC_GALAXY_LOW_PROB accounting (get_low_prob_count) and the FIXED_NODES
# background-mass split that keeps P_gal from going negative in the first place.
# ---------------------------------------------------------------------------

_ALL_METHODS = [
    Nc.DataClusterWLIntegMethod.LNINT,
    Nc.DataClusterWLIntegMethod.FIXED_NODES,
    Nc.DataClusterWLIntegMethod.CUBATURE,
]

# Well-separated galaxies with plausible ellipticities: no galaxy's marginal
# should come anywhere near underflowing.
_HEALTHY_GALAXIES = [
    (0.03, 0.02, 0.60, 0.030, 0.05, -0.02, 0.03),
    (-0.10, 0.15, 0.90, 0.040, -0.04, 0.01, 0.05),
    (0.08, 0.05, 1.50, 0.050, 0.03, -0.01, 0.04),
]

# One galaxy whose observed ellipticity is flatly incompatible with any model
# prediction at this noise level, so its marginal underflows to exactly 0 in
# double precision and the NC_GALAXY_LOW_PROB fallback is substituted.
_UNDERFLOWING_GALAXY = (0.03, 0.02, 0.60, 0.030, 0.9, 0.0, 0.01)


def _build_probe_setup(galaxies, *, z_cl=0.2, sigma_int=0.3, log10m=14.0):
    """A catalog built from explicit (ra, dec, zp, sigma0, e1, e2, std_noise)
    rows, so a test can place a galaxy exactly where it needs one -- fully in
    front of the lens, or far out in the tail of the shape likelihood."""
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hms.param_set_by_name("log10MDelta", log10m)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    hp.param_set_by_name("z", z_cl)
    hp.prepare(cosmo)

    pop_shape = Nc.GalaxyShapePopGauss.new()
    pop_shape.param_set_by_name("sigma", sigma_int)
    pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    obs_z = Nc.GalaxyRedshiftObsGauss.new()

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop_shape, pop_z, obs_z):
        mset.set(model)
    mset.prepare_fparam_map()

    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(ZP_MIN, ZP_MAX)
    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, len(galaxies), cols)
    for i, (ra, dec, zp, sigma0, e1, e2, std_noise) in enumerate(galaxies):
        for col in cols:
            obs.set(col, i, 0.0)
        obs.set("ra", i, ra)
        obs.set("dec", i, dec)
        obs.set("zp", i, zp)
        obs.set("sigma0", i, sigma0)
        obs.set("epsilon_obs_1", i, e1)
        obs.set("epsilon_obs_2", i, e2)
        obs.set("std_noise", i, std_noise)

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)
    dcwlf.set_prec(1.0e-8)

    return dcwlf, mset, redshift_factor, obs


@pytest.mark.parametrize("integ_method", _ALL_METHODS)
def test_low_prob_count_zero_on_healthy_eval(integ_method):
    """On a catalog where every galaxy has a comfortably positive marginal, no
    NC_GALAXY_LOW_PROB substitution happens, under any integration method."""
    dcwlf, mset, _, _ = _build_probe_setup(_HEALTHY_GALAXIES)
    dcwlf.set_integ_method(integ_method)

    m2lnL = dcwlf.m2lnL_val(mset)

    assert dcwlf.get_low_prob_count() == 0
    assert m2lnL < 1.0e5  # nowhere near the 1e6 flat fallback


@pytest.mark.parametrize(
    "integ_method",
    [Nc.DataClusterWLIntegMethod.FIXED_NODES, Nc.DataClusterWLIntegMethod.CUBATURE],
)
def test_low_prob_count_counts_underflowing_galaxies(integ_method):
    """When a galaxy's marginal really does evaluate non-positive -- here by
    honest floating-point underflow, an ellipticity no model can produce at
    this noise level -- the flat fallback is substituted AND counted.

    LNINT is excluded on purpose: it does not substitute, it skips the galaxy
    and stores NaN (see test_low_prob_count_resets_between_methods)."""
    dcwlf, mset, _, _ = _build_probe_setup([_UNDERFLOWING_GALAXY], sigma_int=1.0e-2)
    dcwlf.set_integ_method(integ_method)

    m2lnL = dcwlf.m2lnL_val(mset)

    assert dcwlf.get_low_prob_count() >= 1
    assert m2lnL > 1.0e5  # the flat NC_GALAXY_LOW_PROB wall, ~1e6


def test_low_prob_count_resets_between_methods():
    """get_low_prob_count() describes the LAST evaluation, always.

    LNINT never substitutes NC_GALAXY_LOW_PROB, so its answer is legitimately
    zero -- but it must still clear the counter, or it would keep reporting a
    stale count left behind by an earlier FIXED_NODES/CUBATURE evaluation on
    the same object."""
    dcwlf, mset, _, _ = _build_probe_setup([_UNDERFLOWING_GALAXY], sigma_int=1.0e-2)

    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.CUBATURE)
    dcwlf.m2lnL_val(mset)
    assert dcwlf.get_low_prob_count() >= 1  # counter is now dirty

    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    dcwlf.m2lnL_val(mset)
    assert dcwlf.get_low_prob_count() == 0

    # And back again, so the reset is not a one-way latch.
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)
    dcwlf.m2lnL_val(mset)
    assert dcwlf.get_low_prob_count() >= 1


def test_foreground_only_galaxy_has_no_background_quadrature():
    """A galaxy whose whole photo-z support lies in FRONT of the lens gets no
    background quadrature at all (no NcmIntegralFixed, zero background P(z)
    mass), and its marginal is the pure foreground term p_a * norm.

    The support is what decides this, not the point estimate: at z_cl = 0.2 a
    galaxy with zp = 0.15, sigma0 = 0.02 still integrates over [0.0, 0.334] and
    IS partly behind the lens. zp = 0.10, sigma0 = 0.01 gives [0.018, 0.182],
    which genuinely is not."""
    z_cl = 0.2
    foreground = (0.03, 0.02, 0.10, 0.010, 0.05, -0.02, 0.03)

    dcwlf, mset, redshift_factor, obs = _build_probe_setup([foreground], z_cl=z_cl)

    # Verify the premise rather than trusting it.
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    z_data.read_row(obs, 0)
    z_lo, z_hi = redshift_factor.get_integ_lim(mset, z_data)
    assert z_hi <= z_cl, f"galaxy is not fully in the foreground: [{z_lo}, {z_hi}]"

    results = {}
    for integ_method in _ALL_METHODS:
        dcwlf.set_integ_method(integ_method)
        results[integ_method] = dcwlf.m2lnL_val(mset)
        assert dcwlf.get_low_prob_count() == 0

    # The foreground-only path must reach the same answer as the adaptive
    # schemes, which take no such shortcut.
    assert_allclose(
        results[Nc.DataClusterWLIntegMethod.FIXED_NODES],
        results[Nc.DataClusterWLIntegMethod.LNINT],
        rtol=1.0e-6,
    )
    assert_allclose(
        results[Nc.DataClusterWLIntegMethod.CUBATURE],
        results[Nc.DataClusterWLIntegMethod.LNINT],
        rtol=1.0e-6,
    )


def test_foreground_and_background_galaxies_mix():
    """Both branches of the has_bg split in the same catalog and the same
    prepare() -- the per-galaxy arrays (fixed_bg_nodes, fixed_bg_norm,
    n_total_per_galaxy) must stay correctly aligned when entry 0 has no
    background quadrature and entries 1..n do."""
    galaxies = [(0.03, 0.02, 0.10, 0.010, 0.05, -0.02, 0.03), *_HEALTHY_GALAXIES]

    dcwlf, mset, _, _ = _build_probe_setup(galaxies)

    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    lnint_gal = Ncm.Vector.new(len(galaxies))
    dcwlf.eval_m2lnP_gal(mset, lnint_gal)

    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    # Checked in both fixed-node configurations, since they lay the per-galaxy
    # arrays out differently: with auto-nodes off every background galaxy shares
    # the global grid, while with it on each gets its own size -- the case where
    # a misalignment between fixed_bg_nodes/fixed_bg_norm/n_total_per_galaxy and
    # the foreground entry would actually bite. The tolerances differ only
    # because auto-nodes targets node-reltol (1e-2 by default), not because the
    # claim differs.
    for auto_nodes, rtol in ((False, 1.0e-5), (True, 1.0e-2)):
        dcwlf.set_auto_nodes(auto_nodes)
        fixed_gal = Ncm.Vector.new(len(galaxies))
        dcwlf.eval_m2lnP_gal(mset, fixed_gal)

        assert_allclose(fixed_gal.dup_array(), lnint_gal.dup_array(), rtol=rtol)


def test_auto_node_calibration_warns_once_per_prepare():
    """auto-node calibration that cannot reach node-reltol must warn ONCE per
    prepare(), naming how many galaxies missed and the worst error -- not once
    per galaxy.

    That aggregation is the whole point: ncm_integral_fixed_calibrate() is
    called per galaxy, and a realistic cluster produced thousands of
    near-identical warnings, which is neither readable nor actionable. Run in a
    subprocess because g_warning goes to stderr, not through the Python warnings
    machinery.
    """
    script = (
        "from numcosmo_py import Ncm, Nc\n"
        "Ncm.cfg_init()\n"
        "import sys; sys.path.insert(0, %r)\n"
        "from test_data_cluster_wl_factor import _build_probe_setup, _HEALTHY_GALAXIES\n"
        "dcwlf, mset, _, _ = _build_probe_setup(_HEALTHY_GALAXIES)\n"
        "dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)\n"
        "dcwlf.set_auto_nodes(True)\n"
        "dcwlf.set_node_reltol(1.0e-14)\n"  # unreachable
        "dcwlf.set_max_total_nodes(30)\n"  # within a tiny node budget
        "dcwlf.m2lnL_val(mset)\n" % str(pathlib.Path(__file__).parent)
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "auto-node calibration did not reach" in result.stderr

    # Once per prepare(), not once per galaxy: exactly one aggregate line, and
    # it accounts for all three galaxies rather than reporting one of them.
    assert result.stderr.count("auto-node calibration did not reach") == 1
    assert f"for {len(_HEALTHY_GALAXIES)} of {len(_HEALTHY_GALAXIES)} galaxies" in (
        result.stderr
    )

    # And the per-call warning inside ncm_integral_fixed_calibrate() stayed
    # suppressed -- that is what passing a non-NULL relerr_out buys.
    assert "did not reach reltol" not in result.stderr


# A galaxy whose photo-z support runs far past the lens, so nearly all of its
# P(z) mass is background: the fixed-node quadrature of that mass, W_bg, is then
# competing with the exact full-support normalization, and on a coarse grid it
# OVERSHOOTS it.
_MOSTLY_BACKGROUND_GALAXY = (0.03, 0.02, 2.00, 0.300, 0.05, -0.02, 0.03)


@pytest.mark.parametrize("n_nodes,rule_n", [(2, 1), (4, 2), (6, 3)])
def test_foreground_weight_is_clamped_at_zero(n_nodes, rule_n):
    """The foreground weight (norm - W_bg) is clamped at zero, and must be.

    W_bg is a QUADRATURE estimate of the background P(z) mass while norm is the
    EXACT full-support normalization, so on a coarse grid nothing stops W_bg
    from exceeding norm. Measured on this galaxy at (n_nodes, rule_n) = (2, 1):
    norm = 0.1077 against W_bg = 0.2339, a 117% overshoot. Unclamped that would
    make the foreground term p_a * (norm - W_bg) large and NEGATIVE -- the
    second way (besides a signed correction term) that P_gal could go negative
    and trip the NC_GALAXY_LOW_PROB wall.

    The clamp is therefore load-bearing, not decoration: this test asserts the
    positivity guarantee survives grids far too coarse to resolve the
    integrand, which is exactly the claim the split form is making.
    """
    dcwlf, mset, _, _ = _build_probe_setup([_MOSTLY_BACKGROUND_GALAXY])
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)
    dcwlf.set_auto_nodes(False)
    dcwlf.set_n_nodes(n_nodes)
    dcwlf.set_rule_n(rule_n)

    m2lnL = dcwlf.m2lnL_val(mset)

    assert dcwlf.get_low_prob_count() == 0
    assert np.isfinite(m2lnL)
    assert m2lnL < 1.0e5  # not the flat 1e6 fallback
