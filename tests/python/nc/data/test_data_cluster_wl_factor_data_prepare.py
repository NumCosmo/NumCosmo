#
# test_data_cluster_wl_factor_data_prepare.py
#
# Wed Sep 23 2026
# Copyright  2026  Caio Lima de Oliveira
# <caiolimadeoliveira@pm.me>
#
# test_data_cluster_wl_factor_data_prepare.py
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""nc_data_cluster_wl_factor_data_prepare() prepares galaxies in parallel for
every position, redshift and shape factor, and the result must not depend on
the thread count.

Each run builds everything from scratch -- models, catalog, dataset -- so the
models' lazy first-call state (halo-position rotation, LSST-SRD constants, the
numerically integrated Einasto profile's splines, Cuba's settings) starts cold
and the orchestrator's serial warm-up is exercised every time. The catalog
mixes background galaxies with foreground-only ones, so the warm-up has to
skip past galaxies without a redshift grid.
"""

import os

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

ELLIP_CONV = Nc.GalaxyWLObsEllipConv.TRACE_DET
FRAME = Nc.WLEllipticityFrame.CELESTIAL
Z_CL = 0.2
N_GAL = 40
PARALLEL_THREADS = 4

SHAPE_FACTORS = {
    "var-add": lambda: Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV),
    "cgf": lambda: Nc.GalaxyShapeFactorCGF.new(ELLIP_CONV),
    "laplace": lambda: Nc.GalaxyShapeFactorLaplace.new(ELLIP_CONV),
    "series-lensed": lambda: Nc.GalaxyShapeFactorSeriesLensed.new(ELLIP_CONV, 2),
    "fixed-quad": lambda: Nc.GalaxyShapeFactorFixedQuad.new(ELLIP_CONV),
    "quad": lambda: _quad(),
    "moments-tilt": lambda: Nc.GalaxyShapeFactorMomentsTilt.new(ELLIP_CONV),
    "moments-gauss": lambda: Nc.GalaxyShapeFactorMomentsGauss.new(ELLIP_CONV),
}


def _quad():
    """Quad runs a Cuba cubature per evaluation and the calibration needs
    thousands of them per galaxy: a loose tolerance and a short catalog (see
    N_GAL_OF) keep it a unit test. Accuracy is irrelevant here, only that
    parallel equals serial."""
    gsf = Nc.GalaxyShapeFactorQuad.new(ELLIP_CONV)
    gsf.set_reltol(1.0e-2)
    return gsf


#: Galaxies per shape factor; the rest use N_GAL. Still more than the three
#: foreground-only galaxies plus one warm-up galaxy plus PARALLEL_THREADS.
N_GAL_OF = {"quad": 8}

PROFILES = {
    "nfw": Nc.HaloDensityProfileNFW.new,
    "einasto": Nc.HaloDensityProfileEinasto.new,
}


def _galaxies(n_gal):
    """A fixed synthetic catalog: some galaxies entirely in front of the
    cluster (zp well below Z_CL), the rest behind it."""
    rng = np.random.default_rng(20260923)
    ra = rng.uniform(-0.15, 0.15, n_gal)
    dec = rng.uniform(-0.15, 0.15, n_gal)
    zp = rng.uniform(0.25, 1.5, n_gal)
    zp[:3] = 0.05  # foreground-only first: the warm-up must look past them
    sigma0 = rng.uniform(0.02, 0.05, n_gal)
    e1, e2 = rng.normal(0.0, 0.25, (2, n_gal)).clip(-0.7, 0.7)
    std_noise = rng.uniform(0.02, 0.3, n_gal)
    return ra, dec, zp, sigma0, e1, e2, std_noise


def _gauss_pz_spline(zp, sigma):
    z = np.linspace(max(1.0e-3, zp - 6.0 * sigma), zp + 6.0 * sigma, 64)
    pz = np.exp(-0.5 * ((z - zp) / sigma) ** 2) / (np.sqrt(2.0 * np.pi) * sigma)
    spline = Ncm.SplineCubicNotaknot.new()
    spline.set(
        Ncm.Vector.new_array(z.tolist()), Ncm.Vector.new_array(pz.tolist()), True
    )
    return spline


def _build(shape_name, profile_name, redshift_name, n_gal=None):
    """Builds a fresh FIXED_NODES + auto-nodes dataset; returns (dcwlf, mset)."""
    if n_gal is None:
        n_gal = N_GAL_OF.get(shape_name, N_GAL)
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hms.param_set_by_name("log10MDelta", 14.5)
    dp = PROFILES[profile_name](hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    hp.param_set_by_name("z", Z_CL)
    hp.prepare(cosmo)
    pop_shape = Nc.GalaxyShapePopGauss.new()
    pop_shape.param_set_by_name("sigma", 0.3)

    models = [cosmo, dp, hp, smd, pop_shape]
    if redshift_name == "composed":
        redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(0.0, 5.0)
        models += [
            Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source(),
            Nc.GalaxyRedshiftObsGauss.new(),
        ]
    else:
        redshift_factor = Nc.GalaxyRedshiftFactorSpline.new()

    mset = Ncm.MSet.empty_new()
    for model in models:
        mset.set(model)
    mset.prepare_fparam_map()

    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    shape_factor = SHAPE_FACTORS[shape_name]()

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, n_gal, cols)
    ra, dec, zp, sigma0, e1, e2, std_noise = _galaxies(n_gal)
    values = {
        "ra": ra,
        "dec": dec,
        "z": np.zeros(n_gal),
        "zp": zp,
        "sigma0": sigma0,
        "epsilon_int_1": np.zeros(n_gal),
        "epsilon_int_2": np.zeros(n_gal),
        "epsilon_obs_1": e1,
        "epsilon_obs_2": e2,
        "std_noise": std_noise,
        "c1": np.zeros(n_gal),
        "c2": np.zeros(n_gal),
        "m": np.zeros(n_gal),
    }
    for i in range(n_gal):
        for col in cols:
            obs.set(col, i, float(values[col][i]))
        if redshift_name == "spline":
            # One spline object per row, never shared between rows.
            obs.set_pz(i, _gauss_pz_spline(zp[i], sigma0[i] * (1.0 + zp[i])))

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)
    dcwlf.set_auto_nodes(True)
    dcwlf.set_obs(obs)

    return dcwlf, mset


def _m2lnP(dcwlf, mset):
    m2lnP = Ncm.Vector.new(dcwlf.peek_obs().len())
    dcwlf.eval_m2lnP_gal(mset, m2lnP)
    return np.array(m2lnP.dup_array())


def _node_config(dcwlf):
    node_config = dcwlf.props.node_config
    n_nodes = node_config.get_variant("n-nodes")[1].unpack()
    rule_n = node_config.get_variant("rule-n")[1].unpack()
    return n_nodes, rule_n


def _prepare(shape_name, profile_name, redshift_name, n_threads):
    """Builds a fresh dataset, runs data_prepare() on @n_threads threads and
    returns (node-config n-nodes, node-config rule-n, per-galaxy -2lnP)."""
    dcwlf, mset = _build(shape_name, profile_name, redshift_name)

    Ncm.cfg_set_openmp_nthreads(n_threads)
    try:
        dcwlf.data_prepare(mset)
    finally:
        Ncm.cfg_set_openmp_nthreads(
            int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))
        )

    return (*_node_config(dcwlf), _m2lnP(dcwlf, mset))


def _combinations():
    """Every shape factor on the path that reaches the most shared lazy state
    (Composed + LSST-SRD, numerically integrated Einasto), plus the other
    redshift factor and profile for a cheap shape factor."""
    for shape in SHAPE_FACTORS:
        yield shape, "einasto", "composed"
    yield "var-add", "nfw", "composed"
    yield "var-add", "nfw", "spline"
    yield "var-add", "einasto", "spline"


@pytest.mark.parametrize("shape_name,profile_name,redshift_name", list(_combinations()))
def test_parallel_data_prepare_matches_serial(shape_name, profile_name, redshift_name):
    serial = _prepare(shape_name, profile_name, redshift_name, 1)
    parallel = _prepare(shape_name, profile_name, redshift_name, PARALLEL_THREADS)

    assert serial[0] == parallel[0]
    assert serial[1] == parallel[1]
    assert np.all(np.isfinite(serial[2]))
    np.testing.assert_array_equal(serial[2], parallel[2])


def _reload(dcwlf):
    """A save and load: the node configuration travels as the node-config
    property, the per-galaxy state does not."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    return ser.from_string(ser.to_string(dcwlf, True))


def test_node_config_replay_skips_calibration():
    """A reloaded dataset replays the stored node grid instead of searching
    again, and evaluates bitwise as the dataset that chose it."""
    dcwlf, mset = _build("var-add", "nfw", "composed")
    assert not dcwlf.is_data_prepared()
    dcwlf.data_prepare(mset)
    assert dcwlf.is_data_prepared()
    assert dcwlf.get_calib_count() > 0
    ref = _m2lnP(dcwlf, mset)

    # Saving again before any prepare keeps the loaded configuration.
    dup = _reload(_reload(dcwlf))
    assert dup.props.node_config is not None
    dup.data_prepare(mset)

    assert dup.get_calib_count() == 0
    assert dup.is_data_prepared()
    assert _node_config(dup) == _node_config(dcwlf)
    np.testing.assert_array_equal(_m2lnP(dup, mset), ref)


def test_node_config_mismatch_recalibrates():
    """A stored configuration made under other settings is discarded, and
    the recalibrated dataset is the one a fresh run with those settings
    gives."""
    dcwlf, mset = _build("var-add", "nfw", "composed")
    dcwlf.data_prepare(mset)

    dup = _reload(dcwlf)
    dup.set_node_reltol(1.0e-5)
    dup.data_prepare(mset)
    assert dup.get_calib_count() > 0

    fresh, mset2 = _build("var-add", "nfw", "composed")
    fresh.set_node_reltol(1.0e-5)
    fresh.data_prepare(mset2)

    assert _node_config(dup) == _node_config(fresh)
    np.testing.assert_array_equal(_m2lnP(dup, mset), _m2lnP(fresh, mset2))
    # A tighter tolerance never selects a coarser grid.
    assert sum(_node_config(dup)[0]) >= sum(_node_config(dcwlf)[0])


def test_node_config_cleared_and_prepared_state():
    """Clearing node-config makes the next prepare calibrate; the prepared
    flag follows the method: auto-nodes FIXED_NODES needs a configuration,
    the other methods only the per-galaxy data."""
    dcwlf, mset = _build("var-add", "nfw", "composed")
    dcwlf.data_prepare(mset)

    dup = _reload(dcwlf)
    dup.props.node_config = None
    assert dup.props.node_config is None
    dup.data_prepare(mset)
    assert dup.get_calib_count() > 0
    np.testing.assert_array_equal(_m2lnP(dup, mset), _m2lnP(dcwlf, mset))

    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    assert dcwlf.is_data_prepared()


def test_resample_skips_eager_precomputation():
    """The prepare() that ncm_data_resample() runs only feeds gen(): no
    MomentsTilt table and no auto-node search is spent on it, and the next
    prepare() rebuilds everything from the resampled catalog."""
    dcwlf, mset = _build("moments-tilt", "nfw", "composed")
    gsf = dcwlf.props.shape_factor
    rng = Ncm.RNG.seeded_new(None, 20260923)

    dcwlf.resample(mset, rng)
    assert gsf.get_table_build_count() == 0
    assert dcwlf.get_calib_count() == 0

    # The uncalibrated grid that prepare used is not reported as a node
    # configuration, so a save taken now does not replay it after a load.
    assert not dcwlf.is_data_prepared()
    assert dcwlf.props.node_config is None
    dup = _reload(dcwlf)
    dup.data_prepare(mset)
    assert dup.get_calib_count() > 0

    dcwlf.data_prepare(mset)
    assert gsf.get_table_build_count() > 0
    resampled = _m2lnP(dcwlf, mset)
    assert np.all(np.isfinite(resampled))

    # Same answer as a fresh dataset handed the resampled catalog.
    fresh, mset2 = _build("moments-tilt", "nfw", "composed")
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    fresh.set_obs(ser.dup_obj(dcwlf.peek_obs()))
    fresh.data_prepare(mset2)
    np.testing.assert_array_equal(_m2lnP(fresh, mset2), resampled)


def test_underflowing_shape_likelihood_is_not_the_wall():
    """|epsilon_obs| = 1.7 with a noise of 0.013 is ~1450 nats below the
    support: the linear shape likelihood is exactly 0 at every node, its log
    is finite. FIXED_NODES must return that log, not NC_GALAXY_LOW_PROB."""
    dcwlf, mset = _build("moments-tilt", "nfw", "composed", n_gal=5)
    dcwlf.set_auto_nodes(False)
    obs = dcwlf.peek_obs()
    obs.set("epsilon_obs_1", 4, -0.025)
    obs.set("epsilon_obs_2", 4, -1.697)
    obs.set("std_noise", 4, 0.013)
    dcwlf.set_obs(obs)

    fixed = _m2lnP(dcwlf, mset)
    assert dcwlf.get_low_prob_count() == 0
    assert np.all(np.isfinite(fixed))
    assert 2000.0 < fixed[4] < 1.0e5

    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)
    lnint = _m2lnP(dcwlf, mset)
    assert fixed[4] == pytest.approx(lnint[4], rel=1.0e-4)
