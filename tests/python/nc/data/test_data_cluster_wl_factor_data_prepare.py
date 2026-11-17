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

Each run builds everything from scratch (models, catalog, dataset), so the
models' lazy first-call state (halo-position rotation, LSST-SRD constants, the
numerically integrated Einasto profile's splines) starts
unset and the orchestrator's serial warm-up is exercised every time. The
catalog mixes background galaxies with foreground-only ones, so the warm-up
has to skip galaxies without a redshift grid.
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
    "spectral": lambda: _spectral(),
    "moments-tilt": lambda: Nc.GalaxyShapeFactorMomentsTilt.new(ELLIP_CONV),
    "moments-gauss": lambda: Nc.GalaxyShapeFactorMomentsGauss.new(ELLIP_CONV),
}


def _spectral():
    """Spectral computes a 2D integral per evaluation and the calibration needs
    thousands of them per galaxy: a loose tolerance and a short catalog (see
    N_GAL_OF) keep its run time that of a unit test. The test compares
    parallel with serial, not accuracy."""
    gsf = Nc.GalaxyShapeFactorSpectral.new(ELLIP_CONV)
    gsf.set_reltol(1.0e-3)
    return gsf


#: Galaxies per shape factor; the rest use N_GAL. Still more than the three
#: foreground-only galaxies plus one warm-up galaxy plus PARALLEL_THREADS.
N_GAL_OF = {"spectral": 8}

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
    zp[:3] = 0.05  # foreground-only first: the warm-up must skip them
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


def _build(
    shape_name, profile_name, redshift_name, n_gal=None, zf=100.0, separate_dist=False
):
    """Builds a fresh FIXED_NODES + auto-nodes dataset; returns (dcwlf, mset). The
    surface mass density has a distance of final redshift zf; with separate_dist the
    halo position has its own distance (zf = 100), otherwise it shares that one."""
    if n_gal is None:
        n_gal = N_GAL_OF.get(shape_name, N_GAL)
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(zf)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hms.param_set_by_name("log10MDelta", 14.5)
    dp = PROFILES[profile_name](hms)
    hp = Nc.HaloPosition.new(Nc.Distance.new(100.0) if separate_dist else dist)
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


def _z_top(dcwlf):
    obs = dcwlf.peek_obs()
    return max(obs.peek_pz(i).get_bounds()[1] for i in range(obs.len()))


ZF_WARNING = "above the final redshift zf"


@pytest.mark.parametrize("separate_dist", [False, True])
def test_distance_zf_get_z_max(capfd, separate_dist):
    """The comoving distance is a spline only up to the zf of NcDistance, and above it an
    integral cached at every redshift. The likelihood does not change the distance, a
    model other data may share: it warns when the redshift integration of the galaxies
    goes above zf, and nc_data_cluster_wl_factor_get_z_max() gives the zf to require.
    Raising it after an evaluation, as the warning says and with nothing prepared by
    hand, the models rebuild the distance: the likelihood matches that of a distance
    built with a large zf, also after the cluster moves."""
    zf = 0.5 if separate_dist else 0.55
    ref, mset_ref = _build("var-add", "nfw", "spline", separate_dist=separate_dist)
    dcwlf, mset = _build("var-add", "nfw", "spline", zf=zf, separate_dist=separate_dist)
    dist = mset.peek(Nc.WLSurfaceMassDensity.id()).props.distance

    z_max = dcwlf.get_z_max(mset)
    assert z_max == _z_top(dcwlf) > zf

    capfd.readouterr()
    _m2lnP(ref, mset_ref)
    _m2lnP(dcwlf, mset)
    assert ZF_WARNING in capfd.readouterr().err
    assert dist.props.zf == zf

    dist.require_zf(z_max)
    np.testing.assert_allclose(
        _m2lnP(dcwlf, mset), _m2lnP(ref, mset_ref), rtol=1.0e-8, atol=0.0
    )

    for m in (mset, mset_ref):
        m.peek(Nc.HaloPosition.id()).param_set_by_name("z", 0.3)

    m2lnP = _m2lnP(dcwlf, mset)

    assert ZF_WARNING not in capfd.readouterr().err
    np.testing.assert_allclose(m2lnP, _m2lnP(ref, mset_ref), rtol=1.0e-8, atol=0.0)


def test_spectral_calibration_through_the_line_table():
    """Spectral's calibration evaluates P(eps_obs | z) from the table along the shears of
    each galaxy (line-abstol, default 1e-8) instead of an integral at every probe: the node
    pairs are those of the integrals (line-abstol zero), and so is -2lnP within the
    tolerance of the table."""
    results = {}
    for line_abstol in (0.0, 1.0e-8):
        dcwlf, mset = _build("spectral", "nfw", "composed")
        dcwlf.props.shape_factor.props.line_abstol = line_abstol
        dcwlf.data_prepare(mset)
        assert dcwlf.get_calib_count() > 0
        results[line_abstol] = (*_node_config(dcwlf), _m2lnP(dcwlf, mset))

    direct, table = results[0.0], results[1.0e-8]

    assert table[0] == direct[0]
    assert table[1] == direct[1]
    assert np.all(np.isfinite(table[2]))
    np.testing.assert_allclose(table[2], direct[2], rtol=1.0e-7, atol=1.0e-7)


def _shared_clone(dcwlf, mset):
    """A clone made as NcmFitESMCMC makes its workers': the dataset registered as shared
    after it was prepared, then duplicated (with its own model set) through that
    serializer."""
    dset = Ncm.Dataset.new()
    dset.append_data(dcwlf)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    dset.register_shared(ser)
    clone = ser.dup_obj(dcwlf)
    mset_clone = ser.dup_obj(mset)
    ser.reset(True)
    assert clone.peek_obs() is dcwlf.peek_obs()
    return clone, mset_clone


def test_clone_takes_the_published_grid_and_critical_densities():
    """A clone of a prepared dataset takes its redshift grid and the critical densities
    at its nodes instead of calibrating and computing them, with the same -2lnP. A
    clone at another cosmology computes its own critical densities, which equal those
    of a dataset prepared there, and leaves the original's untouched."""
    dcwlf, mset = _build("var-add", "nfw", "spline")
    dcwlf.data_prepare(mset)
    ref = _m2lnP(dcwlf, mset)

    n_gal = dcwlf.peek_obs().len()
    assert dcwlf.get_crit_count() == n_gal

    clone, mset_clone = _shared_clone(dcwlf, mset)
    np.testing.assert_array_equal(_m2lnP(clone, mset_clone), ref)
    assert clone.get_calib_count() == 0
    # The critical densities were taken from the original, not computed
    assert clone.get_crit_count() == 0

    # Another cosmology in the clone only
    mset_clone.peek(Nc.HICosmo.id()).param_set_by_name("H0", 72.0)
    moved = _m2lnP(clone, mset_clone)
    assert clone.get_crit_count() == n_gal
    np.testing.assert_array_equal(_m2lnP(dcwlf, mset), ref)
    assert not np.array_equal(moved, ref)

    fresh, mset_fresh = _build("var-add", "nfw", "spline")
    fresh.props.node_config = dcwlf.props.node_config
    mset_fresh.peek(Nc.HICosmo.id()).param_set_by_name("H0", 72.0)
    np.testing.assert_array_equal(_m2lnP(fresh, mset_fresh), moved)


def _shared_clones(dcwlf, mset, n):
    """n clones made as NcmFitESMCMC makes its workers', each with its own model set,
    after the dataset was registered as shared (prepared or not)."""
    dset = Ncm.Dataset.new()
    dset.append_data(dcwlf)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    dset.register_shared(ser)
    clones = []
    for _ in range(n):
        clones.append((ser.dup_obj(dcwlf), ser.dup_obj(mset)))
        ser.reset(True)
    return clones


def test_clones_of_an_unprepared_dataset_build_the_grid_once():
    """The clones of a dataset that was not prepared before they were made build the
    redshift grid once: the first calibrates and computes the critical densities, and
    publishes them; the second takes the grid, its pairs of nodes and the critical
    densities, with the same -2lnP bit for bit, as does the original prepared after."""
    dcwlf, mset = _build("var-add", "nfw", "spline")
    n_gal = dcwlf.peek_obs().len()
    (first, mset_first), (second, mset_second) = _shared_clones(dcwlf, mset, 2)

    ref = _m2lnP(first, mset_first)
    assert first.get_calib_count() == n_gal
    assert first.get_crit_count() == n_gal

    np.testing.assert_array_equal(_m2lnP(second, mset_second), ref)
    assert second.get_calib_count() == 0
    assert second.get_crit_count() == 0
    assert _node_config(second) == _node_config(first)

    np.testing.assert_array_equal(_m2lnP(dcwlf, mset), ref)
    assert dcwlf.get_calib_count() == 0


def test_distance_zf_after_model_swap(capfd):
    """A surface mass density swapped in after the first evaluation, with a distance of
    low zf, is warned about too, although the catalog did not change. The warning is
    given once in the process for a redshift and a zf (the zf here is used by no other
    test), so the copies a sampler makes of the data do not repeat it."""
    dcwlf, mset = _build("var-add", "nfw", "spline", separate_dist=True)
    cosmo = mset.peek(Nc.HICosmo.id())
    _m2lnP(dcwlf, mset)

    smd = Nc.WLSurfaceMassDensity.new(Nc.Distance.new(0.6))
    smd.prepare_if_needed(cosmo)
    mset.set(smd)

    capfd.readouterr()
    _m2lnP(dcwlf, mset)
    assert ZF_WARNING in capfd.readouterr().err
    assert smd.props.distance.props.zf == 0.6

    copy = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).dup_obj(dcwlf)
    _m2lnP(copy, mset)
    assert ZF_WARNING not in capfd.readouterr().err


def test_peek_pz_is_usable():
    """The likelihood reads only the knots and values of the P(z), so a catalog keeps
    them as tables; nc_galaxy_wl_obs_peek_pz() still returns a spline that evaluates,
    here from a catalog loaded from its serialization, before and after the likelihood
    read it."""
    dcwlf, mset = _build("var-add", "nfw", "spline")
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    obs = ser.from_string(ser.to_string(dcwlf.peek_obs(), True))
    dcwlf.set_obs(obs)
    _, _, zp, sigma0, *_ = _galaxies(obs.len())

    def peak(i):
        return 1.0 / (np.sqrt(2.0 * np.pi) * sigma0[i] * (1.0 + zp[i]))

    np.testing.assert_allclose(obs.peek_pz(5).eval(zp[5]), peak(5), rtol=1.0e-3)
    _m2lnP(dcwlf, mset)
    np.testing.assert_allclose(obs.peek_pz(6).eval(zp[6]), peak(6), rtol=1.0e-3)


def _reload(dcwlf):
    """A save and load: the node configuration travels as the node-config
    property, the per-galaxy state does not."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    return ser.from_string(ser.to_string(dcwlf, True))


def test_node_config_replay_skips_calibration():
    """A reloaded dataset replays the stored node grid instead of calibrating
    again, and evaluates bitwise identically to the dataset that calibrated
    it."""
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
    MomentsTilt table is built and no auto-node calibration is run for it,
    and the next prepare() rebuilds everything from the resampled catalog."""
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
