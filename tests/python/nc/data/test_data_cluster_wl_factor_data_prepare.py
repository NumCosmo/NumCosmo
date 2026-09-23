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


def _prepare(shape_name, profile_name, redshift_name, n_threads):
    """Builds a fresh dataset, runs data_prepare() on @n_threads threads and
    returns (node-config n-nodes, node-config rule-n, per-galaxy -2lnP)."""
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

    Ncm.cfg_set_openmp_nthreads(n_threads)
    try:
        dcwlf.data_prepare(mset)
    finally:
        Ncm.cfg_set_openmp_nthreads(
            int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))
        )

    node_config = dcwlf.props.node_config
    n_nodes = node_config.get_variant("n-nodes")[1].unpack()
    rule_n = node_config.get_variant("rule-n")[1].unpack()

    m2lnP = Ncm.Vector.new(n_gal)
    dcwlf.eval_m2lnP_gal(mset, m2lnP)

    return n_nodes, rule_n, np.array(m2lnP.dup_array())


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
