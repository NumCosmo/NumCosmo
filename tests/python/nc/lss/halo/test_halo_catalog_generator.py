#!/usr/bin/env python
#
# test_halo_catalog_generator.py
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""Tests for NcHaloCatalogGenerator.

The generator owns the cluster-count sampling pipeline used by
NcDataClusterNCount. These tests check its catalog output shape and that, run
standalone with the same seed, it reproduces the NcDataClusterNCount truth
table (the two share the extracted pipeline).

The reference is the stored seed-0 NcDataClusterNCount catalog at
``data/truth_tables/cluster/ncount_resample_seed0.bin`` (see
test_ncount_resample_truth_table.py). The comparison is tolerance-based so it survives
cross-stack sub-ULP rounding while still catching real regressions.
"""

import math

import numpy as np

from numcosmo_py import Nc, Ncm
from numcosmo_py.catalog import catalog_to_table

Ncm.cfg_init()

AREA = 270 * (math.pi / 180.0) ** 2

TRUTH_TABLE_FILE = "truth_tables/cluster/ncount_resample_seed0.bin"
TRUTH_TABLE_RTOL = 1.0e-7
TRUTH_TABLE_ATOL = 1.0e-12


def _load_truth_table_columns() -> dict[str, np.ndarray]:
    """Load the reference seed-0 catalog as plain arrays keyed by column."""
    path = Ncm.cfg_get_data_filename(TRUTH_TABLE_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    truth = ser.from_binfile(path)
    assert isinstance(truth, Nc.DataClusterNCount)
    return {
        "lnM_obs": np.array(truth.get_lnM_obs().dup_array()),
        "z_obs": np.array(truth.get_z_obs().dup_array()),
        "lnM_true": np.array(truth.get_lnM_true().dup_array()),
        "z_true": np.array(truth.get_z_true().dup_array()),
    }


def _setup():
    """Build cosmology, abundance, models and mset matching the truth-table test."""
    cosmo = Nc.HICosmoDEXcdm(reion=Nc.HIReionCamb(), prim=Nc.HIPrimPowerLaw())

    dist = Nc.Distance.new(2.0)
    psml = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH())
    psml.require_kmin(1.0e-3)
    psml.require_kmax(1.0e3)
    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    # Explicit: the truth table was drawn with the filter at this tolerance.
    psf.set_reltol(1.0e-6)
    psf.set_best_lnr0()

    mulf = Nc.MultiplicityFuncBocquet.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(200.0)
    mulf.set_sim(Nc.MultiplicityFuncBocquetSim.DM)

    hmf = Nc.HaloMassFunction.new(dist, psf, mulf)
    hmf.prepare(cosmo)

    cluster_m = Nc.ClusterMassLnnormal(
        lnMobs_min=math.log(1.0e14), lnMobs_max=math.log(1.0e16)
    )
    cluster_z = Nc.ClusterPhotozGaussGlobal(
        pz_min=0.0, pz_max=0.7, z_bias=0.0, sigma0=0.03
    )
    cad = Nc.ClusterAbundance.new(hmf, None)

    mset = Ncm.MSet.new_array([cosmo, cluster_z, cluster_m])
    for name, value in (
        ("H0", 70.0),
        ("Omegab", 0.05),
        ("Omegac", 0.25),
        ("Omegax", 0.70),
        ("Tgamma0", 2.72),
        ("w", -1.0),
    ):
        cosmo.param_set_by_name(name, value)
    cluster_m.param_set_by_name("bias", 0.0)
    cluster_m.param_set_by_name("sigma", 0.2)

    return cosmo, cad, cluster_z, cluster_m, mset


def test_generate_shape_and_metadata() -> None:
    """The generated catalog is a cluster NcHaloCatalog with the expected columns."""
    cosmo, cad, cluster_z, cluster_m, mset = _setup()
    cad.set_area(AREA)
    cad.prepare(cosmo, cluster_z, cluster_m)

    gen = Nc.HaloCatalogGenerator.new(cad)
    assert gen.peek_abundance() is cad

    rng = Ncm.RNG.seeded_new(None, 0)
    hcat = gen.generate(mset, rng)

    assert isinstance(hcat, Nc.HaloCatalog)
    assert hcat.get_kind() == Nc.HaloCatalogKind.CLUSTER
    # LnNormal mass and global photo-z have observable length 1 and no params.
    assert list(hcat.peek_columns()) == ["z_true", "lnM_true", "z_obs_0", "lnM_obs_0"]
    assert hcat.len() == 4242


def test_generate_with_footprint_adds_positions() -> None:
    """Setting a footprint augments the catalog with in-bounds ra/dec columns."""
    cosmo, cad, cluster_z, cluster_m, mset = _setup()
    cad.set_area(AREA)
    cad.prepare(cosmo, cluster_z, cluster_m)

    footprint = Ncm.SkyFootprintRectangular.new(10.0, 40.0, -5.0, 25.0)
    gen = Nc.HaloCatalogGenerator.new(cad)
    gen.set_footprint(footprint)
    assert gen.peek_footprint() is footprint

    rng = Ncm.RNG.seeded_new(None, 0)
    hcat = gen.generate(mset, rng)
    table = catalog_to_table(hcat)

    assert table.colnames == ["z_true", "lnM_true", "z_obs_0", "lnM_obs_0", "ra", "dec"]
    ra = np.asarray(table["ra"], dtype=np.float64)
    dec = np.asarray(table["dec"], dtype=np.float64)
    assert ra.min() >= 10.0 and ra.max() <= 40.0
    assert dec.min() >= -5.0 and dec.max() <= 25.0
    # Positions are one row per object; interleaving their draws shifts the RNG
    # stream, so the accepted count may differ slightly from the no-footprint run
    # (4242) but stays near it.
    assert len(ra) == hcat.len()
    assert 4100 <= hcat.len() <= 4242


def test_generate_with_radius_adds_r_delta() -> None:
    """Enabling radius output appends an r_Delta column with the SO radius."""
    cosmo, cad, cluster_z, cluster_m, mset = _setup()
    cad.set_area(AREA)
    cad.prepare(cosmo, cluster_z, cluster_m)

    gen = Nc.HaloCatalogGenerator.new(cad)
    gen.set_with_radius(True)
    assert gen.get_with_radius()

    rng = Ncm.RNG.seeded_new(None, 0)
    table = catalog_to_table(gen.generate(mset, rng))

    assert table.colnames == [
        "z_true",
        "lnM_true",
        "z_obs_0",
        "lnM_obs_0",
        "r_Delta",
    ]

    # The mass function uses a critical Delta=200 definition, so
    # M = (4/3) pi Delta rho_crit(z) r^3 with rho_crit(z) = rho_crit0 h^2 E2(z).
    z_true = np.asarray(table["z_true"], dtype=np.float64)
    lnm_true = np.asarray(table["lnM_true"], dtype=np.float64)
    r_delta = np.asarray(table["r_Delta"], dtype=np.float64)

    rho_crit0 = Ncm.C.crit_mass_density_h2_solar_mass_Mpc3() * cosmo.h2()
    e2 = np.array([cosmo.E2(z) for z in z_true])
    delta_rho_bg = 200.0 * rho_crit0 * e2
    expected = np.cbrt(3.0 * np.exp(lnm_true) / (4.0 * np.pi * delta_rho_bg))

    assert r_delta.shape == z_true.shape
    np.testing.assert_allclose(r_delta, expected, rtol=1e-12)


def test_generate_matches_truth_table() -> None:
    """Standalone generation reproduces the NcDataClusterNCount truth table."""
    cosmo, cad, cluster_z, cluster_m, mset = _setup()
    cad.set_area(AREA)
    cad.prepare(cosmo, cluster_z, cluster_m)

    gen = Nc.HaloCatalogGenerator.new(cad)
    rng = Ncm.RNG.seeded_new(None, 0)
    table = catalog_to_table(gen.generate(mset, rng))

    got = {
        "z_true": np.asarray(table["z_true"], dtype=np.float64),
        "lnM_true": np.asarray(table["lnM_true"], dtype=np.float64),
        "z_obs": np.asarray(table["z_obs_0"], dtype=np.float64),
        "lnM_obs": np.asarray(table["lnM_obs_0"], dtype=np.float64),
    }
    truth = _load_truth_table_columns()

    for column, ref in truth.items():
        np.testing.assert_allclose(
            got[column],
            ref,
            rtol=TRUTH_TABLE_RTOL,
            atol=TRUTH_TABLE_ATOL,
            err_msg=column,
        )

    # The sampled (z, lnM) move at 1e-8 with the knots of the splines behind the
    # abundance; the observable draws on top of them are the pipeline's own.
    # lnM_obs - lnM_true is sigma times the draw, and with z_bias = 0 so is
    # (z_obs - z_true) / (1 + z_true): both pinned exactly.
    np.testing.assert_allclose(
        got["lnM_obs"] - got["lnM_true"],
        truth["lnM_obs"] - truth["lnM_true"],
        rtol=0.0,
        atol=1.0e-12,
        err_msg="lnM draw",
    )
    np.testing.assert_allclose(
        (got["z_obs"] - got["z_true"]) / (1.0 + got["z_true"]),
        (truth["z_obs"] - truth["z_true"]) / (1.0 + truth["z_true"]),
        rtol=0.0,
        atol=1.0e-12,
        err_msg="z draw",
    )
