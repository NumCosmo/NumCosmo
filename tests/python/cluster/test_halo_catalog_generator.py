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
standalone with the same seed, it reproduces the pinned NcDataClusterNCount
golden snapshot bit-for-bit (the two share the extracted pipeline).
"""

import math
import hashlib

import numpy as np

from numcosmo_py import Nc, Ncm
from numcosmo_py.catalog import catalog_to_table

Ncm.cfg_init()

AREA = 270 * (math.pi / 180.0) ** 2


def _digest(array: np.ndarray) -> str:
    """Stable short digest of an array's raw bytes."""
    return hashlib.sha256(array.tobytes()).hexdigest()[:16]


def _setup():
    """Build cosmology, abundance, models and mset matching the golden test."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.add_submodel(Nc.HIReionCamb())
    cosmo.add_submodel(Nc.HIPrimPowerLaw())

    dist = Nc.Distance.new(2.0)
    psml = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH())
    psml.require_kmin(1.0e-3)
    psml.require_kmax(1.0e3)
    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
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


def test_generate_matches_golden_snapshot() -> None:
    """Standalone generation reproduces the NcDataClusterNCount golden snapshot."""
    cosmo, cad, cluster_z, cluster_m, mset = _setup()
    cad.set_area(AREA)
    cad.prepare(cosmo, cluster_z, cluster_m)

    gen = Nc.HaloCatalogGenerator.new(cad)
    rng = Ncm.RNG.seeded_new(None, 0)
    table = catalog_to_table(gen.generate(mset, rng))

    z_true = np.asarray(table["z_true"], dtype=np.float64)
    lnM_true = np.asarray(table["lnM_true"], dtype=np.float64)
    z_obs = np.asarray(table["z_obs_0"], dtype=np.float64)
    lnM_obs = np.asarray(table["lnM_obs_0"], dtype=np.float64)

    assert _digest(lnM_obs) == "daffcd2e0548bb36"
    assert _digest(z_obs) == "42aa22695f8e8d7f"
    assert _digest(lnM_true) == "aa0e1d840a11608d"
    assert _digest(z_true) == "8d16a125aa223efa"
