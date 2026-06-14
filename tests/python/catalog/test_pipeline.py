#!/usr/bin/env python
#
# test_pipeline.py
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

"""Tests for the end-to-end MockPipeline orchestrator."""

import math

import numpy as np

from numcosmo_py import Nc, Ncm
from numcosmo_py.catalog import MockCatalogs, MockPipeline

Ncm.cfg_init()

# A small footprint area keeps the realized counts (and hence the test) small.
AREA = 2.0 * (math.pi / 180.0) ** 2


def _setup():
    """Build cosmology, abundance, proxies and mset for the pipeline."""
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

    cad.set_area(AREA)
    cad.prepare(cosmo, cluster_z, cluster_m)

    return cosmo, cad, cluster_z, cluster_m, mset


def _pipeline():
    """A pipeline over the small-area setup with a deterministic central HOD."""
    cosmo, cad, cluster_z, cluster_m, mset = _setup()
    footprint = Ncm.SkyFootprintRectangular.new(10.0, 14.0, -2.0, 2.0)
    hod = Nc.GalaxyHODZheng07.new()
    hod.set_stochastic_central(False)
    return MockPipeline(cad, hod, footprint), mset


def test_generate_returns_linked_tables() -> None:
    """The pipeline returns linked cluster and member astropy tables."""
    pipeline, mset = _pipeline()
    rng = Ncm.RNG.seeded_new(None, 0)

    result = pipeline.generate(mset, rng)

    assert isinstance(result, MockCatalogs)
    assert "cluster_id" in result.clusters.colnames
    assert {"ra", "dec", "r_Delta"}.issubset(result.clusters.colnames)
    assert result.members.colnames == [
        "galaxy_id",
        "parent_id",
        "is_central",
        "ra",
        "dec",
        "z",
    ]
    assert len(result.clusters) > 0
    assert len(result.members) > 0


def test_member_parent_ids_reference_clusters() -> None:
    """Every member links back to a real cluster row."""
    pipeline, mset = _pipeline()
    rng = Ncm.RNG.seeded_new(None, 1)

    result = pipeline.generate(mset, rng)
    cluster_ids = set(np.asarray(result.clusters["cluster_id"]))
    parent_ids = set(np.asarray(result.members["parent_id"]))

    assert parent_ids.issubset(cluster_ids)
    # With a deterministic central each cluster yields at least its central.
    assert parent_ids == cluster_ids


def test_central_sits_on_its_cluster() -> None:
    """Each central member sits exactly on its parent cluster position."""
    pipeline, mset = _pipeline()
    rng = Ncm.RNG.seeded_new(None, 2)

    result = pipeline.generate(mset, rng)
    clusters = result.clusters
    members = result.members

    centrals = members[np.asarray(members["is_central"], dtype=bool)]
    # One central per cluster (deterministic central, every cluster occupied).
    assert len(centrals) == len(clusters)

    by_id = {int(cid): i for i, cid in enumerate(clusters["cluster_id"])}
    for row in centrals:
        i = by_id[int(row["parent_id"])]
        assert math.isclose(float(row["ra"]), float(clusters["ra"][i]), abs_tol=1e-12)
        assert math.isclose(float(row["dec"]), float(clusters["dec"][i]), abs_tol=1e-12)
        assert math.isclose(
            float(row["z"]), float(clusters["z_true"][i]), abs_tol=1e-12
        )
