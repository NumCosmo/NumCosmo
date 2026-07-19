#!/usr/bin/env python
#
# test_cluster_mass_rich_count.py
#
# Fri Jul 10 00:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_cluster_mass_rich_count.py
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

"""Tests on Nc.DataClusterMassRichCount class."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import duplicate_via_serialization

Ncm.cfg_init()


def test_constructor():
    """Test constructor."""
    dmrc = Nc.DataClusterMassRichCount.new()
    assert dmrc is not None
    assert isinstance(dmrc, Nc.DataClusterMassRichCount)

    dmrc2 = dmrc.ref()
    assert dmrc2 == dmrc


@pytest.fixture(name="ascaso")
def fixture_ascaso() -> Nc.ClusterMassAscaso:
    """Fixture for ClusterMassAscaso."""
    ascaso = Nc.ClusterMassAscaso()

    return ascaso


@pytest.fixture(name="n_clusters")
def fixture_n_clusters() -> int:
    """Fixture for number of clusters."""
    return 200


@pytest.fixture(name="cluster_mass_rich_count")
def fixture_cluster_mass_rich_count(
    ascaso: Nc.ClusterMassAscaso, n_clusters: int
) -> Nc.DataClusterMassRichCount:
    """Fixture for DataClusterMassRichCount."""
    dmrc = Nc.DataClusterMassRichCount.new()
    rng = Ncm.RNG.new()

    # Generate random uniform log cluster masses from 10^13M_solar to 10^15M_solar
    log_cluster_masses = np.random.uniform(
        13.0 * np.log(10.0), 15.0 * np.log(10.0), n_clusters
    )

    # Generate random uniform redshift from 0.1 to 0.9
    redshifts = np.random.uniform(0.1, 0.9, n_clusters)

    log_cluster_masses_v = Ncm.Vector.new_array(log_cluster_masses)
    redshifts_v = Ncm.Vector.new_array(redshifts)
    # N is set through resample below, initialize with zeros.
    N_v = Ncm.Vector.new(n_clusters)
    N_v.set_zero()

    dmrc.set_data(log_cluster_masses_v, redshifts_v, N_v)
    mset = Ncm.MSet.new_array([ascaso])

    dmrc.resample(mset, rng)
    N = np.array(dmrc.peek_N().dup_array())

    assert len(N) <= n_clusters

    return dmrc


def test_resample_counts_are_nonnegative_integers(
    cluster_mass_rich_count: Nc.DataClusterMassRichCount,
):
    """Test that resampled richness counts are non-negative integers >= cut."""
    N = np.array(cluster_mass_rich_count.peek_N().dup_array())

    assert np.all(N >= 0.0)
    assert_allclose(N, np.round(N))
    # Default cut is 0.0 -> N_cut = round(exp(0.0)) = 1
    assert np.all(N >= 1.0)


@pytest.fixture(name="fit")
def fixture_fit(
    cluster_mass_rich_count: Nc.DataClusterMassRichCount, ascaso: Nc.ClusterMassAscaso
) -> Ncm.Fit:
    """Fixture for NcmFit object."""
    dset = Ncm.Dataset.new_array([cluster_mass_rich_count])
    likelihood = Ncm.Likelihood.new(dset)
    mset = Ncm.MSet.new_array([ascaso])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    ascaso.param_set_desc("cut", {"fit": False})
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT,
        "ln-neldermead",
        likelihood,
        mset,
        Ncm.FitGradType.NUMDIFF_CENTRAL,
    )
    return fit


def test_data_cluster_mass_rich_count_fit(fit: Ncm.Fit):
    """Test DataClusterMassRichCount fit."""
    mset = fit.peek_mset()
    fparam_len = mset.fparam_len()
    original_params = np.array([mset.fparam_get(i) for i in range(fparam_len)])

    fit.run_restart(Ncm.FitRunMsgs.NONE, 1.0e-2, 0.0)
    new_params = np.array([mset.fparam_get(i) for i in range(fparam_len)])

    # Check that original_params - new_params is close given cov
    diff = original_params - new_params
    assert np.sum(diff**2) < fparam_len


@pytest.mark.omp  # FitMC.set_use_threads(True) exercises the OpenMP-parallel MC path
@pytest.mark.parametrize(
    "mc_type",
    [Ncm.FitMCResampleType.FROM_MODEL, Ncm.FitMCResampleType.BOOTSTRAP_NOMIX],
    ids=["FROM_MODEL", "BOOTSTRAP_NOMIX"],
)
def test_data_cluster_mass_rich_count_bootstrap(
    fit: Ncm.Fit, mc_type: Ncm.FitMCResampleType
):
    """Test DataClusterMassRichCount bootstrap."""
    mset = fit.peek_mset()
    fparam_len = mset.fparam_len()
    original_params = np.array([mset.fparam_get(i) for i in range(fparam_len)])

    mc = Ncm.FitMC.new(fit, mc_type, Ncm.FitRunMsgs.NONE)
    mc.set_use_threads(True)
    mc.start_run()
    mc.run(100)
    mc.end_run()
    mcat = mc.get_catalog()
    assert mcat is not None
    assert isinstance(mcat, Ncm.MSetCatalog)

    cov_m = mcat.get_covar()
    mean_v = mcat.get_mean()
    cov = np.array(cov_m.dup_array()).reshape(fparam_len, fparam_len)

    new_params = np.array(mean_v.dup_array())

    # Check that original_params - new_params is close given cov
    diff = original_params - new_params
    chi2 = np.dot(diff, np.linalg.solve(cov, diff))
    # Bootstrap estimation can be biased
    if mc_type == Ncm.FitMCResampleType.FROM_MODEL:
        assert (
            chi2 < fparam_len * 9.0
        ), "Parameters differ too much from original values"


def test_data_cluster_mass_rich_count_apply_cut(
    cluster_mass_rich_count: Nc.DataClusterMassRichCount, n_clusters: int
):
    """Test DataClusterMassRichCount apply cut."""
    for N_min in range(0, 20, 2):
        cluster_mass_rich_count.apply_cut(N_min)
        N = np.array(cluster_mass_rich_count.peek_N().dup_array())
        assert len(N) <= n_clusters
        assert all(N >= N_min)


def test_serialize_deserialize(cluster_mass_rich_count: Nc.DataClusterMassRichCount):
    """Test serialize and deserialize."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dmrc2 = duplicate_via_serialization(cluster_mass_rich_count, ser)
    assert isinstance(dmrc2, Nc.DataClusterMassRichCount)

    assert dmrc2 is not cluster_mass_rich_count

    assert dmrc2.get_length() == cluster_mass_rich_count.get_length()
    assert dmrc2.get_dof() == cluster_mass_rich_count.get_dof()

    assert_allclose(
        cluster_mass_rich_count.peek_lnM().dup_array(),
        dmrc2.peek_lnM().dup_array(),
    )
    assert_allclose(
        cluster_mass_rich_count.peek_N().dup_array(),
        dmrc2.peek_N().dup_array(),
    )
    assert_allclose(
        cluster_mass_rich_count.peek_z().dup_array(),
        dmrc2.peek_z().dup_array(),
    )
