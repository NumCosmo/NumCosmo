#!/usr/bin/env python
#
# test_py_data_cluster_mass_rich.py
#
# Tue Nov 11 10:15:31 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_cluster_mass_rich.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcmDataClusterMassRich class."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_constructor():
    """Test constructor."""
    cluster_mass_rich = Nc.DataClusterMassRich.new()
    assert cluster_mass_rich is not None
    assert isinstance(cluster_mass_rich, Nc.DataClusterMassRich)

    cluster_mass_rich2 = cluster_mass_rich.ref()
    assert cluster_mass_rich2 == cluster_mass_rich


@pytest.fixture(name="ascaso")
def fixture_ascaso() -> Nc.ClusterMassAscaso:
    """Fixture for ClusterMassAscaso."""
    ascaso = Nc.ClusterMassAscaso()

    return ascaso


@pytest.fixture(name="n_clusters")
def fixture_n_clusters() -> int:
    """Fixture for number of clusters."""
    return 100


@pytest.fixture(name="cluster_mass_rich")
def fixture_cluster_mass_rich(
    ascaso: Nc.ClusterMassAscaso, n_clusters: int
) -> Nc.DataClusterMassRich:
    """Fixture for DataClusterMassRich."""
    dcmr = Nc.DataClusterMassRich.new()
    rng = Ncm.RNG.new()

    # Generate random uniform log cluster masses from 10^13M_solar to 10^15M_solar
    log_cluster_masses = np.random.uniform(
        13.0 * np.log(10.0), 15.0 * np.log(10.0), n_clusters
    )

    # Generate random uniform redshift from 0.1 to 0.9
    redshifts = np.random.uniform(0.1, 0.9, n_clusters)

    log_cluster_masses_v = Ncm.Vector.new_array(log_cluster_masses)
    redshifts_v = Ncm.Vector.new_array(redshifts)
    log_R_v = log_cluster_masses_v.dup()

    dcmr.set_data(log_cluster_masses_v, redshifts_v, log_R_v)
    mset = Ncm.MSet.new_array([ascaso])

    dcmr.resample(mset, rng)
    lnR = np.array(dcmr.peek_lnR().dup_array())

    assert len(lnR) <= n_clusters

    return dcmr


@pytest.fixture(name="fit")
def fixture_fit(
    cluster_mass_rich: Nc.DataClusterMassRich, ascaso: Nc.ClusterMassAscaso
) -> Ncm.Fit:
    """Fixture for NcmFit object."""
    dset = Ncm.Dataset.new_array([cluster_mass_rich])
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


def test_data_cluster_mass_rich_fit(fit: Ncm.Fit):
    """Test DataClusterMassRich fit."""
    mset = fit.peek_mset()
    fparam_len = mset.fparam_len()
    original_params = np.array([mset.fparam_get(i) for i in range(fparam_len)])

    fit.run_restart(Ncm.FitRunMsgs.NONE, 1.0e-2, 0.0)
    new_params = np.array([mset.fparam_get(i) for i in range(fparam_len)])

    # Check that original_params - new_params is close given cov
    diff = original_params - new_params
    assert np.sum(diff**2) < fparam_len


@pytest.mark.parametrize(
    "mc_type",
    [Ncm.FitMCResampleType.FROM_MODEL, Ncm.FitMCResampleType.BOOTSTRAP_NOMIX],
    ids=["FROM_MODEL", "BOOTSTRAP_NOMIX"],
)
def test_data_cluster_mass_rich_bootstrap(fit: Ncm.Fit, mc_type: Ncm.FitMCResampleType):
    """Test DataClusterMassRich bootstrap."""
    mset = fit.peek_mset()
    fparam_len = mset.fparam_len()
    original_params = np.array([mset.fparam_get(i) for i in range(fparam_len)])

    mc = Ncm.FitMC.new(fit, mc_type, Ncm.FitRunMsgs.NONE)
    mc.set_nthreads(2)
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


def test_data_cluster_mass_rich_apply_cut(
    cluster_mass_rich: Nc.DataClusterMassRich, n_clusters: int
):
    """Test DataClusterMassRich apply cut."""
    for cut in np.linspace(0.0, 1.0, 10):
        cluster_mass_rich.apply_cut(cut)
        lnR = np.array(cluster_mass_rich.peek_lnR().dup_array())
        assert len(lnR) <= n_clusters
        assert all(lnR > cut)


def test_serialize_deserialize(cluster_mass_rich: Nc.DataClusterMassRich):
    """Test serialize and deserialize."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    cluster_mass_rich2 = ser.dup_obj(cluster_mass_rich)
    assert isinstance(cluster_mass_rich2, Nc.DataClusterMassRich)

    assert cluster_mass_rich2 is not cluster_mass_rich

    assert cluster_mass_rich2.get_length() == cluster_mass_rich.get_length()
    assert cluster_mass_rich2.get_dof() == cluster_mass_rich.get_dof()

    assert_allclose(
        cluster_mass_rich.peek_lnM().dup_array(),
        cluster_mass_rich2.peek_lnM().dup_array(),
    )
    assert_allclose(
        cluster_mass_rich.peek_lnR().dup_array(),
        cluster_mass_rich2.peek_lnR().dup_array(),
    )
    assert_allclose(
        cluster_mass_rich.peek_z().dup_array(),
        cluster_mass_rich2.peek_z().dup_array(),
    )
