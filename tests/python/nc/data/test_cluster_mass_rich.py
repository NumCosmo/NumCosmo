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
from numcosmo_py.helper import duplicate_via_serialization

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
    # Create sigma_lnR vector (zeros for no observational error in this test)
    sigma_lnR_v = Ncm.Vector.new(n_clusters)
    sigma_lnR_v.set_zero()

    dcmr.set_data(log_cluster_masses_v, redshifts_v, log_R_v, sigma_lnR_v)
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


@pytest.mark.omp  # FitMC.set_use_threads(True) exercises the OpenMP-parallel MC path
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
    cluster_mass_rich2 = duplicate_via_serialization(cluster_mass_rich, ser)
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


def _projection(fprj: float, tau: float = 0.1) -> Nc.ClusterMassProjection:
    """A projection model sharing the log-normal part with the default Ascaso."""
    ascaso = Nc.ClusterMassAscaso()
    projection = Nc.ClusterMassProjection()

    for i in range(ascaso.sparam_len()):
        projection[ascaso.param_name(i)] = ascaso.param_get(i)

    projection["fprj"] = fprj
    projection["tau"] = tau

    return projection


def _dataset(model, lnM, z, lnR):
    """A single-cluster-set dataset built on a given richness model."""
    dcmr = Nc.DataClusterMassRich.new()
    sigma_lnR = Ncm.Vector.new(len(lnM))
    sigma_lnR.set_zero()
    dcmr.set_data(
        Ncm.Vector.new_array(lnM),
        Ncm.Vector.new_array(z),
        Ncm.Vector.new_array(lnR),
        sigma_lnR,
    )

    return Ncm.Dataset.new_array([dcmr]), Ncm.MSet.new_array([model])


def test_is_lognormal():
    """Only a model that adds a component to the richness reports FALSE."""
    assert Nc.ClusterMassAscaso().is_lognormal()
    assert _projection(0.0).is_lognormal()
    assert not _projection(0.2).is_lognormal()


@pytest.fixture(name="sample")
def fixture_sample(ascaso: Nc.ClusterMassAscaso):
    """Masses, redshifts and richnesses above the cut, from the clean relation."""
    rng = np.random.default_rng(7)
    n = 300
    lnM = rng.uniform(14.0 * np.log(10.0), 15.0 * np.log(10.0), n)
    z = rng.uniform(0.2, 0.9, n)
    mu = np.array([ascaso.mu(m, zz) for m, zz in zip(lnM, z)])
    sigma = np.array([ascaso.sigma(m, zz) for m, zz in zip(lnM, z)])
    cut = ascaso.get_cut()
    lnR = np.maximum(rng.normal(mu, sigma), cut + 0.01)

    return lnM.tolist(), z.tolist(), lnR.tolist()


def test_m2lnL_paths_agree(ascaso: Nc.ClusterMassAscaso, sample):
    """With f_prj -> 0 the generic P/intP_bin path must match the closed form."""
    lnM, z, lnR = sample

    dset_a, mset_a = _dataset(ascaso, lnM, z, lnR)
    dset_p, mset_p = _dataset(_projection(1.0e-12), lnM, z, lnR)

    assert_allclose(dset_p.m2lnL_val(mset_p), dset_a.m2lnL_val(mset_a), rtol=1.0e-9)


def test_m2lnL_responds_to_projection(sample):
    """f_prj and tau must reach the likelihood, not just sit in the model."""
    lnM, z, lnR = sample

    def m2lnL(fprj, tau):
        dset, mset = _dataset(_projection(fprj, tau), lnM, z, lnR)

        return dset.m2lnL_val(mset)

    reference = m2lnL(0.0, 0.1)

    assert m2lnL(0.2, 0.1) != pytest.approx(reference)
    assert m2lnL(0.2, 0.5) != pytest.approx(m2lnL(0.2, 0.1))


def test_m2lnL_outlier_stays_finite(ascaso: Nc.ClusterMassAscaso, sample):
    """A cluster deep in the tail costs a large but finite penalty."""
    lnM, z, lnR = sample
    far = list(lnR)
    far[0] = ascaso.mu(lnM[0], z[0]) + 40.0 * ascaso.sigma(lnM[0], z[0])

    dset, mset = _dataset(ascaso, lnM, z, far)

    assert np.isfinite(dset.m2lnL_val(mset))
