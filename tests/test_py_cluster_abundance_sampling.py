#!/usr/bin/env python
#
# test_py_cluster_abundance_sampling.py
#
# Fri May 19 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_cluster_abundance_sampling.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Testing the cluster abundance data object to generate a mock catalog.

The mock catalog is generated by sampling the cosmological model and the
mass-observable relation. The sampling is done using the NcmData object
which is initialized with a ClusterAbundance object and a survey area.
The NcmData object is then initialized using the init_from_sampling
function which samples the cosmological model and the mass-observable
relation and stores the results in the NcmData object. The NcmData object
can then be used to generate a mock catalog by calling the catalog_save
function.
"""

import math
import pytest
import numpy as np
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()


@pytest.fixture(name="cosmo_hmf")
def fixture_cosmo_hmf():
    """Fixture for the halo mass function."""
    cosmo = Nc.HICosmoDEXcdm()
    reion = Nc.HIReionCamb()
    prim = Nc.HIPrimPowerLaw()

    cosmo.add_submodel(reion)
    cosmo.add_submodel(prim)

    dist = Nc.Distance.new(2.0)

    tf = Nc.TransferFuncEH()

    psml = Nc.PowspecMLTransfer.new(tf)
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

    return cosmo, hmf


def test_cluster_abundance_sampling(cosmo_hmf):
    """Test the cluster abundance sampling.

    Example of using the cluster abundance module to generate a mock catalog.
    """
    cosmo, hmf = cosmo_hmf

    lnMobs_min = math.log(1.0e14)
    lnMobs_max = math.log(1.0e16)
    cluster_m = Nc.ClusterMassLnnormal(lnMobs_min=lnMobs_min, lnMobs_max=lnMobs_max)

    z_min = 0.0
    z_max = 0.7
    cluster_z = Nc.ClusterPhotozGaussGlobal(
        pz_min=z_min, pz_max=z_max, z_bias=0.0, sigma0=0.03
    )

    cad = Nc.ClusterAbundance.new(hmf, None)

    ncdata = Nc.DataClusterNCount.new(
        cad, "NcClusterPhotozGaussGlobal", "NcClusterMassLnnormal"
    )

    mset = Ncm.MSet.new_array([cosmo, cluster_z, cluster_m])

    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_by_name("Omegab", 0.05)
    cosmo.param_set_by_name("Omegac", 0.25)
    cosmo.param_set_by_name("Omegax", 0.70)
    cosmo.param_set_by_name("Tgamma0", 2.72)
    cosmo.param_set_by_name("w", -1.0)

    cluster_m.param_set_by_name("bias", 0.0)
    cluster_m.param_set_by_name("sigma", 0.2)

    rng = Ncm.RNG.seeded_new(None, 0)

    ncdata.init_from_sampling(mset, 270 * (math.pi / 180.0) ** 2, rng)

    assert ncdata.has_lnM_true()
    lnM_true = ncdata.get_lnM_true()

    assert ncdata.has_z_true()
    z_true = ncdata.get_z_true()

    lnM_obs = ncdata.get_lnM_obs()
    lnM_obs_params = ncdata.get_lnM_obs_params()

    z_obs = ncdata.get_z_obs()
    z_obs_params = ncdata.get_z_obs_params()

    assert ncdata.get_len() == 4242

    assert ncdata.get_len() == lnM_obs.col_len()
    assert ncdata.get_len() == z_obs.col_len()
    assert ncdata.get_len() == lnM_true.len()
    assert ncdata.get_len() == z_true.len()
    assert lnM_obs.row_len() == 1
    assert z_obs.row_len() == 1
    assert lnM_obs_params is None
    assert z_obs_params is None

    assert np.isfinite(lnM_true.dup_array()).all()
    assert np.isfinite(z_true.dup_array()).all()
    assert np.isfinite(lnM_obs.dup_array()).all()
    assert np.isfinite(z_obs.dup_array()).all()


def test_hmf_volume(cosmo_hmf):
    """Test the halo mass function volume.

    Test the comoving volume of the halo mass function.
    """
    cosmo, hmf = cosmo_hmf

    dist = Nc.Distance.new(2.0)
    dist.prepare(cosmo)

    z = 0.3
    dV = hmf.dv_dzdomega(cosmo, z)

    assert dV > 0.0
    assert np.isfinite(dV)

    dist_dV = dist.comoving_volume_element(cosmo, z)

    assert np.isclose(dV, dist_dV, rtol=1.0e-13, atol=0.0)
