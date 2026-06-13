#!/usr/bin/env python
#
# test_ncount_resample_golden.py
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

"""Golden snapshot guarding NcDataClusterNCount.resample bit-for-bit.

The resample draw order is RNG-sensitive. This test pins the exact output for a
fixed seed so the planned extraction of the sampling pipeline into
NcHaloCatalogGenerator can be verified to preserve behavior bit-for-bit.
"""

import math
import hashlib

import numpy as np

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()


def _digest(array: np.ndarray) -> str:
    """Stable short digest of an array's raw bytes."""
    return hashlib.sha256(array.tobytes()).hexdigest()[:16]


def _resampled_ncount() -> Nc.DataClusterNCount:
    """Build and resample an NcDataClusterNCount deterministically (seed 0)."""
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
    ncdata = Nc.DataClusterNCount.new(
        cad, "NcClusterPhotozGaussGlobal", "NcClusterMassLnnormal"
    )

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

    rng = Ncm.RNG.seeded_new(None, 0)
    ncdata.init_from_sampling(mset, 270 * (math.pi / 180.0) ** 2, rng)

    return ncdata


def test_resample_matches_golden_snapshot() -> None:
    """Seed-0 resample reproduces the pinned output exactly."""
    ncdata = _resampled_ncount()

    assert ncdata.get_len() == 4242

    lnM_obs = np.array(ncdata.get_lnM_obs().dup_array())
    z_obs = np.array(ncdata.get_z_obs().dup_array())
    lnM_true = np.array(ncdata.get_lnM_true().dup_array())
    z_true = np.array(ncdata.get_z_true().dup_array())

    # No observable params for these models.
    assert ncdata.get_lnM_obs_params() is None
    assert ncdata.get_z_obs_params() is None

    assert _digest(lnM_obs) == "daffcd2e0548bb36"
    assert _digest(z_obs) == "42aa22695f8e8d7f"
    assert _digest(lnM_true) == "aa0e1d840a11608d"
    assert _digest(z_true) == "8d16a125aa223efa"
