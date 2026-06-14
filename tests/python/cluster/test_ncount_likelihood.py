#!/usr/bin/env python
#
# test_ncount_likelihood.py
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

"""Small integration test: simulate clusters and evaluate the likelihood.

This exercises the full path from the (generator-backed) NcDataClusterNCount
sampling into the number-counts likelihood, on a small survey area so it stays
fast. It checks that the likelihood evaluates to a finite value at the input
parameters and that it responds to a parameter change.
"""

import math

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

# Small area -> a handful of clusters, keeping the test fast.
AREA = 5.0 * (math.pi / 180.0) ** 2


def _dataset():
    """Simulate a small cluster number-count dataset and its mset."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.add_submodel(Nc.HIReionCamb())
    cosmo.add_submodel(Nc.HIPrimPowerLaw())

    dist = Nc.Distance.new(2.0)
    psml = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH())
    psml.require_kmin(1.0e-3)
    psml.require_kmax(1.0e3)
    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()

    mulf = Nc.MultiplicityFuncTinker.new_full(
        Nc.MultiplicityFuncMassDef.CRITICAL, 500.0
    )
    hmf = Nc.HaloMassFunction.new(dist, psf, mulf)

    # Proven, robust proxy combination (matches the C m2lnL suite): a richness
    # mass proxy and a no-distribution redshift (true z, no z integration).
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    cluster_m = ser.from_string("NcClusterMassAscaso")
    cluster_z = ser.from_string("NcClusterRedshiftNodist{'z-min':<0.1>, 'z-max':<1.0>}")
    cad = Nc.ClusterAbundance.new(hmf, None)

    ncdata = Nc.DataClusterNCount.new(
        cad, "NcClusterRedshiftNodist", "NcClusterMassAscaso"
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

    rng = Ncm.RNG.seeded_new(None, 0)
    ncdata.init_from_sampling(mset, AREA, rng)

    dset = Ncm.Dataset.new()
    dset.append_data(ncdata)

    return dset, mset, ncdata


def test_simulated_dataset_is_nonempty() -> None:
    """Sampling produces a usable number-counts dataset."""
    dset, _, ncdata = _dataset()

    assert dset.get_length() == 1
    assert ncdata.has_lnM_true()
    assert ncdata.has_z_true()
    assert len(ncdata.get_lnM_true().dup_array()) > 0


def test_likelihood_is_finite_at_truth() -> None:
    """The number-counts likelihood evaluates to a finite value."""
    dset, mset, _ = _dataset()

    assert dset.has_m2lnL_val()
    m2lnl = dset.m2lnL_val(mset)
    assert math.isfinite(m2lnl)


def test_likelihood_responds_to_parameter() -> None:
    """Changing a cosmological parameter changes the likelihood value."""
    dset, mset, _ = _dataset()

    m2lnl_truth = dset.m2lnL_val(mset)

    cosmo = mset.peek(Nc.HICosmo.id())
    cosmo.param_set_by_name("Omegac", 0.35)
    m2lnl_shifted = dset.m2lnL_val(mset)

    assert math.isfinite(m2lnl_shifted)
    assert not math.isclose(m2lnl_truth, m2lnl_shifted, rel_tol=1e-6)
