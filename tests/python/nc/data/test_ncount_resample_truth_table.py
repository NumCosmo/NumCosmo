#!/usr/bin/env python
#
# test_ncount_resample_truth_table.py
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

"""Truth table guarding NcDataClusterNCount.resample.

The resample draw order is RNG-sensitive. This test pins the output for a fixed
seed against a stored reference catalog, so the extraction of the sampling
pipeline into NcHaloCatalogGenerator can be verified to preserve behavior.

The reference is the resampled :class:`Nc.DataClusterNCount` itself, serialized
with NumCosmo's GVariant binfile format under
``data/truth_tables/cluster/ncount_resample_seed0.bin``. The comparison is
tolerance-based (not a byte hash) so it survives the sub-ULP drift of the
spline/transcendental evaluations across different libm/GSL/BLAS builds, while
still catching genuine regressions. Regenerate the truth table with::

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    ser.to_binfile(_resampled_ncount(), "<path>")
"""

import math

import numpy as np

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

TRUTH_TABLE_FILE = "truth_tables/cluster/ncount_resample_seed0.bin"
# Cross-stack libm/GSL/BLAS rounding drifts the proxy draws by a few ULP; this
# tolerance absorbs that while still flagging real changes in the draw order.
TRUTH_TABLE_RTOL = 1.0e-7
TRUTH_TABLE_ATOL = 1.0e-12


def _load_truth_table() -> Nc.DataClusterNCount:
    """Load the stored reference NcDataClusterNCount (seed 0)."""
    path = Ncm.cfg_get_data_filename(TRUTH_TABLE_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    table = ser.from_binfile(path)
    assert isinstance(table, Nc.DataClusterNCount)
    return table


def _resampled_ncount() -> Nc.DataClusterNCount:
    """Build and resample an NcDataClusterNCount deterministically (seed 0)."""
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


def test_resample_matches_truth_table() -> None:
    """Seed-0 resample reproduces the stored reference within tolerance."""
    ncdata = _resampled_ncount()
    table = _load_truth_table()

    assert ncdata.get_len() == table.get_len() == 4242

    # No observable params for these models.
    assert ncdata.get_lnM_obs_params() is None
    assert ncdata.get_z_obs_params() is None

    got = {}
    ref = {}
    for getter in ("get_lnM_obs", "get_z_obs", "get_lnM_true", "get_z_true"):
        got[getter] = np.array(getattr(ncdata, getter)().dup_array())
        ref[getter] = np.array(getattr(table, getter)().dup_array())
        np.testing.assert_allclose(
            got[getter],
            ref[getter],
            rtol=TRUTH_TABLE_RTOL,
            atol=TRUTH_TABLE_ATOL,
            err_msg=getter,
        )

    # The sampled (z, lnM) inherit the knots of the splines behind the abundance and move
    # at 1e-8 with them; the observable draws on top of them are the pipeline's own.
    # lnM_obs - lnM_true is sigma times the draw, and with z_bias = 0 so is
    # (z_obs - z_true) / (1 + z_true): both pinned exactly.
    np.testing.assert_allclose(
        got["get_lnM_obs"] - got["get_lnM_true"],
        ref["get_lnM_obs"] - ref["get_lnM_true"],
        rtol=0.0,
        atol=1.0e-12,
        err_msg="lnM draw",
    )
    np.testing.assert_allclose(
        (got["get_z_obs"] - got["get_z_true"]) / (1.0 + got["get_z_true"]),
        (ref["get_z_obs"] - ref["get_z_true"]) / (1.0 + ref["get_z_true"]),
        rtol=0.0,
        atol=1.0e-12,
        err_msg="z draw",
    )
