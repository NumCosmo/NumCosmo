#!/usr/bin/env python
#
# test_ncounts_gauss_varying_sij.py
#
# Thu Aug 13 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_ncounts_gauss_varying_sij.py
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

"""NcDataClusterNCountsGauss with a cosmology-dependent super-sample covariance.

Setting the `ssc-sij` property replaces the fixed `s-matrix` by a #NcXcorSSCSij
recomputed at every prepare, so the whole covariance varies along a chain
instead of being pinned at a fiducial cosmology.

Two invariants matter for the comparison this was built to support, and are
tested here: at the fiducial cosmology the varying setup must reproduce the
fixed one exactly, so the two analyses differ only where they are meant to; and
the resample matrix must stay frozen, so both fit the same mock realization.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

#: Seven top-hat redshift bins and a single mass bin.
Z_EDGES = np.linspace(0.1, 0.8, 8)
LNM_EDGES = [34.0, 35.0]

FIDUCIAL_OMEGAC = 0.2612


def _make_cosmo(omega_c: float = FIDUCIAL_OMEGAC) -> Nc.HICosmo:
    """Build the test cosmology."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("H0", 67.81)
    cosmo.param_set_by_name("Omegac", omega_c)
    cosmo.param_set_by_name("Omegab", 0.0486)
    cosmo.param_set_by_name("w", -1.0)
    cosmo.param_set_by_name("Omegak", 0.0)
    cosmo.add_submodel(Nc.HIReionCamb.new())
    cosmo.add_submodel(Nc.HIPrimPowerLaw.new())

    return cosmo


def _build(varying: bool, use_norma: bool = True):
    """Build a likelihood over the cluster counts, with or without varying S_ij."""
    cosmo = _make_cosmo()

    dist = Nc.Distance.new(3.0)
    ps_ml = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH.new())
    ps_ml.require_kmin(1.0e-6)
    ps_ml.require_kmax(1.0e3)
    psf = Ncm.PowspecFilter.new(ps_ml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()

    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.MEAN)
    mfp = Nc.HaloMassFunction.new(dist, psf, mulf)
    cad = Nc.ClusterAbundance.new(mfp, Nc.HaloBiasPS.new(mfp))

    clusterz = Nc.ClusterRedshiftNodist(z_max=Z_EDGES[-1], z_min=Z_EDGES[0])
    clusterm = Nc.ClusterMassNodist(lnM_max=LNM_EDGES[-1], lnM_min=LNM_EDGES[0])
    mset = Ncm.MSet.new_array([cosmo, clusterm, clusterz])

    z_edges = Ncm.Vector.new_array(Z_EDGES.tolist())

    data = Nc.DataClusterNCountsGauss.new(cad)
    data.set_z_obs(z_edges)
    data.set_lnM_obs(Ncm.Vector.new_array(LNM_EDGES))
    data.set_has_ssc(True)
    data.use_norma(use_norma)

    # One fiducial matrix, used both as the fixed covariance and as the frozen
    # resample covariance, so the mock is the same in both setups.
    s_fiducial = Nc.XcorSSCSij.new(dist, ps_ml, z_edges).eval(cosmo)
    data.set_s_matrix(s_fiducial)
    data.set_resample_s_matrix(s_fiducial)

    if varying:
        data.set_ssc_sij(Nc.XcorSSCSij.new(dist, ps_ml, z_edges))

    dset = Ncm.Dataset.new()
    dset.append_data(data)
    likelihood = Ncm.Likelihood.new(dset)

    data.resample(mset, Ncm.RNG.seeded_new(None, 1234))

    return cosmo, mset, data, likelihood


def _matrix_to_np(matrix: Ncm.Matrix) -> np.ndarray:
    """Convert an #NcmMatrix into a numpy array."""
    return np.array(
        [
            [matrix.get(i, j) for j in range(matrix.ncols())]
            for i in range(matrix.nrows())
        ]
    )


def test_matches_fixed_matrix_at_fiducial() -> None:
    """At the fiducial cosmology the varying setup reproduces the fixed one.

    The fixed matrix was computed by the same calculator at the same cosmology,
    so any difference here would be a wiring bug rather than physics.
    """
    _, mset_fixed, _, lh_fixed = _build(varying=False)
    _, mset_varying, _, lh_varying = _build(varying=True)

    assert_allclose(
        lh_varying.m2lnL_val(mset_varying),
        lh_fixed.m2lnL_val(mset_fixed),
        rtol=1.0e-12,
    )


def test_differs_away_from_fiducial() -> None:
    """Away from the fiducial cosmology the two likelihoods separate."""
    cosmo_fixed, mset_fixed, _, lh_fixed = _build(varying=False)
    cosmo_varying, mset_varying, _, lh_varying = _build(varying=True)

    for omega_c in (0.20, 0.32):
        cosmo_fixed.param_set_by_name("Omegac", omega_c)
        cosmo_varying.param_set_by_name("Omegac", omega_c)

        assert lh_varying.m2lnL_val(mset_varying) != lh_fixed.m2lnL_val(mset_fixed)


def test_sij_tracks_cosmology() -> None:
    """The fitting S_ij moves when the cosmology moves."""
    cosmo, mset, data, likelihood = _build(varying=True)
    ssc_sij = data.get_ssc_sij()

    likelihood.m2lnL_val(mset)
    at_fiducial = _matrix_to_np(ssc_sij.peek_matrix())

    cosmo.param_set_by_name("Omegac", 0.20)
    likelihood.m2lnL_val(mset)
    moved = _matrix_to_np(ssc_sij.peek_matrix())

    assert not np.allclose(at_fiducial, moved)


def test_resample_matrix_stays_frozen() -> None:
    """Moving the cosmology leaves the resample matrix untouched.

    Both the fixed and the varying analyses must fit the same mock, so only the
    fitting covariance may follow the cosmology.
    """
    cosmo, mset, data, likelihood = _build(varying=True)

    before = _matrix_to_np(data.get_resample_s_matrix())

    cosmo.param_set_by_name("Omegac", 0.20)
    likelihood.m2lnL_val(mset)

    after = _matrix_to_np(data.get_resample_s_matrix())

    assert_allclose(after, before, rtol=0.0, atol=0.0)


def test_fixed_s_matrix_untouched() -> None:
    """The fixed s-matrix is left alone; the calculator does not overwrite it."""
    cosmo, mset, data, likelihood = _build(varying=True)

    before = _matrix_to_np(data.get_s_matrix())

    cosmo.param_set_by_name("Omegac", 0.20)
    likelihood.m2lnL_val(mset)

    assert_allclose(_matrix_to_np(data.get_s_matrix()), before, rtol=0.0, atol=0.0)


def test_get_ssc_sij_round_trips() -> None:
    """The calculator set on the data comes back out, and can be unset."""
    _, _, data, _ = _build(varying=True)

    assert isinstance(data.get_ssc_sij(), Nc.XcorSSCSij)

    data.set_ssc_sij(None)

    assert data.get_ssc_sij() is None


def test_set_ssc_sij_keeps_the_data_initialized() -> None:
    """Attaching a calculator must not reset the data's init state.

    The property is restored on deserialization like any other, so a setter
    that cleared `init` would leave a loaded experiment uninitialized and trip
    the assertion in ncm_data_prepare().
    """
    _, _, data, _ = _build(varying=True)

    assert data.is_init()

    data.set_ssc_sij(Nc.XcorSSCSij.new(data.get_ssc_sij().props.dist,
                                       data.get_ssc_sij().props.powspec,
                                       Ncm.Vector.new_array(Z_EDGES.tolist())))

    assert data.is_init()


def test_round_trips_through_a_serialized_experiment() -> None:
    """A serialized dataset reloads with a working calculator.

    This is the path the CLI takes: `numcosmo generate` writes the dataset and
    `run mc` / `run mcmc` read it back, so the calculator has to survive the
    trip and still be usable.
    """
    _, mset, data, likelihood = _build(varying=True)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dset = ser.from_string(ser.to_string(Ncm.Dataset.new_array([data]), True))
    reloaded = dset.peek_data(0)

    assert reloaded.is_init()
    assert isinstance(reloaded.get_ssc_sij(), Nc.XcorSSCSij)

    reloaded_lh = Ncm.Likelihood.new(dset)

    assert_allclose(
        reloaded_lh.m2lnL_val(mset), likelihood.m2lnL_val(mset), rtol=1.0e-12
    )


def test_serialization_keeps_the_calculator() -> None:
    """The data object serializes with its calculator attached."""
    _, _, data, _ = _build(varying=True)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dup = ser.from_string(ser.to_string(data, True))

    ssc_sij = dup.get_ssc_sij()

    assert isinstance(ssc_sij, Nc.XcorSSCSij)
    assert ssc_sij.get_nbins() == len(Z_EDGES) - 1


@pytest.mark.parametrize("omega_c", [0.20, 0.23, 0.29, 0.32])
def test_varying_covariance_shifts_m2lnl_by_a_small_amount(omega_c: float) -> None:
    """The S_ij contribution to -2lnL is a small correction, not a rewrite.

    This is the quantity the whole comparison is about: if freezing S_ij were a
    large effect the reference analysis would be wrong, and if it were exactly
    zero there would be nothing to demonstrate.
    """
    cosmo_fixed, mset_fixed, _, lh_fixed = _build(varying=False)
    cosmo_varying, mset_varying, _, lh_varying = _build(varying=True)

    cosmo_fixed.param_set_by_name("Omegac", omega_c)
    cosmo_varying.param_set_by_name("Omegac", omega_c)

    fixed = lh_fixed.m2lnL_val(mset_fixed)
    varying = lh_varying.m2lnL_val(mset_varying)

    assert 0.0 < abs(varying - fixed) < 1.0
