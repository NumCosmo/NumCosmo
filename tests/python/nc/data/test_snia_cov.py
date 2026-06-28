#!/usr/bin/env python
#
# test_py_data_snia_cov.py
#
# Sun Apr 21 11:15:07 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_snia_cov.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcmDataSNIACov class."""

import pytest
import numpy as np
from gi.repository import GLib

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_constructor():
    """Test constructor."""
    snia_cov = Nc.DataSNIACov.new(False, 0)
    assert snia_cov is not None
    assert isinstance(snia_cov, Nc.DataSNIACov)

    snia_cov2 = snia_cov.ref()
    assert snia_cov2 == snia_cov


@pytest.mark.parametrize(
    "cov_id",
    [
        Nc.DataSNIAId.COV_DES_Y5_STATONLY,
        Nc.DataSNIAId.COV_DES_Y5_STAT_SYS,
        Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT,
        Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL,
        Nc.DataSNIAId.COV_PANTHEON,
        Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_STAT,
        Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_SYS_STAT,
        Nc.DataSNIAId.COV_SNLS3_STAT_ONLY,
        Nc.DataSNIAId.COV_SNLS3_SYS_STAT,
    ],
)
def test_constructor_catalog_id(cov_id):
    """Test constructor with catalog id."""
    snia_cov = Nc.DataSNIACov.new_from_cat_id(cov_id, True)
    assert snia_cov is not None
    assert isinstance(snia_cov, Nc.DataSNIACov)

    assert snia_cov.get_dof() > 0
    assert snia_cov.get_length() > 0


@pytest.mark.acceptance
def test_recover_best_fit_PantheonPlus():
    """Test recover best fit PantheonPlus."""
    cov_id = Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_SYS_STAT
    snia_data = Nc.DataSNIACov.new_from_cat_id(cov_id, True)
    snia_data = snia_data.apply_filter_sh0es_z(0.01, True)
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo.param_set_desc("w", {"fit": True})
    cosmo.param_set_desc("Omegac", {"fit": True})
    cosmo.param_set_desc("H0", {"fit": True})

    mset = Ncm.MSet.new_array([cosmo])
    dset = Ncm.Dataset.new()

    dist = Nc.Distance.new(10.0)

    snia_model = Nc.SNIADistCov.new_by_id(dist, cov_id)
    mset.set(snia_model)

    dset.append_data(snia_data)

    likelihood = Ncm.Likelihood.new(dset)
    best_fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT,
        "ln-neldermead",
        likelihood,
        mset,
        Ncm.FitGradType.NUMDIFF_CENTRAL,
    )

    best_fit.run(Ncm.FitRunMsgs.NONE)
    #  Pantheon+ Sh0es best-fit Flat-wcdm (values obtained using their code)
    np.testing.assert_allclose(
        [cosmo["Omegac"], cosmo["w"], cosmo["H0"]],
        [0.253, -0.91, 73.44],
        rtol=1e-2,
        atol=1e-2,
    )


@pytest.mark.acceptance
def test_recover_best_fit_DES_Y5():
    """Test recover best fit DES Y5."""
    cov_id = Nc.DataSNIAId.COV_DES_Y5_STAT_SYS
    snia_data = Nc.DataSNIACov.new_from_cat_id(cov_id, True)
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo.param_set_desc("w", {"fit": True})
    cosmo.param_set_desc("Omegac", {"fit": True})

    mset = Ncm.MSet.new_array([cosmo])
    dset = Ncm.Dataset.new()

    dist = Nc.Distance.new(10.0)

    snia_model = Nc.SNIADistCov.new_by_id(dist, cov_id)
    mset.set(snia_model)

    dset.append_data(snia_data)

    likelihood = Ncm.Likelihood.new(dset)
    best_fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT,
        "ln-neldermead",
        likelihood,
        mset,
        Ncm.FitGradType.NUMDIFF_CENTRAL,
    )

    best_fit.run(Ncm.FitRunMsgs.NONE)
    # DES Y5 best fit Flat-wcdm from arXiv:2401.02929: Omegam = 0.264, w = -0.80 /
    # Omegab = 0.0432
    np.testing.assert_allclose(
        [cosmo["Omegac"], cosmo["w"]], [0.2208, -0.80], rtol=1e-2, atol=1e-2
    )


def test_resample_type_default_is_auto():
    """The resample strategy defaults to AUTO."""
    snia_cov = Nc.DataSNIACov.new(False, 0)
    assert snia_cov.get_resample_type() == Nc.DataSNIACovResample.AUTO


def test_resample_type_from_cov_always_available():
    """FROM_COV is valid regardless of whether the light-curve covariance exists."""
    snia_cov = Nc.DataSNIACov.new_from_cat_id(
        Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_SYS_STAT, False
    )
    assert not snia_cov.has_complete_cov()
    assert snia_cov.set_resample_type(Nc.DataSNIACovResample.FROM_COV)
    assert snia_cov.get_resample_type() == Nc.DataSNIACovResample.FROM_COV


def test_resample_type_lightcurve_requires_complete_cov():
    """FROM_LIGHTCURVE is rejected (with a GError) on a mu-only dataset."""
    snia_cov = Nc.DataSNIACov.new_from_cat_id(
        Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_SYS_STAT, False
    )
    assert not snia_cov.has_complete_cov()
    with pytest.raises(GLib.Error, match="light-curve"):
        snia_cov.set_resample_type(Nc.DataSNIACovResample.FROM_LIGHTCURVE)
    # The rejected request leaves the previous value untouched.
    assert snia_cov.get_resample_type() == Nc.DataSNIACovResample.AUTO


def test_resample_type_lightcurve_allowed_when_complete():
    """FROM_LIGHTCURVE is accepted when the dataset ships the light-curve cov."""
    snia_cov = Nc.DataSNIACov.new_from_cat_id(
        Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL, False
    )
    assert snia_cov.has_complete_cov()
    assert snia_cov.set_resample_type(Nc.DataSNIACovResample.FROM_LIGHTCURVE)
    assert snia_cov.get_resample_type() == Nc.DataSNIACovResample.FROM_LIGHTCURVE


def test_resample_type_survives_serialization():
    """The chosen strategy round-trips through serialization (raw, unvalidated).

    The property channel stores the value verbatim (validation happens in the
    public setter / at resample time), so an experiment file can carry the choice
    even before the catalog is loaded. An empty catalog is used here because a
    full-catalog dup_obj trips an unrelated abs_mag_set length assertion.
    """
    snia_cov = Nc.DataSNIACov.new(False, 0)
    snia_cov.set_property("resample-type", Nc.DataSNIACovResample.FROM_LIGHTCURVE)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dup = ser.dup_obj(snia_cov)
    assert dup.get_resample_type() == Nc.DataSNIACovResample.FROM_LIGHTCURVE


@pytest.mark.parametrize(
    "resample_type",
    [
        Nc.DataSNIACovResample.AUTO,
        Nc.DataSNIACovResample.FROM_COV,
        Nc.DataSNIACovResample.FROM_LIGHTCURVE,
    ],
)
def test_resample_runs_in_each_mode(resample_type):
    """Each resample strategy produces a finite realization on a complete dataset."""
    snia_cov = Nc.DataSNIACov.new_from_cat_id(
        Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL, False
    )
    snia_cov.set_resample_type(resample_type)

    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    dcov = Nc.SNIADistCov.new_by_id(
        Nc.Distance.new(10.0), Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL
    )
    mset = Ncm.MSet.new_array([cosmo, dcov])
    mset.prepare_fparam_map()

    rng = Ncm.RNG.seeded_new(None, 123)
    snia_cov.resample(mset, rng)
    assert np.all(np.isfinite(snia_cov.peek_mag().dup_array()))
