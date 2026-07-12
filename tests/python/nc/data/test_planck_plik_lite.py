#!/usr/bin/env python
#
# test_planck_plik_lite.py
#
# Fri June 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_plik_lite.py
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

"""Tests on the native NcDataPlanckPlikLite likelihood."""

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_lite import (
    PLIK_LITE_TT_RELPATH,
    NBIN_TT,
    find_baseline_file,
    build_plik_lite_tt,
)

Ncm.cfg_init()

_CLIK = find_baseline_file(PLIK_LITE_TT_RELPATH)
needs_data = pytest.mark.skipif(
    _CLIK is None, reason=f"plik_lite baseline data not found ({PLIK_LITE_TT_RELPATH})"
)


def test_type_is_resampleable_gausscov():
    """The native type is an NcmDataGaussCov (so it inherits resample)."""
    plik = Nc.DataPlanckPlikLite()
    assert isinstance(plik, Ncm.DataGaussCov)
    assert isinstance(plik, Ncm.Data)


@needs_data
def test_construct_and_serialize_roundtrip():
    """Build from clik data, serialize the dataset to binary, reload intact."""
    plik = build_plik_lite_tt(_CLIK)  # no Boltzmann: data-only
    assert plik.get_size() == NBIN_TT

    mean = Ncm.DataGaussCov.peek_mean(plik)
    cov = Ncm.DataGaussCov.peek_cov(plik)

    dset = Ncm.Dataset.new_array([plik])
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    variant = ser.to_variant(dset)
    ser.clear_instances()
    dset2 = ser.from_variant(variant)
    plik2 = dset2.peek_data(0)

    assert isinstance(plik2, Nc.DataPlanckPlikLite)
    assert plik2.get_size() == NBIN_TT
    assert plik2.get_property("calib-name") == "A_planck"

    mean2 = Ncm.DataGaussCov.peek_mean(plik2)
    cov2 = Ncm.DataGaussCov.peek_cov(plik2)
    for i in range(NBIN_TT):
        assert mean2.get(i) == mean.get(i)
        for j in range(0, NBIN_TT, 13):
            assert cov2.get(i, j) == cov.get(i, j)


@pytest.mark.app
@needs_data
def test_matches_clik_reference():
    """Native m2lnL equals the clik plik_lite reference on a CLASS cosmology.

    Both likelihoods share one CBE Boltzmann + cosmology + PlanckFI, so they see
    identical theory Cl. The clik reference reproduces the file's check_value, so
    an exact match validates the native reimplementation. Also checks resample.
    """
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import (
        mset_set_parameters,
        Planck18Types,
    )

    cbe = Nc.HIPertBoltzmannCBE.new()
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    ref = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TT_LITE, cbe
    )
    native = build_plik_lite_tt(_CLIK, cbe)

    ref.prepare(mset)
    native.prepare(mset)
    m2_ref = ref.m2lnL_val(mset)
    m2_nat = native.m2lnL_val(mset)

    assert m2_nat == pytest.approx(m2_ref, abs=1.0e-6)

    # resample: the capability the clik wrapper lacks.
    rng = Ncm.RNG.seeded_new(None, 123)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))
