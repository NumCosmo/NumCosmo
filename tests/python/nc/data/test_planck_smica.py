#!/usr/bin/env python
#
# test_planck_smica.py
#
# Fri June 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_smica.py
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

"""Tests on the native NcDataPlanckSmica (plik high-ell SMICA) likelihood."""

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_smica import (
    PLIK_TT_RELPATH,
    PLIK_TTTEEE_RELPATH,
    build_smica_tt,
    build_smica_ttteee,
)
from numcosmo_py.experiments.planck_lite import find_baseline_file

Ncm.cfg_init()

_CLIK = find_baseline_file(PLIK_TT_RELPATH)
needs_data = pytest.mark.skipif(
    _CLIK is None, reason=f"plik baseline data not found ({PLIK_TT_RELPATH})"
)

_CLIK_TTTEEE = find_baseline_file(PLIK_TTTEEE_RELPATH)
needs_data_ttteee = pytest.mark.skipif(
    _CLIK_TTTEEE is None,
    reason=f"plik baseline data not found ({PLIK_TTTEEE_RELPATH})",
)


def test_type_is_resampleable_gauss():
    """The native type is an NcmDataGauss (so it inherits resample).

    Asserted on the type, not an instance: NcDataPlanckSmica is fully specified
    at construction (all structural properties are CONSTRUCT_ONLY and validated
    in constructed()), so it is deliberately not default-constructible.
    """
    assert issubclass(Nc.DataPlanckSmica, Ncm.DataGauss)


@needs_data
def test_construct_and_serialize_roundtrip():
    """Build from the plik clik data, serialize the dataset, reload intact."""
    smica = build_smica_tt(_CLIK)  # data-only (no Boltzmann)
    npt = smica.get_size()
    assert npt == 765

    mean = Ncm.DataGauss.peek_mean(smica)
    dset = Ncm.Dataset.new_array([smica])
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    smica2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(smica2, Nc.DataPlanckSmica)
    assert smica2.get_size() == npt
    mean2 = Ncm.DataGauss.peek_mean(smica2)
    for i in range(0, npt, 37):
        assert mean2.get(i) == mean.get(i)


@pytest.mark.app
@needs_data
def test_matches_clik_reference():
    """Native SMICA m2lnL agrees with the clik plik TT reference on a CLASS cosmology.

    Both share one CBE Boltzmann + cosmology + PlanckFICorTT. The native
    reimplementation reproduces the full R_q assembly (CMB + 10 foreground/
    calibration components) to machine precision.
    """
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    cbe = Nc.HIPertBoltzmannCBE.new()
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    ref = Nc.DataPlanckLKL.full_new_id(Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TT, cbe)
    native = build_smica_tt(_CLIK, cbe)

    ref.prepare(mset)
    native.prepare(mset)
    m2_ref = ref.m2lnL_val(mset)
    m2_nat = native.m2lnL_val(mset)

    assert m2_nat == pytest.approx(m2_ref, rel=1.0e-12)

    rng = Ncm.RNG.seeded_new(None, 42)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))


@needs_data_ttteee
def test_ttteee_construct_and_serialize_roundtrip():
    """Build the TTTEEE (m=6) SMICA from clik data; serialize/reload intact."""
    smica = build_smica_ttteee(_CLIK_TTTEEE)
    npt = smica.get_size()
    assert npt == 2289

    mean = Ncm.DataGauss.peek_mean(smica)
    dset = Ncm.Dataset.new_array([smica])
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    smica2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(smica2, Nc.DataPlanckSmica)
    assert smica2.get_size() == npt
    mean2 = Ncm.DataGauss.peek_mean(smica2)
    for i in range(0, npt, 53):
        assert mean2.get(i) == mean.get(i)


@pytest.mark.app
@needs_data_ttteee
def test_ttteee_matches_clik_reference():
    """Native TTTEEE (m=6) SMICA m2lnL agrees with the clik reference.

    Same NcDataPlanckSmica class as TT, now with the polarization channels and
    components (CMB TT/EE/TE, pwfe galactic dust, EE end-to-end noise, the
    icalTP two-term calibration mixing and totcalP). Machine precision.
    """
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    cbe = Nc.HIPertBoltzmannCBE.new()
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTTTEEE()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TTTEEE, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    ref = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TTTEEE, cbe
    )
    native = build_smica_ttteee(_CLIK_TTTEEE, cbe)

    ref.prepare(mset)
    native.prepare(mset)
    m2_ref = ref.m2lnL_val(mset)
    m2_nat = native.m2lnL_val(mset)

    assert m2_nat == pytest.approx(m2_ref, rel=1.0e-12)

    rng = Ncm.RNG.seeded_new(None, 42)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))
