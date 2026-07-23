#!/usr/bin/env python
#
# test_planck_lensing.py
#
# Thu July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_lensing.py
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

"""Tests on the native NcDataPlanckLensing (Planck 2018 CMB lensing) likelihood."""

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_lensing import (
    LENSING_FULL_RELPATH,
    LENSING_MARGED_RELPATH,
    build_lensing,
)
from numcosmo_py.experiments.planck_lite import find_baseline_file

Ncm.cfg_init()

# (relpath, DataPlanckLKL enum, whether the file needs CMB spectra) per flavor.
_FLAVORS = {
    "marged": (
        LENSING_MARGED_RELPATH,
        Nc.DataPlanckLKLType.BASELINE_18_LENSING_CMB_MARGINALIZED,
    ),
    "full": (
        LENSING_FULL_RELPATH,
        Nc.DataPlanckLKLType.BASELINE_18_LENSING,
    ),
}


def _clik(relpath):
    return find_baseline_file(relpath)


def _mset():
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()
    return mset


def test_type_is_gauss():
    """The native type is an NcmDataGauss subclass."""
    assert issubclass(Nc.DataPlanckLensing, Ncm.DataGauss)


@pytest.mark.parametrize("flavor", list(_FLAVORS))
def test_construct_and_serialize_roundtrip(flavor):
    """Build from the lensing clik data, serialize and reload intact."""
    relpath, _ = _FLAVORS[flavor]
    if _clik(relpath) is None:
        pytest.skip(f"lensing baseline data not found ({relpath})")

    lens = build_lensing(_clik(relpath))  # data-only (no Boltzmann)
    assert lens.get_length() == 9

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dset = Ncm.Dataset.new_array([lens])
    lens2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(lens2, Nc.DataPlanckLensing)
    assert lens2.get_length() == 9


@pytest.mark.app
@pytest.mark.parametrize("flavor", list(_FLAVORS))
def test_matches_clik_reference(flavor):
    """Native lensing m2lnL agrees with the clik lensing reference.

    Both share one CBE Boltzmann + cosmology + PlanckFICorTT (for A_planck).
    The native inverse-covariance Gaussian assembly of the band-power model
    reproduces the clik dts_lensing likelihood to machine precision.
    """
    relpath, enum = _FLAVORS[flavor]
    if _clik(relpath) is None:
        pytest.skip(f"lensing baseline data not found ({relpath})")

    cbe = Nc.HIPertBoltzmannCBE.new()
    mset = _mset()

    ref = Nc.DataPlanckLKL.full_new_id(enum, cbe)
    native = build_lensing(_clik(relpath), cbe)

    ref.prepare(mset)
    native.prepare(mset)

    assert native.m2lnL_val(mset) == pytest.approx(ref.m2lnL_val(mset), rel=1.0e-11)

    rng = Ncm.RNG.seeded_new(None, 42)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))


@pytest.mark.app
def test_prepare_self_configures_boltzmann():
    """Lensing raises the Boltzmann targets/lmax itself on a bare CBE.

    Unlike plik_lite/smica/commander/simall (which rely on the caller to set the
    targets and lmax), NcDataPlanckLensing self-configures in prepare(). Drive a
    fresh CBE that only has the default TT target at a low lmax and check the full
    file lifts PHIPHI + TT/EE/TE up to its own lmax and still evaluates.
    """
    if _clik(LENSING_FULL_RELPATH) is None:
        pytest.skip(f"lensing baseline data not found ({LENSING_FULL_RELPATH})")

    cbe = Nc.HIPertBoltzmannCBE.new()  # default: TT only, lmax ~30
    mset = _mset()

    native = build_lensing(_clik(LENSING_FULL_RELPATH), cbe)
    native.prepare(mset)

    assert cbe.get_PHIPHI_lmax() >= 2500
    assert cbe.get_TT_lmax() >= 2500
    assert cbe.get_EE_lmax() >= 2500
    assert cbe.get_TE_lmax() >= 2500
    assert np.isfinite(native.m2lnL_val(mset))
