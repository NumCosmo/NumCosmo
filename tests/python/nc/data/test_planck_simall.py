#!/usr/bin/env python
#
# test_planck_simall.py
#
# Wed July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_simall.py
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

"""Tests on the native NcDataPlanckSimall (low-ell SimAll) likelihood."""

import pytest

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_simall import (
    SIMALL_EE_RELPATH,
    SIMALL_BB_RELPATH,
    SIMALL_EEBB_RELPATH,
    build_simall,
)
from numcosmo_py.experiments.planck_lite import find_baseline_file

Ncm.cfg_init()

_EE = find_baseline_file(SIMALL_EE_RELPATH)
needs_ee = pytest.mark.skipif(_EE is None, reason="simall EE data not found")


def test_type_is_ncmdata():
    """The native type is an NcmData subclass."""
    assert issubclass(Nc.DataPlanckSimall, Ncm.Data)


@needs_ee
def test_construct_and_serialize_roundtrip():
    """Build from the simall EE clik data, serialize and reload intact."""
    simall = build_simall(_EE)
    assert simall.get_length() == 28

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dset = Ncm.Dataset.new_array([simall])
    simall2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(simall2, Nc.DataPlanckSimall)
    assert simall2.get_length() == 28


@pytest.mark.app
@pytest.mark.parametrize(
    "relpath,enum,length",
    [
        (SIMALL_EE_RELPATH, "BASELINE_18_LOWL_EE", 28),
        (SIMALL_BB_RELPATH, "BASELINE_18_LOWL_BB", 28),
        (SIMALL_EEBB_RELPATH, "BASELINE_18_LOWL_EEBB", 56),
    ],
)
def test_matches_clik_reference(relpath, enum, length):
    """Native simall m2lnL matches the clik low-ell reference (bit-exact).

    The tabulated log-probability with floor indexing reproduces clik exactly.
    """
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    clik = find_baseline_file(relpath)
    if clik is None:
        pytest.skip(f"simall data not found ({relpath})")

    cbe = Nc.HIPertBoltzmannCBE.new()
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    ref = Nc.DataPlanckLKL.full_new_id(getattr(Nc.DataPlanckLKLType, enum), cbe)
    native = build_simall(clik, cbe)
    assert native.get_length() == length

    ref.prepare(mset)
    native.prepare(mset)
    assert native.m2lnL_val(mset) == pytest.approx(ref.m2lnL_val(mset), rel=1.0e-13)
