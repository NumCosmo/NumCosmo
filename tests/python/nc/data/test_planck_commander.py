#!/usr/bin/env python
#
# test_planck_commander.py
#
# Wed July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_commander.py
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

"""Tests on the native NcDataPlanckCommander (low-ell gibbs_gauss) likelihood."""

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_commander import COMMANDER_RELPATH, build_commander
from numcosmo_py.experiments.planck_lite import find_baseline_file

Ncm.cfg_init()

_CLIK = find_baseline_file(COMMANDER_RELPATH)
needs_data = pytest.mark.skipif(
    _CLIK is None, reason=f"commander baseline data not found ({COMMANDER_RELPATH})"
)


def test_type_is_ncmdata():
    """The native type is an NcmData subclass."""
    assert issubclass(Nc.DataPlanckCommander, Ncm.Data)


@needs_data
def test_construct_and_serialize_roundtrip():
    """Build from the commander clik data, serialize and reload intact."""
    cmd = build_commander(_CLIK)  # data-only (no Boltzmann)
    assert cmd.get_length() == 28

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dset = Ncm.Dataset.new_array([cmd])
    cmd2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(cmd2, Nc.DataPlanckCommander)
    assert cmd2.get_length() == 28


@pytest.mark.app
@needs_data
def test_matches_clik_reference():
    """Native commander m2lnL agrees with the clik low-ell TT reference.

    Both share one CBE Boltzmann + cosmology + PlanckFICorTT (for A_planck).
    The native Gaussianized-Blackwell-Rao assembly (spline transform, Gaussian
    in x-space, spline Jacobian and offset) reproduces clik to ~1e-7.
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

    ref = Nc.DataPlanckLKL.full_new_id(Nc.DataPlanckLKLType.BASELINE_18_LOWL_TT, cbe)
    native = build_commander(_CLIK, cbe)

    ref.prepare(mset)
    native.prepare(mset)
    m2_ref = ref.m2lnL_val(mset)
    m2_nat = native.m2lnL_val(mset)

    assert m2_nat == pytest.approx(m2_ref, rel=1.0e-6)

    rng = Ncm.RNG.seeded_new(None, 42)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))


@pytest.mark.app
@needs_data
def test_clik_pi_compat_bit_identical():
    """With clik_pi_compat, the native m2lnL matches clik bit-for-bit.

    clik's gibbs uses a single-precision pi in the Dl = Cl*l(l+1)/2pi conversion;
    reproducing that removes the ~2e-6 residual of the (more accurate) default.
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

    ref = Nc.DataPlanckLKL.full_new_id(Nc.DataPlanckLKLType.BASELINE_18_LOWL_TT, cbe)
    native = build_commander(_CLIK, cbe, clik_pi_compat=True)

    ref.prepare(mset)
    native.prepare(mset)
    assert native.m2lnL_val(mset) == ref.m2lnL_val(mset)
