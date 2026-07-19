#!/usr/bin/env python
#
# test_py_hireion.py
#
# Sun Sep 08 14:33:22 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_hireion.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numreion is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numreion is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Tests on NcHICosmoidem2 reionlogical model."""

import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc, GLib

Ncm.cfg_init()


def test_nc_hireion_camb():
    """Test NcHICosmoReionCamb initialization."""
    reion = Nc.HIReionCamb()
    assert reion is not None
    assert isinstance(reion, Nc.HIReionCamb)
    assert isinstance(reion, Nc.HIReion)
    assert isinstance(reion, Ncm.Model)


def test_nc_hireion_camb_z_to_tau():
    """Test NcHICosmoReionCamb z_to_tau."""
    reion = Nc.HIReionCamb()
    # reion must already be attached to its host cosmology before it can be
    # reparametrized -- z_to_tau() no longer takes a cosmo argument, it
    # reaches its own host via ncm_model_peek_host().
    cosmo = Nc.HICosmoDEXcdm(reion=reion)
    assert reion is not None
    assert isinstance(reion, Nc.HIReionCamb)
    assert isinstance(reion, Nc.HIReion)
    assert isinstance(reion, Ncm.Model)

    reion.z_to_tau()

    _ = reion["tau_reion"]
    z_re = reion.orig_param_get_by_name("z_re")
    reion["tau_reion"] = 0.123

    assert reion.orig_param_get_by_name("z_re") != z_re
    assert_allclose(reion["tau_reion"], 0.123)


def test_nc_hireion_camb_z_to_tau_round_trip():
    """Setting z_re, converting to tau_reion and back must recover z_re."""
    reion = Nc.HIReionCamb()
    _ = Nc.HICosmoDEXcdm(reion=reion)

    reion.orig_param_set_by_name("z_re", 9.5)
    z_re = reion.orig_param_get_by_name("z_re")

    reion.z_to_tau()
    tau = reion["tau_reion"]

    reion["tau_reion"] = tau
    assert_allclose(reion.orig_param_get_by_name("z_re"), z_re, rtol=1.0e-6)


def test_nc_hireion_camb_z_to_tau_cross_check():
    """The reparam's tau/z_re conversion must match the direct cosmo API."""
    reion = Nc.HIReionCamb()
    cosmo = Nc.HICosmoDEXcdm(reion=reion)

    tau_direct = reion.get_tau(cosmo)

    reion.z_to_tau()
    assert_allclose(reion["tau_reion"], tau_direct, rtol=1.0e-6)

    z_re_direct = reion.calc_z_from_tau(cosmo, tau_direct)
    reion["tau_reion"] = tau_direct
    assert_allclose(reion.orig_param_get_by_name("z_re"), z_re_direct, rtol=1.0e-6)


def test_nc_hireion_camb_z_to_tau_no_host():
    """z_to_tau()/set_z_from_tau() must raise, not crash, without a host."""
    reion = Nc.HIReionCamb()

    with pytest.raises(GLib.Error, match="not attached to a host cosmology"):
        reion.z_to_tau()

    with pytest.raises(GLib.Error, match="not attached to a host cosmology"):
        reion.set_z_from_tau(0.06)


def test_nc_hireion_camb_z_to_tau_reflects_current_host():
    """The reparam must always use reion's current host, not a stale snapshot.

    This is the actual historical bug: the old implementation smuggled a
    caller-supplied `cosmo` reference through a construct-only property, so
    the reparam could silently keep converting against a cosmology other
    than the one reion is really paired with. Now `ncm_model_peek_host()`
    always resolves to the one, permanently-bound host, with no separate
    reference for a caller to get out of sync.
    """
    reion = Nc.HIReionCamb()
    cosmo = Nc.HICosmoDEXcdm(reion=reion)

    reion.z_to_tau()
    reion["tau_reion"] = 0.06
    z_re_before = reion.orig_param_get_by_name("z_re")

    # Changing the host's own background parameters after the reparam is
    # attached must change the tau <-> z_re conversion -- proving it isn't
    # cached from attach time.
    cosmo.param_set_by_name("H0", cosmo.param_get_by_name("H0") * 1.1)
    reion["tau_reion"] = 0.06
    z_re_after = reion.orig_param_get_by_name("z_re")

    assert z_re_before != z_re_after
    assert_allclose(reion["tau_reion"], 0.06, rtol=1.0e-6)


def test_nc_hireion_camb_mset_batch_update_ordering():
    """Regression test for the NcmMSet batch fparam-set ordering fix.

    A submodel reparam that reads its host (like the tau reparam) is only
    correct if the host's own parameters are already committed by the
    time it reads them. `NcmMSet.fparams_set_vector()` sets every free
    parameter's raw value across *every* model first, then fires each
    model's update (hence any attached reparam) -- but that second pass
    must process the host before its submodels, or the submodel's reparam
    reads stale host state. This sets cosmo's H0 and reion's tau_reion in
    a single batch call and checks the resulting z_re reflects the *new*
    H0, not the pre-batch one.
    """
    reion = Nc.HIReionCamb()
    cosmo = Nc.HICosmoDEXcdm(reion=reion)
    reion.z_to_tau()

    found, h0_pid = cosmo.param_index_from_name("H0")
    assert found
    cosmo.param_set_ftype(h0_pid, Ncm.ParamType.FREE)

    found, tau_pid = reion.param_index_from_name("tau_reion")
    assert found
    reion.param_set_ftype(tau_pid, Ncm.ParamType.FREE)

    mset = Ncm.MSet.new_array([cosmo])
    mset.prepare_fparam_map()
    assert mset.fparams_len() == 2

    h0_new = cosmo.param_get_by_name("H0") * 1.1
    tau_target = 0.06

    x = Ncm.Vector.new(2)
    x.set(0, h0_new)
    x.set(1, tau_target)
    mset.fparams_set_vector(x)

    assert_allclose(cosmo.param_get_by_name("H0"), h0_new, rtol=1.0e-12)
    assert_allclose(reion["tau_reion"], tau_target, rtol=1.0e-6)

    # z_re must have been computed against the NEW H0 (cosmo, already
    # updated above), not whatever H0 was before this batch call.
    z_re_actual = reion.orig_param_get_by_name("z_re")
    z_re_expected_new_h0 = reion.calc_z_from_tau(cosmo, tau_target)

    assert_allclose(z_re_actual, z_re_expected_new_h0, rtol=1.0e-10)
