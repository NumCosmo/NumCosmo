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

from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

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
    cosmo = Nc.HICosmoDEXcdm()
    assert reion is not None
    assert isinstance(reion, Nc.HIReionCamb)
    assert isinstance(reion, Nc.HIReion)
    assert isinstance(reion, Ncm.Model)

    reion.z_to_tau(cosmo)

    _ = reion["tau_reion"]
    z_re = reion.orig_param_get_by_name("z_re")
    reion["tau_reion"] = 0.123

    assert reion.orig_param_get_by_name("z_re") != z_re
    assert_allclose(reion["tau_reion"], 0.123)
