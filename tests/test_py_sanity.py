#!/usr/bin/env python
#
# test_py_sanity.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sanity.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests for the Python bindings sanity."""

from numcosmo_py import Ncm
from numcosmo_py import Nc

import math


def test_cfg_init():
    """Test NumCosmo initialization."""
    Ncm.cfg_init()


def test_distances():
    """Test distances."""

    Ncm.cfg_init()

    cosmo = Nc.HICosmo.new_from_name(Nc.HICosmo, "NcHICosmoDEXcdm{'massnu-length':<1>}")
    cosmo.set_reparam(Nc.HICosmoDEReparamCMB.new(cosmo.len()))

    dist = Nc.Distance.new(2.0)

    cosmo.orig_param_set(Nc.HICosmoDESParams.H0, 70.00)
    cosmo.orig_param_set(Nc.HICosmoDESParams.OMEGA_C, 0.25)
    cosmo.orig_param_set(Nc.HICosmoDESParams.OMEGA_X, 0.70)
    cosmo.orig_param_set(Nc.HICosmoDESParams.T_GAMMA0, 2.72)
    cosmo.orig_param_set(Nc.HICosmoDESParams.OMEGA_B, 0.05)
    cosmo.orig_param_set(Nc.HICosmoDEXCDMSParams.W, -1.10)

    cosmo.orig_vparam_set(Nc.HICosmoDEVParams.M, 0, 0.06)
    dist.prepare(cosmo)

    N = 100
    RH_Mpc = cosmo.RH_Mpc()

    for i in range(N):
        z = 1.0 / (N - 1.0) * i
        Dc = dist.comoving(cosmo, z)
        dc = RH_Mpc * Dc

        assert math.isfinite(dc)
