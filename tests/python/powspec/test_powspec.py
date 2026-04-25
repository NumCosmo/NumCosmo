#!/usr/bin/env python
#
# test_powspec.py
#
# Wed Feb 14 13:44:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_powspec.py
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

"""Unit tests for NumCosmo power-spectra without CCL dependency."""

# pylint: disable=redefined-outer-name

import pytest

import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

pytestmark = pytest.mark.powspec

Ncm.cfg_init()


@pytest.fixture
def nc_cosmo() -> Nc.HICosmo:
    """Create a simple NumCosmo cosmology object."""
    cosmo = Nc.HICosmoDEXcdm()
    reion = Nc.HIReionCamb.new()
    prim = Nc.HIPrimPowerLaw.new()

    cosmo.add_submodel(reion)
    cosmo.add_submodel(prim)

    # Set default parameters
    cosmo.props.H0 = 70.0
    cosmo.props.Omegab = 0.05
    cosmo.props.Omegac = 0.25
    cosmo.props.Omegax = 0.70
    cosmo.props.Tgamma0 = 2.7255
    cosmo.props.w = -1.0

    prim.props.n_SA = 0.96
    prim.props.ln10e10ASA = 3.0

    return cosmo


def test_powspec_class(nc_cosmo: Nc.HICosmo) -> None:
    """Test the power spectrum class."""
    ps_ml = Nc.PowspecMLCBE.new()
    cbe = ps_ml.peek_cbe()
    cbe.use_ppf(True)

    ps_ml.prepare(nc_cosmo)

    ik_min = ps_ml.get_intern_k_min()
    ik_max = ps_ml.get_intern_k_max()
    kmin = ps_ml.get_kmin()
    kmax = ps_ml.get_kmax()

    assert ik_min > 0
    assert ik_max > ik_min

    k_a = np.geomspace(kmin, kmax, 100)

    ps_a = np.array([ps_ml.eval(nc_cosmo, 0.0, k) for k in k_a])

    assert np.all(np.isfinite(ps_a))


def test_powspec_class_deriv_z(nc_cosmo: Nc.HICosmo) -> None:
    """Test the power spectrum class."""
    ps_ml = Nc.PowspecMLCBE.new()
    cbe = ps_ml.peek_cbe()
    cbe.use_ppf(True)

    ps_ml.prepare(nc_cosmo)

    ik_min = ps_ml.get_intern_k_min()
    ik_max = ps_ml.get_intern_k_max()
    kmin = ps_ml.get_kmin()
    kmax = ps_ml.get_kmax()

    assert kmin < ik_min
    assert kmax > ik_max
    assert ik_min > 0
    assert ik_max > ik_min

    k_a = np.geomspace(kmin, kmax, 100)

    diff = Ncm.Diff.new()

    for k in k_a:
        dps, _ = diff.rc_d1_1_to_1(0.2, lambda z, k0: ps_ml.eval(nc_cosmo, z, k0), k)
        assert_allclose(dps, ps_ml.deriv_z(nc_cosmo, 0.2, k), atol=0.0, rtol=1.0e-10)
