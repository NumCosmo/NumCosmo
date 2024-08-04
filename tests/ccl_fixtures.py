#!/usr/bin/env python
#
# ccl_fixtures.py
#
# Thu Aug 01 11:44:12 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_powspec.py
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

"""Unit tests for NumCosmo powwer-spectra."""

from itertools import product
import pytest

import numpy as np

import pyccl

from numcosmo_py import Ncm
from numcosmo_py.ccl.nc_ccl import CCLParams

Ncm.cfg_init()


@pytest.fixture(name="k_a")
def fixture_k_a():
    """Fixture for k array."""
    return np.geomspace(1.0e-5, 1.0e1, 1000)


@pytest.fixture(name="z_a")
def fixture_z_a():
    """Fixture for z array."""
    return np.linspace(0.0, 5.0, 1000)


@pytest.fixture(
    name="ccl_cosmo_eh_linear",
    params=product([False, True], range(3)),
    ids=lambda x: f"high_prec={x[0]},index={x[1]}",
)
def fixture_ccl_cosmo_eh_linear(request) -> pyccl.Cosmology:
    """Fixture for CCL Cosmology."""
    Omega_c = 0.25
    Omega_b = 0.05
    Omega_k = 0.0
    h = 0.7
    n_s = 0.96
    Neff = 0.0
    sigma8 = 0.9

    Omega_v_vals = np.array([0.7, 0.65, 0.75])
    w0_vals = np.array([-1.0, -0.9, -1.1])
    wa_vals = np.array([0.0, 0.1, -0.1])

    high_prec, index = request.param

    Omega_k = 1.0 - Omega_c - Omega_b - Omega_v_vals[index]

    if high_prec:
        CCLParams.set_high_prec_params()
    else:
        CCLParams.set_default_params()

    ccl_cosmo = pyccl.Cosmology(
        Omega_c=Omega_c,
        Omega_b=Omega_b,
        Neff=Neff,
        h=h,
        sigma8=sigma8,
        n_s=n_s,
        Omega_k=Omega_k,
        w0=w0_vals[index],
        wa=wa_vals[index],
        transfer_function="eisenstein_hu",
        matter_power_spectrum="linear",
    )

    ccl_cosmo.high_precision = high_prec
    return ccl_cosmo


# CCL Halofit is too slow for high precision
@pytest.fixture(
    name="ccl_cosmo_eh_halofit",
    params=product([False], range(3)),
    ids=lambda x: f"high_prec={x[0]},index={x[1]}",
)
def fixture_ccl_cosmo_eh_halofit(request) -> pyccl.Cosmology:
    """Fixture for CCL Cosmology."""
    Omega_c = 0.25
    Omega_b = 0.05
    Omega_k = 0.0
    h = 0.7
    n_s = 0.96
    Neff = 0.0
    sigma8_vals = [0.9, 0.7, 0.6]

    Omega_v = 0.7
    w0 = -1.0
    wa = 0.0

    high_prec, index = request.param

    Omega_k = 1.0 - Omega_c - Omega_b - Omega_v

    if high_prec:
        CCLParams.set_high_prec_params()
    else:
        CCLParams.set_default_params()

    ccl_cosmo = pyccl.Cosmology(
        Omega_c=Omega_c,
        Omega_b=Omega_b,
        Neff=Neff,
        h=h,
        sigma8=sigma8_vals[index],
        n_s=n_s,
        Omega_k=Omega_k,
        w0=w0,
        wa=wa,
        transfer_function="eisenstein_hu",
        matter_power_spectrum="halofit",
    )

    ccl_cosmo.high_precision = high_prec
    return ccl_cosmo
