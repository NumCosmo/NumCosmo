#!/usr/bin/env python
#
# test_py_ccl_background.py
#
# Thu Aug 01 14:07:10 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_background.py
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

import pytest
from numpy.testing import assert_allclose
import numpy as np
import pyccl

import numcosmo_py.cosmology as ncpy
from numcosmo_py import Ncm
from numcosmo_py.ccl.comparison import (
    compare_power_spectrum_linear,
    compare_power_spectrum_nonlinear,
    compare_sigma_r,
)
from .fixtures_ccl import (  # pylint: disable=unused-import # noqa: F401
    fixture_k_a,
    fixture_z_a,
    fixture_ccl_cosmo_eh_linear,
    fixture_ccl_cosmo_eh_halofit,
    fixture_nc_cosmo_eh_linear,
    fixture_nc_cosmo_eh_halofit,
)

Ncm.cfg_init()

Z_ARRAY = np.linspace(0.0, 5.0, 10)
Z_IDS = [f"z={z:.2f}" for z in Z_ARRAY]


@pytest.mark.parametrize("z", Z_ARRAY, ids=Z_IDS)
def test_power_spectrum_linear(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology, z
) -> None:
    """Compare NumCosmo and CCL linear power spectrum."""
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-8
        if z > 2.0:
            rtol = 1.0e-7
    else:
        rtol = 1.0e-4

    k_test = np.geomspace(5.0e-5, 1.0e3, 1000)

    cmp = compare_power_spectrum_linear(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, k_test, z
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)


@pytest.mark.parametrize("z", Z_ARRAY, ids=Z_IDS)
def test_power_spectrum_nonlinear(
    ccl_cosmo_eh_halofit: pyccl.Cosmology, nc_cosmo_eh_halofit: ncpy.Cosmology, z
) -> None:
    """Compare NumCosmo and CCL nonlinear power spectrum."""
    if ccl_cosmo_eh_halofit.high_precision:
        rtol = 1.0e-19
    else:
        rtol = 1.0e-4
        if z > 3.0:
            rtol = 1.0e-3
        if z > 4.0:
            rtol = 1.0e-2

    k_test = np.geomspace(5.0e-5, 1.0e3, 1000)

    cmp = compare_power_spectrum_nonlinear(
        ccl_cosmo_eh_halofit, nc_cosmo_eh_halofit, k_test, z
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)


@pytest.mark.parametrize("z", Z_ARRAY, ids=Z_IDS)
def test_sigma_r(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology, z
) -> None:
    """Compare NumCosmo and CCL sigma_r."""
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-7
    else:
        rtol = 1.0e-5

    r_test = np.geomspace(5.0e-2, 1.0e2, 1000)

    cmp = compare_sigma_r(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, r_test, z)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)
