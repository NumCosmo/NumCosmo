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

from itertools import product
import time
import pytest

import numpy as np
from numpy.testing import assert_allclose

import pyccl

from numcosmo_py import Ncm, Nc
from numcosmo_py.ccl.nc_ccl import create_nc_obj, CCLParams

from .ccl_fixtures import (  # pylint: disable=unused-import # noqa: F401
    fixture_k_a,
    fixture_z_a,
    fixture_ccl_cosmo_eh_linear,
    fixture_ccl_cosmo_eh_halofit,
)

Ncm.cfg_init()


def test_cmb1_auto(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL transfer functions."""

    cosmo, dist, _, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    dist.prepare(cosmo)

    assert_allclose(cosmo.H0(), ccl_cosmo_eh_linear["h"] * 100, rtol=1e-15)
    assert_allclose(cosmo.Omega_b0(), ccl_cosmo_eh_linear["Omega_b"], rtol=1e-15)
    assert_allclose(cosmo.Omega_c0(), ccl_cosmo_eh_linear["Omega_c"], rtol=1e-15)
    assert_allclose(cosmo.Omega_k0(), ccl_cosmo_eh_linear["Omega_k"], rtol=1e-14)
    assert_allclose(
        cosmo.param_get_by_name("w0"), ccl_cosmo_eh_linear["w0"], rtol=1e-15
    )
    assert_allclose(
        cosmo.param_get_by_name("w1"), ccl_cosmo_eh_linear["wa"], rtol=1e-15
    )
