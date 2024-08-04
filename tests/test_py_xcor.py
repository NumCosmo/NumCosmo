#!/usr/bin/env python
#
# test_py_xcor.py
#
# Thu Aug 01 11:45:10 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor.py
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

import numpy as np
from numpy.testing import assert_allclose

import pyccl

from numcosmo_py import Ncm, Nc
from numcosmo_py.ccl.nc_ccl import create_nc_obj

from .ccl_fixtures import (  # pylint: disable=unused-import # noqa: F401
    fixture_k_a,
    fixture_z_a,
    fixture_ccl_cosmo_eh_linear,
    fixture_ccl_cosmo_eh_halofit,
)

Ncm.cfg_init()


def test_cmb1_auto(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 8.0e-12
    else:
        reltol_target = 8.0e-2

    lmax = 3000
    ells = np.arange(2, lmax + 1)
    z_lss = 1090.0

    cmb1 = pyccl.CMBLensingTracer(ccl_cosmo_eh_linear, z_source=z_lss, n_samples=100000)
    assert cmb1 is not None

    ccl_cl_cmblens_auto1 = pyccl.angular_cl(
        ccl_cosmo_eh_linear, cmb1, cmb1, ells, l_limber=-1
    )

    assert ccl_cl_cmblens_auto1 is not None
    assert all(np.isfinite(ccl_cl_cmblens_auto1))
    assert all(ccl_cl_cmblens_auto1 >= 0.0)

    cosmo, dist, ps_lin, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    recomb = Nc.RecombSeager()

    Nl_lensing = Ncm.Vector.new_array(np.zeros(lmax + 1).tolist())
    nc_cmblens_1 = Nc.XcorLimberKernelCMBLensing.new(dist, recomb, Nl_lensing)
    mset = Ncm.MSet.empty_new()
    mset.set(cosmo)
    mset.push(nc_cmblens_1)
    xcor = Nc.Xcor.new(dist, ps_lin, Nc.XcorLimberMethod.GSL)
    vp_cmblens_auto1 = Ncm.Vector.new(lmax + 1 - 2)
    xcor.prepare(cosmo)
    nc_cmblens_1.prepare(cosmo)
    xcor.limber(nc_cmblens_1, nc_cmblens_1, cosmo, 2, lmax, vp_cmblens_auto1)

    assert_allclose(
        ccl_cl_cmblens_auto1, vp_cmblens_auto1.dup_array(), rtol=reltol_target, atol=0.0
    )
