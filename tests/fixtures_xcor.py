#!/usr/bin/env python
#
# fixtures_xcor.py
#
# Thu Aug 01 11:44:12 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# fixtures_xcor.py
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

"""Fixtures for XCor tests."""

import pytest

import numpy as np

import pyccl

import numcosmo_py.cosmology as ncpy
from numcosmo_py import Ncm, Nc
from .fixtures_ccl import (  # pylint: disable=unused-import # noqa: F401
    fixture_ccl_cosmo_eh_linear,
    fixture_nc_cosmo_eh_linear,
)

Ncm.cfg_init()


@pytest.fixture(name="ccl_cmb_lens")
def fixture_ccl_cmb_lens(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
) -> pyccl.CMBLensingTracer:
    """Fixture for CCL CMB lensing tracer."""
    z_lss = nc_cosmo_eh_linear.dist.decoupling_redshift(nc_cosmo_eh_linear.cosmo)

    ccl_cmb_lens = pyccl.CMBLensingTracer(
        ccl_cosmo_eh_linear, z_source=z_lss, n_samples=10_000
    )

    return ccl_cmb_lens


@pytest.fixture(name="nc_cmb_lens")
def fixture_nc_cmb_lens(
    nc_cosmo_eh_linear: ncpy.Cosmology,
) -> Nc.XcorLimberKernelCMBLensing:
    """Fixture for NumCosmo CMB lensing tracer."""
    lmax = 3000
    nc_cmb_lens = Nc.XcorLimberKernelCMBLensing.new(
        nc_cosmo_eh_linear.dist,
        Nc.RecombSeager(),
        Ncm.Vector.new_array(np.arange(lmax + 1).tolist()),
    )

    return nc_cmb_lens


@pytest.fixture(name="ccl_cmb_isw")
def fixture_ccl_cmb_isw(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
) -> pyccl.CMBLensingTracer:
    """Fixture for CCL CMB ISW tracer."""
    ccl_cmb_isw = pyccl.ISWTracer(ccl_cosmo_eh_linear, n_chi=10_000)
    return ccl_cmb_isw


@pytest.fixture(name="nc_cmb_isw")
def fixture_nc_cmb_isw(
    nc_cosmo_eh_linear: ncpy.Cosmology,
) -> Nc.XcorLimberKernelCMBISW:
    """Fixture for NumCosmo CMB ISW tracer."""
    lmax = 3000
    nc_cmb_isw = Nc.XcorLimberKernelCMBISW.new(
        nc_cosmo_eh_linear.dist,
        nc_cosmo_eh_linear.ps_ml,
        nc_cosmo_eh_linear.recomb,
        Ncm.Vector.new_array(np.arange(lmax + 1).tolist()),
    )

    return nc_cmb_isw


GAL_Z_CENTERS = np.linspace(0.3, 1.2, 4)
GAL_Z_SIGMA = 0.02
GAL_Z_IDS = [f"z={z:.2f}" for z in GAL_Z_CENTERS]

SRC_GAL_Z_CENTERS = np.linspace(0.5, 1.6, 4)
SRC_GAL_Z_SIGMA = 0.02
SRC_GAL_Z_IDS = [f"z={z:.2f}" for z in SRC_GAL_Z_CENTERS]


@pytest.fixture(name="nc_gal", params=GAL_Z_CENTERS, ids=GAL_Z_IDS)
def fixture_nc_gal(
    nc_cosmo_eh_linear: ncpy.Cosmology,
    request,
) -> Nc.XcorLimberKernelGal:
    """Fixture for NumCosmo galaxy tracer."""
    mu = request.param
    sigma = GAL_Z_SIGMA
    z_a = np.linspace(0.0, 2.0, 2000)
    nz_a = np.exp(-((z_a - mu) ** 2) / sigma**2 / 2.0) / np.sqrt(2.0 * np.pi * sigma**2)

    z_v = Ncm.Vector.new_array(z_a.tolist())
    nz_v = Ncm.Vector.new_array(nz_a.tolist())
    dndz = Ncm.SplineCubicNotaknot.new_full(z_v, nz_v, True)

    nc_gal = Nc.XcorLimberKernelGal.new(
        0.2, 2.0, 100, 1.234, dndz, nc_cosmo_eh_linear.dist, False
    )

    return nc_gal


@pytest.fixture(name="ccl_gal")
def fixture_ccl_gal(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_gal: Nc.XcorLimberKernelGal,
) -> pyccl.NumberCountsTracer:
    """Fixture for CCL galaxy tracer."""
    dndz: Ncm.Spline = nc_gal.props.dndz
    assert isinstance(dndz, Ncm.Spline)

    ccl_gal = pyccl.NumberCountsTracer(
        ccl_cosmo_eh_linear,
        has_rsd=False,
        dndz=(
            np.array(dndz.peek_xv().dup_array()),
            np.array(dndz.peek_yv().dup_array()),
        ),
    )
    return ccl_gal


@pytest.fixture(name="nc_weak_lensing", params=SRC_GAL_Z_CENTERS, ids=SRC_GAL_Z_IDS)
def fixture_nc_weak_lensing(
    nc_cosmo_eh_linear: ncpy.Cosmology,
    request,
) -> Nc.XcorLimberKernelWeakLensing:
    """Fixture for NumCosmo weak lensing tracer."""
    mu = request.param
    sigma = SRC_GAL_Z_SIGMA
    z_a = np.linspace(0.0, 2.0, 20_000)
    nz_a = np.exp(-((z_a - mu) ** 2) / sigma**2 / 2.0) / np.sqrt(2.0 * np.pi * sigma**2)

    z_v = Ncm.Vector.new_array(z_a.tolist())
    nz_v = Ncm.Vector.new_array(nz_a.tolist())
    dndz = Ncm.SplineCubicNotaknot.new_full(z_v, nz_v, True)

    nc_wl = Nc.XcorLimberKernelWeakLensing.new(
        0.0, 2.0, dndz, 3.0, 7.0, nc_cosmo_eh_linear.dist
    )

    return nc_wl


@pytest.fixture(name="ccl_weak_lensing")
def fixture_ccl_weak_lensing(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_weak_lensing: Nc.XcorLimberKernelWeakLensing,
) -> pyccl.WeakLensingTracer:
    """Fixture for CCL weak lensing tracer."""
    dndz: Ncm.Spline = nc_weak_lensing.props.dndz
    assert isinstance(dndz, Ncm.Spline)

    ccl_wl = pyccl.WeakLensingTracer(
        ccl_cosmo_eh_linear,
        dndz=(
            np.array(dndz.peek_xv().dup_array()),
            np.array(dndz.peek_yv().dup_array()),
        ),
        n_samples=10_000,
    )
    return ccl_wl