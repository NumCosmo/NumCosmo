#!/usr/bin/env python
#
# fixtures_ccl.py
#
# Thu Aug 01 11:44:12 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# fixtures_ccl.py
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

"""Fixtures for CCL tests."""

from itertools import product
import pytest

import numpy as np

pytest.importorskip("pyccl")
# flake8: noqa: E402
# pylint: disable=wrong-import-position
import pyccl

from numcosmo_py import Ncm, Nc
import numcosmo_py.cosmology as ncpy
from numcosmo_py.ccl.nc_ccl import CCLParams, create_nc_obj


@pytest.fixture(name="k_a")
def fixture_k_a():
    """Fixture for k array."""
    return np.geomspace(1.0e-6, 1.0e2, 1000)


@pytest.fixture(name="z_a")
def fixture_z_a():
    """Fixture for z array."""
    return np.linspace(0.0, 5.0, 1000)


@pytest.fixture(name="z_high_a")
def fixture_z_high_a():
    """Fixture for z array with high redshifts."""
    return np.geomspace(0.01, 1100.0, 2000)


@pytest.fixture(
    name="ccl_cosmo_eh_linear",
    params=product([False, True], range(3)),
    ids=lambda x: f"high_prec={x[0]},index={x[1]}",
    scope="module",
)
def fixture_ccl_cosmo_eh_linear(request) -> pyccl.Cosmology:
    """Fixture for CCL Cosmology."""
    h, Omega_c, Omega_b, Omega_k, Neff = 0.7, 0.25, 0.05, 0.0, 3.046
    sigma8, n_s = 0.9, 0.96

    Omega_v_vals = np.array([0.7, 0.65, 0.75])
    w0_vals = np.array([-1.0, -0.9, -1.1])
    wa_vals = np.array([0.0, 0.1, -0.1])
    m_nu = [[0.0, 0.0, 0.0], [0.04, 0.0, 0.0], [0.02, 0.02, 0.0]]

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
        m_nu=m_nu[index],
        transfer_function="eisenstein_hu",
        matter_power_spectrum="linear",
    )

    ccl_cosmo.high_precision = high_prec
    return ccl_cosmo


# CCL Halofit is too slow for high precision
@pytest.fixture(
    name="ccl_cosmo_eh_halofit",
    params=product([False], range(3)),
    ids=lambda x: f"high_prec={x[0]}, index={x[1]}",
    scope="module",
)
def fixture_ccl_cosmo_eh_halofit(request) -> pyccl.Cosmology:
    """Fixture for CCL Cosmology."""
    Omega_c_vals = [0.25, 0.24, 0.26]
    Omega_b = 0.05
    Omega_k = 0.0
    h = 0.7
    n_s = 0.96
    Neff = 3.046
    sigma8_vals = [0.9, 0.7, 0.6]

    Omega_v = 0.7
    w0 = -1.0
    wa = 0.0
    m_nu = [[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.2, 0.2, 0.0]]

    high_prec, index = request.param

    Omega_k = 1.0 - Omega_c_vals[index] - Omega_b - Omega_v

    if high_prec:
        CCLParams.set_high_prec_params()
    else:
        CCLParams.set_default_params()

    ccl_cosmo = pyccl.Cosmology(
        Omega_c=Omega_c_vals[index],
        Omega_b=Omega_b,
        Neff=Neff,
        h=h,
        sigma8=sigma8_vals[index],
        n_s=n_s,
        Omega_k=Omega_k,
        w0=w0,
        wa=wa,
        m_nu=m_nu[index],
        transfer_function="eisenstein_hu",
        matter_power_spectrum="halofit",
    )

    ccl_cosmo.high_precision = high_prec
    return ccl_cosmo


@pytest.fixture(name="nc_cosmo_eh_linear", scope="module")
def fixture_nc_cosmo_eh_linear(ccl_cosmo_eh_linear) -> ncpy.Cosmology:
    """Fixture for CCL and NumCosmo Cosmology."""
    return create_nc_obj(ccl_cosmo_eh_linear, dist_z_max=2000.0)


@pytest.fixture(name="nc_cosmo_eh_halofit", scope="module")
def fixture_nc_cosmo_eh_halofit(ccl_cosmo_eh_halofit) -> ncpy.Cosmology:
    """Fixture for CCL and NumCosmo Cosmology."""
    return create_nc_obj(ccl_cosmo_eh_halofit, dist_z_max=2000.0)


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


@pytest.fixture(name="ccl_cmb_isw")
def fixture_ccl_cmb_isw(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
) -> pyccl.CMBLensingTracer:
    """Fixture for CCL CMB ISW tracer."""
    ccl_cmb_isw = pyccl.ISWTracer(ccl_cosmo_eh_linear, n_chi=10_000)
    return ccl_cmb_isw


@pytest.fixture(name="ccl_tsz")
def fixture_ccl_tsz(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
) -> pyccl.tSZTracer:
    """Fixture for CCL tSZ tracer."""
    ccl_tsz = pyccl.tSZTracer(ccl_cosmo_eh_linear, z_max=6.0, n_chi=10_000)
    return ccl_tsz


@pytest.fixture(name="ccl_gal")
def fixture_ccl_gal(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_gal: Nc.XcorKernelGal,
) -> pyccl.NumberCountsTracer:
    """Fixture for CCL galaxy tracer."""
    dndz: Ncm.Spline = nc_gal.props.dndz
    assert isinstance(dndz, Ncm.Spline)

    bias = nc_gal.orig_vparam_get(Nc.XcorKernelGalVParams.BIAS, 0)
    magbias = nc_gal.orig_param_get(Nc.XcorKernelGalSParams.MAG_BIAS)
    z_array = np.array(dndz.peek_xv().dup_array())
    ccl_gal = pyccl.NumberCountsTracer(
        ccl_cosmo_eh_linear,
        has_rsd=False,
        dndz=(
            z_array,
            np.array(dndz.peek_yv().dup_array()),
        ),
        bias=(z_array, np.ones_like(z_array) * bias) if bias != 0.0 else None,
        mag_bias=(z_array, np.ones_like(z_array) * magbias) if magbias != 0.0 else None,
        n_samples=len(z_array),
    )
    return ccl_gal


@pytest.fixture(name="ccl_weak_lensing")
def fixture_ccl_weak_lensing(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_weak_lensing: Nc.XcorKernelWeakLensing,
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
