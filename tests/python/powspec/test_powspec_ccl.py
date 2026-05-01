#!/usr/bin/env python
#
# test_powspec_ccl.py
#
# Wed Feb 14 13:44:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_powspec_ccl.py
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

"""Unit tests for NumCosmo power-spectra with CCL comparison."""

import pytest

import numpy as np
from numpy.testing import assert_allclose

pytest.importorskip("getdist")
pytest.importorskip("pyccl")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

import pyccl

pytestmark = [pytest.mark.ccl, pytest.mark.powspec]

import numcosmo_py.cosmology as ncpy
from numcosmo_py import Ncm
from numcosmo_py.helper import npa_to_seq
from numcosmo_py.ccl.nc_ccl import create_nc_obj
from numcosmo_py.ccl.comparison import (
    compare_power_spectrum_linear,
    compare_power_spectrum_nonlinear,
    compare_sigma_r,
)

pytest_plugins = ["python.fixtures_ccl"]

Ncm.cfg_init()


def _test_powspec_transfer_any(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    k_a: np.ndarray,
    z_a: np.ndarray,
    reltol_target: float,
) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    cosmo = nc_cosmo_eh_linear.cosmo
    ps_ml = nc_cosmo_eh_linear.ps_ml

    k_vec = Ncm.Vector.new_array(npa_to_seq(k_a))
    Pk_vec = Ncm.Vector.new(k_vec.len())

    for z in z_a:
        a_i = 1.0 / (1.0 + z)
        pk_ccl = pyccl.linear_matter_power(ccl_cosmo_eh_linear, k_a, a_i)

        ps_ml.eval_vec(cosmo, z, k_vec, Pk_vec)
        pk_nc = Pk_vec.dup_array()

        assert_allclose(pk_nc, pk_ccl, rtol=reltol_target)


def test_powspec_transfer_lowz(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    k_a: np.ndarray,
    z_a: np.ndarray,
) -> None:
    """Compare NumCosmo and CCL transfer functions for low redshifts."""
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-7
        if ccl_cosmo_eh_linear["m_nu"] != 0.0:
            reltol_target = 2.0e-3
    else:
        reltol_target = 1.0e-3
        if ccl_cosmo_eh_linear["m_nu"] != 0.0:
            reltol_target = 2.0e-3
    _test_powspec_transfer_any(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, k_a, z_a, reltol_target
    )


def test_powspec_transfer_highz(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    k_a: np.ndarray,
    z_high_a: np.ndarray,
) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-4
        if ccl_cosmo_eh_linear["m_nu"] != 0.0:
            reltol_target = 4.0e-1
    else:
        reltol_target = 1.0e-2
        if ccl_cosmo_eh_linear["m_nu"] != 0.0:
            reltol_target = 2.0e-1
    _test_powspec_transfer_any(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, k_a, z_high_a, reltol_target
    )


def test_powspec_halofit(
    ccl_cosmo_eh_halofit: pyccl.Cosmology,
    nc_cosmo_eh_halofit: ncpy.Cosmology,
    k_a: np.ndarray,
    z_a: np.ndarray,
) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    cosmo = nc_cosmo_eh_halofit.cosmo
    ps_mnl = nc_cosmo_eh_halofit.ps_mnl
    if ccl_cosmo_eh_halofit.high_precision:
        reltol_target: float = 6.0e-3
    else:
        reltol_target = 5.0e-2

    k_vec = Ncm.Vector.new_array(npa_to_seq(k_a))
    Pk_vec = Ncm.Vector.new(k_vec.len())

    for z in z_a:
        a_i = 1.0 / (1.0 + z)
        ccl_cosmo_eh_halofit.compute_nonlin_power()
        pk_ccl = pyccl.nonlin_matter_power(ccl_cosmo_eh_halofit, k_a, a_i)

        ps_mnl.eval_vec(cosmo, z, k_vec, Pk_vec)
        pk_nc = Pk_vec.dup_array()

        assert_allclose(pk_nc, pk_ccl, rtol=reltol_target)


@pytest.mark.parametrize("sigma8", [0.01, 0.05, 0.1, 0.2, 0.5])
def test_powspec_halofit_linear_universe(
    sigma8, k_a: np.ndarray, z_a: np.ndarray
) -> None:
    """Test NumCosmo for an linear universe (very small sigma8)."""
    # Linear universe, halofit power spectrum
    ccl_cosmo = pyccl.Cosmology(
        Omega_c=0.25,
        Omega_b=0.05,
        Neff=3.046,
        h=0.7,
        sigma8=sigma8,
        n_s=0.96,
        Omega_k=0.0,
        w0=-1.0,
        wa=0.0,
        transfer_function="eisenstein_hu",
        matter_power_spectrum="halofit",
    )

    nc_cosmo = create_nc_obj(ccl_cosmo, ps_nln_z_max=1.0)
    cosmo = nc_cosmo.cosmo
    ps_mnl = nc_cosmo.ps_mnl

    k_vec = Ncm.Vector.new_array(npa_to_seq(k_a))
    Pk_vec = Ncm.Vector.new(k_vec.len())

    step = 10
    for z in z_a[::step]:
        ps_mnl.eval_vec(cosmo, z, k_vec, Pk_vec)
        pk_nc = Pk_vec.dup_array()
        assert all(np.isfinite(pk_nc))

        for i, k in enumerate(k_a[::step]):
            Pk = ps_mnl.eval(cosmo, z, k)
            assert_allclose(Pk_vec.get(i * step), Pk)


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
