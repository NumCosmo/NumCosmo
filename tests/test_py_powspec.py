#!/usr/bin/env python
#
# test_py_powspec.py
#
# Wed Feb 14 13:44:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
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
import time
import pytest

import numpy as np
from numpy.testing import assert_allclose

import pyccl

from numcosmo_py import Ncm, Nc
from numcosmo_py.ccl.nc_ccl import create_nc_obj, CCLParams

Ncm.cfg_init()


@pytest.fixture(name="k_a")
def k_a_fixture():
    """Fixture for k array."""
    return np.geomspace(1.0e-5, 1.0e1, 1000)


@pytest.fixture(name="z_a")
def z_a_fixture():
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


def test_powspec_transfer(
    ccl_cosmo_eh_linear: pyccl.Cosmology, k_a: np.ndarray, z_a: np.ndarray
) -> None:
    """Compare NumCosmo and CCL transfer functions."""

    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-7
    else:
        reltol_target = 1.0e-3

    nc_elapse = 0.0
    ccl_elapse = 0.0

    t0 = time.time()
    cosmo: Nc.HICosmo
    ps_lin: Nc.PowspecML
    cosmo, _, ps_lin, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    t1 = time.time()

    nc_elapse += t1 - t0

    k_vec = Ncm.Vector.new_array(k_a.tolist())
    Pk_vec = Ncm.Vector.new(k_vec.len())

    for z in z_a:
        a_i = 1.0 / (1.0 + z)
        t0 = time.time()
        pk_ccl = pyccl.linear_matter_power(ccl_cosmo_eh_linear, k_a, a_i)
        t1 = time.time()

        ps_lin.eval_vec(cosmo, z, k_vec, Pk_vec)
        pk_nc = Pk_vec.dup_array()
        t2 = time.time()

        nc_elapse += t2 - t1
        ccl_elapse += t1 - t0

        assert_allclose(pk_nc, pk_ccl, rtol=reltol_target)

    print(f"# CCL time: {ccl_elapse:.3f} s, NumCosmo time: {nc_elapse:.3f} s")


def test_powspec_halofit(
    ccl_cosmo_eh_halofit: pyccl.Cosmology, k_a: np.ndarray, z_a: np.ndarray
) -> None:
    """Compare NumCosmo and CCL transfer functions."""

    cosmo, _, _, ps_mnl, _ = create_nc_obj(ccl_cosmo_eh_halofit)

    if ccl_cosmo_eh_halofit.high_precision:
        reltol_target: float = 6.0e-3
    else:
        reltol_target = 5.0e-2

    nc_elapse = 0.0
    ccl_elapse = 0.0

    t0 = time.time()
    ps_mnl.prepare(cosmo)
    t1 = time.time()

    nc_elapse += t1 - t0

    k_vec = Ncm.Vector.new_array(k_a.tolist())
    Pk_vec = Ncm.Vector.new(k_vec.len())

    for z in z_a:
        a_i = 1.0 / (1.0 + z)
        t0 = time.time()
        ccl_cosmo_eh_halofit.compute_nonlin_power()
        pk_ccl = pyccl.nonlin_matter_power(ccl_cosmo_eh_halofit, k_a, a_i)
        t1 = time.time()

        ps_mnl.eval_vec(cosmo, z, k_vec, Pk_vec)
        pk_nc = Pk_vec.dup_array()
        t2 = time.time()

        nc_elapse += t2 - t1
        ccl_elapse += t1 - t0
        assert_allclose(pk_nc, pk_ccl, rtol=reltol_target)

    print(f"# CCL time: {ccl_elapse:.3f} s, NumCosmo time: {nc_elapse:.3f} s")


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

    cosmo, _, _, ps_mnl, _ = create_nc_obj(ccl_cosmo, ps_nln_z_max=1.0)

    k_vec = Ncm.Vector.new_array(k_a.tolist())
    Pk_vec = Ncm.Vector.new(k_vec.len())

    step = 10
    for z in z_a[::step]:
        ps_mnl.eval_vec(cosmo, z, k_vec, Pk_vec)
        pk_nc = Pk_vec.dup_array()
        assert all(np.isfinite(pk_nc))

        for i, k in enumerate(k_a[::step]):
            Pk = ps_mnl.eval(cosmo, z, k)
            assert_allclose(Pk_vec.get(i * step), Pk)
