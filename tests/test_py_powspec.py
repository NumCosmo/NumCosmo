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

import pytest

import numpy as np
from numpy.testing import assert_allclose

import pyccl

import numcosmo_py.cosmology as ncpy
from numcosmo_py import Ncm, Nc
from numcosmo_py.ccl.nc_ccl import create_nc_obj

from .fixtures_ccl import (  # pylint: disable=unused-import # noqa: F401
    fixture_k_a,
    fixture_z_a,
    fixture_z_high_a,
    fixture_ccl_cosmo_eh_linear,
    fixture_ccl_cosmo_eh_halofit,
    fixture_nc_cosmo_eh_linear,
    fixture_nc_cosmo_eh_halofit,
)

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

    k_vec = Ncm.Vector.new_array(k_a.tolist())
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
    else:
        reltol_target = 1.0e-3
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
    else:
        reltol_target = 1.0e-2
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

    k_vec = Ncm.Vector.new_array(k_a.tolist())
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


def test_powspec_class(nc_cosmo_eh_linear: ncpy.Cosmology) -> None:
    """Test the power spectrum class."""
    ps_ml = Nc.PowspecMLCBE.new()
    cbe = ps_ml.peek_cbe()
    cbe.use_ppf(True)

    ps_ml.prepare(nc_cosmo_eh_linear.cosmo)

    ik_min = ps_ml.get_intern_k_min()
    ik_max = ps_ml.get_intern_k_max()
    kmin = ps_ml.get_kmin()
    kmax = ps_ml.get_kmax()

    assert ik_min > 0
    assert ik_max > ik_min

    k_a = np.geomspace(kmin, kmax, 100)

    ps_a = np.array([ps_ml.eval(nc_cosmo_eh_linear.cosmo, 0.0, k) for k in k_a])

    assert np.all(np.isfinite(ps_a))


def test_powspec_class_deriv_z(nc_cosmo_eh_linear: ncpy.Cosmology) -> None:
    """Test the power spectrum class."""
    ps_ml = Nc.PowspecMLCBE.new()
    cbe = ps_ml.peek_cbe()
    cbe.use_ppf(True)

    ps_ml.prepare(nc_cosmo_eh_linear.cosmo)

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
        dps, _ = diff.rc_d1_1_to_1(
            0.2, lambda z, k0: ps_ml.eval(nc_cosmo_eh_linear.cosmo, z, k0), k
        )
        assert_allclose(
            dps, ps_ml.deriv_z(nc_cosmo_eh_linear.cosmo, 0.2, k), atol=0.0, rtol=1e-11
        )
