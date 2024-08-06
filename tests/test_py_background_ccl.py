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

from numpy.testing import assert_allclose

import numpy as np
import pyccl

from numcosmo_py import Ncm
from numcosmo_py.ccl.nc_ccl import create_nc_obj

from .ccl_fixtures import (  # pylint: disable=unused-import # noqa: F401
    fixture_k_a,
    fixture_z_a,
    fixture_ccl_cosmo_eh_linear,
    fixture_ccl_cosmo_eh_halofit,
)

Ncm.cfg_init()


def test_background_params(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    cosmo, _, _, _, _ = create_nc_obj(ccl_cosmo_eh_linear)

    reltol = 1.0e-15

    assert_allclose(cosmo.H0(), ccl_cosmo_eh_linear["h"] * 100, rtol=reltol)
    assert_allclose(cosmo.Omega_b0(), ccl_cosmo_eh_linear["Omega_b"], rtol=reltol)
    assert_allclose(cosmo.Omega_c0(), ccl_cosmo_eh_linear["Omega_c"], rtol=reltol)
    assert_allclose(cosmo.Omega_k0(), ccl_cosmo_eh_linear["Omega_k"], rtol=reltol)
    assert_allclose(cosmo.T_gamma0(), ccl_cosmo_eh_linear["T_CMB"], rtol=reltol)
    assert_allclose(cosmo.Omega_g0(), ccl_cosmo_eh_linear["Omega_g"], rtol=reltol)
    assert_allclose(cosmo.Neff(), ccl_cosmo_eh_linear["Neff"], rtol=reltol)
    assert_allclose(
        cosmo.param_get_by_name("w0"), ccl_cosmo_eh_linear["w0"], rtol=reltol
    )
    assert_allclose(
        cosmo.param_get_by_name("w1"), ccl_cosmo_eh_linear["wa"], rtol=reltol
    )


def test_background_functions(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL background functions."""
    cosmo, _, _, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-15
    else:
        rtol = 1.0e-15

    z_test = np.linspace(0.0, 5.0, 2000)[1:]
    a_test = 1.0 / (1.0 + z_test)

    E2 = np.array([cosmo.E2(z) for z in z_test])

    assert_allclose(
        np.array([cosmo.E2Omega_m(z) for z in z_test]) / E2,
        ccl_cosmo_eh_linear.omega_x(a_test, "matter"),
        rtol=rtol,
    )

    assert_allclose(
        np.array([cosmo.E2Omega_g(z) for z in z_test]) / E2,
        ccl_cosmo_eh_linear.omega_x(a_test, "radiation"),
        rtol=rtol,
    )

    assert_allclose(
        np.array([cosmo.E2Omega_k(z) for z in z_test]) / E2,
        ccl_cosmo_eh_linear.omega_x(a_test, "curvature"),
        rtol=rtol,
    )

    assert_allclose(
        np.array([cosmo.E2Omega_nu(z) for z in z_test]) / E2,
        ccl_cosmo_eh_linear.omega_x(a_test, "neutrinos_rel"),
        rtol=rtol,
    )

    assert_allclose(
        np.array([cosmo.E2Omega_de(z) for z in z_test]) / E2,
        ccl_cosmo_eh_linear.omega_x(a_test, "dark_energy"),
        rtol=rtol * 2.0,
    )


def test_background_h_over_h0(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL background functions."""
    cosmo, _, _, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-9
    else:
        rtol = 1.0e-5

    z_test = np.linspace(0.0, 5.0, 2000)[1:]
    a_test = 1.0 / (1.0 + z_test)

    E2 = np.array([cosmo.E2(z) for z in z_test])

    assert_allclose(E2, ccl_cosmo_eh_linear.h_over_h0(a_test) ** 2, rtol=rtol)


def test_background_distances(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL distances."""
    cosmo, dist, _, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    dist.prepare(cosmo)
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-10
    else:
        rtol = 1.0e-6

    z_test = np.linspace(0.0, 5.0, 2000)[1:]
    a_test = 1.0 / (1.0 + z_test)
    RH_Mpc = cosmo.RH_Mpc()

    assert_allclose(
        [dist.comoving(cosmo, z) * RH_Mpc for z in z_test],
        ccl_cosmo_eh_linear.comoving_radial_distance(a_test),
        rtol=rtol,
    )

    assert_allclose(
        [dist.angular_diameter(cosmo, z) * RH_Mpc for z in z_test],
        ccl_cosmo_eh_linear.angular_diameter_distance(a_test),
        rtol=rtol,
    )

    assert_allclose(
        [dist.luminosity(cosmo, z) * RH_Mpc for z in z_test],
        ccl_cosmo_eh_linear.luminosity_distance(a_test),
        rtol=rtol,
    )

    assert_allclose(
        [dist.dmodulus(cosmo, z) + +5 * np.log10(RH_Mpc) for z in z_test],
        ccl_cosmo_eh_linear.distance_modulus(a_test),
        rtol=rtol,
    )

    assert_allclose(
        [dist.comoving_volume_element(cosmo, z) * RH_Mpc**3 for z in z_test],
        a_test**2 * ccl_cosmo_eh_linear.comoving_volume_element(a_test),
        rtol=rtol * 10.0,
    )


def test_background_scale_factor_of_chi(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL scale factor of comoving distance."""
    cosmo, dist, _, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    dist.compute_inv_comoving(True)
    dist.prepare(cosmo)
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-9
    else:
        rtol = 1.0e-6

    RH_Mpc = cosmo.RH_Mpc()

    c_max = dist.comoving(cosmo, 15.0)
    comoving = np.geomspace(1.0e-4, c_max, 4000)[1:]
    z_array = [dist.inv_comoving(cosmo, c) for c in comoving]

    assert_allclose(
        1.0 / (1.0 + np.array(z_array)),
        ccl_cosmo_eh_linear.scale_factor_of_chi(comoving * RH_Mpc),
        rtol=rtol,
    )


def test_background_growth_lowz(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL growth factor."""
    cosmo, dist, ps_lin, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    dist.prepare(cosmo)
    if ccl_cosmo_eh_linear.high_precision:
        reltol_growth = 1.0e-10
        reltol_rate = 1.0e-8
    else:
        reltol_growth = 1.0e-6
        reltol_rate = 1.0e-6

    k_pivot = 1.0
    z_test = np.linspace(0.0, 5.0, 2000)[1:]
    a_test = 1.0 / (1.0 + z_test)

    gf = ps_lin.peek_gf()
    gf.prepare(cosmo)

    nc_D_a = np.array([gf.eval(cosmo, z) for z in z_test])
    nc_D2_a = np.array(
        [
            ps_lin.eval(cosmo, z, k_pivot) / ps_lin.eval(cosmo, 0.0, k_pivot)
            for z in z_test
        ]
    )
    nc_dD_a = np.array([gf.eval_deriv(cosmo, z) for z in z_test])
    nc_dD2_a = np.array(
        [
            ps_lin.deriv_z(cosmo, z, k_pivot) / ps_lin.eval(cosmo, 0.0, k_pivot)
            for z in z_test
        ]
    )
    nc_f_a = -(1.0 + z_test) * nc_dD_a / nc_D_a

    assert_allclose(nc_D_a**2, nc_D2_a, rtol=1.0e-15, atol=0.0)
    assert_allclose(nc_dD2_a / (2.0 * nc_D_a), nc_dD_a, rtol=1.0e-15, atol=0.0)

    assert_allclose(
        nc_D_a, ccl_cosmo_eh_linear.growth_factor(a_test), rtol=reltol_growth, atol=0.0
    )

    assert_allclose(
        nc_f_a, ccl_cosmo_eh_linear.growth_rate(a_test), rtol=reltol_rate, atol=0.0
    )


def test_background_growth_highz(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL growth factor."""
    cosmo, dist, ps_lin, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    dist.prepare(cosmo)
    if ccl_cosmo_eh_linear.high_precision:
        reltol_growth = 1.0e-4
        reltol_rate = 1.0e-4
    else:
        reltol_growth = 1.0e-3
        reltol_rate = 1.0e-2

    z_test = np.geomspace(0.001, 1100.0, 4000)[1:]
    a_test = 1.0 / (1.0 + z_test)

    gf = ps_lin.peek_gf()
    gf.prepare(cosmo)

    nc_D_a = np.array([gf.eval(cosmo, z) for z in z_test])
    nc_f_a = (
        -(1.0 + z_test) * np.array([gf.eval_deriv(cosmo, z) for z in z_test]) / nc_D_a
    )
    assert_allclose(
        nc_D_a, ccl_cosmo_eh_linear.growth_factor(a_test), rtol=reltol_growth, atol=0.0
    )

    assert_allclose(
        nc_f_a, ccl_cosmo_eh_linear.growth_rate(a_test), rtol=reltol_rate, atol=0.0
    )
