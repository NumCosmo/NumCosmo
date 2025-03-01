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
import matplotlib.pyplot as plt
import pyccl

import numcosmo_py.cosmology as ncpy
from numcosmo_py import Ncm
from numcosmo_py.ccl.comparison import (
    compare_Hubble,
    compare_Omega_m,
    compare_Omega_g,
    compare_Omega_k,
    compare_Omega_nu,
    compare_Omega_mnu,
    compare_Omega_de,
    compare_distance_comoving,
    compare_distance_transverse,
    compare_distance_angular_diameter,
    compare_distance_luminosity,
    compare_distance_lookback_time,
    compare_distance_modulus,
    compare_distance_comoving_volume,
    compare_scale_factor,
    compare_growth_factor,
    compare_growth_rate,
    compare_Sigma_crit,
    CompareFunc1d,
)
from numcosmo_py.plotting.tools import latex_float
from .fixtures_ccl import (  # pylint: disable=unused-import # noqa: F401
    fixture_k_a,
    fixture_z_a,
    fixture_ccl_cosmo_eh_linear,
    fixture_ccl_cosmo_eh_halofit,
    fixture_nc_cosmo_eh_linear,
    fixture_nc_cosmo_eh_halofit,
)

Ncm.cfg_init()


def test_compare_1d() -> None:
    """Test CompareFunc1d class."""
    x = np.linspace(0.0, 2.0, 100, dtype=np.float64)
    y1 = 0.4 + 0.6 * x**2
    y2 = y1 + np.random.normal(0.0, 1.0e-6, y1.shape)
    compare = CompareFunc1d(x, y1, y2)
    assert compare.rel_diff_mean < 1.0e-5
    assert compare.rel_diff_median < 1.0e-5
    assert compare.rel_diff_max < 1.0e-5
    assert compare.rel_diff_min < 1.0e-5
    assert compare.rel_diff_std < 1.0e-5
    assert np.all(compare.abs_diff < 1.0e-5)
    assert compare.model == "unnamed"
    assert compare.name1 == "1"
    assert compare.name2 == "2"
    assert compare.x_symbol == "x"
    assert compare.y_symbol == "y"
    assert compare.x_unit is None
    assert compare.y_unit is None
    assert compare.x_label == "$x$"
    assert compare.y_label == "$y$"
    assert compare.summary_row() == [
        "unnamed",
        "$y$",
        f"${latex_float(compare.rel_diff_min)}$",
        f"${latex_float(compare.rel_diff_max)}$",
        f"${latex_float(compare.rel_diff_mean)}$",
        f"${latex_float(compare.rel_diff_std)}$",
    ]
    assert CompareFunc1d.table_header() == [
        "Model",
        "Quantity",
        "$\\Delta P_{\\text{min}}/P$",
        "$\\Delta P_{\\text{max}}/P$",
        "$\\overline{\\Delta P / P}$",
        "$\\sigma_{\\Delta P / P}$",
    ]


def test_compare_1d_full() -> None:
    """Test CompareFunc1d class."""
    model = "Bob"
    name1 = "Alice"
    name2 = "Fred"
    x_symbol = "x^x"
    y_symbol = "y^y"
    x_unit = "m"
    y_unit = "kg"
    x = np.linspace(0.0, 2.0, 100, dtype=np.float64)
    y1 = 0.4 + 0.6 * x**2
    y2 = y1 + np.random.normal(0.0, 1.0e-6, y1.shape)
    compare = CompareFunc1d(
        x,
        y1,
        y2,
        model=model,
        name1=name1,
        name2=name2,
        x_symbol=x_symbol,
        y_symbol=y_symbol,
        x_unit=x_unit,
        y_unit=y_unit,
    )
    assert compare.model == model
    assert compare.name1 == name1
    assert compare.name2 == name2
    assert compare.x_symbol == x_symbol
    assert compare.y_symbol == y_symbol
    assert compare.x_unit == x_unit
    assert compare.y_unit == y_unit

    assert compare.x_label == f"${x_symbol}$ [{x_unit}]"
    assert compare.y_label == f"${y_symbol}$ [{y_unit}]"

    for use_g in [True, False]:
        for precision in [1, 2, 3, 4, 5]:
            rel_diff_min_str = latex_float(
                compare.rel_diff_min, convert_g=use_g, precision=precision
            )
            rel_diff_max_str = latex_float(
                compare.rel_diff_max, convert_g=use_g, precision=precision
            )
            rel_diff_mean_str = latex_float(
                compare.rel_diff_mean, convert_g=use_g, precision=precision
            )
            rel_diff_std_str = latex_float(
                compare.rel_diff_std, convert_g=use_g, precision=precision
            )
            assert compare.summary_row(convert_g=use_g, precision=precision) == [
                model,
                f"${y_symbol}$ [{y_unit}]",
                f"${rel_diff_min_str}$",
                f"${rel_diff_max_str}$",
                f"${rel_diff_mean_str}$",
                f"${rel_diff_std_str}$",
            ]


def test_compare_1d_plot() -> None:
    """Test CompareFunc1d class."""
    x = np.linspace(0.0, 2.0, 100, dtype=np.float64)
    y1 = 0.4 + 0.6 * x**2
    y2 = y1 + np.random.normal(0.0, 1.0e-6, y1.shape)
    compare = CompareFunc1d(x, y1, y2)
    _, axs = plt.subplots(1, 2)
    xscale = "log"
    yscale = "linear"
    color = "black"
    lw = 1.2
    compare.plot(axs, xscale=xscale, yscale=yscale, color=color, lw=lw)

    # Check axis scales
    assert axs[0].get_xscale() == xscale
    assert axs[0].get_yscale() == yscale

    assert axs[1].get_xscale() == xscale
    assert axs[1].get_yscale() == "log"  # Always use log scale

    # Check that a line was added
    for i in range(2):
        lines = axs[i].get_lines()
        assert len(lines) > 0  # At least one line must be present

        # Check properties of the first line
        line = lines[0]
        assert line.get_color() == color
        assert line.get_linewidth() == lw


def test_background_params(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology
) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    cosmo = nc_cosmo_eh_linear.cosmo

    reltol = 1.0e-7
    # CCL conversion of Neff creates a small error
    reltol_Neff = 1.0e-8

    assert_allclose(cosmo.H0(), ccl_cosmo_eh_linear["h"] * 100, rtol=reltol)
    assert_allclose(cosmo.Omega_b0(), ccl_cosmo_eh_linear["Omega_b"], rtol=reltol)
    assert_allclose(cosmo.Omega_c0(), ccl_cosmo_eh_linear["Omega_c"], rtol=reltol)
    assert_allclose(cosmo.Omega_k0(), ccl_cosmo_eh_linear["Omega_k"], rtol=reltol)
    assert_allclose(cosmo.T_gamma0(), ccl_cosmo_eh_linear["T_CMB"], rtol=reltol)
    assert_allclose(cosmo.Omega_g0(), ccl_cosmo_eh_linear["Omega_g"], rtol=reltol)
    assert_allclose(cosmo.Neff(), ccl_cosmo_eh_linear["Neff"], rtol=reltol_Neff)
    assert_allclose(
        cosmo.param_get_by_name("w0"), ccl_cosmo_eh_linear["w0"], rtol=reltol
    )
    assert_allclose(
        cosmo.param_get_by_name("w1"), ccl_cosmo_eh_linear["wa"], rtol=reltol
    )


def test_background_functions(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology
) -> None:
    """Compare NumCosmo and CCL background functions."""
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-9
        mnu_rtol = 1.0e-6
    else:
        rtol = 1.0e-9
        mnu_rtol = 1.0e-6

    z_test = np.linspace(0.0, 5.0, 2000)[1:]

    cmp = compare_Omega_m(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_Omega_g(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_Omega_k(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_Omega_nu(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_Omega_mnu(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=mnu_rtol)

    cmp = compare_Omega_de(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)


def test_background_h_over_h0(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology
) -> None:
    """Compare NumCosmo and CCL background functions."""
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-9
    else:
        rtol = 1.0e-5

    z_test = np.linspace(0.0, 5.0, 2000)[1:]
    cmp = compare_Hubble(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)


def test_background_distances(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology
) -> None:
    """Compare NumCosmo and CCL distances."""
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-10
        cvol_rtol = 1.0e-9
    else:
        rtol = 1.0e-6
        cvol_rtol = 1.0e-5

    z_test = np.linspace(0.0, 5.0, 2000)[1:]

    cmp = compare_distance_comoving(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_distance_transverse(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_distance_angular_diameter(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_distance_luminosity(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_distance_lookback_time(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_distance_modulus(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)

    cmp = compare_distance_comoving_volume(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=cvol_rtol)


def test_background_scale_factor_of_chi(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology
) -> None:
    """Compare NumCosmo and CCL scale factor of comoving distance."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-9
    else:
        rtol = 1.0e-6

    RH_Mpc = cosmo.RH_Mpc()

    c_max = dist.comoving(cosmo, 15.0)
    comoving = np.geomspace(1.0e-4, c_max, 4000)[1:]
    z_array = np.array([dist.inv_comoving(cosmo, c) for c in comoving])
    a_array = 1.0 / (1.0 + z_array)

    cmp = compare_scale_factor(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, comoving * RH_Mpc
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)
    assert_allclose(a_array, cmp.y1, rtol=rtol)


def test_background_growth_lowz(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology
) -> None:
    """Compare NumCosmo and CCL growth factor."""
    cosmo = nc_cosmo_eh_linear.cosmo
    ps_ml = nc_cosmo_eh_linear.ps_ml
    if ccl_cosmo_eh_linear.high_precision:
        reltol_growth = 1.0e-4
        reltol_rate = 1.0e-4
    else:
        reltol_growth = 1.0e-4
        reltol_rate = 1.0e-4

    k_pivot = 1.0
    z_test = np.linspace(0.0, 5.0, 2000)[1:]
    a_test = 1.0 / (1.0 + z_test)

    nc_D2_a = np.array(
        [
            ps_ml.eval(cosmo, z, k_pivot) / ps_ml.eval(cosmo, 0.0, k_pivot)
            for z in z_test
        ]
    )
    nc_dD2_a = np.array(
        [
            ps_ml.deriv_z(cosmo, z, k_pivot) / ps_ml.eval(cosmo, 0.0, k_pivot)
            for z in z_test
        ]
    )

    cmp_D = compare_growth_factor(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    cmp_f = compare_growth_rate(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)

    assert_allclose(cmp_D.y2**2, nc_D2_a, rtol=1.0e-15, atol=0.0)
    assert_allclose(
        nc_dD2_a / (2.0 * cmp_D.y2),
        -cmp_D.y2 * cmp_f.y2 * a_test,
        rtol=1.0e-15,
        atol=0.0,
    )

    assert_allclose(cmp_D.y1, cmp_D.y2, rtol=reltol_growth, atol=0.0)
    assert_allclose(cmp_f.y1, cmp_f.y2, rtol=reltol_rate, atol=0.0)


def test_background_growth_highz(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology
) -> None:
    """Compare NumCosmo and CCL growth factor."""
    if ccl_cosmo_eh_linear.high_precision:
        reltol_growth = 1.0e-1
        reltol_rate = 1.0e-1
    else:
        reltol_growth = 1.0e-1
        reltol_rate = 1.0e-1

    z_test = np.geomspace(0.001, 1100.0, 4000)[1:]

    cmp = compare_growth_factor(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=reltol_growth, atol=0.0)

    cmp = compare_growth_rate(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test)
    assert_allclose(cmp.y1, cmp.y2, rtol=reltol_rate, atol=0.0)


def test_background_Sigma_crit(
    ccl_cosmo_eh_linear: pyccl.Cosmology, nc_cosmo_eh_linear: ncpy.Cosmology
) -> None:
    """Compare NumCosmo and CCL critical density."""
    if ccl_cosmo_eh_linear.high_precision:
        rtol = 1.0e-10
    else:
        rtol = 1.0e-7

    z_test = np.linspace(1.0, 5.0, 2000)[1:]
    z_l = 0.5

    cmp = compare_Sigma_crit(ccl_cosmo_eh_linear, nc_cosmo_eh_linear, z_test, zl=z_l)
    assert_allclose(cmp.y1, cmp.y2, rtol=rtol)
