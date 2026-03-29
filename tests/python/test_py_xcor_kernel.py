#
# test_py_xcor_kernel.py
#
# Tue Feb 18 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor_kernel.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Unit tests for NcXcorKernel.

Tests all kernel types using parametrized kernel_case fixture. Kernels are collected
with all relevant configurations (e.g., galaxy with/without magbias).
"""

import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology

pytestmark = pytest.mark.xcor

Ncm.cfg_init()
pytest_plugins = ["python.fixtures_xcor"]


def test_limber_vs_limber_z(kernel: Nc.XcorKernel, cosmology: Cosmology) -> None:
    """Test that limber and limber_z give consistent results."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    ps_ml = cosmology.ps_ml
    kernel.set_l_limber(0)  # Treat all ells as limber
    kernel.prepare(cosmo)

    # Test ell values
    ell_array = np.array([10, 50, 100, 500, 1000])

    # Get vectorized evaluation functions
    for ell in ell_array:
        limber_func = kernel.get_eval(cosmo, ell)
        nu = ell + 0.5

        # Evaluate both methods
        kmin, kmax = limber_func.get_range()

        k_array = np.geomspace(kmin, kmax, 100)
        for k in k_array:
            xi = nu / k
            k_Mpc = k / cosmo.RH_Mpc()
            z = dist.inv_comoving(cosmo, xi)
            pk = ps_ml.eval(cosmo, z, k_Mpc)
            result = limber_func.eval_array(k)[0]
            result_z = (
                kernel.eval_limber_z_full(cosmo, z, dist, ell)
                * np.sqrt(pk)
                * np.sqrt(np.pi / 2.0 / nu)
                / k
            )
            max_result = np.max(np.abs(result))
            assert_allclose(
                result,
                result_z,
                rtol=1.0e-3,
                atol=1.0e-3 * max_result,
                err_msg=(
                    f"Limber and limber_z differ at ell={ell}, k={k:.3e}, z={z:.3f}"
                ),
            )


def test_limber_vs_limber_z_vectorized(
    kernel: Nc.XcorKernel, cosmology: Cosmology
) -> None:
    """Test that limber and limber_z give consistent results.

    Validates the vectorized evaluation method by comparing it against the direct
    limber_z_full computation across a range of multipoles."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    ps_ml = cosmology.ps_ml
    kernel.set_l_limber(0)  # Treat all ells as limber
    kernel.prepare(cosmo)

    # Test ell range
    ell_start = 100
    ell_end = ell_start + 8  # 9 elements total

    # Get vectorized evaluation function for ell range
    limber_func = kernel.get_eval_vectorized(cosmo, ell_start, ell_end)
    kmin, kmax = limber_func.get_range()

    # Test at multiple k values
    k_array = np.geomspace(kmin, kmax, 50)
    for k in k_array:
        # Get all 9 results at once
        results = limber_func.eval_array(k)

        # Compare each ell individually
        for i, ell in enumerate(range(ell_start, ell_end + 1)):
            nu = ell + 0.5
            xi = nu / k
            k_Mpc = k / cosmo.RH_Mpc()
            z = dist.inv_comoving(cosmo, xi)
            pk = ps_ml.eval(cosmo, z, k_Mpc)
            result_z = (
                kernel.eval_limber_z_full(cosmo, z, dist, ell)
                * np.sqrt(pk)
                * np.sqrt(np.pi / 2.0 / nu)
                / k
            )
            assert_allclose(
                results[i],
                result_z,
                rtol=1e-13,
                atol=0.0,
                err_msg=(
                    f"Vectorized limber and limber_z differ at "
                    f"ell={ell}, k={k:.3e}, z={z:.3f}"
                ),
            )


def test_limber_vs_non_limber(
    kernel_case: tuple[str, Nc.XcorKernel], cosmology: Cosmology
) -> None:
    """Test Limber approximation accuracy against exact non-Limber calculation.

    Validates that the Limber approximation converges to the exact (non-Limber) result
    at moderate to high multipoles (ell = 200, 500, 800). The test evaluates both
    methods over a k-range and verifies that the maximum relative difference stays
    within kernel-specific tolerances.

    The Limber approximation accuracy exhibits characteristic behavior:

    - Most accurate in the high-k region at moderate to high ell
    - Degrades toward lower k values where the flat-sky approximation breaks down
    - Improves with increasing ell (higher multipoles)

    Different kernel types have different tolerance thresholds based on their physical
    characteristics, ranging from ~1e-4 for weak lensing kernels to ~0.3 for more
    complex kernels like tSZ.
    """
    kernel_id, kernel = kernel_case
    cosmo = cosmology.cosmo
    kernel_tol = {
        "kernel_cmb_lens": {200: 5.0e-2, 500: 3.0e-2, 800: 2.8e-2},
        "kernel_cmb_isw": {200: 4.0, 500: 2.5, 800: 2.0},
        "kernel_tsz": {200: 0.3, 500: 0.3, 800: 0.3},
        "kernel_gal_bin0": {200: 4.0e-3, 500: 6.0e-4, 800: 3.0e-4},
        "kernel_gal_bin1": {200: 8.0e-3, 500: 2.0e-3, 800: 5.1e-4},
        "kernel_gal_bin2": {200: 2.0e-2, 500: 3.0e-3, 800: 9.6e-4},
        "kernel_gal_bin3": {200: 3.0e-2, 500: 4.0e-3, 800: 2.0e-3},
        "kernel_gal_bin4": {200: 4.0e-2, 500: 6.0e-3, 800: 3.0e-3},
        "kernel_wl_bin0": {200: 8.0e-4, 500: 2.0e-4, 800: 6.0e-5},
        "kernel_wl_bin1": {200: 2.0e-3, 500: 2.0e-4, 800: 8.0e-5},
        "kernel_wl_bin2": {200: 2.0e-3, 500: 3.0e-4, 800: 1.0e-4},
        "kernel_wl_bin3": {200: 2.0e-3, 500: 3.0e-4, 800: 2.0e-4},
        "kernel_wl_bin4": {200: 7.0e-4, 500: 2.0e-4, 800: 5.0e-5},
    }
    # Test at moderate to high ell where limber should be accurate
    ell_array = np.array([200, 500, 800])
    default_rtol = 1.0e-13

    for ell in ell_array:
        # Get limber result
        rtol = kernel_tol.get(kernel_id, {}).get(ell, default_rtol)
        kernel.set_l_limber(0)
        kernel.prepare(cosmo)
        limber_func = kernel.get_eval(cosmo, ell)
        kmin_limber, kmax_limber = limber_func.get_range()

        # Get non-limber result (this is the "exact" reference)
        kernel.set_l_limber(-1)
        kernel.prepare(cosmo)

        for _ in range(1):
            non_limber_func = kernel.get_eval(cosmo, ell)

        kmin_non_limber, kmax_non_limber = non_limber_func.get_range()
        kmin = max(kmin_limber, kmin_non_limber)
        kmax = min(kmax_limber, kmax_non_limber)

        # Generate 1000 k values for statistical validation
        k_array = np.geomspace(kmin, kmax, 1000)

        limber_result = np.array([limber_func.eval_array(k)[0] for k in k_array])
        non_limber_result = np.array(
            [non_limber_func.eval_array(k)[0] for k in k_array]
        )
        max_non_limber = np.max(np.abs(non_limber_result))
        rel_diff = np.abs((limber_result - non_limber_result) / max_non_limber)

        assert np.max(rel_diff) < rtol, (
            f"Limber and non-Limber differ for {kernel_id} at ell={ell}: "
            f"max_rel_diff={np.max(rel_diff):.6e} (rtol={rtol:.6e})"
        )


def test_k_projection_limber_vs_non_limber(
    kernel_case: tuple[str, Nc.XcorKernel], cosmology: Cosmology
) -> None:
    """Test k-space projection integral agreement between Limber and non-Limber.

    Validates that the k-space projection integral I = ∫ k³ K²(k) d(ln k) computed
    using the Limber approximation agrees with the exact non-Limber calculation at
    moderate to high multipoles. This tests whether the Limber approximation preserves
    the integrated power, which is the physically relevant quantity for cosmological
    observables like angular power spectra.

    The test evaluates both methods over a trimmed k-range (excluding boundary regions)
    and verifies that the integrated quantities agree within specified tolerance. This
    provides a stronger validation than pointwise comparison, as small local errors can
    cancel when integrated, demonstrating that the Limber approximation is suitable for
    computing physical observables even where pointwise differences exist.
    """
    kernel_id, kernel = kernel_case
    cosmo = cosmology.cosmo

    # Test at moderate to high ell where Limber should be accurate
    ell_array = np.array([100, 500, 800])

    kernel_tol = {
        "kernel_cmb_lens": {100: 1.0e-5, 500: 6.0e-5, 800: 8.0e-5},
        "kernel_cmb_isw": {100: 2.0, 500: 0.5, 800: 0.5},
        "kernel_tsz": {100: 2.0e-4, 500: 2.0e-4, 800: 8.0e-5},
        "kernel_gal_bin0": {100: 1.0e-3, 500: 6.0e-5, 800: 3.0e-5},
        "kernel_gal_bin1": {100: 2.0e-3, 500: 2.0e-4, 800: 6.0e-5},
        "kernel_gal_bin2": {100: 3.0e-3, 500: 3.0e-4, 800: 1.0e-4},
        "kernel_gal_bin3": {100: 4.0e-3, 500: 4.0e-4, 800: 2.0e-4},
        "kernel_gal_bin4": {100: 7.0e-3, 500: 7.1e-4, 800: 4.0e-4},
        "kernel_wl_bin0": {100: 4.0e-4, 500: 2.0e-5, 800: 7.0e-6},
        "kernel_wl_bin1": {100: 2.0e-4, 500: 8.0e-6, 800: 8.0e-6},
        "kernel_wl_bin2": {100: 7.0e-5, 500: 4.0e-6, 800: 2.0e-6},
        "kernel_wl_bin3": {100: 4.0e-5, 500: 8.0e-7, 800: 8.0e-7},
        "kernel_wl_bin4": {100: 1.0e-4, 500: 6.0e-6, 800: 4.0e-6},
    }
    default_rtol = 1.0e-13

    for ell in ell_array:
        # Get limber result
        rtol = kernel_tol.get(kernel_id, {}).get(ell, default_rtol)
        kernel.set_l_limber(0)
        kernel.prepare(cosmo)
        limber_func = kernel.get_eval(cosmo, ell)
        kmin_limber, kmax_limber = limber_func.get_range()

        # Get non-limber result (this is the "exact" reference)
        kernel.set_l_limber(-1)
        kernel.prepare(cosmo)
        non_limber_func = kernel.get_eval(cosmo, ell)
        kmin_non_limber, kmax_non_limber = non_limber_func.get_range()

        kmin = max(kmin_limber, kmin_non_limber)
        kmax = min(kmax_limber, kmax_non_limber)

        # Use 2000 k points for accurate integration
        k_array = np.geomspace(kmin, kmax, 5000)
        lnk = np.log(k_array)

        # Evaluate kernels
        limber_result = np.array([limber_func.eval_array(k)[0] for k in k_array])
        non_limber_result = np.array(
            [non_limber_func.eval_array(k)[0] for k in k_array]
        )

        # Compute k-space projection integrals: I = ∫ k^3 K^2(k) d(ln k)
        integrand_limber = k_array**3 * limber_result**2
        integrand_non_limber = k_array**3 * non_limber_result**2

        I_limber = np.trapezoid(integrand_limber, x=lnk)
        I_exact = np.trapezoid(integrand_non_limber, x=lnk)

        # Check relative agreement
        rel_diff = np.abs((I_limber - I_exact) / I_exact)

        assert rel_diff < rtol, (
            f"k-space projection integral differs for {kernel_id} at ell={ell}: "
            f"I_limber={I_limber:.6e}, I_exact={I_exact:.6e}, "
            f"rel_diff={rel_diff:.6e} (rtol={rtol:.6e})"
        )
