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

# Pinned to one worker under --dist loadgroup: this file is one of the
# xcor lane's memory peaks, and an xdist worker is its own session, so
# without this its cost is paid once per worker rather than once.
pytestmark = [pytest.mark.xcor, pytest.mark.xdist_group("kernel")]

Ncm.cfg_init()
pytest_plugins = ["python.fixtures_xcor"]

RELTOL = 1.0e-7
PEAK_EPSILON = 1.0e-4
# The projection tests compare the Limber approximation's error in the integrated power,
# a few 1e-6, which closures at PEAK_EPSILON move by about 3e-6. At 1e-6 they move it by
# under 1e-7 (lensing bins, ell 500 and 804, against closures at 1e-8), so the comparison
# is the physics.
CONVERGED_TOL = 1.0e-6

# The tests below that build non-Limber closures run the Levin ODE over
# x = k chi, and its cost -- both time and the spectral order the panels must
# carry -- is linear in the top of that range. The kernels' own domain rule
# stops around k = 147 Mpc^-1 with the fixture power spectrum, which for a
# multipole block of 16 means panels holding ~2e5 oscillations and gigabytes of
# operator storage, none of which makes the comparison below any sharper: it
# checks Limber against non-Limber pointwise in k, so a shorter domain compares
# fewer points, each unchanged. Capping the power spectrum at 10 Mpc^-1 for
# these tests costs 12x less memory and 5x less time at the same knot density
# (2711 knots against 2739). The unbounded range belongs in a stress test, not
# in a unit test.
NON_LIMBER_KMAX = 10.0


def _k2_integrals(
    integrand: Nc.XcorKernelIntegrand, k_min: float, k_max: float
) -> np.ndarray:
    """Return int_{k_min}^{k_max} k^2 W_r(k)^2 dk for every row of a closure.

    The closure is piecewise polynomial, and NcmSpectralPanelsProduct integrates the
    square of a panel set exactly. Row r lives in x = k / s_r, so the panels are
    restricted to [k_min / s_r, k_max / s_r] before its product is taken.
    """
    sp = integrand.peek_panels()
    n_rows = integrand.get_len()
    scales = Ncm.Vector.new_array([integrand.get_scale(r) for r in range(n_rows)])
    spectral = Ncm.Spectral.new()
    product = Ncm.SpectralPanelsProduct.new(2)
    result = Ncm.Vector.new(n_rows)
    out = np.zeros(n_rows)

    for r in range(n_rows):
        lo, hi = k_min / scales.get(r), k_max / scales.get(r)
        restricted = Ncm.SpectralPanels.new(sp.get_fdim())

        for i in range(sp.get_n_panels()):
            _, a, b = sp.get_panel(i)
            a_r, b_r = max(a, lo), min(b, hi)

            if a_r < b_r:
                inside, coeffs = sp.restrict(spectral, a_r, b_r)
                assert inside
                restricted.append(a_r, b_r, coeffs)

        if restricted.get_n_panels() > 0:
            product.integrate(restricted, scales, restricted, scales, result, None)
            out[r] = result.get(r)

    return out


@pytest.fixture(name="capped_powspec")
def fixture_capped_powspec(cosmology: Cosmology):
    """Cap the fixture power spectrum for one test, then restore it.

    Function-scoped on purpose: the power spectrum is shared with every other
    xcor test through the cached cosmology, so the cap must not outlive the
    test that asks for it.
    """
    ps_ml = cosmology.ps_ml
    original_kmax = ps_ml.get_kmax()

    ps_ml.set_kmax(NON_LIMBER_KMAX)
    ps_ml.prepare(cosmology.cosmo)

    yield ps_ml

    ps_ml.set_kmax(original_kmax)
    ps_ml.prepare(cosmology.cosmo)


def test_limber_vs_limber_z(kernel: Nc.XcorKernel, cosmology: Cosmology) -> None:
    """Test that limber and limber_z give consistent results."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    ps_ml = cosmology.ps_ml
    kernel.set_l_limber(0)  # Treat all ells as limber
    kernel.prepare(cosmo)

    # Test ell values
    ell_array = np.array([10, 50, 100, 500, 1000])

    max_xi = dist.max_comoving_distance()

    # Get vectorized evaluation functions
    for ell in ell_array:
        limber_func = kernel.get_eval(cosmo, ell)
        nu = ell + 0.5

        # Evaluate both methods
        kmin, kmax = limber_func.get_range()
        kmin = max(kmin, nu / max_xi)  # Ensure kmin does not exceed max_xi limit
        k_array = np.geomspace(kmin, kmax, 100)

        results_list = []
        results_z_list = []
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
            results_list.append(result)
            results_z_list.append(result_z)

        results = np.array(results_list)
        results_z = np.array(results_z_list)
        max_result = np.max(np.abs(results))
        assert_allclose(
            results,
            results_z,
            rtol=1.0e-3,
            atol=1.0e-3 * max_result,
            err_msg=f"Limber and limber_z differ at ell={ell}",
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
    kernel.set_reltol(RELTOL)
    kernel.set_peak_epsilon(PEAK_EPSILON)

    # Test ell range
    ell_start = 100
    ell_end = ell_start + 8  # 9 elements total

    # Get vectorized evaluation function for ell range
    limber_func = kernel.get_eval_vectorized(cosmo, ell_start, ell_end)
    kmin, kmax = limber_func.get_range()

    # Test at multiple k values
    max_xi = dist.max_comoving_distance()
    k_array = np.geomspace(kmin, kmax, 100)

    # Get kernel's valid z range
    zmin_kernel, zmax_kernel, _ = kernel.get_z_range()

    # Build full result arrays first
    results_vectorized_list = []
    results_reference_list = []

    for k in k_array:
        # Get all 9 results at once
        results = limber_func.eval_array(k)
        results_vectorized_list.append(results)

        # Compute reference values for each ell
        row_reference = []
        for ell in range(ell_start, ell_end + 1):
            nu = ell + 0.5
            xi = nu / k
            # Skip if xi is too close to or beyond the maximum comoving distance.
            # Use a 2% safety margin since the kernel base class extends zf to 1e10,
            # making values near max_xi numerically unstable when inverting.
            if xi >= 0.98 * max_xi:
                row_reference.append(np.nan)
                continue
            k_Mpc = k / cosmo.RH_Mpc()
            z = dist.inv_comoving(cosmo, xi)
            # Skip if z is outside the kernel's valid range (e.g., tSZ has zmax=6.0)
            if z < zmin_kernel or z > zmax_kernel:
                row_reference.append(np.nan)
                continue
            pk = ps_ml.eval(cosmo, z, k_Mpc)
            result_z = (
                kernel.eval_limber_z_full(cosmo, z, dist, ell)
                * np.sqrt(pk)
                * np.sqrt(np.pi / 2.0 / nu)
                / k
            )
            row_reference.append(result_z)
        results_reference_list.append(row_reference)

    # Convert to arrays
    results_vectorized = np.array(results_vectorized_list)
    results_reference = np.array(results_reference_list)

    # Create mask for valid (non-NaN) reference values
    valid_mask = ~np.isnan(results_reference)

    # Compute max and atol using only valid values
    max_result = np.max(np.abs(results_reference[valid_mask]))
    peak_epsilon = kernel.get_peak_epsilon()
    atol = max_result * peak_epsilon

    # Compare arrays (only valid elements)
    assert_allclose(
        results_vectorized[valid_mask],
        results_reference[valid_mask],
        rtol=1.0e-5,
        atol=atol,
        err_msg="Vectorized limber and limber_z differ",
    )


def test_limber_vs_non_limber_vectorized(
    kernel_case: tuple[str, Nc.XcorKernel],
    cosmology: Cosmology,
    capped_powspec: Ncm.Powspec,
) -> None:
    """Test Limber vs non-Limber using vectorized evaluation with blocks of ells.

    Validates that vectorized Limber approximation agrees with non-Limber calculation
    when evaluating multiple multipoles simultaneously. Tests blocks of 8 and 16 ells
    at moderate to high multipole ranges.
    """
    kernel_id, kernel = kernel_case
    cosmo = cosmology.cosmo
    kernel_tol = {
        "kernel_cmb_lens": {200: 5.0e-2, 500: 3.0e-2, 800: 2.8e-2},
        "kernel_cmb_isw": {200: 4.0, 500: 2.5, 800: 2.0},
        "kernel_tsz": {200: 0.3, 500: 0.3, 800: 0.3},
        "kernel_cluster_tophat": {200: 0.63, 500: 0.6, 800: 0.6},
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

    # Test with blocks of 8 and 16 ells starting at different positions
    test_configs = [
        (200, 8),  # 8 ells starting at 200
        (500, 8),  # 8 ells starting at 500
        (800, 16),  # 16 ells starting at 800
    ]

    default_rtol = 1.0e-13
    kernel.set_reltol(RELTOL)
    kernel.set_peak_epsilon(PEAK_EPSILON)

    for ell_start, n_ells in test_configs:
        ell_end = ell_start + n_ells - 1
        rtol = kernel_tol.get(kernel_id, {}).get(ell_start, default_rtol)

        # Get limber result with vectorized evaluation
        kernel.set_l_limber(0)
        kernel.prepare(cosmo)
        limber_func = kernel.get_eval_vectorized(cosmo, ell_start, ell_end)
        kmin_limber, kmax_limber = limber_func.get_range()

        # Get non-limber result with vectorized evaluation
        kernel.set_l_limber(-1)
        kernel.prepare(cosmo)
        non_limber_func = kernel.get_eval_vectorized(cosmo, ell_start, ell_end)
        kmin_non_limber, kmax_non_limber = non_limber_func.get_range()

        kmin = max(kmin_limber, kmin_non_limber)
        kmax = min(kmax_limber, kmax_non_limber)

        # Generate 500 k values for validation
        k_array = np.geomspace(kmin, kmax, 500)

        # Evaluate both methods - each returns array of n_ells values per k
        limber_results = np.array([limber_func.eval_array(k) for k in k_array])
        non_limber_results = np.array([non_limber_func.eval_array(k) for k in k_array])

        # Compare across all k and ell values
        max_non_limber = np.max(np.abs(non_limber_results))
        rel_diff = np.abs((limber_results - non_limber_results) / max_non_limber)
        max_rel_diff = np.max(rel_diff)

        assert max_rel_diff < rtol, (
            f"Vectorized Limber and non-Limber differ for {kernel_id} "
            f"at ell_range=[{ell_start}, {ell_end}] (n_ells={n_ells}): "
            f"max_rel_diff={max_rel_diff:.6e} (rtol={rtol:.6e})"
        )


def test_k_projection_limber_vs_non_limber_vectorized(
    kernel_case: tuple[str, Nc.XcorKernel],
    cosmology: Cosmology,
    capped_powspec: Ncm.Powspec,
) -> None:
    """Test k-space projection with vectorized evaluation for blocks of ells.

    Validates that k-space projection integrals I = int k^3 K^2(k) d(ln k) computed
    using vectorized Limber approximation agree with non-Limber calculation when
    evaluating multiple multipoles simultaneously. Tests blocks of 8 and 16 ells.
    """
    kernel_id, kernel = kernel_case
    cosmo = cosmology.cosmo

    kernel_tol = {
        "kernel_cmb_lens": {100: 2.0e-5, 500: 6.0e-5, 800: 8.0e-5},
        "kernel_cmb_isw": {100: 2.0, 500: 0.5, 800: 0.5},
        "kernel_tsz": {100: 2.0e-4, 500: 2.0e-4, 800: 8.0e-5},
        "kernel_cluster_tophat": {100: 4.0e-2, 500: 2.0e-2, 800: 1.0e-2},
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

    # Test with blocks of 8 and 16 ells
    test_configs = [
        (100, 8),  # 8 ells starting at 100
        (500, 8),  # 8 ells starting at 500
        (800, 16),  # 16 ells starting at 800
    ]

    default_rtol = 1.0e-13
    kernel.set_reltol(CONVERGED_TOL)
    kernel.set_peak_epsilon(CONVERGED_TOL)

    for ell_start, n_ells in test_configs:
        ell_end = ell_start + n_ells - 1
        rtol = kernel_tol.get(kernel_id, {}).get(ell_start, default_rtol)

        # Get limber result with vectorized evaluation
        kernel.set_l_limber(0)
        kernel.prepare(cosmo)
        limber_func = kernel.get_eval_vectorized(cosmo, ell_start, ell_end)
        kmin_limber, kmax_limber = limber_func.get_range()

        # Get non-limber result with vectorized evaluation
        kernel.set_l_limber(-1)
        kernel.prepare(cosmo)
        non_limber_func = kernel.get_eval_vectorized(cosmo, ell_start, ell_end)
        kmin_non_limber, kmax_non_limber = non_limber_func.get_range()

        kmin = max(kmin_limber, kmin_non_limber)
        kmax = min(kmax_limber, kmax_non_limber)

        # I = int k^3 K^2(k) d(ln k) = int k^2 K^2 dk per multipole, exact on the panels
        I_limber_block = _k2_integrals(limber_func, kmin, kmax)
        I_exact_block = _k2_integrals(non_limber_func, kmin, kmax)

        for i in range(n_ells):
            ell = ell_start + i
            I_limber = I_limber_block[i]
            I_exact = I_exact_block[i]

            # Check relative agreement for this ell
            rel_diff = np.abs((I_limber - I_exact) / I_exact)

            assert rel_diff < rtol, (
                f"k-space projection integral differs for {kernel_id} "
                f"at ell={ell} (from vectorized range [{ell_start}, {ell_end}]): "
                f"I_limber={I_limber:.6e}, I_exact={I_exact:.6e}, "
                f"rel_diff={rel_diff:.6e} (rtol={rtol:.6e})"
            )
