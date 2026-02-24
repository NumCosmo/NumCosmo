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

from typing import Callable
from functools import cache
import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology

Ncm.cfg_init()


@cache
def get_cosmology() -> Cosmology:
    """Create a default cosmology for testing."""
    return Cosmology.default()


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Create a default cosmology for testing."""
    return get_cosmology()


@pytest.fixture(name="integrator", scope="module")
def fixture_integrator() -> Ncm.SBesselIntegrator:
    """Create a spherical Bessel integrator."""
    return Ncm.SBesselIntegratorLevin(max_order=2**16, lmin=0, lmax=2000)


def _create_galaxy_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator, domagbias: bool = False
) -> Nc.XcorKernelGal:
    """Helper to create a galaxy kernel."""
    dn_dz = Ncm.SplineCubicNotaknot()
    z_array = np.linspace(0.1, 2.0, 50)
    dndz_array = np.exp(-((z_array - 0.8) ** 2) / (2.0 * 0.3**2))
    dndz_array /= np.trapezoid(dndz_array, z_array)

    xv = Ncm.Vector.new_array(z_array.tolist())
    yv = Ncm.Vector.new_array(dndz_array.tolist())
    dn_dz.set(xv, yv, True)

    gal = Nc.XcorKernelGal(
        dndz=dn_dz,
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        integrator=integrator,
        domagbias=domagbias,
    )
    gal["bparam_0"] = 1.5
    return gal


def _create_cmb_lensing_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelCMBLensing:
    """Helper to create a CMB lensing kernel."""
    lmax = 1000
    ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
    return Nc.XcorKernelCMBLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=ell_vec,
        integrator=integrator,
    )


def _create_cmb_isw_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelCMBISW:
    """Helper to create a CMB ISW kernel."""
    lmax = 1000
    ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
    return Nc.XcorKernelCMBISW(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=ell_vec,
        integrator=integrator,
    )


def _create_tsz_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKerneltSZ:
    """Helper to create a tSZ kernel."""
    return Nc.XcorKerneltSZ(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        zmax=6.0,
        integrator=integrator,
    )


def _create_weak_lensing_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelWeakLensing:
    """Helper to create a weak lensing kernel."""
    dn_dz = Ncm.SplineCubicNotaknot()
    z_array = np.linspace(0.1, 3.0, 50)
    dndz_array = (z_array / 1.0) ** 2 * np.exp(-((z_array / 1.0) ** 1.5))
    dndz_array /= np.trapezoid(dndz_array, z_array)

    xv = Ncm.Vector.new_array(z_array.tolist())
    yv = Ncm.Vector.new_array(dndz_array.tolist())
    dn_dz.set(xv, yv, True)

    return Nc.XcorKernelWeakLensing(
        dndz=dn_dz,
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        integrator=integrator,
    )


def _collect_all_kernels(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> list[tuple[str, Nc.XcorKernel]]:
    """Collect all kernel types with all configurations.

    Returns a list of tuples: (kernel_id, kernel)
    """
    cosmo = cosmology.cosmo
    all_kernels: list[tuple[str, Nc.XcorKernel]] = []

    # Define all kernel types with configurations
    kernels_to_create: list[tuple[str, Callable[[], Nc.XcorKernel]]] = [
        ("Gal", lambda: _create_galaxy_kernel(cosmology, integrator, domagbias=False)),
        (
            "Gal:MagBias",
            lambda: _create_galaxy_kernel(cosmology, integrator, domagbias=True),
        ),
        ("CMBLensing", lambda: _create_cmb_lensing_kernel(cosmology, integrator)),
        ("CMBISW", lambda: _create_cmb_isw_kernel(cosmology, integrator)),
        ("tSZ", lambda: _create_tsz_kernel(cosmology, integrator)),
        ("WeakLensing", lambda: _create_weak_lensing_kernel(cosmology, integrator)),
    ]

    # Create each kernel
    for kernel_id, kernel_factory in kernels_to_create:
        kernel = kernel_factory()
        kernel.prepare(cosmo)
        all_kernels.append((kernel_id, kernel))

    return all_kernels


@cache
def get_kernel_cases():
    """Get all kernel cases for parametrization."""
    cosmology = get_cosmology()
    integrator = Ncm.SBesselIntegratorLevin.new(0, 2000)
    return _collect_all_kernels(cosmology, integrator)


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Dynamically generate tests for all kernel types."""
    if "kernel_case" in metafunc.fixturenames:
        cases = get_kernel_cases()

        # Exclude specific kernels from certain tests
        if metafunc.function.__name__ == "test_limber_vs_non_limber":
            # CMBISW has pathological Limber behavior and will be tested separately
            filtered_cases = [
                (kernel_id, kernel)
                for kernel_id, kernel in cases
                if kernel_id not in ["CMBISW"]
            ]
            ids = [kernel_id for kernel_id, _ in filtered_cases]
            metafunc.parametrize("kernel_case", filtered_cases, ids=ids)
        else:
            ids = [kernel_id for kernel_id, _ in cases]
            metafunc.parametrize("kernel_case", cases, ids=ids)


def test_limber_vs_limber_z(
    kernel_case: tuple[str, Nc.XcorKernel], cosmology: Cosmology
) -> None:
    """Test that limber and limber_z give consistent results."""
    _, kernel = kernel_case
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
            assert_allclose(
                result,
                result_z,
                rtol=1e-13,
                atol=0.0,
                err_msg=(
                    f"Limber and limber_z differ at ell={ell}, k={k:.3e}, z={z:.3f}"
                ),
            )


def test_limber_vs_limber_z_vectorized(
    kernel_case: tuple[str, Nc.XcorKernel], cosmology: Cosmology
) -> None:
    """Test that limber and limber_z give consistent results using vectorized evaluation."""
    _, kernel = kernel_case
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
    at moderate to high multipoles (ell = 100, 200, 500). The test evaluates both
    methods over a trimmed k-range (excluding boundary regions) and verifies:

    1. Best agreement is within specified tolerance (rtol ~ 1e-6)
    2. At least 50% of test points pass a relaxed tolerance (1e-3)
    3. Any failures are exclusively at edge regions where boundary effects dominate

    The Limber approximation accuracy exhibits characteristic behavior:
    - Most accurate in the central k-region at moderate to high ell
    - Degrades near integration range boundaries due to truncation effects
    - Improves with increasing ell (higher multipoles)

    By confirming failures concentrate exclusively at edges, this validates that the
    Limber approximation is well-behaved and reliable in the physically relevant bulk
    region where cosmological observables are measured.
    """
    _, kernel = kernel_case

    cosmo = cosmology.cosmo

    # Test at moderate to high ell where limber should be accurate
    ell_array = np.array([100, 200, 500])
    ell_rtol = np.array([1.0e-6, 1.0e-6, 1.0e-6])
    best_agreements = []

    for ell, rtol in zip(ell_array, ell_rtol):
        # Get limber result
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

        # Trim boundaries: skip first 15% and last 20% (in log space)
        log_kmin = np.log(kmin)
        log_kmax = np.log(kmax)
        log_range = log_kmax - log_kmin
        k_test_min = np.exp(log_kmin + 0.15 * log_range)
        k_test_max = np.exp(log_kmax - 0.20 * log_range)

        # Test at 15 k values in the reliable region
        k_array = np.geomspace(k_test_min, k_test_max, 1000)

        limber_result = np.array([limber_func.eval_array(k)[0] for k in k_array])
        non_limber_result = np.array(
            [non_limber_func.eval_array(k)[0] for k in k_array]
        )
        rel_diff = (limber_result - non_limber_result) / non_limber_result
        abs_rel_error = np.abs(rel_diff)

        best_agreements = np.min(abs_rel_error)
        assert best_agreements < rtol

        passing = abs_rel_error < rtol * 1.0e3
        passing_fraction = np.mean(passing) * 100
        assert passing_fraction > 50.0

        # Verify that failures are concentrated at edges where Limber approximation is
        # known to degrade. The Limber approximation has characteristic behavior:
        #
        # - Most accurate in the central k-region at moderate to high ell
        # - Degrades near boundaries due to integration range truncation effects
        #
        # This test confirms that if failures occur, they're exclusively at edges,
        # validating that the Limber approximation is well-behaved in the bulk region.
        if not np.all(passing):
            failing_lnk = np.log(k_array[~passing])
            # Define edge region as 18% from each boundary (in log-space)
            edge_threshold = 0.18 * (np.log(k_test_max) - np.log(k_test_min))
            edge_mask = (failing_lnk < np.log(k_test_min) + edge_threshold) | (
                failing_lnk > np.log(k_test_max) - edge_threshold
            )
            edge_failing_fraction = np.mean(edge_mask) * 100
            assert edge_failing_fraction == 100.0, (
                f"Non-passing points should be exclusively at "
                f"edges where Limber degrades. "
                f"Edge failing fraction: "
                f"{edge_failing_fraction:.1f}% (expected 100.0%)"
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
    ell_array = np.array([100, 200, 500])
    ell_rtol = np.array([1.0e-3, 1.0e-3, 1.0e-4])

    for ell, rtol in zip(ell_array, ell_rtol):
        # Get limber result
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

        # Trim boundaries: skip first 15% and last 20% (in log space)
        log_kmin = np.log(kmin)
        log_kmax = np.log(kmax)
        log_range = log_kmax - log_kmin
        k_test_min = np.exp(log_kmin + 0.0 * log_range)
        k_test_max = np.exp(log_kmax - 0.25 * log_range)

        # Use 2000 k points for accurate integration
        k_array = np.geomspace(k_test_min, k_test_max, 5000)
        lnk = np.log(k_array)

        # Evaluate kernels
        limber_result = np.array([limber_func.eval_array(k)[0] for k in k_array])
        non_limber_result = np.array(
            [non_limber_func.eval_array(k)[0] for k in k_array]
        )

        # Compute k-space projection integrals: I = ∫ k^3 K^2(k) d(ln k)
        integrand_limber = k_array**3 * limber_result**2
        integrand_non_limber = k_array**3 * non_limber_result**2

        # Print the maximum of the integrands for debugging, print both value and
        # index, also print the place they most diverge
        reldiff = np.abs(
            (integrand_limber - integrand_non_limber) / integrand_non_limber
        )
        max_index = np.argmax(reldiff)
        max_k = k_array[max_index]
        max_diff = reldiff[max_index]
        print()
        print(
            f"Max diff for {kernel_id} at ell={ell}: index {max_index}, "
            f"max_k={max_k:.6e}, max_diff={max_diff:.6e}"
        )
        # Print the following 10 knots after the max index to see if the divergence is localized
        for i in range(max(max_index - 10, 0), min(max_index + 10, len(k_array))):
            print(
                f"  k={k_array[i]:.6e}, limber={limber_result[i]:.6e}, "
                f"non_limber={non_limber_result[i]:.6e}, "
                f"reldiff={reldiff[i]:.6e}"
            )

        I_limber = np.trapezoid(integrand_limber, x=lnk)
        I_exact = np.trapezoid(integrand_non_limber, x=lnk)

        # Check relative agreement
        rel_diff = np.abs((I_limber - I_exact) / I_exact)

        assert rel_diff < rtol, (
            f"k-space projection integral differs for {kernel_id} at ell={ell}: "
            f"I_limber={I_limber:.6e}, I_exact={I_exact:.6e}, "
            f"rel_diff={rel_diff:.6e} (rtol={rtol:.6e})"
        )
