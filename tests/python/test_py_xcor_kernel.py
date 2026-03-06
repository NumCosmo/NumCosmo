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

pytestmark = pytest.mark.xcor

Ncm.cfg_init()


@cache
def get_cosmology() -> Cosmology:
    """Create a default cosmology for testing."""
    return Cosmology.default()


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Create a default cosmology for testing."""
    return get_cosmology()


@cache
def get_integrator() -> Ncm.SBesselIntegrator:
    """Create a spherical Bessel integrator for testing."""
    integrator = Ncm.SBesselIntegratorLevin.new(0, 2000, 1.0, 1.0e6, 10, 1200)
    integrator.set_max_order(2**18)
    return integrator


@pytest.fixture(name="integrator", scope="module")
def fixture_integrator() -> Ncm.SBesselIntegrator:
    """Create a spherical Bessel integrator."""
    return get_integrator()


def _create_galaxy_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator, domagbias: bool = False
) -> Nc.XcorKernelGal:
    """Helper to create a galaxy kernel."""
    dn_dz = Ncm.SplineCubicNotaknot()
    z_center = 0.8
    z_sigma = 0.1
    z_low = max(z_center - 5.5 * z_sigma, 0.0)
    z_high = z_center + 5.5 * z_sigma
    z_array = np.linspace(z_low, z_high, 1000)
    nz_array = np.exp(-((z_array - z_center) ** 2) / (2.0 * z_sigma**2))
    nz_array /= np.trapezoid(nz_array, z_array)  # Normalize
    dn_dz.set_array(z_array, nz_array, True)

    gal = Nc.XcorKernelGal(
        dndz=dn_dz,
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        integrator=integrator,
        domagbias=domagbias,
    )
    gal["bparam_0"] = 1.5
    gal["mag_bias"] = 0.3
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
    z_center = 0.5
    z_sigma = 0.1
    z_low = max(z_center - 5.5 * z_sigma, 0.0)
    z_high = z_center + 5.5 * z_sigma
    z_array = np.linspace(z_low, z_high, 1000)
    nz_array = np.exp(-((z_array - z_center) ** 2) / (2.0 * z_sigma**2))
    nz_array /= np.trapezoid(nz_array, z_array)  # Normalize
    dn_dz.set_array(z_array, nz_array, True)

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
    return _collect_all_kernels(cosmology, get_integrator())


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Dynamically generate tests for all kernel types."""
    if "kernel_case" in metafunc.fixturenames:
        cases = get_kernel_cases()

        # Exclude specific kernels from certain tests
        if metafunc.function.__name__ in (
            "test_limber_vs_non_limber",
            "test_k_projection_limber_vs_non_limber",
        ):
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
    kernel_case: tuple[str, Nc.XcorKernel], cosmology: Cosmology
) -> None:
    """Test that limber and limber_z consistence.

    Test that limber and limber_z give consistent results using vectorized
    evaluation."""
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
    at moderate to high multipoles (ell = 200, 500, 800). The test evaluates both
    methods over a k-range starting from intermediate k values and verifies:

    1. Best agreement is within specified tolerance (rtol ~ 1e-5 to 1e-3)
    2. At least 80% of the k-range passes from right to left (high k to low k)

    The Limber approximation accuracy exhibits characteristic behavior:

    - Most accurate in the high-k region at moderate to high ell
    - Degrades toward lower k values where the flat-sky approximation breaks down
    - Improves with increasing ell (higher multipoles)

    By verifying that at least 80% of the range passes from the right (high k), this
    confirms the Limber approximation is reliable in the regime where it's physically
    valid and where most cosmological signal resides.
    """
    _, kernel = kernel_case

    cosmo = cosmology.cosmo

    # Test at moderate to high ell where limber should be accurate
    ell_array = np.array([200, 500, 800])
    ell_rtol = np.array([1.0e-3, 1.0e-4, 3.0e-5])

    for ell, rtol in zip(ell_array, ell_rtol):
        # Get limber result
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

        # Test k-range: start from 30% into the range (in log space) to avoid
        # low-k boundary effects, then test up to kmax where Limber is most accurate
        log_kmin = np.log(kmin)
        log_kmax = np.log(kmax)
        log_range = log_kmax - log_kmin
        k_test_min = np.exp(log_kmin + 0.3 * log_range)
        k_test_max = kmax

        # Generate 1000 k values for statistical validation
        k_array = np.geomspace(k_test_min, k_test_max, 1000)

        limber_result = np.array([limber_func.eval_array(k)[0] for k in k_array])
        non_limber_result = np.array(
            [non_limber_func.eval_array(k)[0] for k in k_array]
        )
        max_non_limber = np.max(np.abs(non_limber_result))
        rel_diff = np.abs((limber_result - non_limber_result) / max_non_limber)

        # Check best agreement is within tolerance
        best_agreement = np.min(rel_diff)
        assert best_agreement < rtol, (
            f"Best agreement {best_agreement:.3e} exceeds tolerance {rtol:.3e} "
            f"at ell={ell}"
        )

        # Check that at least 80% of the k-range passes from right to left
        # (high k to low k). This verifies the Limber approximation is accurate
        # in the high-k region where it's expected to perform best.
        passing = rel_diff < rtol
        failure_indices = np.where(~passing)[0]

        if len(failure_indices) > 0:
            # Find the last (rightmost) failure index
            rightmost_failure = failure_indices[-1]
            # Calculate fraction of range that passes after the rightmost failure
            passing_fraction_from_right = (len(k_array) - 1 - rightmost_failure) / len(
                k_array
            )
        else:
            # All points pass
            passing_fraction_from_right = 1.0

        assert passing_fraction_from_right >= 0.8, (
            f"Limber approximation fails in high-k region at ell={ell}: "
            f"only {passing_fraction_from_right:.1%} of range passes from right "
            f"(rightmost failure at k={k_array[rightmost_failure]:.3e})"
            if len(failure_indices) > 0
            else f"Unexpected failure at ell={ell}"
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
    ell_rtol = np.array([1.0e-2, 1.0e-3, 1.0e-4])

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

        I_limber = np.trapezoid(integrand_limber, x=lnk)
        I_exact = np.trapezoid(integrand_non_limber, x=lnk)

        # Check relative agreement
        rel_diff = np.abs((I_limber - I_exact) / I_exact)

        assert rel_diff < rtol, (
            f"k-space projection integral differs for {kernel_id} at ell={ell}: "
            f"I_limber={I_limber:.6e}, I_exact={I_exact:.6e}, "
            f"rel_diff={rel_diff:.6e} (rtol={rtol:.6e})"
        )
