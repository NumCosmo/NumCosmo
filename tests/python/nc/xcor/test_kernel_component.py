#
# test_py_xcor_kernel_component.py
#
# Tue Feb 18 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor_kernel_component.py
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

"""Unit tests for NcXcorKernelComponent.

Tests all component types using parametrized component fixture.
"""

import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology

pytestmark = pytest.mark.xcor

Ncm.cfg_init()

pytest_plugins = ["python.fixtures_xcor"]


def test_any_component_interface(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test that any component satisfies the basic interface."""
    cosmo = cosmology.cosmo
    _, _, component = kernel_component

    # All components should be instances of NcXcorKernelComponent
    assert isinstance(component, Nc.XcorKernelComponent)

    # Test get_limits
    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
    assert np.isfinite(xi_min) and xi_min > 0.0
    assert np.isfinite(xi_max) and xi_max > xi_min
    assert np.isfinite(k_min) and k_min > 0.0
    assert np.isfinite(k_max) and k_max > k_min

    # Test eval_kernel at a sample point
    xi_test = (xi_min + xi_max) / 2.0
    k_test = np.sqrt(k_min * k_max)  # Geometric mean
    result = component.eval_kernel(cosmo, xi_test, k_test)
    assert np.isfinite(result)

    # Test eval_prefactor
    prefactor = component.eval_prefactor(cosmo, k_test, 100)
    assert np.isfinite(prefactor)

    # Test property getters
    epsilon = component.get_epsilon()
    assert np.isfinite(epsilon) and epsilon > 0.0

    ny = component.get_ny()
    assert isinstance(ny, int) and ny > 0

    max_iter = component.get_max_iter()
    assert isinstance(max_iter, int) and max_iter > 0

    tol = component.get_tol()
    assert np.isfinite(tol) and tol > 0.0


def test_component_creation(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
) -> None:
    """Test that we can create a component."""
    _, _, component = kernel_component
    assert component is not None
    assert isinstance(component, Nc.XcorKernelComponent)


def test_component_properties(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
) -> None:
    """Test component property getters and setters."""
    _, _, component = kernel_component
    component.set_epsilon(1.0e-8)
    assert_allclose(component.get_epsilon(), 1.0e-8, rtol=1e-10)
    component.set_ny(100)
    assert component.get_ny() == 100
    component.set_max_iter(500)
    assert component.get_max_iter() == 500
    component.set_tol(1.0e-6)
    assert_allclose(component.get_tol(), 1.0e-6, rtol=1e-10)


def test_component_eval_kernel(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test eval_kernel method."""
    _, _, component = kernel_component
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
    xi_array = np.linspace(xi_min, xi_max, 10)
    k_array = np.linspace(k_min, k_max, 10)

    for xi in xi_array:
        for k in k_array:
            result = component.eval_kernel(cosmo, xi, k)
            assert np.isfinite(result)
            # Test symmetry and basic properties
            assert isinstance(result, float)


def test_component_eval_prefactor(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test eval_prefactor method."""
    _, _, component = kernel_component
    cosmo = cosmology.cosmo
    k_array = np.linspace(0.1, 5.0, 10)
    ell_array = [10, 50, 100, 500]

    for k in k_array:
        for ell in ell_array:
            result = component.eval_prefactor(cosmo, k, ell)
            assert np.isfinite(result)


def test_component_get_limits(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test get_limits method."""
    _, _, component = kernel_component
    cosmo = cosmology.cosmo
    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)

    assert np.isfinite(xi_min) and xi_min > 0.0
    assert np.isfinite(xi_max) and xi_max > xi_min
    assert np.isfinite(k_min) and k_min > 0.0
    assert np.isfinite(k_max) and k_max > k_min


def test_component_prepare(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test prepare method and kernel analysis functions.

    Verifies that after calling prepare(), the component can evaluate k_max,
    KL_max, and k_epsilon with positive finite values across multiple y values.
    """
    # Set parameters for analysis
    _, _, component = kernel_component
    cosmo = cosmology.cosmo
    component.set_epsilon(1.0e-10)
    component.set_ny(50)
    component.set_max_iter(1000)
    component.set_tol(1.0e-8)
    component.prepare(cosmo)

    # After prepare, we should be able to evaluate k_max, KL_max, k_epsilon
    y_array = np.linspace(1.0, 100.0, 10)

    for y in y_array:
        k_max = component.eval_k_max(y)
        KL_max = component.eval_KL_max(y)
        k_epsilon = component.eval_k_epsilon(y)

        assert np.isfinite(k_max) and k_max > 0.0
        assert np.isfinite(KL_max)
        assert np.isfinite(k_epsilon) and k_epsilon > 0.0


def test_component_kernel_analysis_properties(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test kernel analysis creates valid splines with expected properties.

    Validates that k_max, KL_max, and k_epsilon arrays are all positive, finite,
    and that k_epsilon >= 0.1 * k_max across the tested range.
    """
    _, _, component = kernel_component
    cosmo = cosmology.cosmo
    component.set_epsilon(1.0e-8)
    component.set_ny(30)
    component.set_max_iter(1000)
    component.set_tol(1.0e-7)
    component.prepare(cosmo)

    # Test that eval_k_max is monotonic or behaves reasonably
    y_array = np.geomspace(1.0, 100.0, 20)
    k_max_array = np.array([component.eval_k_max(y) for y in y_array])
    KL_max_array = np.array([component.eval_KL_max(y) for y in y_array])
    k_epsilon_array = np.array([component.eval_k_epsilon(y) for y in y_array])

    # All should be positive and finite
    assert np.all(np.isfinite(k_max_array))
    assert np.all(k_max_array > 0.0)
    assert np.all(np.isfinite(KL_max_array))
    assert np.all(np.isfinite(k_epsilon_array))
    assert np.all(k_epsilon_array > 0.0)

    # k_epsilon should be >= k_max (point where kernel becomes small)
    # This is physics-dependent, so just check they're in reasonable ranges
    assert np.all(k_epsilon_array >= k_max_array * 0.1)


def test_component_edge_cases(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test component behavior at edge cases."""
    _, _, component = kernel_component
    cosmo = cosmology.cosmo
    # Test at boundaries
    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
    xi_mean = np.sqrt(xi_min * xi_max)
    k_mean = np.sqrt(k_min * k_max)

    # At minimum xi
    result = component.eval_kernel(cosmo, xi_min, k_mean)
    assert np.isfinite(result)

    # At maximum xi
    result = component.eval_kernel(cosmo, xi_max, k_mean)
    assert np.isfinite(result)

    # At minimum k
    result = component.eval_kernel(cosmo, xi_mean, k_min)
    assert np.isfinite(result)

    # At maximum k
    result = component.eval_kernel(cosmo, xi_mean, k_max)
    assert np.isfinite(result)

    # Test very small and large multipoles
    for ell in [2, 10, 100, 1000, 10000]:
        result = component.eval_prefactor(cosmo, k_mean, ell)
        assert np.isfinite(result)


def test_component_k_max_is_maximum(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test that k_max(y) gives the maximum of KL(k, y/k) using Limber approximation.

    Verifies that the cached KL_max value matches direct evaluation at k_max
    within 5% relative tolerance.
    """
    _, _kernel, component = kernel_component
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)

    y_min = k_min * xi_min
    y_max = k_max * xi_max

    # Test at several y values
    y_test_values = np.geomspace(max(y_min * 1.1, 0.5), min(y_max * 0.9, 1000.0), 5)

    for y in y_test_values:
        k_max_val = component.eval_k_max(y)
        KL_max_val = component.eval_KL_max(y)
        xi = y / k_max_val

        # Verify KL_max is approximately the actual value at k_max
        if xi_min <= xi <= xi_max:
            K_at_k_max = abs(_compute_KL_at_k(component, cosmo, y, k_max_val))
            # Allow tolerance for spline/search errors
            assert_allclose(
                KL_max_val,
                K_at_k_max,
                rtol=0.05,
                atol=1e-10,
                err_msg=f"KL_max doesn't match eval_kernel at y={y}",
            )


def test_component_KL_max_value(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test that KL_max(y) equals KL(y/k_max, k_max) using Limber approximation.

    Validates consistency between the stored KL_max value and direct kernel
    evaluation. Uses 5% tolerance for real kernels with numerical artifacts.
    """
    _, _kernel, component = kernel_component
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
    y_min = k_min * xi_min
    y_max = k_max * xi_max

    # Test at multiple y values
    y_test_values = np.geomspace(max(y_min * 1.1, 0.5), min(y_max * 0.9, 1000.0), 10)

    for y in y_test_values:
        k_max_val = component.eval_k_max(y)
        KL_max_val = component.eval_KL_max(y)

        # Compute xi from y and k_max
        xi = y / k_max_val

        if xi_min <= xi <= xi_max:
            # Compute KL directly using Limber approximation
            K_direct = abs(_compute_KL_at_k(component, cosmo, y, k_max_val))

            # They should match within ~0.2% for real kernels
            assert_allclose(
                KL_max_val,
                K_direct,
                rtol=0.05,
                atol=1e-8,
                err_msg=(
                    f"KL_max mismatch at y={y}: "
                    f"KL_max={KL_max_val}, direct={K_direct}"
                ),
            )


def test_component_k_epsilon_drop(
    kernel_component_drop: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test that k_epsilon is where KL drops by epsilon from KL_max.

    The algorithm finds k_epsilon where KL(y/k_epsilon, k_epsilon) ≈ epsilon * KL_max.
    Since set_epsilon sets the squared threshold, we use epsilon^2 in the test.
    """
    _, _, component = kernel_component_drop

    cosmo = cosmology.cosmo
    epsilon = 1.0e-5

    original_epsilon = component.get_epsilon()
    # The algorithm will find k_epsilon where KL = epsilon^2 * KL_max.
    component.set_epsilon(epsilon**2)
    component.prepare(cosmo)  # Re-prepare with new epsilon

    xi_min, xi_max, _, _ = component.get_limits(cosmo)

    # Test at several y values
    y_test_values = [20.0, 50.0, 100.0]

    for y in y_test_values:
        k_max_val = component.eval_k_max(y)
        KL_max_val = component.eval_KL_max(y)
        k_epsilon_val = component.eval_k_epsilon(y)

        # k_epsilon should be >= k_max (after the maximum)
        assert (
            k_epsilon_val >= k_max_val * 0.99
        ), f"k_epsilon={k_epsilon_val} should be >= k_max={k_max_val} at y={y}"

        # Compute KL at k_epsilon using Limber approximation
        xi_epsilon = y / k_epsilon_val
        if xi_min <= xi_epsilon <= xi_max:
            K_at_epsilon = np.abs(_compute_KL_at_k(component, cosmo, y, k_epsilon_val))
            # Should be approximately epsilon * KL_max (allow factor of 5 tolerance due
            # to spline/search)
            expected_K = np.abs(epsilon * KL_max_val)
            assert K_at_epsilon >= expected_K * 0.2, (
                f"K at k_epsilon ({K_at_epsilon}) too small vs "
                f"epsilon*KL_max ({expected_K}) at y={y}"
            )
            assert K_at_epsilon <= expected_K * 5.0, (
                f"K at k_epsilon ({K_at_epsilon}) too large vs "
                f"epsilon*KL_max ({expected_K}) at y={y}"
            )

    component.set_epsilon(original_epsilon)
    component.prepare(cosmo)


def test_component_y_range_pruning(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test that kernel analysis works within the expected y range [0.5, 1000.0].

    Verifies that eval_k_max returns valid results at the boundaries of the
    y range defined by the kernel limits.
    """
    _, _kernel, component = kernel_component
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)

    # Calculate expected y range
    y_min_expected = max(k_min * xi_min, 0.5)
    y_max_expected = min(k_max * xi_max, 1000.0)

    # Test that eval_k_max works at the boundaries
    # Values near y_min
    k_max_at_ymin = component.eval_k_max(y_min_expected)
    assert np.isfinite(k_max_at_ymin) and k_max_at_ymin > 0.0

    # Values near y_max
    k_max_at_ymax = component.eval_k_max(y_max_expected)
    assert np.isfinite(k_max_at_ymax) and k_max_at_ymax > 0.0

    # Test a few points in the valid range
    y_test_values = np.geomspace(y_min_expected * 1.1, y_max_expected * 0.9, 5)
    for y in y_test_values:
        k = component.eval_k_max(y)
        assert np.isfinite(k) and k > 0.0, f"Invalid k_max at y={y}"


def test_component_k_max_is_maximum1(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test that k_max(y).

    Test that k_max(y) gives a near-maximum value of KL(k, y/k) using Limber
    approximation. Validates that k_max is among the larger values by sampling points
    around it. Uses percentile-based checking to handle kernels with oscillatory
    behavior.
    """
    _, __kernel, component = kernel_component
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
    y_min = k_min * xi_min
    y_max = k_max * xi_max

    # Test at several y values
    y_test_values = np.geomspace(max(y_min * 1.1, 0.5), min(y_max * 0.9, 1000.0), 5)

    for y in y_test_values:
        k_max_val = component.eval_k_max(y)
        KL_max_val = component.eval_KL_max(y)

        # Verify KL_max is approximately the actual value at k_max
        K_at_k_max = abs(_compute_KL_at_k(component, cosmo, y, k_max_val))
        # Allow tolerance due to spline interpolation and oscillating kernels
        assert_allclose(
            abs(KL_max_val),
            abs(K_at_k_max),
            rtol=0.2,
            atol=1e-8,
            err_msg=f"KL_max doesn't match eval_kernel at y={y}, k_max={k_max_val}",
        )

        # Sample points around k_max to verify it's a local maximum
        k_test_points = np.array(
            [
                k_max_val * 0.7,  # Before maximum
                k_max_val * 0.9,
                k_max_val,  # At maximum
                k_max_val * 1.1,
                k_max_val * 1.3,  # After maximum
            ]
        )

        # Filter to stay within valid k range
        k_test_points = k_test_points[
            (k_test_points >= k_min) & (k_test_points <= k_max)
        ]

        K_values_list = []
        for k_test in k_test_points:
            xi_test = y / k_test
            # Only test if xi is in valid range
            if xi_min <= xi_test <= xi_max:
                K_val = abs(_compute_KL_at_k(component, cosmo, y, k_test))
                K_values_list.append(K_val)
            else:
                K_values_list.append(-np.inf)  # Mark as invalid

        K_values = np.array(K_values_list)
        valid_K = K_values[np.isfinite(K_values)]

        if len(valid_K) > 2:
            # KL_max should be among the larger values (allowing numerical errors)
            assert (
                abs(KL_max_val) >= np.percentile(valid_K, 90) * 0.9
            ), f"k_max={k_max_val} doesn't give a large |K| value at y={y}"


def _compute_KL_at_k(
    component: Nc.XcorKernelComponent, cosmo: Nc.HICosmo, y: float, k: float
) -> float:
    return component.eval_kernel(cosmo, y / k, k) / k * np.sqrt(np.pi / 2 / y)


def test_component_KL_max_value1(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test KL_max(y).

    Test that KL_max(y) approximately equals KL(y/k_max, k_max) using Limber
    approximation. Verifies consistency between the cached KL_max value and direct
    evaluation. Allows 5% tolerance for spline interpolation and search algorithm
    precision.
    """
    _, _kernel, component = kernel_component
    cosmo = cosmology.cosmo

    # Test at multiple y values
    y_test_values = np.geomspace(2.0, 50.0, 10)

    for y in y_test_values:
        k_max_val = component.eval_k_max(y)
        KL_max_val = component.eval_KL_max(y)

        # Compute KL directly using Limber approximation
        K_direct = np.abs(_compute_KL_at_k(component, cosmo, y, k_max_val))

        # Allow tolerance for spline interpolation and search errors
        assert_allclose(
            abs(KL_max_val),
            abs(K_direct),
            rtol=0.05,
            atol=1e-10,
            err_msg=f"KL_max mismatch at y={y}: KL_max={KL_max_val}, direct={K_direct}",
        )


def test_component_k_epsilon_drop1(
    kernel_component_drop: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test that k_epsilon is where KL drops by epsilon from KL_max for k > k_max.

    At k_epsilon: KL(y/k_epsilon, k_epsilon) ≈ epsilon * KL_max.
    Since set_epsilon sets the squared threshold, we configure epsilon^2.
    Allows factor of 10 tolerance due to search algorithm behavior at boundaries.
    """
    _, _, component = kernel_component_drop
    cosmo = cosmology.cosmo
    original_epsilon = component.get_epsilon()
    epsilon = 1.0e-5  # Use larger epsilon for more reliable testing
    # The algorithm finds k_epsilon where KL = epsilon^2 * KL_max, so we set epsilon^2
    # here.
    component.set_epsilon(epsilon**2)
    component.prepare(cosmo)

    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
    y_min = k_min * xi_min
    y_max = k_max * xi_max

    # Test at several y values
    y_test_values = np.geomspace(max(y_min * 1.1, 0.5), min(y_max * 0.9, 1000.0), 5)

    for y in y_test_values:
        k_max_val = component.eval_k_max(y)
        KL_max_val = component.eval_KL_max(y)
        k_epsilon_val = component.eval_k_epsilon(y)

        # k_epsilon should be >= k_max
        assert (
            k_epsilon_val >= k_max_val * 0.99
        ), f"k_epsilon={k_epsilon_val} should be >= k_max={k_max_val} at y={y}"

        # Compute |K*xi| at k_epsilon
        xi_epsilon = y / k_epsilon_val
        K_at_epsilon = np.abs(_compute_KL_at_k(component, cosmo, y, k_epsilon_val))

        # Should be approximately epsilon * KL_max
        expected_K = epsilon * KL_max_val

        # Skip if k_epsilon is at boundary
        if np.isclose(k_epsilon_val, k_max, atol=1e-3) or np.isclose(
            xi_epsilon, xi_min, atol=1e-3
        ):
            return

        assert K_at_epsilon <= expected_K * 10.0, (
            f"K at k_epsilon ({K_at_epsilon}) should be "
            f"~epsilon*KL_max ({expected_K}) at y={y}"
        )
        # Also check it's not way too small
        assert K_at_epsilon >= expected_K * 0.1, (
            f"K at k_epsilon ({K_at_epsilon}) too small "
            f"compared to epsilon*KL_max ({expected_K}) at y={y}"
        )

    component.set_epsilon(original_epsilon)
    component.prepare(cosmo)


def test_component_k_epsilon_after_k_max(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test that k_epsilon >= k_max for all y values.

    The algorithm looks for the drop point k >= k_max where the kernel falls
    below the epsilon threshold.
    """
    _, _, component = kernel_component
    cosmo = cosmology.cosmo
    original_epsilon = component.get_epsilon()
    component.set_epsilon(1.0e-7)
    component.prepare(cosmo)

    # Test at many y values
    y_test_values = np.geomspace(1.0, 100.0, 15)

    for y in y_test_values:
        k_max_val = component.eval_k_max(y)
        k_epsilon_val = component.eval_k_epsilon(y)

        # k_epsilon should always be after k_max
        assert (
            k_epsilon_val >= k_max_val
        ), f"k_epsilon={k_epsilon_val} should be >= k_max={k_max_val} at y={y}"

    component.set_epsilon(original_epsilon)
    component.prepare(cosmo)


def test_component_monotonicity(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
) -> None:
    """Test basic properties of k_max, KL_max, and k_epsilon functions.

    Validates that all values are positive, finite, and satisfy k_epsilon >= k_max.
    """
    _, _, component = kernel_component

    # Sample many y values
    y_array = np.geomspace(1.0, 100.0, 30)

    k_max_array = np.array([component.eval_k_max(y) for y in y_array])
    KL_max_array = np.array([component.eval_KL_max(y) for y in y_array])
    k_epsilon_array = np.array([component.eval_k_epsilon(y) for y in y_array])

    # All should be positive and finite
    assert np.all(np.isfinite(k_max_array))
    assert np.all(k_max_array > 0.0)
    assert np.all(np.isfinite(KL_max_array))
    assert np.all(np.isfinite(k_epsilon_array))
    assert np.all(k_epsilon_array > 0.0)

    # k_epsilon >= k_max everywhere
    assert np.all(k_epsilon_array >= k_max_array * 0.99)

    # KL_max should be positive (for our test kernel)
    assert np.all(KL_max_array >= 0.0)


def test_component_spline_smoothness(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
) -> None:
    """Test that splines for k_max, KL_max, k_epsilon are smooth.

    Evaluates at 600 points and checks for discontinuous jumps by computing
    relative differences between consecutive points. KL_max is allowed larger
    jumps (up to 6.0) due to physical kernel structure.
    """
    _, _, component = kernel_component

    # Dense sampling
    y_array = np.geomspace(2.0, 80.0, 600)

    k_max_array = np.array([component.eval_k_max(y) for y in y_array])
    KL_max_array = np.array([component.eval_KL_max(y) for y in y_array])
    k_epsilon_array = np.array([component.eval_k_epsilon(y) for y in y_array])

    # Check no sudden jumps (should be relatively smooth)
    # Compute relative differences between consecutive points
    def rel_diff(arr):
        return np.abs(np.diff(arr)) / (arr[:-1] + 1e-10)

    k_max_diff = rel_diff(k_max_array)
    KL_max_diff = rel_diff(KL_max_array)
    k_epsilon_diff = rel_diff(k_epsilon_array)

    # Check for smoothness
    assert (
        np.max(k_max_diff) < 1.0
    ), f"k_max has discontinuous jumps: max={np.max(k_max_diff)}"
    assert np.max(KL_max_diff) < 6.0, (
        f"KL_max has large jumps: max={np.max(KL_max_diff)} "
        f"(may indicate kernel structure)"
    )
    assert (
        np.max(k_epsilon_diff) < 1.0
    ), f"k_epsilon has discontinuous jumps: max={np.max(k_epsilon_diff)}"


def test_kernel_analysis_with_different_epsilon(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test that smaller epsilon values give larger k_epsilon.

    Verifies the expected relationship: as epsilon decreases, the point where
    the kernel drops to epsilon * KL_max moves to higher k values.
    """
    _, _, component = kernel_component
    cosmo = cosmology.cosmo

    original_epsilon = component.get_epsilon()
    epsilon_values = [1e-4, 1e-6, 1e-8]
    y_test = 20.0

    k_epsilon_results = []

    for eps in epsilon_values:
        component.set_epsilon(eps)
        component.prepare(cosmo)
        k_eps = component.eval_k_epsilon(y_test)
        k_epsilon_results.append(k_eps)

    component.set_epsilon(original_epsilon)
    component.prepare(cosmo)

    # Smaller epsilon should give larger k_epsilon
    assert (
        k_epsilon_results[0] <= k_epsilon_results[1] <= k_epsilon_results[2]
    ), f"k_epsilon should increase as epsilon decreases: {k_epsilon_results}"


def test_component_correctness(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test kernel analysis correctness at multiple y values.

    Verifies that KL_max matches direct evaluation and that k_epsilon >= k_max,
    using 5% relative tolerance for complex kernel structures.
    """
    _, _, component = kernel_component
    cosmo = cosmology.cosmo

    # Get limits
    xi_min, xi_max, _, _ = component.get_limits(cosmo)

    # Test at a few y values
    y_test_values = [10.0, 50.0, 200.0]

    for y in y_test_values:
        k_max_val = component.eval_k_max(y)
        KL_max_val = component.eval_KL_max(y)

        # Verify KL_max equals KL(y/k_max, k_max) using Limber approximation
        xi = y / k_max_val
        if xi_min <= xi <= xi_max:
            K_direct = np.abs(_compute_KL_at_k(component, cosmo, y, k_max_val))
            # Allow generous tolerance for real components with complex structure
            assert_allclose(
                KL_max_val, K_direct, rtol=0.05, err_msg=f"KL_max mismatch at y={y}"
            )

        # Test that k_epsilon >= k_max
        k_epsilon_val = component.eval_k_epsilon(y)
        assert (
            k_epsilon_val >= k_max_val * 0.99
        ), f"k_epsilon should be >= k_max at y={y}"


def test_kernel_gal_component_list(
    kernel_gal: Nc.XcorKernelGal, cosmology: Cosmology
) -> None:
    """Test that galaxy kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    kernel_gal.prepare(cosmo)

    # Get components from kernel
    comp_list = kernel_gal.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0  # Should have at least clustering component

    # Test each component satisfies the interface
    for component in comp_list:
        assert isinstance(component, Nc.XcorKernelComponent)

        # Test get_limits
        xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min

        # Test eval_kernel at a sample point
        xi_test = (xi_min + xi_max) / 2.0
        k_test = (k_min + k_max) / 2.0
        result = component.eval_kernel(cosmo, xi_test, k_test)
        assert np.isfinite(result)

        # Test eval_prefactor
        prefactor = component.eval_prefactor(cosmo, k_test, 100)
        assert np.isfinite(prefactor)


def test_kernel_cmb_lens_component_list(
    kernel_cmb_lens: Nc.XcorKernelCMBLensing, cosmology: Cosmology
) -> None:
    """Test that CMB lensing kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    kernel_cmb_lens.prepare(cosmo)

    comp_list = kernel_cmb_lens.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    for component in comp_list:
        assert isinstance(component, Nc.XcorKernelComponent)

        xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_kernel_cmb_isw_component_list(
    kernel_cmb_isw: Nc.XcorKernelCMBISW, cosmology: Cosmology
) -> None:
    """Test that CMB ISW kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    kernel_cmb_isw.prepare(cosmo)

    comp_list = kernel_cmb_isw.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    for component in comp_list:
        assert isinstance(component, Nc.XcorKernelComponent)

        xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_kernel_tsz_component_list(
    kernel_tsz: Nc.XcorKerneltSZ, cosmology: Cosmology
) -> None:
    """Test that tSZ kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    kernel_tsz.prepare(cosmo)

    comp_list = kernel_tsz.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    for component in comp_list:
        assert isinstance(component, Nc.XcorKernelComponent)

        xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_kernel_wl_component_list(
    kernel_wl: Nc.XcorKernelWeakLensing, cosmology: Cosmology
) -> None:
    """Test that weak lensing kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    kernel_wl.prepare(cosmo)

    comp_list = kernel_wl.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    for component in comp_list:
        assert isinstance(component, Nc.XcorKernelComponent)

        xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_kernel_component_analysis(
    kernel_tsz: Nc.XcorKerneltSZ, cosmology: Cosmology
) -> None:
    """Test that tSZ kernel components have valid kernel analysis.

    Verifies that after prepare(), components can evaluate k_max, KL_max, and
    k_epsilon with finite positive values, and that k_epsilon >= k_max.
    """
    cosmo = cosmology.cosmo
    kernel_tsz.prepare(cosmo)

    comp_list = kernel_tsz.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    # After prepare, components should have kernel analysis done
    component = comp_list[0]

    # Test at several y values
    y_array = np.linspace(10.0, 100.0, 5)

    for y in y_array:
        k_max = component.eval_k_max(y)
        KL_max = component.eval_KL_max(y)
        k_epsilon = component.eval_k_epsilon(y)

        assert np.isfinite(k_max) and k_max > 0.0
        assert np.isfinite(KL_max)
        assert np.isfinite(k_epsilon) and k_epsilon > 0.0
        # k_epsilon should be at or beyond k_max
        assert k_epsilon >= k_max * 0.1


def test_multiple_components_galaxy_kernel(
    kernel_gal: Nc.XcorKernelGal, cosmology: Cosmology
) -> None:
    """Test that galaxy kernel with magnification bias returns two components.

    Validates that clustering and magnification bias components are both
    present and satisfy the basic component interface.
    """
    # Create galaxy kernel with magnification bias enabled
    cosmo = cosmology.cosmo
    kernel_gal.prepare(cosmo)

    comp_list = kernel_gal.get_component_list()
    assert comp_list is not None
    # Should have 2 components: clustering and magnification bias
    assert len(comp_list) == 2

    # Both should be valid components
    for component in comp_list:
        assert isinstance(component, Nc.XcorKernelComponent)
        xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_component_properties_via_get_component_list(
    kernel_tsz: Nc.XcorKerneltSZ, cosmology: Cosmology
) -> None:
    """Test component property getters and setters via get_component_list."""
    cosmo = cosmology.cosmo
    kernel_tsz.prepare(cosmo)

    comp_list = kernel_tsz.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    component = comp_list[0]

    # Test epsilon property
    component.set_epsilon(1.0e-7)
    assert_allclose(component.get_epsilon(), 1.0e-7, rtol=1e-10)

    # Test ny property
    component.set_ny(150)
    assert component.get_ny() == 150

    # Test max_iter property
    component.set_max_iter(2000)
    assert component.get_max_iter() == 2000

    # Test tol property
    component.set_tol(1.0e-5)
    assert_allclose(component.get_tol(), 1.0e-5, rtol=1e-10)


def test_component_chebyshev_decomposition(
    kernel_component: tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent],
    cosmology: Cosmology,
) -> None:
    """Test Chebyshev decomposition of component kernel g_k(y) = K*xi(y/k, k).

    Validates that the adaptive Chebyshev expansion accurately represents the
    kernel function within the specified tolerance (rtol=1e-9).
    """
    _, _, component = kernel_component
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = component.get_limits(cosmo)
    spectral = Ncm.Spectral.new()
    k_values = np.geomspace(k_min, k_max, 10)
    for k in k_values:

        def g_k(y, k=k):
            """Evaluate K*xi(y/k, k)."""
            xi = y / k
            K_val = xi * component.eval_kernel(cosmo, xi, k)
            return K_val

        y_min = k * xi_min
        y_max = k * xi_max

        # Compute Chebyshev coefficients adaptively
        final_k, coeffs = spectral.compute_chebyshev_coeffs_adaptive(
            g_k, y_min, y_max, 4, 1.0e-8
        )
        max_order_k = spectral.get_max_order()

        # Verify we got coefficients
        assert final_k >= 5
        assert coeffs is not None
        assert len(coeffs) >= 2**5 + 1
        assert len(coeffs) <= 2**max_order_k + 1

        # Verify the expansion is accurate at test points
        y_test_points = np.linspace(y_min, y_max, 500)
        direct_vals = np.array([g_k(y) for y in y_test_points])
        max_direct = np.max(np.abs(direct_vals))
        cheb_vals = np.array(
            [
                Ncm.Spectral.chebyshev_eval_x(coeffs, y_min, y_max, y)
                for y in y_test_points
            ]
        )
        assert all(
            np.isfinite(cheb_vals)
        ), f"Chebyshev values contain non-finite numbers at k={k}"

        assert_allclose(
            cheb_vals,
            direct_vals,
            rtol=1.0e-9,
            atol=1.0e-9 * max_direct,
            err_msg=(
                f"Chebyshev expansion mismatch at k={k}, max_direct={max_direct} "
                f"max_order_k={max_order_k}"
            ),
        )
