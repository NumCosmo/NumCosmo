#
# test_spectral_weighted.py
#
# Tue Apr 29 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_spectral_weighted.py
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

"""Tests for NcmSpectral weighted Chebyshev coefficient computation for integrals."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

pytest.importorskip("scipy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from scipy.integrate import quad

from numcosmo_py import Ncm


class TestSpectralWeighted:
    """Tests for NcmSpectral weighted coefficient methods for integral computation."""

    @pytest.fixture(name="spectral")
    def fixture_spectral(self) -> Ncm.Spectral:
        """Create a spectral object."""
        return Ncm.Spectral.new()

    def test_weighted_constant_integral(self, spectral: Ncm.Spectral) -> None:
        """Test that weighted coefficients give correct integral for constant function.

        For f(x) = c (constant), the integral over [a,b] is c*(b-a).
        In the weighted API, coeffs[0] already includes the interval scaling,
        so the integral is obtained as π * coeffs[0].
        """

        def f_constant(_user_data, _x):
            return 2.5

        a, b = -1.0, 1.0
        k_min = 2
        tol = 1.0e-12

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_constant, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # The first coefficient times π should equal the integral
        expected_integral = 2.5 * (b - a)
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-9, atol=1.0e-10)

    def test_weighted_linear_integral(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients for linear function.

        For f(x) = x on [-1, 1], integral is 0.
        For f(x) = x on [0, 2], integral is 2.
        """

        def f_linear(_user_data, x):
            return x

        # Test on symmetric interval [-1, 1]: integral should be 0
        a1, b1 = -1.0, 1.0
        k_min = 2
        tol = 1.0e-12

        _k1, coeffs1_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_linear, a1, b1, k_min, tol, None
        )
        coeffs1 = np.array(coeffs1_list)

        # Integral of x from -1 to 1 is 0
        computed_integral1 = coeffs1[0] * np.pi
        assert_allclose(computed_integral1, 0.0, rtol=1.0e-9, atol=1.0e-10)

        # Test on [0, 2]: integral should be 2
        a2, b2 = 0.0, 2.0

        _k2, coeffs2_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_linear, a2, b2, k_min, tol, None
        )
        coeffs2 = np.array(coeffs2_list)

        # Integral of x from 0 to 2 is x^2/2 |_0^2 = 2
        computed_integral2 = coeffs2[0] * np.pi
        assert_allclose(computed_integral2, 2.0, rtol=1.0e-9, atol=1.0e-10)

    def test_weighted_quadratic_integral(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients for quadratic function."""

        def f_quadratic(_user_data, x):
            return x * x

        a, b = -1.0, 1.0
        k_min = 2
        tol = 1.0e-12

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_quadratic, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Integral of x^2 from -1 to 1 is x^3/3 |_{-1}^1 = 2/3
        expected_integral = 2.0 / 3.0
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-9, atol=1.0e-10)

    def test_weighted_gaussian_integral(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients for Gaussian function."""

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        a, b = -3.0, 3.0
        k_min = 3
        tol = 1.0e-10

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_gaussian, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Compute reference integral using scipy
        expected_integral, _error = quad(lambda x: np.exp(-x * x), a, b)

        # First coefficient times π should match the integral
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-9, atol=1.0e-10)

    def test_weighted_sin_integral(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients for sine function."""

        def f_sin(_user_data, x):
            return np.sin(x)

        # Test on [0, pi]: integral is 2
        a, b = 0.0, np.pi
        k_min = 3
        tol = 1.0e-10

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_sin, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Integral of sin(x) from 0 to pi is 2
        expected_integral = 2.0
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-9, atol=1.0e-10)

    def test_weighted_exp_integral(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients for exponential function."""

        def f_exp(_user_data, x):
            return np.exp(x)

        a, b = 0.0, 1.0
        k_min = 3
        tol = 1.0e-10

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_exp, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Integral of exp(x) from 0 to 1 is e - 1
        expected_integral = np.exp(1.0) - 1.0
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-9, atol=1.0e-10)

    def test_weighted_vs_scipy_integration(self, spectral: Ncm.Spectral) -> None:
        """Test weighted first coefficient matches scipy quad for various functions."""

        test_cases = [
            (lambda _ud, x: 1.0, -1.0, 1.0, "constant"),
            (lambda _ud, x: x**3, -2.0, 2.0, "cubic"),
            (lambda _ud, x: np.exp(-0.5 * x * x), -2.0, 2.0, "gaussian"),
            (lambda _ud, x: np.sin(2.0 * x), 0.0, np.pi, "sine"),
            (lambda _ud, x: 1.0 / (1.0 + x * x), -1.0, 1.0, "rational"),
            (lambda _ud, x: np.cos(3.0 * x) * np.exp(-x * x), -2.0, 2.0, "product"),
        ]

        k_min = 4
        tol = 1.0e-10

        for func, a, b, name in test_cases:
            _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
                func, a, b, k_min, tol, None
            )
            coeffs = np.array(coeffs_list)

            # Compute reference integral
            expected_integral, _error = quad(lambda x, f=func: f(None, x), a, b)

            # First coefficient times π should match the integral
            computed_integral = coeffs[0] * np.pi
            assert_allclose(
                computed_integral,
                expected_integral,
                rtol=1.0e-8,
                atol=1.0e-9,
                err_msg=f"Failed for {name}",
            )

    def test_weighted_different_intervals(self, spectral: Ncm.Spectral) -> None:
        """Test that weighted coefficients work correctly on different intervals."""

        def f(_user_data, x):
            return np.exp(-x * x)

        intervals = [(-1.0, 1.0), (0.0, 2.0), (-3.0, 3.0), (0.5, 2.5)]

        k_min = 4
        tol = 1.0e-10

        for a, b in intervals:
            _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
                f, a, b, k_min, tol, None
            )
            coeffs = np.array(coeffs_list)

            # Compute reference integral
            expected_integral, _error = quad(lambda x: np.exp(-x * x), a, b)

            # First coefficient times π should match the integral
            computed_integral = coeffs[0] * np.pi
            assert_allclose(
                computed_integral,
                expected_integral,
                rtol=1.0e-8,
                atol=1.0e-9,
                err_msg=f"Failed for interval [{a}, {b}]",
            )

    def test_weighted_adaptive_convergence(self, spectral: Ncm.Spectral) -> None:
        """Test that weighted adaptive method converges properly."""

        def f_sharp(_user_data, x):
            """Sharp Gaussian peak."""
            return np.exp(-10.0 * (x - 0.5) ** 2)

        a, b = -1.0, 2.0
        k_min = 4
        tol = 1.0e-10

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_sharp, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Compute reference integral
        expected_integral, _error = quad(lambda x: f_sharp(None, x), a, b)

        # Should achieve good accuracy
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_weighted_returns_correct_k(self, spectral: Ncm.Spectral) -> None:
        """Test that weighted adaptive returns correct k value."""

        def f(_user_data, x):
            return np.exp(-x * x)

        a, b = -2.0, 2.0
        k_min = 3
        tol = 1.0e-10

        k_final, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        N = len(coeffs)

        # Verify N = 2^k + 1
        expected_N = (1 << k_final) + 1
        assert (
            N == expected_N
        ), f"Length mismatch: N={N}, expected 2^{k_final}+1={expected_N}"

        # k should be reasonable
        assert k_min <= k_final <= spectral.get_max_order()

    def test_weighted_vs_unweighted_coefficients(self, spectral: Ncm.Spectral) -> None:
        """Test that weighted and unweighted coefficients are different.

        The weighted version computes coefficients of F(x)*sqrt(1-t^2), so they
        should differ from the standard coefficients of F(x).
        """

        def f(_user_data, x):
            return np.exp(-x * x)

        a, b = -1.0, 1.0
        k_min = 4
        tol = 1.0e-10

        # Compute weighted coefficients
        _k_w, coeffs_w_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f, a, b, k_min, tol, None
        )
        coeffs_w = np.array(coeffs_w_list)

        # Compute standard coefficients at same order
        N_w = len(coeffs_w)
        _k_std, coeffs_std_list = spectral.compute_chebyshev_coeffs_adaptive(
            f, a, b, k_min, tol, None
        )
        coeffs_std = np.array(coeffs_std_list)

        # If orders match, verify coefficients are different
        if len(coeffs_std) == N_w:
            # They should NOT be equal (different weight functions)
            # Check that at least some coefficients differ significantly
            max_diff = np.max(np.abs(coeffs_w - coeffs_std))
            assert (
                max_diff > 1.0e-6
            ), "Weighted and unweighted coefficients should differ"

        # But the first weighted coefficient should still approximate the integral
        expected_integral, _error = quad(lambda x: f(None, x), a, b)
        computed_integral = coeffs_w[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_weighted_polynomial_exact(self, spectral: Ncm.Spectral) -> None:
        """Test that weighted method is exact for low-degree polynomials.

        For polynomials of degree <= N-1, both standard and weighted Chebyshev
        approximations should be exact (up to numerical precision).
        """

        def f_poly3(_user_data, x):
            """Cubic polynomial."""
            return 1.0 + 2.0 * x - 0.5 * x * x + 0.3 * x * x * x

        a, b = -1.0, 1.0
        k_min = 3  # N = 9, enough for degree 3
        tol = 1.0e-12

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_poly3, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Compute analytical integral: x + x^2 - x^3/6 + 3x^4/40 from -1 to 1
        # = (1 + 1 - 1/6 + 3/40) - (-1 + 1 + 1/6 + 3/40) = 2 - 1/3 = 5/3
        expected_integral = 1.0 * 2.0 + 0.0 - 0.5 * (2.0 / 3.0) + 0.0
        # int(1) = 2, int(x) = 0, int(x^2) = 2/3, int(x^3) = 0

        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-9, atol=1.0e-10)

    def test_weighted_oscillatory_integral(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients for oscillatory function."""

        def f_osc(_user_data, x):
            return np.sin(5.0 * x) * np.exp(-x * x)

        a, b = -2.0, 2.0
        k_min = 5
        tol = 1.0e-10

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_osc, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Compute reference integral
        expected_integral, _error = quad(lambda x: f_osc(None, x), a, b)

        # First coefficient times π should match the integral
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_weighted_sharp_peak_integral(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients for function with sharp peak."""

        def f_peak(_user_data, x):
            """Very sharp Gaussian."""
            return np.exp(-20.0 * (x - 0.3) ** 2)

        a, b = -1.0, 1.0
        k_min = 5
        tol = 1.0e-9

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_peak, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Compute reference integral
        expected_integral, _error = quad(lambda x: f_peak(None, x), a, b)

        # Should still be accurate
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-7, atol=1.0e-8)

    def test_weighted_multiple_refinements(self, spectral: Ncm.Spectral) -> None:
        """Test that weighted method can handle multiple refinement levels."""

        def f(_user_data, x):
            return np.cos(10.0 * x)

        a, b = 0.0, np.pi
        k_min = 3
        tol = 1.0e-9

        k_final, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Should need several refinements for cos(10x)
        assert k_final > k_min, "Should require refinement for oscillatory function"

        # Check integral accuracy
        expected_integral, _error = quad(lambda x: np.cos(10.0 * x), a, b)
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_weighted_boundary_values(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients with specific boundary values."""

        def f_boundary(_user_data, x):
            """Function that is 1 at boundaries, 0 at center."""
            return 1.0 - 4.0 * (x - 0.5) ** 2

        a, b = 0.0, 1.0
        k_min = 3
        tol = 1.0e-10

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_boundary, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Compute reference integral
        # int(1 - 4(x-0.5)^2) from 0 to 1 = int(1 - 4x^2 + 4x - 1) = int(4x - 4x^2)
        # = 2x^2 - 4x^3/3 |_0^1 = 2 - 4/3 = 2/3
        expected_integral = 2.0 / 3.0
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-9, atol=1.0e-10)

    def test_weighted_asymmetric_function(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients for asymmetric function."""

        def f_asym(_user_data, x):
            """Asymmetric function."""
            return x * np.exp(-x * x) if x > 0 else -x * np.exp(-x * x)

        a, b = -2.0, 2.0
        k_min = 4
        tol = 1.0e-10

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_asym, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # Compute reference integral
        expected_integral, _error = quad(lambda x: f_asym(None, x), a, b)

        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_weighted_comparison_with_standard_eval(
        self, spectral: Ncm.Spectral
    ) -> None:
        """Test that weighted coefficients correctly represent F(x)*sqrt(1-t^2).

        The weighted coefficients should represent the function F(x(t))*sqrt(1-t^2),
        not just F(x(t)). We can verify this by evaluating the Chebyshev expansion
        at various points.
        """

        def f(_user_data, x):
            return np.exp(-x * x)

        a, b = -1.0, 1.0
        k_min = 5
        tol = 1.0e-10

        _k, coeffs_w_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f, a, b, k_min, tol, None
        )
        coeffs_w = np.array(coeffs_w_list)

        # Evaluate the weighted Chebyshev expansion at several points
        test_t_values = np.array([0.0, 0.3, 0.5, 0.7, 0.9])

        for t in test_t_values:
            # Evaluate weighted Chebyshev expansion at t
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs_w, t)

            # Convert t to x in [a,b]
            x = Ncm.Spectral.t_to_x(a, b, t)

            # Expected value: F(x) * sqrt(1-t^2)
            expected = f(None, x) * np.sqrt(1.0 - t * t)

            # Allow some tolerance due to approximation
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-7,
                atol=1.0e-8,
                err_msg=f"Failed at t={t}",
            )

    def test_weighted_max_order_limit(self, spectral: Ncm.Spectral) -> None:
        """Test that weighted method respects max_order limit."""

        def f(_user_data, x):
            return np.exp(-x * x)

        a, b = -2.0, 2.0
        max_order = spectral.get_max_order()

        # Request very tight tolerance that might require exceeding max_order
        k_min = 2
        tol = 1.0e-14

        k_final, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # k_final should not exceed max_order
        assert k_final <= max_order, f"k_final={k_final} exceeds max_order={max_order}"

        # Integral should still be reasonable
        expected_integral, _error = quad(lambda x: np.exp(-x * x), a, b)
        computed_integral = coeffs[0] * np.pi
        assert_allclose(computed_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_weighted_zero_function(self, spectral: Ncm.Spectral) -> None:
        """Test weighted coefficients for zero function."""

        def f_zero(_user_data, _x):
            return 0.0

        a, b = -1.0, 1.0
        k_min = 2
        tol = 1.0e-12

        _k, coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_zero, a, b, k_min, tol, None
        )
        coeffs = np.array(coeffs_list)

        # All coefficients should be zero (or very close to zero)
        assert_allclose(coeffs, 0.0, rtol=1.0e-12, atol=1.0e-13)

    def test_weighted_scaled_function(self, spectral: Ncm.Spectral) -> None:
        """Test that scaling the function scales the integral linearly."""

        def f(_user_data, x):
            return np.exp(-x * x)

        a, b = -1.5, 1.5
        k_min = 4
        tol = 1.0e-10

        # Compute integral of f
        _k1, coeffs1_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f, a, b, k_min, tol, None
        )
        coeffs1 = np.array(coeffs1_list)

        # Compute integral of 3*f
        def f_scaled(_user_data, x):
            return 3.0 * np.exp(-x * x)

        _k2, coeffs2_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_scaled, a, b, k_min, tol, None
        )
        coeffs2 = np.array(coeffs2_list)

        # First coefficient of scaled function should be 3x the original
        # (both need the same (b-a)/2 factor, which cancels in the ratio)
        assert_allclose(coeffs2[0], 3.0 * coeffs1[0], rtol=1.0e-9, atol=1.0e-10)


class TestSpectralProductIntegral:
    """Tests for computing integrals of products.

    For two functions F(x) and G(x), the integral of their product can be computed as:
    ∫_a^b F(x)G(x)dx = π·a_0·b_0 + (π/2)·Σ_{k>0} a_k·b_k
    where a_k are weighted coefficients of F and b_k are standard coefficients of G.
    """

    @pytest.fixture(name="spectral")
    def fixture_spectral(self) -> Ncm.Spectral:
        """Create a spectral object."""
        return Ncm.Spectral.new()

    def _compute_product_integral_chebyshev(
        self, coeffs_weighted: np.ndarray, coeffs_standard: np.ndarray
    ) -> float:
        """Compute product integral using Chebyshev coefficients.

        Args:
            coeffs_weighted: Weighted Chebyshev coefficients (a_k)
            coeffs_standard: Standard Chebyshev coefficients (b_k)

        Returns:
            The integral ∫F(x)G(x)dx computed as π·a_0·b_0 + (π/2)·Σ_{k>0} a_k·b_k
        """
        # Ensure both arrays have the same length
        min_len = min(len(coeffs_weighted), len(coeffs_standard))
        a_coeffs = coeffs_weighted[:min_len]
        b_coeffs = coeffs_standard[:min_len]

        # First term: π·a_0·b_0
        result = np.pi * a_coeffs[0] * b_coeffs[0]

        # Sum term: (π/2)·Σ_{k>0} a_k·b_k
        if min_len > 1:
            result += (np.pi / 2.0) * np.sum(a_coeffs[1:] * b_coeffs[1:])

        return result

    def test_product_gaussian_polynomial(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: Gaussian × polynomial."""

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        def g_poly(_user_data, x):
            return 1.0 + x + 0.5 * x * x

        # Use interval that captures most of Gaussian support
        a, b = -3.0, 3.0
        k_min = 5
        tol = 1.0e-10

        # Compute weighted coefficients for Gaussian
        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_gaussian, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        # Compute standard coefficients for polynomial
        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_poly, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        # Compute integral using Chebyshev formula
        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        # Compute reference integral using scipy
        expected_integral, _error = quad(
            lambda x: f_gaussian(None, x) * g_poly(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_product_gaussian_exponential(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: Gaussian × exponential."""

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        def g_exp(_user_data, x):
            return np.exp(0.5 * x)

        a, b = -3.0, 3.0
        k_min = 5
        tol = 1.0e-10

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_gaussian, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_exp, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        expected_integral, _error = quad(
            lambda x: f_gaussian(None, x) * g_exp(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_product_gaussian_sine(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: Gaussian × sine."""

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        def g_sin(_user_data, x):
            return np.sin(2.0 * x)

        a, b = -3.0, 3.0
        k_min = 5
        tol = 1.0e-10

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_gaussian, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_sin, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        expected_integral, _error = quad(
            lambda x: f_gaussian(None, x) * g_sin(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_product_gaussian_cosine(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: Gaussian × cosine."""

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        def g_cos(_user_data, x):
            return np.cos(3.0 * x)

        a, b = -3.0, 3.0
        k_min = 5
        tol = 1.0e-10

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_gaussian, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_cos, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        expected_integral, _error = quad(
            lambda x: f_gaussian(None, x) * g_cos(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_product_gaussian_rational(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: Gaussian × rational function."""

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        def g_rational(_user_data, x):
            return 1.0 / (1.0 + 0.5 * x * x)

        a, b = -3.0, 3.0
        k_min = 5
        tol = 1.0e-10

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_gaussian, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_rational, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        expected_integral, _error = quad(
            lambda x: f_gaussian(None, x) * g_rational(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_product_gaussian_power(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: Gaussian × power function."""

        def f_gaussian(_user_data, x):
            return np.exp(-0.5 * x * x)

        def g_power(_user_data, x):
            return x * x * x * x  # x^4

        a, b = -2.5, 2.5
        k_min = 5
        tol = 1.0e-10

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_gaussian, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_power, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        expected_integral, _error = quad(
            lambda x: f_gaussian(None, x) * g_power(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_product_gaussian_oscillatory(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: Gaussian × oscillatory function."""

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        def g_osc(_user_data, x):
            return np.sin(5.0 * x) * np.cos(3.0 * x)

        a, b = -3.0, 3.0
        k_min = 6
        tol = 1.0e-10

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_gaussian, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_osc, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        expected_integral, _error = quad(
            lambda x: f_gaussian(None, x) * g_osc(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-7, atol=1.0e-8)

    def test_product_constant_functions(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: constant × constant."""

        def f_const(_user_data, _x):
            return 2.0

        def g_const(_user_data, _x):
            return 3.0

        a, b = -1.0, 1.0
        k_min = 2
        tol = 1.0e-12

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_const, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_const, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        # Expected: 2.0 * 3.0 * (b - a) = 12.0
        expected_integral = 2.0 * 3.0 * (b - a)

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_product_linear_quadratic(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: linear × quadratic."""

        def f_linear(_user_data, x):
            return 2.0 * x + 1.0

        def g_quad(_user_data, x):
            return x * x - 0.5

        a, b = -1.0, 1.0
        k_min = 3
        tol = 1.0e-12

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_linear, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_quad, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        # Analytical: (2x+1)(x^2-0.5) = 2x^3 + x^2 - x - 0.5
        # ∫_{-1}^{1} (2x^3 + x^2 - x - 0.5)dx = [x^4/2 + x^3/3 - x^2/2 - x/2]_{-1}^{1}
        # = (1/2 + 1/3 - 1/2 - 1/2) - (1/2 - 1/3 - 1/2 + 1/2)
        # = (1/3 - 1/2) - (1/2 - 1/3) = 2/3 - 1 = -1/3
        expected_integral, _error = quad(
            lambda x: f_linear(None, x) * g_quad(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)

    def test_product_symmetric_functions(self, spectral: Ncm.Spectral) -> None:
        """Test product integral: even function × odd function = 0."""

        def f_even(_user_data, x):
            return np.exp(-x * x)  # Even

        def g_odd(_user_data, x):
            return x * x * x  # x^3 is odd

        a, b = -2.0, 2.0
        k_min = 5
        tol = 1.0e-10

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_even, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_odd, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        # Integral of even × odd on symmetric interval should be ~0
        expected_integral, _error = quad(
            lambda x: f_even(None, x) * g_odd(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)
        assert abs(cheb_integral) < 1.0e-7  # Should be very close to zero

    def test_product_different_intervals(self, spectral: Ncm.Spectral) -> None:
        """Test product integral on various intervals."""

        def f_gaussian(_user_data, x):
            return np.exp(-((x - 1.0) ** 2))  # Shifted Gaussian

        def g_linear(_user_data, x):
            return 1.0 + 0.5 * x

        test_intervals = [(-2.0, 4.0), (0.0, 3.0), (-1.5, 2.5)]

        for a, b in test_intervals:
            k_min = 5
            tol = 1.0e-10

            _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
                f_gaussian, a, b, k_min, tol, None
            )
            a_coeffs = np.array(a_coeffs_list)

            _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
                g_linear, a, b, k_min, tol, None
            )
            b_coeffs = np.array(b_coeffs_list)

            cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

            expected_integral, _error = quad(
                lambda x: f_gaussian(None, x) * g_linear(None, x), a, b
            )

            assert_allclose(
                cheb_integral,
                expected_integral,
                rtol=1.0e-8,
                atol=1.0e-9,
                err_msg=f"Failed for interval [{a}, {b}]",
            )

    def test_product_narrow_gaussian(self, spectral: Ncm.Spectral) -> None:
        """Test product integral with narrow Gaussian."""

        def f_narrow_gaussian(_user_data, x):
            return np.exp(-10.0 * x * x)  # Narrow Gaussian

        def g_poly(_user_data, x):
            return 1.0 + x + x * x

        a, b = -2.0, 2.0
        k_min = 6
        tol = 1.0e-10

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_narrow_gaussian, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_poly, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        expected_integral, _error = quad(
            lambda x: f_narrow_gaussian(None, x) * g_poly(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-7, atol=1.0e-8)

    def test_product_wide_gaussian(self, spectral: Ncm.Spectral) -> None:
        """Test product integral with wide Gaussian (most support captured)."""

        def f_wide_gaussian(_user_data, x):
            return np.exp(-0.2 * x * x)  # Wide Gaussian

        def g_exp(_user_data, x):
            return np.exp(-0.5 * x)

        # Wide interval to capture most of the Gaussian support
        a, b = -5.0, 5.0
        k_min = 5
        tol = 1.0e-10

        _k1, a_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive_weighted(
            f_wide_gaussian, a, b, k_min, tol, None
        )
        a_coeffs = np.array(a_coeffs_list)

        _k2, b_coeffs_list = spectral.compute_chebyshev_coeffs_adaptive(
            g_exp, a, b, k_min, tol, None
        )
        b_coeffs = np.array(b_coeffs_list)

        cheb_integral = self._compute_product_integral_chebyshev(a_coeffs, b_coeffs)

        expected_integral, _error = quad(
            lambda x: f_wide_gaussian(None, x) * g_exp(None, x), a, b
        )

        assert_allclose(cheb_integral, expected_integral, rtol=1.0e-8, atol=1.0e-9)
