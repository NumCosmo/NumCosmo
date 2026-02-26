#
# test_py_spectral.py
#
# Tue Feb 04 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_spectral.py
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

"""Tests for NcmSpectral."""

import pytest
import numpy as np
from numpy.testing import assert_allclose
from numpy.polynomial.chebyshev import chebval
from scipy.special import iv, eval_gegenbauer  # pylint: disable=no-name-in-module
from scipy.fft import dct

from numcosmo_py import Ncm

ATOL = 1.0e-14


class TestSpectral:
    """Tests for NcmSpectral spectral methods."""

    @pytest.fixture(name="spectral")
    def fixture_spectral(self) -> Ncm.Spectral:
        """Create an ODE spectral with l=1."""
        return Ncm.Spectral.new()

    def test_chebT_to_gegenbauer_alpha1_basic(self) -> None:
        """Test Chebyshev to Gegenbauer alpha=1 conversion."""
        N = 16

        # T_0 = C_0^(1)
        c = np.zeros(N)
        c[0] = 1.0

        g = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha1(c))

        expected = np.zeros(N)
        expected[0] = 1.0
        assert_allclose(g, expected, rtol=1.0e-14, atol=ATOL)

        # T_1 = (1/2) C_1^(1)
        c = np.zeros(N)
        c[1] = 1.0

        g = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha1(c))

        expected = np.zeros(N)
        expected[1] = 0.5
        assert_allclose(g, expected, rtol=1.0e-14, atol=ATOL)

    def test_chebT_to_gegenbauer_alpha2_basic(self) -> None:
        """Test Chebyshev to Gegenbauer alpha=2 conversion."""
        N = 16

        # Test T_0 conversion
        c = np.zeros(N)
        c[0] = 1.0

        g = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(c))

        # For T_0, g_0 should have contribution 1/2 * c_0 + c_0/(2*1) = 1.0
        assert g[0] > 0.9

    def test_gegenbauer_alpha1_eval(self) -> None:
        """Test Gegenbauer alpha=1 evaluation."""
        N = 8

        # C_0^(1)(x) = 1, C_1^(1)(x) = 2x
        c = np.zeros(N)
        c[0] = 1.0
        c[1] = 1.0

        # At x=0: result should be 1 + 1*0 = 1
        result = Ncm.Spectral.gegenbauer_alpha1_eval(c, 0.0)
        assert_allclose(result, 1.0, rtol=1.0e-14)

        # At x=0.5: result should be 1 + 1*2*0.5 = 2
        result = Ncm.Spectral.gegenbauer_alpha1_eval(c, 0.5)
        assert_allclose(result, 2.0, rtol=1.0e-14)

        # At x=1: C_n^(1)(1) = n+1, so 1*1 + 1*2 = 3
        result = Ncm.Spectral.gegenbauer_alpha1_eval(c, 1.0)
        assert_allclose(result, 3.0, rtol=1.0e-14)

        # At x=-1: C_n^(1)(-1) = (n+1)*(-1)^n, so 1*1 + 1*2*(-1) = -1
        result = Ncm.Spectral.gegenbauer_alpha1_eval(c, -1.0)
        assert_allclose(result, -1.0, rtol=1.0e-14)

    def test_gegenbauer_alpha2_eval(self) -> None:
        """Test Gegenbauer alpha=2 evaluation."""
        N = 8

        # C_0^(2)(x) = 1, C_1^(2)(x) = 4x
        c = np.zeros(N)
        c[0] = 1.0
        c[1] = 1.0

        # At x=0: result should be 1 + 1*0 = 1
        result = Ncm.Spectral.gegenbauer_alpha2_eval(c, 0.0)
        assert_allclose(result, 1.0, rtol=1.0e-14)

        # At x=0.5: result should be 1 + 1*4*0.5 = 3
        result = Ncm.Spectral.gegenbauer_alpha2_eval(c, 0.5)
        assert_allclose(result, 3.0, rtol=1.0e-14)

    def test_chebyshev_eval_basic(self) -> None:
        """Test Chebyshev evaluation with simple polynomials."""
        N = 8

        # T_0(x) = 1
        a = np.zeros(N)
        a[0] = 1.0

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_eval(a, x)
            assert_allclose(result, 1.0, rtol=1.0e-14)

        # T_1(x) = x
        a = np.zeros(N)
        a[1] = 1.0

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_eval(a, x)
            assert_allclose(result, x, rtol=1.0e-14)

        # T_2(x) = 2x^2 - 1
        a = np.zeros(N)
        a[2] = 1.0

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_eval(a, x)
            expected = 2.0 * x * x - 1.0
            assert_allclose(result, expected, rtol=1.0e-14)

    def test_chebyshev_deriv_basic(self) -> None:
        """Test Chebyshev derivative with simple polynomials."""
        N = 8

        # Derivative of constant is 0
        a = np.zeros(N)
        a[0] = 1.0

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_deriv(a, x)
            assert_allclose(result, 0.0, rtol=1.0e-14, atol=ATOL)

        # Derivative of T_1(x) = x is 1
        a = np.zeros(N)
        a[1] = 1.0

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_deriv(a, x)
            assert_allclose(result, 1.0, rtol=1.0e-14)

        # Derivative of T_2(x) = 2x^2 - 1 is 4x
        a = np.zeros(N)
        a[2] = 1.0

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_deriv(a, x)
            expected = 4.0 * x
            assert_allclose(result, expected, rtol=1.0e-14)

    def test_spectral_serialization(self) -> None:
        """Test that NcmSpectral can be serialized and deserialized."""
        # Create a spectral with custom max_order
        max_order = 12
        spectral = Ncm.Spectral.new_with_max_order(max_order)

        # Serialize and deserialize
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        spectral_dup = ser.dup_obj(spectral)

        # Check that it's a different object
        assert spectral_dup is not None
        assert spectral_dup is not spectral
        assert isinstance(spectral_dup, Ncm.Spectral)

        # Check that max_order is preserved
        assert spectral.get_max_order() == max_order
        assert spectral_dup.get_max_order() == max_order

        # Test that both produce the same results with adaptive computation
        def test_func(_user_data, x):
            """Test function: Gaussian-like."""
            return np.exp(-x * x)

        a, b = -3.0, 3.0
        k_min = 2
        tol = 1.0e-10

        coeffs_orig = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                test_func, a, b, k_min, tol, None
            )
        )
        coeffs_dup = np.array(
            spectral_dup.compute_chebyshev_coeffs_adaptive(
                test_func, a, b, k_min, tol, None
            )
        )

        # Both should produce identical results
        assert len(coeffs_orig) == len(coeffs_dup)
        assert_allclose(coeffs_orig, coeffs_dup, rtol=1.0e-14, atol=ATOL)

        # Verify the approximation works for both
        test_points = np.linspace(a, b, 50)
        for x in test_points:
            val_orig = Ncm.Spectral.chebyshev_eval_x(coeffs_orig, a, b, x)
            val_dup = Ncm.Spectral.chebyshev_eval_x(coeffs_dup, a, b, x)
            expected = test_func(None, x)

            assert_allclose(val_orig, val_dup, rtol=1.0e-14, atol=ATOL)
            assert_allclose(val_orig, expected, rtol=1.0e-9, atol=1.0e-9)

    def test_compute_chebyshev_coeffs_constant(self, spectral: Ncm.Spectral) -> None:
        """Test computing Chebyshev coefficients for constant function."""
        N = 16

        def f_constant(_user_data, _x):
            return 1.0

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_constant, -1.0, 1.0, N, None)
        )

        # f(x) = 1 = T_0(x), so c_0 = 1, rest = 0
        assert_allclose(coeffs[0], 1.0, rtol=1.0e-12)
        assert_allclose(coeffs[1:], 0.0, rtol=1.0e-12, atol=1.0e-12)

    def test_compute_chebyshev_coeffs_linear(self, spectral: Ncm.Spectral) -> None:
        """Test computing Chebyshev coefficients for linear function."""
        N = 16

        def f_linear(_user_data, x):
            return x

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_linear, -1.0, 1.0, N, None)
        )

        # f(x) = x = T_1(x), so c_1 = 1, rest = 0
        assert_allclose(coeffs[0], 0.0, rtol=1.0e-12, atol=1.0e-12)
        assert_allclose(coeffs[1], 1.0, rtol=1.0e-12)
        assert_allclose(coeffs[2:], 0.0, rtol=1.0e-12, atol=1.0e-12)

    def test_compute_chebyshev_coeffs_quadratic(self, spectral: Ncm.Spectral) -> None:
        """Test computing Chebyshev coefficients for quadratic function."""
        N = 16

        def f_quadratic(_user_data, x):
            return x * x

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_quadratic, -1.0, 1.0, N, None)
        )

        # f(x) = x^2 = (T_0 + T_2)/2, so c_0 = 1/2, c_2 = 1/2
        assert_allclose(coeffs[0], 0.5, rtol=1.0e-11)
        assert_allclose(coeffs[1], 0.0, rtol=1.0e-12, atol=1.0e-12)
        assert_allclose(coeffs[2], 0.5, rtol=1.0e-11)
        assert_allclose(coeffs[3:], 0.0, rtol=1.0e-12, atol=1.0e-12)

    def test_compute_chebyshev_coeffs_gaussian(self, spectral: Ncm.Spectral) -> None:
        """Test computing Chebyshev coefficients for Gaussian function."""
        N = 64

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_gaussian, -1.0, 1.0, N, None)
        )

        # Verify by evaluating at test points
        test_points = np.array([0.0, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs, x)
            expected = np.exp(-x * x)
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-10,
                err_msg=f"Gaussian evaluation failed at x={x}",
            )

    def test_adaptive_vs_fixed_constant(self, spectral: Ncm.Spectral) -> None:
        """Test adaptive refinement vs fixed order for constant function."""

        def f_constant(_user_data, _x):
            return 1.0

        # Adaptive should converge quickly
        coeffs_adaptive = spectral.compute_chebyshev_coeffs_adaptive(
            f_constant, -1.0, 1.0, 2, 1e-12, None
        )
        N_adaptive = len(coeffs_adaptive)

        # Fixed order
        coeffs_fixed = np.array(
            spectral.compute_chebyshev_coeffs(f_constant, -1.0, 1.0, N_adaptive, None)
        )

        # Should match
        assert_allclose(
            np.array(coeffs_adaptive),
            coeffs_fixed,
            rtol=1e-12,
        )

    def test_adaptive_vs_fixed_polynomial(self, spectral: Ncm.Spectral) -> None:
        """Test adaptive refinement vs fixed order for polynomial."""

        def f_poly(_user_data, x):
            return 1.0 + 2.0 * x + 3.0 * x * x

        coeffs_adaptive = spectral.compute_chebyshev_coeffs_adaptive(
            f_poly, -1.0, 1.0, 3, 1e-12, None
        )
        N_adaptive = len(coeffs_adaptive)

        coeffs_fixed = np.array(
            spectral.compute_chebyshev_coeffs(f_poly, -1.0, 1.0, N_adaptive, None)
        )

        assert_allclose(
            np.array(coeffs_adaptive),
            coeffs_fixed,
            rtol=1e-11,
        )

    def test_adaptive_gaussian_convergence(self, spectral: Ncm.Spectral) -> None:
        """Test adaptive refinement converges for Gaussian."""

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        coeffs = spectral.compute_chebyshev_coeffs_adaptive(
            f_gaussian, -1.0, 1.0, 3, 1e-10, None
        )

        # Verify accuracy at test points
        test_points = np.linspace(-1.0, 1.0, 20)
        for x in test_points:
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs, x)
            expected = np.exp(-x * x)
            assert_allclose(eval_result, expected, rtol=1e-9, atol=1e-10)

    def test_adaptive_oscillatory_function(self, spectral: Ncm.Spectral) -> None:
        """Test adaptive refinement on oscillatory function."""

        def f_osc(_user_data, x):
            return np.sin(5.0 * x) * np.exp(-x * x)

        coeffs = spectral.compute_chebyshev_coeffs_adaptive(
            f_osc, -1.0, 1.0, 4, 1e-8, None
        )

        # Verify accuracy
        test_points = np.linspace(-1.0, 1.0, 30)
        for x in test_points:
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs, x)
            expected = np.sin(5.0 * x) * np.exp(-x * x)
            assert_allclose(eval_result, expected, rtol=1e-7, atol=1e-8)

    def test_adaptive_analytical_comparison(self, spectral: Ncm.Spectral) -> None:
        """Test adaptive coefficients match analytical solution for simple functions."""

        # Test 1: Constant function
        def f_const(_user_data, _x):
            return 2.5

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                f_const, -1.0, 1.0, 2, 1e-12, None
            )
        )
        assert_allclose(coeffs[0], 2.5, rtol=1e-12)
        assert_allclose(coeffs[1:], 0.0, rtol=1e-12, atol=1e-12)

        # Test 2: Linear function f(x) = 3x
        def f_linear(_user_data, x):
            return 3.0 * x

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                f_linear, -1.0, 1.0, 2, 1e-12, None
            )
        )
        assert_allclose(coeffs[0], 0.0, rtol=1e-12, atol=1e-12)
        assert_allclose(coeffs[1], 3.0, rtol=1e-12)
        assert_allclose(coeffs[2:], 0.0, rtol=1e-12, atol=1e-12)

    def test_adaptive_nested_nodes(self, spectral: Ncm.Spectral) -> None:
        """Test that adaptive refinement properly reuses nested nodes."""

        call_count = [0]

        def f_counted(_user_data, x):
            call_count[0] += 1
            return np.exp(-x * x)

        # Start at k=2 (N=5), refine to k=3 (N=9)
        # Should evaluate 5 nodes initially, then 4 new odd nodes
        coeffs = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                f_counted, -1.0, 1.0, 2, 1e-6, None
            )
        )

        N_final = coeffs.size
        # Should be less than 2*N_final due to node reuse
        assert call_count[0] < 2 * N_final

    def test_adaptive_interval_scaling(self, spectral: Ncm.Spectral) -> None:
        """Test adaptive refinement on different intervals."""

        def f_test(_user_data, x):
            return np.exp(-((x - 2.0) ** 2))

        # Test on [0, 4] instead of [-1, 1]
        coeffs = spectral.compute_chebyshev_coeffs_adaptive(
            f_test, 0.0, 4.0, 3, 1e-8, None
        )

        # Map test points from [0, 4] to [-1, 1] for Chebyshev evaluation
        test_points_orig = np.linspace(0.0, 4.0, 20)
        for x_orig in test_points_orig:
            # Map to [-1, 1]
            x_cheb = 2.0 * (x_orig - 0.0) / (4.0 - 0.0) - 1.0
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs, x_cheb)
            expected = np.exp(-((x_orig - 2.0) ** 2))
            assert_allclose(eval_result, expected, rtol=1e-7, atol=1e-8)

    def test_adaptive_max_order_limit(self) -> None:
        """Test that adaptive refinement respects max_order."""
        spectral = Ncm.Spectral.new_with_max_order(5)  # max N = 2^5+1 = 33

        def f_hard(_user_data, x):
            # Very oscillatory function that won't converge easily
            return np.sin(20.0 * x)

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                f_hard, -1.0, 1.0, 2, 1e-12, None
            )
        )

        # Should stop at max_order
        assert coeffs.size <= 33

    def test_adaptive_max_order_limit_higher(self) -> None:
        """Test that adaptive refinement respects max_order."""
        spectral = Ncm.Spectral.new_with_max_order(16)

        def f_hard(_user_data, x):
            # Very oscillatory function that won't converge easily
            return np.exp(-x * x) * np.sin(40.0 * x)

        def f_harder(_user_data, x):
            # Very oscillatory function that won't converge easily
            return np.exp(-x * x) * np.sin(600.0 * x)

        coeffs40 = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                f_hard, -1.0, 1.0, 2, 1e-12, None
            )
        )
        coeffs600 = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                f_harder, -1.0, 1.0, 2, 1e-12, None
            )
        )

        # Should stop at max_order
        assert coeffs40.size <= 2**16 + 1
        assert coeffs600.size <= 2**16 + 1
        # Should require more points for higher frequency
        assert coeffs40.size < coeffs600.size

    def test_adaptive_performance_simple(self) -> None:
        """Performance test for adaptive refinement with simple function."""
        spectral = Ncm.Spectral.new_with_max_order(16)

        def f_simple(_user_data, x):
            return np.exp(-x * x)

        for _ in range(100):
            coeffs = np.array(
                spectral.compute_chebyshev_coeffs_adaptive(
                    f_simple, -1.0, 1.0, 3, 1e-15, None
                )
            )
            assert coeffs.size > 0

    def test_adaptive_performance_complex(self) -> None:
        """Performance test for adaptive refinement with complex function."""
        spectral = Ncm.Spectral.new_with_max_order(16)

        def f_complex(_user_data, x):
            return np.exp(-x * x) * np.sin(10.0 * x) * np.cos(5.0 * x)

        for _ in range(100):
            coeffs = np.array(
                spectral.compute_chebyshev_coeffs_adaptive(
                    f_complex, -1.0, 1.0, 3, 1e-15, None
                )
            )
            assert coeffs.size > 0

    def test_adaptive_spectral_decay(self, spectral: Ncm.Spectral) -> None:
        """Test that converged adaptive solution shows spectral decay."""

        def f_smooth(_user_data, x):
            return 1.0 / (1.0 + 25.0 * x * x)  # Runge function

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                f_smooth, -1.0, 1.0, 4, 1e-8, None
            )
        )

        # Check that tail coefficients are small
        N = len(coeffs)
        tail_size = min(5, N // 10)
        max_tail = np.max(np.abs(coeffs[-tail_size:]))
        max_head = np.max(np.abs(coeffs[:tail_size]))

        # Tail should be much smaller than head
        assert max_tail < 1e-6 * max_head

    def test_adaptive_vs_fixed_accuracy(self, spectral: Ncm.Spectral) -> None:
        """Test adaptive achieves same accuracy as fixed with fewer/equal points."""

        def f_test(_user_data, x):
            return np.tanh(3.0 * x)

        # Adaptive with tolerance
        coeffs_adaptive = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                f_test, -1.0, 1.0, 3, 1e-10, None
            )
        )
        N_adaptive = len(coeffs_adaptive)

        # Fixed with same N
        coeffs_fixed = np.array(
            spectral.compute_chebyshev_coeffs(f_test, -1.0, 1.0, N_adaptive, None)
        )

        # Both should give similar accuracy
        test_points = np.linspace(-0.9, 0.9, 15)
        for x in test_points:
            result_adaptive = Ncm.Spectral.chebyshev_eval(coeffs_adaptive, x)
            result_fixed = Ncm.Spectral.chebyshev_eval(coeffs_fixed, x)
            expected = np.tanh(3.0 * x)

            assert_allclose(result_adaptive, expected, rtol=1e-9, atol=1e-14)
            assert_allclose(result_fixed, expected, rtol=1e-9, atol=1e-14)
            # Results should be very close to each other
            assert_allclose(result_adaptive, result_fixed, rtol=1e-11, atol=1e-14)

    def test_gaussian_derivative(self, spectral: Ncm.Spectral) -> None:
        """Test derivative of Gaussian function."""
        N = 64

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_gaussian, -1.0, 1.0, N, None)
        )

        # Derivative of exp(-x^2) is -2x*exp(-x^2)
        test_points = np.array([0.0, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            deriv_result = Ncm.Spectral.chebyshev_deriv(coeffs, x)
            expected = -2.0 * x * np.exp(-x * x)
            assert_allclose(
                deriv_result,
                expected,
                rtol=1.0e-9,
                atol=ATOL,
                err_msg=f"Gaussian derivative at x={x} doesn't match",
            )

    def test_rational_function(self, spectral: Ncm.Spectral) -> None:
        """Test with rational function x^2/(1+x^2)^4."""
        N = 64

        def f_rational(_user_data, x):
            return (x * x) / ((1.0 + x * x) ** 4)

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_rational, -1.0, 1.0, N, None)
        )

        # Verify by evaluating at test points
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs, x)
            expected = (x * x) / ((1.0 + x * x) ** 4)
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-10,
                atol=ATOL,
                err_msg=f"Rational function evaluation at x={x} doesn't match",
            )

    def test_rational_derivative(self, spectral: Ncm.Spectral) -> None:
        """Test derivative of rational function."""
        N = 64

        def f_rational(_user_data, x):
            return (x * x) / ((1.0 + x * x) ** 4)

        # Analytical derivative: f'(x) = 2x(1-3x^2)/(1+x^2)^5
        def f_prime(x):
            return 2.0 * x * (1.0 - 3.0 * x * x) / ((1.0 + x * x) ** 5)

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_rational, -1.0, 1.0, N, None)
        )
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            deriv_result = Ncm.Spectral.chebyshev_deriv(coeffs, x)
            expected = f_prime(x)
            assert_allclose(
                deriv_result,
                expected,
                rtol=1.0e-9,
                atol=ATOL,
                err_msg=f"Rational derivative at x={x} doesn't match",
            )

    def test_chebyshev_eval_endpoints(self, spectral: Ncm.Spectral) -> None:
        """Test that Chebyshev evaluation handles endpoints correctly."""
        N = 32

        # Use a simple smooth function
        def f_test(_user_data, x):
            return 1.0 + x + x * x

        coeffs = np.array(spectral.compute_chebyshev_coeffs(f_test, -1.0, 1.0, N, None))

        # Test at endpoints
        result_m1 = Ncm.Spectral.chebyshev_eval(coeffs, -1.0)
        expected_m1 = 1.0 - 1.0 + 1.0
        assert_allclose(result_m1, expected_m1, rtol=1.0e-12)

        result_p1 = Ncm.Spectral.chebyshev_eval(coeffs, 1.0)
        expected_p1 = 1.0 + 1.0 + 1.0
        assert_allclose(result_p1, expected_p1, rtol=1.0e-12)

    def test_chebyshev_deriv_endpoints(self, spectral: Ncm.Spectral) -> None:
        """Test that Chebyshev derivative handles endpoints correctly."""
        N = 32

        # Use f(x) = x^3, f'(x) = 3x^2
        def f_cubic(_user_data, x):
            return x * x * x

        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_cubic, -1.0, 1.0, N, None)
        )

        # Test derivative at endpoints
        result_m1 = Ncm.Spectral.chebyshev_deriv(coeffs, -1.0)
        expected_m1 = 3.0 * 1.0
        assert_allclose(result_m1, expected_m1, rtol=1.0e-11)

        result_p1 = Ncm.Spectral.chebyshev_deriv(coeffs, 1.0)
        expected_p1 = 3.0 * 1.0
        assert_allclose(result_p1, expected_p1, rtol=1.0e-11)

    def test_different_intervals(self, spectral: Ncm.Spectral) -> None:
        """Test computing Chebyshev coefficients on different intervals."""
        N = 32

        def f_test(_user_data, x):
            return x * x

        # Test on [0, 1]
        coeffs = np.array(spectral.compute_chebyshev_coeffs(f_test, 0.0, 1.0, N, None))

        # Map test points from [-1, 1] to [0, 1]
        for xi in [-0.5, 0.0, 0.5]:
            # xi in [-1, 1] maps to x in [0, 1]
            x = 0.5 + 0.5 * xi
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs, xi)
            expected = x * x
            assert_allclose(eval_result, expected, rtol=1.0e-11)

    def test_spectral_cache(self, spectral: Ncm.Spectral) -> None:
        """Test that spectral caching works correctly."""

        def f_test(_user_data, x):
            return np.sin(x)

        # Compute with N=32
        coeffs1 = np.array(
            spectral.compute_chebyshev_coeffs(f_test, -1.0, 1.0, 32, None)
        )

        # Compute again with N=32 (should use cache)
        coeffs2 = np.array(
            spectral.compute_chebyshev_coeffs(f_test, -1.0, 1.0, 32, None)
        )

        # Compute with different N=64 (should reallocate)
        coeffs3 = np.array(
            spectral.compute_chebyshev_coeffs(f_test, -1.0, 1.0, 64, None)
        )

        # Verify results are correct
        for coeffs, N in [(coeffs1, 32), (coeffs2, 32), (coeffs3, 64)]:
            assert coeffs.shape[0] == N
            # Test evaluation
            result = Ncm.Spectral.chebyshev_eval(coeffs, 0.5)
            expected = np.sin(0.5 * 0.5 * 2.0)  # Map to [-1,1]
            assert_allclose(result, expected, rtol=1.0e-10)

    @staticmethod
    def get_analytical_chebyshev_coeffs(power: int, N: int) -> np.ndarray:
        """
        Get analytical Chebyshev coefficients for x^power.

        Uses the recurrence relation for x^n in terms of Chebyshev polynomials.
        x^n can be expressed as a linear combination of T_k(x).
        """
        coeffs = np.zeros(N)

        if power == 0:
            # x^0 = 1 = T_0
            coeffs[0] = 1.0
        elif power == 1:
            # x = T_1
            if N > 1:
                coeffs[1] = 1.0
        elif power == 2:
            # x^2 = (T_0 + T_2) / 2
            coeffs[0] = 1.0 / 2.0
            if N > 2:
                coeffs[2] = 1.0 / 2.0
        elif power == 3:
            # x^3 = (3*T_1 + T_3) / 4
            if N > 1:
                coeffs[1] = 3.0 / 4.0
            if N > 3:
                coeffs[3] = 1.0 / 4.0
        elif power == 4:
            # x^4 = (3*T_0 + 4*T_2 + T_4) / 8
            coeffs[0] = 3.0 / 8.0
            if N > 2:
                coeffs[2] = 4.0 / 8.0
            if N > 4:
                coeffs[4] = 1.0 / 8.0
        elif power == 5:
            # x^5 = (10*T_1 + 5*T_3 + T_5) / 16
            if N > 1:
                coeffs[1] = 10.0 / 16.0
            if N > 3:
                coeffs[3] = 5.0 / 16.0
            if N > 5:
                coeffs[5] = 1.0 / 16.0
        elif power == 6:
            # x^6 = (10*T_0 + 15*T_2 + 6*T_4 + T_6) / 32
            coeffs[0] = 10.0 / 32.0
            if N > 2:
                coeffs[2] = 15.0 / 32.0
            if N > 4:
                coeffs[4] = 6.0 / 32.0
            if N > 6:
                coeffs[6] = 1.0 / 32.0
        elif power == 7:
            # x^7 = (35*T_1 + 21*T_3 + 7*T_5 + T_7) / 64
            if N > 1:
                coeffs[1] = 35.0 / 64.0
            if N > 3:
                coeffs[3] = 21.0 / 64.0
            if N > 5:
                coeffs[5] = 7.0 / 64.0
            if N > 7:
                coeffs[7] = 1.0 / 64.0
        elif power == 8:
            # x^8 = (35*T_0 + 56*T_2 + 28*T_4 + 8*T_6 + T_8) / 128
            coeffs[0] = 35.0 / 128.0
            if N > 2:
                coeffs[2] = 56.0 / 128.0
            if N > 4:
                coeffs[4] = 28.0 / 128.0
            if N > 6:
                coeffs[6] = 8.0 / 128.0
            if N > 8:
                coeffs[8] = 1.0 / 128.0
        elif power == 9:
            # x^9 = (126*T_1 + 84*T_3 + 36*T_5 + 9*T_7 + T_9) / 256
            if N > 1:
                coeffs[1] = 126.0 / 256.0
            if N > 3:
                coeffs[3] = 84.0 / 256.0
            if N > 5:
                coeffs[5] = 36.0 / 256.0
            if N > 7:
                coeffs[7] = 9.0 / 256.0
            if N > 9:
                coeffs[9] = 1.0 / 256.0
        elif power == 10:
            # x^10 = (126*T_0 + 210*T_2 + 120*T_4 + 45*T_6 + 10*T_8 + T_10) / 512
            coeffs[0] = 126.0 / 512.0
            if N > 2:
                coeffs[2] = 210.0 / 512.0
            if N > 4:
                coeffs[4] = 120.0 / 512.0
            if N > 6:
                coeffs[6] = 45.0 / 512.0
            if N > 8:
                coeffs[8] = 10.0 / 512.0
            if N > 10:
                coeffs[10] = 1.0 / 512.0

        return coeffs

    @pytest.mark.parametrize("power", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    @pytest.mark.parametrize("N", [16, 32, 64, 128])
    def test_chebyshev_coeffs_polynomials(
        self, spectral: Ncm.Spectral, power: int, N: int
    ) -> None:
        """
        Test Chebyshev coefficient computation for x^power.

        Computes coefficients using FFT and compares with analytical results.
        Tests various FFT sizes to ensure robustness.
        """

        # Define function f(x) = x^power
        def f_power(_user_data: None, x: float) -> float:
            return x**power

        # Compute Chebyshev coefficients
        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_power, -1.0, 1.0, N, None)
        )

        # Get analytical coefficients
        analytical = self.get_analytical_chebyshev_coeffs(power, N)

        # Compare - only non-zero analytical coefficients should match
        # (FFT might have numerical noise in other coefficients)
        for i in range(min(power + 1, N)):
            if abs(analytical[i]) > 1e-14:
                assert_allclose(
                    coeffs[i],
                    analytical[i],
                    rtol=1.0e-10,
                    atol=0.0,
                    err_msg=f"Mismatch at coefficient {i} for x^{power} with N={N}",
                )

        # Higher order coefficients should be near zero
        for i in range(power + 1, N):
            assert (
                abs(coeffs[i]) < 1e-10
            ), f"Coefficient {i} should be ~0 for x^{power}, got {coeffs[i]}"

    @staticmethod
    def get_analytical_exp_coeffs(N: int) -> np.ndarray:
        """
        Get analytical Chebyshev coefficients for exp(x) on [-1,1].

        The coefficients are: a_k = I_k(1) for k=0, and a_k = 2*I_k(1) for k>0,
        where I_k is the modified Bessel function of the first kind.
        """
        coeffs = np.zeros(N)
        coeffs[0] = iv(0, 1.0)
        for k in range(1, N):
            coeffs[k] = 2.0 * iv(k, 1.0)
        return coeffs

    @staticmethod
    def get_analytical_cos_coeffs(N: int) -> np.ndarray:
        """
        Get analytical Chebyshev coefficients for cos(x) on [-1,1].

        Since cos(x) is an even function, only even Chebyshev polynomials contribute.
        We use the identity: cos(x) = sum_{k=0}^infinity (-1)^k * x^(2k) / (2k)!
        and convert to Chebyshev basis numerically with high precision.
        """
        N_eval = max(N * 4, 256)  # Oversample for accuracy
        x = np.cos(np.pi * np.arange(N_eval) / (N_eval - 1))
        f = np.cos(x)

        coeffs_full = dct(f, type=1) / (N_eval - 1)
        coeffs_full[0] /= 2.0
        coeffs_full[-1] /= 2.0

        # Return truncated to N
        return coeffs_full[:N]

    @staticmethod
    def get_analytical_rational_coeffs(N: int) -> np.ndarray:
        """
        Get analytical Chebyshev coefficients for 1/(2-x) on [-1,1].

        This is a geometric series: 1/(2-x) = sum_{k=0}^infinity x^k/2^(k+1).
        The Chebyshev coefficients can be computed from the polynomial coefficients.
        For simplicity, we compute exactly for small N.
        """
        # 1/(2-x) = 1/2 * 1/(1-x/2) = 1/2 * sum (x/2)^k
        # This is more complex to convert exactly, so we'll use high-precision numerical
        # evaluation at Chebyshev points and DCT (same as the implementation)
        # For testing purposes, we'll compute numerically with high N

        N_eval = max(N * 4, 256)  # Oversample for accuracy
        x = np.cos(np.pi * np.arange(N_eval) / (N_eval - 1))
        f = 1.0 / (2.0 - x)

        coeffs_full = dct(f, type=1) / (N_eval - 1)
        coeffs_full[0] /= 2.0
        coeffs_full[-1] /= 2.0

        # Return truncated to N
        return coeffs_full[:N]

    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_chebyshev_coeffs_exp(self, spectral: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev coefficient computation for exp(x).

        Tests with exponential function which has known coefficients
        in terms of modified Bessel functions.
        """

        # Define function f(x) = exp(x)
        def f_exp(_user_data: None, x: float) -> float:
            return np.exp(x)

        # Compute Chebyshev coefficients
        coeffs = np.array(spectral.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, N, None))

        # Get analytical coefficients
        analytical = self.get_analytical_exp_coeffs(N)

        # Compare - exponential has exponentially decaying coefficients
        for i in range(min(20, N)):  # Check first 20 coefficients
            assert_allclose(
                coeffs[i],
                analytical[i],
                rtol=1.0e-10,
                atol=ATOL,
                err_msg=f"Mismatch at coefficient {i} for exp(x) with N={N}",
            )

    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_chebyshev_coeffs_cos(self, spectral: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev coefficient computation for cos(x).

        Tests with cosine function which has known coefficients
        in terms of Bessel functions.
        """

        # Define function f(x) = cos(x)
        def f_cos(_user_data: None, x: float) -> float:
            return np.cos(x)

        # Compute Chebyshev coefficients
        coeffs = np.array(spectral.compute_chebyshev_coeffs(f_cos, -1.0, 1.0, N, None))

        # Get analytical coefficients
        analytical = self.get_analytical_cos_coeffs(N)

        # Compare - cosine has exponentially decaying coefficients
        for i in range(min(20, N)):  # Check first 20 coefficients
            assert_allclose(
                coeffs[i],
                analytical[i],
                rtol=1.0e-10,
                atol=ATOL,
                err_msg=f"Mismatch at coefficient {i} for cos(x) with N={N}",
            )

    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_chebyshev_coeffs_rational(self, spectral: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev coefficient computation for 1/(2-x).

        Tests with rational function which has a singularity outside [-1,1].
        """

        # Define function f(x) = 1/(2-x)
        def f_rational(_user_data: None, x: float) -> float:
            return 1.0 / (2.0 - x)

        # Compute Chebyshev coefficients
        coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_rational, -1.0, 1.0, N, None)
        )

        # Get analytical coefficients (computed numerically with high precision)
        analytical = self.get_analytical_rational_coeffs(N)

        # Compare
        for i in range(min(20, N)):
            assert_allclose(
                coeffs[i],
                analytical[i],
                rtol=1.0e-9,
                atol=ATOL,
                err_msg=f"Mismatch at coefficient {i} for 1/(2-x) with N={N}",
            )

    @pytest.mark.parametrize("N", [16, 32, 64])
    def test_chebyshev_eval(self, N: int) -> None:
        """
        Test Chebyshev polynomial evaluation.

        Compares against numpy's chebval implementation.
        """
        # Generate random Chebyshev coefficients
        rng = np.random.default_rng(42)
        coeffs = rng.uniform(-1.0, 1.0, N)

        # Test at multiple points
        x_points = np.linspace(-1.0, 1.0, 1024)

        for x in x_points:
            # Evaluate using our implementation
            result = Ncm.Spectral.chebyshev_eval(coeffs, x)

            # Evaluate using numpy
            expected = chebval(x, coeffs)

            assert_allclose(
                result,
                expected,
                rtol=1.0e-12,
                atol=ATOL,
                err_msg=f"Chebyshev eval mismatch at x={x} with N={N}",
            )

    @pytest.mark.parametrize("power", [1, 2, 3, 4, 5])
    @pytest.mark.parametrize("N", [16, 32])
    def test_chebyshev_deriv_polynomials(self, power: int, N: int) -> None:
        """
        Test Chebyshev derivative evaluation for polynomials.

        For f(x) = x^power, the derivative is f'(x) = power * x^(power-1).
        We compute Chebyshev coefficients for x^power and test the derivative.
        """
        # Get Chebyshev coefficients for x^power
        coeffs = self.get_analytical_chebyshev_coeffs(power, N)

        # Test derivative at multiple points
        x_points = np.linspace(-0.95, 0.95, 19)

        for x in x_points:
            # Evaluate derivative using our implementation
            result = Ncm.Spectral.chebyshev_deriv(coeffs, x)

            # Expected derivative: d/dx(x^power) = power * x^(power-1)
            expected = power * x ** (power - 1) if power > 0 else 0.0

            assert_allclose(
                result,
                expected,
                rtol=1.0e-10,
                atol=ATOL,
                err_msg=f"Chebyshev deriv mismatch for x^{power} at x={x} with N={N}",
            )

    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_exp(self, spectral: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev derivative for exp(x).

        The derivative of exp(x) is exp(x) itself.
        """

        # Define function f(x) = exp(x)
        def f_exp(_user_data: None, x: float) -> float:
            return np.exp(x)

        # Compute Chebyshev coefficients
        coeffs = np.array(spectral.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, N, None))

        # Test derivative at multiple points
        x_points = np.linspace(-0.95, 0.95, 19)

        for x in x_points:
            # Evaluate derivative
            result = Ncm.Spectral.chebyshev_deriv(coeffs, x)

            # Expected: d/dx(exp(x)) = exp(x)
            expected = np.exp(x)

            assert_allclose(
                result,
                expected,
                rtol=1.0e-9,
                atol=ATOL,
                err_msg=f"Chebyshev deriv mismatch for exp(x) at x={x} with N={N}",
            )

    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_sin(self, spectral: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev derivative for sin(x).

        The derivative of sin(x) is cos(x).
        """

        # Define function f(x) = sin(x)
        def f_sin(_user_data: None, x: float) -> float:
            return np.sin(x)

        # Compute Chebyshev coefficients
        coeffs = np.array(spectral.compute_chebyshev_coeffs(f_sin, -1.0, 1.0, N, None))

        # Test derivative at multiple points
        x_points = np.linspace(-0.95, 0.95, 19)

        for x in x_points:
            # Evaluate derivative
            result = Ncm.Spectral.chebyshev_deriv(coeffs, x)

            # Expected: d/dx(sin(x)) = cos(x)
            expected = np.cos(x)

            assert_allclose(
                result,
                expected,
                rtol=1.0e-9,
                atol=ATOL,
                err_msg=f"Chebyshev deriv mismatch for sin(x) at x={x} with N={N}",
            )

    def test_chebyshev_deriv_edge_cases(self) -> None:
        """
        Test Chebyshev derivative edge cases.

        Tests single coefficients and constant functions.
        """
        # Test single coefficient (constant function)
        # f(x) = 5.0, so f'(x) = 0
        const = np.array([5.0])

        result = Ncm.Spectral.chebyshev_deriv(const, 0.5)
        assert result == 0.0, "Derivative of constant should be 0"

        result = Ncm.Spectral.chebyshev_deriv(const, 0.0)
        assert result == 0.0, "Derivative of constant should be 0 at x=0"

        result = Ncm.Spectral.chebyshev_deriv(const, 1.0)
        assert result == 0.0, "Derivative of constant should be 0 at x=1"

        result = Ncm.Spectral.chebyshev_deriv(const, -1.0)
        assert result == 0.0, "Derivative of constant should be 0 at x=-1"

        # Test two coefficients (linear function: a0 + a1*x)
        # f(x) = 2.0 + 3.0*x, so f'(x) = 3.0
        linear = np.array([2.0, 3.0])

        x_test = 0.5
        result = Ncm.Spectral.chebyshev_deriv(linear, x_test)
        assert_allclose(result, 3.0, rtol=1.0e-14, atol=ATOL)

        # Test at x=0
        result = Ncm.Spectral.chebyshev_deriv(linear, 0.0)
        assert_allclose(result, 3.0, rtol=1.0e-14, atol=ATOL)

        # Test at x=1
        result = Ncm.Spectral.chebyshev_deriv(linear, 1.0)
        assert_allclose(result, 3.0, rtol=1.0e-14, atol=ATOL)

        # Test at x=-1
        result = Ncm.Spectral.chebyshev_deriv(linear, -1.0)
        assert_allclose(result, 3.0, rtol=1.0e-14, atol=ATOL)

    @pytest.mark.parametrize("power", [0, 1, 2, 3, 4, 5])
    @pytest.mark.parametrize("N", [32])
    def test_chebyshev_deriv_boundary(self, power: int, N: int) -> None:
        """
        Test Chebyshev derivative at boundary points x=+/-1.

        Verifies that derivatives evaluate correctly at the endpoints.
        For f(x) = x^power, f'(x) = power * x^(power-1).
        """
        # Get Chebyshev coefficients for x^power
        coeffs = self.get_analytical_chebyshev_coeffs(power, N)
        coeffs = np.array(coeffs)

        # Test at x = 1: f'(1) = power * 1^(power-1) = power
        result_p1 = Ncm.Spectral.chebyshev_deriv(coeffs, 1.0)
        expected_p1 = power if power > 0 else 0.0
        assert_allclose(
            result_p1,
            expected_p1,
            rtol=1.0e-10,
            atol=ATOL,
            err_msg=f"Derivative of x^{power} at x=1 failed",
        )

        # Test at x = -1: f'(-1) = power * (-1)^(power-1)
        result_m1 = Ncm.Spectral.chebyshev_deriv(coeffs, -1.0)
        expected_m1 = power * ((-1.0) ** (power - 1)) if power > 0 else 0.0
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=1.0e-10,
            atol=ATOL,
            err_msg=f"Derivative of x^{power} at x=-1 failed",
        )

    @pytest.mark.parametrize("n", [1, 2, 5, 10, 20])
    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_single_mode(self, n: int, N: int) -> None:
        """
        Test derivative of a single Chebyshev mode T_n(x).

        f(x) = T_n(x)
        f'(x) = n * U_{n-1}(x)
        """
        coeffs = np.zeros(N)
        coeffs[n] = 1.0

        x_points = np.linspace(-0.9, 0.9, 17)

        for x in x_points:
            result = Ncm.Spectral.chebyshev_deriv(coeffs, x)

            theta = np.arccos(x)
            expected = n * np.sin(n * theta) / np.sin(theta)

            assert_allclose(
                result,
                expected,
                rtol=1.0e-10,
                atol=1.0e-13,
                err_msg=f"Derivative of T_{n} failed at x={x}",
            )

    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_alternating_coeffs(self, N: int) -> None:
        """
        Test derivative with alternating high-frequency coefficients.
        """
        coeffs = np.array([(-1.0) ** k for k in range(N)], dtype=float)

        x_points = np.linspace(-0.8, 0.8, 11)
        h = 1.0e-7

        for x in x_points:
            d_numeric = (
                Ncm.Spectral.chebyshev_eval(coeffs, x + h)
                - Ncm.Spectral.chebyshev_eval(coeffs, x - h)
            ) / (2.0 * h)

            d_cheb = Ncm.Spectral.chebyshev_deriv(coeffs, x)

            assert_allclose(
                d_cheb,
                d_numeric,
                rtol=1.0e-6,
                atol=ATOL,
                err_msg=f"Alternating coeff derivative mismatch at x={x}",
            )

    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_random_coeffs(self, N: int) -> None:
        """
        Test derivative for random Chebyshev coefficients using finite differences.
        """
        rng = np.random.default_rng(1234)
        coeffs = rng.normal(size=N)

        x_points = np.linspace(-0.7, 0.7, 9)
        h = 1.0e-7

        for x in x_points:
            f_plus = Ncm.Spectral.chebyshev_eval(coeffs, x + h)
            f_minus = Ncm.Spectral.chebyshev_eval(coeffs, x - h)
            d_numeric = (f_plus - f_minus) / (2.0 * h)

            d_cheb = Ncm.Spectral.chebyshev_deriv(coeffs, x)

            assert_allclose(
                d_cheb,
                d_numeric,
                rtol=1.0e-6,
                atol=ATOL,
                err_msg=f"Random coeff derivative mismatch at x={x}",
            )

    @pytest.mark.parametrize("n", [1, 2, 3, 5, 10, 20])
    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_Tn_endpoints(self, n: int, N: int) -> None:
        """
        Exact endpoint test for the derivative of T_n(x).

        T_n'(1)  = n^2
        T_n'(-1) = (-1)^(n-1) * n^2
        """
        coeffs = np.zeros(N)
        coeffs[n] = 1.0

        # x = 1
        result_p1 = Ncm.Spectral.chebyshev_deriv(coeffs, 1.0)
        assert_allclose(
            result_p1,
            n * n,
            rtol=0.0,
            atol=ATOL,
            err_msg=f"T_{n}'(1) failed",
        )

        # x = -1
        result_m1 = Ncm.Spectral.chebyshev_deriv(coeffs, -1.0)
        expected_m1 = ((-1.0) ** (n - 1)) * n * n
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=0.0,
            atol=ATOL,
            err_msg=f"T_{n}'(-1) failed",
        )

    @pytest.mark.parametrize("power", [1, 2, 3, 4, 5])
    @pytest.mark.parametrize("N", [32])
    def test_chebT_to_gegenbauer_alpha1(self, power: int, N: int) -> None:
        """
        Test conversion from Chebyshev to Gegenbauer (alpha=1) basis.

        For a polynomial x^power, we convert to Chebyshev basis, then to
        Gegenbauer basis, and verify the polynomial evaluates correctly.
        """
        # Get Chebyshev coefficients for x^power
        cheb_coeffs = self.get_analytical_chebyshev_coeffs(power, N)

        # Convert to Gegenbauer alpha=1 (U_n Chebyshev of second kind)
        gegen = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha1(cheb_coeffs))

        # Test evaluation at multiple points
        x_points = np.linspace(-0.9, 0.9, 11)

        for x in x_points:
            # Evaluate using Chebyshev basis
            cheb_val = Ncm.Spectral.chebyshev_eval(cheb_coeffs, x)

            # Evaluate using Gegenbauer alpha=1 basis
            gegen_val = Ncm.Spectral.gegenbauer_alpha1_eval(gegen, x)

            # Both should give x^power
            expected = x**power

            assert_allclose(
                cheb_val,
                expected,
                rtol=1.0e-10,
                atol=ATOL,
                err_msg=f"Chebyshev eval failed for x^{power} at x={x}",
            )

            assert_allclose(
                gegen_val,
                expected,
                rtol=1.0e-10,
                atol=ATOL,
                err_msg=f"Gegenbauer alpha=1 eval failed for x^{power} at x={x}",
            )

    @pytest.mark.parametrize("power", [1, 2, 3, 4, 5])
    @pytest.mark.parametrize("N", [32])
    def test_chebT_to_gegenbauer_alpha2(self, power: int, N: int) -> None:
        """
        Test conversion from Chebyshev to Gegenbauer (alpha=2) basis.

        For a polynomial x^power, we convert to Chebyshev basis, then to
        Gegenbauer basis, and verify the polynomial evaluates correctly.
        """
        # Get Chebyshev coefficients for x^power
        cheb_coeffs = self.get_analytical_chebyshev_coeffs(power, N)

        # Convert to Gegenbauer alpha=2
        gegen_coeffs = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_coeffs))

        # Test evaluation at multiple points
        x_points = np.linspace(-0.9, 0.9, 11)

        for x in x_points:
            # Evaluate using Chebyshev basis
            cheb_val = Ncm.Spectral.chebyshev_eval(cheb_coeffs, x)

            # Evaluate using Gegenbauer alpha=2 basis
            gegen_val = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_coeffs, x)

            # Both should give x^power
            expected = x**power

            assert_allclose(
                cheb_val,
                expected,
                rtol=1.0e-10,
                atol=ATOL,
                err_msg=f"Chebyshev eval failed for x^{power} at x={x}",
            )

            assert_allclose(
                gegen_val,
                expected,
                rtol=1.0e-10,
                atol=ATOL,
                err_msg=f"Gegenbauer alpha=2 eval failed for x^{power} at x={x}",
            )

    @pytest.mark.parametrize("N", [16, 32])
    @pytest.mark.parametrize("alpha_val", [1, 2])
    def test_gegenbauer_eval_against_scipy(self, N: int, alpha_val: int) -> None:
        """
        Test Gegenbauer polynomial evaluation against scipy.

        Tests that our evaluation matches scipy's eval_gegenbauer for
        individual Gegenbauer polynomials.
        """
        x_points = np.linspace(-0.95, 0.95, 11)

        # Test individual Gegenbauer polynomials C_n^(alpha)
        for n in range(min(10, N)):
            # Create coefficients for pure C_n^(alpha) polynomial
            coeffs = np.zeros(N)
            coeffs[n] = 1.0

            for x in x_points:
                # Evaluate using our implementation
                if alpha_val == 1:
                    result = Ncm.Spectral.gegenbauer_alpha1_eval(coeffs, x)
                else:
                    result = Ncm.Spectral.gegenbauer_alpha2_eval(coeffs, x)

                # Evaluate using scipy
                expected = eval_gegenbauer(n, alpha_val, x)

                assert_allclose(
                    result,
                    expected,
                    rtol=1.0e-10,
                    atol=ATOL,
                    err_msg=f"Gegenbauer C_{n}^({alpha_val}) mismatch at x={x}",
                )

    @pytest.mark.parametrize("N", [32, 64])
    def test_gegenbauer_roundtrip(self, spectral: Ncm.Spectral, N: int) -> None:
        """
        Test round-trip conversion: function -> Chebyshev -> Gegenbauer -> evaluation.

        Uses exp(x) as test function to verify the entire pipeline.
        """

        def f_exp(_user_data: None, x: float) -> float:
            return np.exp(x)

        # Compute Chebyshev coefficients
        coeffs = np.array(spectral.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, N, None))

        # Convert to Gegenbauer alpha=1

        gegen1_coeffs = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha1(coeffs))

        # Convert to Gegenbauer alpha=2
        gegen2_coeffs = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(coeffs))
        # Test evaluation at multiple points
        x_points = np.linspace(-0.8, 0.8, 17)

        for x in x_points:
            expected = np.exp(x)

            # Evaluate with Chebyshev
            cheb_val = Ncm.Spectral.chebyshev_eval(coeffs, x)

            # Evaluate with Gegenbauer alpha=1
            gegen1_val = Ncm.Spectral.gegenbauer_alpha1_eval(gegen1_coeffs, x)

            # Evaluate with Gegenbauer alpha=2
            gegen2_val = Ncm.Spectral.gegenbauer_alpha2_eval(gegen2_coeffs, x)

            # All should match exp(x)
            assert_allclose(
                cheb_val,
                expected,
                rtol=1.0e-8,
                atol=ATOL,
                err_msg=f"Chebyshev roundtrip failed for exp(x) at x={x}",
            )

            assert_allclose(
                gegen1_val,
                expected,
                rtol=1.0e-8,
                atol=ATOL,
                err_msg=f"Gegenbauer alpha=1 roundtrip failed for exp(x) at x={x}",
            )

            assert_allclose(
                gegen2_val,
                expected,
                rtol=1.0e-8,
                atol=ATOL,
                err_msg=f"Gegenbauer alpha=2 roundtrip failed for exp(x) at x={x}",
            )

    @pytest.mark.parametrize("power", [0, 1, 2, 3, 4, 5])
    @pytest.mark.parametrize("N", [32])
    def test_chebyshev_eval_boundary(self, power: int, N: int) -> None:
        """
        Test Chebyshev evaluation at boundary points x=+/-1.

        Verifies that x^power evaluates correctly at the endpoints.
        """
        # Get Chebyshev coefficients for x^power
        cheb_coeffs = self.get_analytical_chebyshev_coeffs(power, N)

        # Test at x = 1
        result_p1 = Ncm.Spectral.chebyshev_eval(cheb_coeffs, 1.0)
        expected_p1 = 1.0**power
        assert_allclose(
            result_p1,
            expected_p1,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg=f"Chebyshev eval of x^{power} at x=1 failed",
        )

        # Test at x = -1
        result_m1 = Ncm.Spectral.chebyshev_eval(cheb_coeffs, -1.0)
        expected_m1 = (-1.0) ** power
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg=f"Chebyshev eval of x^{power} at x=-1 failed",
        )

    @pytest.mark.parametrize("power", [0, 1, 2, 3, 4, 5])
    @pytest.mark.parametrize("N", [32])
    def test_gegenbauer_alpha1_eval_boundary(self, power: int, N: int) -> None:
        """
        Test Gegenbauer alpha=1 evaluation at boundary points x=+/-1.

        Verifies that x^power evaluates correctly at the endpoints.
        """
        # Get Chebyshev coefficients for x^power
        cheb_coeffs = self.get_analytical_chebyshev_coeffs(power, N)

        # Convert to Gegenbauer alpha=1
        gegen_coeffs = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha1(cheb_coeffs))

        # Test at x = 1
        result_p1 = Ncm.Spectral.gegenbauer_alpha1_eval(gegen_coeffs, 1.0)
        expected_p1 = 1.0**power
        assert_allclose(
            result_p1,
            expected_p1,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg=f"Gegenbauer alpha=1 eval of x^{power} at x=1 failed",
        )

        # Test at x = -1
        result_m1 = Ncm.Spectral.gegenbauer_alpha1_eval(gegen_coeffs, -1.0)
        expected_m1 = (-1.0) ** power
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg=f"Gegenbauer alpha=1 eval of x^{power} at x=-1 failed",
        )

    @pytest.mark.parametrize("power", [0, 1, 2, 3, 4, 5])
    @pytest.mark.parametrize("N", [32])
    def test_gegenbauer_alpha2_eval_boundary(self, power: int, N: int) -> None:
        """
        Test Gegenbauer alpha=2 evaluation at boundary points x=+/-1.

        Verifies that x^power evaluates correctly at the endpoints.
        """
        # Get Chebyshev coefficients for x^power
        cheb_coeffs = self.get_analytical_chebyshev_coeffs(power, N)

        # Convert to Gegenbauer alpha=2
        gegen_coeffs = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_coeffs))

        # Test at x = 1
        result_p1 = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_coeffs, 1.0)
        expected_p1 = 1.0**power
        assert_allclose(
            result_p1,
            expected_p1,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg=f"Gegenbauer alpha=2 eval of x^{power} at x=1 failed",
        )

        # Test at x = -1
        result_m1 = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_coeffs, -1.0)
        expected_m1 = (-1.0) ** power
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg=f"Gegenbauer alpha=2 eval of x^{power} at x=-1 failed",
        )

    def test_proj_matrix_identity(self) -> None:
        """Test that projection matrix correctly converts x^2 to Gegenbauer."""
        N = 32

        # Get Chebyshev coefficients for x^2
        cheb_coeffs = self.get_analytical_chebyshev_coeffs(2, N)

        # Get projection matrix
        proj_mat = Ncm.Spectral.get_proj_matrix(N)
        proj_np = proj_mat.to_numpy()

        # Apply matrix: result should be Gegenbauer coeffs of x^2
        gegen_coeffs = proj_np @ cheb_coeffs

        # Verify by converting using the function
        gegen_coeffs_expected = np.array(
            Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_coeffs)
        )

        assert_allclose(
            gegen_coeffs,
            gegen_coeffs_expected,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg="Projection matrix doesn't match conversion function",
        )

    def test_x_operator_on_x(self) -> None:
        """Test that x operator applied to x gives x^2."""
        N = 32

        # Get Chebyshev coefficients for x (input)
        cheb_x = self.get_analytical_chebyshev_coeffs(1, N)

        # Get x operator matrix
        x_mat = Ncm.Spectral.get_x_matrix(N)
        x_np = x_mat.to_numpy()

        # Apply: x * x = x^2
        gegen_result = x_np @ cheb_x

        # Expected: Gegenbauer coefficients of x^2
        cheb_x2 = self.get_analytical_chebyshev_coeffs(2, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_x2))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg="x operator on x doesn't give x^2",
        )

    def test_x_operator_on_x2(self) -> None:
        """Test that x operator applied to x^2 gives x^3."""
        N = 32

        # Get Chebyshev coefficients for x^2 (input)
        cheb_x2 = self.get_analytical_chebyshev_coeffs(2, N)

        # Get x operator matrix
        x_mat = Ncm.Spectral.get_x_matrix(N)
        x_np = x_mat.to_numpy()

        # Apply: x * x^2 = x^3
        gegen_result = x_np @ cheb_x2

        # Expected: Gegenbauer coefficients of x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_x3))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg="x operator on x^2 doesn't give x^3",
        )

    def test_x2_operator_on_x(self) -> None:
        """Test that x^2 operator applied to x gives x^3."""
        N = 32

        # Get Chebyshev coefficients for x (input)
        cheb_x = self.get_analytical_chebyshev_coeffs(1, N)

        # Get x^2 operator matrix
        x2_mat = Ncm.Spectral.get_x2_matrix(N)
        x2_np = x2_mat.to_numpy()

        # Apply: x^2 * x = x^3
        gegen_result = x2_np @ cheb_x

        # Expected: Gegenbauer coefficients of x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_x3))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg="x^2 operator on x doesn't give x^3",
        )

    def test_x2_operator_on_x2(self) -> None:
        """Test that x^2 operator applied to x^2 gives x^4."""
        N = 32

        # Get Chebyshev coefficients for x^2 (input)
        cheb_x2 = self.get_analytical_chebyshev_coeffs(2, N)

        # Get x^2 operator matrix
        x2_mat = Ncm.Spectral.get_x2_matrix(N)
        x2_np = x2_mat.to_numpy()

        # Apply: x^2 * x^2 = x^4
        gegen_result = x2_np @ cheb_x2

        # Expected: Gegenbauer coefficients of x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_x4))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg="x^2 operator on x^2 doesn't give x^4",
        )

    def test_d_operator_on_x2(self) -> None:
        """Test that derivative operator applied to x^2 gives 2x."""
        N = 32

        # Get Chebyshev coefficients for x^2
        cheb_x2 = self.get_analytical_chebyshev_coeffs(2, N)

        # Get derivative operator matrix
        d_mat = Ncm.Spectral.get_d_matrix(N)
        d_np = d_mat.to_numpy()

        # Apply: d/dx(x^2) = 2x
        gegen_result = d_np @ cheb_x2

        # Expected: Gegenbauer coefficients of 2x
        cheb_2x = 2.0 * self.get_analytical_chebyshev_coeffs(1, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_2x))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg="Derivative of x^2 doesn't give 2x",
        )

    def test_d_operator_on_x3(self) -> None:
        """Test that derivative operator applied to x^3 gives 3x^2."""
        N = 32

        # Get Chebyshev coefficients for x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Get derivative operator matrix
        d_mat = Ncm.Spectral.get_d_matrix(N)
        d_np = d_mat.to_numpy()

        # Apply: d/dx(x^3) = 3x^2
        gegen_result = d_np @ cheb_x3

        # Expected: Gegenbauer coefficients of 3x^2
        cheb_3x2 = 3.0 * self.get_analytical_chebyshev_coeffs(2, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_3x2))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg="Derivative of x^3 doesn't give 3x^2",
        )

    def test_d2_operator_on_x3(self) -> None:
        """Test that second derivative operator applied to x^3 gives 6x."""
        N = 32

        # Get Chebyshev coefficients for x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Get second derivative operator matrix
        d2_mat = Ncm.Spectral.get_d2_matrix(N)
        d2_np = d2_mat.to_numpy()

        # Apply: d^2/dx^2(x^3) = 6x
        gegen_result = d2_np @ cheb_x3

        # Expected: Gegenbauer coefficients of 6x
        cheb_6x = 6.0 * self.get_analytical_chebyshev_coeffs(1, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_6x))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg="Second derivative of x^3 doesn't give 6x",
        )

    def test_d2_operator_on_x4(self) -> None:
        """Test that second derivative operator applied to x^4 gives 12x^2."""
        N = 32

        # Get Chebyshev coefficients for x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Get second derivative operator matrix
        d2_mat = Ncm.Spectral.get_d2_matrix(N)
        d2_np = d2_mat.to_numpy()

        # Apply: d^2/dx^2(x^4) = 12x^2
        gegen_result = d2_np @ cheb_x4

        # Expected: Gegenbauer coefficients of 12x^2
        cheb_12x2 = 12.0 * self.get_analytical_chebyshev_coeffs(2, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_12x2))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=ATOL,
            err_msg="Second derivative of x^4 doesn't give 12x^2",
        )

    def test_x_d_operator_on_x3(self) -> None:
        """Test that x*d/dx operator applied to x^3 gives 3x^3."""
        N = 32

        # Get Chebyshev coefficients for x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Get x*d operator matrix
        x_d_mat = Ncm.Spectral.get_x_d_matrix(N)
        x_d_np = x_d_mat.to_numpy()

        # Apply: x * d/dx(x^3) = x * 3x^2 = 3x^3
        gegen_result = x_d_np @ cheb_x3

        # Expected: Gegenbauer coefficients of 3x^3
        cheb_3x3 = 3.0 * self.get_analytical_chebyshev_coeffs(3, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_3x3))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-11,
            atol=ATOL,
            err_msg="x*d operator on x^3 doesn't give 3x^3",
        )

    def test_x_d2_operator_on_x4(self) -> None:
        """Test that x*d^2/dx^2 operator applied to x^4 gives 12x^3."""
        N = 32

        # Get Chebyshev coefficients for x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Get x*d^2 operator matrix
        x_d2_mat = Ncm.Spectral.get_x_d2_matrix(N)
        x_d2_np = x_d2_mat.to_numpy()

        # Apply: x * d^2/dx^2(x^4) = x * 12x^2 = 12x^3
        gegen_result = x_d2_np @ cheb_x4

        # Expected: Gegenbauer coefficients of 12x^3
        cheb_12x3 = 12.0 * self.get_analytical_chebyshev_coeffs(3, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_12x3))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-11,
            atol=ATOL,
            err_msg="x*d^2 operator on x^4 doesn't give 12x^3",
        )

    def test_x2_d2_operator_on_x4(self) -> None:
        """Test that x^2*d^2/dx^2 operator applied to x^4 gives 12x^4."""
        N = 32

        # Get Chebyshev coefficients for x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Get x^2*d^2 operator matrix
        x2_d2_mat = Ncm.Spectral.get_x2_d2_matrix(N)
        x2_d2_np = x2_d2_mat.to_numpy()

        # Apply: x^2 * d^2/dx^2(x^4) = x^2 * 12x^2 = 12x^4
        gegen_result = x2_d2_np @ cheb_x4

        # Expected: Gegenbauer coefficients of 12x^4
        cheb_12x4 = 12.0 * self.get_analytical_chebyshev_coeffs(4, N)
        gegen_expected = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_12x4))

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-11,
            atol=ATOL,
            err_msg="x^2*d^2 operator on x^4 doesn't give 12x^4",
        )

    def test_d_operator_evaluation(self) -> None:
        """Test d/dx operator by evaluating result at points."""
        N = 32

        # Start with x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Apply d operator: d/dx(x^3) = 3x^2
        d_mat = Ncm.Spectral.get_d_matrix(N)
        d_np = d_mat.to_numpy()
        gegen_result = d_np @ cheb_x3

        # Evaluate at test points and compare with 3*x^2
        test_points = np.array([0.0, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            eval_result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_result, x)
            expected = 3.0 * x**2
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-12,
                atol=ATOL,
                err_msg=f"d/dx(x^3) evaluation at x={x} doesn't match 3x^2",
            )

    def test_d2_operator_evaluation(self) -> None:
        """Test d^2/dx^2 operator by evaluating result at points."""
        N = 32

        # Start with x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Apply d^2 operator: d^2/dx^2(x^4) = 12x^2
        d2_mat = Ncm.Spectral.get_d2_matrix(N)
        d2_np = d2_mat.to_numpy()
        gegen_result = d2_np @ cheb_x4

        # Evaluate at test points and compare with 12*x^2
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            eval_result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_result, x)
            expected = 12.0 * x**2
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-12,
                atol=ATOL,
                err_msg=f"d^2/dx^2(x^4) evaluation at x={x} doesn't match 12x^2",
            )

    def test_x_d_operator_evaluation(self) -> None:
        """Test x*d/dx operator by evaluating result at points."""
        N = 32

        # Start with x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Apply x*d operator: x*d/dx(x^4) = x*4x^3 = 4x^4
        x_d_mat = Ncm.Spectral.get_x_d_matrix(N)
        x_d_np = x_d_mat.to_numpy()
        gegen_result = x_d_np @ cheb_x4

        # Evaluate at test points and compare with 4*x^4
        test_points = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            eval_result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_result, x)
            expected = 4.0 * x**4
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-11,
                atol=ATOL,
                err_msg=f"x*d/dx(x^4) evaluation at x={x} doesn't match 4x^4",
            )

    def test_x2_d2_operator_evaluation(self) -> None:
        """Test x^2*d^2/dx^2 operator by evaluating result at points."""
        N = 32

        # Start with x^5
        cheb_x5 = self.get_analytical_chebyshev_coeffs(5, N)

        # Apply x^2*d^2 operator: x^2*d^2/dx^2(x^5) = x^2*20x^3 = 20x^5
        x2_d2_mat = Ncm.Spectral.get_x2_d2_matrix(N)
        x2_d2_np = x2_d2_mat.to_numpy()
        gegen_result = x2_d2_np @ cheb_x5

        # Evaluate at test points and compare with 20*x^5
        test_points = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            eval_result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_result, x)
            expected = 20.0 * x**5
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-11,
                atol=ATOL,
                err_msg=f"x^2*d^2/dx^2(x^5) evaluation at x={x} doesn't match 20x^5",
            )

    def test_x_operator_evaluation(self) -> None:
        """Test x operator by evaluating result at points."""
        N = 32

        # Start with x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Apply x operator: x*x^3 = x^4
        x_mat = Ncm.Spectral.get_x_matrix(N)
        x_np = x_mat.to_numpy()
        gegen_result = x_np @ cheb_x3

        # Evaluate at test points and compare with x^4
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            eval_result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_result, x)
            expected = x**4
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-12,
                atol=ATOL,
                err_msg=f"x*x^3 evaluation at x={x} doesn't match x^4",
            )

    def test_x2_operator_evaluation(self) -> None:
        """Test x^2 operator by evaluating result at points."""
        N = 32

        # Start with x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Apply x^2 operator: x^2*x^3 = x^5
        x2_mat = Ncm.Spectral.get_x2_matrix(N)
        x2_np = x2_mat.to_numpy()
        gegen_result = x2_np @ cheb_x3

        # Evaluate at test points and compare with x^5
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            eval_result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_result, x)
            expected = x**5
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-12,
                atol=ATOL,
                err_msg=f"x^2*x^3 evaluation at x={x} doesn't match x^5",
            )

    def test_gaussian_operators(self, spectral: Ncm.Spectral) -> None:
        """Test operators on exp(-x^2), a non-polynomial function."""
        N = 64

        # Define f(x) = exp(-x^2)
        def f_gaussian(_user_data: None, x: float) -> float:
            return np.exp(-(x**2))

        # Compute Chebyshev coefficients for exp(-x^2)
        cheb_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_gaussian, -1.0, 1.0, N, None)
        )

        # Get operator matrices
        proj_mat = Ncm.Spectral.get_proj_matrix(N)
        x_mat = Ncm.Spectral.get_x_matrix(N)
        x2_mat = Ncm.Spectral.get_x2_matrix(N)
        d_mat = Ncm.Spectral.get_d_matrix(N)
        x_d_mat = Ncm.Spectral.get_x_d_matrix(N)
        d2_mat = Ncm.Spectral.get_d2_matrix(N)
        x_d2_mat = Ncm.Spectral.get_x_d2_matrix(N)
        x2_d2_mat = Ncm.Spectral.get_x2_d2_matrix(N)

        # Convert to numpy
        proj_np = proj_mat.to_numpy()
        x_np = x_mat.to_numpy()
        x2_np = x2_mat.to_numpy()
        d_np = d_mat.to_numpy()
        x_d_np = x_d_mat.to_numpy()
        d2_np = d2_mat.to_numpy()
        x_d2_np = x_d2_mat.to_numpy()
        x2_d2_np = x2_d2_mat.to_numpy()

        # Test points for evaluation
        test_points = np.array([0.0, 0.3, 0.5, 0.7, 0.9])

        # Test 1: Projection gives same function
        gegen_proj = proj_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_proj, x)
            expected = np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=0.0,
                err_msg=f"Projection of exp(-x^2) at x={x} failed",
            )

        # Test 2: x operator: x*exp(-x^2)
        gegen_x = x_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x, x)
            expected = x * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=ATOL,
                err_msg=f"x*exp(-x^2) at x={x} failed",
            )

        # Test 3: x^2 operator: x^2*exp(-x^2)
        gegen_x2 = x2_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x2, x)
            expected = x**2 * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=ATOL,
                err_msg=f"x^2*exp(-x^2) at x={x} failed",
            )

        # Test 4: d operator: d/dx(exp(-x^2)) = -2x*exp(-x^2)
        gegen_d = d_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_d, x)
            expected = -2.0 * x * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-7,
                atol=ATOL,
                err_msg=f"d/dx(exp(-x^2)) at x={x} failed",
            )

        # Test 5: d^2 operator: d^2/dx^2(exp(-x^2)) = (4x^2 - 2)*exp(-x^2)
        gegen_d2 = d2_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_d2, x)
            expected = (4.0 * x**2 - 2.0) * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=ATOL,
                err_msg=f"d^2/dx^2(exp(-x^2)) at x={x} failed",
            )

        # Test 6: x*d operator: x*d/dx(exp(-x^2)) = x*(-2x*exp(-x^2)) = -2x^2*exp(-x^2)
        gegen_x_d = x_d_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x_d, x)
            expected = -2.0 * x**2 * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-7,
                atol=ATOL,
                err_msg=f"x*d/dx(exp(-x^2)) at x={x} failed",
            )

        # Test 7: x*d^2 operator: x*d^2/dx^2(exp(-x^2)) = x*(4x^2 - 2)*exp(-x^2)
        gegen_x_d2 = x_d2_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x_d2, x)
            expected = x * (4.0 * x**2 - 2.0) * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=ATOL,
                err_msg=f"x*d^2/dx^2(exp(-x^2)) at x={x} failed",
            )

        # Test 8: x^2*d^2 operator: x^2*d^2/dx^2(exp(-x^2)) = x^2*(4x^2 - 2)*exp(-x^2)
        gegen_x2_d2 = x2_d2_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x2_d2, x)
            expected = x**2 * (4.0 * x**2 - 2.0) * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-10,
                atol=ATOL,
                err_msg=f"x^2*d^2/dx^2(exp(-x^2)) at x={x} failed",
            )

    def test_rational_operators(self, spectral: Ncm.Spectral) -> None:
        """Test operators on x^2/(1+x^2)^4, a rational function."""
        N = 64

        # Define f(x) = x^2/(1+x^2)^4
        def f_rational(_user_data: None, x: float) -> float:
            return x**2 / (1.0 + x**2) ** 4

        # Analytical derivative: f'(x) = 2x(1-3x^2)/(1+x^2)^5
        def f_prime(x: float) -> float:
            return 2.0 * x * (1.0 - 3.0 * x**2) / (1.0 + x**2) ** 5

        # Analytical second derivative: f''(x) = (2 - 36x^2 + 42x^4)/(1+x^2)^6
        def f_double_prime(x: float) -> float:
            return (2.0 - 36.0 * x**2 + 42.0 * x**4) / (1.0 + x**2) ** 6

        # Compute Chebyshev coefficients
        cheb_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_rational, -1.0, 1.0, N, None)
        )

        # Get operator matrices
        proj_mat = Ncm.Spectral.get_proj_matrix(N)
        x_mat = Ncm.Spectral.get_x_matrix(N)
        x2_mat = Ncm.Spectral.get_x2_matrix(N)
        d_mat = Ncm.Spectral.get_d_matrix(N)
        x_d_mat = Ncm.Spectral.get_x_d_matrix(N)
        d2_mat = Ncm.Spectral.get_d2_matrix(N)
        x_d2_mat = Ncm.Spectral.get_x_d2_matrix(N)
        x2_d2_mat = Ncm.Spectral.get_x2_d2_matrix(N)

        # Convert to numpy
        proj_np = proj_mat.to_numpy()
        x_np = x_mat.to_numpy()
        x2_np = x2_mat.to_numpy()
        d_np = d_mat.to_numpy()
        x_d_np = x_d_mat.to_numpy()
        d2_np = d2_mat.to_numpy()
        x_d2_np = x_d2_mat.to_numpy()
        x2_d2_np = x2_d2_mat.to_numpy()

        # Test points for evaluation
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

        # Test 1: Projection gives same function
        gegen_proj = proj_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_proj, x)
            expected = x**2 / (1.0 + x**2) ** 4
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=1.0e-10,
                err_msg=f"Projection of x^2/(1+x^2)^4 at x={x} failed",
            )

        # Test 2: x operator: x * x^2/(1+x^2)^4 = x^3/(1+x^2)^4
        gegen_x = x_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x, x)
            expected = x**3 / (1.0 + x**2) ** 4
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=1.0e-10,
                err_msg=f"x*f(x) at x={x} failed",
            )

        # Test 3: x^2 operator: x^2 * x^2/(1+x^2)^4 = x^4/(1+x^2)^4
        gegen_x2 = x2_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x2, x)
            expected = x**4 / (1.0 + x**2) ** 4
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=1.0e-10,
                err_msg=f"x^2*f(x) at x={x} failed",
            )

        # Test 4: d operator: d/dx of f(x)
        gegen_d = d_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_d, x)
            expected = f_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-7,
                atol=1.0e-9,
                err_msg=f"d/dx(f(x)) at x={x} failed",
            )

        # Test 5: d^2 operator: d^2/dx^2 of f(x)
        gegen_d2 = d2_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_d2, x)
            expected = f_double_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=1.0e-8,
                err_msg=f"d^2/dx^2(f(x)) at x={x} failed",
            )

        # Test 6: x*d operator: x*f'(x)
        gegen_x_d = x_d_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x_d, x)
            expected = x * f_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-7,
                atol=1.0e-9,
                err_msg=f"x*d/dx(f(x)) at x={x} failed",
            )

        # Test 7: x*d^2 operator: x*f''(x)
        gegen_x_d2 = x_d2_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x_d2, x)
            expected = x * f_double_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=1.0e-8,
                err_msg=f"x*d^2/dx^2(f(x)) at x={x} failed",
            )

        # Test 8: x^2*d^2 operator: x^2*f''(x)
        gegen_x2_d2 = x2_d2_np @ cheb_coeffs
        for x in test_points:
            result = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_x2_d2, x)
            expected = x**2 * f_double_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=1.0e-8,
                err_msg=f"x^2*d^2/dx^2(f(x)) at x={x} failed",
            )

    def test_x_to_t_conversion(self) -> None:
        """Test x to t conversion function."""
        a = -2.0
        b = 5.0

        # Test endpoints
        t_at_a = Ncm.Spectral.x_to_t(a, b, a)
        assert_allclose(t_at_a, -1.0, rtol=1e-14)

        t_at_b = Ncm.Spectral.x_to_t(a, b, b)
        assert_allclose(t_at_b, 1.0, rtol=1e-14)

        # Test midpoint
        x_mid = (a + b) / 2.0
        t_at_mid = Ncm.Spectral.x_to_t(a, b, x_mid)
        assert_allclose(t_at_mid, 0.0, rtol=1e-14)

        # Test several intermediate points
        x_values = np.linspace(a, b, 10)
        for x in x_values:
            t = Ncm.Spectral.x_to_t(a, b, x)
            assert -1.0 <= t <= 1.0

    def test_t_to_x_conversion(self) -> None:
        """Test t to x conversion function."""
        a = -2.0
        b = 5.0

        # Test endpoints
        x_at_minus1 = Ncm.Spectral.t_to_x(a, b, -1.0)
        assert_allclose(x_at_minus1, a, rtol=1e-14)

        x_at_plus1 = Ncm.Spectral.t_to_x(a, b, 1.0)
        assert_allclose(x_at_plus1, b, rtol=1e-14)

        # Test midpoint
        x_at_0 = Ncm.Spectral.t_to_x(a, b, 0.0)
        assert_allclose(x_at_0, (a + b) / 2.0, rtol=1e-14)

        # Test several intermediate points
        t_values = np.linspace(-1.0, 1.0, 10)
        for t in t_values:
            x = Ncm.Spectral.t_to_x(a, b, t)
            assert a <= x <= b

    def test_conversion_round_trip(self) -> None:
        """Test that x -> t -> x conversion is identity."""
        a = -3.5
        b = 7.2

        x_values = np.linspace(a, b, 20)
        for x in x_values:
            t = Ncm.Spectral.x_to_t(a, b, x)
            x_back = Ncm.Spectral.t_to_x(a, b, t)
            assert_allclose(x_back, x, rtol=1e-14)

        # Also test t -> x -> t
        t_values = np.linspace(-1.0, 1.0, 20)
        for t in t_values:
            x = Ncm.Spectral.t_to_x(a, b, t)
            t_back = Ncm.Spectral.x_to_t(a, b, x)
            assert_allclose(t_back, t, rtol=1e-14)

    def test_gegenbauer_alpha1_eval_x(self) -> None:
        """Test Gegenbauer alpha=1 evaluation with x in [a, b]."""
        a = 2.0
        b = 8.0

        # Set up coefficients for C_0^(1) + C_1^(1)
        c = np.zeros(8)
        c[0] = 1.0
        c[1] = 1.0

        # Test at several points in [a, b]
        x_values = np.linspace(a, b, 10)
        for x in x_values:
            # Convert x to t
            t = Ncm.Spectral.x_to_t(a, b, x)

            # Evaluate using both methods
            result_t = Ncm.Spectral.gegenbauer_alpha1_eval(c, t)
            result_x = Ncm.Spectral.gegenbauer_alpha1_eval_x(c, a, b, x)

            # They should match
            assert_allclose(result_x, result_t, rtol=1e-14)

    def test_gegenbauer_alpha2_eval_x(self) -> None:
        """Test Gegenbauer alpha=2 evaluation with x in [a, b]."""
        a = -1.0
        b = 3.0

        # Set up coefficients for C_0^(2) + C_1^(2)
        c = np.zeros(8)
        c[0] = 1.0
        c[1] = 1.0

        # Test at several points in [a, b]
        x_values = np.linspace(a, b, 10)
        for x in x_values:
            # Convert x to t
            t = Ncm.Spectral.x_to_t(a, b, x)

            # Evaluate using both methods
            result_t = Ncm.Spectral.gegenbauer_alpha2_eval(c, t)
            result_x = Ncm.Spectral.gegenbauer_alpha2_eval_x(c, a, b, x)

            # They should match
            assert_allclose(result_x, result_t, rtol=1e-14)

    def test_chebyshev_eval_x(self) -> None:
        """Test Chebyshev evaluation with x in [a, b]."""
        a = 0.5
        b = 4.5

        # Set up Chebyshev coefficients for a simple function
        coeffs = np.array([1.0, 2.0, 1.5, 0.5, 0.2])

        # Test at several points in [a, b]
        x_values = np.linspace(a, b, 15)
        for x in x_values:
            # Convert x to t
            t = Ncm.Spectral.x_to_t(a, b, x)

            # Evaluate using both methods
            result_t = Ncm.Spectral.chebyshev_eval(coeffs, t)
            result_x = Ncm.Spectral.chebyshev_eval_x(coeffs, a, b, x)

            # They should match
            assert_allclose(result_x, result_t, rtol=1e-14)

    def test_chebyshev_deriv_x(self, spectral: Ncm.Spectral) -> None:
        """Test Chebyshev derivative evaluation with x in [a, b]."""
        a = 1.0
        b = 5.0

        # Define a function f(x) = x^3 - 2*x^2 + x + 3
        def f(_user_data, x):
            return x**3 - 2.0 * x**2 + x + 3.0

        # Its derivative is f'(x) = 3*x^2 - 4*x + 1
        def f_prime(x):
            return 3.0 * x**2 - 4.0 * x + 1.0

        # Compute Chebyshev coefficients on [a, b]
        coeffs = spectral.compute_chebyshev_coeffs_adaptive(f, a, b, 3, 1e-12, None)

        # Test derivative at several points
        x_values = np.linspace(a + 0.1, b - 0.1, 10)
        for x in x_values:
            # Compute derivative using the _x function
            df_dx = Ncm.Spectral.chebyshev_deriv_x(coeffs, a, b, x)

            # Compare with analytical derivative
            expected = f_prime(x)
            assert_allclose(
                df_dx, expected, rtol=1e-10, err_msg=f"Derivative mismatch at x={x}"
            )

    def test_chebyshev_coeffs_on_interval(self, spectral: Ncm.Spectral) -> None:
        """Test computing Chebyshev coefficients and evaluating on a custom interval."""
        a = -5.0
        b = 3.0

        # Define f(x) = sin(x) on [a, b]
        def f(_user_data, x):
            return np.sin(x)

        # Compute coefficients
        coeffs = spectral.compute_chebyshev_coeffs_adaptive(f, a, b, 4, 1e-10, None)

        # Evaluate at several points using both t and x methods
        x_test = np.linspace(a, b, 20)
        for x in x_test:
            # Using _x method
            result_x = Ncm.Spectral.chebyshev_eval_x(coeffs, a, b, x)

            # Using manual conversion
            t = Ncm.Spectral.x_to_t(a, b, x)
            result_t = Ncm.Spectral.chebyshev_eval(coeffs, t)

            # Direct evaluation
            expected = f(None, x)

            assert_allclose(result_x, result_t, rtol=1e-14)
            assert_allclose(result_x, expected, rtol=1e-8, atol=1e-10)
