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


class TestSpectral:
    """Tests for NcmSpectral spectral methods."""

    @staticmethod
    def vector_to_numpy(ncm_vec: Ncm.Vector) -> np.ndarray:
        """Convert NcmVector to numpy array."""
        return np.array(ncm_vec.dup_array())

    @pytest.fixture
    def solver(self) -> Ncm.Spectral:
        """Create an ODE solver with l=1."""
        return Ncm.Spectral.new()

    def test_chebT_to_gegenbauer_alpha1_basic(self) -> None:
        """Test Chebyshev to Gegenbauer alpha=1 conversion."""
        N = 16

        # T_0 = C_0^(1)
        c = Ncm.Vector.new(N)
        c.set_zero()
        c.set(0, 1.0)
        g = Ncm.Vector.new(N)
        g.set_zero()

        Ncm.Spectral.chebT_to_gegenbauer_alpha1(c, g)
        g_np = self.vector_to_numpy(g)

        expected = np.zeros(N)
        expected[0] = 1.0
        assert_allclose(g_np, expected, rtol=1.0e-14, atol=1.0e-14)

        # T_1 = (1/2) C_1^(1)
        c = Ncm.Vector.new(N)
        c.set_zero()
        c.set(1, 1.0)
        g = Ncm.Vector.new(N)
        g.set_zero()

        Ncm.Spectral.chebT_to_gegenbauer_alpha1(c, g)
        g_np = self.vector_to_numpy(g)

        expected = np.zeros(N)
        expected[1] = 0.5
        assert_allclose(g_np, expected, rtol=1.0e-14, atol=1.0e-14)

    def test_chebT_to_gegenbauer_alpha2_basic(self) -> None:
        """Test Chebyshev to Gegenbauer alpha=2 conversion."""
        N = 16

        # Test T_0 conversion
        c = Ncm.Vector.new(N)
        c.set_zero()
        c.set(0, 1.0)
        g = Ncm.Vector.new(N)
        g.set_zero()

        Ncm.Spectral.chebT_to_gegenbauer_alpha2(c, g)
        g_np = self.vector_to_numpy(g)

        # For T_0, g_0 should have contribution 1/2 * c_0 + c_0/(2*1) = 1.0
        assert g_np[0] > 0.9

    def test_gegenbauer_alpha1_eval(self) -> None:
        """Test Gegenbauer alpha=1 evaluation."""
        N = 8

        # C_0^(1)(x) = 1, C_1^(1)(x) = 2x
        c = Ncm.Vector.new(N)
        c.set_zero()
        c.set(0, 1.0)
        c.set(1, 1.0)

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
        c = Ncm.Vector.new(N)
        c.set_zero()
        c.set(0, 1.0)
        c.set(1, 1.0)

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
        a = Ncm.Vector.new(N)
        a.set_zero()
        a.set(0, 1.0)

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_eval(a, x)
            assert_allclose(result, 1.0, rtol=1.0e-14)

        # T_1(x) = x
        a = Ncm.Vector.new(N)
        a.set_zero()
        a.set(1, 1.0)

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_eval(a, x)
            assert_allclose(result, x, rtol=1.0e-14)

        # T_2(x) = 2x^2 - 1
        a = Ncm.Vector.new(N)
        a.set_zero()
        a.set(2, 1.0)

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_eval(a, x)
            expected = 2.0 * x * x - 1.0
            assert_allclose(result, expected, rtol=1.0e-14)

    def test_chebyshev_deriv_basic(self) -> None:
        """Test Chebyshev derivative with simple polynomials."""
        N = 8

        # Derivative of constant is 0
        a = Ncm.Vector.new(N)
        a.set_zero()
        a.set(0, 1.0)

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_deriv(a, x)
            assert_allclose(result, 0.0, rtol=1.0e-14, atol=1.0e-14)

        # Derivative of T_1(x) = x is 1
        a = Ncm.Vector.new(N)
        a.set_zero()
        a.set(1, 1.0)

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_deriv(a, x)
            assert_allclose(result, 1.0, rtol=1.0e-14)

        # Derivative of T_2(x) = 2x^2 - 1 is 4x
        a = Ncm.Vector.new(N)
        a.set_zero()
        a.set(2, 1.0)

        for x in [0.0, 0.5, 1.0, -1.0]:
            result = Ncm.Spectral.chebyshev_deriv(a, x)
            expected = 4.0 * x
            assert_allclose(result, expected, rtol=1.0e-14)

    def test_compute_chebyshev_coeffs_constant(self) -> None:
        """Test computing Chebyshev coefficients for constant function."""
        N = 16
        spectral = Ncm.Spectral.new()

        def f_constant(_user_data, _x):
            return 1.0

        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_constant, -1.0, 1.0, coeffs, None)
        coeffs_np = self.vector_to_numpy(coeffs)

        # f(x) = 1 = T_0(x), so c_0 = 1, rest = 0
        assert_allclose(coeffs_np[0], 1.0, rtol=1.0e-12)
        assert_allclose(coeffs_np[1:], 0.0, rtol=1.0e-12, atol=1.0e-12)

    def test_compute_chebyshev_coeffs_linear(self) -> None:
        """Test computing Chebyshev coefficients for linear function."""
        N = 16
        spectral = Ncm.Spectral.new()

        def f_linear(_user_data, x):
            return x

        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_linear, -1.0, 1.0, coeffs, None)
        coeffs_np = self.vector_to_numpy(coeffs)

        # f(x) = x = T_1(x), so c_1 = 1, rest = 0
        assert_allclose(coeffs_np[0], 0.0, rtol=1.0e-12, atol=1.0e-12)
        assert_allclose(coeffs_np[1], 1.0, rtol=1.0e-12)
        assert_allclose(coeffs_np[2:], 0.0, rtol=1.0e-12, atol=1.0e-12)

    def test_compute_chebyshev_coeffs_quadratic(self) -> None:
        """Test computing Chebyshev coefficients for quadratic function."""
        N = 16
        spectral = Ncm.Spectral.new()

        def f_quadratic(_user_data, x):
            return x * x

        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_quadratic, -1.0, 1.0, coeffs, None)
        coeffs_np = self.vector_to_numpy(coeffs)

        # f(x) = x^2 = (T_0 + T_2)/2, so c_0 = 1/2, c_2 = 1/2
        assert_allclose(coeffs_np[0], 0.5, rtol=1.0e-11)
        assert_allclose(coeffs_np[1], 0.0, rtol=1.0e-12, atol=1.0e-12)
        assert_allclose(coeffs_np[2], 0.5, rtol=1.0e-11)
        assert_allclose(coeffs_np[3:], 0.0, rtol=1.0e-12, atol=1.0e-12)

    def test_compute_chebyshev_coeffs_gaussian(self) -> None:
        """Test computing Chebyshev coefficients for Gaussian function."""
        N = 64
        spectral = Ncm.Spectral.new()

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_gaussian, -1.0, 1.0, coeffs, None)

        # Verify by evaluating at test points
        test_points = np.array([0.0, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs, x)
            expected = np.exp(-x * x)
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-10,
                err_msg=f"Gaussian evaluation at x={x} doesn't match",
            )

    def test_gaussian_derivative(self) -> None:
        """Test derivative of Gaussian function."""
        N = 64
        spectral = Ncm.Spectral.new()

        def f_gaussian(_user_data, x):
            return np.exp(-x * x)

        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_gaussian, -1.0, 1.0, coeffs, None)

        # Derivative of exp(-x^2) is -2x*exp(-x^2)
        test_points = np.array([0.0, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            deriv_result = Ncm.Spectral.chebyshev_deriv(coeffs, x)
            expected = -2.0 * x * np.exp(-x * x)
            assert_allclose(
                deriv_result,
                expected,
                rtol=1.0e-9,
                atol=1.0e-15,
                err_msg=f"Gaussian derivative at x={x} doesn't match",
            )

    def test_rational_function(self) -> None:
        """Test with rational function x^2/(1+x^2)^4."""
        N = 64
        spectral = Ncm.Spectral.new()

        def f_rational(_user_data, x):
            return (x * x) / ((1.0 + x * x) ** 4)

        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_rational, -1.0, 1.0, coeffs, None)

        # Verify by evaluating at test points
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs, x)
            expected = (x * x) / ((1.0 + x * x) ** 4)
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-10,
                atol=1.0e-15,
                err_msg=f"Rational function evaluation at x={x} doesn't match",
            )

    def test_rational_derivative(self) -> None:
        """Test derivative of rational function."""
        N = 64
        spectral = Ncm.Spectral.new()

        def f_rational(_user_data, x):
            return (x * x) / ((1.0 + x * x) ** 4)

        # Analytical derivative: f'(x) = 2x(1-3x^2)/(1+x^2)^5
        def f_prime(x):
            return 2.0 * x * (1.0 - 3.0 * x * x) / ((1.0 + x * x) ** 5)

        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_rational, -1.0, 1.0, coeffs, None)

        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            deriv_result = Ncm.Spectral.chebyshev_deriv(coeffs, x)
            expected = f_prime(x)
            assert_allclose(
                deriv_result,
                expected,
                rtol=1.0e-9,
                atol=1.0e-15,
                err_msg=f"Rational derivative at x={x} doesn't match",
            )

    def test_chebyshev_eval_endpoints(self) -> None:
        """Test that Chebyshev evaluation handles endpoints correctly."""
        N = 32
        spectral = Ncm.Spectral.new()

        # Use a simple smooth function
        def f_test(_user_data, x):
            return 1.0 + x + x * x

        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_test, -1.0, 1.0, coeffs, None)

        # Test at endpoints
        result_m1 = Ncm.Spectral.chebyshev_eval(coeffs, -1.0)
        expected_m1 = 1.0 - 1.0 + 1.0
        assert_allclose(result_m1, expected_m1, rtol=1.0e-12)

        result_p1 = Ncm.Spectral.chebyshev_eval(coeffs, 1.0)
        expected_p1 = 1.0 + 1.0 + 1.0
        assert_allclose(result_p1, expected_p1, rtol=1.0e-12)

    def test_chebyshev_deriv_endpoints(self) -> None:
        """Test that Chebyshev derivative handles endpoints correctly."""
        N = 32
        spectral = Ncm.Spectral.new()

        # Use f(x) = x^3, f'(x) = 3x^2
        def f_cubic(_user_data, x):
            return x * x * x

        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_cubic, -1.0, 1.0, coeffs, None)

        # Test derivative at endpoints
        result_m1 = Ncm.Spectral.chebyshev_deriv(coeffs, -1.0)
        expected_m1 = 3.0 * 1.0
        assert_allclose(result_m1, expected_m1, rtol=1.0e-11)

        result_p1 = Ncm.Spectral.chebyshev_deriv(coeffs, 1.0)
        expected_p1 = 3.0 * 1.0
        assert_allclose(result_p1, expected_p1, rtol=1.0e-11)

    def test_different_intervals(self) -> None:
        """Test computing Chebyshev coefficients on different intervals."""
        N = 32
        spectral = Ncm.Spectral.new()

        def f_test(_user_data, x):
            return x * x

        # Test on [0, 1]
        coeffs = Ncm.Vector.new(N)
        spectral.compute_chebyshev_coeffs(f_test, 0.0, 1.0, coeffs, None)

        # Map test points from [-1, 1] to [0, 1]
        for xi in [-0.5, 0.0, 0.5]:
            # xi in [-1, 1] maps to x in [0, 1]
            x = 0.5 + 0.5 * xi
            eval_result = Ncm.Spectral.chebyshev_eval(coeffs, xi)
            expected = x * x
            assert_allclose(eval_result, expected, rtol=1.0e-11)

    def test_spectral_cache(self) -> None:
        """Test that spectral caching works correctly."""
        spectral = Ncm.Spectral.new()

        def f_test(_user_data, x):
            return np.sin(x)

        # Compute with N=32
        coeffs1 = Ncm.Vector.new(32)
        spectral.compute_chebyshev_coeffs(f_test, -1.0, 1.0, coeffs1, None)

        # Compute again with N=32 (should use cache)
        coeffs2 = Ncm.Vector.new(32)
        spectral.compute_chebyshev_coeffs(f_test, -1.0, 1.0, coeffs2, None)

        # Compute with different N=64 (should reallocate)
        coeffs3 = Ncm.Vector.new(64)
        spectral.compute_chebyshev_coeffs(f_test, -1.0, 1.0, coeffs3, None)

        # Verify results are correct
        for coeffs, N in [(coeffs1, 32), (coeffs2, 32), (coeffs3, 64)]:
            assert coeffs.len() == N
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
        self, solver: Ncm.Spectral, power: int, N: int
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
        coeffs_vec = Ncm.Vector.new(N)
        solver.compute_chebyshev_coeffs(f_power, -1.0, 1.0, coeffs_vec, None)
        coeffs = np.array([coeffs_vec.get(i) for i in range(N)])

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
    def test_chebyshev_coeffs_exp(self, solver: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev coefficient computation for exp(x).

        Tests with exponential function which has known coefficients
        in terms of modified Bessel functions.
        """

        # Define function f(x) = exp(x)
        def f_exp(_user_data: None, x: float) -> float:
            return np.exp(x)

        # Compute Chebyshev coefficients
        coeffs_vec = Ncm.Vector.new(N)
        solver.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, coeffs_vec, None)
        coeffs = np.array([coeffs_vec.get(i) for i in range(N)])

        # Get analytical coefficients
        analytical = self.get_analytical_exp_coeffs(N)

        # Compare - exponential has exponentially decaying coefficients
        for i in range(min(20, N)):  # Check first 20 coefficients
            assert_allclose(
                coeffs[i],
                analytical[i],
                rtol=1.0e-10,
                atol=1.0e-14,
                err_msg=f"Mismatch at coefficient {i} for exp(x) with N={N}",
            )

    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_chebyshev_coeffs_cos(self, solver: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev coefficient computation for cos(x).

        Tests with cosine function which has known coefficients
        in terms of Bessel functions.
        """

        # Define function f(x) = cos(x)
        def f_cos(_user_data: None, x: float) -> float:
            return np.cos(x)

        # Compute Chebyshev coefficients
        coeffs_vec = Ncm.Vector.new(N)
        solver.compute_chebyshev_coeffs(f_cos, -1.0, 1.0, coeffs_vec, None)
        coeffs = np.array([coeffs_vec.get(i) for i in range(N)])

        # Get analytical coefficients
        analytical = self.get_analytical_cos_coeffs(N)

        # Compare - cosine has exponentially decaying coefficients
        for i in range(min(20, N)):  # Check first 20 coefficients
            assert_allclose(
                coeffs[i],
                analytical[i],
                rtol=1.0e-10,
                atol=1.0e-14,
                err_msg=f"Mismatch at coefficient {i} for cos(x) with N={N}",
            )

    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_chebyshev_coeffs_rational(self, solver: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev coefficient computation for 1/(2-x).

        Tests with rational function which has a singularity outside [-1,1].
        """

        # Define function f(x) = 1/(2-x)
        def f_rational(_user_data: None, x: float) -> float:
            return 1.0 / (2.0 - x)

        # Compute Chebyshev coefficients
        coeffs_vec = Ncm.Vector.new(N)
        solver.compute_chebyshev_coeffs(f_rational, -1.0, 1.0, coeffs_vec, None)
        coeffs = np.array([coeffs_vec.get(i) for i in range(N)])

        # Get analytical coefficients (computed numerically with high precision)
        analytical = self.get_analytical_rational_coeffs(N)

        # Compare
        for i in range(min(20, N)):
            assert_allclose(
                coeffs[i],
                analytical[i],
                rtol=1.0e-9,
                atol=1.0e-14,
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
        coeffs_vec = Ncm.Vector.new_array(coeffs.tolist())

        # Test at multiple points
        x_points = np.linspace(-1.0, 1.0, 1024)

        for x in x_points:
            # Evaluate using our implementation
            result = Ncm.Spectral.chebyshev_eval(coeffs_vec, x)

            # Evaluate using numpy
            expected = chebval(x, coeffs)

            assert_allclose(
                result,
                expected,
                rtol=1.0e-12,
                atol=1.0e-14,
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
        coeffs_vec = Ncm.Vector.new_array(coeffs.tolist())

        # Test derivative at multiple points
        x_points = np.linspace(-0.95, 0.95, 19)

        for x in x_points:
            # Evaluate derivative using our implementation
            result = Ncm.Spectral.chebyshev_deriv(coeffs_vec, x)

            # Expected derivative: d/dx(x^power) = power * x^(power-1)
            expected = power * x ** (power - 1) if power > 0 else 0.0

            assert_allclose(
                result,
                expected,
                rtol=1.0e-10,
                atol=1.0e-12,
                err_msg=f"Chebyshev deriv mismatch for x^{power} at x={x} with N={N}",
            )

    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_exp(self, solver: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev derivative for exp(x).

        The derivative of exp(x) is exp(x) itself.
        """

        # Define function f(x) = exp(x)
        def f_exp(_user_data: None, x: float) -> float:
            return np.exp(x)

        # Compute Chebyshev coefficients
        coeffs_vec = Ncm.Vector.new(N)
        solver.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, coeffs_vec, None)

        # Test derivative at multiple points
        x_points = np.linspace(-0.95, 0.95, 19)

        for x in x_points:
            # Evaluate derivative
            result = Ncm.Spectral.chebyshev_deriv(coeffs_vec, x)

            # Expected: d/dx(exp(x)) = exp(x)
            expected = np.exp(x)

            assert_allclose(
                result,
                expected,
                rtol=1.0e-9,
                atol=1.0e-12,
                err_msg=f"Chebyshev deriv mismatch for exp(x) at x={x} with N={N}",
            )

    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_sin(self, solver: Ncm.Spectral, N: int) -> None:
        """
        Test Chebyshev derivative for sin(x).

        The derivative of sin(x) is cos(x).
        """

        # Define function f(x) = sin(x)
        def f_sin(_user_data: None, x: float) -> float:
            return np.sin(x)

        # Compute Chebyshev coefficients
        coeffs_vec = Ncm.Vector.new(N)
        solver.compute_chebyshev_coeffs(f_sin, -1.0, 1.0, coeffs_vec, None)

        # Test derivative at multiple points
        x_points = np.linspace(-0.95, 0.95, 19)

        for x in x_points:
            # Evaluate derivative
            result = Ncm.Spectral.chebyshev_deriv(coeffs_vec, x)

            # Expected: d/dx(sin(x)) = cos(x)
            expected = np.cos(x)

            assert_allclose(
                result,
                expected,
                rtol=1.0e-9,
                atol=1.0e-12,
                err_msg=f"Chebyshev deriv mismatch for sin(x) at x={x} with N={N}",
            )

    def test_chebyshev_deriv_edge_cases(self) -> None:
        """
        Test Chebyshev derivative edge cases.

        Tests single coefficients and constant functions.
        """
        # Test single coefficient (constant function)
        # f(x) = 5.0, so f'(x) = 0
        const_vec = Ncm.Vector.new(1)
        const_vec.set(0, 5.0)
        result = Ncm.Spectral.chebyshev_deriv(const_vec, 0.5)
        assert result == 0.0, "Derivative of constant should be 0"

        result = Ncm.Spectral.chebyshev_deriv(const_vec, 0.0)
        assert result == 0.0, "Derivative of constant should be 0 at x=0"

        result = Ncm.Spectral.chebyshev_deriv(const_vec, 1.0)
        assert result == 0.0, "Derivative of constant should be 0 at x=1"

        result = Ncm.Spectral.chebyshev_deriv(const_vec, -1.0)
        assert result == 0.0, "Derivative of constant should be 0 at x=-1"

        # Test two coefficients (linear function: a0 + a1*x)
        # f(x) = 2.0 + 3.0*x, so f'(x) = 3.0
        linear_vec = Ncm.Vector.new(2)
        linear_vec.set(0, 2.0)
        linear_vec.set(1, 3.0)

        x_test = 0.5
        result = Ncm.Spectral.chebyshev_deriv(linear_vec, x_test)
        assert_allclose(result, 3.0, rtol=1.0e-14, atol=1.0e-14)

        # Test at x=0
        result = Ncm.Spectral.chebyshev_deriv(linear_vec, 0.0)
        assert_allclose(result, 3.0, rtol=1.0e-14, atol=1.0e-14)

        # Test at x=1
        result = Ncm.Spectral.chebyshev_deriv(linear_vec, 1.0)
        assert_allclose(result, 3.0, rtol=1.0e-14, atol=1.0e-14)

        # Test at x=-1
        result = Ncm.Spectral.chebyshev_deriv(linear_vec, -1.0)
        assert_allclose(result, 3.0, rtol=1.0e-14, atol=1.0e-14)

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
        coeffs_vec = Ncm.Vector.new_array(coeffs.tolist())

        # Test at x = 1: f'(1) = power * 1^(power-1) = power
        result_p1 = Ncm.Spectral.chebyshev_deriv(coeffs_vec, 1.0)
        expected_p1 = power if power > 0 else 0.0
        assert_allclose(
            result_p1,
            expected_p1,
            rtol=1.0e-10,
            atol=1.0e-12,
            err_msg=f"Derivative of x^{power} at x=1 failed",
        )

        # Test at x = -1: f'(-1) = power * (-1)^(power-1)
        result_m1 = Ncm.Spectral.chebyshev_deriv(coeffs_vec, -1.0)
        expected_m1 = power * ((-1.0) ** (power - 1)) if power > 0 else 0.0
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=1.0e-10,
            atol=1.0e-12,
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
        coeffs_vec = Ncm.Vector.new_array(coeffs.tolist())

        x_points = np.linspace(-0.9, 0.9, 17)

        for x in x_points:
            result = Ncm.Spectral.chebyshev_deriv(coeffs_vec, x)

            theta = np.arccos(x)
            expected = n * np.sin(n * theta) / np.sin(theta)

            assert_allclose(
                result,
                expected,
                rtol=1.0e-10,
                atol=1.0e-12,
                err_msg=f"Derivative of T_{n} failed at x={x}",
            )

    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_alternating_coeffs(self, N: int) -> None:
        """
        Test derivative with alternating high-frequency coefficients.
        """
        coeffs = np.array([(-1.0) ** k for k in range(N)], dtype=float)
        coeffs_vec = Ncm.Vector.new_array(coeffs.tolist())

        x_points = np.linspace(-0.8, 0.8, 11)
        h = 1.0e-7

        for x in x_points:
            d_numeric = (
                Ncm.Spectral.chebyshev_eval(coeffs_vec, x + h)
                - Ncm.Spectral.chebyshev_eval(coeffs_vec, x - h)
            ) / (2.0 * h)

            d_cheb = Ncm.Spectral.chebyshev_deriv(coeffs_vec, x)

            assert_allclose(
                d_cheb,
                d_numeric,
                rtol=1.0e-6,
                atol=1.0e-8,
                err_msg=f"Alternating coeff derivative mismatch at x={x}",
            )

    @pytest.mark.parametrize("N", [32, 64])
    def test_chebyshev_deriv_random_coeffs(self, N: int) -> None:
        """
        Test derivative for random Chebyshev coefficients using finite differences.
        """
        rng = np.random.default_rng(1234)
        coeffs = rng.normal(size=N)
        coeffs_vec = Ncm.Vector.new_array(coeffs.tolist())

        x_points = np.linspace(-0.7, 0.7, 9)
        h = 1.0e-7

        for x in x_points:
            f_plus = Ncm.Spectral.chebyshev_eval(coeffs_vec, x + h)
            f_minus = Ncm.Spectral.chebyshev_eval(coeffs_vec, x - h)
            d_numeric = (f_plus - f_minus) / (2.0 * h)

            d_cheb = Ncm.Spectral.chebyshev_deriv(coeffs_vec, x)

            assert_allclose(
                d_cheb,
                d_numeric,
                rtol=1.0e-6,
                atol=1.0e-8,
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
        coeffs_vec = Ncm.Vector.new_array(coeffs.tolist())

        # x = 1
        result_p1 = Ncm.Spectral.chebyshev_deriv(coeffs_vec, 1.0)
        assert_allclose(
            result_p1,
            n * n,
            rtol=0.0,
            atol=1.0e-12,
            err_msg=f"T_{n}'(1) failed",
        )

        # x = -1
        result_m1 = Ncm.Spectral.chebyshev_deriv(coeffs_vec, -1.0)
        expected_m1 = ((-1.0) ** (n - 1)) * n * n
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=0.0,
            atol=1.0e-12,
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
        cheb_vec = Ncm.Vector.new_array(cheb_coeffs.tolist())
        gegen_vec = Ncm.Vector.new(N)

        # Convert to Gegenbauer alpha=1 (U_n Chebyshev of second kind)
        Ncm.Spectral.chebT_to_gegenbauer_alpha1(cheb_vec, gegen_vec)

        # Test evaluation at multiple points
        x_points = np.linspace(-0.9, 0.9, 11)

        for x in x_points:
            # Evaluate using Chebyshev basis
            cheb_val = Ncm.Spectral.chebyshev_eval(cheb_vec, x)

            # Evaluate using Gegenbauer alpha=1 basis
            gegen_val = Ncm.Spectral.gegenbauer_alpha1_eval(gegen_vec, x)

            # Both should give x^power
            expected = x**power

            assert_allclose(
                cheb_val,
                expected,
                rtol=1.0e-10,
                atol=1.0e-12,
                err_msg=f"Chebyshev eval failed for x^{power} at x={x}",
            )

            assert_allclose(
                gegen_val,
                expected,
                rtol=1.0e-10,
                atol=1.0e-12,
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
        cheb_vec = Ncm.Vector.new_array(cheb_coeffs.tolist())
        gegen_vec = Ncm.Vector.new(N)

        # Convert to Gegenbauer alpha=2
        Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_vec, gegen_vec)

        # Test evaluation at multiple points
        x_points = np.linspace(-0.9, 0.9, 11)

        for x in x_points:
            # Evaluate using Chebyshev basis
            cheb_val = Ncm.Spectral.chebyshev_eval(cheb_vec, x)

            # Evaluate using Gegenbauer alpha=2 basis
            gegen_val = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_vec, x)

            # Both should give x^power
            expected = x**power

            assert_allclose(
                cheb_val,
                expected,
                rtol=1.0e-10,
                atol=1.0e-12,
                err_msg=f"Chebyshev eval failed for x^{power} at x={x}",
            )

            assert_allclose(
                gegen_val,
                expected,
                rtol=1.0e-10,
                atol=1.0e-12,
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
            coeffs_vec = Ncm.Vector.new_array(coeffs.tolist())

            for x in x_points:
                # Evaluate using our implementation
                if alpha_val == 1:
                    result = Ncm.Spectral.gegenbauer_alpha1_eval(coeffs_vec, x)
                else:
                    result = Ncm.Spectral.gegenbauer_alpha2_eval(coeffs_vec, x)

                # Evaluate using scipy
                expected = eval_gegenbauer(n, alpha_val, x)

                assert_allclose(
                    result,
                    expected,
                    rtol=1.0e-10,
                    atol=1.0e-12,
                    err_msg=f"Gegenbauer C_{n}^({alpha_val}) mismatch at x={x}",
                )

    @pytest.mark.parametrize("N", [32, 64])
    def test_gegenbauer_roundtrip(self, solver: Ncm.Spectral, N: int) -> None:
        """
        Test round-trip conversion: function -> Chebyshev -> Gegenbauer -> evaluation.

        Uses exp(x) as test function to verify the entire pipeline.
        """

        def f_exp(_user_data: None, x: float) -> float:
            return np.exp(x)

        # Compute Chebyshev coefficients
        coeffs_vec = Ncm.Vector.new(N)
        solver.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, coeffs_vec, None)

        # Convert to Gegenbauer alpha=1
        gegen1_vec = Ncm.Vector.new(N)
        Ncm.Spectral.chebT_to_gegenbauer_alpha1(coeffs_vec, gegen1_vec)

        # Convert to Gegenbauer alpha=2
        gegen2_vec = Ncm.Vector.new(N)
        Ncm.Spectral.chebT_to_gegenbauer_alpha2(coeffs_vec, gegen2_vec)

        # Test evaluation at multiple points
        x_points = np.linspace(-0.8, 0.8, 17)

        for x in x_points:
            expected = np.exp(x)

            # Evaluate with Chebyshev
            cheb_val = Ncm.Spectral.chebyshev_eval(coeffs_vec, x)

            # Evaluate with Gegenbauer alpha=1
            gegen1_val = Ncm.Spectral.gegenbauer_alpha1_eval(gegen1_vec, x)

            # Evaluate with Gegenbauer alpha=2
            gegen2_val = Ncm.Spectral.gegenbauer_alpha2_eval(gegen2_vec, x)

            # All should match exp(x)
            assert_allclose(
                cheb_val,
                expected,
                rtol=1.0e-8,
                atol=1.0e-10,
                err_msg=f"Chebyshev roundtrip failed for exp(x) at x={x}",
            )

            assert_allclose(
                gegen1_val,
                expected,
                rtol=1.0e-8,
                atol=1.0e-10,
                err_msg=f"Gegenbauer alpha=1 roundtrip failed for exp(x) at x={x}",
            )

            assert_allclose(
                gegen2_val,
                expected,
                rtol=1.0e-8,
                atol=1.0e-10,
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
        cheb_vec = Ncm.Vector.new_array(cheb_coeffs.tolist())

        # Test at x = 1
        result_p1 = Ncm.Spectral.chebyshev_eval(cheb_vec, 1.0)
        expected_p1 = 1.0**power
        assert_allclose(
            result_p1,
            expected_p1,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg=f"Chebyshev eval of x^{power} at x=1 failed",
        )

        # Test at x = -1
        result_m1 = Ncm.Spectral.chebyshev_eval(cheb_vec, -1.0)
        expected_m1 = (-1.0) ** power
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=1.0e-12,
            atol=1.0e-14,
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
        cheb_vec = Ncm.Vector.new_array(cheb_coeffs.tolist())
        gegen_vec = Ncm.Vector.new(N)

        # Convert to Gegenbauer alpha=1
        Ncm.Spectral.chebT_to_gegenbauer_alpha1(cheb_vec, gegen_vec)

        # Test at x = 1
        result_p1 = Ncm.Spectral.gegenbauer_alpha1_eval(gegen_vec, 1.0)
        expected_p1 = 1.0**power
        assert_allclose(
            result_p1,
            expected_p1,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg=f"Gegenbauer alpha=1 eval of x^{power} at x=1 failed",
        )

        # Test at x = -1
        result_m1 = Ncm.Spectral.gegenbauer_alpha1_eval(gegen_vec, -1.0)
        expected_m1 = (-1.0) ** power
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=1.0e-12,
            atol=1.0e-14,
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
        cheb_vec = Ncm.Vector.new_array(cheb_coeffs.tolist())
        gegen_vec = Ncm.Vector.new(N)

        # Convert to Gegenbauer alpha=2
        Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb_vec, gegen_vec)

        # Test at x = 1
        result_p1 = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_vec, 1.0)
        expected_p1 = 1.0**power
        assert_allclose(
            result_p1,
            expected_p1,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg=f"Gegenbauer alpha=2 eval of x^{power} at x=1 failed",
        )

        # Test at x = -1
        result_m1 = Ncm.Spectral.gegenbauer_alpha2_eval(gegen_vec, -1.0)
        expected_m1 = (-1.0) ** power
        assert_allclose(
            result_m1,
            expected_m1,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg=f"Gegenbauer alpha=2 eval of x^{power} at x=-1 failed",
        )
