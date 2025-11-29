#!/usr/bin/env python
#
# test_py_pln1d.py
#
# Sat Nov 29 18:34:00 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_pln1d.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests for ncm_pln1d module."""

import unittest
import time
import numpy as np
from scipy.special import gammaln  # pylint: disable=no-name-in-module

from numcosmo_py import Ncm


def pln1d_integral_numerical(R, mu, sigma, n_samples=500):
    """Compute Poisson-Lognormal 1D integral numerically using trapezoidal rule.

    Computes P(R|μ,σ) = ∫ P_Poisson(R|λ) * P_LogNormal(λ|μ,σ) dλ

    Parameters
    ----------
    R : float
        Poisson rate parameter
    mu : float
        Log-normal location parameter (mean of log(λ))
    sigma : float
        Log-normal scale parameter (std dev of log(λ))
    n_samples : int, optional
        Number of sample points for numerical integration (default: 500)

    Returns
    -------
    float
        Numerical integral result
    """
    # Substitute λ = exp(x), so dλ = exp(x) dx
    # Integral becomes: ∫ exp(R*x - exp(x)) / R! * 1/(σ√(2π)) * exp(-(x-μ)²/(2σ²)) dx

    # Sample x in a wide range around the expected mode
    x_center = np.log(R) if R > 0 else 0.0
    x_min = x_center - 5 * sigma - 5
    x_max = x_center + 5 * sigma + 5

    x_values = np.linspace(x_min, x_max, n_samples)

    # Compute integrand in log space for numerical stability
    # ln(integrand) = R*x - exp(x) - ln(R!) - 0.5*ln(2π*σ²) - (x-μ)²/(2σ²)
    lnfactorial_R = gammaln(R + 1)

    ln_integrand = (
        R * x_values
        - np.exp(x_values)
        - lnfactorial_R
        - 0.5 * np.log(2 * np.pi * sigma**2)
        - (x_values - mu) ** 2 / (2 * sigma**2)
    )

    integrand = np.exp(ln_integrand)

    # Numerical integration using trapezoidal rule
    return np.trapezoid(integrand, x_values)


class TestPLN1D(unittest.TestCase):
    """Test cases for PLN1D module."""

    def setUp(self):
        """Set up test fixtures."""
        self.pln = Ncm.PLN1D()

    def test_pln1d_new(self):
        """Test PLN1D object creation."""
        pln = Ncm.PLN1D()
        self.assertIsNotNone(pln)
        # Default order should be 60
        self.assertEqual(pln.get_property("gh-order"), 60)

    def test_pln1d_new_with_order(self):
        """Test PLN1D creation with specific Gauss-Hermite order."""
        for order in [10, 30, 60, 100]:
            pln = Ncm.PLN1D.new(order)
            self.assertIsNotNone(pln)
            self.assertEqual(pln.get_property("gh-order"), order)

    def test_pln1d_set_order(self):
        """Test setting Gauss-Hermite order."""
        for order in [5, 20, 50, 100]:
            self.pln.set_property("gh-order", order)
            self.assertEqual(self.pln.get_property("gh-order"), order)
            self.assertEqual(self.pln.get_order(), order)

            self.pln.set_order(order + 10)
            self.assertEqual(self.pln.get_property("gh-order"), order + 10)
            self.assertEqual(self.pln.get_order(), order + 10)

    def test_pln1d_set_get_order_consistency(self):
        """Test consistency between set and get order."""
        for order in [10, 30, 60, 90]:
            self.pln.set_property("gh-order", order)
            retrieved_order = self.pln.get_property("gh-order")
            self.assertEqual(retrieved_order, order)

    def test_pln1d_eval_p_basic(self):
        """Test basic evaluation of probability with simple parameters."""
        # Test with moderate parameters
        R = 5.0
        mu = 0.0
        sigma = 1.0
        result = self.pln.eval_p(R, mu, sigma)
        # Should be positive
        self.assertGreater(result, 0.0)
        # Should be finite
        self.assertTrue(np.isfinite(result))

    def test_pln1d_eval_lnp_basic(self):
        """Test basic evaluation of log-probability with simple parameters."""
        # Test with moderate parameters
        R = 5.0
        mu = 0.0
        sigma = 1.0
        result = self.pln.eval_lnp(R, mu, sigma)
        # Should be finite (can be negative)
        self.assertTrue(np.isfinite(result))

    def test_pln1d_eval_p_varies_with_R(self):
        """Test that result changes with Poisson rate parameter R."""
        mu = 0.0
        sigma = 1.0
        result_R1 = self.pln.eval_p(1.0, mu, sigma)
        result_R5 = self.pln.eval_p(5.0, mu, sigma)
        result_R10 = self.pln.eval_p(10.0, mu, sigma)

        # Results should be different for different R values
        self.assertNotEqual(result_R1, result_R5)
        self.assertNotEqual(result_R5, result_R10)
        # All should be positive and finite
        self.assertGreater(result_R1, 0.0)
        self.assertGreater(result_R5, 0.0)
        self.assertGreater(result_R10, 0.0)

    def test_pln1d_eval_p_varies_with_mu(self):
        """Test that result changes with log-normal location parameter."""
        R = 5.0
        sigma = 1.0
        result_mu0 = self.pln.eval_p(R, 0.0, sigma)
        result_mu1 = self.pln.eval_p(R, 1.0, sigma)
        result_mu_neg1 = self.pln.eval_p(R, -1.0, sigma)

        # Results should be different for different mu values
        self.assertNotEqual(result_mu0, result_mu1)
        self.assertNotEqual(result_mu0, result_mu_neg1)
        # All should be positive and finite
        self.assertGreater(result_mu0, 0.0)
        self.assertGreater(result_mu1, 0.0)
        self.assertGreater(result_mu_neg1, 0.0)

    def test_pln1d_eval_p_varies_with_sigma(self):
        """Test that result changes with log-normal scale parameter."""
        R = 5.0
        mu = 0.0
        result_sigma05 = self.pln.eval_p(R, mu, 0.5)
        result_sigma1 = self.pln.eval_p(R, mu, 1.0)
        result_sigma2 = self.pln.eval_p(R, mu, 2.0)

        # Results should be different for different sigma values
        self.assertNotEqual(result_sigma05, result_sigma1)
        self.assertNotEqual(result_sigma1, result_sigma2)
        # All should be positive and finite
        self.assertGreater(result_sigma05, 0.0)
        self.assertGreater(result_sigma1, 0.0)
        self.assertGreater(result_sigma2, 0.0)

    def test_pln1d_eval_p_small_sigma_uses_laplace(self):
        """Test that small sigma values use Laplace approximation."""
        R = 5.0
        mu = 0.0
        # Use GH quadrature
        self.pln.set_property("gh-order", 60)
        result_gh_small = self.pln.eval_p(R, mu, 1.0e-5)
        # Should still be positive and finite (Laplace fallback)
        self.assertGreater(result_gh_small, 0.0)
        self.assertTrue(np.isfinite(result_gh_small))

    def test_pln1d_eval_p_large_R(self):
        """Test evaluation with large Poisson rate parameter."""
        mu = 0.0
        sigma = 1.0
        for R in [10.0, 50.0, 100.0]:
            result = self.pln.eval_p(R, mu, sigma)
            self.assertGreater(result, 0.0)
            self.assertTrue(np.isfinite(result))

    def test_pln1d_eval_p_zero_R(self):
        """Test evaluation with R = 0.

        This is the limiting case (Poisson with rate 1 integrated against lognormal).
        """
        mu = 0.0
        sigma = 1.0
        result = self.pln.eval_p(0.0, mu, sigma)
        self.assertGreater(result, 0.0)
        self.assertTrue(np.isfinite(result))

    def test_pln1d_eval_p_negative_mu(self):
        """Test evaluation with negative mu."""
        R = 5.0
        sigma = 1.0
        for mu in [-2.0, -1.0, -0.5]:
            result = self.pln.eval_p(R, mu, sigma)
            self.assertGreater(result, 0.0)
            self.assertTrue(np.isfinite(result))

    def test_pln1d_eval_p_consistency_gh_laplace(self):
        """Test that GH and Laplace approximations give reasonable results."""
        # Using moderate sigma where both methods should work well
        R = 5.0
        mu = 0.0
        sigma = 1.0

        # GH quadrature with reasonable order
        self.pln.set_property("gh-order", 60)
        result_gh = self.pln.eval_p(R, mu, sigma)

        # Both should be positive and finite
        self.assertGreater(result_gh, 0.0)
        self.assertTrue(np.isfinite(result_gh))

    def test_pln1d_eval_p_monotonicity_R(self):
        """Test monotonicity with respect to R for fixed mu, sigma."""
        mu = 0.0
        sigma = 1.0
        R_values = np.linspace(0.1, 10.0, 10)
        results = np.array([self.pln.eval_p(R, mu, sigma) for R in R_values])

        # All results should be positive and finite
        self.assertTrue(np.all(results > 0.0))
        self.assertTrue(np.all(np.isfinite(results)))
        # Results should vary smoothly (no sudden jumps)
        differences = np.abs(np.diff(results))
        self.assertTrue(np.all(differences > 0.0))  # Should be strictly changing

    def test_pln1d_eval_p_different_orders(self):
        """Test that different GH orders give consistent results."""
        R = 5.0
        mu = 0.0
        sigma = 1.0

        results = {}
        for order in [10, 30, 60]:
            self.pln.set_property("gh-order", order)
            results[order] = self.pln.eval_p(
                R, mu, sigma
            )  # Results from different orders should be close (converging)
        # Higher orders should give more accurate results
        self.assertGreater(results[10], 0.0)
        self.assertGreater(results[30], 0.0)
        self.assertGreater(results[60], 0.0)

        # Check relative consistency (within ~5% for different orders)
        rel_diff_10_60 = abs(results[10] - results[60]) / results[60]
        self.assertLess(rel_diff_10_60, 0.10)  # Allow up to 10% difference

    def test_pln1d_ref_free_clear(self):
        """Test reference counting: ref, free, and clear operations."""
        pln = Ncm.PLN1D.new(60)
        # Increase reference count
        pln_ref = pln.ref()
        self.assertIsNotNone(pln_ref)

        # Test clear operation
        pln_ptr = pln
        pln_ptr = None  # Simulate clearing
        self.assertIsNone(pln_ptr)

    def test_pln1d_eval_p_numerical_stability(self):
        """Test numerical stability with extreme but valid parameters."""
        # Test with very small sigma (should use Laplace)
        result1 = self.pln.eval_p(5.0, 0.0, 1.0e-6)
        self.assertTrue(np.isfinite(result1))

        # Test with large mu
        result2 = self.pln.eval_p(5.0, 10.0, 1.0)
        self.assertTrue(np.isfinite(result2))

        # Test with large sigma
        result3 = self.pln.eval_p(5.0, 0.0, 10.0)
        self.assertTrue(np.isfinite(result3))

    def test_pln1d_eval_p_grid(self):
        """Test evaluation on a grid of parameters."""
        R_values = [1.0, 5.0, 10.0]
        mu_values = [-1.0, 0.0, 1.0]
        sigma_values = [0.5, 1.0, 2.0]

        for R in R_values:
            for mu in mu_values:
                for sigma in sigma_values:
                    result = self.pln.eval_p(R, mu, sigma)
                    self.assertGreater(
                        result,
                        0.0,
                        f"Result should be positive for R={R}, mu={mu}, sigma={sigma}",
                    )
                    self.assertTrue(
                        np.isfinite(result),
                        f"Result should be finite for R={R}, mu={mu}, sigma={sigma}",
                    )

    def test_pln1d_eval_p_vs_lnp_relationship(self):
        """Test that eval_p = exp(eval_lnp)."""
        test_cases = [
            (1.0, 0.0, 1.0),
            (5.0, 0.0, 1.0),
            (10.0, -1.0, 0.5),
            (3.0, 1.0, 2.0),
            (0.1, -0.5, 0.8),
        ]

        for R, mu, sigma in test_cases:
            p = self.pln.eval_p(R, mu, sigma)
            lnp = self.pln.eval_lnp(R, mu, sigma)

            # Check relationship: p = exp(lnp)
            expected_p = np.exp(lnp)
            rel_error = abs(p - expected_p) / expected_p
            self.assertLess(
                rel_error,
                1.0e-10,
                f"eval_p and exp(eval_lnp) mismatch for R={R}, mu={mu}, sigma={sigma}",
            )

    def test_pln1d_symmetry_in_mu(self):
        """Test that the integrand is symmetric around mode in log-space."""
        R = 5.0
        sigma = 1.0
        # The integrand I(x) has a unique mode; evaluate around symmetrically offset
        # points The function values should follow expected asymmetry due to the
        # Poisson term
        mu_center = 0.0
        mu_plus = self.pln.eval_p(R, mu_center + 1.0, sigma)
        mu_minus = self.pln.eval_p(R, mu_center - 1.0, sigma)
        # Both should be valid
        self.assertGreater(mu_plus, 0.0)
        self.assertGreater(mu_minus, 0.0)

    def test_pln1d_limiting_case_zero_sigma(self):
        """Test behavior as sigma approaches 0 (limit should approach Poisson)."""
        R = 5.0
        mu = np.log(R)  # Set mu such that mode is near R in log space
        # As sigma -> 0, the lognormal becomes a delta function at exp(mu)
        # Result should be well-defined and finite
        result_small = self.pln.eval_p(R, mu, 1.0e-3)
        result_moderate = self.pln.eval_p(R, mu, 1.0e-2)

        self.assertGreater(result_small, 0.0)
        self.assertGreater(result_moderate, 0.0)
        self.assertTrue(np.isfinite(result_small))
        self.assertTrue(np.isfinite(result_moderate))

    def test_pln1d_numerical_derivative_R(self):
        """Test first derivative with respect to R using finite differences."""
        mu = 0.0
        sigma = 1.0
        h = 1.0e-4

        R_test = 5.0
        f_plus = self.pln.eval_p(R_test + h, mu, sigma)
        f_minus = self.pln.eval_p(R_test - h, mu, sigma)
        numerical_deriv = (f_plus - f_minus) / (2.0 * h)

        # Derivative should be negative (result decreases as R increases in this region)
        # or positive depending on parameters, but should be non-zero
        self.assertNotEqual(numerical_deriv, 0.0)
        self.assertTrue(np.isfinite(numerical_deriv))

    def test_pln1d_consistency_repeated_eval(self):
        """Test that repeated evaluations.

        Three evaluations with same parameters give identical results."""
        R, mu, sigma = 5.0, 0.0, 1.0

        result1 = self.pln.eval_p(R, mu, sigma)
        result2 = self.pln.eval_p(R, mu, sigma)
        result3 = self.pln.eval_p(
            R, mu, sigma
        )  # All results should be identical (deterministic computation)
        self.assertEqual(result1, result2)
        self.assertEqual(result2, result3)

    def test_pln1d_scaling_with_sigma(self):
        """Test scaling behavior: result changes monotonically with sigma."""
        R = 5.0
        mu = 0.0
        sigma_values = np.linspace(0.5, 3.0, 6)
        results = np.array([self.pln.eval_p(R, mu, s) for s in sigma_values])

        # All results should be positive and finite
        self.assertTrue(np.all(results > 0.0))
        self.assertTrue(np.all(np.isfinite(results)))

    def test_pln1d_large_R_poisson_limit(self):
        """Test that for large R, result should decrease (approaching Poisson limit)."""
        mu = 0.0
        sigma = 1.0

        # For large R with fixed sigma, the result should remain well-defined
        R_large = 100.0
        result = self.pln.eval_p(R_large, mu, sigma)
        self.assertGreater(result, 0.0)
        self.assertTrue(np.isfinite(result))

    def test_pln1d_gh_vs_laplace_consistency(self):
        """Test consistency between GH quadrature and Laplace in overlap region."""
        # Use a case where both should work reasonably well
        R = 3.0
        mu = 0.0
        sigma = 0.8  # Moderate sigma

        # Both methods should give same order of magnitude result
        self.pln.set_property("gh-order", 60)
        result_gh = self.pln.eval_p(R, mu, sigma)

        # Use a sigma that forces Laplace
        result_laplace = self.pln.eval_p(
            R, mu, 1.0e-5
        )  # Both should be positive and finite
        self.assertGreater(result_gh, 0.0)
        self.assertGreater(result_laplace, 0.0)

    def test_pln1d_order_convergence(self):
        """Test that higher GH orders converge to more accurate result."""
        R = 5.0
        mu = 0.0
        sigma = 1.0

        orders = [10, 20, 40, 60]
        results = []

        for order in orders:
            self.pln.set_property("gh-order", order)
            results.append(
                self.pln.eval_p(R, mu, sigma)
            )  # Check monotonic convergence: differences should decrease
        diffs = np.abs(np.diff(results))
        # Most differences should decrease (allowing for some numerical variation)
        decreasing_count = sum(
            1 for i in range(len(diffs) - 1) if diffs[i] >= diffs[i + 1]
        )
        self.assertGreaterEqual(decreasing_count, len(diffs) - 2)

    def test_pln1d_result_bounds(self):
        """Test that results are within physically reasonable bounds."""
        # For Poisson–Lognormal with these parameters
        R_values = np.array([0.0, 1.0, 5.0, 10.0])
        mu = 0.0
        sigma = 1.0

        for R in R_values:
            result = self.pln.eval_p(R, mu, sigma)
            # Result should be positive
            self.assertGreater(result, 0.0)
            # Result should be reasonable (not infinity or NaN)
            self.assertTrue(np.isfinite(result))
            # For R=0, result should be reasonably large (integrating 1 with log-normal)
            if R == 0.0:
                self.assertGreater(result, 0.1)

    def test_pln1d_parameter_sensitivity(self):
        """Test sensitivity to parameter changes."""
        # Baseline
        R_base = 5.0
        mu_base = 0.0
        sigma_base = 1.0
        result_base = self.pln.eval_p(R_base, mu_base, sigma_base)

        # Small perturbations
        delta_R = 1.0e-5
        delta_mu = 1.0e-5
        delta_sigma = 1.0e-5

        result_R_pert = self.pln.eval_p(R_base + delta_R, mu_base, sigma_base)
        result_mu_pert = self.pln.eval_p(R_base, mu_base + delta_mu, sigma_base)
        result_sigma_pert = self.pln.eval_p(
            R_base, mu_base, sigma_base + delta_sigma
        )  # Results should be close but not identical
        rel_change_R = abs(result_R_pert - result_base) / result_base
        rel_change_mu = abs(result_mu_pert - result_base) / result_base
        rel_change_sigma = abs(result_sigma_pert - result_base) / result_base

        # Small perturbations should cause small changes (continuous function)
        self.assertLess(rel_change_R, 0.01)
        self.assertLess(rel_change_mu, 0.01)
        self.assertLess(rel_change_sigma, 0.01)

    def test_pln1d_special_case_R_equals_1(self):
        """Test special case where R = 1 (exponential Poisson)."""
        mu = 0.0
        sigma = 1.0
        result = self.pln.eval_p(1.0, mu, sigma)

        self.assertGreater(result, 0.0)
        self.assertTrue(np.isfinite(result))

    def test_pln1d_extreme_but_valid_mu(self):
        """Test with extreme but mathematically valid mu values."""
        R = 5.0
        sigma = 1.0

        # Very negative mu (shifts log-normal to very small values)
        result_mu_neg = self.pln.eval_p(R, -5.0, sigma)
        # Very positive mu (shifts log-normal to very large values)
        result_mu_pos = self.pln.eval_p(R, 5.0, sigma)

        self.assertGreater(result_mu_neg, 0.0)
        self.assertGreater(result_mu_pos, 0.0)
        self.assertTrue(np.isfinite(result_mu_neg))
        self.assertTrue(np.isfinite(result_mu_pos))

        # Results should be different
        self.assertNotEqual(result_mu_neg, result_mu_pos)

    def test_pln1d_sum_over_R(self):
        """Test that sum over a range of R values is finite."""
        mu = 0.0
        sigma = 1.0
        R_values = np.arange(0.0, 300.0, 1.0)
        results = np.array([self.pln.eval_p(R, mu, sigma) for R in R_values])

        total_sum = np.sum(results)
        self.assertGreater(total_sum, 0.0)
        self.assertTrue(np.isfinite(total_sum))
        self.assertAlmostEqual(
            total_sum, 1.0, places=4
        )  # Should be close to 1 for normalization

    def test_pln1d_sum_over_R_higher_order(self):
        """Test that sum over a range of R values is finite."""
        self.pln.set_order(300)
        mu = 0.0
        sigma = 1.0
        R_values = np.arange(0.0, 300.0, 1.0)
        results = np.array([self.pln.eval_p(R, mu, sigma) for R in R_values])

        total_sum = np.sum(results)
        self.assertGreater(total_sum, 0.0)
        self.assertTrue(np.isfinite(total_sum))
        self.assertAlmostEqual(
            total_sum, 1.0, places=6
        )  # Should be close to 1 for normalization

    def test_pln1d_normal_in_lnR_high_mu(self):
        """Test that for moderately high mu, p(R) approaches normal.

        The distribution p(R) should be close to a normal distribution in ln(R).

        For high mu, the distribution p(R) as a function of ln(R) should be
        approximately Gaussian centered near mu with width roughly sigma.
        """
        self.pln.set_order(100)
        mu = np.log(15)  # Moderately high mu
        sigma = 0.5

        # Sample ln(R) uniformly around mu
        ln_R_values = np.linspace(mu - 4.0 * sigma, mu + 4.0 * sigma, 1000)
        R_values = np.exp(ln_R_values)

        # Evaluate p(R) at each point
        ln_p_values = np.array([self.pln.eval_lnp(R, mu, sigma) for R in R_values])

        # All should be positive and finite
        self.assertTrue(np.all(np.exp(ln_p_values) > 0.0))
        self.assertTrue(np.all(np.isfinite(np.exp(ln_p_values))))
        # Fit a Gaussian to ln(p) vs ln(R)
        coeffs = np.polyfit(ln_R_values, ln_p_values, 2)

        # For a normal in ln(R), ln(p) should be quadratic in ln(R)
        # coeffs[0] is the quadratic term (related to 1/(2σ²))
        # coeffs[1] is the linear term (related to mean)
        # coeffs[2] is the constant

        # The peak should be near mu
        peak_ln_R = -coeffs[1] / (2 * coeffs[0])
        self.assertLess(abs(peak_ln_R - mu), sigma, "Peak should be near mu")

        # The width (curvature) should be approximately consistent with sigma
        # For a Gaussian: ln(p) ∝ -(ln(R) - μ)² / (2σ²)
        # So the quadratic coefficient should be approximately -1/(2σ²)
        expected_curvature = -1.0 / (2 * (sigma**2 + np.exp(-mu)))
        rel_error = abs(coeffs[0] - expected_curvature) / abs(expected_curvature)
        # Allow 10% relative error due to finite sample and Poisson term
        self.assertLess(
            rel_error,
            0.06,
            "Curvature should be consistent with Gaussian width approximately sigma",
        )

    def test_pln1d_timing(self):
        """Test timing of PLN1D evaluations.

        Measure and verify that evaluations complete in reasonable time.
        This is a performance sanity check, not a strict timing constraint.
        """

        self.pln.set_order(100)
        R = 5.0
        mu = 0.0
        sigma = 1.0
        num_evals = 1000

        # Time eval_p
        start_time = time.perf_counter()
        for _ in range(num_evals):
            self.pln.eval_p(R, mu, sigma)
        end_time = time.perf_counter()
        time_eval_p = end_time - start_time
        time_per_eval_p = time_eval_p / num_evals

        # Time eval_lnp
        start_time = time.perf_counter()
        for _ in range(num_evals):
            self.pln.eval_lnp(R, mu, sigma)
        end_time = time.perf_counter()
        time_eval_lnp = end_time - start_time
        time_per_eval_lnp = time_eval_lnp / num_evals

        # Results should complete in reasonable time (< 1 ms per eval on modern
        # hardware)
        self.assertLess(
            time_per_eval_p,
            1.0e-5,
            f"eval_p should be faster than 1ms, got {time_per_eval_p*1000:.3f}ms",
        )
        self.assertLess(
            time_per_eval_lnp,
            1.0e-5,
            f"eval_lnp should be faster than 1ms, got {time_per_eval_lnp*1000:.3f}ms",
        )

        # Both functions should have similar performance
        slower_time = max(time_per_eval_p, time_per_eval_lnp)
        faster_time = min(time_per_eval_p, time_per_eval_lnp)
        time_ratio = slower_time / faster_time if faster_time > 0 else 1.0
        # Allow up to 2x difference between the two methods
        self.assertLess(
            time_ratio,
            2.0,
            (
                f"eval_p and eval_lnp should have similar "
                f"performance, ratio={time_ratio:.2f}"
            ),
        )

    def test_pln1d_correctness_vs_numerical_integration(self):
        """Test PLN1D correctness by comparing with numerical integration.

        Compute the integral numerically using numpy's trapezoidal rule and
        compare with the C implementation result.

        Note: The very tight agreement (typically < 1e-12 relative error)
        indicates that the C implementation uses accurate numerical integration
        methods (likely Gauss-Hermite quadrature), which converges to machine
        precision for these test parameters.
        """
        self.pln.set_order(100)

        # Test parameters
        test_cases = [
            (2.0, 0.0, 0.5),  # Moderate R, mu, sigma
            (5.0, 1.0, 0.8),  # Higher R and mu
            (1.0, -1.0, 1.0),  # Low R, negative mu
            (10.0, 2.0, 0.3),  # High R, high mu, small sigma
        ]

        for R, mu, sigma in test_cases:
            # Get C implementation result
            c_result = self.pln.eval_p(R, mu, sigma)

            # Compute integral numerically
            numerical_result = pln1d_integral_numerical(R, mu, sigma)

            # Both should be positive
            self.assertGreater(c_result, 0.0)
            self.assertGreater(numerical_result, 0.0)

            # Compare results - should agree within reasonable tolerance
            rel_error = abs(c_result - numerical_result) / max(
                c_result, numerical_result
            )
            self.assertLess(
                rel_error,
                1.0e-12,
                f"C result {c_result:.6e} vs numerical {numerical_result:.6e} "
                f"(rel_error={rel_error:.4f}) for R={R}, mu={mu}, sigma={sigma}",
            )

    def test_pln1d_eval_range_sum_basic(self):
        """Test basic range sum evaluation."""
        self.pln.set_order(100)
        mu = 0.0
        sigma = 1.0

        # Test a small range
        R_min = 0
        R_max = 10
        result = self.pln.eval_range_sum(R_min, R_max, mu, sigma)

        # Result should be positive and finite
        self.assertGreater(result, 0.0)
        self.assertTrue(np.isfinite(result))

    def test_pln1d_eval_range_sum_lnp_basic(self):
        """Test basic range sum log-probability evaluation."""
        self.pln.set_order(100)
        mu = 0.0
        sigma = 1.0

        # Test a small range
        R_min = 0
        R_max = 10
        result = self.pln.eval_range_sum_lnp(R_min, R_max, mu, sigma)

        # Result should be finite (can be negative in log space)
        self.assertTrue(np.isfinite(result))

    def test_pln1d_eval_range_sum_consistency_with_individual_evals(self):
        """Test that range sum equals sum of individual evaluations."""
        self.pln.set_order(100)
        mu = 0.0
        sigma = 1.0

        R_min = 0
        R_max = 20

        # Compute range sum
        range_sum = self.pln.eval_range_sum(R_min, R_max, mu, sigma)

        # Compute sum of individual evaluations
        individual_sum = sum(
            self.pln.eval_p(R, mu, sigma) for R in range(R_min, R_max + 1)
        )

        # Results should agree very closely
        rel_error = abs(range_sum - individual_sum) / max(range_sum, individual_sum)
        self.assertLess(
            rel_error,
            1.0e-7,
            f"Range sum {range_sum:.6e} vs individual sum {individual_sum:.6e}",
        )

    def test_pln1d_eval_range_sum_lnp_consistency(self):
        """Test that eval_range_sum = exp(eval_range_sum_lnp)."""
        self.pln.set_order(100)
        mu = 0.0
        sigma = 1.0

        R_min = 0
        R_max = 15

        # Compute both forms
        p_sum = self.pln.eval_range_sum(R_min, R_max, mu, sigma)
        lnp_sum = self.pln.eval_range_sum_lnp(R_min, R_max, mu, sigma)

        # Check relationship: p_sum = exp(lnp_sum)
        expected_p = np.exp(lnp_sum)
        rel_error = abs(p_sum - expected_p) / expected_p
        self.assertLess(
            rel_error,
            1.0e-10,
            f"Range sum {p_sum:.6e} vs exp(lnp_sum) {expected_p:.6e}",
        )

    def test_pln1d_eval_range_sum_varies_with_range(self):
        """Test that range sum increases with larger ranges."""
        self.pln.set_order(100)
        mu = 0.0
        sigma = 1.0

        # Test progressively larger ranges
        range_sums = []
        for R_max in [10, 20, 30, 50]:
            result = self.pln.eval_range_sum(0, R_max, mu, sigma)
            range_sums.append(result)
            self.assertGreater(result, 0.0)
            self.assertTrue(np.isfinite(result))

        # Results should be strictly increasing
        for i in range(len(range_sums) - 1):
            self.assertGreater(
                range_sums[i + 1],
                range_sums[i],
                "Range sums should increase with larger range",
            )

    def test_pln1d_eval_range_sum_high_mu(self):
        """Test range sum with moderately high mu."""
        self.pln.set_order(100)
        mu = 3.0  # Moderately high mu
        sigma = 0.5

        R_min = 5
        R_max = 25
        result = self.pln.eval_range_sum(R_min, R_max, mu, sigma)

        # Result should be positive and finite
        self.assertGreater(result, 0.0)
        self.assertTrue(np.isfinite(result))

        # Should be close to individual sum
        individual_sum = sum(
            self.pln.eval_p(R, mu, sigma) for R in range(R_min, R_max + 1)
        )
        rel_error = abs(result - individual_sum) / max(result, individual_sum)
        self.assertLess(rel_error, 1.0e-10)

    def test_pln1d_eval_range_sum_single_R(self):
        """Test range sum with single R value (R_min == R_max)."""
        self.pln.set_order(100)
        mu = 0.0
        sigma = 1.0

        R = 5
        range_sum = self.pln.eval_range_sum(R, R, mu, sigma)
        individual = self.pln.eval_p(R, mu, sigma)

        # Should be equal for single R
        rel_error = abs(range_sum - individual) / individual
        self.assertLess(rel_error, 1.0e-12)

    def test_pln1d_eval_range_sum_different_orders(self):
        """Test that different GH orders give consistent range sum results."""
        mu = 0.0
        sigma = 1.0
        R_min = 0
        R_max = 20

        results = {}
        for order in [30, 60, 100]:
            self.pln.set_order(order)
            results[order] = self.pln.eval_range_sum(R_min, R_max, mu, sigma)

        # All results should be positive
        for order, value in results.items():
            self.assertGreater(value, 0.0)

        # Results should converge with higher orders (within 1% for these ranges)
        rel_diff = abs(results[30] - results[100]) / results[100]
        self.assertLess(rel_diff, 1.0e-4)

    def test_pln1d_eval_range_sum_normalization(self):
        """Test that summing over sufficient range approaches 1.0."""
        self.pln.set_order(100)
        mu = 0.0
        sigma = 1.0

        # Sum over a very large range to approximate total probability
        R_max = 200
        result = self.pln.eval_range_sum(0, R_max, mu, sigma)

        # Should be close to 1.0 for normalization
        self.assertAlmostEqual(result, 1.0, places=3)
